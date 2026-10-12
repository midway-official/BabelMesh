#include "babelmesh/cad.hpp"
#include "babelmesh/mesh_conversion.hpp"
#include "babelmesh/mesh_optimizer.hpp"
#include <Interface_Static.hxx>
#include <STEPControl_Controller.hxx>
#include <STEPControl_Writer.hxx>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
namespace babelmesh {
namespace {
using Json = nlohmann::json;
std::string stamp() {
  return std::to_string(
      std::chrono::steady_clock::now().time_since_epoch().count());
}
void jsonFile(const std::filesystem::path &path, const Json &value) {
  std::ofstream stream(path);
  stream << value.dump(2) << '\n';
  stream.flush();
  if (!stream)
    throw std::runtime_error("cannot write conversion report");
}
void bundle(const MixedMesh &m, const std::filesystem::path &directory,
            bool overwrite) {
  const std::vector<std::string> names = {"mesh.bsm", "mesh.vtu", "mesh.json",
                                          "mesh_report.json"};
  for (const auto &name : names)
    if (std::filesystem::exists(directory / name) && !overwrite)
      throw std::runtime_error("output exists; use --overwrite");
  std::filesystem::create_directories(directory);
  const auto temp = directory / (".staging-" + stamp());
  std::filesystem::create_directory(temp);
  try {
    writeMixedBabelSimV3(m, temp / "mesh.bsm");
    writeMixedVtu(m, temp / "mesh.vtu");
    writeMixedCheckpoint(m, temp / "mesh.json");
    jsonFile(temp / "mesh_report.json", mixedMeshReport(m));
    for (const auto &name : names)
      std::filesystem::rename(temp / name, directory / name);
    std::filesystem::remove(temp);
  } catch (...) {
    std::filesystem::remove_all(temp);
    throw;
  }
}
} // namespace
int runConvertMesh(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: babelmesh convert INPUT --to FORMAT --out FILE "
                 "[--from FORMAT] [--length-scale FACTOR] [--overwrite] "
                 "[--checkpoint FILE]\nformats: babelsim-v3, checkpoint, "
                 "fluent-ascii, cgns, cgns-poly, vtu\n";
    return 2;
  }
  try {
    std::string from = "auto", to;
    std::filesystem::path output, checkpoint;
    bool overwrite = false;
    double lengthScale = 1;
    for (int i = 2; i < argc; ++i) {
      const std::string flag = argv[i];
      if ((flag == "--from" || flag == "--to" || flag == "--format" ||
           flag == "--out" || flag == "--length-scale" ||
           flag == "--checkpoint") &&
          i + 1 < argc) {
        const std::string value = argv[++i];
        if (flag == "--from")
          from = value;
        else if (flag == "--to" || flag == "--format")
          to = value;
        else if (flag == "--out")
          output = value;
        else if (flag == "--checkpoint")
          checkpoint = value;
        else
          lengthScale = std::stod(value);
      } else if (flag == "--overwrite")
        overwrite = true;
      else
        throw std::runtime_error("unknown conversion argument: " + flag);
    }
    if (output.empty() || to.empty() || !std::isfinite(lengthScale) ||
        lengthScale <= 0)
      throw std::runtime_error(
          "conversion requires --to, --out and positive finite length scale");
    const auto reportFile =
        std::filesystem::path(output.string() + ".report.json");
    std::set<std::filesystem::path> files = {
        std::filesystem::weakly_canonical(std::filesystem::absolute(argv[1]))};
    for (const auto &file : {output, reportFile, checkpoint}) {
      if (file.empty())
        continue;
      if (!files
               .insert(std::filesystem::weakly_canonical(
                   std::filesystem::absolute(file)))
               .second)
        throw std::runtime_error("conversion input, output, checkpoint and "
                                 "report must be distinct files");
      if (std::filesystem::exists(file) && !overwrite)
        throw std::runtime_error("output exists; use --overwrite: " +
                                 file.string());
    }
    auto mesh = readExchangeMesh(argv[1], from);
    for (auto &point : mesh.points_m)
      for (auto &c : point)
        c *= lengthScale;
    mesh.extrusion_thickness_m *= lengthScale;
    if (mesh.metadata.contains("domain_bounds_constraint")) {
      auto &bounds = mesh.metadata["domain_bounds_constraint"];
      for (const auto *key : {"minimum_m", "maximum_m"})
        for (auto &coordinate : bounds.at(key))
          coordinate = coordinate.get<double>() * lengthScale;
      if (lengthScale != 1)
        bounds["explicit_length_scale"] =
            bounds.value("explicit_length_scale", 1.0) * lengthScale;
    }
    validateMixedMesh(mesh);
    const auto token = stamp();
    const auto temporary =
        std::filesystem::path(output.string() + ".staging-" + token);
    const auto reportTemp =
        std::filesystem::path(reportFile.string() + ".staging-" + token);
    const auto checkpointTemp =
        std::filesystem::path(checkpoint.string() + ".staging-" + token);
    std::filesystem::create_directories(
        std::filesystem::absolute(output).parent_path());
    try {
      writeExchangeMesh(mesh, temporary, to);
      auto report = mixedMeshReport(mesh);
      report["conversion"] = {
          {"source", std::filesystem::absolute(argv[1]).string()},
          {"destination", std::filesystem::absolute(output).string()},
          {"destination_format", to},
          {"length_scale", lengthScale},
          {"added_cells", 0},
          {"commercial_application_import_tested", false},
          {"scope", "volume mesh geometry/connectivity and boundary groups; "
                    "solver settings and fields are not converted"}};
      jsonFile(reportTemp, report);
      if (!checkpoint.empty()) {
        std::filesystem::create_directories(
            std::filesystem::absolute(checkpoint).parent_path());
        writeMixedCheckpoint(mesh, checkpointTemp);
      }
      std::filesystem::rename(temporary, output);
      std::filesystem::rename(reportTemp, reportFile);
      if (!checkpoint.empty())
        std::filesystem::rename(checkpointTemp, checkpoint);
    } catch (...) {
      for (const auto &file : {temporary, reportTemp, checkpointTemp})
        std::filesystem::remove(file);
      throw;
    }
    std::cout << "Converted " << mesh.cells.size()
              << " cells, 0 added: " << output << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "mesh conversion failed: " << e.what() << '\n';
    return 1;
  }
}
int runOptimizeMesh(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: babelmesh optimize INPUT --out DIRECTORY "
                 "[--iterations N] [--merge-transition] [--overwrite]\n";
    return 2;
  }
  try {
    std::filesystem::path out;
    bool overwrite = false;
    Json options = Json::object();
    for (int i = 2; i < argc; ++i) {
      const std::string flag = argv[i];
      if (flag == "--out" && i + 1 < argc)
        out = argv[++i];
      else if (flag == "--iterations" && i + 1 < argc)
        options["iterations"] = std::stoi(argv[++i]);
      else if (flag == "--overwrite")
        overwrite = true;
      else if (flag == "--merge-transition")
        options["merge_transition_cells"] = true;
      else
        throw std::runtime_error("unknown optimizer argument");
    }
    if (out.empty())
      throw std::runtime_error("optimizer requires --out DIRECTORY");
    validateOptimization(options);
    auto mesh = readExchangeMesh(argv[1]);
    optimizeMixedMesh(mesh, options);
    bundle(mesh, out, overwrite);
    std::cout << "Optimized mesh saved in " << out << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "mesh optimization failed: " << e.what() << '\n';
    return 1;
  }
}
int runExportCad(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: babelmesh export-cad PROJECT.json --out FILE.step "
                 "[--entity ID] [--overwrite]\n";
    return 2;
  }
  try {
    std::filesystem::path out;
    std::string entity;
    bool overwrite = false;
    for (int i = 2; i < argc; ++i) {
      const std::string flag = argv[i];
      if (flag == "--out" && i + 1 < argc)
        out = argv[++i];
      else if (flag == "--entity" && i + 1 < argc)
        entity = argv[++i];
      else if (flag == "--overwrite")
        overwrite = true;
      else
        throw std::runtime_error("unknown CAD export argument");
    }
    if (out.empty())
      throw std::runtime_error("CAD export requires --out FILE.step");
    if (std::filesystem::exists(out) && !overwrite)
      throw std::runtime_error("STEP output exists; use --overwrite");
    const auto model = buildCadModel(argv[1]);
    const auto shape = entity.empty() ? model.fluidDomain()
                                      : model.entities().at(entity).shape;
    const auto temp = out.string() + ".staging-" + stamp();
    std::filesystem::create_directories(
        std::filesystem::absolute(out).parent_path());
    try {
      if (!STEPControl_Controller::Init())
        throw std::runtime_error("cannot initialize STEP transfer controller");
      const auto previousCascade =
          std::string(Interface_Static::CVal("xstep.cascade.unit"));
      const auto previousWrite =
          std::string(Interface_Static::CVal("write.step.unit"));
      struct RestoreStepUnits {
        std::string cascade, write;
        ~RestoreStepUnits() {
          Interface_Static::SetCVal("xstep.cascade.unit", cascade.c_str());
          Interface_Static::SetCVal("write.step.unit", write.c_str());
        }
      } restore{previousCascade, previousWrite};
      // B-Rep coordinates are already metres in the BabelMesh CAD contract.
      // Both the source and destination units must be set before the writer
      // constructs its transfer model; otherwise OCCT assumes millimetres.
      if (!Interface_Static::SetCVal("xstep.cascade.unit", "M") ||
          !Interface_Static::SetCVal("write.step.unit", "M"))
        throw std::runtime_error("cannot configure metre STEP transfer units");
      STEPControl_Writer writer;
      if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone ||
          writer.Write(temp.c_str()) != IFSelect_RetDone)
        throw std::runtime_error("OCCT STEP export failed");
      std::filesystem::rename(temp, out);
    } catch (...) {
      std::filesystem::remove(temp);
      throw;
    }
    std::cout << "STEP saved (meters): " << out << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "CAD export failed: " << e.what() << '\n';
    return 1;
  }
}
} // namespace babelmesh
