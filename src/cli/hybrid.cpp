#include "babelmesh/hybrid_mesher.hpp"
#include "babelmesh/mesh_conversion.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace babelmesh {
int runGenerateHybrid(int argc, char **argv) {
  try {
    if (argc < 2) {
      std::cerr << "usage: babelmesh generate PROJECT.json [--out DIRECTORY] "
                   "[--overwrite]\n";
      return 2;
    }
    const auto projectFile = std::filesystem::absolute(argv[1]);
    nlohmann::json project;
    {
      std::ifstream in(projectFile);
      if (!in)
        throw std::runtime_error("cannot open project");
      in >> project;
    }
    std::filesystem::path out;
    for (int i = 2; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--out" && i + 1 < argc)
        out = argv[++i];
      else if (option == "--overwrite")
        project["output"]["overwrite"] = true;
      else {
        std::cerr << "unknown hybrid generate argument: " << option << '\n';
        return 2;
      }
    }
    project = normalizeHybridProject(project);
    validateHybridProject(project);
    if (out.empty()) {
      out = project.at("output").at("directory").get<std::string>();
      if (out.is_relative())
        out = projectFile.parent_path() / out;
    }
    out = std::filesystem::absolute(out);
    const bool overwrite = project.at("output").at("overwrite").get<bool>();
    const std::vector<std::string> artifacts = {
        "mesh.bsm",     "mesh.vtu", "mesh.json", "mesh_report.json",
        "project.json", "mesh.msh", "mesh.cgns"};
    for (const auto &name : artifacts)
      if (std::filesystem::exists(out / name) && !overwrite)
        throw std::runtime_error("output exists and overwrite is false: " +
                                 (out / name).string());
    const auto start = std::chrono::steady_clock::now();
    const auto cad = buildCadModel(projectFile);
    auto mesh = generateHybridMesh(project, cad);
    auto report = mixedMeshReport(mesh);
    report["project_file"] = projectFile.string();
    report["elapsed_seconds"] =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    report["status"] =
        report.at("quality").at("face_quality_thresholds_passed").get<bool>()
            ? "topology-and-face-quality-passed"
            : "topology-passed-quality-diagnostic";
    std::filesystem::create_directories(out);
    // All expensive generation and validation precedes publishing any result.
    const auto staging =
        out /
        (".staging-" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(staging);
    try {
      for (const auto &format : project.at("output").at("formats")) {
        if (format == "babelsim-v3")
          writeMixedBabelSimV3(mesh, staging / "mesh.bsm");
        else if (format == "vtu")
          writeMixedVtu(mesh, staging / "mesh.vtu");
        else if (format == "fluent-ascii")
          writeFluentAscii(mesh, staging / "mesh.msh");
        else if (format == "cgns" || format == "cgns-poly")
          writeCgns(mesh, staging / "mesh.cgns", format == "cgns-poly");
      }
      writeMixedCheckpoint(mesh, staging / "mesh.json");
      {
        std::ofstream stream(staging / "mesh_report.json");
        stream << report.dump(2) << '\n';
        if (!stream)
          throw std::runtime_error("cannot write report");
      }
      {
        std::ofstream stream(staging / "project.json");
        stream << project.dump(2) << '\n';
        if (!stream)
          throw std::runtime_error("cannot write project snapshot");
      }
      for (const auto &name : artifacts)
        if (std::filesystem::exists(staging / name))
          std::filesystem::rename(staging / name, out / name);
      std::filesystem::remove(staging);
    } catch (...) {
      std::filesystem::remove_all(staging);
      throw;
    }
    std::cout << "Generated " << mesh.cells.size() << " hybrid cells in " << out
              << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "hybrid generate failed: " << e.what() << '\n';
    return 1;
  } catch (const std::string &e) {
    std::cerr << "hybrid backend failed: " << e << '\n';
    return 1;
  } catch (...) {
    std::cerr << "hybrid backend failed with an unknown exception\n";
    return 1;
  }
}
int runCheckMixed(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: babelmesh check MESH.json\n";
    return 2;
  }
  try {
    const auto mesh = readMixedCheckpoint(argv[1]);
    std::cout << mixedMeshReport(mesh).dump(2) << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "mesh check failed: " << e.what() << '\n';
    return 1;
  }
}
} // namespace babelmesh
