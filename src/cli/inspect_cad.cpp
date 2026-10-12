#include "babelmesh/cad.hpp"

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace babelmesh {
namespace {

std::string surfaceVtu(const CadModel& model) {
    std::ostringstream points, connectivity, offsets, types, entityIds, faceIds;
    points.imbue(std::locale::classic());
    points << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::size_t totalPoints = 0;
    std::size_t totalTriangles = 0;
    for (const auto& pair : model.entities()) {
        totalPoints += pair.second.surface_mesh.vertices.size();
        totalTriangles += pair.second.surface_mesh.triangles.size();
    }
    std::size_t pointBase = 0;
    std::size_t triangleOffset = 0;
    std::size_t entityIndex = 0;
    for (const auto& pair : model.entities()) {
        for (const auto& point : pair.second.surface_mesh.vertices) {
            points << point[0] << ' ' << point[1] << ' ' << point[2] << ' ';
        }
        for (const auto& triangle : pair.second.surface_mesh.triangles) {
            for (const auto vertex : triangle.vertices) connectivity << (pointBase + vertex) << ' ';
            triangleOffset += 3;
            offsets << triangleOffset << ' ';
            types << "5 ";
            entityIds << entityIndex << ' ';
            faceIds << triangle.face_index << ' ';
        }
        pointBase += pair.second.surface_mesh.vertices.size();
        ++entityIndex;
    }
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           << " <UnstructuredGrid>\n  <Piece NumberOfPoints=\"" << totalPoints << "\" NumberOfCells=\"" << totalTriangles << "\">\n"
           << "   <PointData/>\n   <CellData>\n"
           << "    <DataArray type=\"UInt64\" Name=\"cad_entity_id\" format=\"ascii\">" << entityIds.str() << "</DataArray>\n"
           << "    <DataArray type=\"UInt64\" Name=\"cad_face_id\" format=\"ascii\">" << faceIds.str() << "</DataArray>\n"
           << "   </CellData>\n   <Points>\n"
           << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << points.str() << "</DataArray>\n"
           << "   </Points>\n   <Cells>\n"
           << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << connectivity.str() << "</DataArray>\n"
           << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << offsets.str() << "</DataArray>\n"
           << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << types.str() << "</DataArray>\n"
           << "   </Cells>\n  </Piece>\n </UnstructuredGrid>\n</VTKFile>\n";
    return output.str();
}

void writeAtomically(const std::filesystem::path& path, const std::string& contents, const std::string& suffix) {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto temporaryPath = path.string() + ".tmp." + std::to_string(::getpid()) + "." + std::to_string(stamp) + suffix;
    {
        std::ofstream output(temporaryPath, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("cannot create temporary output: " + temporaryPath);
        output.imbue(std::locale::classic());
        output << contents;
        output.flush();
        if (!output) throw std::runtime_error("failed while writing temporary output: " + temporaryPath);
    }
    std::error_code error;
    std::filesystem::rename(temporaryPath, path, error);
    if (error) {
        std::filesystem::remove(temporaryPath);
        throw std::runtime_error("cannot publish output file: " + error.message());
    }
}

}  // namespace

int runInspectCad(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: babelmesh inspect-cad PROJECT.json [--out DIRECTORY]\n";
            return 2;
        }
        const std::filesystem::path projectPath = std::filesystem::absolute(argv[1]);
        std::optional<std::filesystem::path> requestedOutput;
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--out" && i + 1 < argc) requestedOutput = argv[++i];
            else {
                std::cerr << "unknown inspect-cad argument: " << option << "\n";
                return 2;
            }
        }
        std::ifstream projectStream(projectPath);
        if (!projectStream) throw std::runtime_error("cannot open project file: " + projectPath.string());
        nlohmann::json project;
        projectStream >> project;
        const bool overwrite = project.value("output", nlohmann::json::object()).value("overwrite", false);
        std::filesystem::path outputDirectory;
        if (requestedOutput) outputDirectory = *requestedOutput;
        else {
            const std::string configured = project.at("output").at("directory").get<std::string>();
            outputDirectory = configured;
            if (outputDirectory.is_relative()) outputDirectory = projectPath.parent_path() / outputDirectory;
        }
        std::filesystem::create_directories(outputDirectory);
        const auto reportPath = outputDirectory / "geometry_report.json";
        const auto surfacePath = outputDirectory / "cad_surface.vtu";
        if ((std::filesystem::exists(reportPath) || std::filesystem::exists(surfacePath) || std::filesystem::exists(outputDirectory / "cad_preview.json")) && !overwrite) {
            throw std::runtime_error("CAD output already exists and output.overwrite is false: " + outputDirectory.string());
        }

        const CadModel model = buildCadModel(projectPath);
        const auto report = geometryReport(model, projectPath);
        writeAtomically(surfacePath, surfaceVtu(model), ".vtu");
        writeAtomically(reportPath, report.dump(2) + "\n", ".json");
        writeAtomically(outputDirectory/"cad_preview.json",cadPreview(model).dump()+"\n",".json");
        std::cout << "CAD inspection passed: " << model.entities().size() << " entities, "
                  << model.groups().size() << " physical groups\n"
                  << "report: " << reportPath.string() << '\n'
                  << "surface: " << surfacePath.string() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "inspect-cad failed: " << e.what() << '\n';
        return 1;
    }
}

}  // namespace babelmesh
