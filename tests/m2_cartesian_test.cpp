#include "babelmesh/cartesian_mesh.hpp"

#include "babelsim/mesh_io.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

nlohmann::json project() {
    return {
        {"schema_version", 1},
        {"units", {{"length", "mm"}, {"output_length", "m"}}},
        {"cad", {{"imports", nlohmann::json::array()}, {"operations", nlohmann::json::array()}, {"groups", nlohmann::json::array()}}},
        {"fluid_domain", {{"mode", "external"}, {"outer_bounds", {{"minimum", {0, 0, 0}}, {"maximum", {2000, 3000, 4000}}}},
            {"solids", {"body"}}, {"seeds", {{100, 100, 100}}}}},
        {"background", {{"bounds", {{"minimum", {0, 0, 0}}, {"maximum", {2000, 3000, 4000}}}},
            {"cell_counts", {2, 3, 4}}, {"uniform_cells", true}}},
        {"sizing", {{"minimum_size", 0.01}, {"maximum_size", 1000.0}, {"max_level", 0}, {"sources", nlohmann::json::array()}}},
        {"boundary_layers", {{"enabled", false}, {"targets", nlohmann::json::array()}}},
        {"quality", {{"profile", "strict-cfd-v1"}, {"mode", "strict"}}},
        {"resources", {{"max_cells", 10000}, {"max_memory_mb", 256}, {"threads", 1}}},
        {"output", {{"directory", "unused"}, {"formats", {"babelsim-v3"}}, {"overwrite", true}}}
    };
}

std::array<double, 3> cellCentre(const babelmesh::HexMesh& mesh, std::int32_t cell) {
    std::array<double, 3> centre{};
    for (const auto point : mesh.cells.at(static_cast<std::size_t>(cell))) {
        const auto& xyz = mesh.points_m.at(static_cast<std::size_t>(point));
        for (std::size_t axis = 0; axis < 3; ++axis) centre[axis] += xyz[axis] / 8.0;
    }
    return centre;
}

void checkOwnerOutwardFaceRings(const babelmesh::HexMesh& mesh) {
    for (const auto& face : mesh.faces) {
        std::array<double, 3> area{}, centre{};
        for (std::size_t i = 0; i < face.vertices.size(); ++i) {
            const auto& a = mesh.points_m.at(static_cast<std::size_t>(face.vertices[i]));
            const auto& b = mesh.points_m.at(static_cast<std::size_t>(face.vertices[(i + 1U) % face.vertices.size()]));
            for (std::size_t axis = 0; axis < 3; ++axis) centre[axis] += a[axis] / static_cast<double>(face.vertices.size());
            area[0] += 0.5 * (a[1] * b[2] - a[2] * b[1]);
            area[1] += 0.5 * (a[2] * b[0] - a[0] * b[2]);
            area[2] += 0.5 * (a[0] * b[1] - a[1] * b[0]);
        }
        const auto ownerCentre = cellCentre(mesh, face.owner);
        double ownerDot = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis) ownerDot += area[axis] * (centre[axis] - ownerCentre[axis]);
        require(ownerDot > 0.0, "generated face area vector must point out of its owner cell");
        if (face.neighbour >= 0) {
            const auto neighbourCentre = cellCentre(mesh, face.neighbour);
            double neighbourDot = 0.0;
            for (std::size_t axis = 0; axis < 3; ++axis) neighbourDot += area[axis] * (neighbourCentre[axis] - ownerCentre[axis]);
            require(neighbourDot > 0.0 && face.patch == -1, "internal-face owner/neighbour orientation or patch is invalid");
        } else {
            require(face.patch >= 0, "boundary face must reference a physical patch");
        }
    }
}

void roundTrip(const babelmesh::HexMesh& mesh, const std::filesystem::path& file,
               std::int32_t cells, std::int32_t faces, std::int32_t vertices, double volume) {
    babelmesh::writeBabelSimV3(mesh, file);
    const babelsim::Mesh readBack = babelsim::readMeshFile(file);
    readBack.validate();
    require(readBack.cellCount() == cells, "BabelSim reader cell count mismatch");
    require(readBack.faceCount() == faces, "BabelSim reader face count mismatch");
    require(readBack.vertexCount() == vertices, "BabelSim reader vertex count mismatch");
    require(readBack.patchCount() == 6, "BabelSim reader patch count mismatch");
    require(readBack.patchName(0) == "xmin" && readBack.patchName(5) == "zmax", "BabelSim patch names did not round-trip");
    for (babelsim::Index cell = 0; cell < readBack.cellCount(); ++cell) {
        require(std::abs(readBack.cellVolume(cell) - volume) < 1.0e-12, "reader reported an incorrect cell volume");
    }
}

}  // namespace

int main() {
    try {
        for (const auto& invalid : std::vector<nlohmann::json>{0.0, -0.001, "0.001", nullptr}) {
            auto invalidProject = project();
            invalidProject["quality"]["max_surface_deviation_m"] = invalid;
            bool rejected = false;
            try { (void)babelmesh::makeUniformCartesianMesh(invalidProject); }
            catch (const std::exception& error) {
                rejected = std::string(error.what()).find("max_surface_deviation_m") != std::string::npos;
            }
            require(rejected, "invalid explicit geometry tolerance must be rejected before generation");
        }
        const auto generated = babelmesh::makeUniformCartesianMesh(project());
        require(generated.cells.size() == 24U, "expected 24 generated hex cells");
        require(generated.points_m.size() == 60U, "expected 60 unique grid vertices");
        require(generated.faces.size() == 98U, "expected 98 unique owner-neighbour faces");
        require(generated.patches.size() == 6U, "expected six exterior patches");

        checkOwnerOutwardFaceRings(generated);
        const auto temp = std::filesystem::temp_directory_path();
        const auto file = temp / "babelmesh-m2-roundtrip.bsm";
        roundTrip(generated, file, 24, 98, 60, 1.0);

        auto oneCellProject = project();
        oneCellProject["background"]["bounds"]["maximum"] = {1000, 1000, 1000};
        oneCellProject["background"]["cell_counts"] = {1, 1, 1};
        oneCellProject["fluid_domain"]["outer_bounds"]["maximum"] = {1000, 1000, 1000};
        const auto oneCell = babelmesh::makeUniformCartesianMesh(oneCellProject);
        checkOwnerOutwardFaceRings(oneCell);
        roundTrip(oneCell, temp / "babelmesh-m2-one-cell.bsm", 1, 6, 8, 1.0);

        auto twoByTwoProject = project();
        twoByTwoProject["background"]["cell_counts"] = {2, 2, 2};
        twoByTwoProject["sizing"]["maximum_size"] = 2000.0;
        const auto twoByTwo = babelmesh::makeUniformCartesianMesh(twoByTwoProject);
        checkOwnerOutwardFaceRings(twoByTwo);
        roundTrip(twoByTwo, temp / "babelmesh-m2-two-by-two.bsm", 8, 36, 27, 3.0);

        auto badMesh = generated;
        badMesh.faces.front().owner = 999;
        const auto badFile = temp / "babelmesh-m2-invalid-owner.bsm";
        babelmesh::writeBabelSimV3(badMesh, badFile);
        bool invalidRejected = false;
        try { (void)babelsim::readMeshFile(badFile); }
        catch (const std::exception&) { invalidRejected = true; }
        require(invalidRejected, "BabelSim reader must reject an invalid owner index");

        for (const auto& path : {file, temp / "babelmesh-m2-one-cell.bsm", temp / "babelmesh-m2-two-by-two.bsm", badFile}) {
            std::filesystem::remove(path);
        }
        std::cout << "M2 Cartesian acceptance passed: 1/8/24 hex cells, outward owner rings, metric normalization, invalid-face rejection, and BabelSim v3 reader round-trip\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "M2 Cartesian acceptance failed: " << error.what() << '\n';
        return 1;
    }
}
