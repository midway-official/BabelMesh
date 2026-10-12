#include "babelmesh/cad.hpp"
#include "babelmesh/fluid_domain.hpp"
#include "babelmesh/polyhedral_audit.hpp"
#include "babelsim/mesh_io.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

using Json = nlohmann::json;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Json externalBoxProject() {
    return {
        {"schema_version", 1}, {"units", {{"length", "m"}, {"output_length", "m"}}},
        {"cad", {
            {"imports", Json::array()},
            {"operations", Json::array({{{"id", "solid"}, {"operation", "box"}, {"inputs", Json::array()},
                {"parameters", {{"minimum", {1, 1, 1}}, {"maximum", {3, 3, 3}}}}}})},
            {"groups", Json::array()}
        }},
        {"fluid_domain", {{"mode", "external"}, {"outer_bounds", {{"minimum", {0, 0, 0}}, {"maximum", {4, 4, 4}}}},
            {"solids", {"solid"}}, {"seeds", {{0.5, 0.5, 0.5}}}}},
        {"background", {{"bounds", {{"minimum", {0, 0, 0}}, {"maximum", {4, 4, 4}}}},
            {"cell_counts", {4, 4, 4}}, {"uniform_cells", true}}},
        {"sizing", {{"minimum_size", 0.1}, {"maximum_size", 1.0}, {"surface_size", 0.5}, {"max_level", 4},
            {"balance_face_neighbors", 1}, {"sources", Json::array()}}},
        {"boundary_layers", {{"enabled", false}, {"targets", Json::array()}, {"strict_coverage", true}}},
        {"quality", {{"profile", "strict-cfd-v1"}, {"mode", "strict"}}},
        {"resources", {{"max_cells", 100000}, {"max_memory_mb", 256}, {"threads", 1}}},
        {"output", {{"directory", "unused"}, {"formats", {"babelsim-v3", "vtu"}}, {"overwrite", true}}}
    };
}

Json disconnectedInternalProject() {
    Json project = externalBoxProject();
    project["cad"]["operations"] = Json::array({
        {{"id", "fluid_a"}, {"operation", "box"}, {"inputs", Json::array()},
            {"parameters", {{"minimum", {0, 0, 0}}, {"maximum", {1, 1, 1}}}}},
        {{"id", "fluid_b"}, {"operation", "box"}, {"inputs", Json::array()},
            {"parameters", {{"minimum", {3, 3, 3}}, {"maximum", {4, 4, 4}}}}}
    });
    project["fluid_domain"] = {{"mode", "internal"}, {"fluid_solids", {"fluid_a", "fluid_b"}},
        {"seeds", {{0.5, 0.5, 0.5}}}, {"caps", Json::array()}};
    return project;
}

Json sphereExternalProject() {
    Json project = externalBoxProject();
    project["cad"]["operations"] = Json::array({
        {{"id", "sphere"}, {"operation", "sphere"}, {"inputs", Json::array()},
            {"parameters", {{"center", {1, 1, 1}}, {"radius", 0.7}}}}
    });
    project["fluid_domain"] = {{"mode", "external"},
        {"outer_bounds", {{"minimum", {0, 0, 0}}, {"maximum", {2, 2, 2}}}},
        {"solids", {"sphere"}}, {"seeds", {{0.1, 0.1, 0.1}}}};
    project["background"] = {{"bounds", {{"minimum", {0, 0, 0}}, {"maximum", {2, 2, 2}}}},
        {"cell_counts", {8, 8, 8}}, {"uniform_cells", true}};
    project["sizing"]["minimum_size"] = 0.015625;
    project["sizing"]["maximum_size"] = 0.25;
    project["sizing"]["surface_size"] = 0.015625;
    project["quality"]["mode"] = "strict";
    project["resources"]["max_cells"] = 500000;
    project["resources"]["max_memory_mb"] = 2048;
    return project;
}

Json rotatedBoxExternalProject() {
    Json project = sphereExternalProject();
    constexpr double angle = 0.52359877559829887308;
    project["cad"]["operations"] = Json::array({
        {{"id", "unrotated"}, {"operation", "box"}, {"inputs", Json::array()},
            {"parameters", {{"minimum", {0.5, 0.5, 0.5}}, {"maximum", {1.5, 1.5, 1.5}}}}},
        {{"id", "rotated"}, {"operation", "transform"}, {"inputs", {"unrotated"}},
            {"parameters", {{"rotation", {0.0, 0.0, angle}},
                {"translation", {0.6339745962155614, -0.3660254037844386, 0.0}}}}}
    });
    project["fluid_domain"]["solids"] = {"rotated"};
    return project;
}

}  // namespace

int main() {
    try {
        const auto projectPath = std::filesystem::temp_directory_path() / "babelmesh-m4-domain-project.json";
        {
            std::ofstream output(projectPath);
            output << externalBoxProject().dump(2) << '\n';
        }
        const auto cad = babelmesh::buildCadModel(projectPath);
        const auto inside = babelmesh::classifyFluidCell(cad, {0.1, 0.1, 0.1, 0.8, 0.8, 0.8});
        const auto outside = babelmesh::classifyFluidCell(cad, {1.2, 1.2, 1.2, 1.8, 1.8, 1.8});
        const auto intersect = babelmesh::classifyFluidCell(cad, {0.8, 1.2, 1.2, 1.2, 1.8, 1.8});
        require(inside.relation == babelmesh::CellRelation::inside,
            "cell fully inside the fluid volume must classify as inside");
        require(outside.relation == babelmesh::CellRelation::outside,
            "cell fully inside a subtracted solid must classify as outside the fluid domain");
        require(intersect.relation == babelmesh::CellRelation::intersect && intersect.candidate_surface_triangles > 0U &&
                intersect.intersecting_surface_triangles > 0U,
            "cell crossed by a CAD boundary triangle must classify as intersect with query evidence");
        const auto projected = babelmesh::projectPointToFluidBoundary(cad, {0.9, 2.0, 2.0});
        require(std::abs(projected.point_m[0] - 1.0) <= 1.0e-9 &&
                std::abs(projected.point_m[1] - 2.0) <= 1.0e-9 &&
                std::abs(projected.point_m[2] - 2.0) <= 1.0e-9 &&
                std::abs(projected.distance_m - 0.1) <= 1.0e-9,
            "CAD boundary projection must return the exact closest trimmed B-Rep face point");
        const auto background = babelmesh::makeUniformCartesianMesh(externalBoxProject());
        auto fluid = babelmesh::extractSeedConnectedFluid(background, externalBoxProject(), cad);
        require(fluid.cad_clipping_applied && fluid.cad_removed_cells == 8U,
            "seed-connected extraction must remove cells classified inside the obstacle solid; removed=" +
                std::to_string(fluid.cad_removed_cells) + ", retained=" + std::to_string(fluid.cells.size()));
        require(fluid.cells.size() == 56U && std::abs(fluid.cad_mesh_volume_relative_error) <= 1.0e-12,
            "grid-aligned external box extraction must preserve exact fluid volume");
        const auto cadPatch = std::find_if(fluid.patches.begin(), fluid.patches.end(),
            [](const auto& patch) { return patch.first == "cad_boundary"; });
        require(cadPatch != fluid.patches.end(), "exposed obstacle faces must receive a CAD boundary patch");
        const auto cadPatchIndex = static_cast<std::int32_t>(std::distance(fluid.patches.begin(), cadPatch));
        const auto cadBoundaryFaces = static_cast<std::size_t>(std::count_if(fluid.faces.begin(), fluid.faces.end(),
            [&](const auto& face) { return face.neighbour < 0 && face.patch == cadPatchIndex; }));
        require(cadBoundaryFaces == 24U, "grid-aligned obstacle wall should expose six faces with four tiles each");
        const auto refinedFluid = babelmesh::makeAdaptiveCartesianMesh(externalBoxProject(), &cad, true);
        require(refinedFluid.cad_clipping_applied && refinedFluid.cad_surface_refinement_applied &&
                refinedFluid.cad_surface_target_size_m == 0.5,
            "M4 generation path must refine CAD-intersecting cells to the configured surface size before extraction");
        require(refinedFluid.cells.size() > 56U &&
                std::abs(refinedFluid.cad_mesh_volume_relative_error) <= 1.0e-12,
            "surface refinement must preserve exact volume for a grid-aligned solid while adding boundary resolution; cells=" +
                std::to_string(refinedFluid.cells.size()) + ", relative error=" +
                std::to_string(refinedFluid.cad_mesh_volume_relative_error) + ", CAD volume=" +
                std::to_string(refinedFluid.cad_fluid_domain_volume_m3) + ", cell volume=" +
                std::to_string(std::accumulate(refinedFluid.cell_volumes_m3.begin(), refinedFluid.cell_volumes_m3.end(), 0.0)) +
                ", collapsed faces=" + std::to_string(refinedFluid.cad_collapsed_faces_removed) +
                ", deformed vertices=" + std::to_string(refinedFluid.cad_deformed_vertices) +
                ", max deformation=" + std::to_string(refinedFluid.cad_max_deformation_m));
        const auto meshFile = std::filesystem::temp_directory_path() / "babelmesh-m4-extracted-domain.bsm";
        babelmesh::writeBabelSimV3(fluid, meshFile);
        const auto readBack = babelsim::readMeshFile(meshFile);
        readBack.validate();
        require(readBack.cellCount() == 56 && readBack.faceCount() == static_cast<babelsim::Index>(fluid.faces.size()),
            "fluid extraction must preserve the BABELSIM_MESH 3 reader contract");
        std::filesystem::remove(meshFile);

        auto disconnectedProject = disconnectedInternalProject();
        const auto disconnectedPath = std::filesystem::temp_directory_path() / "babelmesh-m4-disconnected-project.json";
        {
            std::ofstream output(disconnectedPath);
            output << disconnectedProject.dump(2) << '\n';
        }
        const auto disconnectedCad = babelmesh::buildCadModel(disconnectedPath);
        const auto selectedComponent = babelmesh::extractSeedConnectedFluid(
            babelmesh::makeUniformCartesianMesh(disconnectedProject), disconnectedProject, disconnectedCad);
        require(selectedComponent.cells.size() == 1U && selectedComponent.cad_removed_cells == 63U,
            "seed-driven flood fill must discard an unseeded disconnected internal fluid component");
        require(selectedComponent.cad_selected_components == 1U &&
                std::abs(selectedComponent.cad_fluid_domain_volume_m3 - 1.0) <= 1.0e-12 &&
                std::abs(selectedComponent.cad_mesh_volume_relative_error) <= 1.0e-12,
            "volume comparison must use the selected seed-connected CAD component");
        std::filesystem::remove(disconnectedPath);

        auto sphereProject = sphereExternalProject();
        const auto spherePath = std::filesystem::temp_directory_path() / "babelmesh-m4-sphere-project.json";
        {
            std::ofstream output(spherePath);
            output << sphereProject.dump(2) << '\n';
        }
        const auto sphereCad = babelmesh::buildCadModel(spherePath);
        auto coarseStrictSphere = sphereProject;
        coarseStrictSphere["sizing"]["surface_size"] = 0.25;
        bool rejectedCoarseSurface = false;
        try { (void)babelmesh::makeAdaptiveCartesianMesh(coarseStrictSphere, &sphereCad, true); }
        catch (const std::exception& error) {
            std::cerr << "coarse strict sphere rejected: " << error.what() << '\n';
            rejectedCoarseSurface = std::string(error.what()).find("CAD surface snapping failed strict geometry checks") != std::string::npos;
        }
        require(rejectedCoarseSurface,
            "strict M4 must reject a coarse sphere staircase that exceeds the locked CAD surface deviation threshold");

        Json unspecifiedSurface = sphereProject;
        unspecifiedSurface["sizing"].erase("surface_size");
        bool rejectedUnspecifiedSurface = false;
        try { (void)babelmesh::makeAdaptiveCartesianMesh(unspecifiedSurface, &sphereCad, true); }
        catch (const std::exception& error) {
            rejectedUnspecifiedSurface = std::string(error.what()).find("requires an explicit sizing.surface_size") != std::string::npos;
        }
        require(rejectedUnspecifiedSurface,
            "CAD boundary generation must not silently inherit minimum_size as the surface resolution");

        Json unreachableSurface = sphereProject;
        unreachableSurface["sizing"]["minimum_size"] = 0.001;
        unreachableSurface["sizing"]["surface_size"] = 0.001;
        bool rejectedUnreachableSurface = false;
        try { (void)babelmesh::makeAdaptiveCartesianMesh(unreachableSurface, &sphereCad, true); }
        catch (const std::exception& error) {
            rejectedUnreachableSurface = std::string(error.what()).find("cannot be met before sizing.max_level") != std::string::npos;
        }
        require(rejectedUnreachableSurface,
            "CAD boundary generation must reject resolution targets that the configured max_level cannot attain");
        // Explicitly authorized 3 mm sphere geometry/quality tradeoff (2026-09-24).
        // Keep the default coarse-sphere rejection above and all face/cell
        // quality thresholds below unchanged; planar cases use their defaults.
        sphereProject["quality"]["max_surface_deviation_m"] = 0.003;
        const auto snappedSphere = babelmesh::makeAdaptiveCartesianMesh(sphereProject, &sphereCad, true);
        require(snappedSphere.cad_clipping_applied && snappedSphere.cad_surface_refinement_applied &&
                snappedSphere.cad_snapping_applied && snappedSphere.cad_snap_fraction > 0.0 &&
                snappedSphere.cad_surface_deviation_m <= snappedSphere.cad_surface_deviation_tolerance_m,
            "strict curved CAD case must fit inside its explicitly configured geometry tolerance");
        require(snappedSphere.cad_surface_deviation_overridden &&
                snappedSphere.cad_surface_deviation_tolerance_m == 0.003,
            "sphere geometry tradeoff must be explicitly reported as a 3 mm override");
        require(snappedSphere.cad_surface_target_size_m == 0.015625,
            "sphere boundary cells must honor the configured maximum isotropic edge length");
        require(std::all_of(snappedSphere.cell_volumes_m3.begin(), snappedSphere.cell_volumes_m3.end(),
                [](double volume) { return std::isfinite(volume) && volume > 0.0; }),
            "CAD snap rollback must retain positive finite polyhedral cell volumes");
        const auto sphereQuality=babelmesh::auditMeshQuality(snappedSphere);
        require(sphereQuality.minimum_regular_hex_scaled_jacobian>=0.1,
            "sphere regular-hex scaled Jacobian is below strict-cfd-v1 threshold: "+
                std::to_string(sphereQuality.minimum_regular_hex_scaled_jacobian));
        require(sphereQuality.maximum_non_orthogonality_degrees<=75.0 &&
                sphereQuality.p99_non_orthogonality_degrees<=60.0,
            "sphere non-orthogonality exceeds strict-cfd-v1 max/P99 gates");
        require(sphereQuality.maximum_internal_skewness<=0.5 &&
                sphereQuality.p99_internal_skewness<=0.25,
            "sphere internal skewness exceeds strict-cfd-v1 max/P99 gates");
        require(snappedSphere.cad_bidirectional_surface_audit &&
                snappedSphere.cad_mesh_to_surface_max_m<=0.001 &&
                snappedSphere.cad_surface_to_mesh_max_m<=0.001,
            "sphere must pass sampled mesh-to-CAD and CAD-to-mesh geometry checks within the configured 1 mm tolerance");
        const auto sphereMeshFile = std::filesystem::temp_directory_path() / "babelmesh-m4-sphere-candidate.bsm";
        babelmesh::writeBabelSimV3(snappedSphere, sphereMeshFile);
        const auto sphereReadBack = babelsim::readMeshFile(sphereMeshFile);
        sphereReadBack.validate();
        std::filesystem::remove(sphereMeshFile);
        std::filesystem::remove(spherePath);

        auto rotatedProject = rotatedBoxExternalProject();
        const auto rotatedPath = std::filesystem::temp_directory_path() / "babelmesh-m4-rotated-box-project.json";
        {
            std::ofstream output(rotatedPath);
            output << rotatedProject.dump(2) << '\n';
        }
        const auto rotatedCad = babelmesh::buildCadModel(rotatedPath);
        const auto rotatedMesh = babelmesh::makeAdaptiveCartesianMesh(rotatedProject, &rotatedCad, true);
        require(rotatedMesh.cad_snapping_applied &&
                rotatedMesh.cad_surface_deviation_m <= rotatedMesh.cad_surface_deviation_tolerance_m,
            "30-degree rotated box must pass strict CAD boundary snapping and sampled geometry tolerance");
        const auto rotatedMeshFile = std::filesystem::temp_directory_path() / "babelmesh-m4-rotated-box.bsm";
        babelmesh::writeBabelSimV3(rotatedMesh, rotatedMeshFile);
        babelsim::readMeshFile(rotatedMeshFile).validate();
        std::filesystem::remove(rotatedMeshFile);
        std::filesystem::remove(rotatedPath);
        require(babelmesh::toString(inside.relation) == "inside" &&
                babelmesh::toString(outside.relation) == "outside" &&
                babelmesh::toString(intersect.relation) == "intersect",
            "cell relation diagnostics must use stable names");
        std::filesystem::remove(projectPath);
        std::cout << "M4 fluid-domain foundation passed: B-Rep classification, seed-connected extraction, exact CAD projection and BabelSim readback; sphere snap fraction="
                  << snappedSphere.cad_snap_fraction << ", sampled deviation=" << snappedSphere.cad_surface_deviation_m
                  << " m, relative volume error=" << snappedSphere.cad_mesh_volume_relative_error
                  << ", quality max_nonorth=" << sphereQuality.maximum_non_orthogonality_degrees
                  << " p99_nonorth=" << sphereQuality.p99_non_orthogonality_degrees
                  << " max_skew=" << sphereQuality.maximum_internal_skewness
                  << " p99_skew=" << sphereQuality.p99_internal_skewness
                  << " min_hex_jacobian=" << sphereQuality.minimum_regular_hex_scaled_jacobian << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "M4 fluid-domain classification failed: " << error.what() << '\n';
        return 1;
    }
}
