#include "babelmesh/cartesian_mesh.hpp"
#include "babelmesh/cad.hpp"
#include "babelsim/mesh_io.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <set>

namespace {

using Json = nlohmann::json;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Json baseProject(std::array<double, 3> maximum, std::array<int, 3> counts) {
    return {
        {"schema_version", 1}, {"units", {{"length", "m"}, {"output_length", "m"}}},
        {"cad", {{"imports", Json::array()}, {"operations", Json::array()}, {"groups", Json::array()}}},
        {"fluid_domain", {{"mode", "external"}, {"outer_bounds", {{"minimum", {0, 0, 0}}, {"maximum", maximum}}},
            {"solids", {"unused_in_background_test"}}, {"seeds", {{0.1, 0.1, 0.1}}}}},
        {"background", {{"bounds", {{"minimum", {0, 0, 0}}, {"maximum", maximum}}}, {"cell_counts", counts}, {"uniform_cells", true}}},
        {"sizing", {{"minimum_size", 0.125}, {"maximum_size", 1.0}, {"max_level", 3}, {"balance_face_neighbors", 1}, {"sources", Json::array()}}},
        {"boundary_layers", {{"enabled", false}, {"targets", Json::array()}}},
        {"quality", {{"profile", "strict-cfd-v1"}, {"mode", "strict"}}},
        {"resources", {{"max_cells", 10000}, {"max_memory_mb", 512}, {"threads", 1}}},
        {"output", {{"directory", "unused"}, {"formats", {"babelsim-v3", "vtu"}}, {"overwrite", true}}}
    };
}

Json cadSizingProject() {
    auto project = baseProject({1.0, 1.0, 1.0}, {2, 2, 2});
    project["cad"]["operations"] = Json::array({
        {{"id", "sizing_sphere"}, {"operation", "sphere"}, {"inputs", Json::array()},
            {"parameters", {{"center", {0.25, 0.5, 0.5}}, {"radius", 0.1}}}},
        {{"id", "plate_a"}, {"operation", "box"}, {"inputs", Json::array()},
            {"parameters", {{"minimum", {0.45, 0.45, 0.4}}, {"maximum", {0.55, 0.55, 0.45}}}}},
        {{"id", "plate_b"}, {"operation", "box"}, {"inputs", Json::array()},
            {"parameters", {{"minimum", {0.45, 0.45, 0.55}}, {"maximum", {0.55, 0.55, 0.6}}}}}
    });
    const auto group = [](const std::string& name, const std::string& operation) {
        return Json{{"name", name}, {"kind", "wall"},
            {"selector", {{"operation_id", operation}, {"entity_type", "face"}}}};
    };
    project["cad"]["groups"] = Json::array({group("sphere_surface", "sizing_sphere"),
        group("plate_a_surface", "plate_a"), group("plate_b_surface", "plate_b"),
        Json{{"name", "plate_a_edges"}, {"kind", "generic"},
            {"selector", {{"operation_id", "plate_a"}, {"entity_type", "edge"}}}}});
    project["fluid_domain"]["solids"] = {"sizing_sphere", "plate_a", "plate_b"};
    project["fluid_domain"]["seeds"] = {{0.9, 0.9, 0.9}};
    return project;
}

std::array<double, 6> cellBounds(const babelmesh::HexMesh& mesh, const std::array<std::int32_t, 8>& cell) {
    std::array<double, 6> result{1.0e100, 1.0e100, 1.0e100, -1.0e100, -1.0e100, -1.0e100};
    for (const auto id : cell) for (std::size_t axis = 0; axis < 3; ++axis) {
        const double value = mesh.points_m.at(static_cast<std::size_t>(id))[axis];
        result[axis] = std::min(result[axis], value);
        result[axis + 3U] = std::max(result[axis + 3U], value);
    }
    return result;
}

void checkFaceEdgeCornerBalance(const babelmesh::HexMesh& mesh) {
    for (std::size_t a = 0; a < mesh.cells.size(); ++a) for (std::size_t b = a + 1U; b < mesh.cells.size(); ++b) {
        const auto x = cellBounds(mesh, mesh.cells[a]);
        const auto y = cellBounds(mesh, mesh.cells[b]);
        bool touch = true;
        bool hasBoundary = false;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double overlap = std::min(x[axis + 3U], y[axis + 3U]) - std::max(x[axis], y[axis]);
            const double tol = 1.0e-12;
            if (overlap < -tol) touch = false;
            if (std::abs(overlap) <= tol) hasBoundary = true;
        }
        if (touch && hasBoundary) {
            const int difference = std::abs(mesh.cell_levels[a] - mesh.cell_levels[b]);
            require(difference <= 1, "octree is not 2:1 balanced across a face/edge/corner");
        }
    }
}

void checkOneCoarseFourFineInterface(const babelmesh::HexMesh& mesh) {
    std::set<std::int32_t> fineCells;
    std::set<std::int32_t> coarseCells;
    std::set<std::vector<std::int32_t>> faceKeys;
    double interfaceArea = 0.0;
    for (const auto& face : mesh.faces) {
        auto key = face.vertices;
        std::sort(key.begin(), key.end());
        require(faceKeys.insert(key).second, "adaptive topology contains a duplicate face");
        if (face.neighbour >= 0 && std::abs(mesh.cell_levels[static_cast<std::size_t>(face.owner)] -
                mesh.cell_levels[static_cast<std::size_t>(face.neighbour)]) == 1) {
            const auto ownerLevel = mesh.cell_levels[static_cast<std::size_t>(face.owner)];
            const auto neighbourLevel = mesh.cell_levels[static_cast<std::size_t>(face.neighbour)];
            const auto fine = ownerLevel > neighbourLevel ? face.owner : face.neighbour;
            const auto coarse = ownerLevel > neighbourLevel ? face.neighbour : face.owner;
            bool xOne = true;
            for (const auto vertex : face.vertices) {
                xOne = xOne && std::abs(mesh.points_m[static_cast<std::size_t>(vertex)][0] - 1.0) <= 1.0e-12;
            }
            if (xOne) {
                fineCells.insert(fine);
                coarseCells.insert(coarse);
                double twiceArea = 0.0;
                for (std::size_t i = 0; i < face.vertices.size(); ++i) {
                    const auto& a = mesh.points_m[static_cast<std::size_t>(face.vertices[i])];
                    const auto& b = mesh.points_m[static_cast<std::size_t>(face.vertices[(i + 1U) % face.vertices.size()])];
                    twiceArea += a[1] * b[2] - b[1] * a[2];
                }
                interfaceArea += 0.5 * std::abs(twiceArea);
            }
        }
    }
    require(fineCells.size() == 4U && coarseCells.size() == 1U,
        "one coarse cell must share an interface with exactly four fine cells");
    require(std::abs(interfaceArea - 1.0) <= 1.0e-12,
        "coarse-to-fine shared subfaces must tile the interface exactly once");
    std::vector<bool> usedPoints(mesh.points_m.size(), false);
    for (const auto& face : mesh.faces) for (const auto vertex : face.vertices) {
        usedPoints[static_cast<std::size_t>(vertex)] = true;
    }
    require(std::all_of(usedPoints.begin(), usedPoints.end(), [](bool used) { return used; }),
        "adaptive mesh contains an unreferenced/orphan point");
}

void readerRoundTrip(const babelmesh::HexMesh& mesh, const std::filesystem::path& file, double expectedVolume) {
    babelmesh::writeBabelSimV3(mesh, file);
    const auto readBack = babelsim::readMeshFile(file);
    readBack.validate();
    require(readBack.cellCount() == static_cast<babelsim::Index>(mesh.cells.size()), "BabelSim reader cell count mismatch");
    require(readBack.faceCount() == static_cast<babelsim::Index>(mesh.faces.size()), "BabelSim reader face count mismatch");
    double volume = 0.0;
    for (babelsim::Index cell = 0; cell < readBack.cellCount(); ++cell) volume += readBack.cellVolume(cell);
    require(std::abs(volume - expectedVolume) <= expectedVolume * 1.0e-12, "adaptive mesh total volume mismatch");
}

}  // namespace

int main() {
    try {
        const auto temp = std::filesystem::temp_directory_path();
        auto oneFineFace = baseProject({2.0, 1.0, 1.0}, {2, 1, 1});
        oneFineFace["sizing"]["sources"].push_back({{"id", "left_box"}, {"kind", "box"}, {"target_size", 0.5},
            {"region", {{"bounds", {{"minimum", {0.0, 0.0, 0.0}}, {"maximum", {0.9, 1.0, 1.0}}}}}}});
        const auto refined = babelmesh::makeAdaptiveCartesianMesh(oneFineFace);
        require(refined.cells.size() == 9U, "one refined root cell should produce eight fine cells plus one coarse cell");
        require(refined.faces.size() > 6U * refined.cells.size(), "adaptive interface should be represented by shared face subelements");
        require(refined.source_intersection_counts.at("left_box") > 0U, "box source should report affected cells");
        require(refined.sizing_source_markers.at("left_box").size() == refined.cells.size(),
            "every adaptive cell should carry an explicit source-membership marker");
        require(std::count(refined.sizing_source_markers.at("left_box").begin(),
                    refined.sizing_source_markers.at("left_box").end(), static_cast<std::uint8_t>(1U)) ==
                static_cast<std::ptrdiff_t>(refined.source_intersection_counts.at("left_box")),
            "source markers and reported source intersection count must agree");
        require(std::count(refined.strict_hex.begin(), refined.strict_hex.end(), false) > 0,
            "coarse-to-fine transition cell must be reported as polyhedral");
        require(std::count_if(refined.faces.begin(), refined.faces.end(), [](const auto& face) { return face.vertices.size() == 3U; }) > 0,
            "coarse-face split vertices must propagate into adjacent face edge rings and be triangulated through a shared face center");
        checkFaceEdgeCornerBalance(refined);
        checkOneCoarseFourFineInterface(refined);
        readerRoundTrip(refined, temp / "babelmesh-m3-one-coarse-four-fine.bsm", 2.0);

        const auto vtu = temp / "babelmesh-m3-transition.vtu";
        babelmesh::writeLegacyVtu(refined, vtu);
        std::ifstream vtuInput(vtu);
        const std::string vtuText((std::istreambuf_iterator<char>(vtuInput)), std::istreambuf_iterator<char>());
        require(vtuText.find("42 ") != std::string::npos && vtuText.find("Name=\"faceoffsets\"") != std::string::npos,
            "VTU must use explicit polyhedron face streams for transition cells");
        require(vtuText.find("Name=\"level\"") != std::string::npos &&
                vtuText.find("Name=\"cell_type\"") != std::string::npos &&
                vtuText.find("Name=\"source_left_box\"") != std::string::npos,
            "VTU cell data must expose refinement level, cell type, and per-source membership");

        auto levels = baseProject({2.0, 2.0, 2.0}, {2, 2, 2});
        levels["sizing"]["sources"].push_back({{"id", "corner_thin_band"}, {"kind", "box"}, {"target_size", 0.25},
            {"region", {{"bounds", {{"minimum", {0.02, 0.02, 0.02}}, {"maximum", {0.08, 0.08, 0.08}}}}}}});
        const auto multilevel = babelmesh::makeAdaptiveCartesianMesh(levels);
        require(*std::max_element(multilevel.cell_levels.begin(), multilevel.cell_levels.end()) >= 2,
            "thin box band away from the initial cell center must trigger at least three levels");
        checkFaceEdgeCornerBalance(multilevel);
        readerRoundTrip(multilevel, temp / "babelmesh-m3-multilevel.bsm", 8.0);

        auto rootVertex = baseProject({2.0, 2.0, 2.0}, {2, 2, 2});
        rootVertex["sizing"]["sources"].push_back({{"id", "multi_root_vertex"}, {"kind", "box"}, {"target_size", 0.125},
            {"region", {{"bounds", {{"minimum", {0.98, 0.98, 0.98}}, {"maximum", {1.02, 1.02, 1.02}}}}}}});
        const auto crossRoot = babelmesh::makeAdaptiveCartesianMesh(rootVertex);
        require(*std::max_element(crossRoot.cell_levels.begin(), crossRoot.cell_levels.end()) >= 3,
            "a thin feature at an internal root-grid vertex must refine across multiple roots for at least three levels");
        checkFaceEdgeCornerBalance(crossRoot);
        readerRoundTrip(crossRoot, temp / "babelmesh-m3-root-vertex.bsm", 8.0);

        auto sphereProject = baseProject({1.0, 1.0, 1.0}, {1, 1, 1});
        sphereProject["sizing"]["sources"].push_back({{"id", "off_center_sphere"}, {"kind", "sphere"}, {"target_size", 0.25},
            {"region", {{"center", {0.13, 0.17, 0.21}}, {"radius", 0.035}}}});
        const auto sphereMesh = babelmesh::makeAdaptiveCartesianMesh(sphereProject);
        require(sphereMesh.cells.size() > 1U, "off-center sphere must trigger refinement without a cell-center hit");
        checkFaceEdgeCornerBalance(sphereMesh);
        readerRoundTrip(sphereMesh, temp / "babelmesh-m3-sphere.bsm", 1.0);

        auto cylinderProject = baseProject({1.0, 1.0, 1.0}, {1, 1, 1});
        cylinderProject["sizing"]["sources"].push_back({{"id", "vertical_cylinder"}, {"kind", "cylinder"}, {"target_size", 0.25},
            {"region", {{"center", {0.5, 0.5, 0.5}}, {"axis", {0, 0, 1}}, {"radius", 0.06}, {"length", 0.4}}}});
        const auto cylinderMesh = babelmesh::makeAdaptiveCartesianMesh(cylinderProject);
        require(cylinderMesh.cells.size() > 1U, "finite cylinder must trigger local refinement");
        checkFaceEdgeCornerBalance(cylinderMesh);
        readerRoundTrip(cylinderMesh, temp / "babelmesh-m3-cylinder.bsm", 1.0);

        auto cadProject = cadSizingProject();
        const auto cadProjectFile = temp / "babelmesh-m3-cad-sizing-project.json";
        {
            std::ofstream output(cadProjectFile);
            output << cadProject.dump(2) << '\n';
        }
        const auto cadModel = babelmesh::buildCadModel(cadProjectFile);
        cadProject["sizing"]["sources"].push_back({{"id", "cad_surface"}, {"kind", "cad_group"},
            {"target_size", 0.25}, {"region", {{"cad_group", "sphere_surface"}}}});
        const auto cadSurfaceMesh = babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel);
        require(cadSurfaceMesh.source_intersection_counts.at("cad_surface") > 0U,
            "CAD group sizing must find cells intersecting an exact OCCT face group");
        checkFaceEdgeCornerBalance(cadSurfaceMesh);
        readerRoundTrip(cadSurfaceMesh, temp / "babelmesh-m3-cad-surface.bsm", 1.0);

        cadProject["sizing"]["sources"] = Json::array({{{"id", "cad_distance"}, {"kind", "distance"},
            {"target_size", 0.25}, {"region", {{"cad_group", "sphere_surface"}, {"distance", 0.2}}}}});
        const auto cadDistanceMesh = babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel);
        require(cadDistanceMesh.source_intersection_counts.at("cad_distance") > cadSurfaceMesh.source_intersection_counts.at("cad_surface"),
            "positive distance sizing should include a wider band than exact CAD surface sizing");
        checkFaceEdgeCornerBalance(cadDistanceMesh);
        readerRoundTrip(cadDistanceMesh, temp / "babelmesh-m3-cad-distance.bsm", 1.0);

        cadProject["sizing"]["sources"] = Json::array({{{"id", "cad_edge_distance"}, {"kind", "distance"},
            {"target_size", 0.25}, {"region", {{"cad_group", "plate_a_edges"}, {"distance", 0.1}}}}});
        const auto cadEdgeMesh = babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel);
        require(cadEdgeMesh.source_intersection_counts.at("cad_edge_distance") > 0U,
            "CAD edge groups should support exact distance-based local sizing");

        cadProject["sizing"]["sources"] = Json::array({{{"id", "planar_curvature"}, {"kind", "curvature"},
            {"target_size", 0.25}, {"region", {{"cad_group", "plate_a_surface"}, {"max_radius", 0.5}}}}});
        bool rejectedPlanarCurvature = false;
        try { (void)babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel); }
        catch (const std::exception& error) {
            rejectedPlanarCurvature = std::string(error.what()).find("no selected CAD face with radius of curvature") != std::string::npos;
        }
        require(rejectedPlanarCurvature, "curvature sizing must not silently accept an entirely planar CAD group");

        cadProject["sizing"]["sources"] = Json::array({{{"id", "curved_surface"}, {"kind", "curvature"},
            {"target_size", 0.25}, {"region", {{"cad_group", "sphere_surface"}, {"distance", 0.05}, {"max_radius", 1.0}}}}});
        const auto curvatureMesh = babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel);
        require(curvatureMesh.source_intersection_counts.at("curved_surface") > 0U,
            "curvature sizing should select high-curvature OCCT faces and refine their neighborhood");
        require(curvatureMesh.sizing_source_details.at("curved_surface").at("curvature_sampling") ==
                "3x3 UV samples restricted to trimmed CAD faces",
            "mesh report metadata must state the curvature sampling method and criterion");
        require(curvatureMesh.sizing_source_details.at("curved_surface").at("curvature_max_radius_m") == 1.0,
            "curvature selector radius must be reported separately from its target cell size");
        checkFaceEdgeCornerBalance(curvatureMesh);

        cadProject["sizing"]["sources"] = Json::array({{{"id", "plate_gap"}, {"kind", "gap"},
            {"target_size", 0.25}, {"region", {{"cad_group", "plate_a_surface"},
                {"other_cad_group", "plate_b_surface"}, {"distance", 0.12}}}}});
        const auto gapMesh = babelmesh::makeAdaptiveCartesianMesh(cadProject, &cadModel);
        require(gapMesh.source_intersection_counts.at("plate_gap") > 0U,
            "gap sizing should refine cells lying within the requested distance band of both CAD groups");
        checkFaceEdgeCornerBalance(gapMesh);

        auto overBudget = oneFineFace;
        overBudget["resources"]["max_cells"] = 8;
        bool rejected = false;
        try { (void)babelmesh::makeAdaptiveCartesianMesh(overBudget); }
        catch (const std::exception& error) { rejected = std::string(error.what()).find("max_cells") != std::string::npos; }
        require(rejected, "refinement over max_cells must fail with a budget diagnostic");

        for (const auto& file : {temp / "babelmesh-m3-one-coarse-four-fine.bsm", temp / "babelmesh-m3-transition.vtu",
                 temp / "babelmesh-m3-multilevel.bsm", temp / "babelmesh-m3-sphere.bsm", temp / "babelmesh-m3-cylinder.bsm",
                 temp / "babelmesh-m3-cad-surface.bsm", temp / "babelmesh-m3-cad-distance.bsm",
                 temp / "babelmesh-m3-cad-sizing-project.json"}) std::filesystem::remove(file);
        std::cout << "M3 adaptive acceptance passed: analytic and CAD sizing, curvature/gap bands, thin-band intersection, 2:1 face-edge-corner balance, shared transition faces, budget rejection, and BabelSim readback\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "M3 adaptive acceptance failed: " << error.what() << '\n';
        return 1;
    }
}
