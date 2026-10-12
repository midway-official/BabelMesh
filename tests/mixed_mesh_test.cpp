#include "babelmesh/hybrid_mesher.hpp"
#include "babelmesh/mesh_optimizer.hpp"
#include "babelmesh/mixed_mesh.hpp"
#include "babelsim/mesh_io.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
using namespace babelmesh;
void require(bool condition, const char *why) {
  if (!condition)
    throw std::runtime_error(why);
}
BoundaryMap boundaryFor(const MixedMesh &m) {
  BoundaryMap b;
  for (const auto &cell : m.cells)
    for (const auto &face : cellFaceRings(cell))
      b.emplace(faceKey(face), 0);
  return b;
}
template <class Function> void rejects(Function &&f, const char *why) {
  try {
    f();
  } catch (const std::exception &) {
    return;
  }
  throw std::runtime_error(why);
}
} // namespace
int main(int argc, char **argv) {
  try {
    using namespace babelmesh;
    MixedMesh m;
    m.patches = {{"wall", "wall"}};
    m.points_m = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0},    {0, 1, 0}, {0, 0, 1},
                  {1, 0, 1}, {1, 1, 1}, {0, 1, 1},    {2, 0, 0}, {3, 0, 0},
                  {3, 1, 0}, {2, 1, 0}, {2.5, .5, 1}, {4, 0, 0}, {5, 0, 0},
                  {4, 1, 0}, {4, 0, 1}, {6, 0, 0},    {7, 0, 0}, {6, 1, 0},
                  {6, 0, 1}, {7, 0, 1}, {6, 1, 1}};
    m.cells = {
        {CellType::hexahedron, CellRegion::cartesian, {0, 1, 2, 3, 4, 5, 6, 7}},
        {CellType::pyramid, CellRegion::transition, {8, 9, 10, 11, 12}},
        {CellType::tetrahedron, CellRegion::transition, {13, 14, 15, 16}},
        {CellType::prism,
         CellRegion::boundary_layer,
         {17, 18, 19, 20, 21, 22}}};
    buildMixedFaces(m, boundaryFor(m));
    require(std::abs(m.volumes_m3[0] - 1) < 1e-13, "hex volume");
    require(std::abs(m.volumes_m3[1] - 1.0 / 3) < 1e-13, "pyramid volume");
    require(std::abs(m.volumes_m3[2] - 1.0 / 6) < 1e-13, "tet volume");
    require(std::abs(m.volumes_m3[3] - .5) < 1e-13, "prism volume");
    const auto directory =
        std::filesystem::temp_directory_path() / "babelmesh-mixed-contract";
    std::filesystem::create_directories(directory);
    writeMixedBabelSimV3(m, directory / "mesh.bsm");
    writeMixedCheckpoint(m, directory / "mesh.json");
    const auto solver = babelsim::readMeshFile(directory / "mesh.bsm");
    solver.validate();
    require(solver.cellCount() == 4, "real BabelSim reader cell count");
    for (int c = 0; c < 4; ++c)
      require(std::abs(solver.cellVolume(c) -
                       m.volumes_m3[static_cast<std::size_t>(c)]) < 1e-13,
              "reader and producer volume differ");
    const auto roundtrip = readMixedCheckpoint(directory / "mesh.json");
    require(roundtrip.faces.size() == m.faces.size(), "checkpoint round-trip");
    auto bounded = m;
    bounded.metadata["domain_bounds_constraint"] = {{"minimum_m", {0, 0, 0}},
                                                    {"maximum_m", {7, 1, 1}}};
    validateMixedMesh(bounded);
    require(bounded.metadata.at("bounds_audit").at("maximum_excess_m") == 0,
            "closed domain envelope rejected a boundary vertex");
    bounded.points_m[0][0] = -1e-12;
    rejects([&] { validateMixedMesh(bounded); },
            "out-of-bounds node accepted with a geometric tolerance");
    rejects([&] { optimizeMixedMesh(bounded); },
            "optimizer accepted an expanded domain");
    auto duplicate = m;
    duplicate.faces.push_back(duplicate.faces[0]);
    rejects([&] { validateMixedMesh(duplicate); }, "duplicate face accepted");
    auto open = m;
    open.faces.pop_back();
    rejects([&] { validateMixedMesh(open); }, "open cell accepted");
    auto inverted = m;
    std::swap(inverted.cells[2].vertices[1], inverted.cells[2].vertices[2]);
    rejects([&] { buildMixedFaces(inverted, boundaryFor(inverted)); },
            "inverted tet accepted");
    rejects(
        [&] {
          auto missing = boundaryFor(m);
          missing.erase(missing.begin());
          buildMixedFaces(m, missing);
        },
        "unassigned boundary accepted");
    auto nonmanifold = m;
    nonmanifold.cells.push_back(m.cells[0]);
    nonmanifold.cells.push_back(m.cells[0]);
    rejects([&] { buildMixedFaces(nonmanifold, boundaryFor(nonmanifold)); },
            "overlapping cells accepted");
    auto invalidPlanar = m;
    invalidPlanar.generation_dimension = 2;
    invalidPlanar.extrusion_axis = 2;
    invalidPlanar.extrusion_thickness_m = 1;
    rejects([&] { validateMixedMesh(invalidPlanar); },
            "2D mode accepted tetrahedra or pyramids");
    MixedMesh poly;
    poly.patches = {{"wall", "wall"}};
    poly.points_m = {{0, 0, 0},     {1, 0, 0},    {1.3, .7, 0}, {.5, 1.2, 0},
                     {-.2, .6, 0},  {0, 0, .1},   {1, 0, .1},   {1.3, .7, .1},
                     {.5, 1.2, .1}, {-.2, .6, .1}};
    std::vector<FaceKey> rings = {{4, 3, 2, 1, 0}, {5, 6, 7, 8, 9}};
    for (std::int32_t i = 0; i < 5; ++i)
      rings.push_back({i, (i + 1) % 5, (i + 1) % 5 + 5, i + 5});
    poly.cells.emplace_back(CellType::polyhedron, CellRegion::cartesian,
                            FaceKey{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, rings);
    poly.generation_dimension = 2;
    poly.extrusion_axis = 2;
    poly.extrusion_thickness_m = .1;
    buildMixedFaces(poly, boundaryFor(poly));
    writeMixedBabelSimV3(poly, directory / "poly.bsm");
    const auto polySolver = babelsim::readMeshFile(directory / "poly.bsm");
    polySolver.validate();
    require(std::abs(polySolver.cellVolume(0) - poly.volumes_m3[0]) < 1e-13,
            "real polygonal prism volume differs in BabelSim");
    auto multilayer = poly;
    multilayer.points_m[5][2] = .05;
    rejects([&] { validateMixedMesh(multilayer); },
            "2D accepted more than two planes");
    auto mortar = poly;
    mortar.points_m = {{0, 0, 0},  {.5, 0, 0},  {1, 0, 0},  {1, .5, 0},
                       {1, 1, 0},  {.5, 1, 0},  {0, 1, 0},  {0, .5, 0},
                       {0, 0, .1}, {.5, 0, .1}, {1, 0, .1}, {1, .5, .1},
                       {1, 1, .1}, {.5, 1, .1}, {0, 1, .1}, {0, .5, .1}};
    std::vector<FaceKey> mortarFaces = {{7, 6, 5, 4, 3, 2, 1, 0},
                                        {8, 9, 10, 11, 12, 13, 14, 15}};
    for (std::int32_t i = 0; i < 8; ++i)
      mortarFaces.push_back({i, (i + 1) % 8, (i + 1) % 8 + 8, i + 8});
    mortar.cells = {
        MixedCell{CellType::polyhedron,
                  CellRegion::cartesian,
                  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                  mortarFaces}};
    buildMixedFaces(mortar, boundaryFor(mortar));
    writeMixedBabelSimV3(mortar, directory / "mortar.bsm");
    const auto mortarSolver = babelsim::readMeshFile(directory / "mortar.bsm");
    mortarSolver.validate();
    require(std::abs(mortarSolver.cellVolume(0) - .1) < 1e-13,
            "BabelSim rejected hanging-node cap polygons");
    // The shared triangular face has analytic normalized skewness > .25.
    // Removing that face through a convex union improves FVM geometry while
    // retaining the exact physical boundary and the eight Cartesian cells.
    MixedMesh cleanup;
    cleanup.patches = {{"wall", "wall"}};
    cleanup.points_m = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {.98, .01, .5}, {.98, .01, -.5}};
    cleanup.cells = {
        {CellType::tetrahedron, CellRegion::transition, {0, 1, 2, 3}},
        {CellType::tetrahedron, CellRegion::transition, {0, 2, 1, 4}}};
    for (std::int32_t i = 0; i < 8; ++i) {
      const auto start = static_cast<std::int32_t>(cleanup.points_m.size());
      for (std::size_t v = 0; v < 8; ++v) {
        auto p = m.points_m[v];
        p[0] += 2 + 2 * i;
        cleanup.points_m.push_back(p);
      }
      cleanup.cells.push_back({CellType::hexahedron,
                               CellRegion::cartesian,
                               {start, start + 1, start + 2, start + 3,
                                start + 4, start + 5, start + 6, start + 7}});
    }
    std::map<FaceKey, int> occurrences;
    for (const auto &c : cleanup.cells)
      for (const auto &f : cellFaceRings(c))
        ++occurrences[faceKey(f)];
    BoundaryMap outside;
    for (const auto &item : occurrences)
      if (item.second == 1)
        outside.emplace(item.first, 0);
    buildMixedFaces(cleanup, outside);
    const double initialVolume = mixedMeshReport(cleanup).at("volume_m3");
    require(!mixedMeshReport(cleanup)
                 .at("quality")
                 .at("face_quality_thresholds_passed")
                 .get<bool>(),
            "cleanup fixture does not exercise a failing skew percentile");
    optimizeMixedMesh(cleanup, {{"merge_transition_cells", true}});
    const auto cleanupReport = mixedMeshReport(cleanup);
    require(cleanup.cells.size() == 9 &&
                cleanup.metadata.at("optimization")
                        .at("merged_transition_tet_pairs") == 1,
            "convex transition tetrahedron merge was not performed");
    require(cleanupReport.at("quality")
                .at("face_quality_thresholds_passed")
                .get<bool>(),
            "convex merge did not improve the analytic skew fixture");
    require(std::abs(cleanupReport.at("volume_m3").get<double>() -
                     initialVolume) < 1e-12,
            "convex merge changed the fluid volume");
    writeMixedBabelSimV3(cleanup, directory / "cleanup.bsm");
    const auto cleanupSolver =
        babelsim::readMeshFile(directory / "cleanup.bsm");
    cleanupSolver.validate();
    require(cleanupSolver.cellCount() == 9,
            "BabelSim rejected merged topology");
    require(argc > 1, "test requires its source-root argument");
    for (const auto *name : {"three", "four"}) {
      const auto file = std::filesystem::path(argv[1]) / "tests" / "fixtures" /
                        (std::string("topology-") + name + "-tet-star.json");
      auto star = readMixedCheckpoint(file);
      const auto original = star;
      const auto before = mixedMeshReport(star);
      require(!before.at("quality").at("face_quality_thresholds_passed"),
              "tet-star fixture must exercise poor initial quality");
      auto disabled = star;
      optimizeMixedMesh(disabled, {{"topology_swaps", false}});
      require(disabled.cells.size() == star.cells.size() &&
                  disabled.points_m == star.points_m,
              "disabled swaps modified the closed boundary-only fixture");
      optimizeMixedMesh(star);
      const auto after = mixedMeshReport(star);
      const auto &topology =
          star.metadata.at("optimization").at("topology_swaps");
      const bool three = std::string(name) == "three";
      require(topology.at(three ? "accepted_3_to_2" : "accepted_4_to_4") == 1,
              "expected local tetrahedron swap was not accepted");
      require(
          star.cells.size() == (three ? 2U : 4U) &&
              star.points_m == original.points_m,
          "tet swap added cells/nodes or moved constrained boundary points");
      require(std::abs(after.at("volume_m3").get<double>() -
                       before.at("volume_m3").get<double>()) < 1e-12,
              "tet swap changed the cavity volume");
      BoundaryMap oldSkin, newSkin;
      for (const auto &f : original.faces)
        if (f.neighbour < 0)
          oldSkin.emplace(faceKey(f.vertices), f.patch);
      for (const auto &f : star.faces)
        if (f.neighbour < 0)
          newSkin.emplace(faceKey(f.vertices), f.patch);
      require(oldSkin == newSkin,
              "tet swap retriangulated the fixed cavity skin");
      for (const auto *metric :
           {"max_non_orthogonality_degrees", "p99_non_orthogonality_degrees",
            "max_internal_skewness", "p99_internal_skewness"})
        require(after.at("quality").at(metric).get<double>() <=
                    before.at("quality").at(metric).get<double>() + 1e-8,
                "tet swap worsened a global maximum/p99 FVM metric");
      if (!three)
        require(after.at("quality").at("face_quality_thresholds_passed"),
                "4-to-4 fixture did not cross its analytic quality gate");
      writeMixedBabelSimV3(star, directory / "swap.bsm");
      const auto read = babelsim::readMeshFile(directory / "swap.bsm");
      read.validate();
      require(static_cast<std::size_t>(read.cellCount()) == star.cells.size(),
              "BabelSim rejected swapped topology");
      auto layer = original;
      for (auto &cell : layer.cells)
        cell.region = CellRegion::boundary_layer;
      optimizeMixedMesh(layer);
      require(layer.cells.size() == original.cells.size() &&
                  layer.points_m == original.points_m,
              "topology optimizer modified fixed layer cells");
    }
    rejects(
        [&] {
          validateOptimization({{"topology_swaps", "yes"}});
        },
        "optimizer accepted an invalid topology-swaps option type");
    MixedMesh regular;
    regular.patches = {{"wall", "wall"}};
    regular.points_m = {{0, 0, 0},
                        {1, 0, 0},
                        {.5, std::sqrt(3.) / 2, 0},
                        {.5, std::sqrt(3.) / 6, std::sqrt(2. / 3)}};
    regular.cells = {
        {CellType::tetrahedron, CellRegion::transition, {0, 1, 2, 3}}};
    buildMixedFaces(regular, boundaryFor(regular));
    const double ratio = mixedMeshReport(regular)
                             .at("tetrahedron_shape_quality")
                             .at("regions")
                             .at("transition")
                             .at("minimum_mean_ratio");
    require(std::abs(ratio - 1) < 1e-12,
            "mean-ratio metric is not normalized for a regular tetrahedron");
    std::filesystem::remove_all(directory);
    std::cout << "Four real mixed cell types, shell rejection, checkpoint and "
                 "BabelSim reader passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
