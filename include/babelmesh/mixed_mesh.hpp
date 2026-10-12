#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace babelmesh {

// Values are the linear VTK cell identifiers, not flags on a nominal hex.
enum class CellType : std::uint8_t {
  tetrahedron = 10,
  hexahedron = 12,
  prism = 13,
  pyramid = 14,
  polyhedron = 42
};
enum class CellRegion : std::uint8_t {
  cartesian = 0,
  boundary_layer = 1,
  transition = 2
};

struct MixedCell {
  CellType type{CellType::tetrahedron};
  CellRegion region{CellRegion::transition};
  std::vector<std::int32_t> vertices;
  // Explicit rings for conforming Cartesian coarse/fine cells and polygonal
  // pyramids.
  std::vector<std::vector<std::int32_t>> poly_faces;
  MixedCell() = default;
  MixedCell(CellType t, CellRegion r, std::vector<std::int32_t> v,
            std::vector<std::vector<std::int32_t>> f = {})
      : type(t), region(r), vertices(std::move(v)), poly_faces(std::move(f)) {}
};
struct MixedFace {
  std::vector<std::int32_t> vertices;
  std::int32_t owner{-1}, neighbour{-1}, patch{-1};
};
struct MeshPatch {
  std::string name, kind;
};
using FaceKey = std::vector<std::int32_t>;
using BoundaryMap = std::map<FaceKey, std::int32_t>;

struct MixedMesh {
  int generation_dimension{3};
  int extrusion_axis{2};
  double extrusion_thickness_m{};
  std::vector<std::array<double, 3>> points_m;
  std::vector<MixedCell> cells;
  std::vector<MixedFace> faces;
  std::vector<MeshPatch> patches;
  std::vector<double> volumes_m3;
  std::vector<std::array<double, 3>> centres_m;
  nlohmann::json metadata = nlohmann::json::object();
};

FaceKey faceKey(FaceKey vertices);
std::vector<FaceKey> cellFaceRings(const MixedCell &cell);
const char *cellTypeName(CellType type);
// Resolve every interface once; require an exact declaration for every external
// face.
void buildMixedFaces(MixedMesh &mesh, const BoundaryMap &boundary);
void validateMixedMesh(MixedMesh &mesh);
nlohmann::json validateMixedMeshBounds(const MixedMesh &mesh);
nlohmann::json mixedMeshReport(const MixedMesh &mesh);
void writeMixedBabelSimV3(const MixedMesh &mesh,
                          const std::filesystem::path &file);
void writeMixedVtu(const MixedMesh &mesh, const std::filesystem::path &file);
void writeMixedCheckpoint(const MixedMesh &mesh,
                          const std::filesystem::path &file);
MixedMesh readMixedCheckpoint(const std::filesystem::path &file);

} // namespace babelmesh
