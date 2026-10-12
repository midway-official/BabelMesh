#pragma once
#include "babelmesh/cad.hpp"
#include "babelmesh/mixed_mesh.hpp"
namespace babelmesh {
struct BackgroundGrid {
  std::vector<std::array<double, 3>> points;
  std::vector<MixedCell> volumes;
  std::vector<FaceKey> polygons;
  nlohmann::json report;
};
// Face-conforming octree in 3D; plane-only, edge-conforming quadtree in 2D.
BackgroundGrid generateBackgroundGrid(const nlohmann::json &project,
                                      const CadModel &cad);
} // namespace babelmesh
