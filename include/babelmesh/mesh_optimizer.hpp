#pragma once
#include "babelmesh/mixed_mesh.hpp"
namespace babelmesh {
// Preserve Cartesian and boundary-layer geometry. Relocation and optional
// transition-tet swaps/convex merges preserve positive shells, fixed domain
// bounds, and the maximum/p99 FVM quality gates without adding cells.
void validateOptimization(const nlohmann::json &options);
nlohmann::json improveTransitionTetTopology(MixedMesh &mesh,
                                            std::size_t maxTrials);
void optimizeMixedMesh(
    MixedMesh &mesh, const nlohmann::json &options = nlohmann::json::object());
} // namespace babelmesh
