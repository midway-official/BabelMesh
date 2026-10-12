#pragma once
#include "babelmesh/cartesian_mesh.hpp"

namespace babelmesh {
// Necessary export conditions, not a self-intersection or full CFD quality test.
struct PolyhedralAudit {
    std::size_t open_cells{};
    std::size_t area_closure_failures{};
    std::size_t point_mean_halfspace_failures{};
    std::size_t first_open_cell{static_cast<std::size_t>(-1)};
    std::size_t first_halfspace_cell{static_cast<std::size_t>(-1)};
};
struct MeshQualityAudit {
    double minimum_cell_volume{};
    double minimum_regular_hex_scaled_jacobian{};
    double maximum_non_orthogonality_degrees{};
    std::size_t maximum_non_orthogonality_face_index{static_cast<std::size_t>(-1)};
    double p99_non_orthogonality_degrees{};
    double maximum_internal_skewness{};
    double p99_internal_skewness{};
    std::size_t negative_or_zero_volume_cells{};
    std::size_t degenerate_faces{};
    std::size_t regular_hex_cells{};
    std::vector<double> cell_volumes;
    std::vector<double> non_orthogonality;
    std::vector<double> non_orthogonality_by_face;
    std::vector<double> internal_skewness;
    std::vector<double> non_orthogonality_weights;
    std::vector<double> internal_skewness_weights;
};
PolyhedralAudit auditPolyhedralShells(const HexMesh& mesh);
MeshQualityAudit auditMeshQuality(const HexMesh& mesh);
// Insert existing, incident-cell vertices on exactly collinear edges. No point
// movement, merging, or area removal. Return inserted face-ring entries.
std::size_t conformCollinearCellEdges(HexMesh& mesh, double tolerance);
void requirePolyhedralShells(const HexMesh& mesh);
std::size_t consolidateCoplanarFaces(HexMesh& mesh, double planeTolerance, double relativeWarp = 0.0);
}  // namespace babelmesh
