#pragma once
#include "babelmesh/cad.hpp"
#include "babelmesh/cartesian_mesh.hpp"

namespace babelmesh {
struct SurfaceDistanceAudit {
    double mesh_to_cad_max_m{};
    double cad_to_mesh_max_m{};
    double requested_sample_spacing_m{};
    double cad_tessellation_deflection_m{};
    std::size_t mesh_to_cad_samples{};
    std::size_t cad_to_mesh_samples{};
};
// Bidirectional sampled distances, not a certified continuous Hausdorff bound.
// CAD-side samples are projected from the tessellation onto the exact B-Rep.
// The budget applies separately to each direction and exhaustion is an error.
SurfaceDistanceAudit auditSurfaceDistance(const HexMesh& mesh, const CadModel& cad,
    double sampleSpacing, std::size_t maximumSamples = 2000000);
} // namespace babelmesh
