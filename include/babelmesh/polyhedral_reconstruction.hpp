#pragma once

#include "babelmesh/cartesian_mesh.hpp"

#include <cstddef>
#include <limits>

namespace babelmesh {

// Replace low-Jacobian or high-nonorthogonality cells with one polygonal
// pyramid per boundary face. Shared parent faces retain their geometry;
// triangular pyramid interfaces are inserted from an interior point. The
// operation is transactional and returns zero without modifying the mesh when
// selected cells cannot be safely star-decomposed.
std::size_t splitLowQualityHexesIntoPyramids(HexMesh& mesh, double minimumJacobian = 0.1,
    std::size_t maximumCells = std::numeric_limits<std::size_t>::max(),
    double maximumFaceNonOrthogonalityDegrees = std::numeric_limits<double>::infinity());

}  // namespace babelmesh
