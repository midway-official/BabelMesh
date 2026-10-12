#pragma once

#include "babelmesh/cad.hpp"
#include "babelmesh/cartesian_mesh.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <string_view>

namespace babelmesh {

enum class CellRelation {
    inside,
    outside,
    intersect,
    uncertain
};

struct CellClassification {
    CellRelation relation{CellRelation::uncertain};
    std::size_t candidate_surface_triangles{};
    std::size_t intersecting_surface_triangles{};
    std::size_t inside_samples{};
    std::size_t outside_samples{};
    std::size_t on_surface_samples{};
    bool center_inside_fluid{};
    bool center_on_surface{};
};

struct BoundaryProjection {
    std::array<double, 3> point_m{};
    double distance_m{};
    std::size_t face_index{};
};

CellClassification classifyFluidCell(
    const CadModel& cadModel,
    const std::array<double, 6>& cellBounds,
    double tolerance = 1.0e-9);

bool surfaceIntersectsBox(
    const SurfaceMesh& surface,
    const std::array<double, 6>& bounds,
    double tolerance = 1.0e-9);

bool surfaceIntersectsBoxExceptBoundary(
    const SurfaceMesh& surface,
    const std::array<double, 6>& cellBounds,
    const std::array<double, 6>& domainBounds,
    double tolerance = 1.0e-9);

HexMesh extractSeedConnectedFluid(
    HexMesh background,
    const nlohmann::json& project,
    const CadModel& cadModel);

void snapFluidBoundaryToCad(
    HexMesh& mesh,
    const CadModel& cadModel,
    double maximumSurfaceDeviation,
    bool strictQuality,
    std::size_t maximumCells = std::numeric_limits<std::size_t>::max());

std::string_view toString(CellRelation relation) noexcept;
BoundaryProjection projectPointToFluidBoundary(
    const CadModel& cadModel,
    const std::array<double, 3>& point_m,
    double tolerance = 1.0e-9);

}  // namespace babelmesh
