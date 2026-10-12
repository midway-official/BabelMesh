#include "babelmesh/fluid_domain.hpp"
#include "babelmesh/polyhedral_audit.hpp"
#include "babelmesh/polyhedral_reconstruction.hpp"

#include <BRepBndLib.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopAbs_State.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <deque>
#include <functional>
#include <numeric>
#include <queue>
#include <set>
#include <vector>

namespace babelmesh {
namespace {

bool boundsIntersect(const std::array<double, 6>& a, const std::array<double, 6>& b, double tolerance) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (a[axis + 3U] + tolerance < b[axis] || b[axis + 3U] + tolerance < a[axis]) return false;
    }
    return true;
}

bool separatedOnAxis(const std::array<std::array<double, 3>, 3>& triangle,
                     const std::array<double, 3>& half,
                     const std::array<double, 3>& axis,
                     double tolerance) {
    const double axisSquared = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if (axisSquared <= 1.0e-28) return false;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto& vertex : triangle) {
        const double projection = vertex[0] * axis[0] + vertex[1] * axis[1] + vertex[2] * axis[2];
        minimum = std::min(minimum, projection);
        maximum = std::max(maximum, projection);
    }
    const double radius = (half[0] + tolerance) * std::abs(axis[0]) +
        (half[1] + tolerance) * std::abs(axis[1]) + (half[2] + tolerance) * std::abs(axis[2]);
    return minimum > radius || maximum < -radius;
}

bool triangleIntersectsBox(const std::array<std::array<double, 3>, 3>& points,
                           const std::array<double, 6>& bounds, double tolerance) {
    std::array<double, 3> center{}, half{};
    std::array<std::array<double, 3>, 3> triangle{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        center[axis] = 0.5 * (bounds[axis] + bounds[axis + 3U]);
        half[axis] = 0.5 * (bounds[axis + 3U] - bounds[axis]);
        for (std::size_t vertex = 0; vertex < 3; ++vertex) triangle[vertex][axis] = points[vertex][axis] - center[axis];
    }

    constexpr std::array<std::array<double, 3>, 3> boxAxes{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (const auto& axis : boxAxes) if (separatedOnAxis(triangle, half, axis, tolerance)) return false;

    std::array<std::array<double, 3>, 3> edges{};
    for (std::size_t edge = 0; edge < 3; ++edge) {
        const auto& a = triangle[edge];
        const auto& b = triangle[(edge + 1U) % 3U];
        edges[edge] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    }
    for (const auto& edge : edges) for (const auto& axis : boxAxes) {
        const std::array<double, 3> cross{edge[1] * axis[2] - edge[2] * axis[1],
            edge[2] * axis[0] - edge[0] * axis[2], edge[0] * axis[1] - edge[1] * axis[0]};
        if (separatedOnAxis(triangle, half, cross, tolerance)) return false;
    }
    const std::array<double, 3> normal{edges[0][1] * edges[1][2] - edges[0][2] * edges[1][1],
        edges[0][2] * edges[1][0] - edges[0][0] * edges[1][2],
        edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0]};
    return !separatedOnAxis(triangle, half, normal, tolerance);
}

std::vector<std::size_t> surfaceCandidates(const SurfaceMesh& surface,
                                          const std::array<double, 6>& bounds, double tolerance) {
    std::vector<std::size_t> result;
    if (surface.triangles.empty()) return result;
    if (surface.bvh_nodes.empty()) {
        result.reserve(surface.triangles.size());
        for (std::size_t triangle = 0; triangle < surface.triangles.size(); ++triangle) result.push_back(triangle);
        return result;
    }
    std::vector<std::size_t> stack{0U};
    while (!stack.empty()) {
        const auto nodeIndex = stack.back();
        stack.pop_back();
        if (nodeIndex >= surface.bvh_nodes.size()) throw std::runtime_error("fluid surface BVH contains an invalid node index");
        const auto& node = surface.bvh_nodes[nodeIndex];
        if (!boundsIntersect(node.bounds, bounds, tolerance)) continue;
        if (node.count > 0U) {
            const std::size_t end = node.begin + node.count;
            if (end > surface.bvh_triangle_order.size()) throw std::runtime_error("fluid surface BVH leaf has an invalid triangle range");
            for (std::size_t item = node.begin; item < end; ++item) {
                const auto triangleIndex = surface.bvh_triangle_order[item];
                if (triangleIndex >= surface.triangles.size()) throw std::runtime_error("fluid surface BVH references an invalid triangle");
                const auto& triangle = surface.triangles[triangleIndex];
                std::array<double, 6> triangleBounds{1.0e300, 1.0e300, 1.0e300, -1.0e300, -1.0e300, -1.0e300};
                for (const auto vertex : triangle.vertices) {
                    if (vertex >= surface.vertices.size()) throw std::runtime_error("fluid surface triangle references an invalid vertex");
                    for (std::size_t axis = 0; axis < 3; ++axis) {
                        const double coordinate = surface.vertices[vertex][axis];
                        triangleBounds[axis] = std::min(triangleBounds[axis], coordinate);
                        triangleBounds[axis + 3U] = std::max(triangleBounds[axis + 3U], coordinate);
                    }
                }
                if (boundsIntersect(triangleBounds, bounds, tolerance)) result.push_back(triangleIndex);
            }
        } else {
            if (node.left >= surface.bvh_nodes.size() || node.right >= surface.bvh_nodes.size()) {
                throw std::runtime_error("fluid surface BVH internal node has an invalid child index");
            }
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }
    return result;
}

double projectLengthScale(const std::string& unit) {
    if (unit == "m") return 1.0;
    if (unit == "mm") return 1.0e-3;
    if (unit == "cm") return 1.0e-2;
    if (unit == "inch") return 0.0254;
    throw std::runtime_error("unsupported project length unit '" + unit + "'");
}

std::array<double, 6> cellBounds(const HexMesh& mesh, const std::array<std::int32_t, 8>& cell) {
    std::array<double, 6> result{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const auto vertex : cell) {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= mesh.points_m.size()) {
            throw std::runtime_error("fluid extraction cell has an invalid point index");
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            result[axis] = std::min(result[axis], mesh.points_m[static_cast<std::size_t>(vertex)][axis]);
            result[axis + 3U] = std::max(result[axis + 3U], mesh.points_m[static_cast<std::size_t>(vertex)][axis]);
        }
    }
    return result;
}

bool containsPoint(const std::array<double, 6>& bounds, const std::array<double, 3>& point, double tolerance) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (point[axis] < bounds[axis] - tolerance || point[axis] > bounds[axis + 3U] + tolerance) return false;
    }
    return true;
}

double squaredDistance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    double result = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double delta = a[axis] - b[axis];
        result += delta * delta;
    }
    return result;
}

std::array<double, 3> closestPointOnTriangle(const std::array<double, 3>& point,
                                              const std::array<double, 3>& a,
                                              const std::array<double, 3>& b,
                                              const std::array<double, 3>& c) {
    const std::array<double, 3> ab{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const std::array<double, 3> ac{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const std::array<double, 3> ap{point[0] - a[0], point[1] - a[1], point[2] - a[2]};
    const double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
    const double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
    if (d1 <= 0.0 && d2 <= 0.0) return a;
    const std::array<double, 3> bp{point[0] - b[0], point[1] - b[1], point[2] - b[2]};
    const double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
    const double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
    if (d3 >= 0.0 && d4 <= d3) return b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return {a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]};
    }
    const std::array<double, 3> cp{point[0] - c[0], point[1] - c[1], point[2] - c[2]};
    const double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
    const double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
    if (d6 >= 0.0 && d5 <= d6) return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return {a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]};
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const std::array<double, 3> bc{c[0] - b[0], c[1] - b[1], c[2] - b[2]};
        return {b[0] + w * bc[0], b[1] + w * bc[1], b[2] + w * bc[2]};
    }
    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator;
    const double w = vc * denominator;
    return {a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w};
}

double squaredDistanceToBounds(const std::array<double, 3>& point, const std::array<double, 6>& bounds) {
    double result = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double delta = point[axis] < bounds[axis] ? bounds[axis] - point[axis] :
            (point[axis] > bounds[axis + 3U] ? point[axis] - bounds[axis + 3U] : 0.0);
        result += delta * delta;
    }
    return result;
}

std::vector<std::size_t> nearestSurfaceFaces(const SurfaceMesh& surface,
                                             const std::array<double, 3>& point) {
    std::vector<std::size_t> result;
    if (surface.triangles.empty()) return result;
    if (surface.bvh_nodes.empty()) {
        double minimum = std::numeric_limits<double>::infinity();
        std::vector<std::pair<double, std::size_t>> distances;
        distances.reserve(surface.triangles.size());
        for (std::size_t index = 0; index < surface.triangles.size(); ++index) {
            const auto& triangle = surface.triangles[index];
            const auto closest = closestPointOnTriangle(point,
                surface.vertices[triangle.vertices[0]], surface.vertices[triangle.vertices[1]], surface.vertices[triangle.vertices[2]]);
            const double distance = squaredDistance(point, closest);
            minimum = std::min(minimum, distance);
            distances.emplace_back(distance, triangle.face_index);
        }
        const double threshold = minimum + surface.linear_deflection_m * surface.linear_deflection_m;
        for (const auto& candidate : distances) if (candidate.first <= threshold) result.push_back(candidate.second);
    } else {
        using Entry = std::pair<double, std::size_t>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pending;
        pending.emplace(squaredDistanceToBounds(point, surface.bvh_nodes.front().bounds), 0U);
        double minimum = std::numeric_limits<double>::infinity();
        std::vector<std::pair<double, std::size_t>> faceDistances;
        while (!pending.empty()) {
            const auto [boundDistance, nodeIndex] = pending.top();
            pending.pop();
            if (boundDistance > minimum + surface.linear_deflection_m * surface.linear_deflection_m) break;
            const auto& node = surface.bvh_nodes[nodeIndex];
            if (node.count > 0U) {
                const std::size_t end = node.begin + node.count;
                if (end > surface.bvh_triangle_order.size()) throw std::runtime_error("nearest-surface BVH leaf has an invalid triangle range");
                for (std::size_t item = node.begin; item < end; ++item) {
                    const auto triangleIndex = surface.bvh_triangle_order[item];
                    const auto& triangle = surface.triangles.at(triangleIndex);
                    const auto closest = closestPointOnTriangle(point,
                        surface.vertices.at(triangle.vertices[0]), surface.vertices.at(triangle.vertices[1]),
                        surface.vertices.at(triangle.vertices[2]));
                    const double distance = squaredDistance(point, closest);
                    minimum = std::min(minimum, distance);
                    faceDistances.emplace_back(distance, triangle.face_index);
                }
            } else {
                if (node.left >= surface.bvh_nodes.size() || node.right >= surface.bvh_nodes.size()) throw std::runtime_error("nearest-surface BVH has an invalid child reference");
                const double leftDistance = squaredDistanceToBounds(point, surface.bvh_nodes[node.left].bounds);
                const double rightDistance = squaredDistanceToBounds(point, surface.bvh_nodes[node.right].bounds);
                if (leftDistance <= minimum + surface.linear_deflection_m * surface.linear_deflection_m) pending.emplace(leftDistance, node.left);
                if (rightDistance <= minimum + surface.linear_deflection_m * surface.linear_deflection_m) pending.emplace(rightDistance, node.right);
            }
        }
        const double threshold = minimum + surface.linear_deflection_m * surface.linear_deflection_m;
        for (const auto& candidate : faceDistances) if (candidate.first <= threshold) result.push_back(candidate.second);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool faceCenterHasFluidPassage(const CadModel& cadModel, const HexMesh& mesh,
                               const MeshFace& face, double tolerance) {
    if (face.vertices.size() < 3U) return false;
    std::array<double, 3> center{};
    std::array<double, 3> minimum{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    std::array<double, 3> maximum{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const auto vertex : face.vertices) {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= mesh.points_m.size()) throw std::runtime_error("face-center connectivity check found an invalid point index");
        const auto& point = mesh.points_m[static_cast<std::size_t>(vertex)];
        for (std::size_t axis = 0; axis < 3; ++axis) {
            center[axis] += point[axis] / static_cast<double>(face.vertices.size());
            minimum[axis] = std::min(minimum[axis], point[axis]);
            maximum[axis] = std::max(maximum[axis], point[axis]);
        }
    }
    std::size_t normalAxis = 3U;
    double faceScale = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (maximum[axis] - minimum[axis] <= tolerance * 2.0) normalAxis = axis;
        else faceScale = std::max(faceScale, maximum[axis] - minimum[axis]);
    }
    if (normalAxis >= 3U) return false;
    const double epsilon = std::max(tolerance * 16.0, faceScale * 1.0e-6);
    BRepClass3d_SolidClassifier classifier(cadModel.fluidDomain());
    auto stateAt = [&](double sign) {
        auto point = center;
        point[normalAxis] += sign * epsilon;
        classifier.Perform(gp_Pnt(point[0], point[1], point[2]), tolerance);
        return classifier.State();
    };
    try {
        return stateAt(-1.0) == TopAbs_IN && stateAt(1.0) == TopAbs_IN;
    } catch (const Standard_Failure& error) {
        const char* message = error.GetMessageString();
        throw std::runtime_error(std::string("OCCT face passage classification failed: ") + (message ? message : "unknown OCCT failure"));
    }
}

}  // namespace

CellClassification classifyFluidCell(const CadModel& cadModel,
                                     const std::array<double, 6>& cellBounds,
                                     double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) throw std::runtime_error("cell classification tolerance must be finite and nonnegative");
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(cellBounds[axis]) || !std::isfinite(cellBounds[axis + 3U]) || cellBounds[axis + 3U] <= cellBounds[axis]) {
            throw std::runtime_error("cell classification requires finite, positive bounds");
        }
    }
    const auto& surface = cadModel.fluidSurfaceMesh();
    const auto candidates = surfaceCandidates(surface, cellBounds, tolerance);
    CellClassification result;
    result.candidate_surface_triangles = candidates.size();
    for (const auto triangleIndex : candidates) {
        const auto& triangle = surface.triangles[triangleIndex];
        std::array<std::array<double, 3>, 3> points{};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) points[vertex] = surface.vertices[triangle.vertices[vertex]];
        if (triangleIntersectsBox(points, cellBounds, tolerance)) ++result.intersecting_surface_triangles;
    }

    BRepClass3d_SolidClassifier classifier(cadModel.fluidDomain());
    std::array<std::array<double, 3>, 9> samples{};
    for (std::size_t corner = 0; corner < 8; ++corner) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            samples[corner][axis] = ((corner & (1U << axis)) != 0U) ? cellBounds[axis + 3U] : cellBounds[axis];
        }
    }
    for (std::size_t axis = 0; axis < 3; ++axis) samples[8][axis] = 0.5 * (cellBounds[axis] + cellBounds[axis + 3U]);
    auto classifySample = [&](std::size_t index) {
        const auto& sample = samples[index];
        classifier.Perform(gp_Pnt(sample[0], sample[1], sample[2]), tolerance);
        const TopAbs_State state = classifier.State();
        if (state == TopAbs_IN) ++result.inside_samples;
        else if (state == TopAbs_OUT) ++result.outside_samples;
        else if (state == TopAbs_ON) ++result.on_surface_samples;
        if (index == 8U) {
            result.center_inside_fluid = state == TopAbs_IN;
            result.center_on_surface = state == TopAbs_ON;
        }
        return state;
    };
    TopAbs_State centerState = TopAbs_UNKNOWN;
    try {
        centerState = classifySample(8U);
        if (result.intersecting_surface_triangles == 0U && !candidates.empty()) {
            for (std::size_t index = 0; index < 8U; ++index) (void)classifySample(index);
        }
    } catch (const Standard_Failure& error) {
        const char* message = error.GetMessageString();
        throw std::runtime_error(std::string("OCCT fluid-domain point classification failed: ") + (message ? message : "unknown OCCT failure"));
    }

    if (result.intersecting_surface_triangles > 0U ||
        (result.inside_samples > 0U && result.outside_samples > 0U)) {
        result.relation = CellRelation::intersect;
    } else if (result.on_surface_samples > 0U || centerState == TopAbs_UNKNOWN) {
        result.relation = CellRelation::uncertain;
    } else if (result.inside_samples > 0U && result.outside_samples == 0U) {
        result.relation = CellRelation::inside;
    } else if (result.outside_samples > 0U && result.inside_samples == 0U) {
        result.relation = CellRelation::outside;
    }
    return result;
}

std::string_view toString(CellRelation relation) noexcept {
    switch (relation) {
        case CellRelation::inside: return "inside";
        case CellRelation::outside: return "outside";
        case CellRelation::intersect: return "intersect";
        case CellRelation::uncertain: return "uncertain";
    }
    return "uncertain";
}

bool surfaceIntersectsBox(const SurfaceMesh& surface,
                          const std::array<double, 6>& bounds,
                          double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) throw std::runtime_error("surface intersection tolerance must be finite and nonnegative");
    const auto candidates = surfaceCandidates(surface, bounds, tolerance);
    for (const auto triangleIndex : candidates) {
        const auto& triangle = surface.triangles[triangleIndex];
        std::array<std::array<double, 3>, 3> points{};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) points[vertex] = surface.vertices[triangle.vertices[vertex]];
        if (triangleIntersectsBox(points, bounds, tolerance)) return true;
    }
    return false;
}

bool surfaceIntersectsBoxExceptBoundary(const SurfaceMesh& surface,
                                       const std::array<double, 6>& cellBounds,
                                       const std::array<double, 6>& domainBounds,
                                       double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) throw std::runtime_error("surface intersection tolerance must be finite and nonnegative");
    const auto candidates = surfaceCandidates(surface, cellBounds, tolerance);
    for (const auto triangleIndex : candidates) {
        const auto& triangle = surface.triangles[triangleIndex];
        std::array<std::array<double, 3>, 3> points{};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) points[vertex] = surface.vertices[triangle.vertices[vertex]];
        if (!triangleIntersectsBox(points, cellBounds, tolerance)) continue;
        bool liesOnDomainBoundary = false;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const bool onMinimum = std::all_of(points.begin(), points.end(), [&](const auto& point) {
                return std::abs(point[axis] - domainBounds[axis]) <= tolerance;
            });
            const bool onMaximum = std::all_of(points.begin(), points.end(), [&](const auto& point) {
                return std::abs(point[axis] - domainBounds[axis + 3U]) <= tolerance;
            });
            liesOnDomainBoundary = liesOnDomainBoundary || onMinimum || onMaximum;
        }
        if (!liesOnDomainBoundary) return true;
    }
    return false;
}

HexMesh extractSeedConnectedFluid(HexMesh background,
                                  const nlohmann::json& project,
                                  const CadModel& cadModel) {
    if (background.cells.empty() || background.faces.empty()) throw std::runtime_error("fluid extraction requires a nonempty background mesh");
    if (background.cells.size() != background.cell_levels.size() || background.cells.size() != background.strict_hex.size()) {
        throw std::runtime_error("fluid extraction received inconsistent per-cell arrays");
    }
    const double scale = projectLengthScale(project.at("units").at("length").get<std::string>());
    const double tolerance = std::max(1.0e-10, 1.0e-9 * (background.bounds_m[3] - background.bounds_m[0]));
    const std::size_t cellCount = background.cells.size();
    std::vector<std::array<double, 6>> bounds(cellCount);
    std::vector<CellClassification> classifications(cellCount);
    std::vector<bool> retained(cellCount, false);
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
        bounds[cell] = cellBounds(background, background.cells[cell]);
        classifications[cell] = classifyFluidCell(cadModel, bounds[cell], tolerance);
        if (classifications[cell].relation == CellRelation::intersect) ++background.cad_intersect_cells;
        if (classifications[cell].relation == CellRelation::uncertain) ++background.cad_uncertain_cells;
        retained[cell] = classifications[cell].center_inside_fluid;
    }

    const auto& seeds = project.at("fluid_domain").at("seeds");
    if (!seeds.is_array() || seeds.empty()) throw std::runtime_error("fluid extraction requires at least one configured fluid seed");
    std::vector<std::size_t> seedCells;
    std::vector<std::array<double, 3>> seedPoints;
    for (const auto& seedValue : seeds) {
        if (!seedValue.is_array() || seedValue.size() != 3U) throw std::runtime_error("fluid seed must have three coordinates");
        std::array<double, 3> seed{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            seed[axis] = seedValue[axis].get<double>() * scale;
            if (!std::isfinite(seed[axis])) throw std::runtime_error("fluid seed contains a non-finite coordinate");
        }
        std::size_t selected = cellCount;
        double bestDistance = std::numeric_limits<double>::infinity();
        for (std::size_t cell = 0; cell < cellCount; ++cell) {
            if (!containsPoint(bounds[cell], seed, tolerance)) continue;
            const std::array<double, 3> center{0.5 * (bounds[cell][0] + bounds[cell][3]),
                0.5 * (bounds[cell][1] + bounds[cell][4]), 0.5 * (bounds[cell][2] + bounds[cell][5])};
            const double distance = squaredDistance(center, seed);
            const double preference = classifications[cell].center_inside_fluid ? 0.0 : 1.0e100;
            if (preference + distance < bestDistance) {
                selected = cell;
                bestDistance = preference + distance;
            }
        }
        if (selected == cellCount) {
            throw std::runtime_error("fluid seed lies outside the generated background mesh; seed= (" +
                std::to_string(seed[0]) + ", " + std::to_string(seed[1]) + ", " + std::to_string(seed[2]) + ")");
        }
        if (!retained[selected]) {
            retained[selected] = true;
            ++background.cad_seed_forced_cells;
        }
        seedCells.push_back(selected);
        seedPoints.push_back(seed);
    }

    std::vector<std::vector<std::size_t>> adjacency(cellCount);
    for (const auto& face : background.faces) {
        if (face.neighbour < 0) continue;
        if (face.owner < 0 || static_cast<std::size_t>(face.owner) >= cellCount ||
            static_cast<std::size_t>(face.neighbour) >= cellCount) {
            throw std::runtime_error("background face has an invalid cell reference during fluid extraction");
        }
        const auto owner = static_cast<std::size_t>(face.owner);
        const auto neighbour = static_cast<std::size_t>(face.neighbour);
        if (!retained[owner] || !retained[neighbour]) continue;
        if (faceCenterHasFluidPassage(cadModel, background, face, tolerance)) {
            adjacency[owner].push_back(neighbour);
            adjacency[neighbour].push_back(owner);
        }
    }

    std::vector<bool> connected(cellCount, false);
    std::deque<std::size_t> queue;
    for (const auto seedCell : seedCells) {
        if (!connected[seedCell]) {
            connected[seedCell] = true;
            queue.push_back(seedCell);
        }
    }
    while (!queue.empty()) {
        const auto current = queue.front();
        queue.pop_front();
        for (const auto neighbour : adjacency[current]) if (!connected[neighbour]) {
            connected[neighbour] = true;
            queue.push_back(neighbour);
        }
    }
    retained = std::move(connected);
    const std::size_t retainedCount = static_cast<std::size_t>(std::count(retained.begin(), retained.end(), true));
    if (retainedCount == 0U) throw std::runtime_error("fluid seed flood fill selected no cells");

    std::vector<std::int32_t> oldToNew(cellCount, -1);
    std::vector<std::array<std::int32_t, 8>> newCells;
    std::vector<std::int32_t> newLevels;
    std::vector<bool> newStrictHex;
    std::vector<double> newVolumes, newAspects;
    newCells.reserve(retainedCount); newLevels.reserve(retainedCount); newStrictHex.reserve(retainedCount);
    for (std::size_t cell = 0; cell < cellCount; ++cell) if (retained[cell]) {
        oldToNew[cell] = static_cast<std::int32_t>(newCells.size());
        newCells.push_back(background.cells[cell]);
        newLevels.push_back(background.cell_levels[cell]);
        newStrictHex.push_back(background.strict_hex[cell]);
        if (cell < background.cell_volumes_m3.size()) newVolumes.push_back(background.cell_volumes_m3[cell]);
        if (cell < background.cell_aspect_ratios.size()) newAspects.push_back(background.cell_aspect_ratios[cell]);
    }

    auto cadPatch = std::find_if(background.patches.begin(), background.patches.end(),
        [](const auto& patch) { return patch.first == "cad_boundary"; });
    if (cadPatch == background.patches.end()) {
        background.patches.emplace_back("cad_boundary", "wall");
        cadPatch = std::prev(background.patches.end());
    }
    const auto cadPatchIndex = static_cast<std::int32_t>(std::distance(background.patches.begin(), cadPatch));
    std::vector<MeshFace> newFaces;
    newFaces.reserve(background.faces.size());
    for (auto face : background.faces) {
        const auto oldOwner = static_cast<std::size_t>(face.owner);
        const bool ownerRetained = oldOwner < cellCount && retained[oldOwner];
        const bool neighbourExists = face.neighbour >= 0;
        const auto oldNeighbour = neighbourExists ? static_cast<std::size_t>(face.neighbour) : cellCount;
        const bool neighbourRetained = neighbourExists && oldNeighbour < cellCount && retained[oldNeighbour];
        if (!ownerRetained && !neighbourRetained) continue;
        if (ownerRetained && neighbourRetained) {
            face.owner = oldToNew[oldOwner];
            face.neighbour = oldToNew[oldNeighbour];
            face.patch = -1;
        } else if (ownerRetained) {
            face.owner = oldToNew[oldOwner];
            face.neighbour = -1;
            if (neighbourExists) face.patch = cadPatchIndex;
        } else {
            face.owner = oldToNew[oldNeighbour];
            face.neighbour = -1;
            face.patch = cadPatchIndex;
            std::reverse(face.vertices.begin(), face.vertices.end());
        }
        newFaces.push_back(std::move(face));
    }

    for (auto& source : background.sizing_source_markers) {
        if (source.second.size() != cellCount) throw std::runtime_error("sizing source markers do not match background cell count");
        std::vector<std::uint8_t> filtered;
        filtered.reserve(retainedCount);
        std::size_t marked = 0;
        for (std::size_t cell = 0; cell < cellCount; ++cell) if (retained[cell]) {
            filtered.push_back(source.second[cell]);
            if (source.second[cell] != 0U) ++marked;
        }
        source.second = std::move(filtered);
        background.source_intersection_counts[source.first] = marked;
        const auto detail = background.sizing_source_details.find(source.first);
        if (detail != background.sizing_source_details.end()) detail->second["marked_cell_count"] = marked;
    }

    background.cells = std::move(newCells);
    background.cell_levels = std::move(newLevels);
    background.strict_hex = std::move(newStrictHex);
    background.cell_volumes_m3 = std::move(newVolumes);
    background.cell_aspect_ratios = std::move(newAspects);
    background.faces = std::move(newFaces);
    std::vector<bool> usedPoints(background.points_m.size(), false);
    for (const auto& face : background.faces) for (const auto vertex : face.vertices) {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= usedPoints.size()) throw std::runtime_error("extracted face has an invalid point index");
        usedPoints[static_cast<std::size_t>(vertex)] = true;
    }
    std::vector<std::int32_t> pointRemap(background.points_m.size(), -1);
    std::vector<std::array<double, 3>> compactPoints;
    for (std::size_t point = 0; point < usedPoints.size(); ++point) if (usedPoints[point]) {
        pointRemap[point] = static_cast<std::int32_t>(compactPoints.size());
        compactPoints.push_back(background.points_m[point]);
    }
    for (auto& face : background.faces) for (auto& vertex : face.vertices) vertex = pointRemap[static_cast<std::size_t>(vertex)];
    for (auto& cell : background.cells) for (auto& vertex : cell) vertex = pointRemap[static_cast<std::size_t>(vertex)];
    background.points_m = std::move(compactPoints);

    std::vector<TopoDS_Shape> domainComponents;
    if (cadModel.fluidDomain().ShapeType() == TopAbs_SOLID) domainComponents.push_back(cadModel.fluidDomain());
    else for (TopExp_Explorer explorer(cadModel.fluidDomain(), TopAbs_SOLID); explorer.More(); explorer.Next()) {
        domainComponents.push_back(explorer.Current());
    }
    if (domainComponents.empty()) domainComponents.push_back(cadModel.fluidDomain());
    for (const auto& component : domainComponents) {
        bool selectedBySeed = false;
        BRepClass3d_SolidClassifier componentClassifier(component);
        for (const auto& seed : seedPoints) {
            componentClassifier.Perform(gp_Pnt(seed[0], seed[1], seed[2]), tolerance);
            if (componentClassifier.State() == TopAbs_IN || componentClassifier.State() == TopAbs_ON) {
                selectedBySeed = true;
                break;
            }
        }
        if (!selectedBySeed) continue;
        GProp_GProps volumeProperties;
        BRepGProp::VolumeProperties(component, volumeProperties);
        background.cad_fluid_domain_volume_m3 += std::abs(volumeProperties.Mass());
        ++background.cad_selected_components;
    }
    const double meshVolume = std::accumulate(background.cell_volumes_m3.begin(), background.cell_volumes_m3.end(), 0.0);
    if (background.cad_fluid_domain_volume_m3 <= 0.0) throw std::runtime_error("CAD fluid domain has a nonpositive volume");
    background.cad_mesh_volume_relative_error = std::abs(meshVolume - background.cad_fluid_domain_volume_m3) /
        background.cad_fluid_domain_volume_m3;
    background.cad_removed_cells = cellCount - retainedCount;
    background.cad_clipping_applied = true;
    return background;
}

void snapFluidBoundaryToCad(HexMesh& mesh,
                            const CadModel& cadModel,
                            double maximumSurfaceDeviation,
                            bool strictQuality,
                            std::size_t maximumCells) {
    if (!std::isfinite(maximumSurfaceDeviation) || maximumSurfaceDeviation < 0.0) {
        throw std::runtime_error("maximum CAD surface deviation must be finite and nonnegative");
    }
    mesh.cad_surface_deviation_tolerance_m = maximumSurfaceDeviation;
    const auto patch = std::find_if(mesh.patches.begin(), mesh.patches.end(),
        [](const auto& item) { return item.first == "cad_boundary"; });
    if (patch == mesh.patches.end()) return;
    const auto patchIndex = static_cast<std::int32_t>(std::distance(mesh.patches.begin(), patch));
    std::vector<TopoDS_Face> cadFaces;
    for (TopExp_Explorer explorer(cadModel.fluidDomain(), TopAbs_FACE); explorer.More(); explorer.Next()) {
        cadFaces.push_back(TopoDS::Face(explorer.Current()));
    }
    if (cadFaces.empty()) throw std::runtime_error("CAD snapping requires at least one trimmed fluid-boundary face");
    struct CadFeatureEdge {
        TopoDS_Edge edge;
        std::vector<std::size_t> incidentFaces;
        Bnd_Box bounds;
    };
    std::vector<CadFeatureEdge> cadFeatureEdges;
    for (std::size_t faceIndex = 0; faceIndex < cadFaces.size(); ++faceIndex) {
        for (TopExp_Explorer explorer(cadFaces[faceIndex], TopAbs_EDGE); explorer.More(); explorer.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(explorer.Current());
            auto found = std::find_if(cadFeatureEdges.begin(), cadFeatureEdges.end(), [&](const CadFeatureEdge& item) {
                return item.edge.IsSame(edge);
            });
            if (found == cadFeatureEdges.end()) {
                cadFeatureEdges.push_back({edge, {faceIndex}, Bnd_Box{}});
            } else {
                if (std::find(found->incidentFaces.begin(), found->incidentFaces.end(), faceIndex) == found->incidentFaces.end()) {
                    found->incidentFaces.push_back(faceIndex);
                }
            }
        }
    }
    cadFeatureEdges.erase(std::remove_if(cadFeatureEdges.begin(), cadFeatureEdges.end(), [](const CadFeatureEdge& edge) {
        return edge.incidentFaces.size() < 2U;
    }), cadFeatureEdges.end());
    for (auto& edge : cadFeatureEdges) std::sort(edge.incidentFaces.begin(), edge.incidentFaces.end());
    auto projectToCadFace = [&](const std::array<double, 3>& point, std::size_t faceIndex) {
        if (faceIndex >= cadFaces.size()) throw std::runtime_error("CAD face projection received an invalid face index");
        const TopoDS_Vertex query = BRepBuilderAPI_MakeVertex(gp_Pnt(point[0], point[1], point[2])).Vertex();
        BRepExtrema_DistShapeShape distance(query, cadFaces[faceIndex]);
        distance.Perform();
        if (!distance.IsDone() || distance.NbSolution() == 0) {
            throw std::runtime_error("OCCT failed to project a boundary point to its associated CAD face");
        }
        const gp_Pnt projected = distance.PointOnShape2(1);
        return BoundaryProjection{{projected.X(), projected.Y(), projected.Z()}, distance.Value(), faceIndex};
    };
    std::set<std::int32_t> boundaryVertices;
    for (const auto& face : mesh.faces) if (face.neighbour < 0 && face.patch == patchIndex) {
        boundaryVertices.insert(face.vertices.begin(), face.vertices.end());
    }
    if (boundaryVertices.empty()) return;

    struct TargetPoint {
        std::int32_t vertex{};
        std::array<double, 3> original{};
        std::array<double, 3> target{};
        std::size_t faceIndex{};
    };
    std::vector<TargetPoint> targets;
    targets.reserve(boundaryVertices.size());
    for (const auto vertex : boundaryVertices) {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= mesh.points_m.size()) throw std::runtime_error("CAD snapping found an invalid boundary point index");
        const auto point = mesh.points_m[static_cast<std::size_t>(vertex)];
        const auto projection = projectPointToFluidBoundary(cadModel, point);
        targets.push_back({vertex, point, projection.point_m, projection.face_index});
    }

    // A staircase wall polygon can straddle a CAD crease even when every one
    // of its vertices has been projected to a trimmed CAD face.  Split such a
    // polygon at the associated CAD edge so the emitted boundary is piecewise
    // conformal to the CAD faces instead of bridging the crease with a warped
    // quad.  The edge point is also inserted into every incident mesh face so
    // the resulting polyhedral cells remain closed.
    const std::size_t noFace = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> cadFaceAtVertex(mesh.points_m.size(), noFace);
    std::vector<std::int32_t> targetAtVertex(mesh.points_m.size(), -1);
    for (std::size_t target = 0; target < targets.size(); ++target) {
        const auto vertex = static_cast<std::size_t>(targets[target].vertex);
        cadFaceAtVertex[vertex] = targets[target].faceIndex;
        targetAtVertex[vertex] = static_cast<std::int32_t>(target);
    }
    std::vector<std::set<std::size_t>> supportedFaces(mesh.points_m.size());
    std::vector<std::size_t> supportAtFace(mesh.faces.size(), noFace);
    std::vector<double> captureRadius(mesh.points_m.size(), 0.0);
    const bool planarDomain = std::all_of(cadFaces.begin(), cadFaces.end(), [](const auto& face) {
        return BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane;
    });
    for (std::size_t meshFace = 0; meshFace < mesh.faces.size(); ++meshFace) {
        const auto& face = mesh.faces[meshFace];
        if (face.neighbour >= 0 || face.patch != patchIndex || face.vertices.empty()) continue;
        std::map<std::size_t, std::size_t> counts;
        for (const auto vertex : face.vertices) {
            if (vertex < 0 || static_cast<std::size_t>(vertex) >= cadFaceAtVertex.size()) continue;
            const auto cadFace = cadFaceAtVertex[static_cast<std::size_t>(vertex)];
            if (cadFace != noFace) ++counts[cadFace];
        }
        if (counts.empty()) continue;
        auto support = std::max_element(counts.begin(), counts.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; })->first;
        if (planarDomain) {
            // A nearest-point tie at a sharp edge can select a perpendicular
            // CAD face. Associate the wall FACE first, using its outward
            // orientation, then constrain its vertices to shared features.
            std::array<double, 3> center{}, normal{};
            const auto& origin = mesh.points_m[static_cast<std::size_t>(face.vertices.front())];
            for (std::size_t j = 0; j < face.vertices.size(); ++j) {
                const auto v = face.vertices[j];
                const auto w = face.vertices[(j + 1U) % face.vertices.size()];
                const auto& a = mesh.points_m[static_cast<std::size_t>(v)];
                const auto& b = mesh.points_m[static_cast<std::size_t>(w)];
                const double length = std::sqrt(squaredDistance(a, b));
                captureRadius[static_cast<std::size_t>(v)] = std::max(captureRadius[static_cast<std::size_t>(v)], 2.0 * length);
                captureRadius[static_cast<std::size_t>(w)] = std::max(captureRadius[static_cast<std::size_t>(w)], 2.0 * length);
                for (std::size_t k = 0; k < 3; ++k) {
                    center[k] += a[k] / static_cast<double>(face.vertices.size());
                    normal[k] += (a[(k+1U)%3U]-origin[(k+1U)%3U])*(b[(k+2U)%3U]-origin[(k+2U)%3U]) -
                        (a[(k+2U)%3U]-origin[(k+2U)%3U])*(b[(k+1U)%3U]-origin[(k+1U)%3U]);
                }
            }
            double closest = std::numeric_limits<double>::infinity();
            for (std::size_t candidate = 0; candidate < cadFaces.size(); ++candidate) {
                auto n = BRepAdaptor_Surface(cadFaces[candidate]).Plane().Axis().Direction();
                if (cadFaces[candidate].Orientation() == TopAbs_REVERSED) n.Reverse();
                const double alignment = normal[0]*n.X() + normal[1]*n.Y() + normal[2]*n.Z();
                if (alignment <= 1.0e-10 * std::hypot(normal[0], normal[1], normal[2])) continue;
                const auto projection = projectToCadFace(center, candidate);
                if (projection.distance_m < closest) { closest = projection.distance_m; support = candidate; }
            }
            if (!std::isfinite(closest)) throw std::runtime_error("CAD planar wall has no orientation-compatible support face");
            supportAtFace[meshFace] = support;
        }
        for (const auto vertex : face.vertices) if (vertex >= 0 && static_cast<std::size_t>(vertex) < supportedFaces.size()) {
            supportedFaces[static_cast<std::size_t>(vertex)].insert(support);
        }
    }
    // Planar supports have an orientation-aware face association and a local
    // edge-length capture radius. The legacy curved-face path still uses its
    // conservative capture bound until general curved feature constraints are
    // implemented. Neither capture bound changes the final CAD error gate.
    const double featureSnapLimit = std::max(maximumSurfaceDeviation, 1.0e-12);
    for (std::size_t vertex = 0; vertex < supportedFaces.size(); ++vertex) {
        if (planarDomain && supportedFaces[vertex].size() == 1U && targetAtVertex[vertex] >= 0) {
            auto& target = targets[static_cast<std::size_t>(targetAtVertex[vertex])];
            const auto projection = projectToCadFace(target.original, *supportedFaces[vertex].begin());
            target.target = projection.point_m;
            target.faceIndex = projection.face_index;
            cadFaceAtVertex[vertex] = projection.face_index;
        }
        if (supportedFaces[vertex].size() < 2U || targetAtVertex[vertex] < 0) continue;
        const auto targetIndex = static_cast<std::size_t>(targetAtVertex[vertex]);
        const auto& support = supportedFaces[vertex];
        const TopoDS_Face& firstFace = cadFaces[*support.begin()];
        const gp_Pnt queryPoint(targets[targetIndex].original[0], targets[targetIndex].original[1], targets[targetIndex].original[2]);
        const TopoDS_Vertex query = BRepBuilderAPI_MakeVertex(queryPoint).Vertex();
        BoundaryProjection best;
        best.distance_m = std::numeric_limits<double>::infinity();
        try {
            for (TopExp_Explorer edgeExplorer(firstFace, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next()) {
                const TopoDS_Edge edge = TopoDS::Edge(edgeExplorer.Current());
                bool common = true;
                for (auto other = std::next(support.begin()); other != support.end() && common; ++other) {
                    bool present = false;
                    for (TopExp_Explorer candidate(cadFaces[*other], TopAbs_EDGE); candidate.More(); candidate.Next()) {
                        if (edge.IsSame(TopoDS::Edge(candidate.Current()))) { present = true; break; }
                    }
                    common = present;
                }
                if (!common) continue;
                BRepExtrema_DistShapeShape distance(query, edge);
                distance.Perform();
                if (!distance.IsDone()) throw std::runtime_error("OCCT failed to project a boundary vertex to a shared CAD edge");
                if (distance.Value() < best.distance_m) {
                    const gp_Pnt projected = distance.PointOnShape2(1);
                    best = {{projected.X(), projected.Y(), projected.Z()}, distance.Value(), *support.begin()};
                }
            }
            if (!std::isfinite(best.distance_m)) {
                for (TopExp_Explorer vertexExplorer(firstFace, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next()) {
                    const TopoDS_Vertex featureVertex = TopoDS::Vertex(vertexExplorer.Current());
                    bool common = true;
                    for (auto other = std::next(support.begin()); other != support.end() && common; ++other) {
                        bool present = false;
                        for (TopExp_Explorer candidate(cadFaces[*other], TopAbs_VERTEX); candidate.More(); candidate.Next()) {
                            if (featureVertex.IsSame(TopoDS::Vertex(candidate.Current()))) { present = true; break; }
                        }
                        common = present;
                    }
                    if (!common) continue;
                    BRepExtrema_DistShapeShape distance(query, featureVertex);
                    distance.Perform();
                    if (!distance.IsDone()) throw std::runtime_error("OCCT failed to project a boundary vertex to a shared CAD vertex");
                    if (distance.Value() < best.distance_m) {
                        const gp_Pnt projected = distance.PointOnShape2(1);
                        best = {{projected.X(), projected.Y(), projected.Z()}, distance.Value(), *support.begin()};
                    }
                }
            }
        } catch (const Standard_Failure& error) {
            const char* message = error.GetMessageString();
            throw std::runtime_error(std::string("OCCT CAD feature association failed: ") + (message ? message : "unknown OCCT failure"));
        }
        const double allowedCapture = planarDomain ? captureRadius[vertex] : featureSnapLimit;
        if (planarDomain && (!std::isfinite(best.distance_m) || best.distance_m > allowedCapture)) {
            throw std::runtime_error("CAD planar feature capture requires local topology reconstruction at vertex " + std::to_string(vertex));
        }
        if (std::isfinite(best.distance_m) && best.distance_m <= allowedCapture) {
            targets[targetIndex].target = best.point_m;
            targets[targetIndex].faceIndex = best.face_index;
            cadFaceAtVertex[vertex] = best.face_index;
        }
    }
    using EdgeKey = std::pair<std::int32_t, std::int32_t>;
    auto edgeKey = [](std::int32_t a, std::int32_t b) -> EdgeKey {
        return a < b ? EdgeKey{a, b} : EdgeKey{b, a};
    };
    const double domainDiagonal = std::hypot(mesh.bounds_m[3] - mesh.bounds_m[0],
        mesh.bounds_m[4] - mesh.bounds_m[1], mesh.bounds_m[5] - mesh.bounds_m[2]);
    const double degenerateAreaTolerance = 1.0e-13 * std::max(1.0, domainDiagonal * domainDiagonal);
    // Only weld numerical coincidences.  A local edge may collapse at a CAD
    // feature, but a merely short edge is still physical mesh topology.
    const double featureMergeTolerance = std::max(1.0e-12,
        64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, domainDiagonal));
    std::vector<std::array<double, 3>> pointOrigins = mesh.points_m;
    std::map<EdgeKey, std::int32_t> featurePointByEdge;
    auto featurePointOnEdge = [&](std::int32_t startVertex, std::int32_t endVertex,
                                  std::size_t startFace, std::size_t endFace) -> std::int32_t {
        const EdgeKey key = edgeKey(startVertex, endVertex);
        const auto existing = featurePointByEdge.find(key);
        if (existing != featurePointByEdge.end()) return existing->second;
        if (startVertex < 0 || endVertex < 0 || static_cast<std::size_t>(startVertex) >= targetAtVertex.size() ||
            static_cast<std::size_t>(endVertex) >= targetAtVertex.size()) {
            throw std::runtime_error("CAD feature split encountered an invalid mesh edge");
        }
        const auto startTargetIndex = targetAtVertex[static_cast<std::size_t>(startVertex)];
        const auto endTargetIndex = targetAtVertex[static_cast<std::size_t>(endVertex)];
        if (startTargetIndex < 0 || endTargetIndex < 0) {
            throw std::runtime_error("CAD feature split edge endpoint has no boundary projection association");
        }
        const TargetPoint startTarget = targets[static_cast<std::size_t>(startTargetIndex)];
        const TargetPoint endTarget = targets[static_cast<std::size_t>(endTargetIndex)];
        if (startFace == endFace) {
            throw std::runtime_error("CAD feature split was requested across a mesh edge with one CAD face association");
        }
        const gp_Pnt sourceStart(startTarget.original[0], startTarget.original[1], startTarget.original[2]);
        const gp_Pnt sourceEnd(endTarget.original[0], endTarget.original[1], endTarget.original[2]);
        const TopoDS_Edge sourceEdge = BRepBuilderAPI_MakeEdge(sourceStart, sourceEnd).Edge();
        const auto& startCadFace = cadFaces[startFace];
        const auto& endCadFace = cadFaces[endFace];
        // A newly inserted feature node may move from the Cartesian edge to
        // the true B-Rep crease by at most one local surface spacing.  Final
        // CAD deviation, positive cell volume, and BabelSim topology remain
        // independent acceptance gates.
        const double featureIntersectionTolerance = std::max(featureSnapLimit, mesh.cad_surface_target_size_m);
        double bestDistance = std::numeric_limits<double>::infinity();
        double fraction = 0.5;
        gp_Pnt featurePoint;
        bool hasCommonFeatureEdge = false;
        try {
            for (TopExp_Explorer edgeExplorer(startCadFace, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next()) {
                const TopoDS_Edge cadEdge = TopoDS::Edge(edgeExplorer.Current());
                bool common = false;
                for (TopExp_Explorer otherEdgeExplorer(endCadFace, TopAbs_EDGE); otherEdgeExplorer.More(); otherEdgeExplorer.Next()) {
                    if (cadEdge.IsSame(TopoDS::Edge(otherEdgeExplorer.Current()))) { common = true; break; }
                }
                if (!common) continue;
                hasCommonFeatureEdge = true;
                BRepExtrema_DistShapeShape distance(sourceEdge, cadEdge);
                distance.Perform();
                if (!distance.IsDone()) throw std::runtime_error("OCCT failed to intersect a mesh edge with its shared CAD feature edge");
                const std::array<double, 3> direction{endTarget.original[0] - startTarget.original[0],
                    endTarget.original[1] - startTarget.original[1], endTarget.original[2] - startTarget.original[2]};
                const double lengthSquared = squaredDistance(startTarget.original, endTarget.original);
                if (!(lengthSquared > 0.0)) throw std::runtime_error("CAD feature split encountered a zero-length mesh edge");
                for (int solution = 1; solution <= distance.NbSolution(); ++solution) {
                    const gp_Pnt onMesh = distance.PointOnShape1(solution);
                    const gp_Pnt onCad = distance.PointOnShape2(solution);
                    const std::array<double, 3> pointOnMesh{onMesh.X(), onMesh.Y(), onMesh.Z()};
                    const double t = std::clamp(((pointOnMesh[0] - startTarget.original[0]) * direction[0] +
                        (pointOnMesh[1] - startTarget.original[1]) * direction[1] +
                        (pointOnMesh[2] - startTarget.original[2]) * direction[2]) / lengthSquared, 0.0, 1.0);
                    if (distance.Value() < bestDistance ||
                        (distance.Value() == bestDistance && std::abs(t - 0.5) < std::abs(fraction - 0.5))) {
                        bestDistance = distance.Value();
                        fraction = t;
                        featurePoint = onCad;
                    }
                }
            }
        } catch (const Standard_Failure& error) {
            const char* message = error.GetMessageString();
            throw std::runtime_error(std::string("OCCT CAD feature-edge intersection failed: ") + (message ? message : "unknown OCCT failure"));
        }
        if (!hasCommonFeatureEdge) {
            throw std::runtime_error("CAD feature split requested for two trimmed faces without a shared topological edge");
        }
        if (!std::isfinite(bestDistance) || bestDistance > featureIntersectionTolerance) {
            throw std::runtime_error("CAD feature edge does not intersect its associated mesh edge within tolerance; distance=" +
                std::to_string(bestDistance) + " m, mesh edge=" + std::to_string(startVertex) + "-" +
                std::to_string(endVertex) + ", CAD faces=" + std::to_string(startFace) + "-" + std::to_string(endFace) +
                ", allowed=" + std::to_string(featureIntersectionTolerance) + " m");
        }
        const BoundaryProjection projection{{featurePoint.X(), featurePoint.Y(), featurePoint.Z()}, bestDistance, startFace};
        const double mergeToleranceSquared = featureMergeTolerance * featureMergeTolerance;
        for (const auto& target : targets) {
            if ((target.vertex == startVertex || target.vertex == endVertex) &&
                squaredDistance(target.target, projection.point_m) <= mergeToleranceSquared &&
                (target.faceIndex == startFace || target.faceIndex == endFace)) {
                featurePointByEdge.emplace(key, target.vertex);
                return target.vertex;
            }
        }
        if (mesh.points_m.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("CAD feature split exceeds the supported point index range");
        }
        const auto original = std::array<double, 3>{
            startTarget.original[0] + fraction * (endTarget.original[0] - startTarget.original[0]),
            startTarget.original[1] + fraction * (endTarget.original[1] - startTarget.original[1]),
            startTarget.original[2] + fraction * (endTarget.original[2] - startTarget.original[2])};
        const auto vertex = static_cast<std::int32_t>(mesh.points_m.size());
        mesh.points_m.push_back(original);
        pointOrigins.push_back(original);
        targets.push_back({vertex, original, projection.point_m, projection.face_index});
        cadFaceAtVertex.push_back(projection.face_index);
        targetAtVertex.push_back(static_cast<std::int32_t>(targets.size() - 1U));
        featurePointByEdge.emplace(key, vertex);
        return vertex;
    };

    const std::size_t originalFaceCount = mesh.faces.size();
    std::vector<bool> faceWasSplit(originalFaceCount, false);
    std::vector<std::vector<MeshFace>> splitFaces(originalFaceCount);
    for (std::size_t faceIndex = 0; faceIndex < originalFaceCount; ++faceIndex) {
        const auto& face = mesh.faces[faceIndex];
        if (face.neighbour >= 0 || face.patch != patchIndex || face.vertices.size() < 3U) continue;
        std::vector<std::size_t> associatedFaces;
        associatedFaces.reserve(face.vertices.size());
        for (const auto vertex : face.vertices) {
            if (vertex < 0 || static_cast<std::size_t>(vertex) >= cadFaceAtVertex.size() ||
                cadFaceAtVertex[static_cast<std::size_t>(vertex)] == noFace) {
                throw std::runtime_error("CAD wall face has a vertex without a trimmed-face association");
            }
            associatedFaces.push_back(cadFaceAtVertex[static_cast<std::size_t>(vertex)]);
        }
        std::map<std::size_t, std::size_t> faceCounts;
        for (const auto associatedFace : associatedFaces) ++faceCounts[associatedFace];
        const auto majority = std::max_element(faceCounts.begin(), faceCounts.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });
        if (majority == faceCounts.end()) continue;
        const auto supportFace = supportAtFace[faceIndex] != noFace ? supportAtFace[faceIndex] : majority->first;
        const double associationTolerance = std::max(1.0e-10, 1.0e-8 * domainDiagonal);
        for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex) {
            if (associatedFaces[vertex] == supportFace) continue;
            const auto targetIndex = targetAtVertex[static_cast<std::size_t>(face.vertices[vertex])];
            if (targetIndex < 0) continue;
            const auto& target = targets[static_cast<std::size_t>(targetIndex)];
            const auto supportProjection = projectToCadFace(target.target, supportFace);
            if (supportProjection.distance_m <= associationTolerance) associatedFaces[vertex] = supportFace;
        }
        std::size_t start = face.vertices.size();
        for (std::size_t vertex = 0; vertex < associatedFaces.size(); ++vertex) {
            const std::size_t previous = (vertex + associatedFaces.size() - 1U) % associatedFaces.size();
            if (associatedFaces[vertex] != associatedFaces[previous]) { start = vertex; break; }
        }
        if (start == face.vertices.size()) continue;
        std::size_t consumed = 0U;
        std::vector<std::int32_t> interiorRing;
        while (consumed < face.vertices.size()) {
            const std::size_t runFace = associatedFaces[start];
            std::size_t end = start;
            std::size_t runLength = 1U;
            while (consumed + runLength < face.vertices.size() &&
                associatedFaces[(end + 1U) % face.vertices.size()] == runFace) {
                end = (end + 1U) % face.vertices.size();
                ++runLength;
            }
            const std::size_t previous = (start + face.vertices.size() - 1U) % face.vertices.size();
            const std::size_t next = (end + 1U) % face.vertices.size();
            const auto startCut = featurePointOnEdge(face.vertices[previous], face.vertices[start],
                associatedFaces[previous], associatedFaces[start]);
            const auto endCut = featurePointOnEdge(face.vertices[end], face.vertices[next],
                associatedFaces[end], associatedFaces[next]);
            interiorRing.push_back(startCut);
            MeshFace piece = face;
            piece.vertices.clear();
            piece.vertices.push_back(startCut);
            for (std::size_t vertex = start;; vertex = (vertex + 1U) % face.vertices.size()) {
                piece.vertices.push_back(face.vertices[vertex]);
                if (vertex == end) break;
            }
            piece.vertices.push_back(endCut);
            splitFaces[faceIndex].push_back(std::move(piece));
            consumed += runLength;
            start = next;
        }
        if (interiorRing.size() > 2U) {
            // Wedges alone leave the middle uncovered when three or more CAD
            // face runs meet one Cartesian wall polygon.  This central region
            // closes the partition and shares every feature-cut edge.
            MeshFace centerPiece = face;
            centerPiece.vertices = std::move(interiorRing);
            splitFaces[faceIndex].push_back(std::move(centerPiece));
        }
        faceWasSplit[faceIndex] = true;
        if (face.owner >= 0 && static_cast<std::size_t>(face.owner) < mesh.strict_hex.size()) {
            mesh.strict_hex[static_cast<std::size_t>(face.owner)] = false;
        }
    }

    auto signedFaceVolume = [&](const MeshFace& face) {
        if (face.vertices.size() < 3U) return 0.0L;
        const auto& origin = mesh.points_m.at(static_cast<std::size_t>(face.vertices.front()));
        long double result = 0.0L;
        for (std::size_t vertex = 1; vertex + 1U < face.vertices.size(); ++vertex) {
            const auto& a = mesh.points_m.at(static_cast<std::size_t>(face.vertices[vertex]));
            const auto& b = mesh.points_m.at(static_cast<std::size_t>(face.vertices[vertex + 1U]));
            const long double cx = static_cast<long double>(a[1]) * b[2] - static_cast<long double>(a[2]) * b[1];
            const long double cy = static_cast<long double>(a[2]) * b[0] - static_cast<long double>(a[0]) * b[2];
            const long double cz = static_cast<long double>(a[0]) * b[1] - static_cast<long double>(a[1]) * b[0];
            result += (static_cast<long double>(origin[0]) * cx + static_cast<long double>(origin[1]) * cy +
                static_cast<long double>(origin[2]) * cz) / 6.0L;
        }
        return result;
    };
    const long double splitVolumeTolerance = 1.0e-12L * std::max(1.0, domainDiagonal * domainDiagonal * domainDiagonal);
    for (std::size_t faceIndex = 0; faceIndex < originalFaceCount; ++faceIndex) if (faceWasSplit[faceIndex]) {
        const long double originalContribution = signedFaceVolume(mesh.faces[faceIndex]);
        long double splitContribution = 0.0L;
        for (const auto& piece : splitFaces[faceIndex]) splitContribution += signedFaceVolume(piece);
        if (std::abs(splitContribution - originalContribution) > splitVolumeTolerance) {
            throw std::runtime_error("CAD feature face partition changed signed volume at boundary face " +
                std::to_string(faceIndex) + "; pieces=" + std::to_string(splitFaces[faceIndex].size()) +
                ", original=" + std::to_string(static_cast<double>(originalContribution)) +
                ", split=" + std::to_string(static_cast<double>(splitContribution)) +
                " m^3; feature cuts do not form a complete partition of the original face");
        }
    }

    if (!featurePointByEdge.empty()) {
        std::vector<MeshFace> conformingFaces;
        conformingFaces.reserve(mesh.faces.size() + featurePointByEdge.size());
        auto removeConsecutiveDuplicateVertices = [](std::vector<std::int32_t>& ring) {
            std::vector<std::int32_t> compact;
            compact.reserve(ring.size());
            for (const auto vertex : ring) if (compact.empty() || compact.back() != vertex) compact.push_back(vertex);
            if (compact.size() > 1U && compact.front() == compact.back()) compact.pop_back();
            ring = std::move(compact);
        };
        for (std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
            if (faceWasSplit[faceIndex]) {
                for (auto& splitFace : splitFaces[faceIndex]) {
                    removeConsecutiveDuplicateVertices(splitFace.vertices);
                    if (splitFace.vertices.size() >= 3U) conformingFaces.push_back(std::move(splitFace));
                }
                continue;
            }
            MeshFace face = std::move(mesh.faces[faceIndex]);
            std::vector<std::int32_t> expanded;
            expanded.reserve(face.vertices.size() + 1U);
            bool insertedFeaturePoint = false;
            for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex) {
                const auto startVertex = face.vertices[vertex];
                const auto endVertex = face.vertices[(vertex + 1U) % face.vertices.size()];
                expanded.push_back(startVertex);
                const auto feature = featurePointByEdge.find(edgeKey(startVertex, endVertex));
                if (feature != featurePointByEdge.end()) {
                    expanded.push_back(feature->second);
                    insertedFeaturePoint = true;
                }
            }
            face.vertices = std::move(expanded);
            removeConsecutiveDuplicateVertices(face.vertices);
            if (face.vertices.size() < 3U) throw std::runtime_error("CAD feature split collapsed a non-feature mesh face");
            if (insertedFeaturePoint) {
                if (face.owner >= 0 && static_cast<std::size_t>(face.owner) < mesh.strict_hex.size()) {
                    mesh.strict_hex[static_cast<std::size_t>(face.owner)] = false;
                }
                if (face.neighbour >= 0 && static_cast<std::size_t>(face.neighbour) < mesh.strict_hex.size()) {
                    mesh.strict_hex[static_cast<std::size_t>(face.neighbour)] = false;
                }
            }
            conformingFaces.push_back(std::move(face));
        }
        mesh.faces = std::move(conformingFaces);
    }

    // On smooth curved CAD patches, redistribute projected surface vertices
    // tangentially before applying their displacement to the volume mesh.
    // Nodes supported by multiple CAD faces stay fixed on their shared feature.
    // Every proposed node remains constrained to its associated B-Rep face.
    if (!planarDomain) {
        std::vector<std::int32_t> targetForVertex(mesh.points_m.size(), -1);
        for (std::size_t target = 0; target < targets.size(); ++target) {
            const auto vertex = static_cast<std::size_t>(targets[target].vertex);
            if (vertex < targetForVertex.size()) targetForVertex[vertex] = static_cast<std::int32_t>(target);
        }
        std::vector<std::set<std::size_t>> neighbors(targets.size());
        std::vector<bool> movableTarget(targets.size(),false);
        for(std::size_t i=0;i<targets.size();++i){
            const auto vertex=static_cast<std::size_t>(targets[i].vertex);
            movableTarget[i]=targets[i].faceIndex<cadFaces.size() &&
                BRepAdaptor_Surface(cadFaces[targets[i].faceIndex]).GetType()!=GeomAbs_Plane &&
                (vertex>=supportedFaces.size()||supportedFaces[vertex].size()<=1U);
        }
        for (const auto& face : mesh.faces) {
            if (face.neighbour >= 0 || face.patch != patchIndex) continue;
            for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
                const auto a = face.vertices[edge];
                const auto b = face.vertices[(edge + 1U) % face.vertices.size()];
                if (a < 0 || b < 0 || static_cast<std::size_t>(a) >= targetForVertex.size() ||
                    static_cast<std::size_t>(b) >= targetForVertex.size()) continue;
                const auto ta = targetForVertex[static_cast<std::size_t>(a)];
                const auto tb = targetForVertex[static_cast<std::size_t>(b)];
                if (ta < 0 || tb < 0 || ta == tb) continue;
                const auto ia = static_cast<std::size_t>(ta), ib = static_cast<std::size_t>(tb);
                if (targets[ia].faceIndex != targets[ib].faceIndex) continue;
                neighbors[ia].insert(ib);
                neighbors[ib].insert(ia);
            }
        }
        std::vector<std::array<double, 3>> current(targets.size());
        for (std::size_t i = 0; i < targets.size(); ++i) current[i] = targets[i].target;
        constexpr std::size_t smoothingPasses = 4U;
        constexpr double relaxation = 0.25;
        constexpr double maximumStepPerPass = 0.05;
        for (std::size_t pass = 0; pass < smoothingPasses; ++pass) {
            auto next = current;
            double passMaximumMove = 0.0;
            for (std::size_t i = 0; i < current.size(); ++i) {
                if (!movableTarget[i] || neighbors[i].size() < 3U) continue;
                std::array<double, 3> average{};
                double localEdgeLength = 0.0;
                for (const auto neighbor : neighbors[i]) {
                    for (std::size_t axis = 0; axis < 3; ++axis) average[axis] += current[neighbor][axis];
                    localEdgeLength += std::sqrt(squaredDistance(current[i], current[neighbor]));
                }
                for (auto& coordinate : average) coordinate /= static_cast<double>(neighbors[i].size());
                localEdgeLength /= static_cast<double>(neighbors[i].size());
                const auto smoothProjection = projectToCadFace(average, targets[i].faceIndex);
                std::array<double, 3> proposal{};
                for (std::size_t axis = 0; axis < 3; ++axis)
                    proposal[axis] = current[i][axis] + relaxation * (smoothProjection.point_m[axis] - current[i][axis]);
                const double stepLength = std::sqrt(squaredDistance(proposal, current[i]));
                const double maximumStep = maximumStepPerPass * localEdgeLength;
                if (stepLength > maximumStep && stepLength > 0.0) {
                    const double scale = maximumStep / stepLength;
                    for (std::size_t axis = 0; axis < 3; ++axis)
                        proposal[axis] = current[i][axis] + scale * (proposal[axis] - current[i][axis]);
                }
                const auto constrained = projectToCadFace(proposal, targets[i].faceIndex);
                next[i] = constrained.point_m;
                passMaximumMove = std::max(passMaximumMove, std::sqrt(squaredDistance(next[i], current[i])));
            }
            if (passMaximumMove <= std::max(1.0e-13, mesh.cad_surface_target_size_m * 1.0e-10)) break;
            current = std::move(next);
            ++mesh.cad_tangential_smoothing_iterations;
            mesh.cad_tangential_smoothing_max_move_m = std::max(
                mesh.cad_tangential_smoothing_max_move_m, passMaximumMove);
        }
        for (std::size_t i = 0; i < targets.size(); ++i) targets[i].target = current[i];
    }

    struct DeformationPoint {
        std::int32_t vertex{};
        std::array<double, 3> original{};
        std::array<double, 3> target{};
    };
    std::vector<DeformationPoint> propagatedTargets;
    std::vector<bool> dirichlet(mesh.points_m.size(), false);
    std::vector<std::array<double, 3>> displacement(mesh.points_m.size());
    for (const auto& target : targets) {
        const auto vertex = static_cast<std::size_t>(target.vertex);
        dirichlet[vertex] = true;
        for (std::size_t axis = 0; axis < 3; ++axis) displacement[vertex][axis] = target.target[axis] - target.original[axis];
    }
    std::vector<bool> externalBoundaryVertex(mesh.points_m.size(), false);
    for (const auto& face : mesh.faces) if (face.neighbour < 0 && face.patch != patchIndex) {
        for (const auto vertex : face.vertices) externalBoundaryVertex[static_cast<std::size_t>(vertex)] = true;
    }
    for (std::size_t vertex = 0; vertex < externalBoundaryVertex.size(); ++vertex) if (externalBoundaryVertex[vertex]) {
        if (dirichlet[vertex] && std::hypot(displacement[vertex][0], displacement[vertex][1], displacement[vertex][2]) >
            featureMergeTolerance) {
            throw std::runtime_error("CAD boundary vertex also belongs to a fixed outer-domain patch");
        }
        dirichlet[vertex] = true;
        displacement[vertex] = {};
    }
    std::vector<std::vector<std::pair<std::int32_t, double>>> pointAdjacency(mesh.points_m.size());
    std::set<EdgeKey> meshEdges;
    for (const auto& face : mesh.faces) for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
        const auto a = face.vertices[edge];
        const auto b = face.vertices[(edge + 1U) % face.vertices.size()];
        const EdgeKey key = edgeKey(a, b);
        if (!meshEdges.insert(key).second) continue;
        const auto& pa = mesh.points_m[static_cast<std::size_t>(a)];
        const auto& pb = mesh.points_m[static_cast<std::size_t>(b)];
        const double lengthSquared = squaredDistance(pa, pb);
        if (!std::isfinite(lengthSquared) || lengthSquared <= featureMergeTolerance * featureMergeTolerance) continue;
        const double weight = 1.0 / lengthSquared;
        pointAdjacency[static_cast<std::size_t>(a)].emplace_back(b, weight);
        pointAdjacency[static_cast<std::size_t>(b)].emplace_back(a, weight);
    }
    std::vector<std::array<double, 3>> nextDisplacement = displacement;
    const double deformationTolerance = std::max(1.0e-12, mesh.cad_surface_target_size_m * 1.0e-8);
    constexpr std::size_t maxDeformationIterations = 300U;
    for (std::size_t iteration = 0; iteration < maxDeformationIterations; ++iteration) {
        double maximumChange = 0.0;
        for (std::size_t vertex = 0; vertex < displacement.size(); ++vertex) {
            if (dirichlet[vertex] || pointAdjacency[vertex].empty()) continue;
            std::array<double, 3> average{};
            double totalWeight = 0.0;
            for (const auto& [neighbour, weight] : pointAdjacency[vertex]) {
                totalWeight += weight;
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    average[axis] += weight * displacement[static_cast<std::size_t>(neighbour)][axis];
                }
            }
            if (!(totalWeight > 0.0)) continue;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                average[axis] /= totalWeight;
                const double updated = displacement[vertex][axis] + 0.8 * (average[axis] - displacement[vertex][axis]);
                maximumChange = std::max(maximumChange, std::abs(updated - displacement[vertex][axis]));
                nextDisplacement[vertex][axis] = updated;
            }
        }
        displacement.swap(nextDisplacement);
        mesh.cad_deformation_iterations = iteration + 1U;
        if (maximumChange <= deformationTolerance) break;
    }
    for (std::size_t vertex = 0; vertex < displacement.size(); ++vertex) {
        if (dirichlet[vertex]) continue;
        const double magnitude = std::hypot(displacement[vertex][0], displacement[vertex][1], displacement[vertex][2]);
        if (magnitude <= deformationTolerance) continue;
        if (vertex >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("CAD core deformation exceeds the supported point index range");
        }
        auto target = mesh.points_m[vertex];
        for (std::size_t axis = 0; axis < 3; ++axis) target[axis] += displacement[vertex][axis];
        propagatedTargets.push_back({static_cast<std::int32_t>(vertex), mesh.points_m[vertex], target});
        ++mesh.cad_deformed_vertices;
        mesh.cad_max_deformation_m = std::max(mesh.cad_max_deformation_m, magnitude);
    }

    auto cellVolumes = [&]() {
        std::vector<double> volumes(mesh.cells.size(), 0.0);
        for (const auto& face : mesh.faces) {
            if (face.vertices.size() < 3U) throw std::runtime_error("CAD snapping encountered a face with fewer than three vertices");
            if (face.owner < 0 || static_cast<std::size_t>(face.owner) >= volumes.size()) throw std::runtime_error("CAD snapping encountered an invalid face owner");
            const auto& p0 = mesh.points_m.at(static_cast<std::size_t>(face.vertices[0]));
            long double contribution = 0.0L;
            for (std::size_t i = 1; i + 1U < face.vertices.size(); ++i) {
                const auto& a = mesh.points_m.at(static_cast<std::size_t>(face.vertices[i]));
                const auto& b = mesh.points_m.at(static_cast<std::size_t>(face.vertices[i + 1U]));
                const long double cx = static_cast<long double>(a[1]) * b[2] - static_cast<long double>(a[2]) * b[1];
                const long double cy = static_cast<long double>(a[2]) * b[0] - static_cast<long double>(a[0]) * b[2];
                const long double cz = static_cast<long double>(a[0]) * b[1] - static_cast<long double>(a[1]) * b[0];
                contribution += (static_cast<long double>(p0[0]) * cx + static_cast<long double>(p0[1]) * cy +
                    static_cast<long double>(p0[2]) * cz) / 6.0L;
            }
            volumes[static_cast<std::size_t>(face.owner)] += static_cast<double>(contribution);
            if (face.neighbour >= 0) {
                if (static_cast<std::size_t>(face.neighbour) >= volumes.size()) throw std::runtime_error("CAD snapping encountered an invalid face neighbour");
                volumes[static_cast<std::size_t>(face.neighbour)] -= static_cast<double>(contribution);
            }
        }
        return volumes;
    };

    struct DeviationCheck {
        double maximum{};
        std::size_t samples{};
        std::size_t meshFace{};
        std::size_t cadFace{};
        std::array<double, 3> point{};
    };
    auto maximumDeviation = [&]() {
        DeviationCheck check;
        for (std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
            const auto& face = mesh.faces[faceIndex];
            if (face.neighbour >= 0 || face.patch != patchIndex) continue;
            if (face.vertices.size() < 3U) continue;
            std::array<double, 3> center{};
            for (const auto vertex : face.vertices) for (std::size_t axis = 0; axis < 3; ++axis) {
                center[axis] += mesh.points_m[static_cast<std::size_t>(vertex)][axis] / static_cast<double>(face.vertices.size());
            }
            const auto centerProjection = projectPointToFluidBoundary(cadModel, center);
            if (centerProjection.distance_m > check.maximum) {
                check.maximum = centerProjection.distance_m;
                check.meshFace = faceIndex;
                check.cadFace = centerProjection.face_index;
                check.point = center;
            }
            ++check.samples;
            if (strictQuality && check.maximum > maximumSurfaceDeviation) return check;
            for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
                const auto& a = mesh.points_m[static_cast<std::size_t>(face.vertices[edge])];
                const auto& b = mesh.points_m[static_cast<std::size_t>(face.vertices[(edge + 1U) % face.vertices.size()])];
                const std::array<double, 3> midpoint{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
                const auto midpointProjection = projectPointToFluidBoundary(cadModel, midpoint);
                if (midpointProjection.distance_m > check.maximum) {
                    check.maximum = midpointProjection.distance_m;
                    check.meshFace = faceIndex;
                    check.cadFace = midpointProjection.face_index;
                    check.point = midpoint;
                }
                ++check.samples;
                if (strictQuality && check.maximum > maximumSurfaceDeviation) return check;
            }
        }
        return check;
    };

    using Vector3 = std::array<double, 3>;
    std::vector<std::int32_t> targetForVertex(mesh.points_m.size(), -1);
    std::vector<bool> movableQualityTarget(targets.size(),false);
    for (std::size_t target = 0; target < targets.size(); ++target) {
        const auto vertex = static_cast<std::size_t>(targets[target].vertex);
        if (vertex < targetForVertex.size()) targetForVertex[vertex] = static_cast<std::int32_t>(target);
        movableQualityTarget[target]=targets[target].faceIndex<cadFaces.size() &&
            BRepAdaptor_Surface(cadFaces[targets[target].faceIndex]).GetType()!=GeomAbs_Plane &&
            (vertex>=supportedFaces.size()||supportedFaces[vertex].size()<=1U);
    }
    const auto movableQualityTargetCount=static_cast<std::size_t>(std::count(
        movableQualityTarget.begin(),movableQualityTarget.end(),true));
    std::vector<Vector3> qualityOffsets(mesh.points_m.size());
    auto crossProduct = [](const Vector3& a, const Vector3& b) -> Vector3 {
        return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
    };
    auto vectorDot = [](const Vector3& a, const Vector3& b) {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    };
    auto vectorSubtract = [](const Vector3& a,const Vector3& b)->Vector3 {
        return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
    };
    auto jacobianPenalty = [&](std::vector<Vector3>* gradient) {
        using BadCorner=std::pair<double,std::pair<std::size_t,std::size_t>>;
        std::vector<BadCorner> badCorners;
        for (std::size_t cell = 0; cell < mesh.cells.size() && cell < mesh.strict_hex.size(); ++cell) {
            if (!mesh.strict_hex[cell]) continue;
            const auto& ids = mesh.cells[cell];
            for (std::size_t corner = 0; corner < 8U; ++corner) {
                const auto p = mesh.points_m[static_cast<std::size_t>(ids[corner])];
                Vector3 u=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^1U])],p);
                Vector3 v=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^3U])],p);
                Vector3 w=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^4U])],p);
                const double su=((corner&1U)^((corner>>1U)&1U))?-1.0:1.0;
                const double sv=(corner&2U)?-1.0:1.0;
                const double sw=(corner&4U)?-1.0:1.0;
                for(std::size_t axis=0;axis<3;++axis){u[axis]*=su;v[axis]*=sv;w[axis]*=sw;}
                const double lu2=vectorDot(u,u),lv2=vectorDot(v,v),lw2=vectorDot(w,w);
                const double denominator=std::sqrt(lu2*lv2*lw2);
                const double value=denominator>0?vectorDot(u,crossProduct(v,w))/denominator:-1.0;
                const double deficit=0.1-value;
                if(!(deficit>0.0))continue;
                badCorners.push_back({deficit,{cell,corner}});
            }
        }
        if(badCorners.empty()){
            if(gradient)gradient->assign(targets.size(),Vector3{});
            return 0.0;
        }
        constexpr double temperature=0.005;
        double maximumDeficit=0.0;
        for(const auto& item:badCorners)maximumDeficit=std::max(maximumDeficit,item.first);
        double partition=0.0;
        for(const auto& item:badCorners)partition+=std::exp((item.first-maximumDeficit)/temperature);
        const double energy=maximumDeficit+temperature*std::log(partition);
        if(!gradient)return energy;
        gradient->assign(mesh.points_m.size(),Vector3{});
        for(const auto& item:badCorners){
            const auto cell=item.second.first,corner=item.second.second;
            const auto& ids=mesh.cells[cell];
            const auto p=mesh.points_m[static_cast<std::size_t>(ids[corner])];
            Vector3 u=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^1U])],p);
            Vector3 v=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^3U])],p);
            Vector3 w=vectorSubtract(mesh.points_m[static_cast<std::size_t>(ids[corner^4U])],p);
            const double su=((corner&1U)^((corner>>1U)&1U))?-1.0:1.0;
            const double sv=(corner&2U)?-1.0:1.0;
            const double sw=(corner&4U)?-1.0:1.0;
            for(std::size_t axis=0;axis<3;++axis){u[axis]*=su;v[axis]*=sv;w[axis]*=sw;}
            const double lu2=vectorDot(u,u),lv2=vectorDot(v,v),lw2=vectorDot(w,w);
            const double denominator=std::sqrt(lu2*lv2*lw2);
            if(!(denominator>0.0)||lu2<=0||lv2<=0||lw2<=0)continue;
            const double value=vectorDot(u,crossProduct(v,w))/denominator;
            const auto uv=crossProduct(v,w),vw=crossProduct(w,u),wu=crossProduct(u,v);
            Vector3 gu{},gv{},gw{};
            for(std::size_t axis=0;axis<3;++axis){
                gu[axis]=uv[axis]/denominator-value*u[axis]/lu2;
                gv[axis]=vw[axis]/denominator-value*v[axis]/lv2;
                gw[axis]=wu[axis]/denominator-value*w[axis]/lw2;
            }
            const double weight=std::exp((item.first-maximumDeficit)/temperature)/partition;
            const std::array<std::size_t,4> localVertices{corner,corner^1U,corner^3U,corner^4U};
            for(std::size_t axis=0;axis<3;++axis){
                const std::array<double,4> derivatives{
                    -(su*gu[axis]+sv*gv[axis]+sw*gw[axis]),su*gu[axis],sv*gv[axis],sw*gw[axis]};
                for(std::size_t n=0;n<localVertices.size();++n){
                    const auto vertex=static_cast<std::size_t>(ids[localVertices[n]]);
                    if(vertex>=targetForVertex.size())continue;
                    (*gradient)[vertex][axis]+=weight*derivatives[n];
                }
            }
        }
        return energy;
    };
    auto optimizeSurfaceJacobian = [&](double boundaryFraction,DeviationCheck& acceptedOptimizationDeviation) {
        if (planarDomain || !(boundaryFraction>0.0)) return false;
        std::vector<Vector3> gradient;
        const double before=jacobianPenalty(&gradient);
        if(!(before>1.0e-16))return false;
        ++mesh.cad_quality_optimization_attempts;
        if(mesh.cad_quality_optimization_attempts==1U)
            mesh.cad_quality_optimization_objective_before=before;
        mesh.cad_quality_optimization_objective_after=before;
        const auto beforeQuality=auditMeshQuality(mesh);
        const auto oldTargets=targets;
        const auto oldPointPositions=mesh.points_m;
        const auto oldQualityOffsets=qualityOffsets;
        constexpr std::array<double,5> lineSearch{{1.0,0.5,0.25,0.125,0.0625}};
        const double nominalStep=0.025*mesh.cad_surface_target_size_m;
        for(const double lineScale:lineSearch){
            targets=oldTargets;
            mesh.points_m=oldPointPositions;
            qualityOffsets=oldQualityOffsets;
            auto trialOffsets=oldQualityOffsets;
            bool moved=false;
            for(std::size_t vertex=0;vertex<gradient.size();++vertex){
                const double norm=std::sqrt(vectorDot(gradient[vertex],gradient[vertex]));
                if(!(norm>1.0e-14)||!std::isfinite(norm))continue;
                const auto targetIndex=vertex<targetForVertex.size()?targetForVertex[vertex]:-1;
                if(targetIndex>=0){
                    const auto i=static_cast<std::size_t>(targetIndex);
                    if(!movableQualityTarget[i])continue;
                    Vector3 proposal=targets[i].target;
                    for(std::size_t axis=0;axis<3;++axis)
                        proposal[axis]+=lineScale*nominalStep/boundaryFraction*gradient[vertex][axis]/norm;
                    const auto projection=projectToCadFace(proposal,targets[i].faceIndex);
                    const double actualMove=boundaryFraction*std::sqrt(squaredDistance(projection.point_m,targets[i].target));
                    if(!(actualMove>mesh.cad_surface_target_size_m*1.0e-10))continue;
                    ++mesh.cad_quality_optimization_proposed_vertices;
                    targets[i].target=projection.point_m;
                    auto& point=mesh.points_m[vertex];
                    for(std::size_t axis=0;axis<3;++axis)
                        point[axis]=targets[i].original[axis]+boundaryFraction*(targets[i].target[axis]-targets[i].original[axis]);
                }else{
                    if(vertex<externalBoundaryVertex.size()&&externalBoundaryVertex[vertex])continue;
                    ++mesh.cad_quality_optimization_proposed_vertices;
                    for(std::size_t axis=0;axis<3;++axis){
                        const double delta=lineScale*nominalStep*gradient[vertex][axis]/norm;
                        mesh.points_m[vertex][axis]+=delta;
                        trialOffsets[vertex][axis]+=delta;
                    }
                }
                moved=true;
            }
            if(!moved)break;
            const double after=jacobianPenalty(nullptr);
            mesh.cad_quality_optimization_objective_after=after;
            const auto candidateQuality=auditMeshQuality(mesh);
            const bool noFaceQualityRegression=
                candidateQuality.maximum_non_orthogonality_degrees<=beforeQuality.maximum_non_orthogonality_degrees+0.25 &&
                candidateQuality.p99_non_orthogonality_degrees<=beforeQuality.p99_non_orthogonality_degrees+0.05 &&
                candidateQuality.maximum_internal_skewness<=beforeQuality.maximum_internal_skewness+0.005 &&
                candidateQuality.p99_internal_skewness<=beforeQuality.p99_internal_skewness+0.005;
            const bool jacobianNotRegressed=
                candidateQuality.minimum_regular_hex_scaled_jacobian>=beforeQuality.minimum_regular_hex_scaled_jacobian-1.0e-8;
            const bool strictTopologyValid=candidateQuality.negative_or_zero_volume_cells==0&&candidateQuality.degenerate_faces==0;
            const bool objectiveDecreased=after<before-std::max(1.0e-14,before*1.0e-8);
            const bool objectiveProgress=
                candidateQuality.minimum_regular_hex_scaled_jacobian>beforeQuality.minimum_regular_hex_scaled_jacobian+1.0e-8 ||
                candidateQuality.p99_non_orthogonality_degrees<beforeQuality.p99_non_orthogonality_degrees-0.01 ||
                candidateQuality.p99_internal_skewness<beforeQuality.p99_internal_skewness-1.0e-5;
            if(!strictTopologyValid||!noFaceQualityRegression||!jacobianNotRegressed||!objectiveDecreased||!objectiveProgress)continue;
            const auto deviation=maximumDeviation();
            if(deviation.maximum>maximumSurfaceDeviation)continue;
            qualityOffsets=std::move(trialOffsets);
            acceptedOptimizationDeviation=deviation;
            ++mesh.cad_quality_optimization_iterations;
            if(mesh.cad_quality_optimization_iterations==1U)mesh.cad_quality_optimization_objective_before=before;
            return true;
        }
        targets=oldTargets;
        mesh.points_m=oldPointPositions;
        qualityOffsets=oldQualityOffsets;
        return false;
    };

    const auto originalVolumes = cellVolumes();
    const auto unconsolidatedFaces=mesh.faces;
    const auto originalStrictHex=mesh.strict_hex;
    std::vector<double> acceptedVolumes;
    double acceptedFraction = 0.0;
    double acceptedCoreDeformationFraction = 0.0;
    double acceptedDeviation = std::numeric_limits<double>::infinity();
    std::size_t deviationSamples = 0;
    DeviationCheck latestDeviation;
    std::map<double,DeviationCheck> deviationCache;
    std::size_t firstInvalidVolumeCell = std::numeric_limits<std::size_t>::max();
    double firstInvalidVolume = std::numeric_limits<double>::quiet_NaN();
    std::size_t qualityRollbackCount=0;
    std::string latestQualityFailure;
    double lastGeometryPassingFraction=0.0;
    bool geometryFractionSearchComplete=false;
    bool jacobianOptimizationStopped=false;
    // An explicitly enlarged geometric tolerance must be allowed to trade
    // more normal displacement for quality; stopping at 88% projection would
    // otherwise leave much of the authorized geometric budget unexplored.
    const int maximumSnapAttempts = strictQuality && mesh.cad_surface_deviation_overridden ? 51 : 13;
    for (int attempt = 0; attempt < maximumSnapAttempts; ++attempt) {
        if (attempt >= 13 && geometryFractionSearchComplete) break;
        mesh.faces=unconsolidatedFaces;
        mesh.strict_hex=originalStrictHex;
        mesh.cad_consolidated_faces=0;
        double fraction = strictQuality ? std::max(0.0,1.0-0.01*attempt) : std::ldexp(1.0, -attempt);
        if(strictQuality && geometryFractionSearchComplete)fraction=lastGeometryPassingFraction;
        const double coreFraction = strictQuality ? std::ldexp(1.0,-std::min(attempt,12)) : fraction;
        for (const auto& target : propagatedTargets) {
            auto& point = mesh.points_m[static_cast<std::size_t>(target.vertex)];
            for (std::size_t axis = 0; axis < 3; ++axis) point[axis] = target.original[axis] + coreFraction * (target.target[axis] - target.original[axis]);
        }
        for (const auto& target : targets) {
            auto& point = mesh.points_m[static_cast<std::size_t>(target.vertex)];
            const double snapFraction = fraction;
            for (std::size_t axis = 0; axis < 3; ++axis) point[axis] = target.original[axis] + snapFraction * (target.target[axis] - target.original[axis]);
        }
        for(std::size_t vertex=0;vertex<qualityOffsets.size();++vertex){
            if(targetForVertex[vertex]>=0)continue;
            for(std::size_t axis=0;axis<3;++axis)mesh.points_m[vertex][axis]+=qualityOffsets[vertex][axis];
        }
            auto volumes = cellVolumes();
        bool volumeValid = volumes.size() == mesh.cells.size();
        for (std::size_t cell = 0; cell < volumes.size() && volumeValid; ++cell) {
            const double reference = cell < mesh.cell_volumes_m3.size() ? mesh.cell_volumes_m3[cell] : originalVolumes[cell];
            const double lower = std::max(1.0e-30, std::abs(reference) * 1.0e-12);
            volumeValid = std::isfinite(volumes[cell]) && volumes[cell] > lower;
            if (!volumeValid && firstInvalidVolumeCell == std::numeric_limits<std::size_t>::max()) {
                firstInvalidVolumeCell = cell;
                firstInvalidVolume = volumes[cell];
            }
        }
        if (!volumeValid) continue;
        if(strictQuality) {
            mesh.cad_consolidated_faces=consolidateCoplanarFaces(mesh,
                std::max(1.0e-11*domainDiagonal,featureMergeTolerance),0.05);
            volumes=cellVolumes();
            volumeValid=volumes.size()==mesh.cells.size();
            for(std::size_t cell=0;cell<volumes.size() && volumeValid;++cell){
                const double reference=cell<mesh.cell_volumes_m3.size()?mesh.cell_volumes_m3[cell]:originalVolumes[cell];
                const double lower=std::max(1.0e-30,std::abs(reference)*1.0e-12);
                volumeValid=std::isfinite(volumes[cell]) && volumes[cell]>lower;
            }
            if(!volumeValid)continue;
            auto deviationResult=deviationCache.find(fraction);
            if(deviationResult==deviationCache.end())
                deviationResult=deviationCache.emplace(fraction,maximumDeviation()).first;
            auto deviation=deviationResult->second;
            latestDeviation=deviation;
            acceptedDeviation=deviation.maximum;
            deviationSamples=deviation.samples;
            if(deviation.maximum>maximumSurfaceDeviation){
                // Once a geometry-safe boundary fraction has been found,
                // keep it fixed and spend remaining trials on the independent
                // core deformation fraction. Otherwise shrinking both at
                // once can make geometry fail before core-quality recovery is
                // explored.
                if(lastGeometryPassingFraction>0.0)geometryFractionSearchComplete=true;
                continue;
            }
            lastGeometryPassingFraction=fraction;
            if(mesh.cad_quality_optimization_iterations<3U&&!jacobianOptimizationStopped){
                DeviationCheck optimizedDeviation=deviation;
                if(optimizeSurfaceJacobian(fraction,optimizedDeviation)){
                    deviationCache.clear();
                    deviation=optimizedDeviation;
                    latestDeviation=deviation;
                    acceptedDeviation=deviation.maximum;
                    deviationSamples=deviation.samples;
                }else jacobianOptimizationStopped=true;
            }
            if(deviation.maximum>maximumSurfaceDeviation){
                geometryFractionSearchComplete=true;
                continue;
            }
            HexMesh reconstructed;
            HexMesh* qualityMesh=&mesh;
            const auto baseQuality=auditMeshQuality(mesh);
            const bool baseQualityValid=baseQuality.negative_or_zero_volume_cells==0 && baseQuality.degenerate_faces==0 &&
                baseQuality.minimum_regular_hex_scaled_jacobian>=0.1 &&
                baseQuality.maximum_non_orthogonality_degrees<=75.0 && baseQuality.p99_non_orthogonality_degrees<=60.0 &&
                baseQuality.maximum_internal_skewness<=0.5 && baseQuality.p99_internal_skewness<=0.25;
            std::size_t reconstructedParents=0;
            if(!baseQualityValid) {
                reconstructed=mesh;
                reconstructedParents=splitLowQualityHexesIntoPyramids(reconstructed,0.1,maximumCells,75.0);
                if(reconstructedParents>0U)qualityMesh=&reconstructed;
            }
            const auto quality=auditMeshQuality(*qualityMesh);
            const bool qualityValid=quality.negative_or_zero_volume_cells==0 && quality.degenerate_faces==0 &&
                quality.minimum_regular_hex_scaled_jacobian>=0.1 &&
                quality.maximum_non_orthogonality_degrees<=75.0 && quality.p99_non_orthogonality_degrees<=60.0 &&
                quality.maximum_internal_skewness<=0.5 && quality.p99_internal_skewness<=0.25;
            if(qualityValid && reconstructedParents>0U) {
                reconstructed.cell_volumes_m3=quality.cell_volumes;
                mesh=std::move(reconstructed);
                volumes=quality.cell_volumes;
            }
            if(!qualityValid){
                ++qualityRollbackCount;
                latestQualityFailure="min_scaled_jacobian="+std::to_string(quality.minimum_regular_hex_scaled_jacobian)+
                    ", reconstructed_parent_cells="+std::to_string(reconstructedParents)+
                    ", pre_reconstruction_min_scaled_jacobian="+std::to_string(baseQuality.minimum_regular_hex_scaled_jacobian)+
                    ", pre_reconstruction_nonorth="+std::to_string(baseQuality.maximum_non_orthogonality_degrees)+"/"+
                    std::to_string(baseQuality.p99_non_orthogonality_degrees)+
                    ", pre_reconstruction_skew="+std::to_string(baseQuality.maximum_internal_skewness)+"/"+
                    std::to_string(baseQuality.p99_internal_skewness)+
                    ", surface_quality_optimization_iterations="+
                    std::to_string(mesh.cad_quality_optimization_iterations)+
                    ", surface_quality_optimization_attempts="+
                    std::to_string(mesh.cad_quality_optimization_attempts)+
                    ", CAD_faces="+std::to_string(cadFaces.size())+
                    ", movable_CAD_targets="+std::to_string(movableQualityTargetCount)+
                    ", surface_quality_optimization_proposed_vertices="+
                    std::to_string(mesh.cad_quality_optimization_proposed_vertices)+
                    ", surface_quality_objective="+
                    std::to_string(mesh.cad_quality_optimization_objective_before)+"->"+
                    std::to_string(mesh.cad_quality_optimization_objective_after)+
                    ", boundary_fraction="+std::to_string(fraction)+", core_fraction="+std::to_string(coreFraction)+
                    ", max/p99 nonorth="+std::to_string(quality.maximum_non_orthogonality_degrees)+"/"+
                    std::to_string(quality.p99_non_orthogonality_degrees)+", max/p99 skew="+
                    std::to_string(quality.maximum_internal_skewness)+"/"+
                    std::to_string(quality.p99_internal_skewness);
                if(quality.maximum_non_orthogonality_face_index<qualityMesh->faces.size()) {
                    const auto& worstFace=qualityMesh->faces[quality.maximum_non_orthogonality_face_index];
                    latestQualityFailure+=", worst_nonorth_face="+
                        std::to_string(quality.maximum_non_orthogonality_face_index)+
                        "(owner="+std::to_string(worstFace.owner)+",neighbour="+
                        std::to_string(worstFace.neighbour)+",patch="+std::to_string(worstFace.patch)+",ring=";
                    for(const auto vertex:worstFace.vertices) {
                        const auto& point=qualityMesh->points_m[static_cast<std::size_t>(vertex)];
                        latestQualityFailure+=std::to_string(vertex)+":"+
                            std::to_string(point[0])+","+std::to_string(point[1])+","+
                            std::to_string(point[2])+";";
                    }
                    latestQualityFailure+=")";
                }
                continue;
            }
        }
        const auto deviation = strictQuality ? latestDeviation : maximumDeviation();
        if(!strictQuality)latestDeviation=deviation;
        if (strictQuality && deviation.maximum > maximumSurfaceDeviation) {
            acceptedDeviation = deviation.maximum;
            deviationSamples = deviation.samples;
            continue;
        }
        acceptedFraction = fraction;
        acceptedCoreDeformationFraction = coreFraction;
        acceptedVolumes = std::move(volumes);
        acceptedDeviation = deviation.maximum;
        deviationSamples = deviation.samples;
        break;
    }
    if (acceptedFraction == 0.0) {
        std::string failureDetails;
        if (firstInvalidVolumeCell != std::numeric_limits<std::size_t>::max()) {
            failureDetails = ", nonpositive cell volume at cell " + std::to_string(firstInvalidVolumeCell) +
                " (signed volume=" + std::to_string(firstInvalidVolume) + " m^3)";
            std::set<std::int32_t> diagnosticCellVertices;
            for (const auto& face : mesh.faces) if (face.owner == static_cast<std::int32_t>(firstInvalidVolumeCell) ||
                face.neighbour == static_cast<std::int32_t>(firstInvalidVolumeCell)) {
                diagnosticCellVertices.insert(face.vertices.begin(), face.vertices.end());
            }
            std::array<double, 3> diagnosticCenter{};
            for (const auto vertex : diagnosticCellVertices) for (std::size_t axis = 0; axis < 3; ++axis) {
                diagnosticCenter[axis] += mesh.points_m[static_cast<std::size_t>(vertex)][axis] /
                    static_cast<double>(diagnosticCellVertices.size());
            }
            std::array<double, 3> closureArea{};
            long double localCellVolume = 0.0L;
            for (std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
                const auto& face = mesh.faces[faceIndex];
                if (face.owner != static_cast<std::int32_t>(firstInvalidVolumeCell) &&
                    face.neighbour != static_cast<std::int32_t>(firstInvalidVolumeCell)) continue;
                failureDetails += " {face=" + std::to_string(faceIndex) + ",patch=" + std::to_string(face.patch) + ",ring=";
                double faceContribution = 0.0;
                if (face.vertices.size() >= 3U) {
                    const auto& p0 = mesh.points_m[static_cast<std::size_t>(face.vertices[0])];
                    long double contribution = 0.0L;
                    for (std::size_t vertex = 1; vertex + 1U < face.vertices.size(); ++vertex) {
                        const auto& a = mesh.points_m[static_cast<std::size_t>(face.vertices[vertex])];
                        const auto& b = mesh.points_m[static_cast<std::size_t>(face.vertices[vertex + 1U])];
                        const long double cx = static_cast<long double>(a[1]) * b[2] - static_cast<long double>(a[2]) * b[1];
                        const long double cy = static_cast<long double>(a[2]) * b[0] - static_cast<long double>(a[0]) * b[2];
                        const long double cz = static_cast<long double>(a[0]) * b[1] - static_cast<long double>(a[1]) * b[0];
                        contribution += (static_cast<long double>(p0[0]) * cx + static_cast<long double>(p0[1]) * cy +
                            static_cast<long double>(p0[2]) * cz) / 6.0L;
                        const std::array<double, 3> da{a[0] - p0[0], a[1] - p0[1], a[2] - p0[2]};
                        const std::array<double, 3> db{b[0] - p0[0], b[1] - p0[1], b[2] - p0[2]};
                        const std::array<double, 3> cross{da[1] * db[2] - da[2] * db[1],
                            da[2] * db[0] - da[0] * db[2], da[0] * db[1] - da[1] * db[0]};
                        const double sign = face.owner == static_cast<std::int32_t>(firstInvalidVolumeCell) ? 1.0 : -1.0;
                        for (std::size_t axis = 0; axis < 3; ++axis) closureArea[axis] += sign * 0.5 * cross[axis];
                        const std::array<double, 3> relative{p0[0] - diagnosticCenter[0],
                            p0[1] - diagnosticCenter[1], p0[2] - diagnosticCenter[2]};
                        localCellVolume += static_cast<long double>(sign) *
                            (relative[0] * cross[0] + relative[1] * cross[1] + relative[2] * cross[2]) / 6.0L;
                    }
                    faceContribution = static_cast<double>(contribution) *
                        (face.owner == static_cast<std::int32_t>(firstInvalidVolumeCell) ? 1.0 : -1.0);
                }
                failureDetails += "signed-face-volume=" + std::to_string(faceContribution) + ",ring=";
                for (const auto vertex : face.vertices) {
                    const auto& point = mesh.points_m[static_cast<std::size_t>(vertex)];
                    failureDetails += std::to_string(vertex) + ":(" + std::to_string(point[0]) + "," +
                        std::to_string(point[1]) + "," + std::to_string(point[2]) + ");";
                }
                failureDetails += "}";
            }
            failureDetails += ", closure-area=(" + std::to_string(closureArea[0]) + "," +
                std::to_string(closureArea[1]) + "," + std::to_string(closureArea[2]) +
                "), local-origin volume=" + std::to_string(static_cast<double>(localCellVolume)) + " m^3";
        }
        if(qualityRollbackCount)failureDetails += ", mesh-quality rollback count="+
            std::to_string(qualityRollbackCount)+", last rejected metrics="+latestQualityFailure;
        for (const auto& target : targets) mesh.points_m[static_cast<std::size_t>(target.vertex)] = target.original;
        std::string offendingGeometry;
        if (latestDeviation.meshFace < mesh.faces.size()) {
            const auto& face = mesh.faces[latestDeviation.meshFace];
            for (const auto vertex : face.vertices) {
                if (vertex < 0 || static_cast<std::size_t>(vertex) >= targetAtVertex.size()) continue;
                const auto targetIndex = targetAtVertex[static_cast<std::size_t>(vertex)];
                if (targetIndex < 0) continue;
                const auto& target = targets[static_cast<std::size_t>(targetIndex)];
                offendingGeometry += " [point=" + std::to_string(vertex) + ",CAD face=" +
                    std::to_string(target.faceIndex) + ",target=(" + std::to_string(target.target[0]) + "," +
                    std::to_string(target.target[1]) + "," + std::to_string(target.target[2]) + ")]";
            }
        }
        const bool geometryRejected=std::isfinite(acceptedDeviation) && acceptedDeviation>maximumSurfaceDeviation;
        throw std::runtime_error(std::string(geometryRejected ?
            "CAD surface snapping failed strict geometry checks" :
            "CAD surface snapping failed strict positive-volume/full-quality checks")+"; target vertices="+
            std::to_string(targets.size()) + ", first violating sample deviation=" + std::to_string(acceptedDeviation) +
            " m, allowed=" + std::to_string(maximumSurfaceDeviation) + " m, offending wall face=" +
            std::to_string(latestDeviation.meshFace) + ", CAD face=" + std::to_string(latestDeviation.cadFace) +
            ", sample=(" + std::to_string(latestDeviation.point[0]) + "," +
            std::to_string(latestDeviation.point[1]) + "," + std::to_string(latestDeviation.point[2]) + ")" +
            failureDetails +
            offendingGeometry);
    }
    mesh.cad_quality_rollback_count=qualityRollbackCount;

    // Snapping independent Cartesian boundary vertices onto the same CAD
    // feature can create different subdivisions of that one CAD edge.  Bring
    // every face edge on a shared B-Rep edge to the same global node chain
    // before dropping zero-area feature patches.  Without this pass, a
    // collapsed line patch may be the only thing connecting otherwise valid
    // neighboring cell-face rings.
    using FeatureNode = std::pair<double, std::int32_t>;
    std::vector<std::vector<std::pair<std::size_t, double>>> featureEdgesAtVertex(mesh.points_m.size());
    const double featureEdgeTolerance = std::max(featureMergeTolerance * 8.0, 1.0e-11 * domainDiagonal);
    for (auto& edge : cadFeatureEdges) {
        BRepBndLib::Add(edge.edge, edge.bounds, Standard_False);
        if (!edge.bounds.IsVoid()) edge.bounds.Enlarge(featureEdgeTolerance);
    }
    try {
        for (const auto& target : targets) {
            const auto vertex = static_cast<std::size_t>(target.vertex);
            const gp_Pnt targetPoint(target.target[0], target.target[1], target.target[2]);
            const TopoDS_Vertex query = BRepBuilderAPI_MakeVertex(targetPoint).Vertex();
            for (std::size_t edgeIndex = 0; edgeIndex < cadFeatureEdges.size(); ++edgeIndex) {
                const auto& featureEdge = cadFeatureEdges[edgeIndex];
                if (featureEdge.bounds.IsVoid() || featureEdge.bounds.IsOut(targetPoint)) continue;
                BRepExtrema_DistShapeShape distance(query, featureEdge.edge);
                distance.Perform();
                if (!distance.IsDone()) throw std::runtime_error("OCCT failed to associate a CAD feature-edge mesh node");
                if (distance.Value() > featureEdgeTolerance || distance.NbSolution() == 0) continue;
                Standard_Real parameter = 0.0;
                if (distance.SupportTypeShape2(1) == BRepExtrema_IsOnEdge) {
                    distance.ParOnEdgeS2(1, parameter);
                } else if (distance.SupportTypeShape2(1) == BRepExtrema_IsVertex) {
                    const TopoDS_Vertex edgeVertex = TopoDS::Vertex(distance.SupportOnShape2(1));
                    parameter = BRep_Tool::Parameter(edgeVertex, featureEdge.edge);
                } else {
                    continue;
                }
                if (!std::isfinite(parameter)) throw std::runtime_error("CAD feature-edge projection returned a non-finite curve parameter");
                featureEdgesAtVertex[vertex].emplace_back(edgeIndex, parameter);
            }
        }
    } catch (const Standard_Failure& error) {
        const char* message = error.GetMessageString();
        throw std::runtime_error(std::string("OCCT CAD feature-edge node association failed: ") +
            (message ? message : "unknown OCCT failure"));
    }
    const double featureParameterTolerance = 1.0e-12;
    std::vector<std::vector<FeatureNode>> nodesOnFeatureEdge(cadFeatureEdges.size());
    for (std::size_t vertex = 0; vertex < featureEdgesAtVertex.size(); ++vertex) {
        for (const auto& [edgeIndex, parameter] : featureEdgesAtVertex[vertex]) {
            nodesOnFeatureEdge[edgeIndex].emplace_back(parameter, static_cast<std::int32_t>(vertex));
        }
    }
    for (auto& nodes : nodesOnFeatureEdge) {
        std::sort(nodes.begin(), nodes.end(), [](const FeatureNode& a, const FeatureNode& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
    }
    for (auto& face : mesh.faces) {
        if (face.neighbour < 0 || face.vertices.size() < 2U) continue;
        const std::set<std::int32_t> existingFaceVertices(face.vertices.begin(), face.vertices.end());
        std::vector<std::int32_t> expanded;
        expanded.reserve(face.vertices.size() + 2U);
        for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
            const auto first = face.vertices[edge];
            const auto second = face.vertices[(edge + 1U) % face.vertices.size()];
            expanded.push_back(first);
            if (first < 0 || second < 0 || static_cast<std::size_t>(first) >= featureEdgesAtVertex.size() ||
                static_cast<std::size_t>(second) >= featureEdgesAtVertex.size()) {
                throw std::runtime_error("CAD feature-edge synchronization encountered an invalid face vertex");
            }
            std::set<std::size_t> commonFeatureEdges;
            for (const auto& [firstEdge, firstParameter] : featureEdgesAtVertex[static_cast<std::size_t>(first)]) {
                for (const auto& [secondEdge, secondParameter] : featureEdgesAtVertex[static_cast<std::size_t>(second)]) {
                    if (firstEdge != secondEdge) continue;
                    commonFeatureEdges.insert(firstEdge);
                }
            }
            std::set<std::int32_t> insertedNodes;
            for (const auto commonEdge : commonFeatureEdges) {
                const auto firstMembership = std::find_if(featureEdgesAtVertex[static_cast<std::size_t>(first)].begin(),
                    featureEdgesAtVertex[static_cast<std::size_t>(first)].end(), [&](const auto& item) { return item.first == commonEdge; });
                const auto secondMembership = std::find_if(featureEdgesAtVertex[static_cast<std::size_t>(second)].begin(),
                    featureEdgesAtVertex[static_cast<std::size_t>(second)].end(), [&](const auto& item) { return item.first == commonEdge; });
                const double lower = std::min(firstMembership->second, secondMembership->second);
                const double upper = std::max(firstMembership->second, secondMembership->second);
                if (upper - lower <= featureParameterTolerance) continue;
                const bool increasing = firstMembership->second < secondMembership->second;
                std::vector<FeatureNode> edgeInsertions;
                for (const auto& node : nodesOnFeatureEdge[commonEdge]) {
                    if (node.second != first && node.second != second && existingFaceVertices.count(node.second) == 0U &&
                        node.first > lower + featureParameterTolerance &&
                        node.first < upper - featureParameterTolerance && insertedNodes.insert(node.second).second) {
                        edgeInsertions.push_back(node);
                    }
                }
                std::sort(edgeInsertions.begin(), edgeInsertions.end(), [increasing](const FeatureNode& a, const FeatureNode& b) {
                    if (a.first != b.first) return increasing ? a.first < b.first : a.first > b.first;
                    return a.second < b.second;
                });
                for (const auto& node : edgeInsertions) if (expanded.back() != node.second) expanded.push_back(node.second);
            }
        }
        if (expanded.size() != face.vertices.size()) {
            face.vertices = std::move(expanded);
            if (face.owner >= 0 && static_cast<std::size_t>(face.owner) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.owner)] = false;
            }
            if (face.neighbour >= 0 && static_cast<std::size_t>(face.neighbour) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.neighbour)] = false;
            }
        }
    }

    // Canonicalize nodes that were already topologically coincident before
    // snapping, plus independently created feature cuts that land at the same
    // CAD vertex.  Distant Cartesian vertices can both project to one corner;
    // welding those would collapse a whole boundary face and hide a bad snap.
    std::vector<std::int32_t> pointRemap(mesh.points_m.size(), -1);
    std::map<std::tuple<std::int64_t, std::int64_t, std::int64_t>, std::vector<std::int32_t>> pointBuckets;
    const double mergeToleranceSquared = featureMergeTolerance * featureMergeTolerance;
    std::set<EdgeKey> topologicalEdges;
    for (const auto& face : mesh.faces) for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
        topologicalEdges.insert(edgeKey(face.vertices[edge], face.vertices[(edge + 1U) % face.vertices.size()]));
    }
    for (std::size_t point = 0; point < mesh.points_m.size(); ++point) {
        std::array<std::int64_t, 3> bucket{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const long double coordinate = std::floor(static_cast<long double>(mesh.points_m[point][axis]) / featureMergeTolerance);
            if (coordinate < static_cast<long double>(std::numeric_limits<std::int64_t>::min() + 2) ||
                coordinate > static_cast<long double>(std::numeric_limits<std::int64_t>::max() - 2)) {
                throw std::runtime_error("CAD feature point canonicalization exceeds the supported coordinate range");
            }
            bucket[axis] = static_cast<std::int64_t>(coordinate);
        }
        std::int32_t canonical = static_cast<std::int32_t>(point);
        for (std::int64_t dx = -1; dx <= 1 && canonical == static_cast<std::int32_t>(point); ++dx) {
            for (std::int64_t dy = -1; dy <= 1 && canonical == static_cast<std::int32_t>(point); ++dy) {
                for (std::int64_t dz = -1; dz <= 1 && canonical == static_cast<std::int32_t>(point); ++dz) {
                    const auto found = pointBuckets.find({bucket[0] + dx, bucket[1] + dy, bucket[2] + dz});
                    if (found == pointBuckets.end()) continue;
                    for (const auto candidate : found->second) {
                        bool shareFeatureEdge = false;
                        for (const auto& [pointEdge, parameter] : featureEdgesAtVertex[point]) {
                            (void)parameter;
                            if (std::any_of(featureEdgesAtVertex[static_cast<std::size_t>(candidate)].begin(),
                                featureEdgesAtVertex[static_cast<std::size_t>(candidate)].end(),
                                [&](const auto& item) { return item.first == pointEdge; })) {
                                shareFeatureEdge = true;
                                break;
                            }
                        }
                        if (squaredDistance(mesh.points_m[point], mesh.points_m[static_cast<std::size_t>(candidate)]) <= mergeToleranceSquared &&
                            (squaredDistance(pointOrigins[point], pointOrigins[static_cast<std::size_t>(candidate)]) <= mergeToleranceSquared ||
                             topologicalEdges.find(edgeKey(static_cast<std::int32_t>(point), candidate)) != topologicalEdges.end() ||
                             shareFeatureEdge)) {
                        canonical = pointRemap[static_cast<std::size_t>(candidate)];
                        break;
                    }
                    }
                }
            }
        }
        pointRemap[point] = canonical;
        if (canonical == static_cast<std::int32_t>(point)) {
            pointBuckets[{bucket[0], bucket[1], bucket[2]}].push_back(static_cast<std::int32_t>(point));
        }
    }
    // A snap can collapse the endpoints of a real mesh edge even when the
    // endpoints were far apart before projection.  Weld such an edge as one
    // global topological operation so every incident face sees the same ID.
    std::vector<std::int32_t> weldParent(mesh.points_m.size());
    std::iota(weldParent.begin(), weldParent.end(), 0);
    auto weldRoot = [&](std::int32_t vertex) {
        std::int32_t root = vertex;
        while (weldParent[static_cast<std::size_t>(root)] != root) {
            root = weldParent[static_cast<std::size_t>(root)];
        }
        while (weldParent[static_cast<std::size_t>(vertex)] != vertex) {
            const auto next = weldParent[static_cast<std::size_t>(vertex)];
            weldParent[static_cast<std::size_t>(vertex)] = root;
            vertex = next;
        }
        return root;
    };
    for (const auto& face : mesh.faces) for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
        const auto first = face.vertices[edge];
        const auto second = face.vertices[(edge + 1U) % face.vertices.size()];
        if (first < 0 || second < 0 || static_cast<std::size_t>(first) >= pointRemap.size() ||
            static_cast<std::size_t>(second) >= pointRemap.size()) {
            throw std::runtime_error("CAD feature edge welding encountered an invalid face vertex");
        }
        if (squaredDistance(mesh.points_m[static_cast<std::size_t>(first)],
                mesh.points_m[static_cast<std::size_t>(second)]) > mergeToleranceSquared) continue;
        const auto a = weldRoot(pointRemap[static_cast<std::size_t>(first)]);
        const auto b = weldRoot(pointRemap[static_cast<std::size_t>(second)]);
        if (a != b) weldParent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
    }
    for (auto& canonical : pointRemap) canonical = weldRoot(canonical);
    std::vector<MeshFace> canonicalFaces;
    canonicalFaces.reserve(mesh.faces.size());
    std::vector<MeshFace> collapsedFaces;
    auto compactRing = [](std::vector<std::int32_t>& ring, std::int32_t owner, std::int32_t neighbour,
                          const std::vector<std::int32_t>& sourceRing) {
        std::vector<std::int32_t> compact;
        compact.reserve(ring.size());
        for (const auto vertex : ring) {
            if (compact.empty() || vertex != compact.back()) compact.push_back(vertex);
        }
        if (compact.size() > 1U && compact.front() == compact.back()) compact.pop_back();
        std::set<std::int32_t> unique(compact.begin(), compact.end());
        if (unique.size() != compact.size()) {
            std::string details = "CAD feature point welding produced a non-simple face ring (owner " +
                std::to_string(owner) + ", neighbour " + std::to_string(neighbour) + ", vertices";
            for (const auto vertex : compact) details += " " + std::to_string(vertex);
            details += "; source";
            for (const auto vertex : sourceRing) details += " " + std::to_string(vertex);
            details += ")";
            throw std::runtime_error(details);
        }
        ring = std::move(compact);
    };
    for (auto face : mesh.faces) {
        bool changed = false;
        const auto sourceRing = face.vertices;
        for (auto& vertex : face.vertices) {
            if (vertex < 0 || static_cast<std::size_t>(vertex) >= pointRemap.size()) {
                throw std::runtime_error("CAD snapping canonicalization encountered an invalid face vertex");
            }
            const auto mapped = pointRemap[static_cast<std::size_t>(vertex)];
            changed = changed || mapped != vertex;
            vertex = mapped;
        }
        const std::size_t originalRingSize = face.vertices.size();
        compactRing(face.vertices, face.owner, face.neighbour, sourceRing);
        changed = changed || face.vertices.size() != originalRingSize;
        if (changed) {
            if (face.owner >= 0 && static_cast<std::size_t>(face.owner) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.owner)] = false;
            }
            if (face.neighbour >= 0 && static_cast<std::size_t>(face.neighbour) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.neighbour)] = false;
            }
        }
        double triangulatedArea = 0.0;
        if (face.vertices.size() >= 3U) {
            const auto& origin = mesh.points_m[static_cast<std::size_t>(face.vertices[0])];
            for (std::size_t vertex = 1; vertex + 1U < face.vertices.size(); ++vertex) {
                const auto& a = mesh.points_m[static_cast<std::size_t>(face.vertices[vertex])];
                const auto& b = mesh.points_m[static_cast<std::size_t>(face.vertices[vertex + 1U])];
                const std::array<double, 3> oa{a[0] - origin[0], a[1] - origin[1], a[2] - origin[2]};
                const std::array<double, 3> ob{b[0] - origin[0], b[1] - origin[1], b[2] - origin[2]};
                const std::array<double, 3> area{oa[1] * ob[2] - oa[2] * ob[1],
                    oa[2] * ob[0] - oa[0] * ob[2], oa[0] * ob[1] - oa[1] * ob[0]};
                triangulatedArea += 0.5 * std::hypot(area[0], area[1], area[2]);
            }
        }
        if (face.vertices.size() < 3U || triangulatedArea <= degenerateAreaTolerance) {
            ++mesh.cad_collapsed_faces_removed;
            collapsedFaces.push_back(std::move(face));
            continue;
        }
        // BabelSim's native reader uses a stable fan for convex rings.  Put
        // the fan root at the strongest corner so a legitimate collinear
        // refinement vertex along another edge cannot become a zero-area fan
        // triangle.  This is a cyclic rotation and preserves face orientation.
        if (face.vertices.size() > 3U) {
            std::size_t root = 0;
            double strongestCorner = -1.0;
            for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex) {
                const auto& a = mesh.points_m[static_cast<std::size_t>(face.vertices[(vertex + face.vertices.size() - 1U) % face.vertices.size()])];
                const auto& b = mesh.points_m[static_cast<std::size_t>(face.vertices[vertex])];
                const auto& c = mesh.points_m[static_cast<std::size_t>(face.vertices[(vertex + 1U) % face.vertices.size()])];
                const std::array<double, 3> ab{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
                const std::array<double, 3> cb{c[0] - b[0], c[1] - b[1], c[2] - b[2]};
                const std::array<double, 3> corner{ab[1] * cb[2] - ab[2] * cb[1],
                    ab[2] * cb[0] - ab[0] * cb[2], ab[0] * cb[1] - ab[1] * cb[0]};
                const double strength = std::hypot(corner[0], corner[1], corner[2]);
                if (strength > strongestCorner) { strongestCorner = strength; root = vertex; }
            }
            std::rotate(face.vertices.begin(), face.vertices.begin() + static_cast<std::ptrdiff_t>(root), face.vertices.end());
        }
        canonicalFaces.push_back(std::move(face));
    }
    mesh.faces = std::move(canonicalFaces);
    // A zero-area CAD feature patch can be the only face that carries a
    // subdivision point along an otherwise valid boundary edge. Removing the
    // collapsed patch without transferring that point leaves an open cell
    // shell. Propagate the collapsed ring's collinear vertices into every
    // surviving incident face that still has the corresponding unsplit edge.
    for (const auto& collapsed : collapsedFaces) {
        if (collapsed.vertices.size() < 3U) continue;
        std::set<EdgeKey> collapsedEdges;
        for (std::size_t edge = 0; edge < collapsed.vertices.size(); ++edge) {
            collapsedEdges.insert(edgeKey(collapsed.vertices[edge],
                collapsed.vertices[(edge + 1U) % collapsed.vertices.size()]));
        }
        for (auto& face : mesh.faces) {
            const bool sharesOwner = collapsed.owner >= 0 &&
                (face.owner == collapsed.owner || face.neighbour == collapsed.owner);
            const bool sharesNeighbour = collapsed.neighbour >= 0 &&
                (face.owner == collapsed.neighbour || face.neighbour == collapsed.neighbour);
            if (!sharesOwner && !sharesNeighbour) continue;
            std::vector<std::int32_t> expanded;
            expanded.reserve(face.vertices.size() + collapsed.vertices.size());
            bool repaired = false;
            for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
                const auto first = face.vertices[edge];
                const auto second = face.vertices[(edge + 1U) % face.vertices.size()];
                expanded.push_back(first);
                if (collapsedEdges.count(edgeKey(first, second)) == 0U) continue;
                const auto& a = mesh.points_m[static_cast<std::size_t>(first)];
                const auto& b = mesh.points_m[static_cast<std::size_t>(second)];
                const std::array<double, 3> direction{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                const double lengthSquared = direction[0] * direction[0] + direction[1] * direction[1] +
                    direction[2] * direction[2];
                if (!(lengthSquared > featureMergeTolerance * featureMergeTolerance)) continue;
                std::vector<std::pair<double, std::int32_t>> insertions;
                for (const auto vertex : collapsed.vertices) {
                    if (vertex == first || vertex == second ||
                        std::find(face.vertices.begin(), face.vertices.end(), vertex) != face.vertices.end()) continue;
                    const auto& point = mesh.points_m[static_cast<std::size_t>(vertex)];
                    const std::array<double, 3> offset{point[0] - a[0], point[1] - a[1], point[2] - a[2]};
                    const double parameter = (offset[0] * direction[0] + offset[1] * direction[1] +
                        offset[2] * direction[2]) / lengthSquared;
                    if (parameter <= 1.0e-12 || parameter >= 1.0 - 1.0e-12) continue;
                    const std::array<double, 3> residual{offset[0] - parameter * direction[0],
                        offset[1] - parameter * direction[1], offset[2] - parameter * direction[2]};
                    const double residualSquared = residual[0] * residual[0] + residual[1] * residual[1] +
                        residual[2] * residual[2];
                    if (residualSquared <= featureEdgeTolerance * featureEdgeTolerance) {
                        insertions.emplace_back(parameter, vertex);
                    }
                }
                std::sort(insertions.begin(), insertions.end(), [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) return lhs.first < rhs.first;
                    return lhs.second < rhs.second;
                });
                double previousParameter = -1.0;
                for (const auto& [parameter, vertex] : insertions) {
                    if (parameter - previousParameter <= 1.0e-12) continue;
                    if (expanded.back() != vertex) expanded.push_back(vertex);
                    repaired = true;
                    previousParameter = parameter;
                }
            }
            if (!repaired) continue;
            face.vertices = std::move(expanded);
            ++mesh.cad_collapsed_face_edge_repairs;
            if (face.owner >= 0 && static_cast<std::size_t>(face.owner) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.owner)] = false;
            }
            if (face.neighbour >= 0 && static_cast<std::size_t>(face.neighbour) < mesh.strict_hex.size()) {
                mesh.strict_hex[static_cast<std::size_t>(face.neighbour)] = false;
            }
        }
    }
    for (auto& cell : mesh.cells) {
        std::set<std::int32_t> unique;
        for (auto& vertex : cell) {
            if (vertex < 0 || static_cast<std::size_t>(vertex) >= pointRemap.size()) {
                throw std::runtime_error("CAD snapping canonicalization encountered an invalid cell vertex");
            }
            vertex = pointRemap[static_cast<std::size_t>(vertex)];
            unique.insert(vertex);
        }
        if (unique.size() < cell.size()) {
            const auto cellIndex = static_cast<std::size_t>(&cell - mesh.cells.data());
            if (cellIndex < mesh.strict_hex.size()) mesh.strict_hex[cellIndex] = false;
        }
    }
    std::vector<bool> usedPoints(mesh.points_m.size(), false);
    for (const auto& face : mesh.faces) for (const auto vertex : face.vertices) usedPoints[static_cast<std::size_t>(vertex)] = true;
    for (const auto& cell : mesh.cells) for (const auto vertex : cell) usedPoints[static_cast<std::size_t>(vertex)] = true;
    std::vector<std::int32_t> compactPointRemap(mesh.points_m.size(), -1);
    std::vector<std::array<double, 3>> compactPoints;
    for (std::size_t point = 0; point < mesh.points_m.size(); ++point) if (usedPoints[point] && pointRemap[point] == static_cast<std::int32_t>(point)) {
        compactPointRemap[point] = static_cast<std::int32_t>(compactPoints.size());
        compactPoints.push_back(mesh.points_m[point]);
    }
    for (std::size_t point = 0; point < pointRemap.size(); ++point) {
        const auto canonical = pointRemap[point];
        if (canonical >= 0 && static_cast<std::size_t>(canonical) < compactPointRemap.size()) {
            compactPointRemap[point] = compactPointRemap[static_cast<std::size_t>(canonical)];
        }
    }
    for (auto& face : mesh.faces) for (auto& vertex : face.vertices) vertex = compactPointRemap[static_cast<std::size_t>(vertex)];
    for (auto& cell : mesh.cells) for (auto& vertex : cell) vertex = compactPointRemap[static_cast<std::size_t>(vertex)];
    mesh.points_m = std::move(compactPoints);

    // Canonicalization can leave a long edge opposite two collinear subedges.
    // Propagate the existing subdivisions through every shared face before
    // accepting the cell shells; an area-vector sum alone cannot detect this.
    for (int pass = 0; pass < 4; ++pass) {
        const auto inserted = conformCollinearCellEdges(mesh, featureEdgeTolerance);
        mesh.cad_collapsed_face_edge_repairs += inserted;
        if (inserted == 0U) break;
    }
    mesh.cad_consolidated_faces += consolidateCoplanarFaces(mesh, featureEdgeTolerance, 0.05);
    requirePolyhedralShells(mesh);
    const auto finalQuality=auditMeshQuality(mesh);
    if(strictQuality && (finalQuality.negative_or_zero_volume_cells!=0 || finalQuality.degenerate_faces!=0 ||
       finalQuality.minimum_regular_hex_scaled_jacobian<0.1 ||
       finalQuality.maximum_non_orthogonality_degrees>75.0 || finalQuality.p99_non_orthogonality_degrees>60.0 ||
       finalQuality.maximum_internal_skewness>0.5 || finalQuality.p99_internal_skewness>0.25))
        throw std::runtime_error("final CAD feature canonicalization failed strict mesh-quality gate: min_scaled_jacobian="+
            std::to_string(finalQuality.minimum_regular_hex_scaled_jacobian)+", max/p99 nonorth="+
            std::to_string(finalQuality.maximum_non_orthogonality_degrees)+"/"+
            std::to_string(finalQuality.p99_non_orthogonality_degrees)+", max/p99 skew="+
            std::to_string(finalQuality.maximum_internal_skewness)+"/"+
            std::to_string(finalQuality.p99_internal_skewness));
    // Consolidation changes face-interior sample locations. Audit the final
    // exported boundary rather than retaining the pre-consolidation samples.
    const auto finalDeviation = maximumDeviation();
    if (strictQuality && finalDeviation.maximum > maximumSurfaceDeviation)
        throw std::runtime_error("consolidated CAD boundary exceeds surface deviation at face " +
            std::to_string(finalDeviation.meshFace));
    acceptedDeviation = finalDeviation.maximum;
    deviationSamples = finalDeviation.samples;

    auto finalVolumes = cellVolumes();
    for (std::size_t cell = 0; cell < finalVolumes.size(); ++cell) {
        const double reference = cell < acceptedVolumes.size() ? acceptedVolumes[cell] : 0.0;
        const double lower = std::max(1.0e-30, std::abs(reference) * 1.0e-12);
        if (!std::isfinite(finalVolumes[cell]) || finalVolumes[cell] <= lower) {
            throw std::runtime_error("CAD feature canonicalization produced a nonpositive polyhedral cell volume at cell " +
                std::to_string(cell));
        }
    }
    mesh.cell_volumes_m3 = std::move(finalVolumes);
    mesh.cad_snapping_applied = acceptedFraction > 0.0;
    mesh.cad_snap_fraction = acceptedFraction;
    mesh.cad_core_deformation_fraction = acceptedCoreDeformationFraction;
    mesh.cad_surface_deviation_m = acceptedDeviation;
    mesh.cad_snapping_vertices = targets.size();
    std::set<std::size_t> projectedFaces;
    for (const auto& target : targets) projectedFaces.insert(target.faceIndex);
    mesh.cad_projection_faces = projectedFaces.size();
    mesh.cad_surface_deviation_samples = deviationSamples;
    if (mesh.cad_fluid_domain_volume_m3 > 0.0) {
        const double volume = std::accumulate(mesh.cell_volumes_m3.begin(), mesh.cell_volumes_m3.end(), 0.0);
        mesh.cad_mesh_volume_relative_error = std::abs(volume - mesh.cad_fluid_domain_volume_m3) / mesh.cad_fluid_domain_volume_m3;
    }
    (void)deviationSamples;
}

BoundaryProjection projectPointToFluidBoundary(const CadModel& cadModel,
                                               const std::array<double, 3>& point_m,
                                               double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) throw std::runtime_error("boundary projection tolerance must be finite and nonnegative");
    for (const auto coordinate : point_m) if (!std::isfinite(coordinate)) {
        throw std::runtime_error("boundary projection point must be finite");
    }
    const TopoDS_Vertex query = BRepBuilderAPI_MakeVertex(gp_Pnt(point_m[0], point_m[1], point_m[2])).Vertex();
    BoundaryProjection result;
    result.distance_m = std::numeric_limits<double>::infinity();
    std::vector<TopoDS_Face> faces;
    for (TopExp_Explorer explorer(cadModel.fluidDomain(), TopAbs_FACE); explorer.More(); explorer.Next()) {
        faces.push_back(TopoDS::Face(explorer.Current()));
    }
    if (faces.empty()) throw std::runtime_error("fluid domain has no boundary faces");
    auto candidates = nearestSurfaceFaces(cadModel.fluidSurfaceMesh(), point_m);
    if (candidates.empty()) {
        candidates.resize(faces.size());
        std::iota(candidates.begin(), candidates.end(), 0U);
    }
    try {
        for (const auto faceIndex : candidates) {
            if (faceIndex >= faces.size()) throw std::runtime_error("fluid surface triangulation references an invalid CAD face");
            const TopoDS_Face& face = faces[faceIndex];
            Bnd_Box bounds;
            BRepBndLib::Add(face, bounds, Standard_False);
            if (bounds.IsVoid()) continue;
            Standard_Real xmin = 0.0, ymin = 0.0, zmin = 0.0, xmax = 0.0, ymax = 0.0, zmax = 0.0;
            bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
            const double dx = point_m[0] < xmin ? xmin - point_m[0] : (point_m[0] > xmax ? point_m[0] - xmax : 0.0);
            const double dy = point_m[1] < ymin ? ymin - point_m[1] : (point_m[1] > ymax ? point_m[1] - ymax : 0.0);
            const double dz = point_m[2] < zmin ? zmin - point_m[2] : (point_m[2] > zmax ? point_m[2] - zmax : 0.0);
            if (std::sqrt(dx * dx + dy * dy + dz * dz) > result.distance_m + tolerance) continue;
            BRepExtrema_DistShapeShape distance(query, face);
            distance.Perform();
            if (!distance.IsDone()) throw std::runtime_error("OCCT failed to project a mesh point to a CAD face");
            if (distance.Value() < result.distance_m) {
                const gp_Pnt projected = distance.PointOnShape2(1);
                result.point_m = {projected.X(), projected.Y(), projected.Z()};
                result.distance_m = distance.Value();
                result.face_index = faceIndex;
            }
        }
    } catch (const Standard_Failure& error) {
        const char* message = error.GetMessageString();
        throw std::runtime_error(std::string("OCCT fluid-boundary projection failed: ") + (message ? message : "unknown OCCT failure"));
    }
    if (!std::isfinite(result.distance_m)) throw std::runtime_error("fluid domain has no projectable boundary face");
    return result;
}

}  // namespace babelmesh
