#include "babelmesh/cartesian_mesh.hpp"
#include "babelmesh/polyhedral_audit.hpp"

#include "babelmesh/cad.hpp"
#include "babelmesh/fluid_domain.hpp"
#include "babelmesh/surface_distance.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>

namespace babelmesh {
namespace {

using Json = nlohmann::json;

double lengthScale(const std::string& unit) {
    if (unit == "m") return 1.0;
    if (unit == "mm") return 1.0e-3;
    if (unit == "cm") return 1.0e-2;
    if (unit == "inch") return 0.0254;
    throw std::runtime_error("unsupported project length unit '" + unit + "'");
}

void atomicWrite(const std::filesystem::path& file, const std::function<void(std::ostream&)>& writer) {
    std::filesystem::create_directories(file.parent_path());
    const auto temporary = file.string() + ".tmp";
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("cannot create output file: " + temporary);
        writer(output);
        output.flush();
        if (!output) throw std::runtime_error("failed while writing output file: " + temporary);
    }
    std::error_code error;
    std::filesystem::rename(temporary, file, error);
    if (error) {
        std::filesystem::remove(temporary);
        throw std::runtime_error("cannot atomically replace output file '" + file.string() + "': " + error.message());
    }
}

std::int32_t cellId(const HexMesh& mesh, std::int32_t i, std::int32_t j, std::int32_t k) {
    const auto nx = mesh.cell_counts[0];
    const auto ny = mesh.cell_counts[1];
    const std::int64_t id = static_cast<std::int64_t>(i) + static_cast<std::int64_t>(nx) *
        (static_cast<std::int64_t>(j) + static_cast<std::int64_t>(ny) * k);
    if (id > std::numeric_limits<std::int32_t>::max()) throw std::runtime_error("Cartesian cell index exceeds BabelSim's 32-bit limit");
    return static_cast<std::int32_t>(id);
}

std::int32_t pointId(const HexMesh& mesh, std::int32_t i, std::int32_t j, std::int32_t k) {
    const auto nx = mesh.cell_counts[0] + 1;
    const auto ny = mesh.cell_counts[1] + 1;
    const std::int64_t id = static_cast<std::int64_t>(i) + static_cast<std::int64_t>(nx) *
        (static_cast<std::int64_t>(j) + static_cast<std::int64_t>(ny) * k);
    if (id > std::numeric_limits<std::int32_t>::max()) throw std::runtime_error("Cartesian point index exceeds BabelSim's 32-bit limit");
    return static_cast<std::int32_t>(id);
}

long double estimateCartesianStorageBytes(const HexMesh& mesh) {
    const long double nx = mesh.cell_counts[0], ny = mesh.cell_counts[1], nz = mesh.cell_counts[2];
    const long double pointCount = (nx + 1.0L) * (ny + 1.0L) * (nz + 1.0L);
    const long double cellCount = nx * ny * nz;
    const long double faceCount = (nx + 1.0L) * ny * nz + nx * (ny + 1.0L) * nz + nx * ny * (nz + 1.0L);
    return pointCount * sizeof(std::array<double, 3>) + cellCount * sizeof(std::array<std::int32_t, 8>) +
        faceCount * (sizeof(MeshFace) + 4.0L * sizeof(std::int32_t));
}

void requireObject(const Json& value, const std::string& where) {
    if (!value.is_object()) throw std::runtime_error(where + " must be an object");
}

void requireExactKeys(const Json& value, const std::vector<std::string>& allowed, const std::string& where) {
    requireObject(value, where);
    for (auto it = value.begin(); it != value.end(); ++it) {
        if (std::find(allowed.begin(), allowed.end(), it.key()) == allowed.end()) {
            throw std::runtime_error(where + " has unknown key '" + it.key() + "'");
        }
    }
}

void emitFace(HexMesh& mesh, std::vector<std::int32_t> ring,
              std::int32_t owner, std::int32_t neighbour, std::int32_t patch) {
    if (neighbour < 0) {
        if (patch < 0 || static_cast<std::size_t>(patch) >= mesh.patches.size()) {
            throw std::runtime_error("Cartesian boundary face has an invalid patch index");
        }
    } else if (patch != -1) {
        throw std::runtime_error("Cartesian internal face must not have a boundary patch");
    }
    mesh.faces.push_back({std::move(ring), owner, neighbour, patch});
}

void makeXFaces(HexMesh& mesh) {
    const auto nx = mesh.cell_counts[0], ny = mesh.cell_counts[1], nz = mesh.cell_counts[2];
    for (std::int32_t i = 0; i <= nx; ++i) for (std::int32_t k = 0; k < nz; ++k) for (std::int32_t j = 0; j < ny; ++j) {
        std::vector<std::int32_t> ring{pointId(mesh, i, j, k), pointId(mesh, i, j + 1, k),
            pointId(mesh, i, j + 1, k + 1), pointId(mesh, i, j, k + 1)};
        std::int32_t owner, neighbour = -1, patch = -1;
        if (i == 0) { std::reverse(ring.begin(), ring.end()); owner = cellId(mesh, 0, j, k); patch = 0; }
        else if (i == nx) { owner = cellId(mesh, nx - 1, j, k); patch = 1; }
        else { owner = cellId(mesh, i - 1, j, k); neighbour = cellId(mesh, i, j, k); }
        emitFace(mesh, std::move(ring), owner, neighbour, patch);
    }
}

void makeYFaces(HexMesh& mesh) {
    const auto nx = mesh.cell_counts[0], ny = mesh.cell_counts[1], nz = mesh.cell_counts[2];
    for (std::int32_t j = 0; j <= ny; ++j) for (std::int32_t k = 0; k < nz; ++k) for (std::int32_t i = 0; i < nx; ++i) {
        std::vector<std::int32_t> ring{pointId(mesh, i, j, k), pointId(mesh, i, j, k + 1),
            pointId(mesh, i + 1, j, k + 1), pointId(mesh, i + 1, j, k)};
        std::int32_t owner, neighbour = -1, patch = -1;
        if (j == 0) { std::reverse(ring.begin(), ring.end()); owner = cellId(mesh, i, 0, k); patch = 2; }
        else if (j == ny) { owner = cellId(mesh, i, ny - 1, k); patch = 3; }
        else { owner = cellId(mesh, i, j - 1, k); neighbour = cellId(mesh, i, j, k); }
        emitFace(mesh, std::move(ring), owner, neighbour, patch);
    }
}

void makeZFaces(HexMesh& mesh) {
    const auto nx = mesh.cell_counts[0], ny = mesh.cell_counts[1], nz = mesh.cell_counts[2];
    for (std::int32_t k = 0; k <= nz; ++k) for (std::int32_t j = 0; j < ny; ++j) for (std::int32_t i = 0; i < nx; ++i) {
        std::vector<std::int32_t> ring{pointId(mesh, i, j, k), pointId(mesh, i + 1, j, k),
            pointId(mesh, i + 1, j + 1, k), pointId(mesh, i, j + 1, k)};
        std::int32_t owner, neighbour = -1, patch = -1;
        if (k == 0) { std::reverse(ring.begin(), ring.end()); owner = cellId(mesh, i, j, 0); patch = 4; }
        else if (k == nz) { owner = cellId(mesh, i, j, nz - 1); patch = 5; }
        else { owner = cellId(mesh, i, j, k - 1); neighbour = cellId(mesh, i, j, k); }
        emitFace(mesh, std::move(ring), owner, neighbour, patch);
    }
}

using LogicalPoint = std::array<std::int64_t, 3>;

struct OctreeNode {
    LogicalPoint lower{};
    std::int64_t step{};
    std::int32_t level{};
    std::array<std::int32_t, 8> children{{-1, -1, -1, -1, -1, -1, -1, -1}};
    bool leaf{true};
};

class OctreeForest {
public:
    OctreeForest(std::array<std::int32_t, 3> roots, std::int64_t rootStep)
        : m_rootCounts(roots), m_rootStep(rootStep) {
        const auto count = static_cast<std::size_t>(roots[0]) * static_cast<std::size_t>(roots[1]) * static_cast<std::size_t>(roots[2]);
        m_roots.resize(count);
        for (std::int32_t k = 0; k < roots[2]; ++k) for (std::int32_t j = 0; j < roots[1]; ++j) for (std::int32_t i = 0; i < roots[0]; ++i) {
            OctreeNode root;
            root.lower = {static_cast<std::int64_t>(i) * rootStep, static_cast<std::int64_t>(j) * rootStep,
                static_cast<std::int64_t>(k) * rootStep};
            root.step = rootStep;
            m_roots[rootIndex(i, j, k)] = addNode(root);
        }
        m_leafCount = m_roots.size();
    }

    const OctreeNode& node(std::int32_t id) const { return m_nodes.at(static_cast<std::size_t>(id)); }
    OctreeNode& node(std::int32_t id) { return m_nodes.at(static_cast<std::size_t>(id)); }
    std::int64_t rootStep() const { return m_rootStep; }
    std::size_t leafCount() const { return m_leafCount; }
    std::size_t nodeCount() const { return m_nodes.size(); }

    std::vector<std::int32_t> leaves() const {
        std::vector<std::int32_t> result;
        result.reserve(m_leafCount);
        for (const auto root : m_roots) collectLeaves(root, result);
        return result;
    }

    std::int32_t locate(const LogicalPoint& point) const {
        std::array<std::int64_t, 3> rootCoordinate{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (point[axis] < 0 || point[axis] >= static_cast<std::int64_t>(m_rootCounts[axis]) * m_rootStep) return -1;
            rootCoordinate[axis] = point[axis] / m_rootStep;
        }
        std::int32_t current = m_roots[rootIndex(static_cast<std::int32_t>(rootCoordinate[0]),
            static_cast<std::int32_t>(rootCoordinate[1]), static_cast<std::int32_t>(rootCoordinate[2]))];
        while (!node(current).leaf) {
            const OctreeNode& parent = node(current);
            const std::int64_t half = parent.step / 2;
            std::size_t child = 0;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                if (point[axis] >= parent.lower[axis] + half) child |= (1U << axis);
            }
            current = parent.children[child];
            if (current < 0) throw std::runtime_error("octree internal node has a missing child");
        }
        return current;
    }

    std::array<std::int32_t, 8> split(std::int32_t id) {
        OctreeNode& parent = node(id);
        if (!parent.leaf || parent.step < 2) throw std::runtime_error("invalid octree split request");
        const std::int64_t half = parent.step / 2;
        const auto lower = parent.lower;
        const auto level = parent.level + 1;
        std::array<std::int32_t, 8> children{};
        for (std::size_t child = 0; child < children.size(); ++child) {
            OctreeNode item;
            item.lower = lower;
            for (std::size_t axis = 0; axis < 3; ++axis) if ((child & (1U << axis)) != 0U) item.lower[axis] += half;
            item.step = half;
            item.level = level;
            children[child] = addNode(item);
        }
        OctreeNode& updatedParent = node(id);
        updatedParent.children = children;
        updatedParent.leaf = false;
        m_leafCount += 7U;
        return children;
    }

private:
    std::size_t rootIndex(std::int32_t i, std::int32_t j, std::int32_t k) const {
        return static_cast<std::size_t>(i) + static_cast<std::size_t>(m_rootCounts[0]) *
            (static_cast<std::size_t>(j) + static_cast<std::size_t>(m_rootCounts[1]) * static_cast<std::size_t>(k));
    }

    std::int32_t addNode(const OctreeNode& item) {
        if (m_nodes.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("octree node count exceeds BabelSim's 32-bit index limit");
        }
        m_nodes.push_back(item);
        return static_cast<std::int32_t>(m_nodes.size() - 1U);
    }

    void collectLeaves(std::int32_t id, std::vector<std::int32_t>& output) const {
        const OctreeNode& item = node(id);
        if (item.leaf) { output.push_back(id); return; }
        for (const auto child : item.children) collectLeaves(child, output);
    }

    std::array<std::int32_t, 3> m_rootCounts;
    std::int64_t m_rootStep;
    std::vector<std::int32_t> m_roots;
    std::vector<OctreeNode> m_nodes;
    std::size_t m_leafCount{};
};

struct SizingSource {
    std::string id;
    std::string kind;
    double targetSize{};
    double curvatureMaxRadius{};
    std::array<double, 6> bounds{};
    bool hasSphere{};
    std::array<double, 3> sphereCenter{};
    double sphereRadius{};
    bool hasCylinder{};
    std::array<double, 3> cylinderCenter{};
    std::array<double, 3> cylinderAxis{};
    double cylinderHalfLength{};
    double cylinderRadius{};
    double distance{};
    std::string cadGroup;
    std::string otherCadGroup;
    std::vector<TopoDS_Shape> cadShapes;
    std::vector<std::array<double, 6>> cadShapeBounds;
    std::vector<TopoDS_Shape> otherCadShapes;
    std::vector<std::array<double, 6>> otherCadShapeBounds;
};

std::array<double, 6> physicalBounds(const OctreeNode& node, const HexMesh& mesh, std::int64_t rootStep) {
    std::array<double, 6> result{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double rootSize = (mesh.bounds_m[axis + 3U] - mesh.bounds_m[axis]) / mesh.cell_counts[axis];
        result[axis] = mesh.bounds_m[axis] + static_cast<double>(node.lower[axis]) * rootSize / static_cast<double>(rootStep);
        result[axis + 3U] = mesh.bounds_m[axis] + static_cast<double>(node.lower[axis] + node.step) * rootSize / static_cast<double>(rootStep);
    }
    return result;
}

bool boxesIntersect(const std::array<double, 6>& a, const std::array<double, 6>& b) {
    for (std::size_t axis = 0; axis < 3; ++axis) if (a[axis + 3U] < b[axis] || b[axis + 3U] < a[axis]) return false;
    return true;
}

std::array<double, 6> shapeBounds(const TopoDS_Shape& shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box, Standard_False);
    if (box.IsVoid()) throw std::runtime_error("CAD sizing group contains a shape with empty bounds");
    std::array<double, 6> result{};
    box.Get(result[0], result[1], result[2], result[3], result[4], result[5]);
    return result;
}

std::array<double, 6> expandedBounds(const std::array<double, 6>& bounds, double amount) {
    auto result = bounds;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        result[axis] -= amount;
        result[axis + 3U] += amount;
    }
    return result;
}

const PhysicalGroup& physicalGroup(const CadModel* cadModel, const std::string& name,
                                   const std::string& sourceId) {
    if (cadModel == nullptr) {
        throw std::runtime_error("sizing source '" + sourceId + "' references CAD group '" + name +
            "' but no validated CAD model was supplied");
    }
    const auto found = std::find_if(cadModel->groups().begin(), cadModel->groups().end(),
        [&](const PhysicalGroup& group) { return group.name == name; });
    if (found == cadModel->groups().end()) {
        throw std::runtime_error("sizing source '" + sourceId + "' refers to unknown CAD group '" + name + "'");
    }
    return *found;
}

std::vector<TopoDS_Shape> groupShapes(const CadModel* cadModel, const std::string& name,
                                      const std::string& sourceId, bool facesOnly) {
    std::vector<TopoDS_Shape> result;
    const auto& group = physicalGroup(cadModel, name, sourceId);
    for (const auto& shape : group.selected_subshapes) {
        if (!facesOnly) {
            result.push_back(shape);
        } else if (shape.ShapeType() == TopAbs_FACE) {
            result.push_back(shape);
        } else {
            for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
                result.push_back(explorer.Current());
            }
        }
    }
    if (result.empty()) {
        throw std::runtime_error("sizing source '" + sourceId + "' CAD group '" + name +
            "' has no usable selected shapes");
    }
    return result;
}

std::array<double, 6> unionBounds(const std::vector<std::array<double, 6>>& bounds) {
    if (bounds.empty()) throw std::runtime_error("cannot compute bounds of an empty CAD group");
    auto result = bounds.front();
    for (std::size_t index = 1; index < bounds.size(); ++index) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            result[axis] = std::min(result[axis], bounds[index][axis]);
            result[axis + 3U] = std::max(result[axis + 3U], bounds[index][axis + 3U]);
        }
    }
    return result;
}

bool faceHasCurvatureAtRadius(const TopoDS_Shape& shape, double maximumRadius) {
    if (shape.ShapeType() != TopAbs_FACE) return false;
    const TopoDS_Face face = TopoDS::Face(shape);
    BRepAdaptor_Surface surface(face, Standard_True);
    Standard_Real uMin = 0.0, uMax = 0.0, vMin = 0.0, vMax = 0.0;
    BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);
    if (!std::isfinite(uMin) || !std::isfinite(uMax) || !std::isfinite(vMin) || !std::isfinite(vMax) ||
        !(uMax > uMin) || !(vMax > vMin)) return false;
    const std::array<double, 3> fractions{0.25, 0.5, 0.75};
    for (const double fu : fractions) for (const double fv : fractions) {
        const double u = uMin + fu * (uMax - uMin);
        const double v = vMin + fv * (vMax - vMin);
        BRepClass_FaceClassifier classifier(face, gp_Pnt2d(u, v), 1.0e-9);
        if (classifier.State() != TopAbs_IN && classifier.State() != TopAbs_ON) continue;
        BRepLProp_SLProps properties(surface, u, v, 2, 1.0e-10);
        if (!properties.IsCurvatureDefined()) continue;
        const double curvature = std::max(std::abs(properties.MaxCurvature()), std::abs(properties.MinCurvature()));
        if (std::isfinite(curvature) && curvature >= 1.0 / maximumRadius) return true;
    }
    return false;
}

TopoDS_Shape cellBoxShape(const std::array<double, 6>& bounds) {
    return BRepPrimAPI_MakeBox(gp_Pnt(bounds[0], bounds[1], bounds[2]),
        gp_Pnt(bounds[3], bounds[4], bounds[5])).Shape();
}

double distanceToShapes(const TopoDS_Shape& cellShape, const std::array<double, 6>& cellBounds,
                        const std::vector<TopoDS_Shape>& shapes,
                        const std::vector<std::array<double, 6>>& bounds) {
    const auto boxDistance = [&](const std::array<double, 6>& candidate) {
        double squared = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double delta = candidate[axis + 3U] < cellBounds[axis]
                ? cellBounds[axis] - candidate[axis + 3U]
                : (cellBounds[axis + 3U] < candidate[axis]
                    ? candidate[axis] - cellBounds[axis + 3U] : 0.0);
            squared += delta * delta;
        }
        return std::sqrt(squared);
    };
    std::size_t nearestBounds = 0;
    double minimum = boxDistance(bounds.front());
    for (std::size_t index = 1; index < bounds.size(); ++index) {
        const double candidate = boxDistance(bounds[index]);
        if (candidate < minimum) { minimum = candidate; nearestBounds = index; }
    }
    auto exactDistance = [&](std::size_t index) {
        BRepExtrema_DistShapeShape distance(cellShape, shapes[index]);
        distance.Perform();
        if (!distance.IsDone()) throw std::runtime_error("OCCT failed to evaluate distance to a CAD sizing group");
        return distance.Value();
    };
    minimum = exactDistance(nearestBounds);
    for (std::size_t index = 0; index < shapes.size(); ++index) {
        if (index == nearestBounds || boxDistance(bounds[index]) >= minimum) continue;
        minimum = std::min(minimum, exactDistance(index));
    }
    return minimum;
}

bool sourceIntersects(const SizingSource& source, const std::array<double, 6>& cellBounds) {
    if (!boxesIntersect(source.bounds, cellBounds)) return false;
    if (source.hasSphere) {
        double squaredDistance = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double nearest = std::clamp(source.sphereCenter[axis], cellBounds[axis], cellBounds[axis + 3U]);
            const double delta = nearest - source.sphereCenter[axis];
            squaredDistance += delta * delta;
        }
        if (squaredDistance > source.sphereRadius * source.sphereRadius) return false;
    }
    if (source.hasCylinder) {
        const auto distanceSquared = [&](double parameter) {
            double value = 0.0;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double point = source.cylinderCenter[axis] + parameter * source.cylinderAxis[axis];
                const double nearest = std::clamp(point, cellBounds[axis], cellBounds[axis + 3U]);
                const double delta = point - nearest;
                value += delta * delta;
            }
            return value;
        };
        const auto derivative = [&](double parameter) {
            double value = 0.0;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double point = source.cylinderCenter[axis] + parameter * source.cylinderAxis[axis];
                if (point < cellBounds[axis]) value += source.cylinderAxis[axis] * (point - cellBounds[axis]);
                else if (point > cellBounds[axis + 3U]) value += source.cylinderAxis[axis] * (point - cellBounds[axis + 3U]);
            }
            return value;
        };
        double lower = -source.cylinderHalfLength, upper = source.cylinderHalfLength;
        if (derivative(lower) < 0.0 && derivative(upper) > 0.0) {
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double middle = 0.5 * (lower + upper);
                if (derivative(middle) < 0.0) lower = middle;
                else upper = middle;
            }
        } else if (derivative(lower) >= 0.0) upper = lower;
        else lower = upper;
        if (distanceSquared(0.5 * (lower + upper)) > source.cylinderRadius * source.cylinderRadius) return false;
    }
    if (!source.cadShapes.empty()) {
        const auto cellShape = cellBoxShape(cellBounds);
        const double firstDistance = distanceToShapes(cellShape, cellBounds, source.cadShapes, source.cadShapeBounds);
        if (firstDistance > source.distance) return false;
        if (!source.otherCadShapes.empty()) {
            const double secondDistance = distanceToShapes(cellShape, cellBounds,
                source.otherCadShapes, source.otherCadShapeBounds);
            if (secondDistance > source.distance) return false;
        }
    }
    return true;
}

std::vector<SizingSource> parseSizingSources(const Json& project, double scale, const CadModel* cadModel) {
    std::vector<SizingSource> result;
    const double minimumSize = project.at("sizing").at("minimum_size").get<double>() * scale;
    for (const auto& item : project.at("sizing").at("sources")) {
        SizingSource source;
        source.id = item.at("id").get<std::string>();
        source.kind = item.at("kind").get<std::string>();
        source.targetSize = item.at("target_size").get<double>() * scale;
        if (!std::isfinite(source.targetSize) || source.targetSize < minimumSize) {
            throw std::runtime_error("sizing source '" + source.id + "' target_size is below sizing.minimum_size or not finite");
        }
        if (!item.contains("region")) throw std::runtime_error("sizing source '" + source.id + "' needs a region in M3");
        const Json& region = item.at("region");
        if (region.contains("bounds")) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                source.bounds[axis] = region.at("bounds").at("minimum").at(axis).get<double>() * scale;
                source.bounds[axis + 3U] = region.at("bounds").at("maximum").at(axis).get<double>() * scale;
                if (!std::isfinite(source.bounds[axis]) || !std::isfinite(source.bounds[axis + 3U]) ||
                    source.bounds[axis + 3U] < source.bounds[axis]) {
                    throw std::runtime_error("sizing source '" + source.id + "' has invalid bounds");
                }
            }
        }
        if (source.kind == "sphere") {
            for (std::size_t axis = 0; axis < 3; ++axis) source.sphereCenter[axis] = region.at("center").at(axis).get<double>() * scale;
            source.sphereRadius = region.at("radius").get<double>() * scale;
            if (!std::isfinite(source.sphereRadius) || source.sphereRadius <= 0.0) throw std::runtime_error("sizing sphere radius must be finite and positive");
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double lower = source.sphereCenter[axis] - source.sphereRadius;
                const double upper = source.sphereCenter[axis] + source.sphereRadius;
                if (region.contains("bounds")) {
                    source.bounds[axis] = std::min(source.bounds[axis], lower);
                    source.bounds[axis + 3U] = std::max(source.bounds[axis + 3U], upper);
                } else {
                    source.bounds[axis] = lower;
                    source.bounds[axis + 3U] = upper;
                }
            }
            source.hasSphere = true;
        } else if (source.kind == "cylinder") {
            double norm = 0.0;
            for (std::size_t d = 0; d < 3; ++d) { source.cylinderAxis[d] = region.at("axis").at(d).get<double>(); norm += source.cylinderAxis[d] * source.cylinderAxis[d]; }
            norm = std::sqrt(norm);
            if (!std::isfinite(norm) || norm <= 1.0e-15) throw std::runtime_error("sizing cylinder axis must be nonzero");
            for (double& component : source.cylinderAxis) component /= norm;
            source.cylinderRadius = region.at("radius").get<double>() * scale;
            source.cylinderHalfLength = 0.5 * region.at("length").get<double>() * scale;
            if (!std::isfinite(source.cylinderRadius) || source.cylinderRadius <= 0.0 ||
                !std::isfinite(source.cylinderHalfLength) || source.cylinderHalfLength <= 0.0) {
                throw std::runtime_error("sizing cylinder radius and length must be finite and positive");
            }
            for (std::size_t d = 0; d < 3; ++d) {
                source.cylinderCenter[d] = region.at("center").at(d).get<double>() * scale;
                const double extent = std::abs(source.cylinderAxis[d]) * source.cylinderHalfLength +
                    source.cylinderRadius * std::sqrt(std::max(0.0, 1.0 - source.cylinderAxis[d] * source.cylinderAxis[d]));
                const double lower = source.cylinderCenter[d] - extent;
                const double upper = source.cylinderCenter[d] + extent;
                if (region.contains("bounds")) {
                    source.bounds[d] = std::min(source.bounds[d], lower);
                    source.bounds[d + 3U] = std::max(source.bounds[d + 3U], upper);
                } else {
                    source.bounds[d] = lower;
                    source.bounds[d + 3U] = upper;
                }
            }
            source.hasCylinder = true;
        } else if (source.kind == "box" && !region.contains("bounds")) {
            throw std::runtime_error("box sizing source '" + source.id + "' requires region.bounds");
        }

        if (source.kind == "cad_group" || source.kind == "distance" ||
            source.kind == "curvature" || source.kind == "gap") {
            if (!region.contains("cad_group")) {
                throw std::runtime_error("CAD sizing source '" + source.id + "' requires region.cad_group");
            }
            if (source.kind == "distance" || source.kind == "gap") {
                if (!region.contains("distance")) {
                    throw std::runtime_error("sizing source '" + source.id + "' requires a positive region.distance");
                }
                source.distance = region.at("distance").get<double>() * scale;
                if (!std::isfinite(source.distance) || source.distance <= 0.0) {
                    throw std::runtime_error("sizing source '" + source.id + "' region.distance must be finite and positive");
                }
            } else if (source.kind == "cad_group") {
                source.distance = region.value("distance", 0.0) * scale;
            } else {
                source.distance = region.value("distance", item.at("target_size").get<double>()) * scale;
            }
            if (!std::isfinite(source.distance) || source.distance < 0.0) {
                throw std::runtime_error("sizing source '" + source.id + "' region.distance must be finite and nonnegative");
            }

            const bool curvature = source.kind == "curvature";
            source.cadGroup = region.at("cad_group").get<std::string>();
            if (curvature) {
                if (!region.contains("max_radius")) {
                    throw std::runtime_error("curvature sizing source '" + source.id +
                        "' requires region.max_radius; target_size is the cell-size target");
                }
                source.curvatureMaxRadius = region.at("max_radius").get<double>() * scale;
                if (!std::isfinite(source.curvatureMaxRadius) || source.curvatureMaxRadius <= 0.0) {
                    throw std::runtime_error("curvature sizing source '" + source.id + "' region.max_radius must be finite and positive");
                }
            }
            source.cadShapes = groupShapes(cadModel, source.cadGroup, source.id, curvature);
            if (curvature) {
                source.cadShapes.erase(std::remove_if(source.cadShapes.begin(), source.cadShapes.end(),
                    [&](const TopoDS_Shape& shape) { return !faceHasCurvatureAtRadius(shape, source.curvatureMaxRadius); }),
                    source.cadShapes.end());
                if (source.cadShapes.empty()) {
                    throw std::runtime_error("curvature sizing source '" + source.id +
                        "' found no selected CAD face with radius of curvature at or below region.max_radius");
                }
            }
            for (const auto& shape : source.cadShapes) source.cadShapeBounds.push_back(shapeBounds(shape));
            source.bounds = expandedBounds(unionBounds(source.cadShapeBounds), source.distance);

            if (source.kind == "gap") {
                if (!region.contains("other_cad_group")) {
                    throw std::runtime_error("gap sizing source '" + source.id + "' requires region.other_cad_group");
                }
                source.otherCadShapes = groupShapes(cadModel,
                    region.at("other_cad_group").get<std::string>(), source.id, false);
                source.otherCadGroup = region.at("other_cad_group").get<std::string>();
                for (const auto& shape : source.otherCadShapes) source.otherCadShapeBounds.push_back(shapeBounds(shape));
                const auto otherBounds = expandedBounds(unionBounds(source.otherCadShapeBounds), source.distance);
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    source.bounds[axis] = std::max(source.bounds[axis], otherBounds[axis]);
                    source.bounds[axis + 3U] = std::min(source.bounds[axis + 3U], otherBounds[axis + 3U]);
                    if (source.bounds[axis + 3U] < source.bounds[axis]) {
                        throw std::runtime_error("gap sizing source '" + source.id +
                            "' groups have no overlapping distance bands; increase region.distance");
                    }
                }
            }
        } else if (source.kind != "box" && source.kind != "sphere" && source.kind != "cylinder") {
            throw std::runtime_error("unknown sizing source kind '" + source.kind + "'");
        }
        result.push_back(std::move(source));
    }
    return result;
}

std::array<double, 3> logicalToPhysical(const LogicalPoint& point, const HexMesh& mesh, std::int64_t rootStep) {
    std::array<double, 3> result{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double rootSize = (mesh.bounds_m[axis + 3U] - mesh.bounds_m[axis]) / mesh.cell_counts[axis];
        result[axis] = mesh.bounds_m[axis] + static_cast<double>(point[axis]) * rootSize / static_cast<double>(rootStep);
    }
    return result;
}

std::vector<std::int32_t> touchingLeaves(const OctreeForest& tree, const OctreeNode& cell, int axis, int side, bool faceOnly) {
    std::set<std::int32_t> found;
    if (!faceOnly) return {};
    const std::vector<std::int64_t> fractions = cell.step >= 4 && cell.step % 4 == 0 ?
        std::vector<std::int64_t>{1, 2, 3} : std::vector<std::int64_t>{2};
    std::array<int, 2> tangent{};
    int t = 0;
    for (int d = 0; d < 3; ++d) if (d != axis) tangent[static_cast<std::size_t>(t++)] = d;
    for (const auto f0 : fractions) for (const auto f1 : fractions) {
        LogicalPoint probe{};
        for (int d = 0; d < 3; ++d) {
            if (d == axis) probe[static_cast<std::size_t>(d)] = side > 0 ? cell.lower[static_cast<std::size_t>(d)] + cell.step + 1 : cell.lower[static_cast<std::size_t>(d)] - 1;
            else {
                const int tangentIndex = tangent[0] == d ? 0 : 1;
                const auto fraction = tangentIndex == 0 ? f0 : f1;
                probe[static_cast<std::size_t>(d)] = cell.lower[static_cast<std::size_t>(d)] + cell.step * fraction / 4;
            }
        }
        const auto neighbor = tree.locate(probe);
        if (neighbor >= 0) found.insert(neighbor);
    }
    return {found.begin(), found.end()};
}

std::vector<std::int32_t> adjacentLeaves(const OctreeForest& tree, const OctreeNode& cell) {
    std::set<std::int32_t> found;
    const auto sampleFractions = cell.step >= 4 && cell.step % 4 == 0 ? std::vector<std::int64_t>{1, 2, 3} : std::vector<std::int64_t>{2};
    for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        const std::array<int, 3> direction{dx, dy, dz};
        std::vector<int> tangents;
        for (int axis = 0; axis < 3; ++axis) if (direction[static_cast<std::size_t>(axis)] == 0) tangents.push_back(axis);
        std::array<std::int64_t, 3> selections{2, 2, 2};
        const std::size_t combinations = static_cast<std::size_t>(std::pow(static_cast<double>(sampleFractions.size()), static_cast<double>(tangents.size())));
        for (std::size_t combination = 0; combination < combinations; ++combination) {
            std::size_t remainder = combination;
            for (const int axis : tangents) {
                selections[static_cast<std::size_t>(axis)] = sampleFractions[remainder % sampleFractions.size()];
                remainder /= sampleFractions.size();
            }
            LogicalPoint probe{};
            for (std::size_t axis = 0; axis < 3; ++axis) {
                if (direction[axis] < 0) probe[axis] = cell.lower[axis] - 1;
                else if (direction[axis] > 0) probe[axis] = cell.lower[axis] + cell.step + 1;
                else probe[axis] = cell.lower[axis] + cell.step * selections[axis] / 4;
            }
            const auto neighbor = tree.locate(probe);
            if (neighbor >= 0) found.insert(neighbor);
        }
    }
    return {found.begin(), found.end()};
}

std::array<LogicalPoint, 4> positiveFaceCorners(const OctreeNode& cell, int axis, int side,
                                                const std::array<std::int64_t, 3>& low,
                                                const std::array<std::int64_t, 3>& high) {
    const auto plane = side > 0 ? cell.lower[static_cast<std::size_t>(axis)] + cell.step : cell.lower[static_cast<std::size_t>(axis)];
    if (axis == 0) return {{{plane, low[1], low[2]}, {plane, high[1], low[2]}, {plane, high[1], high[2]}, {plane, low[1], high[2]}}};
    if (axis == 1) return {{{low[0], plane, low[2]}, {low[0], plane, high[2]}, {high[0], plane, high[2]}, {high[0], plane, low[2]}}};
    return {{{low[0], low[1], plane}, {high[0], low[1], plane}, {high[0], high[1], plane}, {low[0], high[1], plane}}};
}

}  // namespace

HexMesh makeUniformCartesianMesh(const Json& project) {
    requireExactKeys(project, {"schema_version", "units", "cad", "fluid_domain", "background", "sizing",
        "boundary_layers", "quality", "resources", "output"}, "project");
    if (project.value("schema_version", 0) != 1) throw std::runtime_error("unsupported project schema_version");
    if (project.at("quality").contains("max_surface_deviation_m")) {
        const auto& value = project.at("quality").at("max_surface_deviation_m");
        if (!value.is_number() || !std::isfinite(value.get<double>()) || value.get<double>() <= 0.0)
            throw std::runtime_error("quality.max_surface_deviation_m must be finite and positive (meters)");
    }
    requireExactKeys(project.at("units"), {"length", "output_length"}, "units");
    if (project.at("units").value("output_length", std::string{}) != "m") throw std::runtime_error("output_length must be meters");
    const double scale = lengthScale(project.at("units").at("length").get<std::string>());
    const double maximumCellSize = project.at("sizing").at("maximum_size").get<double>() * scale;
    if (!std::isfinite(maximumCellSize) || maximumCellSize <= 0.0) throw std::runtime_error("sizing.maximum_size must be finite and positive");

    requireExactKeys(project.at("background"), {"bounds", "cell_counts", "uniform_cells"}, "background");
    const Json& background = project.at("background");
    if (!background.value("uniform_cells", true)) throw std::runtime_error("M2 supports only a uniform Cartesian background grid");
    requireExactKeys(background.at("bounds"), {"minimum", "maximum"}, "background.bounds");
    if (!background.at("cell_counts").is_array() || background.at("cell_counts").size() != 3U) {
        throw std::runtime_error("background.cell_counts must contain exactly three positive integers");
    }

    HexMesh mesh;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto rawCount = background.at("cell_counts")[axis];
        if (!rawCount.is_number_integer() || rawCount.get<std::int64_t>() <= 0 ||
            rawCount.get<std::int64_t>() > std::numeric_limits<std::int32_t>::max()) {
            throw std::runtime_error("background.cell_counts entries must fit positive 32-bit integers");
        }
        mesh.cell_counts[axis] = rawCount.get<std::int32_t>();
        const double lower = background.at("bounds").at("minimum").at(axis).get<double>() * scale;
        const double upper = background.at("bounds").at("maximum").at(axis).get<double>() * scale;
        if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower) {
            throw std::runtime_error("background bounds must be finite and have positive extent on each axis");
        }
        mesh.bounds_m[axis] = lower;
        mesh.bounds_m[axis + 3U] = upper;
    }
    const auto maxCells = project.at("resources").at("max_cells").get<std::uint64_t>();
    std::uint64_t cellCount = 1U;
    for (const auto count : mesh.cell_counts) {
        const auto factor = static_cast<std::uint64_t>(count);
        const std::uint64_t limit = std::min(maxCells, static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()));
        if (factor == 0U || cellCount > limit / factor) {
            throw std::runtime_error("uniform Cartesian cell count exceeds resources.max_cells or BabelSim's 32-bit index limit");
        }
        cellCount *= factor;
    }
    const std::uint64_t nxp = static_cast<std::uint64_t>(mesh.cell_counts[0]) + 1U;
    const std::uint64_t nyp = static_cast<std::uint64_t>(mesh.cell_counts[1]) + 1U;
    const std::uint64_t nzp = static_cast<std::uint64_t>(mesh.cell_counts[2]) + 1U;
    const auto pointLimit = static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max());
    if (nxp > pointLimit / nyp || nxp * nyp > pointLimit / nzp) {
        throw std::runtime_error("uniform Cartesian point count exceeds BabelSim's 32-bit index limit");
    }
    const std::uint64_t pointCount = nxp * nyp * nzp;
    const long double nx = mesh.cell_counts[0], ny = mesh.cell_counts[1], nz = mesh.cell_counts[2];
    const long double faceCountEstimate = (nx + 1.0L) * ny * nz + nx * (ny + 1.0L) * nz + nx * ny * (nz + 1.0L);
    if (faceCountEstimate > static_cast<long double>(pointLimit)) {
        throw std::runtime_error("uniform Cartesian face count exceeds BabelSim's 32-bit index limit");
    }
    const std::uint64_t maxMemoryMb = project.at("resources").at("max_memory_mb").get<std::uint64_t>();
    const long double memoryBudget = static_cast<long double>(maxMemoryMb) * 1024.0L * 1024.0L;
    if (estimateCartesianStorageBytes(mesh) > memoryBudget * 0.75L) {
        throw std::runtime_error("estimated uniform Cartesian storage exceeds 75% of resources.max_memory_mb");
    }
    const std::string domainMode = project.at("fluid_domain").at("mode").get<std::string>();
    if (domainMode == "external") {
        const auto expectedOuter = project.at("fluid_domain").at("outer_bounds");
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double domainLower = expectedOuter.at("minimum").at(axis).get<double>() * scale;
            const double domainUpper = expectedOuter.at("maximum").at(axis).get<double>() * scale;
            const double tolerance = std::max(1.0, std::abs(mesh.bounds_m[axis + 3U])) * 1.0e-12;
            if (std::abs(domainLower - mesh.bounds_m[axis]) > tolerance ||
                std::abs(domainUpper - mesh.bounds_m[axis + 3U]) > tolerance) {
                throw std::runtime_error("external-domain Cartesian background bounds must equal fluid_domain.outer_bounds");
            }
        }
    } else if (domainMode != "internal") {
        throw std::runtime_error("fluid_domain.mode must be 'external' or 'internal'");
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double cellSize = (mesh.bounds_m[axis + 3U] - mesh.bounds_m[axis]) / mesh.cell_counts[axis];
        if (cellSize > maximumCellSize * (1.0 + 1.0e-12)) {
            throw std::runtime_error("background.cell_counts produce a cell larger than sizing.maximum_size");
        }
    }
    if (project.at("boundary_layers").at("enabled").get<bool>()) throw std::runtime_error("boundary layers are added in M5 and are not available in M2");
    if (!project.at("sizing").at("sources").empty()) throw std::runtime_error("regional sizing sources are added in M3 and are not available in M2");

    mesh.points_m.reserve(static_cast<std::size_t>(pointCount));
    for (std::int32_t k = 0; k <= mesh.cell_counts[2]; ++k) {
        for (std::int32_t j = 0; j <= mesh.cell_counts[1]; ++j) {
            for (std::int32_t i = 0; i <= mesh.cell_counts[0]; ++i) {
                std::array<double, 3> point{};
                const std::array<std::int32_t, 3> ijk{i, j, k};
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    const double fraction = static_cast<double>(ijk[axis]) / static_cast<double>(mesh.cell_counts[axis]);
                    point[axis] = mesh.bounds_m[axis] + fraction * (mesh.bounds_m[axis + 3U] - mesh.bounds_m[axis]);
                }
                mesh.points_m.push_back(point);
            }
        }
    }
    mesh.cells.reserve(static_cast<std::size_t>(cellCount));
    const double initialDx = (mesh.bounds_m[3] - mesh.bounds_m[0]) / mesh.cell_counts[0];
    const double initialDy = (mesh.bounds_m[4] - mesh.bounds_m[1]) / mesh.cell_counts[1];
    const double initialDz = (mesh.bounds_m[5] - mesh.bounds_m[2]) / mesh.cell_counts[2];
    const double initialCellVolume = initialDx * initialDy * initialDz;
    const double initialAspect = std::max({initialDx, initialDy, initialDz}) / std::min({initialDx, initialDy, initialDz});
    for (std::int32_t k = 0; k < mesh.cell_counts[2]; ++k) for (std::int32_t j = 0; j < mesh.cell_counts[1]; ++j) for (std::int32_t i = 0; i < mesh.cell_counts[0]; ++i) {
        mesh.cells.push_back({pointId(mesh, i, j, k), pointId(mesh, i + 1, j, k), pointId(mesh, i + 1, j + 1, k), pointId(mesh, i, j + 1, k),
            pointId(mesh, i, j, k + 1), pointId(mesh, i + 1, j, k + 1), pointId(mesh, i + 1, j + 1, k + 1), pointId(mesh, i, j + 1, k + 1)});
        mesh.cell_levels.push_back(0);
        mesh.strict_hex.push_back(true);
        mesh.cell_volumes_m3.push_back(initialCellVolume);
        mesh.cell_aspect_ratios.push_back(initialAspect);
    }
    mesh.patches = {{"xmin", "generic"}, {"xmax", "generic"}, {"ymin", "generic"},
        {"ymax", "generic"}, {"zmin", "generic"}, {"zmax", "generic"}};
    makeXFaces(mesh);
    makeYFaces(mesh);
    makeZFaces(mesh);
    return mesh;
}

HexMesh makeAdaptiveCartesianMesh(const Json& project, const CadModel* cadModel, bool refineFluidBoundary) {
    if (refineFluidBoundary && cadModel == nullptr) throw std::runtime_error("CAD surface refinement requires a validated CAD model");
    if (project.at("sizing").at("sources").empty() && !refineFluidBoundary) return makeUniformCartesianMesh(project);
    Json baseProject = project;
    baseProject["sizing"]["sources"] = Json::array();
    HexMesh mesh = makeUniformCartesianMesh(baseProject);
    const auto rootCounts = mesh.cell_counts;
    const auto bounds = mesh.bounds_m;
    mesh = HexMesh{};
    mesh.adaptive = true;
    mesh.cad_surface_refinement_applied = refineFluidBoundary;
    mesh.cell_counts = rootCounts;
    mesh.bounds_m = bounds;

    const std::int32_t maxLevel = project.at("sizing").at("max_level").get<std::int32_t>();
    if (maxLevel < 0 || maxLevel > 20) throw std::runtime_error("sizing.max_level must be in [0, 20]");
    const std::int64_t rootStep = std::int64_t{1} << (maxLevel + 1);
    const double scale = lengthScale(project.at("units").at("length").get<std::string>());
    const auto sources = parseSizingSources(project, scale, cadModel);
    if (sources.empty() && !refineFluidBoundary) return makeUniformCartesianMesh(project);
    const double minimumSize = project.at("sizing").at("minimum_size").get<double>() * scale;
    if (refineFluidBoundary && !project.at("sizing").contains("surface_size")) {
        throw std::runtime_error("CAD boundary meshing requires an explicit sizing.surface_size; it caps the longest physical edge of cells in the boundary band, independent of CAD fit tolerance");
    }
    const double surfaceSizeInput = project.at("sizing").value("surface_size", project.at("sizing").at("maximum_size").get<double>());
    const double surfaceSize = surfaceSizeInput * scale;
    if (refineFluidBoundary && (!std::isfinite(surfaceSize) || surfaceSize < minimumSize || surfaceSize <= 0.0)) {
        throw std::runtime_error("sizing.surface_size must be finite and at least sizing.minimum_size (the hard lower bound on cell size)");
    }
    mesh.cad_surface_target_size_m = refineFluidBoundary ? surfaceSize : 0.0;
    mesh.cad_surface_buffer_layers = project.at("sizing").value("surface_buffer_layers",2);
    if(mesh.cad_surface_buffer_layers<0 || mesh.cad_surface_buffer_layers>8)
        throw std::runtime_error("sizing.surface_buffer_layers must be in [0,8]");
    const std::uint64_t maxCells = project.at("resources").at("max_cells").get<std::uint64_t>();
    const std::uint64_t maxMemoryMb = project.at("resources").at("max_memory_mb").get<std::uint64_t>();
    const long double memoryBudget = static_cast<long double>(maxMemoryMb) * 1024.0L * 1024.0L * 0.75L;
    OctreeForest tree(rootCounts, rootStep);
    auto enforceBudget = [&]() {
        if (tree.leafCount() > maxCells || tree.leafCount() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("adaptive refinement exceeds resources.max_cells or BabelSim's 32-bit cell limit");
        }
        const long double estimated = static_cast<long double>(tree.leafCount()) * 512.0L +
            static_cast<long double>(tree.nodeCount()) * 192.0L;
        if (estimated > memoryBudget) throw std::runtime_error("adaptive refinement exceeds 75% of resources.max_memory_mb");
    };
    enforceBudget();

    std::vector<const SizingSource*> prioritized;
    prioritized.reserve(sources.size());
    for (const auto& source : sources) prioritized.push_back(&source);
    std::stable_sort(prioritized.begin(), prioritized.end(), [](const auto* a, const auto* b) { return a->targetSize < b->targetSize; });
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<std::int32_t> refine;
        for (const auto id : tree.leaves()) {
            const auto& cell = tree.node(id);
            const auto cellBounds = physicalBounds(cell, mesh, rootStep);
            const double cellSize = std::max({cellBounds[3] - cellBounds[0], cellBounds[4] - cellBounds[1], cellBounds[5] - cellBounds[2]});
            for (const auto* source : prioritized) {
                if (!sourceIntersects(*source, cellBounds) || cellSize <= source->targetSize * (1.0 + 1.0e-12)) continue;
                if (0.5 * cellSize < minimumSize * (1.0 - 1.0e-12)) {
                    throw std::runtime_error("sizing source '" + source->id + "' target_size would require cells below sizing.minimum_size");
                }
                if (cell.level >= maxLevel) {
                    throw std::runtime_error("sizing source '" + source->id + "' cannot meet target_size before sizing.max_level");
                }
                refine.insert(id);
                break;
            }
        }
        for (const auto id : refine) { tree.split(id); changed = true; enforceBudget(); }
    }

    if (refineFluidBoundary) {
        bool surfaceChanged = true;
        while (surfaceChanged) {
            surfaceChanged = false;
            std::set<std::int32_t> refine;
            for (const auto id : tree.leaves()) {
                const auto& cell = tree.node(id);
                const auto cellBounds = physicalBounds(cell, mesh, rootStep);
                const double cellEdge = std::max({cellBounds[3] - cellBounds[0], cellBounds[4] - cellBounds[1], cellBounds[5] - cellBounds[2]});
                const bool meetsSurfaceSize = cellEdge <= surfaceSize * (1.0 + 1.0e-12);
                if (meetsSurfaceSize) continue;
                bool surfaceHit = false;
                auto queryBounds=cellBounds;
                for(std::size_t axis=0;axis<3;++axis){
                    queryBounds[axis]-=surfaceSize*mesh.cad_surface_buffer_layers;
                    queryBounds[axis+3]+=surfaceSize*mesh.cad_surface_buffer_layers;
                }
                if (cadModel->domainMode() == "external") {
                    surfaceHit = surfaceIntersectsBoxExceptBoundary(cadModel->fluidSurfaceMesh(),
                        queryBounds, mesh.bounds_m, 1.0e-9 * scale);
                } else {
                    const auto classification = classifyFluidCell(*cadModel, queryBounds, 1.0e-9 * scale);
                    surfaceHit = classification.relation == CellRelation::intersect ||
                        classification.relation == CellRelation::uncertain;
                }
                if (!surfaceHit) continue;
                if (0.5 * cellEdge < minimumSize * (1.0 - 1.0e-12)) {
                    throw std::runtime_error("sizing.surface_size would require cells below sizing.minimum_size; relax one of these conflicting constraints");
                }
                if (cell.level >= maxLevel) {
                    throw std::runtime_error("CAD boundary sizing.surface_size cannot be met before sizing.max_level; increase max_level or coarsen surface_size");
                }
                refine.insert(id);
            }
            for (const auto id : refine) { tree.split(id); surfaceChanged = true; enforceBudget(); }
        }
    }

    bool balanced = false;
    while (!balanced) {
        balanced = true;
        std::set<std::int32_t> refine;
        const auto leaves = tree.leaves();
        for (const auto id : leaves) {
            const auto& cell = tree.node(id);
            for (const auto neighborId : adjacentLeaves(tree, cell)) {
                const auto& neighbor = tree.node(neighborId);
                if (cell.level + 1 < neighbor.level) {
                    const auto cellBounds = physicalBounds(cell, mesh, rootStep);
                    const double cellSize = std::max({cellBounds[3]-cellBounds[0],cellBounds[4]-cellBounds[1],cellBounds[5]-cellBounds[2]});
                    if (0.5 * cellSize < minimumSize * (1.0 - 1.0e-12))
                        throw std::runtime_error("2:1 balance would violate sizing.minimum_size");
                    refine.insert(id);
                } else if (neighbor.level + 1 < cell.level) {
                    const auto neighborBounds = physicalBounds(neighbor, mesh, rootStep);
                    const double neighborSize = std::max({neighborBounds[3]-neighborBounds[0],neighborBounds[4]-neighborBounds[1],neighborBounds[5]-neighborBounds[2]});
                    if (0.5 * neighborSize < minimumSize * (1.0 - 1.0e-12))
                        throw std::runtime_error("2:1 balance would violate sizing.minimum_size");
                    refine.insert(neighborId);
                }
            }
        }
        for (const auto id : refine) {
            if (tree.node(id).level >= maxLevel) throw std::runtime_error("2:1 balance requires refinement beyond sizing.max_level");
            tree.split(id);
            balanced = false;
            enforceBudget();
        }
    }

    std::vector<std::int32_t> leaves = tree.leaves();
    std::sort(leaves.begin(), leaves.end(), [&](std::int32_t a, std::int32_t b) {
        const auto& x = tree.node(a); const auto& y = tree.node(b);
        return std::tie(x.lower[2], x.lower[1], x.lower[0], x.level) < std::tie(y.lower[2], y.lower[1], y.lower[0], y.level);
    });
    std::vector<std::int32_t> nodeToCell(tree.nodeCount(), -1);
    std::vector<OctreeNode> cells;
    cells.reserve(leaves.size());
    for (const auto id : leaves) {
        nodeToCell[static_cast<std::size_t>(id)] = static_cast<std::int32_t>(cells.size());
        cells.push_back(tree.node(id));
    }
    for (const auto& source : sources) {
        std::size_t covered = 0;
        std::vector<std::uint8_t> markers;
        markers.reserve(cells.size());
        for (const auto& cell : cells) {
            const bool intersects = sourceIntersects(source, physicalBounds(cell, mesh, rootStep));
            markers.push_back(intersects ? 1U : 0U);
            if (intersects) ++covered;
        }
        mesh.source_intersection_counts[source.id] = covered;
        mesh.sizing_source_markers.emplace(source.id, std::move(markers));
        Json details{{"kind", source.kind}, {"target_size_m", source.targetSize}, {"marked_cell_count", covered}};
        if (!source.cadGroup.empty()) details["cad_group"] = source.cadGroup;
        if (!source.otherCadGroup.empty()) details["other_cad_group"] = source.otherCadGroup;
        if (source.kind == "distance" || source.kind == "gap" || source.kind == "curvature" || source.kind == "cad_group") {
            details["distance_band_m"] = source.distance;
        }
        if (source.kind == "curvature") {
            details["curvature_sampling"] = "3x3 UV samples restricted to trimmed CAD faces";
            details["curvature_max_radius_m"] = source.curvatureMaxRadius;
            details["curvature_size_criterion"] = "cell target_size is independent of max radius selector";
            details["maximum_selected_curvature_radius_m"] = source.targetSize;
        }
        mesh.sizing_source_details.emplace(source.id, std::move(details));
        if (covered == 0U) throw std::runtime_error("sizing source '" + source.id + "' intersects no background cells");
    }

    std::map<LogicalPoint, std::int32_t> pointIds;
    auto pointFor = [&](const LogicalPoint& point) {
        const auto found = pointIds.find(point);
        if (found != pointIds.end()) return found->second;
        if (mesh.points_m.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("adaptive point count exceeds BabelSim's 32-bit index limit");
        }
        const auto physical = logicalToPhysical(point, mesh, rootStep);
        const auto id = static_cast<std::int32_t>(mesh.points_m.size());
        mesh.points_m.push_back(physical);
        pointIds.emplace(point, id);
        return id;
    };
    mesh.cells.reserve(cells.size());
    mesh.strict_hex.assign(cells.size(), true);
    mesh.cell_levels.reserve(cells.size());
    mesh.cell_volumes_m3.reserve(cells.size());
    mesh.cell_aspect_ratios.reserve(cells.size());
    for (std::size_t ci = 0; ci < cells.size(); ++ci) {
        const auto& cell = cells[ci];
        const auto hi = LogicalPoint{cell.lower[0] + cell.step, cell.lower[1] + cell.step, cell.lower[2] + cell.step};
        const std::array<LogicalPoint, 8> corners{{cell.lower, {hi[0], cell.lower[1], cell.lower[2]},
            {hi[0], hi[1], cell.lower[2]}, {cell.lower[0], hi[1], cell.lower[2]},
            {cell.lower[0], cell.lower[1], hi[2]}, {hi[0], cell.lower[1], hi[2]}, hi, {cell.lower[0], hi[1], hi[2]}}};
        std::array<std::int32_t, 8> vertexIds{};
        for (std::size_t v = 0; v < corners.size(); ++v) vertexIds[v] = pointFor(corners[v]);
        mesh.cells.push_back(vertexIds);
        mesh.cell_levels.push_back(cell.level);
        const auto physical = physicalBounds(cell, mesh, rootStep);
        const double dx = physical[3] - physical[0], dy = physical[4] - physical[1], dz = physical[5] - physical[2];
        mesh.cell_volumes_m3.push_back(dx * dy * dz);
        mesh.cell_aspect_ratios.push_back(std::max({dx, dy, dz}) / std::min({dx, dy, dz}));
    }

    mesh.patches = {{"xmin", "generic"}, {"xmax", "generic"}, {"ymin", "generic"},
        {"ymax", "generic"}, {"zmin", "generic"}, {"zmax", "generic"}};
    struct FaceTile {
        std::array<LogicalPoint, 4> corners{};
        std::vector<LogicalPoint> ring;
        std::int32_t owner{};
        std::int32_t neighbour{-1};
        std::int32_t patch{-1};
    };
    std::vector<FaceTile> tiles;
    const std::array<std::int64_t, 3> domainHigh{static_cast<std::int64_t>(rootCounts[0]) * rootStep,
        static_cast<std::int64_t>(rootCounts[1]) * rootStep, static_cast<std::int64_t>(rootCounts[2]) * rootStep};
    for (std::size_t ci = 0; ci < cells.size(); ++ci) {
        const auto& cell = cells[ci];
        for (int axis = 0; axis < 3; ++axis) for (const int side : {-1, 1}) {
            const auto neighborNodes = touchingLeaves(tree, cell, axis, side, true);
            if (neighborNodes.size() > 1U) mesh.strict_hex[ci] = false;
            if (neighborNodes.empty()) {
                const bool boundary = side < 0 ? cell.lower[static_cast<std::size_t>(axis)] == 0 :
                    cell.lower[static_cast<std::size_t>(axis)] + cell.step == domainHigh[static_cast<std::size_t>(axis)];
                if (!boundary) throw std::runtime_error("adaptive Cartesian topology contains an uncovered interior face");
                if (side > 0) {
                    std::array<std::int64_t, 3> lo = cell.lower;
                    std::array<std::int64_t, 3> hi{cell.lower[0] + cell.step, cell.lower[1] + cell.step, cell.lower[2] + cell.step};
                    FaceTile tile;
                    tile.corners = positiveFaceCorners(cell, axis, side, lo, hi);
                    tile.owner = static_cast<std::int32_t>(ci);
                    tile.patch = 2 * axis + 1;
                    tiles.push_back(std::move(tile));
                } else {
                    std::array<std::int64_t, 3> lo = cell.lower;
                    std::array<std::int64_t, 3> hi{cell.lower[0] + cell.step, cell.lower[1] + cell.step, cell.lower[2] + cell.step};
                    FaceTile tile;
                    tile.corners = positiveFaceCorners(cell, axis, side, lo, hi);
                    std::reverse(tile.corners.begin(), tile.corners.end());
                    tile.owner = static_cast<std::int32_t>(ci);
                    tile.patch = 2 * axis;
                    tiles.push_back(std::move(tile));
                }
                continue;
            }
            if (side < 0) continue;
            for (const auto neighborNode : neighborNodes) {
                const auto& neighbor = tree.node(neighborNode);
                std::array<std::int64_t, 3> lo{}, hi{};
                bool positiveArea = true;
                for (int d = 0; d < 3; ++d) {
                    if (d == axis) { lo[static_cast<std::size_t>(d)] = cell.lower[static_cast<std::size_t>(d)] + cell.step; hi[static_cast<std::size_t>(d)] = lo[static_cast<std::size_t>(d)]; }
                    else {
                        lo[static_cast<std::size_t>(d)] = std::max(cell.lower[static_cast<std::size_t>(d)], neighbor.lower[static_cast<std::size_t>(d)]);
                        hi[static_cast<std::size_t>(d)] = std::min(cell.lower[static_cast<std::size_t>(d)] + cell.step, neighbor.lower[static_cast<std::size_t>(d)] + neighbor.step);
                        positiveArea = positiveArea && hi[static_cast<std::size_t>(d)] > lo[static_cast<std::size_t>(d)];
                    }
                }
                if (!positiveArea) continue;
                const auto neighborCell = nodeToCell.at(static_cast<std::size_t>(neighborNode));
                if (neighborCell < 0) throw std::runtime_error("adaptive face references a non-leaf neighbor");
                FaceTile tile;
                tile.corners = positiveFaceCorners(cell, axis, side, lo, hi);
                tile.owner = static_cast<std::int32_t>(ci);
                tile.neighbour = neighborCell;
                tiles.push_back(std::move(tile));
            }
        }
    }

    using LineKey = std::tuple<int, std::int64_t, std::int64_t>;
    auto lineKey = [](const LogicalPoint& point, int axis) -> LineKey {
        if (axis == 0) return {axis, point[1], point[2]};
        if (axis == 1) return {axis, point[0], point[2]};
        return {axis, point[0], point[1]};
    };
    std::map<LineKey, std::set<std::int64_t>> lineVertices;
    for (const auto& tile : tiles) for (const auto& point : tile.corners) for (int axis = 0; axis < 3; ++axis) {
        lineVertices[lineKey(point, axis)].insert(point[static_cast<std::size_t>(axis)]);
    }
    for (auto& tile : tiles) {
        tile.ring.clear();
        for (std::size_t edge = 0; edge < 4U; ++edge) {
            const auto& a = tile.corners[edge];
            const auto& b = tile.corners[(edge + 1U) % 4U];
            int axis = -1;
            for (int d = 0; d < 3; ++d) if (a[static_cast<std::size_t>(d)] != b[static_cast<std::size_t>(d)]) axis = d;
            if (axis < 0) throw std::runtime_error("adaptive face contains a zero-length edge");
            const auto& candidates = lineVertices.at(lineKey(a, axis));
            std::vector<std::int64_t> ordered;
            for (const auto value : candidates) if (value >= std::min(a[static_cast<std::size_t>(axis)], b[static_cast<std::size_t>(axis)]) &&
                value <= std::max(a[static_cast<std::size_t>(axis)], b[static_cast<std::size_t>(axis)])) ordered.push_back(value);
            if (a[static_cast<std::size_t>(axis)] > b[static_cast<std::size_t>(axis)]) std::reverse(ordered.begin(), ordered.end());
            if (ordered.size() < 2U || ordered.front() != a[static_cast<std::size_t>(axis)] || ordered.back() != b[static_cast<std::size_t>(axis)]) {
                throw std::runtime_error("adaptive face edge overlay lost an endpoint");
            }
            for (std::size_t q = 0; q + 1U < ordered.size(); ++q) {
                LogicalPoint point = a;
                point[static_cast<std::size_t>(axis)] = ordered[q];
                tile.ring.push_back(point);
            }
        }
    }
    for (const auto& tile : tiles) {
        const bool subdividedEdge = tile.ring.size() != 4U;
        if (subdividedEdge) {
            mesh.strict_hex[static_cast<std::size_t>(tile.owner)] = false;
            if (tile.neighbour >= 0) mesh.strict_hex[static_cast<std::size_t>(tile.neighbour)] = false;
        }
        if (!subdividedEdge) {
            std::vector<std::int32_t> ring;
            for (const auto& point : tile.ring) ring.push_back(pointFor(point));
            mesh.faces.push_back({std::move(ring), tile.owner, tile.neighbour, tile.patch});
        } else {
            LogicalPoint center{};
            for (const auto& point : tile.corners) for (std::size_t axis = 0; axis < 3; ++axis) center[axis] += point[axis];
            for (auto& coordinate : center) coordinate /= 4;
            const auto centerId = pointFor(center);
            for (std::size_t i = 0; i < tile.ring.size(); ++i) {
                const auto next = (i + 1U) % tile.ring.size();
                mesh.faces.push_back({{centerId, pointFor(tile.ring[i]), pointFor(tile.ring[next])},
                    tile.owner, tile.neighbour, tile.patch});
            }
        }
    }
    if (refineFluidBoundary) {
        mesh = extractSeedConnectedFluid(std::move(mesh), project, *cadModel);
        const double dx = mesh.bounds_m[3] - mesh.bounds_m[0];
        const double dy = mesh.bounds_m[4] - mesh.bounds_m[1];
        const double dz = mesh.bounds_m[5] - mesh.bounds_m[2];
        const double domainDiagonal = std::hypot(dx, dy, dz);
        const double maximumDeviation = project.at("quality").value("max_surface_deviation_m",
            std::max(domainDiagonal * 1.0e-4, 1.0e-9));
        mesh.cad_surface_deviation_overridden = project.at("quality").contains("max_surface_deviation_m");
        const bool strictQuality = project.at("quality").value("mode", std::string("strict")) == "strict";
        snapFluidBoundaryToCad(mesh, *cadModel, maximumDeviation, strictQuality,
            project.at("resources").at("max_cells").get<std::size_t>());
        if(strictQuality){
            const auto audit=auditSurfaceDistance(mesh,*cadModel,surfaceSize,2000000U);
            if(std::max(audit.mesh_to_cad_max_m,audit.cad_to_mesh_max_m)>maximumDeviation)
                throw std::runtime_error("CAD bidirectional sampled surface error exceeds geometry tolerance: mesh-to-CAD="+
                    std::to_string(audit.mesh_to_cad_max_m)+" m, CAD-to-mesh="+
                    std::to_string(audit.cad_to_mesh_max_m)+" m, allowed="+std::to_string(maximumDeviation)+" m");
            mesh.cad_bidirectional_surface_audit=true;
            mesh.cad_mesh_to_surface_max_m=audit.mesh_to_cad_max_m;
            mesh.cad_surface_to_mesh_max_m=audit.cad_to_mesh_max_m;
            mesh.cad_surface_audit_sample_spacing_m=audit.requested_sample_spacing_m;
            mesh.cad_mesh_to_surface_samples=audit.mesh_to_cad_samples;
            mesh.cad_surface_to_mesh_samples=audit.cad_to_mesh_samples;
        }
        return mesh;
    }
    return mesh;
}

Json meshReport(const HexMesh& mesh, const std::filesystem::path& projectFile) {
    const double cellVolume = mesh.cell_volumes_m3.empty() ?
        ((mesh.bounds_m[3] - mesh.bounds_m[0]) / mesh.cell_counts[0]) *
            ((mesh.bounds_m[4] - mesh.bounds_m[1]) / mesh.cell_counts[1]) *
            ((mesh.bounds_m[5] - mesh.bounds_m[2]) / mesh.cell_counts[2]) : mesh.cell_volumes_m3.front();
    const auto volumeRange = mesh.cell_volumes_m3.empty() ? std::pair<double, double>{cellVolume, cellVolume} :
        std::pair<double, double>{*std::min_element(mesh.cell_volumes_m3.begin(), mesh.cell_volumes_m3.end()),
                                  *std::max_element(mesh.cell_volumes_m3.begin(), mesh.cell_volumes_m3.end())};
    const double maximumAspect = mesh.cell_aspect_ratios.empty() ? 1.0 :
        *std::max_element(mesh.cell_aspect_ratios.begin(), mesh.cell_aspect_ratios.end());
    const auto strictHexCount = mesh.strict_hex.empty() ? mesh.cells.size() :
        static_cast<std::size_t>(std::count(mesh.strict_hex.begin(), mesh.strict_hex.end(), true));
    Json report;
    report["schema_version"] = 1;
    report["project"] = projectFile.string();
    report["generator"] = mesh.adaptive ? "adaptive-cartesian-octree" : "uniform-cartesian-hex";
    report["cad_clipping_applied"] = mesh.cad_clipping_applied;
    report["cad_removed_cells"] = mesh.cad_removed_cells;
    report["cad_intersect_cells"] = mesh.cad_intersect_cells;
    report["cad_uncertain_cells"] = mesh.cad_uncertain_cells;
    report["cad_seed_forced_cells"] = mesh.cad_seed_forced_cells;
    report["cad_selected_components"] = mesh.cad_selected_components;
    report["cad_surface_refinement_applied"] = mesh.cad_surface_refinement_applied;
    report["cad_surface_target_size_m"] = mesh.cad_surface_target_size_m > 0.0 ? Json(mesh.cad_surface_target_size_m) : Json(nullptr);
    report["cad_surface_buffer_layers"] = mesh.cad_surface_buffer_layers;
    report["cad_snapping_applied"] = mesh.cad_snapping_applied;
    report["cad_snap_fraction"] = mesh.cad_surface_refinement_applied ? Json(mesh.cad_snap_fraction) : Json(nullptr);
    report["cad_snapping_vertices"] = mesh.cad_snapping_vertices;
    report["cad_projection_faces"] = mesh.cad_projection_faces;
    report["cad_collapsed_faces_removed_during_feature_canonicalization"] = mesh.cad_collapsed_faces_removed;
    report["cad_collapsed_face_edge_repairs"] = mesh.cad_collapsed_face_edge_repairs;
    report["cad_consolidated_faces"] = mesh.cad_consolidated_faces;
    report["cad_quality_rollback_count"] = mesh.cad_quality_rollback_count;
    report["cad_deformed_vertices"] = mesh.cad_deformed_vertices;
    report["cad_deformation_iterations"] = mesh.cad_deformation_iterations;
    report["cad_tangential_smoothing_iterations"] = mesh.cad_tangential_smoothing_iterations;
    report["cad_tangential_smoothing_max_move_m"] = mesh.cad_surface_refinement_applied ?
        Json(mesh.cad_tangential_smoothing_max_move_m) : Json(nullptr);
    report["cad_quality_optimization_iterations"] = mesh.cad_quality_optimization_iterations;
    report["cad_quality_optimization_attempts"] = mesh.cad_quality_optimization_attempts;
    report["cad_quality_optimization_proposed_vertices"] = mesh.cad_quality_optimization_proposed_vertices;
    report["cad_quality_optimization_objective_before"] = mesh.cad_surface_refinement_applied ?
        Json(mesh.cad_quality_optimization_objective_before) : Json(nullptr);
    report["cad_quality_optimization_objective_after"] = mesh.cad_surface_refinement_applied ?
        Json(mesh.cad_quality_optimization_objective_after) : Json(nullptr);
    report["cad_max_deformation_m"] = mesh.cad_surface_refinement_applied ? Json(mesh.cad_max_deformation_m) : Json(nullptr);
    report["cad_core_deformation_fraction"] = mesh.cad_surface_refinement_applied ? Json(mesh.cad_core_deformation_fraction) : Json(nullptr);
    report["cad_surface_deviation_m"] = mesh.cad_surface_refinement_applied ? Json(mesh.cad_surface_deviation_m) : Json(nullptr);
    report["cad_surface_deviation_tolerance_m"] = mesh.cad_surface_refinement_applied ? Json(mesh.cad_surface_deviation_tolerance_m) : Json(nullptr);
    report["cad_surface_deviation_overridden"] = mesh.cad_surface_deviation_overridden;
    report["cad_surface_deviation_samples"] = mesh.cad_surface_deviation_samples;
    report["cad_bidirectional_surface_audit"] = {{"completed",mesh.cad_bidirectional_surface_audit},
        {"mesh_to_cad_max_m",mesh.cad_bidirectional_surface_audit?Json(mesh.cad_mesh_to_surface_max_m):Json(nullptr)},
        {"cad_to_mesh_max_m",mesh.cad_bidirectional_surface_audit?Json(mesh.cad_surface_to_mesh_max_m):Json(nullptr)},
        {"sample_spacing_m",mesh.cad_bidirectional_surface_audit?Json(mesh.cad_surface_audit_sample_spacing_m):Json(nullptr)},
        {"mesh_to_cad_samples",mesh.cad_mesh_to_surface_samples},{"cad_to_mesh_samples",mesh.cad_surface_to_mesh_samples},
        {"continuous_hausdorff_bound",false}};
    report["cad_fluid_domain_volume_m3"] = mesh.cad_fluid_domain_volume_m3 > 0.0 ? Json(mesh.cad_fluid_domain_volume_m3) : Json(nullptr);
    report["cad_mesh_volume_relative_error"] = mesh.cad_clipping_applied ? Json(mesh.cad_mesh_volume_relative_error) : Json(nullptr);
    report["limitations"] = mesh.cad_clipping_applied ?
        Json::array({mesh.cad_snapping_applied ?
            "CAD boundary snapping passed two-way sampled CAD distance and positive-volume checks; see face and cell quality metrics" :
            "CAD extraction is a staircase candidate; the CAD snap was rolled back in permissive mode or no CAD boundary vertices required movement"}) :
        (mesh.adaptive ?
            Json::array({"M3 refines the Cartesian background only; CAD clipping is scheduled for M4"}) :
            Json::array({"M2 meshes the Cartesian background bounds only; CAD obstacle subtraction is scheduled for M4"}));
    report["vertices"] = mesh.points_m.size();
    report["faces"] = mesh.faces.size();
    report["cells"] = mesh.cells.size();
    report["strict_hex_cells"] = strictHexCount;
    report["polyhedron_cells"] = mesh.cells.size() - strictHexCount;
    const double domainVolume = (mesh.bounds_m[3] - mesh.bounds_m[0]) * (mesh.bounds_m[4] - mesh.bounds_m[1]) * (mesh.bounds_m[5] - mesh.bounds_m[2]);
    report["total_volume_m3"] = mesh.cell_volumes_m3.empty() ? domainVolume :
        std::accumulate(mesh.cell_volumes_m3.begin(), mesh.cell_volumes_m3.end(), 0.0);
    report["maximum_level"] = mesh.cell_levels.empty() ? 0 : *std::max_element(mesh.cell_levels.begin(), mesh.cell_levels.end());
    report["boundary_faces"] = std::count_if(mesh.faces.begin(), mesh.faces.end(), [](const MeshFace& face) { return face.neighbour < 0; });
    report["cell_counts"] = mesh.cell_counts;
    report["bounds_m"] = {{mesh.bounds_m[0], mesh.bounds_m[1], mesh.bounds_m[2]},
        {mesh.bounds_m[3], mesh.bounds_m[4], mesh.bounds_m[5]}};
    const bool constantCellVolume = std::abs(volumeRange.second - volumeRange.first) <= std::max(volumeRange.second, 1.0e-300) * 1.0e-12;
    report["cell_volume_m3"] = constantCellVolume ? Json(cellVolume) : Json(nullptr);
    const long double estimatedBytes = static_cast<long double>(mesh.points_m.size()) * sizeof(std::array<double, 3>) +
        static_cast<long double>(mesh.cells.size()) * sizeof(std::array<std::int32_t, 8>) +
        static_cast<long double>(mesh.faces.size()) * (sizeof(MeshFace) + 8.0L * sizeof(std::int32_t));
    report["estimated_mesh_storage_bytes"] = static_cast<std::uint64_t>(estimatedBytes);
    const auto qualityAudit = auditMeshQuality(mesh);
    auto weightedP99 = [](const std::vector<double>& values, const std::vector<double>& weights) {
        std::vector<std::pair<double,double>> sorted;
        double total=0.0;
        for(std::size_t i=0;i<values.size() && i<weights.size();++i) {
            sorted.emplace_back(values[i],weights[i]); total+=weights[i];
        }
        if(sorted.empty() || !(total>0)) return 0.0;
        std::sort(sorted.begin(),sorted.end(),[](const auto& a,const auto& b){return a.first<b.first;});
        double cumulative=0;
        for(const auto& entry:sorted) { cumulative+=entry.second; if(cumulative>=.99*total) return entry.first; }
        return sorted.back().first;
    };
    report["quality"] = {{"min_cell_volume_m3", qualityAudit.minimum_cell_volume},
        {"max_cell_volume_m3", volumeRange.second},
        {"min_scaled_jacobian", qualityAudit.regular_hex_cells ? Json(qualityAudit.minimum_regular_hex_scaled_jacobian) : Json(nullptr)},
        {"regular_hex_cells_evaluated", qualityAudit.regular_hex_cells},
        {"max_non_orthogonality_degrees", qualityAudit.maximum_non_orthogonality_degrees},
        {"p99_non_orthogonality_degrees_face_count", qualityAudit.p99_non_orthogonality_degrees},
        {"p99_non_orthogonality_degrees_volume_weighted", weightedP99(qualityAudit.non_orthogonality,qualityAudit.non_orthogonality_weights)},
        {"max_internal_skewness", qualityAudit.maximum_internal_skewness},
        {"p99_internal_skewness_face_count", qualityAudit.p99_internal_skewness},
        {"p99_internal_skewness_volume_weighted", weightedP99(qualityAudit.internal_skewness,qualityAudit.internal_skewness_weights)},
        {"max_cell_aspect_ratio", mesh.cad_surface_refinement_applied ? Json(nullptr) : Json(maximumAspect)},
        {"negative_or_zero_volume_cells", qualityAudit.negative_or_zero_volume_cells},
        {"degenerate_faces", qualityAudit.degenerate_faces}, {"full_cfd_profile_evaluated", false}};
    report["patches"] = mesh.patches;
    report["sizing_source_intersection_counts"] = mesh.source_intersection_counts;
    report["sizing_sources"] = mesh.sizing_source_details;
    return report;
}

void writeBabelSimV3(const HexMesh& mesh, const std::filesystem::path& file) {
    atomicWrite(file, [&](std::ostream& output) {
        output << "BABELSIM_MESH 3\nvertices " << mesh.points_m.size() << '\n' << std::setprecision(17);
        for (const auto& point : mesh.points_m) output << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
        output << "faces " << mesh.faces.size() << '\n';
        for (const auto& face : mesh.faces) {
            output << "face " << face.vertices.size() << ' ' << face.owner << ' ' << face.neighbour << ' ' << face.patch;
            for (const auto vertex : face.vertices) output << ' ' << vertex;
            output << '\n';
        }
        output << "patches " << mesh.patches.size() << '\n';
        for (const auto& patch : mesh.patches) output << "patch " << patch.first << ' ' << patch.second << '\n';
        output << "end\n";
    });
}

void writeLegacyVtu(const HexMesh& mesh, const std::filesystem::path& file) {
    atomicWrite(file, [&](std::ostream& output) {
        std::vector<std::vector<std::size_t>> facesByCell(mesh.cells.size());
        for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
            const auto& item = mesh.faces[face];
            facesByCell.at(static_cast<std::size_t>(item.owner)).push_back(face);
            if (item.neighbour >= 0) facesByCell.at(static_cast<std::size_t>(item.neighbour)).push_back(face);
        }
        std::vector<std::int64_t> connectivity, offsets, faceStream, faceOffsets, types;
        connectivity.reserve(mesh.cells.size() * 8U);
        offsets.reserve(mesh.cells.size());
        types.reserve(mesh.cells.size());
        for (std::size_t cellIdValue = 0; cellIdValue < mesh.cells.size(); ++cellIdValue) {
            const bool strictHex = mesh.strict_hex.empty() || mesh.strict_hex.at(cellIdValue);
            std::set<std::int32_t> cellVertices;
            if (strictHex) {
                for (const auto vertex : mesh.cells[cellIdValue]) cellVertices.insert(vertex);
                for (const auto vertex : mesh.cells[cellIdValue]) connectivity.push_back(vertex);
            } else {
                for (const auto faceId : facesByCell[cellIdValue]) {
                    for (const auto vertex : mesh.faces[faceId].vertices) cellVertices.insert(vertex);
                }
                for (const auto vertex : cellVertices) connectivity.push_back(vertex);
            }
            offsets.push_back(static_cast<std::int64_t>(connectivity.size()));
            types.push_back(strictHex ? 12 : 42);
            if (strictHex) {
                faceOffsets.push_back(-1);
            } else {
                faceStream.push_back(static_cast<std::int64_t>(facesByCell[cellIdValue].size()));
                for (const auto faceId : facesByCell[cellIdValue]) {
                    const auto& face = mesh.faces[faceId];
                    std::vector<std::int32_t> ring = face.vertices;
                    if (face.neighbour == static_cast<std::int32_t>(cellIdValue)) std::reverse(ring.begin(), ring.end());
                    faceStream.push_back(static_cast<std::int64_t>(ring.size()));
                    for (const auto vertex : ring) faceStream.push_back(vertex);
                }
                faceOffsets.push_back(static_cast<std::int64_t>(faceStream.size()));
            }
        }
        output << std::setprecision(17)
            << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "<UnstructuredGrid><Piece NumberOfPoints=\"" << mesh.points_m.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">\n"
            << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& point : mesh.points_m) output << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
        output << "</DataArray></Points><Cells><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto item : connectivity) output << item << ' ';
        output << "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        for (const auto item : offsets) output << item << ' ';
        output << "\n</DataArray><DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (const auto item : types) output << item << ' ';
        output << "\n</DataArray><DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";
        for (const auto item : faceStream) output << item << ' ';
        output << "\n</DataArray><DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n";
        for (const auto item : faceOffsets) output << item << ' ';
        output << "\n</DataArray></Cells><CellData Scalars=\"cell_type\">\n"
               << "<DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
        for (const auto level : mesh.cell_levels) output << level << ' ';
        output << "\n</DataArray><DataArray type=\"UInt8\" Name=\"cell_type\" format=\"ascii\">\n";
        for (const auto type : types) output << type << ' ';
        output << "\n</DataArray><DataArray type=\"UInt8\" Name=\"strict_hex\" format=\"ascii\">\n";
        for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
            output << ((mesh.strict_hex.empty() || mesh.strict_hex[cell]) ? 1 : 0) << ' ';
        }
        output << "\n</DataArray>";
        for (const auto& source : mesh.sizing_source_markers) {
            output << "<DataArray type=\"UInt8\" Name=\"source_" << source.first << "\" format=\"ascii\">\n";
            for (const auto marker : source.second) output << static_cast<unsigned int>(marker) << ' ';
            output << "\n</DataArray>";
        }
        output << "</CellData></Piece></UnstructuredGrid></VTKFile>\n";
    });
}

int runGenerate(int argc, char** argv) {
    try {
    if (argc < 2) throw std::runtime_error("usage: babelmesh generate PROJECT.json [--out DIRECTORY] [--cad-staircase]");
    const std::filesystem::path projectFile = std::filesystem::absolute(argv[1]);
    std::filesystem::path outputOverride;
    bool cadStaircase = false;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--out" && i + 1 < argc) outputOverride = argv[++i];
        else if (std::string(argv[i]) == "--cad-staircase") cadStaircase = true;
        else throw std::runtime_error("unknown generate option '" + std::string(argv[i]) + "'");
    }
    std::ifstream input(projectFile);
    if (!input) throw std::runtime_error("cannot open project file: " + projectFile.string());
    Json project;
    input >> project;
    const CadModel cadModel = buildCadModel(projectFile);  // Validate all project sections and CAD geometry before generating any output.
    const HexMesh mesh = cadStaircase ? makeAdaptiveCartesianMesh(project, &cadModel, true) :
        (project.at("sizing").at("sources").empty() ?
            makeUniformCartesianMesh(project) : makeAdaptiveCartesianMesh(project, &cadModel));
    const auto requestedDirectory = outputOverride.empty()
        ? std::filesystem::path(project.at("output").at("directory").get<std::string>()) : outputOverride;
    const auto outputDirectory = requestedDirectory.is_absolute() ? requestedDirectory : projectFile.parent_path() / requestedDirectory;
    const bool overwrite = project.at("output").value("overwrite", false);
    const auto formats = project.at("output").at("formats");
    const auto meshFile = outputDirectory / "mesh.bsm";
    const auto vtuFile = outputDirectory / "mesh.vtu";
    const auto reportFile = outputDirectory / "mesh_report.json";
    if (!overwrite && (std::filesystem::exists(meshFile) || std::filesystem::exists(vtuFile) || std::filesystem::exists(reportFile))) {
        throw std::runtime_error("output already exists; set output.overwrite=true or choose another directory");
    }
    for (const auto& format : formats) {
        const std::string name = format.get<std::string>();
        if (name == "babelsim-v3") writeBabelSimV3(mesh, meshFile);
        else if (name == "vtu") writeLegacyVtu(mesh, vtuFile);
        else throw std::runtime_error("M2 does not support requested output format '" + name + "'");
    }
    atomicWrite(reportFile, [&](std::ostream& output) { output << meshReport(mesh, projectFile).dump(2) << '\n'; });
    const auto strictHexCount = mesh.strict_hex.empty() ? mesh.cells.size() :
        static_cast<std::size_t>(std::count(mesh.strict_hex.begin(), mesh.strict_hex.end(), true));
    std::cout << "generated " << (mesh.adaptive ? "adaptive" : "uniform") << " Cartesian mesh: " << mesh.cells.size()
              << " cells (" << strictHexCount << " strict hex, " << mesh.cells.size() - strictHexCount << " polyhedral transitions), "
              << mesh.faces.size() << " faces, " << mesh.points_m.size() << " vertices\n"
              << (mesh.cad_clipping_applied ?
                  (mesh.cad_snapping_applied ? "seed-connected CAD extraction and feature-aware snapping passed sampled bidirectional surface-distance and mesh-quality checks; see mesh_report.json for audit coverage\n" :
                      "seed-connected CAD staircase candidate generated without accepted boundary snapping\n") :
                  "CAD clipping is not yet applied; use --cad-staircase for the current M4 staircase path\n") << "output: "
              << outputDirectory.string() << '\n';
    return 0;
    } catch (const std::exception& error) {
        std::cerr << "generate failed: " << error.what() << '\n';
        return 1;
    }
}

}  // namespace babelmesh
