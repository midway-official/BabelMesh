#include "babelmesh/mixed_mesh.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace babelmesh {
namespace {
using Vec = Eigen::Vector3d;
Vec point(const MixedMesh &m, std::int32_t id) {
  if (id < 0 || static_cast<std::size_t>(id) >= m.points_m.size())
    throw std::runtime_error("mixed mesh vertex index is out of range");
  const auto &p = m.points_m[static_cast<std::size_t>(id)];
  return {p[0], p[1], p[2]};
}
std::array<double, 3> array(const Vec &p) { return {p.x(), p.y(), p.z()}; }
Vec centre(const MixedMesh &m, const FaceKey &ring) {
  Vec c = Vec::Zero();
  for (auto v : ring)
    c += point(m, v);
  return c / static_cast<double>(ring.size());
}
Vec area(const MixedMesh &m, const FaceKey &ring) {
  const Vec origin = point(m, ring.front());
  Vec a = Vec::Zero();
  for (std::size_t j = 1; j + 1 < ring.size(); ++j)
    a += (point(m, ring[j]) - origin).cross(point(m, ring[j + 1]) - origin) *
         0.5;
  return a;
}
FaceKey canonicalRing(FaceKey ring) {
  std::rotate(ring.begin(), std::min_element(ring.begin(), ring.end()),
              ring.end());
  return ring;
}
FaceKey solverCompatibleRing(const MixedMesh &mesh, FaceKey ring) {
  // BabelSim triangulates a polygon from its first vertex and rejects zero
  // fan triangles. A quadtree mortar polygon can include collinear midpoints;
  // choose a cyclic origin whose fan is nondegenerate without deleting them.
  double best = -1;
  std::size_t origin = 0;
  for (std::size_t start = 0; start < ring.size(); ++start) {
    const Vec a = point(mesh, ring[start]);
    double minimum = std::numeric_limits<double>::max();
    for (std::size_t j = 1; j + 1 < ring.size(); ++j) {
      const Vec u = point(mesh, ring[(start + j) % ring.size()]) - a,
                v = point(mesh, ring[(start + j + 1) % ring.size()]) - a;
      minimum = std::min(
          minimum, u.cross(v).squaredNorm() /
                       std::max(u.squaredNorm() * v.squaredNorm(), 1e-300));
    }
    if (minimum > best) {
      best = minimum;
      origin = start;
    }
  }
  if (best <= 1e-24)
    throw std::runtime_error("polygon cannot form a nondegenerate vertex fan "
                             "for the BabelSim reader");
  std::rotate(ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(origin),
              ring.end());
  return ring;
}
bool opposite(const FaceKey &a, const FaceKey &b) {
  if (a.size() != b.size())
    return false;
  const auto it = std::find(b.begin(), b.end(), a[0]);
  if (it == b.end())
    return false;
  const auto base = static_cast<std::size_t>(it - b.begin());
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[(base + b.size() - i) % b.size()])
      return false;
  return true;
}
Vec faceCentroid(const MixedMesh &m, const FaceKey &ring) {
  const Vec origin = point(m, ring[0]);
  Vec result = Vec::Zero();
  double weight = 0;
  for (std::size_t j = 1; j + 1 < ring.size(); ++j) {
    const Vec a = point(m, ring[j]), b = point(m, ring[j + 1]);
    const double area = (a - origin).cross(b - origin).norm();
    result += area * (origin + a + b) / 3;
    weight += area;
  }
  if (!(weight > 0))
    throw std::runtime_error("zero face area");
  return result / weight;
}
double percentile(std::vector<double> v, double q) {
  if (v.empty())
    return 0.0;
  const auto i =
      static_cast<std::size_t>(std::ceil(q * static_cast<double>(v.size()))) -
      1;
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(i),
                   v.end());
  return v[i];
}
void checkStream(const std::ofstream &out, const std::filesystem::path &p) {
  if (!out)
    throw std::runtime_error("cannot write mesh artifact: " + p.string());
}
} // namespace

FaceKey faceKey(FaceKey v) {
  std::sort(v.begin(), v.end());
  return v;
}
const char *cellTypeName(CellType type) {
  switch (type) {
  case CellType::tetrahedron:
    return "tetrahedron";
  case CellType::hexahedron:
    return "hexahedron";
  case CellType::prism:
    return "prism";
  case CellType::pyramid:
    return "pyramid";
  case CellType::polyhedron:
    return "polyhedron";
  }
  throw std::runtime_error("unsupported mixed cell type");
}
std::vector<FaceKey> cellFaceRings(const MixedCell &c) {
  if (c.type == CellType::polyhedron) {
    if (c.poly_faces.size() < 4 || c.vertices.size() < 4)
      throw std::runtime_error("invalid explicit polyhedron");
    const std::set<std::int32_t> vertices(c.vertices.begin(), c.vertices.end());
    if (vertices.size() != c.vertices.size())
      throw std::runtime_error("duplicate polyhedron vertex");
    for (const auto &f : c.poly_faces) {
      if (f.size() < 3 ||
          std::set<std::int32_t>(f.begin(), f.end()).size() != f.size())
        throw std::runtime_error("invalid polyhedron face");
      for (auto v : f)
        if (!vertices.count(v))
          throw std::runtime_error(
              "polyhedron face vertex is missing from its cell");
    }
    return c.poly_faces;
  }
  std::vector<std::vector<std::size_t>> local;
  std::size_t expected = 0;
  switch (c.type) {
  case CellType::tetrahedron:
    expected = 4;
    local = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    break;
  case CellType::hexahedron:
    expected = 8;
    local = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
             {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    break;
  case CellType::prism:
    expected = 6;
    local = {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    break;
  case CellType::pyramid:
    expected = 5;
    local = {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
    break;
  default:
    throw std::runtime_error("unknown mixed cell type");
  }
  if (c.vertices.size() != expected ||
      std::set<std::int32_t>(c.vertices.begin(), c.vertices.end()).size() !=
          expected)
    throw std::runtime_error("invalid " + std::string(cellTypeName(c.type)) +
                             " connectivity");
  std::vector<FaceKey> result;
  for (const auto &face : local) {
    FaceKey ids;
    for (auto i : face)
      ids.push_back(c.vertices[i]);
    result.push_back(std::move(ids));
  }
  return result;
}

void buildMixedFaces(MixedMesh &mesh, const BoundaryMap &boundary) {
  mesh.faces.clear();
  std::map<FaceKey, std::size_t> lookup;
  for (std::size_t c = 0; c < mesh.cells.size(); ++c) {
    for (auto ring : cellFaceRings(mesh.cells[c])) {
      const auto key = faceKey(ring);
      auto inserted = lookup.emplace(key, mesh.faces.size());
      if (inserted.second)
        mesh.faces.push_back({solverCompatibleRing(mesh, std::move(ring)),
                              static_cast<std::int32_t>(c), -1, -1});
      else {
        auto &f = mesh.faces[inserted.first->second];
        if (f.neighbour >= 0)
          throw std::runtime_error("non-manifold shared face");
        if (!opposite(f.vertices, ring))
          throw std::runtime_error(
              "overlapping cells or inconsistent shared-face orientation");
        f.neighbour = static_cast<std::int32_t>(c);
      }
    }
  }
  std::size_t declarations = 0;
  for (auto &f : mesh.faces) {
    auto it = boundary.find(faceKey(f.vertices));
    if (f.neighbour < 0) {
      if (it == boundary.end())
        throw std::runtime_error("unassigned boundary face; cavity and "
                                 "core/layer interface do not conform");
      f.patch = it->second;
      ++declarations;
    } else if (it != boundary.end())
      throw std::runtime_error("declared boundary face is internal");
  }
  if (declarations != boundary.size())
    throw std::runtime_error(
        "boundary declaration is missing from the final mesh");
  validateMixedMesh(mesh);
}

nlohmann::json validateMixedMeshBounds(const MixedMesh &mesh) {
  using Json = nlohmann::json;
  if (!mesh.metadata.contains("domain_bounds_constraint"))
    return {{"available", false}};
  const auto &bounds = mesh.metadata.at("domain_bounds_constraint");
  const auto low = bounds.at("minimum_m").get<std::array<double, 3>>();
  const auto high = bounds.at("maximum_m").get<std::array<double, 3>>();
  auto actualLow = high, actualHigh = low;
  for (std::size_t axis = 0; axis < 3; ++axis)
    if (!std::isfinite(low[axis]) || !std::isfinite(high[axis]) ||
        low[axis] >= high[axis])
      throw std::runtime_error("invalid locked domain bounds");
  for (std::size_t node = 0; node < mesh.points_m.size(); ++node)
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double x = mesh.points_m[node][axis];
      if (!std::isfinite(x) || x < low[axis] || x > high[axis])
        throw std::runtime_error(
            "mesh point " + std::to_string(node) +
            " exceeds locked background/domain bounds on axis " +
            std::to_string(axis));
      actualLow[axis] = std::min(actualLow[axis], x);
      actualHigh[axis] = std::max(actualHigh[axis], x);
    }
  return Json{{"available", true},
              {"passed", true},
              {"checked_points", mesh.points_m.size()},
              {"minimum_m", actualLow},
              {"maximum_m", actualHigh},
              {"maximum_excess_m", 0.0},
              {"method", "all vertices within fixed closed coordinate bounds; "
                         "convex cells retain the same envelope"}};
}
void validateMixedMesh(MixedMesh &mesh) {
  mesh.metadata["bounds_audit"] = validateMixedMeshBounds(mesh);
  constexpr auto limit =
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
  if (mesh.cells.empty() || mesh.points_m.empty() || mesh.faces.empty() ||
      mesh.patches.empty() || mesh.cells.size() >= limit ||
      mesh.points_m.size() >= limit || mesh.faces.size() >= limit)
    throw std::runtime_error(
        "empty mixed mesh or 32-bit indexing limit exceeded");
  if (mesh.generation_dimension != 2 && mesh.generation_dimension != 3)
    throw std::runtime_error("invalid generation dimension");
  if (mesh.generation_dimension == 2 &&
      (mesh.extrusion_axis < 0 || mesh.extrusion_axis > 2 ||
       !std::isfinite(mesh.extrusion_thickness_m) ||
       mesh.extrusion_thickness_m <= 0))
    throw std::runtime_error("invalid planar extrusion metadata");
  for (const auto &cell : mesh.cells) {
    if (static_cast<int>(cell.region) > 2)
      throw std::runtime_error("invalid cell region");
    if (mesh.generation_dimension == 2 && cell.type != CellType::hexahedron &&
        cell.type != CellType::prism && cell.type != CellType::polyhedron)
      throw std::runtime_error(
          "planar generation allows only extruded hexahedra, triangular prisms "
          "and polygonal prisms; tetrahedra and pyramids are forbidden");
  }
  if (mesh.generation_dimension == 2) {
    const auto axis = static_cast<std::size_t>(mesh.extrusion_axis);
    double low = std::numeric_limits<double>::max(), high = -low;
    for (const auto &p : mesh.points_m) {
      low = std::min(low, p[axis]);
      high = std::max(high, p[axis]);
    }
    const double tolerance =
        std::max(mesh.extrusion_thickness_m * 1e-8,
                 std::max(std::abs(low), std::abs(high)) * 1e-12);
    if (std::abs(high - low - mesh.extrusion_thickness_m) > tolerance)
      throw std::runtime_error(
          "planar mesh span differs from its extrusion thickness");
    for (const auto &p : mesh.points_m)
      if (std::min(std::abs(p[axis] - low), std::abs(p[axis] - high)) >
          tolerance)
        throw std::runtime_error(
            "planar mesh must have exactly one extrusion layer");
  }
  std::set<std::string> patchNames;
  const std::set<std::string> kinds = {"generic", "wall",     "inlet",
                                       "outlet",  "symmetry", "processor"};
  for (const auto &p : mesh.patches) {
    if (p.name.empty() ||
        p.name.find_first_of(" \t\r\n") != std::string::npos ||
        !patchNames.insert(p.name).second || !kinds.count(p.kind))
      throw std::runtime_error("invalid or duplicate mixed mesh patch");
  }
  for (const auto &p : mesh.points_m)
    for (double v : p)
      if (!std::isfinite(v))
        throw std::runtime_error("non-finite mixed mesh coordinate");
  std::vector<std::map<FaceKey, FaceKey>> expected(mesh.cells.size()),
      actual(mesh.cells.size());
  for (std::size_t c = 0; c < mesh.cells.size(); ++c)
    for (const auto &ring : cellFaceRings(mesh.cells[c])) {
      for (auto v : ring)
        (void)point(mesh, v);
      if (!expected[c].emplace(faceKey(ring), canonicalRing(ring)).second)
        throw std::runtime_error("duplicate declared cell face");
    }
  std::set<FaceKey> unique;
  for (const auto &f : mesh.faces) {
    if (f.owner < 0 || static_cast<std::size_t>(f.owner) >= mesh.cells.size() ||
        f.neighbour == f.owner || f.neighbour < -1 ||
        (f.neighbour >= 0 &&
         static_cast<std::size_t>(f.neighbour) >= mesh.cells.size()))
      throw std::runtime_error("invalid face owner/neighbour");
    if ((f.neighbour >= 0 && f.patch != -1) ||
        (f.neighbour < 0 && (f.patch < 0 || static_cast<std::size_t>(f.patch) >=
                                                mesh.patches.size())))
      throw std::runtime_error("invalid internal/boundary patch assignment");
    if (f.vertices.size() < 3 || !unique.insert(faceKey(f.vertices)).second)
      throw std::runtime_error("duplicate or degenerate mesh face");
    actual[static_cast<std::size_t>(f.owner)].emplace(
        faceKey(f.vertices), canonicalRing(f.vertices));
    if (f.neighbour >= 0) {
      auto reverse = f.vertices;
      std::reverse(reverse.begin(), reverse.end());
      actual[static_cast<std::size_t>(f.neighbour)].emplace(
          faceKey(reverse), canonicalRing(reverse));
    }
  }
  mesh.volumes_m3.assign(mesh.cells.size(), 0.0);
  mesh.centres_m.resize(mesh.cells.size());
  for (std::size_t c = 0; c < mesh.cells.size(); ++c) {
    if (actual[c] != expected[c])
      throw std::runtime_error(
          "cell shell differs from its declared connectivity at cell " +
          std::to_string(c));
    const Vec ref = centre(mesh, mesh.cells[c].vertices);
    double volume = 0.0;
    Vec moment = Vec::Zero(), closure = Vec::Zero();
    double surfaceArea = 0.0;
    std::map<std::pair<std::int32_t, std::int32_t>, int> edgeBalance,
        edgeCounts;
    for (const auto &entry : actual[c]) {
      const auto &ring = entry.second;
      const Vec av = area(mesh, ring);
      const double a = av.norm();
      if (!(a > 0.0) || !std::isfinite(a))
        throw std::runtime_error("degenerate face at cell " +
                                 std::to_string(c));
      closure += av;
      surfaceArea += a;
      if (av.dot(centre(mesh, ring) - ref) <= 0.0)
        throw std::runtime_error("inverted or non-convex cell " +
                                 std::to_string(c));
      const Vec p0 = point(mesh, ring[0]);
      for (std::size_t j = 1; j + 1 < ring.size(); ++j) {
        const Vec p1 = point(mesh, ring[j]), p2 = point(mesh, ring[j + 1]);
        const double v = (p0 - ref).dot((p1 - ref).cross(p2 - ref)) / 6.0;
        if (v < -1e-12 * std::pow((p0 - ref).norm() + (p1 - ref).norm() +
                                      (p2 - ref).norm(),
                                  3))
          throw std::runtime_error("negative subvolume in cell " +
                                   std::to_string(c));
        if (v <= 0.0)
          continue; // Collinear vertices retain conforming edge connectivity.
        volume += v;
        moment += v * (ref + p0 + p1 + p2) / 4.0;
      }
      for (std::size_t j = 0; j < ring.size(); ++j) {
        const auto aId = ring[j], bId = ring[(j + 1) % ring.size()];
        edgeBalance[std::minmax(aId, bId)] += aId < bId ? 1 : -1;
        ++edgeCounts[std::minmax(aId, bId)];
      }
    }
    if (closure.norm() > surfaceArea * 1e-9 ||
        std::any_of(edgeCounts.begin(), edgeCounts.end(),
                    [](const auto &e) { return e.second != 2; }) ||
        std::any_of(edgeBalance.begin(), edgeBalance.end(),
                    [](const auto &e) { return e.second != 0; }))
      throw std::runtime_error("open cell shell at cell " + std::to_string(c));
    if (!(volume > 0.0) || !std::isfinite(volume))
      throw std::runtime_error("non-positive cell volume");
    mesh.volumes_m3[c] = volume;
    mesh.centres_m[c] = array(moment / volume);
  }
}

nlohmann::json mixedMeshReport(const MixedMesh &mesh) {
  using Json = nlohmann::json;
  Json counts = {{"hexahedron", 0},
                 {"prism", 0},
                 {"tetrahedron", 0},
                 {"pyramid", 0},
                 {"polyhedron", 0}};
  Json regions = {{"cartesian", 0}, {"boundary_layer", 0}, {"transition", 0}};
  constexpr const char *regionNames[] = {"cartesian", "boundary_layer",
                                         "transition"};
  for (const auto &c : mesh.cells) {
    const auto name = cellTypeName(c.type);
    counts[name] = counts[name].get<std::size_t>() + 1;
    const auto r = static_cast<std::size_t>(c.region);
    if (r >= 3)
      throw std::runtime_error("invalid cell region");
    regions[regionNames[r]] = regions[regionNames[r]].get<std::size_t>() + 1;
  }
  std::vector<double> nonorth, skew;
  std::array<std::vector<double>, 3> regionAngles, regionSkews, meanRatios;
  for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
    const auto &cell = mesh.cells[i];
    if (cell.type != CellType::tetrahedron)
      continue;
    double squaredEdges = 0;
    for (std::size_t a = 0; a < 4; ++a)
      for (std::size_t b = a + 1; b < 4; ++b)
        squaredEdges +=
            (point(mesh, cell.vertices[a]) - point(mesh, cell.vertices[b]))
                .squaredNorm();
    const double meanRatio =
        12 * std::pow(3 * mesh.volumes_m3[i], 2.0 / 3) / squaredEdges;
    meanRatios[static_cast<std::size_t>(cell.region)].push_back(meanRatio);
  }
  std::size_t worst = 0, worstSkew = 0, boundary = 0;
  double maximumSkew = -1;
  double maximum = 0;
  for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
    const auto &f = mesh.faces[i];
    const Vec av = area(mesh, f.vertices), fc = faceCentroid(mesh, f.vertices);
    const auto &oc = mesh.centres_m[static_cast<std::size_t>(f.owner)];
    const Vec owner(oc[0], oc[1], oc[2]);
    Vec d = fc - owner;
    if (f.neighbour >= 0) {
      const auto &nc = mesh.centres_m[static_cast<std::size_t>(f.neighbour)];
      d = Vec(nc[0], nc[1], nc[2]) - owner;
      const double projection = av.dot(d);
      if (projection > 0.0) {
        const Vec hit = owner + d * (av.dot(fc - owner) / projection);
        skew.push_back((fc - hit).norm() / std::sqrt(av.norm()));
      } else
        skew.push_back(std::numeric_limits<double>::max());
      if (skew.back() > maximumSkew) {
        maximumSkew = skew.back();
        worstSkew = i;
      }
    } else
      ++boundary;
    const double cosine =
        std::clamp(av.dot(d) / (av.norm() * d.norm()), -1.0, 1.0);
    const double angle = std::acos(cosine) * 180.0 / 3.14159265358979323846;
    nonorth.push_back(angle);
    const auto ownerRegion = static_cast<std::size_t>(
        mesh.cells[static_cast<std::size_t>(f.owner)].region);
    regionAngles[ownerRegion].push_back(angle);
    if (f.neighbour >= 0) {
      regionSkews[ownerRegion].push_back(skew.back());
      const auto neighbourRegion = static_cast<std::size_t>(
          mesh.cells[static_cast<std::size_t>(f.neighbour)].region);
      if (neighbourRegion != ownerRegion) {
        regionAngles[neighbourRegion].push_back(angle);
        regionSkews[neighbourRegion].push_back(skew.back());
      }
    }
    if (angle > maximum) {
      maximum = angle;
      worst = i;
    }
  }
  Json patchCounts = Json::object();
  for (const auto &f : mesh.faces)
    if (f.patch >= 0) {
      const auto &name = mesh.patches[static_cast<std::size_t>(f.patch)].name;
      patchCounts[name] = patchCounts.value(name, std::size_t{0}) + 1;
    }
  const double maxSkew =
      skew.empty() ? 0 : *std::max_element(skew.begin(), skew.end());
  // A continuous aggregate lets constrained optimization improve bad faces
  // while the maximum and percentile indices stay on unchanged faces.
  double penalty = 0;
  for (double angle : nonorth)
    penalty += std::pow(angle / 60, 4);
  if (!nonorth.empty())
    penalty /= static_cast<double>(nonorth.size());
  double skewPenalty = 0;
  for (double value : skew)
    skewPenalty += std::pow(value / .25, 4);
  if (!skew.empty())
    penalty += skewPenalty / static_cast<double>(skew.size());
  const double p99 = percentile(nonorth, .99), skewP99 = percentile(skew, .99);
  Json regionQuality = Json::object(), shapeQuality = Json::object();
  for (std::size_t i = 0; i < 3; ++i) {
    if (!regionAngles[i].empty())
      regionQuality[regionNames[i]] = {
          {"faces", regionAngles[i].size()},
          {"max_non_orthogonality_degrees",
           *std::max_element(regionAngles[i].begin(), regionAngles[i].end())},
          {"p99_non_orthogonality_degrees", percentile(regionAngles[i], .99)},
          {"max_internal_skewness",
           regionSkews[i].empty() ? 0
                                  : *std::max_element(regionSkews[i].begin(),
                                                      regionSkews[i].end())},
          {"p99_internal_skewness", percentile(regionSkews[i], .99)}};
    if (!meanRatios[i].empty())
      shapeQuality[regionNames[i]] = {
          {"tetrahedra", meanRatios[i].size()},
          {"minimum_mean_ratio",
           *std::min_element(meanRatios[i].begin(), meanRatios[i].end())},
          {"p01_mean_ratio", percentile(meanRatios[i], .01)},
          {"median_mean_ratio", percentile(meanRatios[i], .5)}};
  }
  return {
      {"generator", "cartesian-hybrid"},
      {"generation_dimension", mesh.generation_dimension},
      {"topological_dimension", 3},
      {"extrusion_axis", mesh.extrusion_axis},
      {"extrusion_thickness_m", mesh.extrusion_thickness_m},
      {"cells", mesh.cells.size()},
      {"vertices", mesh.points_m.size()},
      {"faces", mesh.faces.size()},
      {"boundary_faces", boundary},
      {"cell_types", counts},
      {"regions", regions},
      {"patch_faces", patchCounts},
      {"volume_m3",
       std::accumulate(mesh.volumes_m3.begin(), mesh.volumes_m3.end(), 0.0)},
      {"region_face_quality",
       {{"interface_policy", "cross-region faces counted in both regions"},
        {"regions", regionQuality}}},
      {"tetrahedron_shape_quality",
       {{"definition", "12*(3*V)^(2/3)/sum(squared edge lengths); regular=1"},
        {"used_as_fvm_acceptance_gate", false},
        {"regions", shapeQuality}}},
      {"quality",
       {{"skewness_definition", "distance from face centroid to centre-line "
                                "intersection / sqrt(face area)"},
        {"min_cell_volume_m3",
         *std::min_element(mesh.volumes_m3.begin(), mesh.volumes_m3.end())},
        {"negative_or_zero_volume_cells", 0},
        {"open_cells", 0},
        {"unassigned_boundary_faces", 0},
        {"max_non_orthogonality_degrees", maximum},
        {"p99_non_orthogonality_degrees", p99},
        {"max_internal_skewness", maxSkew},
        {"p99_internal_skewness", skewP99},
        {"mean_quartic_face_penalty", penalty},
        {"worst_non_orthogonality_face", worst},
        {"worst_internal_skewness_face", worstSkew},
        {"face_quality_thresholds_passed",
         maximum <= 75 && p99 <= 60 && maxSkew <= .5 && skewP99 <= .25},
        {"full_cfd_profile_evaluated", false}}},
      {"metadata", mesh.metadata}};
}

void writeMixedBabelSimV3(const MixedMesh &mesh,
                          const std::filesystem::path &file) {
  std::ofstream out(file);
  checkStream(out, file);
  out << std::setprecision(17);
  out << "BABELSIM_MESH 3\nvertices " << mesh.points_m.size() << '\n';
  for (const auto &p : mesh.points_m)
    out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
  out << "faces " << mesh.faces.size() << '\n';
  for (const auto &f : mesh.faces) {
    const auto ring = solverCompatibleRing(mesh, f.vertices);
    out << "face " << f.vertices.size() << ' ' << f.owner << ' ' << f.neighbour
        << ' ' << f.patch;
    for (auto v : ring)
      out << ' ' << v;
    out << '\n';
  }
  out << "patches " << mesh.patches.size() << '\n';
  for (const auto &p : mesh.patches)
    out << "patch " << p.name << ' ' << p.kind << '\n';
  out << "end\n";
  out.flush();
  checkStream(out, file);
}
void writeMixedVtu(const MixedMesh &mesh, const std::filesystem::path &file) {
  std::ofstream out(file);
  checkStream(out, file);
  out << std::setprecision(17);
  out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" "
         "version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "<UnstructuredGrid><Piece NumberOfPoints=\"" << mesh.points_m.size()
      << "\" NumberOfCells=\"" << mesh.cells.size() << "\">\n"
      << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";
  for (const auto &p : mesh.points_m)
    out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
  out << "</DataArray></Points><Cells>\n<DataArray type=\"Int32\" "
         "Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto &c : mesh.cells) {
    for (auto v : c.vertices)
      out << v << ' ';
    out << '\n';
  }
  out << "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" "
         "format=\"ascii\">\n";
  std::size_t offset = 0;
  for (const auto &c : mesh.cells) {
    offset += c.vertices.size();
    out << offset << ' ';
  }
  out << "</DataArray><DataArray type=\"UInt8\" Name=\"types\" "
         "format=\"ascii\">\n";
  for (const auto &c : mesh.cells)
    out << static_cast<unsigned int>(c.type) << ' ';
  out << "</DataArray>";
  if (std::any_of(mesh.cells.begin(), mesh.cells.end(), [](const auto &c) {
        return c.type == CellType::polyhedron;
      })) {
    out << "<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";
    std::vector<std::int64_t> offsets;
    std::int64_t cursor = 0;
    for (const auto &cell : mesh.cells) {
      if (cell.type != CellType::polyhedron) {
        offsets.push_back(-1);
        continue;
      }
      out << cell.poly_faces.size() << ' ';
      ++cursor;
      for (const auto &f : cell.poly_faces) {
        out << f.size() << ' ';
        ++cursor;
        for (auto v : f) {
          out << v << ' ';
          ++cursor;
        }
      }
      offsets.push_back(cursor);
    }
    out << "</DataArray><DataArray type=\"Int64\" Name=\"faceoffsets\" "
           "format=\"ascii\">\n";
    for (auto offsetValue : offsets)
      out << offsetValue << ' ';
    out << "</DataArray>";
  }
  out << "</Cells><CellData>\n<DataArray type=\"UInt8\" Name=\"region\" "
         "format=\"ascii\">\n";
  for (const auto &c : mesh.cells)
    out << static_cast<unsigned int>(c.region) << ' ';
  out << "</DataArray><DataArray type=\"Float64\" Name=\"volume_m3\" "
         "format=\"ascii\">\n";
  for (double v : mesh.volumes_m3)
    out << v << ' ';
  out << "</DataArray></CellData></Piece></UnstructuredGrid></VTKFile>\n";
  out.flush();
  checkStream(out, file);
}
void writeMixedCheckpoint(const MixedMesh &m,
                          const std::filesystem::path &file) {
  nlohmann::json j = {{"checkpoint_version", 1},
                      {"generation_dimension", m.generation_dimension},
                      {"extrusion_axis", m.extrusion_axis},
                      {"extrusion_thickness_m", m.extrusion_thickness_m},
                      {"points_m", m.points_m},
                      {"cells", nlohmann::json::array()},
                      {"faces", nlohmann::json::array()},
                      {"patches", nlohmann::json::array()},
                      {"metadata", m.metadata}};
  for (const auto &c : m.cells)
    j["cells"].push_back({{"type", static_cast<int>(c.type)},
                          {"region", static_cast<int>(c.region)},
                          {"vertices", c.vertices},
                          {"poly_faces", c.poly_faces}});
  for (const auto &f : m.faces)
    j["faces"].push_back({{"vertices", f.vertices},
                          {"owner", f.owner},
                          {"neighbour", f.neighbour},
                          {"patch", f.patch}});
  for (const auto &p : m.patches)
    j["patches"].push_back({{"name", p.name}, {"kind", p.kind}});
  std::ofstream out(file);
  checkStream(out, file);
  out << j.dump() << '\n';
  out.flush();
  checkStream(out, file);
}
MixedMesh readMixedCheckpoint(const std::filesystem::path &file) {
  std::ifstream in(file);
  if (!in)
    throw std::runtime_error("cannot read checkpoint: " + file.string());
  nlohmann::json j;
  in >> j;
  if (j.at("checkpoint_version") != 1)
    throw std::runtime_error("unsupported mesh checkpoint version");
  MixedMesh m;
  m.generation_dimension = j.at("generation_dimension");
  m.extrusion_axis = j.at("extrusion_axis");
  m.extrusion_thickness_m = j.at("extrusion_thickness_m");
  m.points_m = j.at("points_m").get<decltype(m.points_m)>();
  m.metadata = j.at("metadata");
  for (const auto &c : j.at("cells"))
    m.cells.push_back({static_cast<CellType>(c.at("type").get<int>()),
                       static_cast<CellRegion>(c.at("region").get<int>()),
                       c.at("vertices").get<FaceKey>(),
                       c.value("poly_faces", std::vector<FaceKey>{})});
  for (const auto &f : j.at("faces"))
    m.faces.push_back({f.at("vertices").get<FaceKey>(), f.at("owner"),
                       f.at("neighbour"), f.at("patch")});
  for (const auto &p : j.at("patches"))
    m.patches.push_back({p.at("name"), p.at("kind")});
  validateMixedMesh(m);
  return m;
}
} // namespace babelmesh
