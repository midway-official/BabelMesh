#include "babelmesh/mesh_optimizer.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace babelmesh {
namespace {
using Vec = Eigen::Vector3d;
using Edge = std::array<std::int32_t, 2>;
using Json = nlohmann::json;
Vec point(const MixedMesh &m, std::int32_t node) {
  const auto &p = m.points_m.at(static_cast<std::size_t>(node));
  return {p[0], p[1], p[2]};
}
struct FaceGeometry {
  double angle{}, skew{}, penalty{};
};
FaceGeometry measure(const MixedMesh &m, const FaceKey &ring, const Vec &owner,
                     const Vec *neighbour) {
  const Vec a = point(m, ring[0]), b = point(m, ring[1]), c = point(m, ring[2]);
  const Vec area = (b - a).cross(c - a) * .5, centre = (a + b + c) / 3;
  const Vec direction = neighbour ? *neighbour - owner : centre - owner;
  const double projection = area.dot(direction);
  if (!(projection > 0) || !(direction.norm() > 0))
    return {180, 1e100, 1e100};
  const double angle =
      std::acos(std::clamp(projection / (area.norm() * direction.norm()), -1.0,
                           1.0)) *
      180 / std::acos(-1.0);
  const double skew =
      neighbour ? (centre - owner -
                   direction * (area.dot(centre - owner) / projection))
                          .norm() /
                      std::sqrt(area.norm())
                : 0;
  return {angle, skew, std::pow(angle / 60, 4) + std::pow(skew / .25, 4)};
}
Vec centre(const MixedMesh &m, std::int32_t cell) {
  const auto &p = m.centres_m.at(static_cast<std::size_t>(cell));
  return {p[0], p[1], p[2]};
}
struct Cavity {
  std::map<FaceKey, FaceKey> shell;
  std::map<FaceKey, Vec> outsideCentres;
};
std::map<FaceKey, FaceKey> shell(const std::vector<MixedCell> &cells) {
  std::map<FaceKey, std::vector<FaceKey>> occurrences;
  for (const auto &cell : cells)
    for (const auto &ring : cellFaceRings(cell))
      occurrences[faceKey(ring)].push_back(ring);
  std::map<FaceKey, FaceKey> result;
  for (const auto &item : occurrences) {
    if (item.second.size() > 2)
      return {};
    if (item.second.size() == 1)
      result.emplace(item.first, item.second.front());
  }
  return result;
}
double localPenalty(const MixedMesh &m, const std::vector<MixedCell> &cells,
                    const Cavity &cavity) {
  struct Entry {
    FaceKey ring;
    Vec owner;
    Vec neighbour;
    bool shared{};
  };
  std::map<FaceKey, Entry> faces;
  for (const auto &cell : cells) {
    Vec centroid = Vec::Zero();
    for (auto node : cell.vertices)
      centroid += point(m, node) / 4;
    for (const auto &ring : cellFaceRings(cell)) {
      auto result =
          faces.emplace(faceKey(ring), Entry{ring, centroid, {}, false});
      if (!result.second) {
        result.first->second.neighbour = centroid;
        result.first->second.shared = true;
      }
    }
  }
  double penalty = 0;
  for (const auto &item : faces) {
    const auto &f = item.second;
    const auto outside = cavity.outsideCentres.find(item.first);
    const Vec *neighbour = f.shared ? &f.neighbour
                           : outside != cavity.outsideCentres.end()
                               ? &outside->second
                               : nullptr;
    penalty += measure(m, f.ring, f.owner, neighbour).penalty;
  }
  return penalty / static_cast<double>(faces.size());
}
MixedCell orientedTet(const MixedMesh &m, FaceKey vertices) {
  const Vec a = point(m, vertices[0]);
  if ((point(m, vertices[1]) - a)
          .cross(point(m, vertices[2]) - a)
          .dot(point(m, vertices[3]) - a) < 0)
    std::swap(vertices[1], vertices[2]);
  return {CellType::tetrahedron, CellRegion::transition, std::move(vertices)};
}
double volume(const MixedMesh &m, const std::vector<MixedCell> &cells) {
  double total = 0;
  for (const auto &cell : cells) {
    const Vec a = point(m, cell.vertices[0]);
    const double v = (point(m, cell.vertices[1]) - a)
                         .cross(point(m, cell.vertices[2]) - a)
                         .dot(point(m, cell.vertices[3]) - a) /
                     6;
    if (!(v > 0))
      return -1;
    total += v;
  }
  return total;
}
bool better(const Json &next, const Json &old) {
  double score = 0, previous = 0;
  for (const auto &item :
       {std::make_pair("max_non_orthogonality_degrees", 75.0),
        std::make_pair("p99_non_orthogonality_degrees", 60.0),
        std::make_pair("max_internal_skewness", .5),
        std::make_pair("p99_internal_skewness", .25)}) {
    const double a = next.at(item.first), b = old.at(item.first);
    if (a > b + 1e-8)
      return false;
    score += std::pow(a / item.second, 4);
    previous += std::pow(b / item.second, 4);
  }
  score += next.at("mean_quartic_face_penalty").get<double>();
  previous += old.at("mean_quartic_face_penalty").get<double>();
  return score < previous * (1 - 1e-7);
}
} // namespace

Json improveTransitionTetTopology(MixedMesh &m, std::size_t maxTrials) {
  // A restricted edge-removal operator following the 3->2 and 4->4 cavity
  // constructions in Klingner/Shewchuk (2007) and Marot/Remacle (2020).
  // Unlike general N-edge removal, these cases never increase the cell count.
  Json quality = mixedMeshReport(m).at("quality");
  std::size_t trials = 0, prefiltered = 0, flips32 = 0, flips44 = 0;
  BoundaryMap boundary;
  for (const auto &face : m.faces)
    if (face.neighbour < 0)
      boundary.emplace(faceKey(face.vertices), face.patch);
  while (trials < maxTrials &&
         !quality.at("face_quality_thresholds_passed").get<bool>()) {
    std::map<Edge, std::vector<std::size_t>> stars;
    std::vector<double> severity(m.cells.size(), 0);
    std::map<FaceKey, const MixedFace *> faceIndex;
    for (const auto &face : m.faces) {
      faceIndex.emplace(faceKey(face.vertices), &face);
      if (face.vertices.size() != 3)
        continue;
      const Vec owner = centre(m, face.owner);
      const Vec neighbour =
          face.neighbour >= 0 ? centre(m, face.neighbour) : Vec{};
      const auto geometry = measure(m, face.vertices, owner,
                                    face.neighbour >= 0 ? &neighbour : nullptr);
      const double bad = std::max(geometry.angle / 60, geometry.skew / .25);
      severity[static_cast<std::size_t>(face.owner)] =
          std::max(severity[static_cast<std::size_t>(face.owner)], bad);
      if (face.neighbour >= 0)
        severity[static_cast<std::size_t>(face.neighbour)] =
            std::max(severity[static_cast<std::size_t>(face.neighbour)], bad);
    }
    // Include every incident cell, so an edge touching a hex, pyramid or layer
    // cannot be misidentified as an entirely removable transition-tet star.
    for (std::size_t i = 0; i < m.cells.size(); ++i)
      for (const auto &ring : cellFaceRings(m.cells[i]))
        for (std::size_t j = 0; j < ring.size(); ++j) {
          Edge edge{ring[j], ring[(j + 1) % ring.size()]};
          if (edge[0] > edge[1])
            std::swap(edge[0], edge[1]);
          auto &list = stars[edge];
          if (list.empty() || list.back() != i)
            list.push_back(i);
        }
    std::vector<std::pair<double, Edge>> order;
    for (const auto &item : stars) {
      if (item.second.size() != 3 && item.second.size() != 4)
        continue;
      double worst = 0;
      bool eligible = true;
      for (auto i : item.second) {
        const auto &cell = m.cells[i];
        eligible &= cell.type == CellType::tetrahedron &&
                    cell.region == CellRegion::transition;
        worst = std::max(worst, severity[i]);
      }
      if (eligible && worst > .8)
        order.emplace_back(worst, item.first);
    }
    std::sort(order.rbegin(), order.rend());
    bool improved = false;
    for (const auto &item : order) {
      const auto &edge = item.second;
      const auto &indices = stars.at(edge);
      std::map<std::int32_t, std::set<std::int32_t>> link;
      std::vector<MixedCell> original;
      for (auto i : indices) {
        original.push_back(m.cells[i]);
        FaceKey other;
        for (auto node : m.cells[i].vertices)
          if (node != edge[0] && node != edge[1])
            other.push_back(node);
        link[other[0]].insert(other[1]);
        link[other[1]].insert(other[0]);
      }
      if (link.size() != indices.size() ||
          std::any_of(link.begin(), link.end(),
                      [](const auto &p) { return p.second.size() != 2; }))
        continue; // An open/boundary star is never removed.
      FaceKey ring{link.begin()->first};
      auto previous = std::int32_t{-1}, current = ring.front();
      for (std::size_t i = 1; i < link.size(); ++i) {
        const auto &adjacent = link.at(current);
        auto next = *adjacent.begin();
        if (next == previous)
          next = *adjacent.rbegin();
        ring.push_back(next);
        previous = current;
        current = next;
      }
      if (std::set<std::int32_t>(ring.begin(), ring.end()).size() !=
              link.size() ||
          !link.at(ring.back()).count(ring.front()))
        continue;
      Cavity cavity;
      cavity.shell = shell(original);
      if (cavity.shell.size() != indices.size() * 2)
        continue;
      for (const auto &f : cavity.shell) {
        const auto &outside = *faceIndex.at(f.first);
        auto cell = outside.neighbour;
        if (!std::binary_search(indices.begin(), indices.end(),
                                static_cast<std::size_t>(outside.owner)))
          cell = outside.owner;
        if (cell >= 0)
          cavity.outsideCentres.emplace(f.first, centre(m, cell));
      }
      const double oldPenalty = localPenalty(m, original, cavity);
      const double oldVolume = volume(m, original);
      const int templates = ring.size() == 3 ? 1 : 2;
      for (int pattern = 0; pattern < templates && trials < maxTrials;
           ++pattern) {
        std::vector<FaceKey> triangles;
        if (ring.size() == 3)
          triangles.push_back(ring);
        else if (pattern == 0)
          triangles = {{ring[0], ring[1], ring[2]},
                       {ring[0], ring[2], ring[3]}};
        else
          triangles = {{ring[0], ring[1], ring[3]},
                       {ring[1], ring[2], ring[3]}};
        std::vector<MixedCell> replacement;
        for (const auto &triangle : triangles)
          for (auto apex : edge) {
            auto vertices = triangle;
            vertices.push_back(apex);
            replacement.push_back(orientedTet(m, std::move(vertices)));
          }
        const auto skin = shell(replacement);
        bool sameSkin = skin.size() == cavity.shell.size();
        for (const auto &f : skin)
          sameSkin &= cavity.shell.count(f.first) != 0;
        const double newVolume = volume(m, replacement);
        if (!sameSkin || newVolume <= 0 ||
            std::abs(newVolume - oldVolume) > oldVolume * 1e-10 ||
            localPenalty(m, replacement, cavity) >= oldPenalty * (1 - 1e-7)) {
          ++prefiltered;
          continue;
        }
        ++trials;
        MixedMesh candidate = m;
        for (std::size_t i = 0; i < replacement.size(); ++i)
          candidate.cells[indices[i]] = replacement[i];
        if (replacement.size() < indices.size())
          candidate.cells.erase(candidate.cells.begin() +
                                static_cast<std::ptrdiff_t>(indices.back()));
        try {
          // Checks shared-face orientation, convex positive cells and locked
          // bounds. The unchanged, oriented cavity skin and volume are also
          // checked above; no CAD/layer/interface facet is retriangulated.
          buildMixedFaces(candidate, boundary);
          const auto proposed = mixedMeshReport(candidate).at("quality");
          if (better(proposed, quality)) {
            m = std::move(candidate);
            quality = proposed;
            if (indices.size() == 3)
              ++flips32;
            else
              ++flips44;
            improved = true;
            break;
          }
        } catch (const std::exception &) {
          // Reject a degenerate/invalid partition without modifying the mesh.
        }
      }
      if (improved || trials == maxTrials)
        break;
    }
    if (!improved)
      break;
  }
  return {{"method", "bounded transition-tet 3-to-2 and 4-to-4 edge removal"},
          {"trials", trials},
          {"max_trials", maxTrials},
          {"prefiltered_candidates", prefiltered},
          {"accepted_3_to_2", flips32},
          {"accepted_4_to_4", flips44},
          {"removed_cells", flips32},
          {"added_nodes", 0},
          {"added_cells", 0},
          {"preserved_cavity_skin", true},
          {"after", quality}};
}
} // namespace babelmesh
