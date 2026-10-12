#include "babelmesh/mesh_optimizer.hpp"
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
namespace babelmesh {
namespace {
using Vec = Eigen::Vector3d;
using Json = nlohmann::json;
Vec point(const MixedMesh &m, std::size_t id) {
  const auto &p = m.points_m.at(id);
  return {p[0], p[1], p[2]};
}
// Exchange meshes need not retain generation-region tags. Recognize only an
// exact axis-aligned box shell, including subdivided coarse/fine box faces.
// A warped hex or an arbitrary convex polyhedron must remain eligible.
bool cartesianBox(const MixedMesh &m, std::size_t index) {
  const auto &cell = m.cells[index];
  if (cell.type != CellType::hexahedron && cell.type != CellType::polyhedron)
    return false;
  Vec low = point(m, static_cast<std::size_t>(cell.vertices.front())),
      high = low;
  double coordinateScale = low.cwiseAbs().maxCoeff();
  for (auto node : cell.vertices) {
    const Vec p = point(m, static_cast<std::size_t>(node));
    low = low.cwiseMin(p);
    high = high.cwiseMax(p);
    coordinateScale = std::max(coordinateScale, p.cwiseAbs().maxCoeff());
  }
  const Vec extent = high - low;
  const double tolerance =
      std::max(extent.maxCoeff() * 1e-12,
               coordinateScale * 64 * std::numeric_limits<double>::epsilon());
  if (extent.minCoeff() <= tolerance)
    return false;
  unsigned corners = 0, sides = 0;
  for (auto node : cell.vertices) {
    const Vec p = point(m, static_cast<std::size_t>(node));
    unsigned corner = 0;
    bool isCorner = true;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(p[axis] - high[axis]) <= tolerance)
        corner |= 1U << axis;
      else if (std::abs(p[axis] - low[axis]) > tolerance)
        isCorner = false;
    }
    if (isCorner)
      corners |= 1U << corner;
  }
  if (corners != 255U)
    return false;
  for (const auto &ring : cellFaceRings(cell)) {
    bool onSide = false;
    for (int axis = 0; axis < 3; ++axis)
      for (int side = 0; side < 2; ++side) {
        const double plane = side == 0 ? low[axis] : high[axis];
        if (std::all_of(ring.begin(), ring.end(), [&](auto node) {
              return std::abs(point(m, static_cast<std::size_t>(node))[axis] -
                              plane) <= tolerance;
            })) {
          onSide = true;
          sides |= 1U << (2 * axis + side);
        }
      }
    if (!onSide)
      return false;
  }
  const double boxVolume = extent.prod();
  return sides == 63U &&
         std::abs(m.volumes_m3[index] - boxVolume) <= boxVolume * 1e-10;
}
double objective(const Json &q) {
  return std::pow(q.at("max_non_orthogonality_degrees").get<double>() / 75, 4) +
         std::pow(q.at("p99_non_orthogonality_degrees").get<double>() / 60, 4) +
         std::pow(q.at("max_internal_skewness").get<double>() / .5, 4) +
         std::pow(q.at("p99_internal_skewness").get<double>() / .25, 4) +
         q.value("mean_quartic_face_penalty", 0.0);
}
bool doesNotWorsen(const Json &a, const Json &b) {
  for (const auto *key :
       {"max_non_orthogonality_degrees", "p99_non_orthogonality_degrees",
        "max_internal_skewness", "p99_internal_skewness"})
    if (a.at(key).get<double>() > b.at(key).get<double>() + 1e-8)
      return false;
  return true;
}
std::size_t mergeTransitionTets(MixedMesh &m, std::size_t budget) {
  std::size_t merged = 0;
  auto quality = mixedMeshReport(m).at("quality");
  while (merged < budget &&
         !quality.at("face_quality_thresholds_passed").get<bool>()) {
    std::vector<std::pair<double, std::size_t>> candidates;
    for (std::size_t i = 0; i < m.faces.size(); ++i) {
      const auto &f = m.faces[i];
      if (f.neighbour < 0 || f.vertices.size() != 3)
        continue;
      const auto &a = m.cells.at(static_cast<std::size_t>(f.owner)),
                 &b = m.cells.at(static_cast<std::size_t>(f.neighbour));
      if (a.type != CellType::tetrahedron || b.type != CellType::tetrahedron ||
          a.region != CellRegion::transition ||
          b.region != CellRegion::transition)
        continue;
      const Vec x = point(m, static_cast<std::size_t>(f.vertices[0])),
                y = point(m, static_cast<std::size_t>(f.vertices[1])),
                z = point(m, static_cast<std::size_t>(f.vertices[2])),
                area = (y - x).cross(z - x) * .5, fc = (x + y + z) / 3;
      const auto &op = m.centres_m.at(static_cast<std::size_t>(f.owner)),
                 &np = m.centres_m.at(static_cast<std::size_t>(f.neighbour));
      const Vec owner(op[0], op[1], op[2]), neighbour(np[0], np[1], np[2]),
          direction = neighbour - owner;
      const double denominator = area.dot(direction);
      if (denominator <= 0)
        continue;
      const Vec hit = owner + direction * (area.dot(fc - owner) / denominator);
      const double skew = (fc - hit).norm() / std::sqrt(area.norm());
      if (skew > .24)
        candidates.push_back({skew, i});
    }
    std::sort(candidates.rbegin(), candidates.rend());
    bool improved = false;
    BoundaryMap boundary;
    for (const auto &f : m.faces)
      if (f.neighbour < 0)
        boundary.emplace(faceKey(f.vertices), f.patch);
    for (const auto &item : candidates) {
      const auto &face = m.faces[item.second];
      const auto first = static_cast<std::size_t>(face.owner),
                 second = static_cast<std::size_t>(face.neighbour);
      std::vector<FaceKey> shell;
      const auto common = faceKey(face.vertices);
      std::set<std::int32_t> nodes;
      for (auto cell : {first, second})
        for (auto ring : cellFaceRings(m.cells[cell]))
          if (faceKey(ring) != common) {
            nodes.insert(ring.begin(), ring.end());
            shell.push_back(std::move(ring));
          }
      MixedMesh candidate = m;
      candidate.cells[first] = {CellType::polyhedron, CellRegion::transition,
                                FaceKey(nodes.begin(), nodes.end()),
                                std::move(shell)};
      candidate.cells.erase(candidate.cells.begin() +
                            static_cast<std::ptrdiff_t>(second));
      try {
        buildMixedFaces(candidate,
                        boundary); // Includes positive convex-shell checks.
        const auto proposed = mixedMeshReport(candidate).at("quality");
        if (doesNotWorsen(proposed, quality) &&
            objective(proposed) < objective(quality) * (1 - 1e-7)) {
          m = std::move(candidate);
          quality = proposed;
          ++merged;
          improved = true;
          break;
        }
      } catch (const std::exception &) {
        // A non-convex union or a failed shared-face contract is discarded.
      }
    }
    if (!improved)
      break;
  }
  return merged;
}
} // namespace
void validateOptimization(const Json &p) {
  if (!p.is_object())
    throw std::runtime_error("optimization must be an object");
  for (auto it = p.begin(); it != p.end(); ++it)
    if (it.key() != "enabled" && it.key() != "iterations" &&
        it.key() != "relaxation" && it.key() != "max_cell_growth_fraction" &&
        it.key() != "topology_swaps" && it.key() != "merge_transition_cells" &&
        it.key() != "max_merge_fraction")
      throw std::runtime_error("unknown optimization option: " + it.key());
  if (p.contains("enabled") && !p.at("enabled").is_boolean())
    throw std::runtime_error("optimization.enabled must be boolean");
  if (p.contains("topology_swaps") && !p.at("topology_swaps").is_boolean())
    throw std::runtime_error("optimization.topology_swaps must be boolean");
  if (p.contains("merge_transition_cells") &&
      !p.at("merge_transition_cells").is_boolean())
    throw std::runtime_error("merge_transition_cells must be boolean");
  const double merge = p.value("max_merge_fraction", .1);
  if (!std::isfinite(merge) || merge < 0 || merge > .1)
    throw std::runtime_error("max_merge_fraction must be in [0,0.1]");
  if (p.contains("iterations") &&
      (!p.at("iterations").is_number_integer() ||
       p.at("iterations").get<int>() < 1 || p.at("iterations").get<int>() > 30))
    throw std::runtime_error(
        "optimization.iterations must be an integer in [1,30]");
  const double relax = p.value("relaxation", .5);
  if (!std::isfinite(relax) || relax <= 0 || relax > 1)
    throw std::runtime_error("optimization.relaxation must be in (0,1]");
  const double growth = p.value("max_cell_growth_fraction", 0.0);
  if (!std::isfinite(growth) || growth < 0 || growth > .05)
    throw std::runtime_error(
        "optimization.max_cell_growth_fraction must be in [0,0.05]");
}
void optimizeMixedMesh(MixedMesh &m, const Json &p) {
  validateOptimization(p);
  if (!m.metadata.contains("domain_bounds_constraint") && !m.points_m.empty()) {
    auto low = m.points_m.front(), high = low;
    for (const auto &point : m.points_m)
      for (std::size_t axis = 0; axis < 3; ++axis) {
        low[axis] = std::min(low[axis], point[axis]);
        high[axis] = std::max(high[axis], point[axis]);
      }
    m.metadata["domain_bounds_constraint"] = {
        {"minimum_m", low},
        {"maximum_m", high},
        {"source", "initial imported mesh envelope"}};
  }
  validateMixedMesh(m);
  const auto before = mixedMeshReport(m).at("quality");
  Json report = {
      {"method",
       "constrained transition-node relocation with quality rollback"},
      {"enabled", p.value("enabled", true)},
      {"cells_before", m.cells.size()},
      {"cells_after", m.cells.size()},
      {"added_cells", 0},
      {"allowed_cell_growth_fraction",
       p.value("max_cell_growth_fraction", 0.0)},
      {"before", before}};
  if (!p.value("enabled", true)) {
    report["after"] = before;
    report["accepted_sweeps"] = 0;
    m.metadata["optimization"] = report;
    return;
  }
  const bool planar = m.generation_dimension == 2;
  const int axis = m.extrusion_axis, u = (axis + 1) % 3, v = (axis + 2) % 3;
  std::vector<std::size_t> group(m.points_m.size());
  std::vector<std::vector<std::size_t>> members;
  std::map<std::array<double, 2>, std::size_t> pairs;
  for (std::size_t i = 0; i < m.points_m.size(); ++i) {
    std::size_t id = members.size();
    if (planar) {
      const auto key =
          std::array<double, 2>{m.points_m[i][static_cast<std::size_t>(u)],
                                m.points_m[i][static_cast<std::size_t>(v)]};
      const auto found = pairs.find(key);
      if (found != pairs.end())
        id = found->second;
      else
        pairs.emplace(key, id);
    }
    if (id == members.size())
      members.emplace_back();
    members[id].push_back(i);
    group[i] = id;
  }
  std::vector<bool> pinned(members.size(), false);
  const bool regionTagsAvailable =
      m.metadata.value("generation_region_tags_available", true);
  // Lost wall-layer tags cannot be reconstructed from a shell. Keep imported
  // geometry unchanged rather than inadvertently changing layer heights.
  if (!regionTagsAvailable)
    std::fill(pinned.begin(), pinned.end(), true);
  std::vector<std::set<std::size_t>> neighbours(members.size());
  std::size_t geometricCartesianCells = 0;
  for (std::size_t c = 0; c < m.cells.size(); ++c) {
    const auto &cell = m.cells[c];
    const bool box = cartesianBox(m, c);
    geometricCartesianCells += box ? 1U : 0U;
    if (cell.region != CellRegion::transition || box)
      for (auto i : cell.vertices)
        pinned[group.at(static_cast<std::size_t>(i))] = true;
  }
  for (const auto &face : m.faces) {
    const bool symmetryCap =
        planar && face.neighbour < 0 &&
        m.patches.at(static_cast<std::size_t>(face.patch)).kind == "symmetry" &&
        (m.patches.at(static_cast<std::size_t>(face.patch)).name == "front" ||
         m.patches.at(static_cast<std::size_t>(face.patch)).name == "back");
    if (face.neighbour < 0 && !symmetryCap)
      for (auto i : face.vertices)
        pinned[group.at(static_cast<std::size_t>(i))] = true;
    for (std::size_t j = 0; j < face.vertices.size(); ++j) {
      const auto a = group.at(static_cast<std::size_t>(face.vertices[j])),
                 b = group.at(static_cast<std::size_t>(
                     face.vertices[(j + 1) % face.vertices.size()]));
      if (a != b) {
        neighbours[a].insert(b);
        neighbours[b].insert(a);
      }
    }
  }
  const auto eligible =
      static_cast<std::size_t>(std::count(pinned.begin(), pinned.end(), false));
  const auto initial = m.points_m;
  Json quality = before;
  double score = objective(quality);
  int accepted = 0, trials = 0, circumcentreSweeps = 0, alignmentSweeps = 0;
  std::size_t guardedGroups = 0;
  for (int iteration = 0; eligible > 0 && iteration < p.value("iterations", 6);
       ++iteration) {
    const auto old = m.points_m;
    std::vector<Vec> delta(members.size(), Vec::Zero());
    for (std::size_t id = 0; id < members.size(); ++id)
      if (!pinned[id] && !neighbours[id].empty()) {
        const Vec current = point(m, members[id].front());
        Vec target = Vec::Zero();
        double weights = 0;
        for (auto other : neighbours[id]) {
          Vec q = point(m, members[other].front());
          if (planar)
            q[axis] = current[axis];
          const double weight = 1 / std::max((q - current).norm(), 1e-20);
          target += q * weight;
          weights += weight;
        }
        delta[id] = target / weights - current;
        if (planar)
          delta[id][axis] = 0;
      }
    const auto graphDelta = delta;
    std::vector<Vec> circumcentreDelta(members.size(), Vec::Zero());
    std::vector<Vec> alignmentDelta(members.size(), Vec::Zero());
    if (!planar) {
      std::vector<Vec> targets(members.size(), Vec::Zero());
      std::vector<double> weights(members.size(), 0);
      for (std::size_t c = 0; c < m.cells.size(); ++c) {
        const auto &cell = m.cells[c];
        if (cell.type != CellType::tetrahedron ||
            cell.region != CellRegion::transition)
          continue;
        const Vec a = point(m, static_cast<std::size_t>(cell.vertices[0]));
        Eigen::Matrix3d matrix;
        Vec rhs;
        for (int j = 0; j < 3; ++j) {
          const Vec edge =
              point(m, static_cast<std::size_t>(
                           cell.vertices[static_cast<std::size_t>(j + 1)])) -
              a;
          matrix.row(j) = 2 * edge.transpose();
          rhs[j] = edge.squaredNorm();
        }
        const Eigen::FullPivLU<Eigen::Matrix3d> lu(matrix);
        if (!lu.isInvertible())
          continue;
        const Vec centre = a + lu.solve(rhs);
        if (!centre.allFinite())
          continue;
        for (auto node : cell.vertices) {
          const auto id = group.at(static_cast<std::size_t>(node));
          targets[id] += m.volumes_m3[c] * centre;
          weights[id] += m.volumes_m3[c];
        }
      }
      for (std::size_t id = 0; id < members.size(); ++id)
        if (!pinned[id] && weights[id] > 0 && !neighbours[id].empty()) {
          const Vec current = point(m, members[id].front());
          Vec move = targets[id] / weights[id] - current;
          double shortest = std::numeric_limits<double>::max();
          for (auto other : neighbours[id])
            shortest = std::min(
                shortest, (point(m, members[other].front()) - current).norm());
          if (move.norm() > shortest * .4)
            move *= shortest * .4 / move.norm();
          circumcentreDelta[id] = move;
        }
      std::vector<double> alignmentWeights(members.size(), 0);
      for (const auto &face : m.faces) {
        if (face.neighbour < 0 || face.vertices.size() != 3)
          continue;
        const Vec a = point(m, static_cast<std::size_t>(face.vertices[0])),
                  b = point(m, static_cast<std::size_t>(face.vertices[1])),
                  c = point(m, static_cast<std::size_t>(face.vertices[2]));
        const Vec area = (b - a).cross(c - a) * .5, fc = (a + b + c) / 3;
        const auto &op = m.centres_m.at(static_cast<std::size_t>(face.owner)),
                   &np =
                       m.centres_m.at(static_cast<std::size_t>(face.neighbour));
        const Vec owner(op[0], op[1], op[2]), neighbour(np[0], np[1], np[2]),
            direction = neighbour - owner;
        const double denominator = area.dot(direction);
        if (denominator <= 0)
          continue;
        const double t =
            std::clamp(area.dot(fc - owner) / denominator, 0.0, 1.0);
        const Vec error = fc - (owner + t * direction);
        const double skew = error.norm() / std::sqrt(area.norm());
        if (skew < .12)
          continue;
        for (const auto &item : {std::make_pair(face.owner, 1 - t),
                                 std::make_pair(face.neighbour, t)}) {
          const auto &cell = m.cells.at(static_cast<std::size_t>(item.first));
          if (cell.type != CellType::tetrahedron)
            continue;
          for (auto node : cell.vertices) {
            const auto id = group.at(static_cast<std::size_t>(node));
            if (pinned[id] ||
                std::find(face.vertices.begin(), face.vertices.end(), node) !=
                    face.vertices.end())
              continue;
            // The opposite tetrahedron vertex contributes one quarter of its
            // volume centroid. Correct the centre-line hit in the face plane.
            const double weight = item.second * skew * skew;
            alignmentDelta[id] += error * (4 * weight);
            alignmentWeights[id] += weight;
          }
        }
      }
      for (std::size_t id = 0; id < members.size(); ++id)
        if (alignmentWeights[id] > 0 && !neighbours[id].empty()) {
          alignmentDelta[id] /= alignmentWeights[id];
          double shortest = std::numeric_limits<double>::max();
          for (auto other : neighbours[id])
            shortest = std::min(shortest, (point(m, members[other].front()) -
                                           point(m, members[id].front()))
                                              .norm());
          if (alignmentDelta[id].norm() > shortest * .25)
            alignmentDelta[id] *= shortest * .25 / alignmentDelta[id].norm();
        }
    }
    bool improved = false;
    // Try a second, volume-weighted tetrahedron circumcentre target when the
    // graph sweep stalls. The same full FVM quality and shell gates apply.
    for (int method = 0; method < (planar ? 1 : 3); ++method) {
      delta = method == 0   ? graphDelta
              : method == 1 ? circumcentreDelta
                            : alignmentDelta;
      for (double factor = p.value("relaxation", .5); factor >= .015;
           factor *= .5) {
        for (int correction = 0; correction < 8; ++correction) {
          ++trials;
          for (std::size_t id = 0; id < members.size(); ++id)
            for (auto node : members[id]) {
              const Vec q = Vec(old[node][0], old[node][1], old[node][2]) +
                            factor * delta[id];
              m.points_m[node] = {q.x(), q.y(), q.z()};
            }
          try {
            validateMixedMesh(m);
            const auto candidate = mixedMeshReport(m).at("quality");
            const double proposed = objective(candidate);
            const bool angleOk =
                candidate.at("max_non_orthogonality_degrees").get<double>() <=
                quality.at("max_non_orthogonality_degrees").get<double>() +
                    1e-8;
            const bool skewOk =
                candidate.at("max_internal_skewness").get<double>() <=
                quality.at("max_internal_skewness").get<double>() + 1e-8;
            const bool percentilesOk =
                candidate.at("p99_non_orthogonality_degrees").get<double>() <=
                    quality.at("p99_non_orthogonality_degrees").get<double>() +
                        1e-8 &&
                candidate.at("p99_internal_skewness").get<double>() <=
                    quality.at("p99_internal_skewness").get<double>() + 1e-8;
            if (proposed < score * (1 - 1e-7) && angleOk && skewOk &&
                percentilesOk) {
              score = proposed;
              quality = candidate;
              improved = true;
              ++accepted;
              if (method == 1)
                ++circumcentreSweeps;
              if (method == 2)
                ++alignmentSweeps;
              break;
            }
            if (angleOk && skewOk)
              break;
            // A few sensitive stars should not prevent all other transition
            // vertices from improving. Hold their complete incident cells at
            // the old coordinates and retry the same bounded sweep.
            bool guarded = false;
            for (const auto &key : {"worst_non_orthogonality_face",
                                    "worst_internal_skewness_face"}) {
              if ((std::string(key) == "worst_non_orthogonality_face" &&
                   angleOk) ||
                  (std::string(key) == "worst_internal_skewness_face" &&
                   skewOk))
                continue;
              const auto &face =
                  m.faces.at(candidate.at(key).get<std::size_t>());
              for (auto cell : {face.owner, face.neighbour})
                if (cell >= 0)
                  for (auto node :
                       m.cells.at(static_cast<std::size_t>(cell)).vertices) {
                    const auto id = group.at(static_cast<std::size_t>(node));
                    if (delta[id].squaredNorm() > 0) {
                      delta[id] = Vec::Zero();
                      ++guardedGroups;
                      guarded = true;
                    }
                  }
            }
            if (!guarded)
              break;
          } catch (const std::exception &) {
            break; // Invalid shells are rolled back before trying a smaller
                   // step.
          }
        }
        if (improved)
          break;
        m.points_m = old;
      }
      if (improved)
        break;
    }
    if (!improved) {
      m.points_m = old;
      validateMixedMesh(m);
      break;
    }
  }
  std::size_t merged = 0;
  // A bounded single-node pass resolves conflicting batch targets. Fixed
  // geometry and all four worst/percentile gates still apply to every move.
  int localAccepted = 0, localTrials = 0;
  if (!planar && regionTagsAvailable &&
      !quality.at("face_quality_thresholds_passed").get<bool>()) {
    struct Candidate {
      double skew;
      std::size_t group;
      Vec move;
    };
    std::vector<Candidate> candidates;
    for (const auto &face : m.faces) {
      if (face.neighbour < 0 || face.vertices.size() != 3)
        continue;
      const Vec a = point(m, static_cast<std::size_t>(face.vertices[0])),
                b = point(m, static_cast<std::size_t>(face.vertices[1])),
                c = point(m, static_cast<std::size_t>(face.vertices[2]));
      const Vec area = (b - a).cross(c - a) * .5, fc = (a + b + c) / 3;
      const auto &op = m.centres_m[static_cast<std::size_t>(face.owner)],
                 &np = m.centres_m[static_cast<std::size_t>(face.neighbour)];
      const Vec owner(op[0], op[1], op[2]), neighbour(np[0], np[1], np[2]),
          direction = neighbour - owner;
      const double denominator = area.dot(direction);
      if (denominator <= 0)
        continue;
      const double t = area.dot(fc - owner) / denominator;
      const Vec error = fc - (owner + t * direction);
      const double skew = error.norm() / std::sqrt(area.norm());
      if (skew <= .2)
        continue;
      for (const auto &item : {std::make_pair(face.owner, 1 - t),
                               std::make_pair(face.neighbour, t)}) {
        const auto &cell = m.cells[static_cast<std::size_t>(item.first)];
        if (cell.type != CellType::tetrahedron ||
            cell.region != CellRegion::transition || item.second <= .05)
          continue;
        for (auto node : cell.vertices) {
          const auto id = group[static_cast<std::size_t>(node)];
          if (pinned[id] ||
              std::find(face.vertices.begin(), face.vertices.end(), node) !=
                  face.vertices.end())
            continue;
          Vec move = error * (4 / item.second);
          double shortest = std::numeric_limits<double>::max();
          for (auto other : neighbours[id])
            shortest = std::min(shortest, (point(m, members[other].front()) -
                                           point(m, members[id].front()))
                                              .norm());
          if (move.norm() > shortest * .2)
            move *= shortest * .2 / move.norm();
          candidates.push_back({skew, id, move});
        }
      }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const auto &a, const auto &b) { return a.skew > b.skew; });
    std::set<std::size_t> attempted;
    for (const auto &candidate : candidates) {
      if (attempted.size() == 32)
        break;
      if (!attempted.insert(candidate.group).second)
        continue;
      const auto old = m.points_m;
      bool improved = false;
      for (double factor = .5; factor >= .0625; factor *= .5) {
        ++localTrials;
        for (auto node : members[candidate.group]) {
          const Vec p(old[node][0], old[node][1], old[node][2]);
          const Vec q = p + candidate.move * factor;
          m.points_m[node] = {q.x(), q.y(), q.z()};
        }
        try {
          validateMixedMesh(m);
          const auto proposed = mixedMeshReport(m).at("quality");
          if (doesNotWorsen(proposed, quality) &&
              objective(proposed) < score * (1 - 1e-7)) {
            score = objective(proposed);
            quality = proposed;
            ++localAccepted;
            improved = true;
            break;
          }
        } catch (const std::exception &) {
          // The next smaller move is evaluated from the original geometry.
        }
      }
      if (!improved)
        m.points_m = old;
    }
    validateMixedMesh(m);
  }
  Json topology = {{"enabled", false}};
  Json postSwapRelocation = {{"enabled", false}};
  if (!planar && regionTagsAvailable && p.value("topology_swaps", true)) {
    topology = improveTransitionTetTopology(
        m,
        static_cast<std::size_t>(std::min(256, p.value("iterations", 6) * 32)));
    topology["enabled"] = true;
    quality = mixedMeshReport(m).at("quality");
    score = objective(quality);
    if (topology.at("accepted_3_to_2").get<std::size_t>() +
                topology.at("accepted_4_to_4").get<std::size_t>() >
            0 &&
        !quality.at("face_quality_thresholds_passed").get<bool>()) {
      // Rebuild the adjacency after swaps before trying a short relocation
      // pass. Explicitly disable topology/merging to bound the schedule and
      // prevent recursion beyond this one smoothing pass.
      auto smoothing = p;
      smoothing["topology_swaps"] = false;
      smoothing["merge_transition_cells"] = false;
      smoothing["iterations"] = std::min(2, p.value("iterations", 6));
      optimizeMixedMesh(m, smoothing);
      postSwapRelocation = m.metadata.at("optimization");
      quality = mixedMeshReport(m).at("quality");
      score = objective(quality);
    }
  }
  if (!planar && regionTagsAvailable &&
      p.value("merge_transition_cells", false)) {
    merged = mergeTransitionTets(
        m, static_cast<std::size_t>(
               std::floor(m.cells.size() * p.value("max_merge_fraction", .1))));
    quality = mixedMeshReport(m).at("quality");
    score = objective(quality);
  }
  std::size_t moved = 0;
  double maxMove = 0;
  for (std::size_t i = 0; i < m.points_m.size(); ++i) {
    const Vec original(initial[i][0], initial[i][1], initial[i][2]);
    const double distance = (point(m, i) - original).norm();
    if (distance > 0) {
      ++moved;
      maxMove = std::max(maxMove, distance);
    }
  }
  report.update(
      {{"eligible_node_groups", eligible},
       {"geometrically_protected_cartesian_cells", geometricCartesianCells},
       {"generation_region_tags_available", regionTagsAvailable},
       {"missing_region_tags_policy", "preserve all vertices and cells"},
       {"cells_after", m.cells.size()},
       {"merged_transition_tet_pairs", merged},
       {"topology_swaps", topology},
       {"post_swap_relocation", postSwapRelocation},
       {"fixed_node_groups", members.size() - eligible},
       {"accepted_sweeps", accepted},
       {"line_search_trials", trials},
       {"circumcentre_sweeps", circumcentreSweeps},
       {"face_alignment_sweeps", alignmentSweeps},
       {"local_face_alignment_moves", localAccepted},
       {"local_face_alignment_trials", localTrials},
       {"guarded_sensitive_node_groups", guardedGroups},
       {"moved_vertices", moved},
       {"max_displacement_m", maxMove},
       {"objective_before", objective(before)},
       {"objective_after", score},
       {"after", quality},
       {"constraints",
        {"fixed tagged or geometrically detected Cartesian vertices",
         "fixed physical-boundary and tagged boundary-layer vertices",
         "paired planar extrusion vertices", "positive convex shells",
         "fixed closed background/domain bounds",
         "max and p99 face-quality metrics do not worsen", "no added cells"}}});
  m.metadata["optimization"] = report;
  std::cout << "Optimizer: " << accepted << " accepted sweeps, " << moved
            << " vertices moved, 0 cells added, objective " << objective(before)
            << " -> " << score << '\n'
            << std::flush;
}
} // namespace babelmesh
