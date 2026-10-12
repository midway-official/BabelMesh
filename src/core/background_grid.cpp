#include "babelmesh/background_grid.hpp"
#include "babelmesh/cartesian_mesh.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
namespace babelmesh {
namespace {
using Json = nlohmann::json;
using I = std::int64_t;
struct Leaf {
  I x{}, y{}, width{};
  int level{};
};
double scaleOf(const Json &p) {
  const auto unit = p.at("units").at("length").get<std::string>();
  return unit == "mm" ? .001 : unit == "cm" ? .01 : unit == "inch" ? .0254 : 1;
}
Json backgroundSources(const Json &p) {
  Json result = Json::array();
  for (auto s : p.at("sizing").at("sources")) {
    const auto phase = s.value("apply_to", "all");
    if (phase != "all" && phase != "background")
      continue;
    s.erase("apply_to");
    result.push_back(s);
  }
  return result;
}
bool intersects(const Json &source, const std::array<double, 3> &lo,
                const std::array<double, 3> &hi, double scale) {
  const auto kind = source.at("kind");
  const auto &r = source.at("region");
  if (kind == "box") {
    for (std::size_t a = 0; a < 3; ++a)
      if (lo[a] > r.at("bounds").at("maximum")[a].get<double>() * scale ||
          hi[a] < r.at("bounds").at("minimum")[a].get<double>() * scale)
        return false;
    return true;
  }
  if (kind == "sphere") {
    double d2 = 0;
    for (std::size_t a = 0; a < 3; ++a) {
      const double c = r.at("center")[a].get<double>() * scale,
                   d = std::max({lo[a] - c, c - hi[a], 0.0});
      d2 += d * d;
    }
    const double radius = r.at("radius").get<double>() * scale;
    return d2 <= radius * radius;
  }
  throw std::runtime_error(
      "2D background refinement supports box and sphere regions; choose "
      "apply_to=background for supported 3D CAD-distance sources");
}
bool neighbours(const Leaf &a, const Leaf &b) {
  return ((a.x + a.width == b.x || b.x + b.width == a.x) &&
          std::max(a.y, b.y) < std::min(a.y + a.width, b.y + b.width)) ||
         ((a.y + a.width == b.y || b.y + b.width == a.y) &&
          std::max(a.x, b.x) < std::min(a.x + a.width, b.x + b.width));
}
} // namespace
BackgroundGrid generateBackgroundGrid(const Json &p, const CadModel &cad) {
  BackgroundGrid result;
  const auto sources = backgroundSources(p);
  const double scale = scaleOf(p);
  if (p.at("meshing").at("dimension") == 3) {
    auto legacy = p;
    legacy["schema_version"] = 1;
    legacy.erase("meshing");
    legacy.erase("optimization");
    legacy["boundary_layers"]["enabled"] = false;
    legacy["sizing"]["sources"] = sources;
    legacy["sizing"].erase("volume_size");
    const auto grid = makeAdaptiveCartesianMesh(legacy, &cad, false);
    result.points = grid.points_m;
    std::vector<std::vector<FaceKey>> rings(grid.cells.size());
    for (const auto &face : grid.faces) {
      rings.at(static_cast<std::size_t>(face.owner)).push_back(face.vertices);
      if (face.neighbour >= 0) {
        auto reverse = face.vertices;
        std::reverse(reverse.begin(), reverse.end());
        rings.at(static_cast<std::size_t>(face.neighbour))
            .push_back(std::move(reverse));
      }
    }
    for (std::size_t i = 0; i < grid.cells.size(); ++i) {
      if (grid.strict_hex.at(i))
        result.volumes.emplace_back(
            CellType::hexahedron, CellRegion::cartesian,
            FaceKey(grid.cells[i].begin(), grid.cells[i].end()));
      else {
        std::set<std::int32_t> vertices;
        for (const auto &ring : rings[i])
          vertices.insert(ring.begin(), ring.end());
        result.volumes.emplace_back(CellType::polyhedron, CellRegion::cartesian,
                                    FaceKey(vertices.begin(), vertices.end()),
                                    std::move(rings[i]));
      }
    }
    result.report = {
        {"strategy",
         sources.empty() ? "uniform-cartesian" : "face-conforming-octree"},
        {"candidate_cells", grid.cells.size()},
        {"maximum_level", grid.cell_levels.empty()
                              ? 0
                              : *std::max_element(grid.cell_levels.begin(),
                                                  grid.cell_levels.end())},
        {"source_intersection_counts", grid.source_intersection_counts}};
    return result;
  }
  const auto axisName = p.at("meshing").at("section_axis").get<std::string>();
  const int axis = axisName == "x"   ? 0
                   : axisName == "y" ? 1
                                     : 2,
            u = (axis + 1) % 3, v = (axis + 2) % 3;
  const int levels = p.at("sizing").at("max_level").get<int>();
  if (levels < 0 || levels > 20)
    throw std::runtime_error("sizing.max_level must be in [0,20]");
  const I factor = I{1} << levels;
  const auto &b = p.at("background");
  const int nx = b.at("cell_counts")[u], ny = b.at("cell_counts")[v];
  const auto maximum = p.at("resources").at("max_cells").get<std::size_t>();
  auto checkBudget = [&](std::size_t count) {
    if (count > maximum)
      throw std::runtime_error("2D background exceeds resources.max_cells");
    if (static_cast<double>(count) * 3072 / (1024 * 1024) >
        p.at("resources").at("max_memory_mb").get<double>())
      throw std::runtime_error(
          "2D background estimate exceeds resources.max_memory_mb");
  };
  checkBudget(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
  std::array<double, 3>
      lo = b.at("bounds").at("minimum").get<std::array<double, 3>>(),
      hi = b.at("bounds").at("maximum").get<std::array<double, 3>>();
  for (std::size_t a = 0; a < 3; ++a) {
    lo[a] *= scale;
    hi[a] *= scale;
  }
  const double dx = (hi[static_cast<std::size_t>(u)] -
                     lo[static_cast<std::size_t>(u)]) /
                    nx,
               dy = (hi[static_cast<std::size_t>(v)] -
                     lo[static_cast<std::size_t>(v)]) /
                    ny,
               coordinate =
                   p.at("meshing").at("section_coordinate").get<double>() *
                   scale;
  std::vector<Leaf> leaves;
  for (int y = 0; y < ny; ++y)
    for (int x = 0; x < nx; ++x)
      leaves.push_back({x * factor, y * factor, factor, 0});
  const double floor = p.at("sizing").at("minimum_size").get<double>() * scale;
  auto split = [&](const std::vector<bool> &marking) {
    std::vector<Leaf> next;
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      const auto &c = leaves[i];
      if (!marking[i]) {
        next.push_back(c);
        continue;
      }
      if (c.level >= levels || std::max(dx, dy) * static_cast<double>(c.width) /
                                       static_cast<double>(factor) * .5 <
                                   floor * (1 - 1e-12))
        throw std::runtime_error("2D refinement cannot meet target or balance "
                                 "without violating minimum_size/max_level");
      const I w = c.width / 2;
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          next.push_back({c.x + k * w, c.y + j * w, w, c.level + 1});
      if (next.size() > maximum)
        throw std::runtime_error("2D refinement exceeds resources.max_cells");
    }
    if (next.size() > maximum)
      throw std::runtime_error("2D refinement exceeds resources.max_cells");
    checkBudget(next.size());
    leaves = std::move(next);
  };
  for (;;) {
    std::vector<bool> marking(leaves.size(), false);
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      const auto &c = leaves[i];
      auto lower = lo, upper = lo;
      lower[static_cast<std::size_t>(axis)] =
          upper[static_cast<std::size_t>(axis)] = coordinate;
      lower[static_cast<std::size_t>(u)] +=
          dx * static_cast<double>(c.x) / static_cast<double>(factor);
      upper[static_cast<std::size_t>(u)] =
          lower[static_cast<std::size_t>(u)] +
          dx * static_cast<double>(c.width) / static_cast<double>(factor);
      lower[static_cast<std::size_t>(v)] +=
          dy * static_cast<double>(c.y) / static_cast<double>(factor);
      upper[static_cast<std::size_t>(v)] =
          lower[static_cast<std::size_t>(v)] +
          dy * static_cast<double>(c.width) / static_cast<double>(factor);
      const double size = std::max(dx, dy) * static_cast<double>(c.width) /
                          static_cast<double>(factor);
      for (const auto &source : sources)
        if (size >
                source.at("target_size").get<double>() * scale * (1 + 1e-12) &&
            intersects(source, lower, upper, scale))
          marking[i] = true;
    }
    if (std::none_of(marking.begin(), marking.end(), [](bool b) { return b; }))
      break;
    split(marking);
  }
  // A sparse line index finds face neighbours without an all-pairs leaf scan.
  for (;;) {
    std::map<I, std::vector<std::size_t>> vertical, horizontal;
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      const auto &c = leaves[i];
      vertical[c.x].push_back(i);
      vertical[c.x + c.width].push_back(i);
      horizontal[c.y].push_back(i);
      horizontal[c.y + c.width].push_back(i);
    }
    std::vector<bool> marking(leaves.size(), false);
    auto balance = [&](auto &lines, bool vert) {
      for (auto &line : lines) {
        auto &ids = line.second;
        std::sort(ids.begin(), ids.end(), [&](auto a, auto b) {
          return vert ? leaves[a].y < leaves[b].y : leaves[a].x < leaves[b].x;
        });
        for (std::size_t i = 0; i < ids.size(); ++i)
          for (std::size_t j = i + 1; j < ids.size(); ++j) {
            const auto &a = leaves[ids[i]];
            const auto &b = leaves[ids[j]];
            if ((vert ? b.y : b.x) >= (vert ? a.y : a.x) + a.width)
              break;
            if (neighbours(a, b) && std::abs(a.level - b.level) > 1)
              marking[a.level < b.level ? ids[i] : ids[j]] = true;
          }
      }
    };
    balance(vertical, true);
    balance(horizontal, false);
    if (std::none_of(marking.begin(), marking.end(), [](bool b) { return b; }))
      break;
    split(marking);
  }
  std::map<std::array<I, 2>, std::int32_t> nodes;
  std::map<I, std::set<I>> xs, ys;
  auto node = [&](I x, I y) {
    const std::array<I, 2> key{x, y};
    auto old = nodes.find(key);
    if (old != nodes.end())
      return old->second;
    auto point = lo;
    point[static_cast<std::size_t>(axis)] = coordinate;
    point[static_cast<std::size_t>(u)] +=
        dx * static_cast<double>(x) / static_cast<double>(factor);
    point[static_cast<std::size_t>(v)] +=
        dy * static_cast<double>(y) / static_cast<double>(factor);
    const auto id = static_cast<std::int32_t>(result.points.size());
    result.points.push_back(point);
    nodes.emplace(key, id);
    xs[x].insert(y);
    ys[y].insert(x);
    return id;
  };
  for (const auto &c : leaves) {
    node(c.x, c.y);
    node(c.x + c.width, c.y);
    node(c.x + c.width, c.y + c.width);
    node(c.x, c.y + c.width);
  }
  int maxLevel = 0;
  for (const auto &c : leaves) {
    FaceKey polygon;
    for (auto it = ys[c.y].lower_bound(c.x);
         it != ys[c.y].end() && *it < c.x + c.width; ++it)
      polygon.push_back(node(*it, c.y));
    for (auto it = xs[c.x + c.width].lower_bound(c.y);
         it != xs[c.x + c.width].end() && *it < c.y + c.width; ++it)
      polygon.push_back(node(c.x + c.width, *it));
    auto &top = ys[c.y + c.width];
    for (auto it = top.upper_bound(c.x + c.width); it != top.begin();) {
      --it;
      if (*it <= c.x)
        break;
      polygon.push_back(node(*it, c.y + c.width));
    }
    auto &left = xs[c.x];
    for (auto it = left.upper_bound(c.y + c.width); it != left.begin();) {
      --it;
      if (*it <= c.y)
        break;
      polygon.push_back(node(c.x, *it));
    }
    result.polygons.push_back(std::move(polygon));
    maxLevel = std::max(maxLevel, c.level);
  }
  result.report = {{"strategy", sources.empty() ? "uniform-cartesian"
                                                : "edge-conforming-quadtree"},
                   {"candidate_cells", leaves.size()},
                   {"maximum_level", maxLevel},
                   {"extrusion_layers", 1},
                   {"balanced_neighbor_level_difference", 1}};
  return result;
}
} // namespace babelmesh
