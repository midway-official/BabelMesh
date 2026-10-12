#include "babelmesh/mesh_conversion.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
namespace babelmesh {
namespace {
FaceKey rotateRing(FaceKey ring) {
  std::rotate(ring.begin(), std::min_element(ring.begin(), ring.end()),
              ring.end());
  return ring;
}
bool matches(const MixedCell &cell, const std::vector<FaceKey> &faces) {
  std::set<FaceKey> want, actual;
  for (const auto &f : faces)
    want.insert(rotateRing(f));
  for (const auto &f : cellFaceRings(cell))
    actual.insert(rotateRing(f));
  return want == actual;
}
MixedCell recognize(std::vector<FaceKey> faces) {
  std::set<std::int32_t> vertices;
  for (const auto &face : faces)
    vertices.insert(face.begin(), face.end());
  MixedCell result{CellType::polyhedron, CellRegion::transition,
                   FaceKey(vertices.begin(), vertices.end()), faces};
  for (const auto &face : faces) {
    auto base = face;
    std::reverse(base.begin(), base.end());
    CellType type = CellType::polyhedron;
    if (vertices.size() == 4 && faces.size() == 4 && base.size() == 3)
      type = CellType::tetrahedron;
    else if (vertices.size() == 5 && faces.size() == 5 && base.size() == 4)
      type = CellType::pyramid;
    else if (vertices.size() == 6 && faces.size() == 5 && base.size() == 3)
      type = CellType::prism;
    else if (vertices.size() == 8 && faces.size() == 6 && base.size() == 4)
      type = CellType::hexahedron;
    if (type == CellType::polyhedron)
      continue;
    FaceKey candidate = base;
    const std::set<std::int32_t> onBase(base.begin(), base.end());
    if (type == CellType::tetrahedron || type == CellType::pyramid) {
      for (auto node : vertices)
        if (!onBase.count(node))
          candidate.push_back(node);
    } else {
      bool okay = true;
      for (auto node : base) {
        std::set<std::int32_t> opposite;
        for (const auto &f : faces)
          for (std::size_t j = 0; j < f.size(); ++j)
            if (f[j] == node) {
              const auto before = f[(j + f.size() - 1) % f.size()],
                         after = f[(j + 1) % f.size()];
              if (!onBase.count(before))
                opposite.insert(before);
              if (!onBase.count(after))
                opposite.insert(after);
            }
        if (opposite.size() != 1) {
          okay = false;
          break;
        }
        candidate.push_back(*opposite.begin());
      }
      if (!okay)
        continue;
    }
    MixedCell standard{type, CellRegion::transition, candidate};
    if (matches(standard, faces))
      return standard;
  }
  return result;
}
// Native BabelSim meshes carry physical patches but no generation metadata.
// Infer a planar rule only when the complete topology is an exact paired slab
// with two opposite symmetry caps; a thin 3D volume is insufficient evidence.
void inferPlanarExtrusion(MixedMesh &m) {
  if (std::any_of(m.cells.begin(), m.cells.end(), [](const auto &cell) {
        return cell.type == CellType::tetrahedron ||
               cell.type == CellType::pyramid;
      }))
    return;
  int front = -1, back = -1;
  for (std::size_t i = 0; i < m.patches.size(); ++i) {
    if (m.patches[i].kind != "symmetry")
      continue;
    if (m.patches[i].name == "front")
      front = static_cast<int>(i);
    if (m.patches[i].name == "back")
      back = static_cast<int>(i);
  }
  if (front < 0 || back < 0)
    return;
  for (int axis = 0; axis < 3; ++axis) {
    double low = m.points_m.front()[static_cast<std::size_t>(axis)], high = low;
    for (const auto &p : m.points_m) {
      low = std::min(low, p[static_cast<std::size_t>(axis)]);
      high = std::max(high, p[static_cast<std::size_t>(axis)]);
    }
    if (high <= low)
      continue;
    // Exact levels also make the optimizer's in-plane pairing unambiguous.
    if (std::any_of(m.points_m.begin(), m.points_m.end(), [&](const auto &p) {
          return p[static_cast<std::size_t>(axis)] != low &&
                 p[static_cast<std::size_t>(axis)] != high;
        }))
      continue;
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    bool paired = true;
    for (const auto &cell : m.cells) {
      std::map<std::array<double, 2>, unsigned> levels;
      for (auto node : cell.vertices) {
        const auto &p = m.points_m[static_cast<std::size_t>(node)];
        const auto key = std::array<double, 2>{p[static_cast<std::size_t>(u)],
                                               p[static_cast<std::size_t>(v)]};
        const unsigned level = p[static_cast<std::size_t>(axis)] == low ? 1U : 2U;
        if (levels[key] & level)
          paired = false;
        levels[key] |= level;
      }
      for (const auto &level : levels)
        paired = paired && level.second == 3U;
    }
    if (!paired)
      continue;
    int frontLevel = -1, backLevel = -1;
    std::size_t capFaces = 0;
    for (const auto &face : m.faces) {
      if (face.neighbour >= 0)
        continue;
      const double plane = m.points_m[static_cast<std::size_t>(face.vertices[0])]
                                        [static_cast<std::size_t>(axis)];
      const bool cap = std::all_of(face.vertices.begin(), face.vertices.end(),
                                  [&](auto node) {
        return m.points_m[static_cast<std::size_t>(node)]
                         [static_cast<std::size_t>(axis)] == plane;
      });
      const bool namedCap = face.patch == front || face.patch == back;
      if (cap != namedCap) {
        paired = false;
        break;
      }
      if (cap) {
        ++capFaces;
        int &label = face.patch == front ? frontLevel : backLevel;
        const int level = plane == low ? 0 : 1;
        if (label >= 0 && label != level)
          paired = false;
        label = level;
      }
    }
    if (paired && capFaces == 2 * m.cells.size() && frontLevel >= 0 &&
        backLevel >= 0 && frontLevel != backLevel) {
      m.generation_dimension = 2;
      m.extrusion_axis = axis;
      m.extrusion_thickness_m = high - low;
      m.metadata["inferred_generation_rule"] =
          "exact paired slab with complete front/back symmetry caps";
      validateMixedMesh(m);
      return;
    }
  }
}
} // namespace
void reconstructFaceCells(MixedMesh &m) {
  std::int32_t maximum = -1;
  for (const auto &f : m.faces)
    maximum = std::max({maximum, f.owner, f.neighbour});
  if (maximum < 0)
    throw std::runtime_error("exchange mesh has no volume cells");
  std::vector<std::vector<FaceKey>> rings(static_cast<std::size_t>(maximum) +
                                          1);
  for (const auto &face : m.faces) {
    if (face.owner < 0)
      throw std::runtime_error("exchange face has no owner");
    rings.at(static_cast<std::size_t>(face.owner)).push_back(face.vertices);
    if (face.neighbour >= 0) {
      auto reverse = face.vertices;
      std::reverse(reverse.begin(), reverse.end());
      rings.at(static_cast<std::size_t>(face.neighbour))
          .push_back(std::move(reverse));
    }
  }
  m.cells.clear();
  for (auto &shell : rings)
    m.cells.push_back(recognize(std::move(shell)));
  validateMixedMesh(m);
}
MixedMesh readMixedBabelSimV3(const std::filesystem::path &file) {
  std::ifstream in(file);
  if (!in)
    throw std::runtime_error("cannot read BabelSim mesh");
  std::string word;
  int version = 0;
  in >> word >> version;
  if (word != "BABELSIM_MESH" || version != 3)
    throw std::runtime_error("expected BABELSIM_MESH 3");
  MixedMesh m;
  std::size_t count = 0;
  auto section = [&](const char *expected) {
    in >> word >> count;
    if (!in || word != expected || count == 0 ||
        count >=
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
      throw std::runtime_error(std::string("invalid BabelSim ") + expected +
                               " section");
  };
  section("vertices");
  m.points_m.resize(count);
  for (auto &p : m.points_m)
    in >> p[0] >> p[1] >> p[2];
  section("faces");
  m.faces.resize(count);
  for (auto &f : m.faces) {
    std::size_t size = 0;
    in >> word >> size >> f.owner >> f.neighbour >> f.patch;
    if (word != "face" || size < 3 || size > m.points_m.size())
      throw std::runtime_error("invalid BabelSim face");
    f.vertices.resize(size);
    for (auto &node : f.vertices)
      in >> node;
  }
  section("patches");
  m.patches.resize(count);
  for (auto &p : m.patches) {
    in >> word >> p.name >> p.kind;
    if (word != "patch")
      throw std::runtime_error("invalid BabelSim patch");
  }
  in >> word;
  if (!in || word != "end")
    throw std::runtime_error("incomplete BabelSim mesh");
  if (in >> word)
    throw std::runtime_error("trailing BabelSim mesh data");
  reconstructFaceCells(m);
  inferPlanarExtrusion(m);
  m.metadata["conversion_source"] = {{"format", "babelsim-v3"},
                                     {"file", file.string()}};
  m.metadata["generation_region_tags_available"] = false;
  return m;
}
MixedMesh readExchangeMesh(const std::filesystem::path &file,
                           const std::string &requested) {
  auto format = requested;
  if (format == "auto") {
    const auto ext = file.extension().string();
    format = ext == ".bsm"                                      ? "babelsim-v3"
             : ext == ".json"                                   ? "checkpoint"
             : ext == ".cgns" || ext == ".cga" || ext == ".cgh" ? "cgns"
             : ext == ".msh" || ext == ".cas" || ext == ".grd"  ? "fluent-ascii"
                                                                : "unknown";
  }
  if (format == "babelsim-v3")
    return readMixedBabelSimV3(file);
  if (format == "checkpoint")
    return readMixedCheckpoint(file);
  if (format == "fluent-ascii")
    return readFluentAscii(file);
  if (format == "cgns" || format == "cgns-poly")
    return readCgns(file);
  throw std::runtime_error(
      "unsupported input format; use BabelSim v3, checkpoint, Fluent ASCII or "
      "CGNS; native .ccm/.sim/.hex/.igg and Fluent HDF5 need a vendor exchange "
      "export");
}
void writeExchangeMesh(const MixedMesh &mesh, const std::filesystem::path &file,
                       const std::string &format) {
  if (format == "babelsim-v3")
    writeMixedBabelSimV3(mesh, file);
  else if (format == "checkpoint")
    writeMixedCheckpoint(mesh, file);
  else if (format == "vtu")
    writeMixedVtu(mesh, file);
  else if (format == "fluent-ascii")
    writeFluentAscii(mesh, file);
  else if (format == "cgns" || format == "cgns-poly")
    writeCgns(mesh, file, format == "cgns-poly");
  else
    throw std::runtime_error("unsupported output exchange format");
}
} // namespace babelmesh
