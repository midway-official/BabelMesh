#include "babelmesh/mesh_conversion.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
namespace babelmesh {
namespace {
using Json = nlohmann::json;
std::string hex(std::size_t value) {
  std::ostringstream out;
  out << std::hex << value;
  return out.str();
}
int code(const std::string &kind) {
  return kind == "wall"        ? 3
         : kind == "inlet"     ? 10
         : kind == "outlet"    ? 5
         : kind == "symmetry"  ? 7
         : kind == "processor" ? 24
                               : 3;
}
std::string kindOf(int bc) {
  return bc == 3                           ? "wall"
         : bc == 4 || bc == 10 || bc == 20 ? "inlet"
         : bc == 5 || bc == 36             ? "outlet"
         : bc == 7                         ? "symmetry"
                                           : "generic";
}
int cellCode(CellType type) {
  switch (type) {
  case CellType::tetrahedron:
    return 2;
  case CellType::hexahedron:
    return 4;
  case CellType::pyramid:
    return 5;
  case CellType::prism:
    return 6;
  case CellType::polyhedron:
    return 7;
  }
  throw std::runtime_error("unsupported Fluent cell type");
}
class Tokens {
  std::istream &input;
  std::string buffered;

public:
  explicit Tokens(std::istream &in) : input(in) {}
  void putback(std::string value) { buffered = std::move(value); }
  std::string next() {
    if (!buffered.empty()) {
      auto result = buffered;
      buffered.clear();
      return result;
    }
    int c;
    do {
      c = input.get();
    } while (c != EOF && std::isspace(c));
    if (c == EOF)
      return {};
    if (c == '(' || c == ')')
      return std::string(1, static_cast<char>(c));
    if (c == 0 || (c < 32 && c != '\t' && c != '\n' && c != '\r'))
      throw std::runtime_error(
          "binary Fluent is unsupported; export ASCII .msh or CGNS");
    std::string result;
    if (c == '"') {
      bool closed = false;
      while ((c = input.get()) != EOF) {
        if (c == '"') {
          closed = true;
          break;
        }
        if (c == '\\') {
          c = input.get();
          if (c == EOF)
            break;
        }
        result.push_back(static_cast<char>(c));
      }
      if (!closed)
        throw std::runtime_error("unterminated Fluent string");
      return result;
    }
    result.push_back(static_cast<char>(c));
    while (input.peek() != EOF && !std::isspace(input.peek()) &&
           input.peek() != '(' && input.peek() != ')')
      result.push_back(static_cast<char>(input.get()));
    return result;
  }
  std::vector<std::string> body() {
    if (next() != "(")
      throw std::runtime_error("expected Fluent section list");
    std::vector<std::string> result;
    int depth = 1;
    while (depth) {
      auto t = next();
      if (t.empty())
        throw std::runtime_error("truncated Fluent section");
      if (t == "(")
        ++depth;
      else if (t == ")")
        --depth;
      else
        result.push_back(std::move(t));
    }
    return result;
  }
};
std::size_t number(const std::string &s, int base = 16) {
  std::size_t used = 0;
  const auto value = std::stoull(s, &used, base);
  if (used != s.size() ||
      value >=
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    throw std::runtime_error("invalid Fluent index");
  return static_cast<std::size_t>(value);
}
} // namespace
void writeFluentAscii(const MixedMesh &m, const std::filesystem::path &file) {
  std::ofstream out(file);
  if (!out)
    throw std::runtime_error("cannot create Fluent file");
  out << std::setprecision(17);
  out << "(0 \"BabelMesh conforming face-based mesh; length unit m\")\n(2 3)\n";
  Json regions = Json::array();
  for (const auto &cell : m.cells)
    regions.push_back(static_cast<int>(cell.region));
  const Json data = {
      {"generation_dimension", m.generation_dimension},
      {"generation_region_tags_available",
       m.metadata.value("generation_region_tags_available", true)},
      {"extrusion_axis", m.extrusion_axis},
      {"extrusion_thickness_m", m.extrusion_thickness_m},
      {"patches", Json::array()},
      {"cell_regions", regions}};
  auto metadata = data;
  if (m.metadata.contains("domain_bounds_constraint"))
    metadata["domain_bounds_constraint"] =
        m.metadata.at("domain_bounds_constraint");
  for (const auto &patch : m.patches)
    metadata["patches"].push_back({patch.name, patch.kind});
  out << "(0 " << Json("BabelMesh:" + metadata.dump()).dump() << ")\n";
  out << "(10 (0 1 " << hex(m.points_m.size()) << " 0 3))\n(10 (1 1 "
      << hex(m.points_m.size()) << " 1 3)(\n";
  for (const auto &p : m.points_m)
    out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
  out << "))\n(12 (0 1 " << hex(m.cells.size()) << " 0))\n(12 (2 1 "
      << hex(m.cells.size()) << " 1 0)(\n";
  for (const auto &cell : m.cells)
    out << hex(static_cast<std::size_t>(cellCode(cell.type))) << ' ';
  out << "\n))\n(13 (0 1 " << hex(m.faces.size()) << " 0))\n";
  std::vector<std::vector<std::size_t>> zones(m.patches.size() + 1);
  for (std::size_t i = 0; i < m.faces.size(); ++i) {
    const auto &face = m.faces[i];
    zones[face.neighbour >= 0 ? 0 : static_cast<std::size_t>(face.patch) + 1]
        .push_back(i);
  }
  std::size_t first = 1;
  for (std::size_t z = 0; z < zones.size(); ++z) {
    if (zones[z].empty())
      continue;
    out << "(13 (" << hex(z + 3) << ' ' << hex(first) << ' '
        << hex(first + zones[z].size() - 1) << ' '
        << hex(static_cast<std::size_t>(z ? code(m.patches[z - 1].kind) : 2))
        << " 0)(\n";
    for (auto f : zones[z]) {
      const auto &face = m.faces[f];
      out << hex(
          face.vertices.size()); // Fluent legacy file normals point towards c0;
                                 // reverse BabelSim's outward owner ring.
      for (auto it = face.vertices.rbegin(); it != face.vertices.rend(); ++it)
        out << ' ' << hex(static_cast<std::size_t>(*it) + 1);
      out << ' ' << hex(static_cast<std::size_t>(face.owner) + 1) << ' '
          << (face.neighbour < 0
                  ? "0"
                  : hex(static_cast<std::size_t>(face.neighbour) + 1))
          << '\n';
    }
    out << "))\n";
    first += zones[z].size();
  }
  out << "(39 (2 fluid fluid)())\n";
  for (std::size_t z = 0; z < zones.size(); ++z)
    if (!zones[z].empty())
      out << "(39 (" << z + 3 << ' '
          << (z ? m.patches[z - 1].kind == "generic"  ? "wall"
                  : m.patches[z - 1].kind == "inlet"  ? "velocity-inlet"
                  : m.patches[z - 1].kind == "outlet" ? "pressure-outlet"
                                                      : m.patches[z - 1].kind
                : "interior")
          << ' ' << Json(z ? m.patches[z - 1].name : "interior").dump()
          << ")())\n";
  out.flush();
  if (!out)
    throw std::runtime_error("Fluent write failed");
}
MixedMesh readFluentAscii(const std::filesystem::path &file) {
  std::ifstream in(file, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot read Fluent mesh");
  Tokens tokens(in);
  MixedMesh m;
  std::map<std::size_t, std::string> names;
  std::map<std::size_t, int> boundaryCode;
  std::map<std::size_t, std::size_t> faceZones;
  std::set<std::size_t> definedNodes, definedFaces, definedCells;
  std::size_t declaredCells = 0, declaredFaces = 0, declaredNodes = 0;
  Json metadata;
  for (;;) {
    const auto opening = tokens.next();
    if (opening.empty())
      break;
    if (opening != "(")
      throw std::runtime_error("expected Fluent section opening");
    const auto section = number(tokens.next(), 10);
    if (section == 0) {
      const auto comment = tokens.next();
      if (comment.rfind("BabelMesh:", 0) == 0)
        metadata = Json::parse(comment.substr(10));
      if (tokens.next() != ")")
        throw std::runtime_error("invalid Fluent comment");
      continue;
    }
    if (section == 2) {
      if (tokens.next() != "3")
        throw std::runtime_error("Fluent 2D planar input needs explicit "
                                 "extrusion; use a 3D one-layer volume export");
      if (tokens.next() != ")")
        throw std::runtime_error("invalid Fluent dimension");
      continue;
    }
    if (section == 2010 || section == 3010 || section == 2012 ||
        section == 3012 || section == 2013 || section == 3013)
      throw std::runtime_error(
          "binary Fluent is unsupported; export ASCII .msh or CGNS");
    if (section == 58 || section == 59 || section == 61)
      throw std::runtime_error("Fluent refinement trees/nonconformal face "
                               "parents require a conforming exchange export");
    if (section == 10 || section == 12 || section == 13) {
      const auto header = tokens.body();
      if (header.size() < 4)
        throw std::runtime_error("short Fluent mesh header");
      const auto zone = number(header[0]), first = number(header[1]),
                 last = number(header[2]);
      if (first != 1 && !zone)
        throw std::runtime_error("unsupported declaration numbering");
      if (last < first || !first)
        throw std::runtime_error("invalid Fluent index range");
      if (!zone) {
        if (section == 10)
          declaredNodes = last;
        else if (section == 12)
          declaredCells = last;
        else
          declaredFaces = last;
        if (tokens.next() != ")")
          throw std::runtime_error("invalid Fluent declaration");
        continue;
      }
      std::vector<std::string> values;
      auto begin = tokens.next();
      const bool homogeneous =
          section == 12 && header.size() == 5 && number(header[4]) > 0;
      if (begin == "(") {
        tokens.putback(begin);
        values = tokens.body();
        if (tokens.next() != ")")
          throw std::runtime_error("Fluent section not closed");
      } else if (begin != ")" || !homogeneous)
        throw std::runtime_error("missing Fluent connectivity body");
      std::size_t cursor = 0;
      auto take = [&]() {
        if (cursor >= values.size())
          throw std::runtime_error("short Fluent section payload");
        return values[cursor++];
      };
      if (section == 10) {
        if (header.size() != 5 || number(header[4]) != 3)
          throw std::runtime_error("Fluent coordinate section must be 3D");
        m.points_m.resize(std::max(m.points_m.size(), last),
                          {std::numeric_limits<double>::quiet_NaN(), 0, 0});
        for (auto id = first; id <= last; ++id) {
          if (!definedNodes.insert(id).second)
            throw std::runtime_error("duplicate Fluent node");
          for (auto &coordinate : m.points_m[id - 1])
            coordinate = std::stod(take());
        }
      } else if (section == 12) {
        const auto type = header.size() == 5 ? number(header[4]) : 0;
        for (auto id = first; id <= last; ++id) {
          if (!definedCells.insert(id).second)
            throw std::runtime_error("duplicate Fluent cell");
          const auto element = type ? type : number(take());
          if (element != 2 && element != 4 && element != 5 && element != 6 &&
              element != 7)
            throw std::runtime_error("unsupported Fluent volume cell type");
        }
      } else {
        if (header.size() != 5)
          throw std::runtime_error("invalid Fluent face header");
        const auto type = number(header[4]);
        boundaryCode[zone] = static_cast<int>(number(header[3]));
        m.faces.resize(std::max(m.faces.size(), last));
        for (auto id = first; id <= last; ++id) {
          if (!definedFaces.insert(id).second)
            throw std::runtime_error("duplicate Fluent face");
          const auto count = type == 0 || type == 5 ? number(take()) : type;
          if (count < 3 || count > m.points_m.size() + declaredNodes)
            throw std::runtime_error("invalid Fluent face vertex count");
          MixedFace f;
          for (std::size_t j = 0; j < count; ++j)
            f.vertices.push_back(static_cast<std::int32_t>(number(take())) - 1);
          const auto c0 = number(take()), c1 = number(take());
          if (!c0 && !c1)
            throw std::runtime_error(
                "Fluent surface-only face has no volume owner");
          if (c0) {
            std::reverse(f.vertices.begin(), f.vertices.end());
            f.owner = static_cast<std::int32_t>(c0) - 1;
            f.neighbour = static_cast<std::int32_t>(c1) - 1;
          } else {
            f.owner = static_cast<std::int32_t>(c1) - 1;
            f.neighbour = -1;
          }
          m.faces[id - 1] = std::move(f);
          faceZones[id - 1] = zone;
        }
      }
      if (cursor != values.size())
        throw std::runtime_error("unexpected Fluent section payload");
      continue;
    }
    if (section == 39 || section == 45) {
      const auto header = tokens.body();
      if (header.size() >= 3)
        names[number(header[0], 10)] = header[2];
      tokens.body();
      if (tokens.next() != ")")
        throw std::runtime_error("invalid Fluent zone declaration");
      continue;
    }
    // Ignore non-mesh case settings while retaining the original topology.
    int depth = 1;
    while (depth) {
      const auto item = tokens.next();
      if (item.empty())
        throw std::runtime_error("truncated Fluent auxiliary section");
      if (item == "(")
        ++depth;
      else if (item == ")")
        --depth;
    }
  }
  if (m.points_m.size() != definedNodes.size() ||
      m.faces.size() != definedFaces.size() ||
      (declaredNodes && declaredNodes != m.points_m.size()) ||
      (declaredFaces && declaredFaces != m.faces.size()))
    throw std::runtime_error("Fluent mesh has missing nodes or faces");
  std::map<std::size_t, std::int32_t> patches;
  for (std::size_t i = 0; i < m.faces.size(); ++i)
    if (m.faces[i].neighbour < 0) {
      const auto zone = faceZones.at(i);
      if (!patches.count(zone)) {
        patches[zone] = static_cast<std::int32_t>(m.patches.size());
        m.patches.push_back({names.count(zone) ? names.at(zone)
                                               : "zone_" + std::to_string(zone),
                             kindOf(boundaryCode.at(zone))});
      }
      m.faces[i].patch = patches.at(zone);
    }
  if (metadata.is_object()) {
    m.generation_dimension = metadata.value("generation_dimension", 3);
    m.extrusion_axis = metadata.value("extrusion_axis", 2);
    m.extrusion_thickness_m = metadata.value("extrusion_thickness_m", 0.0);
    for (auto &patch : m.patches)
      for (const auto &original : metadata.value("patches", Json::array()))
        if (original.at(0) == patch.name)
          patch.kind = original.at(1);
  }
  reconstructFaceCells(m);
  if ((declaredCells && declaredCells != m.cells.size()) ||
      definedCells.size() != m.cells.size())
    throw std::runtime_error(
        "Fluent cell declaration differs from face connectivity");
  if (metadata.is_object() && metadata.contains("cell_regions")) {
    if (metadata.at("cell_regions").size() != m.cells.size())
      throw std::runtime_error("Fluent region metadata count differs");
    for (std::size_t i = 0; i < m.cells.size(); ++i)
      m.cells[i].region =
          static_cast<CellRegion>(metadata.at("cell_regions")[i].get<int>());
  }
  if (metadata.is_object() && metadata.contains("domain_bounds_constraint"))
    m.metadata["domain_bounds_constraint"] =
        metadata.at("domain_bounds_constraint");
  validateMixedMesh(m);
  m.metadata["generation_region_tags_available"] =
      metadata.is_object() && metadata.contains("cell_regions") &&
      metadata.value("generation_region_tags_available", true);
  m.metadata["conversion_source"] = {
      {"format", "fluent-ascii"}, {"file", file.string()}, {"mesh_only", true}};
  return m;
}
} // namespace babelmesh
