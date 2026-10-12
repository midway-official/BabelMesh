#include "babelmesh/mesh_conversion.hpp"
#include <algorithm>
#include <cgnslib.h>
#include <cstring>
#include <map>
#include <set>
#include <stdexcept>
namespace babelmesh {
namespace {
using Json = nlohmann::json;
void checked(int code) {
  if (code != CG_OK)
    throw std::runtime_error(std::string("CGNS: ") + cg_get_error());
}
struct File {
  int id{};
  explicit File(const std::filesystem::path &p, int mode) {
    checked(cg_open(p.string().c_str(), mode, &id));
  }
  ~File() {
    if (id)
      cg_close(id);
  }
  void close() {
    if (id) {
      checked(cg_close(id));
      id = 0;
    }
  }
};
ElementType_t typeOf(CellType t) {
  switch (t) {
  case CellType::tetrahedron:
    return TETRA_4;
  case CellType::hexahedron:
    return HEXA_8;
  case CellType::pyramid:
    return PYRA_5;
  case CellType::prism:
    return PENTA_6;
  default:
    throw std::runtime_error("canonical CGNS cannot represent a polyhedron");
  }
}
CellType cellType(ElementType_t type) {
  switch (type) {
  case TETRA_4:
    return CellType::tetrahedron;
  case HEXA_8:
    return CellType::hexahedron;
  case PYRA_5:
    return CellType::pyramid;
  case PENTA_6:
    return CellType::prism;
  default:
    throw std::runtime_error(
        "CGNS supports first-order tetra/pyramid/prism/hex or NGON/NFACE only");
  }
}
BCType_t bcType(const std::string &s) {
  return s == "wall"       ? BCWall
         : s == "inlet"    ? BCInflow
         : s == "outlet"   ? BCOutflow
         : s == "symmetry" ? BCSymmetryPlane
                           : BCTypeUserDefined;
}
std::string patchKind(BCType_t t) {
  return t == BCWall || t == BCWallViscous || t == BCWallInviscid ||
                 t == BCWallViscousHeatFlux || t == BCWallViscousIsothermal
             ? "wall"
         : t == BCInflow || t == BCInflowSubsonic || t == BCInflowSupersonic
             ? "inlet"
         : t == BCOutflow || t == BCOutflowSubsonic || t == BCOutflowSupersonic
             ? "outlet"
         : t == BCSymmetryPlane ? "symmetry"
                                : "generic";
}
struct Section {
  ElementType_t type{};
  cgsize_t first{}, last{};
  std::vector<cgsize_t> connectivity, offsets, parent;
};
} // namespace
void writeCgns(const MixedMesh &m, const std::filesystem::path &path,
               bool forcePolyhedral) {
  const bool poly =
      forcePolyhedral ||
      std::any_of(m.cells.begin(), m.cells.end(),
                  [](const auto &c) { return c.type == CellType::polyhedron; });
  File file(path, CG_MODE_WRITE);
  const int fn = file.id;
  int B = 0, Z = 0, S = 0;
  checked(cg_base_write(fn, "BabelMesh", 3, 3, &B));
  checked(cg_goto(fn, B, "end"));
  checked(cg_units_write(Kilogram, Meter, Second, Kelvin, Radian));
  checked(cg_dataclass_write(Dimensional));
  const cgsize_t size[] = {static_cast<cgsize_t>(m.points_m.size()),
                           static_cast<cgsize_t>(m.cells.size()), 0};
  checked(cg_zone_write(fn, B, "Fluid", size, Unstructured, &Z));
  constexpr const char *names[] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
  for (std::size_t a = 0; a < 3; ++a) {
    std::vector<double> values;
    values.reserve(m.points_m.size());
    for (const auto &p : m.points_m)
      values.push_back(p[a]);
    int C = 0;
    checked(cg_coord_write(fn, B, Z, RealDouble, names[a], values.data(), &C));
    checked(cg_goto(fn, B, "Zone_t", Z, "GridCoordinates_t", 1, "DataArray_t",
                    C, "end"));
    const float exponents[5] = {0, 1, 0, 0, 0};
    checked(cg_exponents_write(RealSingle, exponents));
  }
  Json metadata = {{"generation_dimension", m.generation_dimension},
                   {"generation_region_tags_available",
                    m.metadata.value("generation_region_tags_available", true)},
                   {"extrusion_axis", m.extrusion_axis},
                   {"extrusion_thickness_m", m.extrusion_thickness_m},
                   {"cell_regions", Json::array()},
                   {"patch_kinds", Json::object()}};
  for (const auto &cell : m.cells)
    metadata["cell_regions"].push_back(static_cast<int>(cell.region));
  for (const auto &patch : m.patches)
    metadata["patch_kinds"][patch.name] = patch.kind;
  if (m.metadata.contains("domain_bounds_constraint"))
    metadata["domain_bounds_constraint"] =
        m.metadata.at("domain_bounds_constraint");
  checked(cg_goto(fn, B, "Zone_t", Z, "end"));
  checked(cg_descriptor_write("BabelMesh", metadata.dump().c_str()));
  std::map<std::size_t, cgsize_t> boundaryIds;
  if (poly) {
    std::vector<cgsize_t> faces, faceOffsets = {0};
    for (std::size_t i = 0; i < m.faces.size(); ++i) {
      for (auto node : m.faces[i].vertices)
        faces.push_back(node + 1);
      faceOffsets.push_back(static_cast<cgsize_t>(faces.size()));
      boundaryIds[i] = static_cast<cgsize_t>(i + 1);
    }
    checked(cg_poly_section_write(fn, B, Z, "Faces", NGON_n, 1,
                                  static_cast<cgsize_t>(m.faces.size()), 0,
                                  faces.data(), faceOffsets.data(), &S));
    std::vector<std::vector<cgsize_t>> incidence(m.cells.size());
    for (std::size_t i = 0; i < m.faces.size(); ++i) {
      const auto &f = m.faces[i];
      incidence.at(static_cast<std::size_t>(f.owner))
          .push_back(static_cast<cgsize_t>(i + 1));
      if (f.neighbour >= 0)
        incidence.at(static_cast<std::size_t>(f.neighbour))
            .push_back(-static_cast<cgsize_t>(i + 1));
    }
    std::vector<cgsize_t> cells, cellOffsets = {0};
    for (const auto &facesOfCell : incidence) {
      cells.insert(cells.end(), facesOfCell.begin(), facesOfCell.end());
      cellOffsets.push_back(static_cast<cgsize_t>(cells.size()));
    }
    checked(cg_poly_section_write(
        fn, B, Z, "Cells", NFACE_n, static_cast<cgsize_t>(m.faces.size() + 1),
        static_cast<cgsize_t>(m.faces.size() + m.cells.size()), 0, cells.data(),
        cellOffsets.data(), &S));
  } else {
    std::vector<cgsize_t> elements, offsets = {0};
    for (const auto &c : m.cells) {
      elements.push_back(typeOf(c.type));
      for (auto node : c.vertices)
        elements.push_back(node + 1);
      offsets.push_back(static_cast<cgsize_t>(elements.size()));
    }
    checked(cg_poly_section_write(fn, B, Z, "Cells", MIXED, 1,
                                  static_cast<cgsize_t>(m.cells.size()), 0,
                                  elements.data(), offsets.data(), &S));
    cgsize_t next = static_cast<cgsize_t>(m.cells.size() + 1);
    for (std::size_t patch = 0; patch < m.patches.size(); ++patch) {
      // Homogeneous boundary sections are understood by older CGNS consumers
      // as well as readers of variable-length MIXED volume sections.
      for (std::size_t arity : {std::size_t{3}, std::size_t{4}}) {
        std::vector<cgsize_t> elements;
        const cgsize_t begin = next;
        for (std::size_t i = 0; i < m.faces.size(); ++i)
          if (m.faces[i].patch == static_cast<std::int32_t>(patch) &&
              m.faces[i].vertices.size() == arity) {
            for (auto node : m.faces[i].vertices)
              elements.push_back(node + 1);
            boundaryIds[i] = next++;
          }
        if (next > begin)
          checked(cg_section_write(fn, B, Z,
                                   ("Boundary_" + std::to_string(patch) +
                                    (arity == 3 ? "_tri" : "_quad"))
                                       .c_str(),
                                   arity == 3 ? TRI_3 : QUAD_4, begin, next - 1,
                                   0, elements.data(), &S));
      }
    }
  }
  for (std::size_t patch = 0; patch < m.patches.size(); ++patch) {
    const auto &p = m.patches[patch];
    if (p.name.size() > 32)
      throw std::runtime_error("CGNS boundary name exceeds 32 characters: " +
                               p.name);
    std::vector<cgsize_t> ids;
    for (std::size_t i = 0; i < m.faces.size(); ++i)
      if (m.faces[i].patch == static_cast<std::int32_t>(patch))
        ids.push_back(boundaryIds.at(i));
    if (ids.empty())
      continue;
    int BC = 0;
    checked(cg_boco_write(fn, B, Z, p.name.c_str(), bcType(p.kind), PointList,
                          static_cast<cgsize_t>(ids.size()), ids.data(), &BC));
    checked(cg_boco_gridlocation_write(fn, B, Z, BC, FaceCenter));
  }
  file.close();
}
MixedMesh readCgns(const std::filesystem::path &path) {
  File file(path, CG_MODE_READ);
  const int fn = file.id;
  int bases = 0;
  checked(cg_nbases(fn, &bases));
  if (bases != 1)
    throw std::runtime_error("CGNS converter requires one base; combine the "
                             "mesh bases in the exporting application");
  char name[33]{};
  int cellDim = 0, physicalDim = 0;
  checked(cg_base_read(fn, 1, name, &cellDim, &physicalDim));
  if (cellDim != 3 || physicalDim != 3)
    throw std::runtime_error("CGNS input must be a 3D volume mesh; planar 2D "
                             "input needs a one-layer extrusion export");
  int zones = 0;
  checked(cg_nzones(fn, 1, &zones));
  if (zones != 1)
    throw std::runtime_error("CGNS converter requires one conforming volume "
                             "zone; multizone/overset/periodic connectivity "
                             "needs an explicit region/interface mapping");
  ZoneType_t zoneType{};
  checked(cg_zone_type(fn, 1, 1, &zoneType));
  if (zoneType != Unstructured)
    throw std::runtime_error("CGNS structured input must be exported as an "
                             "unstructured volume zone");
  cgsize_t size[9]{};
  checked(cg_zone_read(fn, 1, 1, name, size));
  const std::string zoneName(name);
  if (size[0] <= 0 || size[1] <= 0)
    throw std::runtime_error("empty CGNS volume mesh");
  MixedMesh m;
  m.points_m.resize(static_cast<std::size_t>(size[0]));
  double unitScale = 1;
  std::string unit = "unspecified; assumed m";
  // CGNS units inherit from the nearest ancestor. Coordinates can override
  // Zone_t/Base_t units at GridCoordinates_t or DataArray_t.
  auto coordinateUnit = [&]() {
    MassUnits_t mass{};
    LengthUnits_t length{};
    TimeUnits_t time{};
    TemperatureUnits_t temperature{};
    AngleUnits_t angle{};
    const int code = cg_units_read(&mass, &length, &time, &temperature, &angle);
    if (code == CG_NODE_NOT_FOUND)
      return false;
    checked(code);
    switch (length) {
    case Meter:
      unitScale = 1;
      break;
    case Centimeter:
      unitScale = .01;
      break;
    case Millimeter:
      unitScale = .001;
      break;
    case Foot:
      unitScale = .3048;
      break;
    case Inch:
      unitScale = .0254;
      break;
    case LengthUnitsNull:
    case LengthUnitsUserDefined:
      return false;
    default:
      throw std::runtime_error("unsupported CGNS coordinate length unit");
    }
    unit = cg_LengthUnitsName(length);
    return true;
  };
  checked(cg_goto(fn, 1, "end"));
  coordinateUnit();
  checked(cg_goto(fn, 1, "Zone_t", 1, "end"));
  coordinateUnit();
  checked(cg_goto(fn, 1, "Zone_t", 1, "GridCoordinates_t", 1, "end"));
  coordinateUnit();
  const double inheritedScale = unitScale;
  const std::string inheritedUnit = unit;
  int coordinateCount = 0;
  checked(cg_ncoords(fn, 1, 1, &coordinateCount));
  const cgsize_t first = 1, last = size[0];
  constexpr const char *names[] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
  for (std::size_t a = 0; a < 3; ++a) {
    unitScale = inheritedScale;
    unit = inheritedUnit;
    int array = 0;
    for (int c = 1; c <= coordinateCount; ++c) {
      DataType_t type{};
      char label[33]{};
      checked(cg_coord_info(fn, 1, 1, c, &type, label));
      if (std::strcmp(label, names[a]) == 0)
        array = c;
    }
    if (!array)
      throw std::runtime_error("missing CGNS Cartesian coordinate array");
    checked(cg_goto(fn, 1, "Zone_t", 1, "GridCoordinates_t", 1, "DataArray_t",
                    array, "end"));
    coordinateUnit();
    std::vector<double> coordinates(m.points_m.size());
    checked(cg_coord_read(fn, 1, 1, names[a], RealDouble, &first, &last,
                          coordinates.data()));
    for (std::size_t i = 0; i < coordinates.size(); ++i)
      m.points_m[i][a] = coordinates[i] * unitScale;
  }
  Json metadata = Json::object();
  checked(cg_goto(fn, 1, "Zone_t", 1, "end"));
  int descriptors = 0;
  checked(cg_ndescriptors(&descriptors));
  for (int D = 1; D <= descriptors; ++D) {
    char *text = nullptr;
    char label[33]{};
    checked(cg_descriptor_read(D, label, &text));
    if (std::strcmp(label, "BabelMesh") == 0)
      metadata = Json::parse(text);
    cg_free(text);
  }
  int sectionCount = 0;
  checked(cg_nsections(fn, 1, 1, &sectionCount));
  std::map<cgsize_t, FaceKey> faceById;
  std::map<cgsize_t, MixedCell> cellById;
  std::map<cgsize_t, std::vector<cgsize_t>> polyCells;
  for (int S = 1; S <= sectionCount; ++S) {
    Section section;
    int nbndry = 0, parents = 0;
    char label[33]{};
    checked(cg_section_read(fn, 1, 1, S, label, &section.type, &section.first,
                            &section.last, &nbndry, &parents));
    cgsize_t length = 0;
    checked(cg_ElementDataSize(fn, 1, 1, S, &length));
    if (length <= 0 || section.last < section.first)
      throw std::runtime_error("invalid CGNS connectivity section");
    section.connectivity.resize(static_cast<std::size_t>(length));
    const auto count =
        static_cast<std::size_t>(section.last - section.first + 1);
    if (parents)
      section.parent.resize(count * 4);
    if (section.type == NGON_n || section.type == NFACE_n) {
      section.offsets.resize(count + 1);
      checked(cg_poly_elements_read(fn, 1, 1, S, section.connectivity.data(),
                                    section.offsets.data(),
                                    parents ? section.parent.data() : nullptr));
      for (std::size_t i = 0; i < count; ++i) {
        const auto begin = section.offsets.at(i),
                   end = section.offsets.at(i + 1);
        if (begin < 0 || end < begin ||
            static_cast<std::size_t>(end) > section.connectivity.size())
          throw std::runtime_error("invalid CGNS poly offsets");
        std::vector<cgsize_t> ring(section.connectivity.begin() + begin,
                                   section.connectivity.begin() + end);
        const auto id = section.first + static_cast<cgsize_t>(i);
        if (section.type == NGON_n) {
          FaceKey face;
          for (auto node : ring) {
            if (node <= 0 || node > size[0])
              throw std::runtime_error("invalid CGNS node");
            face.push_back(node - 1);
          }
          if (!faceById.emplace(id, std::move(face)).second)
            throw std::runtime_error("duplicate CGNS face id");
        } else if (!polyCells.emplace(id, std::move(ring)).second)
          throw std::runtime_error("duplicate CGNS NFACE id");
      }
      continue;
    }
    if (section.type == MIXED) {
      section.offsets.resize(count + 1);
      checked(cg_poly_elements_read(fn, 1, 1, S, section.connectivity.data(),
                                    section.offsets.data(),
                                    parents ? section.parent.data() : nullptr));
    } else
      checked(cg_elements_read(fn, 1, 1, S, section.connectivity.data(),
                               parents ? section.parent.data() : nullptr));
    std::size_t offset = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const auto type =
          section.type == MIXED
              ? static_cast<ElementType_t>(section.connectivity.at(offset++))
              : section.type;
      int npe = 0;
      checked(cg_npe(type, &npe));
      if (npe <= 0 ||
          offset + static_cast<std::size_t>(npe) > section.connectivity.size())
        throw std::runtime_error("invalid CGNS standard cell connectivity");
      FaceKey ring;
      for (int j = 0; j < npe; ++j) {
        const auto node = section.connectivity.at(offset++);
        if (node <= 0 || node > size[0])
          throw std::runtime_error("invalid CGNS node index");
        ring.push_back(node - 1);
      }
      const auto id = section.first + static_cast<cgsize_t>(i);
      if (type == TRI_3 || type == QUAD_4)
        faceById[id] = ring;
      else {
        const MixedCell cell{cellType(type), CellRegion::transition, ring};
        if (!cellById.emplace(id, cell).second)
          throw std::runtime_error("duplicate CGNS volume id");
      }
    }
    if (offset != section.connectivity.size())
      throw std::runtime_error("trailing CGNS element connectivity");
  }
  if (!polyCells.empty()) {
    if (!cellById.empty())
      throw std::runtime_error("CGNS mixed canonical/NFACE sections require a "
                               "unified cell representation");
    for (const auto &entry : polyCells) {
      std::set<std::int32_t> vertices;
      std::vector<FaceKey> faces;
      for (auto signedId : entry.second) {
        if (signedId == 0 || signedId == std::numeric_limits<cgsize_t>::min())
          throw std::runtime_error("invalid signed CGNS face id");
        auto face = faceById.at(std::abs(signedId));
        if (signedId < 0)
          std::reverse(face.begin(), face.end());
        vertices.insert(face.begin(), face.end());
        faces.push_back(std::move(face));
      }
      cellById.emplace(entry.first,
                       MixedCell{CellType::polyhedron, CellRegion::transition,
                                 FaceKey(vertices.begin(), vertices.end()),
                                 std::move(faces)});
    }
  }
  if (cellById.size() != static_cast<std::size_t>(size[1]))
    throw std::runtime_error("CGNS cell count differs from volume sections; "
                             "NGON without NFACE is unsupported");
  for (const auto &c : cellById)
    m.cells.push_back(c.second);
  BoundaryMap boundary;
  int bcCount = 0;
  checked(cg_nbocos(fn, 1, 1, &bcCount));
  for (int BC = 1; BC <= bcCount; ++BC) {
    BCType_t type{};
    PointSetType_t pointType{};
    cgsize_t npnts = 0, normalSize = 0;
    int normalIndex[3]{}, datasets = 0;
    DataType_t normalType{};
    char label[33]{};
    checked(cg_boco_info(fn, 1, 1, BC, label, &type, &pointType, &npnts,
                         normalIndex, &normalSize, &normalType, &datasets));
    GridLocation_t location = Vertex;
    const int locationCode = cg_boco_gridlocation_read(fn, 1, 1, BC, &location);
    if (locationCode != CG_OK)
      location = Vertex;
    if (location != FaceCenter && pointType != ElementList &&
        pointType != ElementRange)
      throw std::runtime_error("CGNS BC must use face/element ids; vertex BC "
                               "needs explicit face mapping");
    if (pointType != PointList && pointType != PointRange &&
        pointType != ElementList && pointType != ElementRange)
      throw std::runtime_error("unsupported CGNS BC point set");
    std::vector<cgsize_t> ids(static_cast<std::size_t>(npnts));
    std::vector<double> normals(static_cast<std::size_t>(normalSize) * 2);
    checked(cg_boco_read(fn, 1, 1, BC, ids.data(),
                         normalSize ? normals.data() : nullptr));
    if (pointType == PointRange || pointType == ElementRange) {
      if (ids.size() != 2 || ids[1] < ids[0])
        throw std::runtime_error("invalid CGNS BC range");
      const auto end = ids[1], begin = ids[0];
      ids.clear();
      for (auto id = begin; id <= end; ++id)
        ids.push_back(id);
    }
    const auto patch = static_cast<std::int32_t>(m.patches.size());
    std::string kind = patchKind(type);
    if (metadata.is_object() &&
        metadata.value("patch_kinds", Json::object()).contains(label))
      kind = metadata.at("patch_kinds").at(label);
    m.patches.push_back({label, kind});
    for (auto id : ids) {
      const auto face = faceById.find(id);
      if (face == faceById.end())
        throw std::runtime_error("CGNS BC references a missing face");
      if (!boundary.emplace(faceKey(face->second), patch).second)
        throw std::runtime_error(
            "CGNS boundary face belongs to more than one BC");
    }
  }
  if (metadata.is_object()) {
    m.generation_dimension = metadata.value("generation_dimension", 3);
    m.extrusion_axis = metadata.value("extrusion_axis", 2);
    m.extrusion_thickness_m = metadata.value("extrusion_thickness_m", 0.0);
    // BabelMesh's descriptor stores explicit meters, independently of a
    // CGNS coordinate array's inherited length unit.
    if (metadata.contains("domain_bounds_constraint"))
      m.metadata["domain_bounds_constraint"] =
          metadata.at("domain_bounds_constraint");
    if (metadata.contains("cell_regions")) {
      if (metadata.at("cell_regions").size() != m.cells.size())
        throw std::runtime_error("CGNS region metadata count differs");
      for (std::size_t i = 0; i < m.cells.size(); ++i)
        m.cells[i].region =
            static_cast<CellRegion>(metadata.at("cell_regions")[i].get<int>());
    }
  }
  buildMixedFaces(m, boundary);
  if (!polyCells.empty()) {
    const auto regions = metadata.value("cell_regions", Json::array());
    reconstructFaceCells(m);
    if (regions.size() == m.cells.size())
      for (std::size_t i = 0; i < m.cells.size(); ++i)
        m.cells[i].region = static_cast<CellRegion>(regions[i].get<int>());
  }
  m.metadata["generation_region_tags_available"] =
      metadata.is_object() && metadata.contains("cell_regions") &&
      metadata.value("generation_region_tags_available", true);
  m.metadata["conversion_source"] = {{"format", "cgns"},
                                     {"file", path.string()},
                                     {"zone", zoneName},
                                     {"mesh_only", true},
                                     {"coordinate_unit", inheritedUnit},
                                     {"coordinates_normalized_to", "m"}};
  file.close();
  return m;
}
} // namespace babelmesh
