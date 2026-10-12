#include "babelmesh/hybrid_mesher.hpp"
#include "babelmesh/background_grid.hpp"
#include "babelmesh/fluid_domain.hpp"
#include "babelmesh/mesh_optimizer.hpp"
#include "babelmesh/surface_distance.hpp"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <Eigen/Dense>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <gmsh.h>
#include <gp_Pln.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>

namespace babelmesh {
namespace {
using Json = nlohmann::json;
using Vec = Eigen::Vector3d;
using Index = std::int32_t;
struct SkinFace {
  FaceKey vertices;
  Index patch{};
};
struct PlanarCell {
  FaceKey vertices;
  CellRegion region{};
};
struct Loop {
  FaceKey vertices;
  std::vector<Index> patches;
  int parent{-1};
  int depth{};
};
std::mutex gmshMutex;

Vec vec(const std::array<double, 3> &p) { return {p[0], p[1], p[2]}; }
std::array<double, 3> arr(const Vec &p) { return {p.x(), p.y(), p.z()}; }
double scaleOf(const Json &p) {
  const auto unit = p.at("units").at("length").get<std::string>();
  if (unit == "m")
    return 1;
  if (unit == "mm")
    return .001;
  if (unit == "cm")
    return .01;
  if (unit == "inch")
    return .0254;
  throw std::runtime_error("unsupported project length unit");
}
double positive(const Json &object, const char *key) {
  const auto &v = object.at(key);
  if (!v.is_number() || !std::isfinite(v.get<double>()) || v.get<double>() <= 0)
    throw std::runtime_error(std::string(key) + " must be finite and positive");
  return v.get<double>();
}
void keys(const Json &p, const std::set<std::string> &allowed,
          const std::string &where) {
  if (!p.is_object())
    throw std::runtime_error(where + " must be an object");
  for (auto it = p.begin(); it != p.end(); ++it)
    if (!allowed.count(it.key()))
      throw std::runtime_error(where + " has unknown key '" + it.key() + "'");
}
struct DomainQuery {
  std::vector<std::unique_ptr<BRepClass3d_SolidClassifier>> solids;
  explicit DomainQuery(const TopoDS_Shape &s) {
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next())
      solids.push_back(std::make_unique<BRepClass3d_SolidClassifier>(
          TopoDS::Solid(ex.Current())));
  }
  bool inside(const Vec &p, double tolerance = 1e-9) {
    for (auto &s : solids) {
      s->Perform(gp_Pnt(p.x(), p.y(), p.z()), tolerance);
      if (s->State() == TopAbs_IN || s->State() == TopAbs_ON)
        return true;
    }
    return false;
  }
};
struct Session {
  Session() {
    gmsh::initialize(0, nullptr, false);
    gmsh::option::setNumber("General.Terminal", 1);
  }
  ~Session() {
    try {
      gmsh::finalize();
    } catch (...) {
    }
  }
};
Index addPoint(MixedMesh &m, const Vec &p) {
  if (m.points_m.size() >=
      static_cast<std::size_t>(std::numeric_limits<Index>::max()))
    throw std::runtime_error("32-bit point index budget exceeded");
  m.points_m.push_back(arr(p));
  return static_cast<Index>(m.points_m.size() - 1);
}
Index addPatch(MixedMesh &m, const std::string &name, const std::string &kind) {
  for (std::size_t i = 0; i < m.patches.size(); ++i)
    if (m.patches[i].name == name) {
      if (m.patches[i].kind != kind)
        throw std::runtime_error("patch kind conflict for " + name);
      return static_cast<Index>(i);
    }
  m.patches.push_back({name, kind});
  return static_cast<Index>(m.patches.size() - 1);
}
void boundaryFace(BoundaryMap &b, const FaceKey &ring, Index patch) {
  if (!b.emplace(faceKey(ring), patch).second)
    throw std::runtime_error("duplicate physical boundary face");
}
Vec faceCentre(const MixedMesh &m, const FaceKey &ring) {
  Vec c = Vec::Zero();
  for (auto i : ring)
    c += vec(m.points_m.at(static_cast<std::size_t>(i)));
  return c / static_cast<double>(ring.size());
}
double signedTet(const MixedMesh &m, const FaceKey &v) {
  const Vec a = vec(m.points_m.at(static_cast<std::size_t>(v[0])));
  return (vec(m.points_m.at(static_cast<std::size_t>(v[1]))) - a)
             .dot(
                 (vec(m.points_m.at(static_cast<std::size_t>(v[2]))) - a)
                     .cross(vec(m.points_m.at(static_cast<std::size_t>(v[3]))) -
                            a)) /
         6;
}
void addTet(MixedMesh &m, FaceKey v, CellRegion r) {
  if (signedTet(m, v) < 0)
    std::swap(v[1], v[2]);
  if (!(signedTet(m, v) > 0))
    throw std::runtime_error("boundary layer creates a degenerate tetrahedron");
  m.cells.push_back({CellType::tetrahedron, r, std::move(v)});
}

void options(const Json &p, double h, const std::string &phase) {
  gmsh::option::setNumber("General.NumThreads",
                          p.at("resources").at("threads").get<int>());
  gmsh::option::setNumber("Mesh.MeshSizeMin",
                          positive(p.at("sizing"), "minimum_size") *
                              scaleOf(p));
  gmsh::option::setNumber("Mesh.MeshSizeMax", h);
  gmsh::option::setNumber("Mesh.ElementOrder", 1);
  // Sparse node tags identify fixed background/layer nodes across adapters.
  // Automatic renumbering would overwrite retained core coordinates.
  gmsh::option::setNumber("Mesh.Renumber", 0);
  gmsh::option::setNumber("Mesh.Algorithm", 6);
  gmsh::option::setNumber("Mesh.Algorithm3D", 1);
  gmsh::option::setNumber("Mesh.Optimize", 1);
  gmsh::option::setNumber("Mesh.RandomSeed", 1);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  std::vector<double> fields;
  const int uniform = gmsh::model::mesh::field::add("MathEval");
  std::ostringstream constant;
  constant << std::setprecision(17) << h;
  gmsh::model::mesh::field::setString(uniform, "F", constant.str());
  fields.push_back(uniform);
  const double scale = scaleOf(p);
  for (const auto &source : p.at("sizing").at("sources")) {
    const auto apply = source.value("apply_to", "all");
    if (apply != "all" && apply != phase)
      continue;
    const auto kind = source.at("kind").get<std::string>();
    const double target = source.at("target_size").get<double>() * scale;
    const auto &region = source.at("region");
    int field = -1;
    if (kind == "box") {
      field = gmsh::model::mesh::field::add("Box");
      gmsh::model::mesh::field::setNumber(field, "VIn", target);
      gmsh::model::mesh::field::setNumber(field, "VOut", h);
      constexpr const char *low[] = {"XMin", "YMin", "ZMin"};
      constexpr const char *high[] = {"XMax", "YMax", "ZMax"};
      for (std::size_t a = 0; a < 3; ++a) {
        gmsh::model::mesh::field::setNumber(
            field, low[a],
            region.at("bounds").at("minimum")[a].get<double>() * scale);
        gmsh::model::mesh::field::setNumber(
            field, high[a],
            region.at("bounds").at("maximum")[a].get<double>() * scale);
      }
      gmsh::model::mesh::field::setNumber(
          field, "Thickness", source.value("transition_width", 0.0) * scale);
    } else if (kind == "sphere") {
      std::ostringstream formula;
      formula << std::setprecision(17);
      const auto &c = region.at("center");
      const double width =
          std::max(source.value("transition_width", 0.0) * scale, h * .05);
      formula << target << "+(" << h - target << ")*Min(1,Max(0,(Sqrt((x-("
              << c[0].get<double>() * scale << "))^2+(y-("
              << c[1].get<double>() * scale << "))^2+(z-("
              << c[2].get<double>() * scale << "))^2)-"
              << region.at("radius").get<double>() * scale << ")/" << width
              << "))";
      field = gmsh::model::mesh::field::add("MathEval");
      gmsh::model::mesh::field::setString(field, "F", formula.str());
    } else
      throw std::runtime_error("hybrid surface/transition sizing currently "
                               "supports box and sphere regions; source " +
                               source.at("id").get<std::string>() +
                               " is unsupported");
    fields.push_back(field);
  }
  const int combine = gmsh::model::mesh::field::add("Min");
  gmsh::model::mesh::field::setNumbers(combine, "FieldsList", fields);
  gmsh::model::mesh::field::setAsBackgroundMesh(combine);
}
std::map<std::size_t, Index> importNodes(MixedMesh &m) {
  std::vector<std::size_t> tags;
  std::vector<double> xyz, param;
  gmsh::model::mesh::getNodes(tags, xyz, param);
  std::map<std::size_t, Index> map;
  for (std::size_t i = 0; i < tags.size(); ++i)
    map.emplace(tags[i],
                addPoint(m, Vec(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2])));
  return map;
}

Index patchAt(MixedMesh &m, const CadModel &cad, const Vec &p,
              const Json &project, double h) {
  const double scale = scaleOf(project), tol = std::max(1e-8, h * .1);
  double nearest = std::numeric_limits<double>::max();
  const PhysicalGroup *chosen = nullptr;
  const auto vertex =
      BRepBuilderAPI_MakeVertex(gp_Pnt(p.x(), p.y(), p.z())).Vertex();
  for (const auto &group : cad.groups()) {
    double distance = std::numeric_limits<double>::max();
    for (const auto &s : group.selected_subshapes) {
      BRepExtrema_DistShapeShape q(vertex, s);
      q.Perform();
      if (q.IsDone())
        distance = std::min(distance, q.Value());
    }
    if (distance < tol && distance < nearest - 1e-10) {
      nearest = distance;
      chosen = &group;
    } else if (distance < tol && std::abs(distance - nearest) <= 1e-10 &&
               chosen && chosen->name != group.name)
      throw std::runtime_error("overlapping physical groups '" + chosen->name +
                               "' and '" + group.name + "'");
  }
  if (chosen)
    return addPatch(m, chosen->name, chosen->kind);
  const auto &b = project.at("background").at("bounds");
  constexpr const char *names[] = {"xmin", "xmax", "ymin",
                                   "ymax", "zmin", "zmax"};
  for (int axis = 0; axis < 3; ++axis)
    for (int side = 0; side < 2; ++side) {
      const double value = b.at(side ? "maximum" : "minimum")
                               .at(static_cast<std::size_t>(axis))
                               .get<double>() *
                           scale;
      if (std::abs(p[axis] - value) < 1e-8)
        return addPatch(m, names[2 * axis + side], "generic");
    }
  return addPatch(m, "wall", "wall");
}

TopoDS_Shape sectionShape(const TopoDS_Shape &fluid, int axis,
                          double coordinate, const Json &p) {
  Vec origin = Vec::Zero(), normal = Vec::Zero();
  origin[axis] = coordinate;
  normal[axis] = 1;
  const double scale = scaleOf(p);
  double span = 0;
  for (int i = 0; i < 3; ++i)
    span = std::max(
        span,
        std::max(
            std::abs(
                p.at("background").at("bounds").at("minimum")[i].get<double>() *
                scale),
            std::abs(
                p.at("background").at("bounds").at("maximum")[i].get<double>() *
                scale)));
  const double radius = 2 * (span + coordinate * coordinate + 1);
  const auto plane = BRepBuilderAPI_MakeFace(
                         gp_Pln(gp_Pnt(origin.x(), origin.y(), origin.z()),
                                gp_Dir(normal.x(), normal.y(), normal.z())),
                         -radius, radius, -radius, radius)
                         .Face();
  BRepAlgoAPI_Common common(fluid, plane);
  common.Build();
  if (!common.IsDone() || common.Shape().IsNull())
    throw std::runtime_error("CAD planar section failed");
  bool face = false;
  for (TopExp_Explorer ex(common.Shape(), TopAbs_FACE); ex.More(); ex.Next())
    face = true;
  if (!face)
    throw std::runtime_error("selected plane does not cut a fluid area");
  return common.Shape();
}

std::vector<SkinFace> cadSkin(MixedMesh &mesh, const CadModel &cad,
                              const Json &p, double h, DomainQuery &query) {
  gmsh::model::add("cad-skin");
  gmsh::vectorpair entities;
  TopoDS_Shape shape = cad.fluidDomain();
  const auto &cfg = p.at("meshing");
  const double scale = scaleOf(p);
  if (mesh.generation_dimension == 2) {
    const double coordinate =
        cfg.at("section_coordinate").get<double>() * scale;
    shape = sectionShape(shape, mesh.extrusion_axis, coordinate, p);
  }
  // Exchange B-Rep through the public file API to avoid coupling OCCT ABI
  // versions. The CLI owns a per-job temporary path, removed after import.
  const auto stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  const auto brep = std::filesystem::temp_directory_path() /
                    ("babelmesh-skin-" + std::to_string(stamp) + ".brep");
  if (!BRepTools::Write(shape, brep.string().c_str()))
    throw std::runtime_error("cannot serialize CAD for surface mesher");
  try {
    gmsh::model::occ::importShapes(brep.string(), entities, true, "brep");
  } catch (...) {
    std::filesystem::remove(brep);
    throw;
  }
  std::filesystem::remove(brep);
  gmsh::model::occ::synchronize();
  options(p, h, "surface");
  gmsh::model::mesh::generate(mesh.generation_dimension - 1);
  const auto nodes = importNodes(mesh);
  std::vector<SkinFace> result;
  gmsh::vectorpair skin;
  gmsh::model::getEntities(skin, mesh.generation_dimension - 1);
  for (const auto &entity : skin) {
    std::vector<int> types;
    std::vector<std::vector<std::size_t>> ids, connectivity;
    gmsh::model::mesh::getElements(types, ids, connectivity, entity.first,
                                   entity.second);
    for (std::size_t t = 0; t < types.size(); ++t) {
      const int expected = mesh.generation_dimension == 3 ? 2 : 1;
      if (types[t] != expected)
        throw std::runtime_error(
            "surface backend returned non-linear or non-simplex elements");
      const auto size = static_cast<std::size_t>(mesh.generation_dimension);
      for (std::size_t j = 0; j < connectivity[t].size(); j += size) {
        FaceKey ring;
        for (std::size_t k = 0; k < size; ++k)
          ring.push_back(nodes.at(connectivity[t][j + k]));
        const Vec c = faceCentre(mesh, ring);
        const Index patch = patchAt(mesh, cad, c, p, h);
        if (mesh.generation_dimension == 3) {
          const Vec a = vec(mesh.points_m[static_cast<std::size_t>(ring[0])]);
          Vec n =
              (vec(mesh.points_m[static_cast<std::size_t>(ring[1])]) - a)
                  .cross(vec(mesh.points_m[static_cast<std::size_t>(ring[2])]) -
                         a)
                  .normalized();
          const auto projection = projectPointToFluidBoundary(cad, arr(c));
          if (query.inside(vec(projection.point_m) +
                           n * std::max(1e-7, h * 1e-4)))
            std::swap(ring[1], ring[2]);
        }
        result.push_back({std::move(ring), patch});
      }
    }
  }
  return result;
}

std::vector<Loop> loopsFromEdges(const std::vector<SkinFace> &edges) {
  std::map<Index, std::vector<std::size_t>> incidence;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].vertices.size() != 2)
      throw std::runtime_error("invalid planar boundary edge");
    for (auto v : edges[i].vertices)
      incidence[v].push_back(i);
  }
  for (const auto &node : incidence)
    if (node.second.size() != 2)
      throw std::runtime_error(
          "planar boundary is open or non-manifold at node " +
          std::to_string(node.first));
  std::vector<bool> used(edges.size());
  std::vector<Loop> loops;
  for (std::size_t start = 0; start < edges.size(); ++start)
    if (!used[start]) {
      Loop loop;
      Index current = edges[start].vertices[0], first = current;
      std::size_t edge = start;
      do {
        if (used[edge])
          throw std::runtime_error("planar boundary loop repeats an edge");
        used[edge] = true;
        loop.vertices.push_back(current);
        loop.patches.push_back(edges[edge].patch);
        current = edges[edge].vertices[0] == current ? edges[edge].vertices[1]
                                                     : edges[edge].vertices[0];
        const auto &next = incidence.at(current);
        edge = next[0] == edge ? next[1] : next[0];
      } while (current != first);
      if (loop.vertices.size() < 3)
        throw std::runtime_error("degenerate planar loop");
      loops.push_back(std::move(loop));
    }
  return loops;
}
bool inPolygon(const Vec &p, const Loop &loop, const MixedMesh &m, int u,
               int v) {
  bool in = false;
  for (std::size_t i = 0, j = loop.vertices.size() - 1;
       i < loop.vertices.size(); j = i++) {
    const Vec a = vec(m.points_m[static_cast<std::size_t>(loop.vertices[i])]),
              b = vec(m.points_m[static_cast<std::size_t>(loop.vertices[j])]);
    if ((a[v] > p[v]) != (b[v] > p[v]) &&
        p[u] < (b[u] - a[u]) * (p[v] - a[v]) / (b[v] - a[v]) + a[u])
      in = !in;
  }
  return in;
}
void nestLoops(std::vector<Loop> &loops, const MixedMesh &m, int u, int v) {
  std::vector<double> area(loops.size());
  for (std::size_t i = 0; i < loops.size(); ++i)
    for (std::size_t j = 0; j < loops[i].vertices.size(); ++j) {
      const Vec a =
          vec(m.points_m[static_cast<std::size_t>(loops[i].vertices[j])]);
      const Vec b = vec(m.points_m[static_cast<std::size_t>(
          loops[i].vertices[(j + 1) % loops[i].vertices.size()])]);
      area[i] += a[u] * b[v] - b[u] * a[v];
    }
  for (std::size_t i = 0; i < loops.size(); ++i) {
    double best = std::numeric_limits<double>::max();
    const Vec sample =
        vec(m.points_m[static_cast<std::size_t>(loops[i].vertices[0])]);
    for (std::size_t j = 0; j < loops.size(); ++j)
      if (i != j && std::abs(area[j]) > std::abs(area[i]) &&
          std::abs(area[j]) < best && inPolygon(sample, loops[j], m, u, v)) {
        best = std::abs(area[j]);
        loops[i].parent = static_cast<int>(j);
      }
  }
  for (auto &loop : loops) {
    for (int parent = loop.parent; parent >= 0;
         parent = loops[static_cast<std::size_t>(parent)].parent)
      ++loop.depth;
  }
}

SurfaceMesh cavitySurface(const MixedMesh &mesh,
                          const std::vector<SkinFace> &skin) {
  SurfaceMesh result;
  result.vertices = mesh.points_m;
  for (const auto &f : skin) {
    if (mesh.generation_dimension == 3)
      result.triangles.push_back({{static_cast<std::size_t>(f.vertices[0]),
                                   static_cast<std::size_t>(f.vertices[1]),
                                   static_cast<std::size_t>(f.vertices[2])},
                                  0});
    else {
      const auto a = result.vertices.at(
                     static_cast<std::size_t>(f.vertices[0])),
                 b = result.vertices.at(
                     static_cast<std::size_t>(f.vertices[1]));
      auto A = a, B = b;
      A[static_cast<std::size_t>(mesh.extrusion_axis)] -= 1;
      B[static_cast<std::size_t>(mesh.extrusion_axis)] += 1;
      auto C = a, D = b;
      C[static_cast<std::size_t>(mesh.extrusion_axis)] += 1;
      D[static_cast<std::size_t>(mesh.extrusion_axis)] -= 1;
      const auto start = result.vertices.size();
      result.vertices.insert(result.vertices.end(), {A, D, B, C});
      result.triangles.push_back({{start, start + 1, start + 2}, 0});
      result.triangles.push_back({{start, start + 2, start + 3}, 0});
    }
  }
  buildSurfaceIndex(result);
  return result;
}
int backgroundSide(const MixedMesh &mesh, const FaceKey &ring, const Json &p) {
  const auto &bounds = p.at("background").at("bounds");
  const double scale = scaleOf(p);
  for (int axis = 0; axis < 3; ++axis)
    for (int side = 0; side < 2; ++side) {
      const double value =
          bounds.at(side ? "maximum" : "minimum")[axis].get<double>() * scale;
      const double span = (bounds.at("maximum")[axis].get<double>() -
                           bounds.at("minimum")[axis].get<double>()) *
                          scale;
      if (std::all_of(ring.begin(), ring.end(), [&](auto node) {
            return std::abs(mesh.points_m.at(static_cast<std::size_t>(
                                node))[static_cast<std::size_t>(axis)] -
                            value) <= span * 1e-9;
          }))
        return 2 * axis + side;
    }
  return -1;
}
double faceArea(const MixedMesh &mesh, const FaceKey &ring) {
  const Vec centre = faceCentre(mesh, ring);
  Vec area = Vec::Zero();
  for (std::size_t j = 0; j < ring.size(); ++j)
    area += (vec(mesh.points_m.at(static_cast<std::size_t>(ring[j]))) - centre)
                .cross(vec(mesh.points_m.at(static_cast<std::size_t>(
                           ring[(j + 1) % ring.size()]))) -
                       centre);
  return area.norm() * .5;
}
std::vector<SkinFace>
coreCells(MixedMesh &mesh, std::vector<PlanarCell> &planar, const CadModel &cad,
          const Json &p, const std::vector<SkinFace> &skin, DomainQuery &query,
          const BackgroundGrid &grid, bool retainOuter = false) {
  (void)cad;
  const auto surface = cavitySurface(mesh, skin);
  std::vector<SkinFace> fitSkin;
  for (const auto &face : skin)
    if (!retainOuter || backgroundSide(mesh, face.vertices, p) < 0)
      fitSkin.push_back(face);
  const auto collisionSurface =
      retainOuter ? cavitySurface(mesh, fitSkin) : surface;
  const auto loops = mesh.generation_dimension == 2 ? loopsFromEdges(skin)
                                                    : std::vector<Loop>{};
  const int axis = mesh.extrusion_axis, u = (axis + 1) % 3, v = (axis + 2) % 3;
  std::vector<FaceKey> cells;
  std::vector<MixedCell> volumes;
  std::size_t nonFluid = 0, reserved = 0, refinementTransition = 0,
              morphologyRemoved = 0;
  double minimumClearance = std::numeric_limits<double>::max();
  const auto size = mesh.generation_dimension == 3 ? grid.volumes.size()
                                                   : grid.polygons.size();
  for (std::size_t i = 0; i < size; ++i) {
    auto vertices = mesh.generation_dimension == 3 ? grid.volumes[i].vertices
                                                   : grid.polygons[i];
    Vec low = Vec::Constant(std::numeric_limits<double>::max()), high = -low;
    for (auto node : vertices) {
      low = low.cwiseMin(vec(mesh.points_m.at(static_cast<std::size_t>(node))));
      high =
          high.cwiseMax(vec(mesh.points_m.at(static_cast<std::size_t>(node))));
    }
    const Vec centre = (low + high) * .5;
    if (!query.inside(centre)) {
      ++nonFluid;
      continue;
    }
    // A hanging-face Cartesian cell has poor coarse/fine centroid alignment.
    // Reserve it for a conforming hybrid transition instead of keeping a
    // skewed mortar interface in the final 3D Cartesian core.
    if (mesh.generation_dimension == 3 &&
        grid.volumes[i].type == CellType::polyhedron) {
      ++refinementTransition;
      continue;
    }
    bool inside = false;
    if (mesh.generation_dimension == 3)
      inside = pointInsideSurface(surface, arr(centre));
    else
      for (const auto &loop : loops)
        if (inPolygon(centre, loop, mesh, u, v))
          inside = !inside;
    if (!inside) {
      ++reserved;
      continue;
    }
    double step = 0;
    for (int a = 0; a < 3; ++a)
      if (mesh.generation_dimension == 3 || a != axis)
        step = std::max(step, high[a] - low[a]);
    const double clearance =
        step * std::max(.35, p.at("meshing").value("core_buffer", .35));
    std::array<double, 6> box;
    for (int a = 0; a < 3; ++a) {
      box[static_cast<std::size_t>(a)] = low[a] - clearance;
      box[static_cast<std::size_t>(a + 3)] = high[a] + clearance;
    }
    if (mesh.generation_dimension == 2) {
      box[static_cast<std::size_t>(axis)] = centre[axis] - 1e-10;
      box[static_cast<std::size_t>(axis + 3)] = centre[axis] + 1e-10;
    }
    if (!collisionSurface.triangles.empty() &&
        surfaceIntersectsBox(collisionSurface, box, 1e-10)) {
      ++reserved;
      continue;
    }
    minimumClearance = std::min(minimumClearance, clearance);
    cells.push_back(vertices);
    if (mesh.generation_dimension == 3) {
      auto cell = grid.volumes[i];
      cell.vertices = vertices;
      volumes.push_back(std::move(cell));
    }
  }
  mesh.metadata["core_clearance_m"] =
      minimumClearance == std::numeric_limits<double>::max() ? 0
                                                             : minimumClearance;
  // Remove edge-touching core configurations instead of passing non-manifold
  // PLCs to the filler.
  std::vector<bool> keep(cells.size(), true);
  for (;;) {
    std::map<FaceKey, std::vector<std::pair<std::size_t, FaceKey>>> faces;
    for (std::size_t i = 0; i < cells.size(); ++i)
      if (keep[i]) {
        std::vector<FaceKey> rings;
        if (mesh.generation_dimension == 3)
          rings = cellFaceRings(volumes[i]);
        else
          for (std::size_t j = 0; j < cells[i].size(); ++j)
            rings.push_back({cells[i][j], cells[i][(j + 1) % cells[i].size()]});
        for (auto &ring : rings)
          faces[faceKey(ring)].push_back({i, std::move(ring)});
      }
    std::map<FaceKey, std::vector<std::size_t>> incidence;
    for (const auto &f : faces)
      if (f.second.size() == 1) {
        const auto &ring = f.second[0].second;
        if (mesh.generation_dimension == 3)
          for (std::size_t j = 0; j < ring.size(); ++j)
            incidence[faceKey({ring[j], ring[(j + 1) % ring.size()]})]
                .push_back(f.second[0].first);
        else
          for (auto nodeId : ring)
            incidence[{nodeId}].push_back(f.second[0].first);
      }
    std::set<std::size_t> remove;
    for (const auto &item : incidence)
      if (item.second.size() != 2)
        for (auto c : item.second)
          remove.insert(c);
    if (remove.empty()) {
      std::vector<SkinFace> shell;
      for (const auto &f : faces)
        if (f.second.size() == 1)
          shell.push_back({f.second[0].second, -1});
      for (std::size_t i = 0; i < cells.size(); ++i)
        if (keep[i]) {
          if (mesh.generation_dimension == 3)
            mesh.cells.push_back(std::move(volumes[i]));
          else
            planar.push_back({std::move(cells[i]), CellRegion::cartesian});
        }
      mesh.metadata["background_carving"] = {
          {"candidate_cells", size},
          {"removed_nonfluid_cells", nonFluid},
          {"reserved_fit_and_layer_cells", reserved},
          {"reserved_refinement_transition_cells", refinementTransition},
          {"removed_nonmanifold_core_cells", morphologyRemoved},
          {"retained_cartesian_cells", cells.size() - morphologyRemoved},
          {"effective_core_buffer",
           std::max(.35, p.at("meshing").value("core_buffer", .35))}};
      return shell;
    }
    for (auto i : remove) {
      keep[i] = false;
      ++morphologyRemoved;
    }
  }
}

std::map<Index, Json> layerTargets(const MixedMesh &m, const Json &p) {
  std::map<Index, Json> result;
  if (!p.at("boundary_layers").at("enabled").get<bool>())
    return result;
  for (const auto &target : p.at("boundary_layers").at("targets")) {
    const auto name = target.at("group").get<std::string>();
    Index patch = -1;
    for (std::size_t i = 0; i < m.patches.size(); ++i)
      if (m.patches[i].name == name)
        patch = static_cast<Index>(i);
    if (patch < 0)
      throw std::runtime_error(
          "boundary-layer target '" + name +
          "' has no surface in the selected generation domain");
    if (m.patches[static_cast<std::size_t>(patch)].kind != "wall")
      throw std::runtime_error("boundary layers require a wall group");
    if (!result.emplace(patch, target).second)
      throw std::runtime_error("duplicate boundary-layer target " + name);
  }
  return result;
}
double totalHeight(const Json &target, double scale) {
  double h = positive(target, "first_layer") * scale, total = 0;
  const double g = positive(target, "growth_ratio");
  const int n = target.at("count").get<int>();
  for (int i = 0; i < n; ++i) {
    total += h;
    h *= g;
  }
  return total;
}
Json measuredLayerHeights(const std::vector<std::vector<double>> &samples,
                          double first, double growth) {
  Json result = Json::array();
  double expected = first;
  for (std::size_t layer = 0; layer < samples.size(); ++layer) {
    auto values = samples[layer];
    if (values.empty() || *std::min_element(values.begin(), values.end()) <= 0)
      throw std::runtime_error(
          "wall layer has no positive measured normal height");
    std::vector<double> error;
    for (double value : values)
      error.push_back(std::abs(value / expected - 1));
    std::sort(error.begin(), error.end());
    result.push_back(
        {{"layer", layer + 1},
         {"requested_height_m", expected},
         {"measured_min_normal_height_m",
          *std::min_element(values.begin(), values.end())},
         {"measured_max_normal_height_m",
          *std::max_element(values.begin(), values.end())},
         {"p95_relative_height_error",
          error[static_cast<std::size_t>(std::ceil(.95 * error.size())) - 1]},
         {"samples", values.size()}});
    expected *= growth;
  }
  return result;
}

void layers2(MixedMesh &m, std::vector<SkinFace> &skin,
             std::vector<PlanarCell> &cells, const Json &p,
             DomainQuery &query) {
  const auto targets = layerTargets(m, p);
  if (targets.empty())
    return;
  auto loops = loopsFromEdges(skin);
  const int axis = m.extrusion_axis, u = (axis + 1) % 3, v = (axis + 2) % 3;
  const double scale = scaleOf(p);
  std::vector<SkinFace> replacement;
  Json reports = Json::array();
  for (const auto &loop : loops) {
    const auto target = targets.find(loop.patches[0]);
    if (target == targets.end()) {
      if (std::any_of(loop.patches.begin(), loop.patches.end(),
                      [&](auto patch) { return targets.count(patch) != 0; }))
        throw std::runtime_error(
            "2D boundary-layer target must cover a closed boundary loop");
      for (std::size_t j = 0; j < loop.vertices.size(); ++j)
        replacement.push_back(
            {{loop.vertices[j], loop.vertices[(j + 1) % loop.vertices.size()]},
             loop.patches[j]});
      continue;
    }
    if (std::any_of(loop.patches.begin(), loop.patches.end(),
                    [&](auto patch) { return patch != target->first; }))
      throw std::runtime_error(
          "2D boundary-layer target must cover the whole closed loop");
    const std::size_t n = loop.vertices.size();
    std::vector<Vec> normals(n);
    for (std::size_t j = 0; j < n; ++j) {
      const Vec a = vec(m.points_m[static_cast<std::size_t>(loop.vertices[j])]),
                b = vec(m.points_m[static_cast<std::size_t>(
                    loop.vertices[(j + 1) % n])]);
      Vec direction = Vec::Zero();
      direction[u] = -(b[v] - a[v]);
      direction[v] = b[u] - a[u];
      direction.normalize();
      if (!query.inside((a + b) * .5 +
                        direction * positive(p.at("sizing"), "surface_size") *
                            scale * .01))
        direction = -direction;
      normals[j] = direction;
    }
    std::vector<Vec> displacements(n);
    for (std::size_t j = 0; j < n; ++j) {
      Vec d = normals[(j + n - 1) % n] + normals[j];
      const double dot = d.dot(normals[j]);
      if (dot < .1)
        throw std::runtime_error("2D boundary-layer corner is too sharp");
      displacements[j] = d / dot;
    }
    FaceKey lower = loop.vertices;
    double distance = 0, h = positive(target->second, "first_layer") * scale;
    const int count = target->second.at("count").get<int>();
    std::vector<std::vector<double>> heights(static_cast<std::size_t>(count));
    for (int layer = 0; layer < count; ++layer) {
      distance += h;
      FaceKey upper;
      for (std::size_t j = 0; j < n; ++j) {
        const Vec q =
            vec(m.points_m[static_cast<std::size_t>(loop.vertices[j])]) +
            displacements[j] * distance;
        // Along-segment samples also reject crossing an intervening solid or
        // thin gap.
        for (int sample = 1; sample <= 8; ++sample)
          if (!query.inside(
                  vec(m.points_m[static_cast<std::size_t>(loop.vertices[j])]) +
                  displacements[j] * (distance * sample / 8)))
            throw std::runtime_error(
                "2D boundary layer exceeds available fluid clearance");
        upper.push_back(addPoint(m, q));
      }
      for (std::size_t j = 0; j < n; ++j) {
        const auto k = (j + 1) % n;
        for (auto node : {j, k})
          heights[static_cast<std::size_t>(layer)].push_back(
              (vec(m.points_m[static_cast<std::size_t>(upper[node])]) -
               vec(m.points_m[static_cast<std::size_t>(lower[node])]))
                  .dot(normals[j]));
        FaceKey a = {lower[j], lower[k], upper[k]},
                b = {lower[j], upper[k], upper[j]};
        auto orient = [&](FaceKey &tri) {
          const Vec x = vec(m.points_m[static_cast<std::size_t>(tri[0])]),
                    y = vec(m.points_m[static_cast<std::size_t>(tri[1])]),
                    z = vec(m.points_m[static_cast<std::size_t>(tri[2])]);
          const double signedArea =
              (y[u] - x[u]) * (z[v] - x[v]) - (y[v] - x[v]) * (z[u] - x[u]);
          if (std::abs(signedArea) < h * h * 1e-10)
            throw std::runtime_error("2D layer strip is degenerate");
          if (signedArea < 0)
            std::swap(tri[1], tri[2]);
        };
        orient(a);
        orient(b);
        cells.push_back({std::move(a), CellRegion::boundary_layer});
        cells.push_back({std::move(b), CellRegion::boundary_layer});
      }
      lower = std::move(upper);
      h *= target->second.at("growth_ratio").get<double>();
    }
    for (std::size_t j = 0; j < n; ++j)
      replacement.push_back({{lower[j], lower[(j + 1) % n]}, target->first});
    reports.push_back(
        {{"group", m.patches[static_cast<std::size_t>(target->first)].name},
         {"requested_count", count},
         {"generated_count", count},
         {"first_layer_m",
          target->second.at("first_layer").get<double>() * scale},
         {"growth_ratio", target->second.at("growth_ratio")},
         {"total_thickness_m", distance},
         {"height_measurements",
          measuredLayerHeights(heights,
                               positive(target->second, "first_layer") * scale,
                               positive(target->second, "growth_ratio"))},
         {"coverage", 1.0}});
  }
  skin = std::move(replacement);
  m.metadata["boundary_layers"] = reports;
}

void layers3(MixedMesh &m, std::vector<SkinFace> &skin, const Json &p,
             DomainQuery &query) {
  const auto targets = layerTargets(m, p);
  if (targets.empty())
    return;
  const double scale = scaleOf(p);
  std::vector<SkinFace> replacement;
  for (const auto &f : skin)
    if (!targets.count(f.patch))
      replacement.push_back(f);
  Json reports = Json::array();
  std::set<Index> usedNodes;
  for (const auto &target : targets) {
    std::vector<SkinFace> wall;
    std::set<Index> nodes;
    std::map<FaceKey, int> edges;
    for (const auto &f : skin)
      if (f.patch == target.first) {
        wall.push_back(f);
        for (auto n : f.vertices)
          nodes.insert(n);
        for (std::size_t j = 0; j < 3; ++j)
          ++edges[faceKey({f.vertices[j], f.vertices[(j + 1) % 3]})];
      }
    if (std::any_of(edges.begin(), edges.end(),
                    [](const auto &e) { return e.second != 2; }))
      throw std::runtime_error(
          "3D layer target must be a closed wall shell; partial-patch layer "
          "termination needs an explicit geometric partition");
    for (auto n : nodes)
      if (!usedNodes.insert(n).second)
        throw std::runtime_error("boundary-layer targets touch with "
                                 "incompatible layer specifications");
    std::map<Index, Eigen::Matrix3d> lhs;
    std::map<Index, Vec> rhs;
    for (auto n : nodes) {
      lhs[n] = Eigen::Matrix3d::Zero();
      rhs[n] = Vec::Zero();
    }
    for (const auto &f : wall) {
      const Vec a = vec(m.points_m[static_cast<std::size_t>(f.vertices[0])]);
      const Vec normal =
          (vec(m.points_m[static_cast<std::size_t>(f.vertices[1])]) - a)
              .cross(vec(m.points_m[static_cast<std::size_t>(f.vertices[2])]) -
                     a)
              .normalized();
      for (auto n : f.vertices) {
        lhs[n] += normal * normal.transpose();
        rhs[n] -= normal;
      }
    }
    std::map<Index, Vec> displacement;
    for (auto n : nodes) {
      const Eigen::Matrix3d matrix =
          lhs.at(n) + Eigen::Matrix3d::Identity() * 1e-6;
      const Vec d = matrix.ldlt().solve(rhs.at(n));
      if (!d.allFinite() || d.norm() > 4 || d.norm() < .5)
        throw std::runtime_error(
            "3D layer normal is ill-conditioned at a sharp corner");
      displacement[n] = d;
    }
    std::map<Index, Index> lower;
    for (auto n : nodes)
      lower[n] = n;
    double distance = 0, h = positive(target.second, "first_layer") * scale;
    const int count = target.second.at("count").get<int>();
    std::vector<std::vector<double>> heights(static_cast<std::size_t>(count));
    for (int layer = 0; layer < count; ++layer) {
      distance += h;
      std::map<Index, Index> upper;
      for (auto n : nodes) {
        const Vec base = vec(m.points_m[static_cast<std::size_t>(n)]);
        for (int sample = 1; sample <= 8; ++sample)
          if (!query.inside(base +
                            displacement.at(n) * (distance * sample / 8)))
            throw std::runtime_error(
                "3D boundary layer exceeds available fluid clearance");
        upper[n] = addPoint(m, base + displacement.at(n) * distance);
      }
      for (const auto &f : wall) {
        const Vec base =
            vec(m.points_m[static_cast<std::size_t>(f.vertices[0])]);
        const Vec inward =
            -(vec(m.points_m[static_cast<std::size_t>(f.vertices[1])]) - base)
                 .cross(
                     vec(m.points_m[static_cast<std::size_t>(f.vertices[2])]) -
                     base)
                 .normalized();
        for (auto node : f.vertices)
          heights[static_cast<std::size_t>(layer)].push_back(
              (vec(m.points_m[static_cast<std::size_t>(upper.at(node))]) -
               vec(m.points_m[static_cast<std::size_t>(lower.at(node))]))
                  .dot(inward));
        auto ids = f.vertices;
        std::sort(ids.begin(), ids.end());
        const Index a = lower.at(ids[0]), b = lower.at(ids[1]),
                    c = lower.at(ids[2]);
        const Index A = upper.at(ids[0]), B = upper.at(ids[1]),
                    C = upper.at(ids[2]);
        addTet(m, {a, b, c, C}, CellRegion::boundary_layer);
        addTet(m, {a, b, B, C}, CellRegion::boundary_layer);
        addTet(m, {a, A, B, C}, CellRegion::boundary_layer);
      }
      lower = std::move(upper);
      h *= target.second.at("growth_ratio").get<double>();
    }
    for (const auto &f : wall)
      replacement.push_back({{lower.at(f.vertices[0]), lower.at(f.vertices[1]),
                              lower.at(f.vertices[2])},
                             f.patch});
    reports.push_back(
        {{"group", m.patches[static_cast<std::size_t>(target.first)].name},
         {"requested_count", count},
         {"generated_count", count},
         {"first_layer_m",
          target.second.at("first_layer").get<double>() * scale},
         {"growth_ratio", target.second.at("growth_ratio")},
         {"total_thickness_m", distance},
         {"height_measurements",
          measuredLayerHeights(heights,
                               positive(target.second, "first_layer") * scale,
                               positive(target.second, "growth_ratio"))},
         {"coverage", 1.0}});
  }
  skin = std::move(replacement);
  m.metadata["boundary_layers"] = reports;
}

void injectAllNodes(MixedMesh &m, int dim, int tag) {
  std::vector<std::size_t> ids(m.points_m.size());
  std::iota(ids.begin(), ids.end(), std::size_t{1});
  std::vector<double> xyz;
  xyz.reserve(3 * m.points_m.size());
  for (const auto &p : m.points_m)
    for (auto v : p)
      xyz.push_back(v);
  gmsh::model::mesh::addNodes(dim, tag, ids, xyz);
}
void refreshNodes(MixedMesh &m) {
  std::vector<std::size_t> tags;
  std::vector<double> xyz, param;
  gmsh::model::mesh::getNodes(tags, xyz, param);
  const auto size =
      tags.empty() ? 0 : *std::max_element(tags.begin(), tags.end());
  m.points_m.resize(std::max(size, m.points_m.size()),
                    {std::numeric_limits<double>::quiet_NaN(), 0, 0});
  for (std::size_t i = 0; i < tags.size(); ++i)
    m.points_m.at(tags[i] - 1) = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
}
void fill3(MixedMesh &m, const std::vector<SkinFace> &skin,
           const std::vector<SkinFace> &core, const Json &p, double h) {
  if (skin.empty() && core.empty())
    return;
  // Place the transition pyramids explicitly. A bounded apex height avoids the
  // arbitrarily flat pyramids produced by unconstrained automatic recovery.
  std::vector<SkinFace> pyramidShell;
  const double clearance = m.metadata.at("core_clearance_m").get<double>();
  for (const auto &face : core) {
    const Vec centre = faceCentre(m, face.vertices);
    Vec normal = Vec::Zero();
    double minimumEdge = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < face.vertices.size(); ++j) {
      const Vec q = vec(
                    m.points_m.at(static_cast<std::size_t>(face.vertices[j]))),
                r = vec(m.points_m.at(static_cast<std::size_t>(
                    face.vertices[(j + 1) % face.vertices.size()])));
      normal += (q - centre).cross(r - centre);
      minimumEdge = std::min(minimumEdge, (r - q).norm());
    }
    normal.normalize();
    const double height = std::min(.4 * minimumEdge, clearance * .8);
    auto ids = face.vertices;
    const auto apex = addPoint(m, centre + normal * height);
    ids.push_back(apex);
    auto lower = face.vertices;
    std::reverse(lower.begin(), lower.end());
    std::vector<FaceKey> rings = {lower};
    for (std::size_t j = 0; j < face.vertices.size(); ++j)
      rings.push_back({face.vertices[j],
                       face.vertices[(j + 1) % face.vertices.size()], apex});
    if (face.vertices.size() == 3) {
      addTet(m, ids, CellRegion::transition);
    } else if (face.vertices.size() == 4)
      m.cells.emplace_back(CellType::pyramid, CellRegion::transition, ids);
    else
      m.cells.emplace_back(CellType::polyhedron, CellRegion::transition, ids,
                           rings);
    for (std::size_t j = 1; j < rings.size(); ++j)
      pyramidShell.push_back({rings[j], -1});
  }
  gmsh::clear();
  gmsh::model::add("hybrid-cavity");
  options(p, h, "transition");
  const int outer = gmsh::model::addDiscreteEntity(2);
  // Background nodes inside retained cells and unused CAD box nodes must not
  // become points in the cavity Delaunay triangulation.
  std::set<Index> referenced;
  const std::array<const std::vector<SkinFace> *, 2> fillShells = {
      &skin, &pyramidShell};
  for (const auto *shell : fillShells)
    for (const auto &face : *shell)
      referenced.insert(face.vertices.begin(), face.vertices.end());
  std::vector<std::size_t> nodeTags;
  std::vector<double> xyz;
  for (auto node : referenced) {
    nodeTags.push_back(static_cast<std::size_t>(node) + 1);
    for (double coordinate : m.points_m.at(static_cast<std::size_t>(node)))
      xyz.push_back(coordinate);
  }
  gmsh::model::mesh::addNodes(2, outer, nodeTags, xyz);
  std::vector<SkinFace> plc = skin;
  for (auto face : pyramidShell) {
    // Core/pyramid outward faces point into the remaining cavity.
    std::reverse(face.vertices.begin(), face.vertices.end());
    plc.push_back(std::move(face));
  }
  // Build one entity per connected closed shell, then explicitly assign holes
  // to the smallest enclosing outward shell. An opaque disconnected surface
  // cannot tell the volume backend which component is a solid exclusion.
  std::map<FaceKey, std::vector<std::size_t>> edgeFaces;
  for (std::size_t i = 0; i < plc.size(); ++i)
    for (std::size_t j = 0; j < 3; ++j)
      edgeFaces[faceKey({plc[i].vertices[j], plc[i].vertices[(j + 1) % 3]})]
          .push_back(i);
  std::vector<std::vector<std::size_t>> neighbours(plc.size());
  for (const auto &item : edgeFaces) {
    if (item.second.size() != 2)
      throw std::runtime_error(
          "3D carved cavity has an open or non-manifold PLC edge");
    neighbours[item.second[0]].push_back(item.second[1]);
    neighbours[item.second[1]].push_back(item.second[0]);
  }
  std::vector<bool> visited(plc.size());
  std::vector<std::vector<SkinFace>> components;
  std::vector<double> signedVolumes;
  std::vector<int> shellTags;
  for (std::size_t start = 0; start < plc.size(); ++start) {
    if (visited[start])
      continue;
    std::vector<SkinFace> component;
    std::queue<std::size_t> pending;
    pending.push(start);
    visited[start] = true;
    while (!pending.empty()) {
      const auto face = pending.front();
      pending.pop();
      component.push_back(plc[face]);
      for (auto next : neighbours[face])
        if (!visited[next]) {
          visited[next] = true;
          pending.push(next);
        }
    }
    const Vec origin = faceCentre(m, component.front().vertices);
    double volume = 0;
    std::vector<std::size_t> triangles;
    for (const auto &face : component) {
      const auto &ring = face.vertices;
      volume +=
          (vec(m.points_m.at(static_cast<std::size_t>(ring[0]))) - origin)
              .dot((vec(m.points_m.at(static_cast<std::size_t>(ring[1]))) -
                    origin)
                       .cross(vec(m.points_m.at(
                                  static_cast<std::size_t>(ring[2]))) -
                              origin)) /
          6;
      for (auto node : ring)
        triangles.push_back(static_cast<std::size_t>(node) + 1);
    }
    if (!std::isfinite(volume) || volume == 0)
      throw std::runtime_error("3D cavity shell has zero signed volume");
    const int tag =
        components.empty() ? outer : gmsh::model::addDiscreteEntity(2);
    gmsh::model::mesh::addElementsByType(tag, 2, {}, triangles);
    shellTags.push_back(gmsh::model::geo::addSurfaceLoop({tag}));
    signedVolumes.push_back(volume);
    components.push_back(std::move(component));
  }
  std::vector<SurfaceMesh> indexes(components.size());
  std::vector<std::vector<int>> holes(components.size());
  for (std::size_t i = 0; i < components.size(); ++i)
    if (signedVolumes[i] > 0)
      indexes[i] = cavitySurface(m, components[i]);
  for (std::size_t i = 0; i < components.size(); ++i) {
    if (signedVolumes[i] > 0)
      continue;
    const auto point = m.points_m.at(
        static_cast<std::size_t>(components[i].front().vertices[0]));
    std::size_t parent = components.size();
    double smallest = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < components.size(); ++j)
      if (signedVolumes[j] > 0 && signedVolumes[j] < smallest &&
          pointInsideSurface(indexes[j], point)) {
        parent = j;
        smallest = signedVolumes[j];
      }
    if (parent == components.size())
      throw std::runtime_error("3D cavity hole has no enclosing fluid shell");
    holes[parent].push_back(shellTags[i]);
  }
  std::size_t cavities = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
    if (signedVolumes[i] > 0) {
      auto shells = holes[i];
      shells.insert(shells.begin(), shellTags[i]);
      gmsh::model::geo::addVolume(shells);
      ++cavities;
    }
  if (!cavities)
    throw std::runtime_error("3D cavity has no outward fluid shell");
  m.metadata["cavity_components"] = cavities;
  m.metadata["cavity_exclusion_shells"] = components.size() - cavities;
  gmsh::model::geo::synchronize();
  gmsh::option::setNumber("Mesh.MeshOnlyEmpty", 1);
  gmsh::model::mesh::generate(3);
  refreshNodes(m);
  std::vector<int> types;
  std::vector<std::vector<std::size_t>> ids, connectivity;
  gmsh::model::mesh::getElements(types, ids, connectivity, 3);
  for (std::size_t i = 0; i < types.size(); ++i) {
    const int size = types[i] == 4 ? 4 : (types[i] == 7 ? 5 : 0);
    if (!size)
      throw std::runtime_error(
          "3D cavity mesher returned an unsupported cell type");
    for (std::size_t j = 0; j < connectivity[i].size();
         j += static_cast<std::size_t>(size)) {
      FaceKey ring;
      for (int k = 0; k < size; ++k)
        ring.push_back(static_cast<Index>(
            connectivity[i][j + static_cast<std::size_t>(k)] - 1));
      m.cells.push_back(
          {types[i] == 4 ? CellType::tetrahedron : CellType::pyramid,
           CellRegion::transition, std::move(ring)});
    }
  }
}
void fill2(MixedMesh &m, std::vector<PlanarCell> &cells,
           const std::vector<SkinFace> &skin, const std::vector<SkinFace> &core,
           const Json &p, double h) {
  if (skin.empty() && core.empty())
    return;
  std::vector<SkinFace> edges = skin;
  edges.insert(edges.end(), core.begin(), core.end());
  auto loops = loopsFromEdges(edges);
  nestLoops(loops, m, (m.extrusion_axis + 1) % 3, (m.extrusion_axis + 2) % 3);
  const int axis = m.extrusion_axis, u = (axis + 1) % 3, v = (axis + 2) % 3;
  const double coordinate = m.points_m.at(static_cast<std::size_t>(
      skin.front().vertices.front()))[static_cast<std::size_t>(axis)];
  // Discrete curves carry no plane parametrization. Mesh in an explicit local
  // XY plane, then restore the original section axis and exact coordinate.
  for (auto &node : m.points_m)
    node = {node[static_cast<std::size_t>(u)],
            node[static_cast<std::size_t>(v)], 0};
  gmsh::clear();
  gmsh::model::add("planar-cavity");
  options(p, h, "transition");
  std::vector<int> curveLoops;
  for (std::size_t i = 0; i < loops.size(); ++i) {
    const int curve = gmsh::model::addDiscreteEntity(1);
    if (i == 0)
      injectAllNodes(m, 1, curve);
    std::vector<std::size_t> lines;
    const auto &ring = loops[i].vertices;
    for (std::size_t j = 0; j < ring.size(); ++j) {
      lines.push_back(static_cast<std::size_t>(ring[j]) + 1);
      lines.push_back(static_cast<std::size_t>(ring[(j + 1) % ring.size()]) +
                      1);
    }
    gmsh::model::mesh::addElementsByType(curve, 1, {}, lines);
    curveLoops.push_back(gmsh::model::geo::addCurveLoop({curve}));
  }
  for (std::size_t i = 0; i < loops.size(); ++i)
    if (loops[i].depth % 2 == 0) {
      std::vector<int> boundary = {curveLoops[i]};
      for (std::size_t j = 0; j < loops.size(); ++j)
        if (loops[j].parent == static_cast<int>(i))
          boundary.push_back(curveLoops[j]);
      gmsh::model::geo::addPlaneSurface(boundary);
    }
  gmsh::model::geo::synchronize();
  gmsh::option::setNumber("Mesh.MeshOnlyEmpty", 1);
  gmsh::model::mesh::generate(2);
  refreshNodes(m);
  for (auto &node : m.points_m) {
    const auto local = node;
    node[static_cast<std::size_t>(axis)] = coordinate;
    node[static_cast<std::size_t>(u)] = local[0];
    node[static_cast<std::size_t>(v)] = local[1];
  }
  std::vector<int> types;
  std::vector<std::vector<std::size_t>> ids, connectivity;
  gmsh::model::mesh::getElements(types, ids, connectivity, 2);
  for (std::size_t i = 0; i < types.size(); ++i) {
    if (types[i] != 2)
      throw std::runtime_error(
          "planar cavity mesher returned an unsupported cell");
    for (std::size_t j = 0; j < connectivity[i].size(); j += 3) {
      FaceKey ring;
      for (std::size_t k = 0; k < 3; ++k)
        ring.push_back(static_cast<Index>(connectivity[i][j + k] - 1));
      const Vec a = vec(m.points_m[static_cast<std::size_t>(ring[0])]),
                b = vec(m.points_m[static_cast<std::size_t>(ring[1])]),
                c = vec(m.points_m[static_cast<std::size_t>(ring[2])]);
      if ((b[u] - a[u]) * (c[v] - a[v]) - (b[v] - a[v]) * (c[u] - a[u]) < 0)
        std::swap(ring[1], ring[2]);
      cells.push_back({std::move(ring), CellRegion::transition});
    }
  }
}
void extrude2(MixedMesh &m, const std::vector<PlanarCell> &cells,
              const std::vector<SkinFace> &skin, BoundaryMap &boundary) {
  const Index size = static_cast<Index>(m.points_m.size());
  const auto plane = m.points_m;
  const Index front = addPatch(m, "front", "symmetry"),
              back = addPatch(m, "back", "symmetry");
  for (auto &p : m.points_m)
    p[static_cast<std::size_t>(m.extrusion_axis)] -=
        m.extrusion_thickness_m * .5;
  for (auto p : plane) {
    p[static_cast<std::size_t>(m.extrusion_axis)] +=
        m.extrusion_thickness_m * .5;
    m.points_m.push_back(p);
  }
  for (const auto &c : cells) {
    auto ids = c.vertices;
    for (auto n : c.vertices)
      ids.push_back(n + size);
    if (c.vertices.size() == 3 || c.vertices.size() == 4)
      m.cells.emplace_back(c.vertices.size() == 4 ? CellType::hexahedron
                                                  : CellType::prism,
                           c.region, ids);
    else {
      auto low = c.vertices;
      std::reverse(low.begin(), low.end());
      auto high = c.vertices;
      for (auto &n : high)
        n += size;
      std::vector<FaceKey> faces = {low, high};
      for (std::size_t j = 0; j < c.vertices.size(); ++j) {
        const auto a = c.vertices[j],
                   b = c.vertices[(j + 1) % c.vertices.size()];
        faces.push_back({a, b, b + size, a + size});
      }
      m.cells.emplace_back(CellType::polyhedron, c.region, ids, faces);
    }
    auto lower = c.vertices;
    std::reverse(lower.begin(), lower.end());
    boundaryFace(boundary, lower, front);
    auto upper = c.vertices;
    for (auto &n : upper)
      n += size;
    boundaryFace(boundary, upper, back);
  }
  for (const auto &f : skin) {
    const auto a = f.vertices[0], b = f.vertices[1];
    boundaryFace(boundary, {a, b, b + size, a + size}, f.patch);
  }
}
void compactNodes(MixedMesh &m, BoundaryMap &boundary) {
  std::vector<Index> map(m.points_m.size(), -1);
  std::vector<std::array<double, 3>> points;
  for (auto &c : m.cells)
    for (auto &n : c.vertices) {
      if (map.at(static_cast<std::size_t>(n)) < 0) {
        map[static_cast<std::size_t>(n)] = static_cast<Index>(points.size());
        points.push_back(m.points_m[static_cast<std::size_t>(n)]);
      }
      n = map[static_cast<std::size_t>(n)];
    }
  for (auto &cell : m.cells)
    for (auto &face : cell.poly_faces)
      for (auto &n : face)
        n = map.at(static_cast<std::size_t>(n));
  BoundaryMap remapped;
  for (const auto &f : boundary) {
    auto ring = f.first;
    for (auto &n : ring) {
      n = map.at(static_cast<std::size_t>(n));
      if (n < 0)
        throw std::runtime_error("boundary node was lost during filling");
    }
    boundaryFace(remapped, ring, f.second);
  }
  boundary = std::move(remapped);
  m.points_m = std::move(points);
}
} // namespace

Json normalizeHybridProject(Json p) {
  auto &background = p.at("background");
  if (background.contains("cell_size")) {
    if (background.contains("cell_counts"))
      throw std::runtime_error(
          "choose background.cell_counts or background.cell_size, not both");
    const auto &size = background.at("cell_size");
    if (!size.is_array() || size.size() != 3)
      throw std::runtime_error(
          "background.cell_size requires three components");
    Json counts = Json::array();
    const int dimension = p.at("meshing").at("dimension");
    const auto axis = p.at("meshing").value("section_axis", "z");
    const int extrusion = axis == "x" ? 0 : axis == "y" ? 1 : 2;
    for (int a = 0; a < 3; ++a) {
      if (!size[a].is_number() || !std::isfinite(size[a].get<double>()) ||
          size[a].get<double>() <= 0)
        throw std::runtime_error(
            "background.cell_size must be finite and positive");
      const double length =
          background.at("bounds").at("maximum")[a].get<double>() -
          background.at("bounds").at("minimum")[a].get<double>();
      const double n = dimension == 2 && a == extrusion
                           ? 1
                           : std::ceil(length / size[a].get<double>());
      if (n < 1 || n > std::numeric_limits<int>::max())
        throw std::runtime_error(
            "background.cell_size produces an invalid cell count");
      counts.push_back(static_cast<int>(n));
    }
    background["cell_counts"] = counts;
    background.erase("cell_size");
  }
  return p;
}
void validateHybridProject(const Json &p) {
  if (p.at("schema_version") != 2)
    throw std::runtime_error("hybrid meshing requires schema_version 2");
  validateOptimization(p.value("optimization", Json::object()));
  const auto &m = p.at("meshing");
  keys(m,
       {"dimension", "section_axis", "section_coordinate", "thickness",
        "core_buffer"},
       "meshing");
  if (!m.at("dimension").is_number_integer() ||
      (m.at("dimension") != 2 && m.at("dimension") != 3))
    throw std::runtime_error("meshing.dimension must be 2 or 3");
  positive(p.at("sizing"), "surface_size");
  positive(p.at("sizing"), "minimum_size");
  positive(p.at("sizing"), "maximum_size");
  if (positive(p.at("sizing"), "surface_size") <
          positive(p.at("sizing"), "minimum_size") ||
      positive(p.at("sizing"), "maximum_size") <
          positive(p.at("sizing"), "minimum_size"))
    throw std::runtime_error(
        "surface_size/maximum_size must be >= minimum_size");
  if (!p.at("sizing").at("max_level").is_number_integer() ||
      p.at("sizing").at("max_level").get<int>() < 0 ||
      p.at("sizing").at("max_level").get<int>() > 20)
    throw std::runtime_error("max_level must be in [0,20]");
  if (p.at("sizing").value("balance_face_neighbors", 1) != 1)
    throw std::runtime_error("hybrid background requires 2:1 face-neighbor "
                             "grading (balance_face_neighbors=1)");
  if (p.at("sizing").contains("volume_size"))
    positive(p.at("sizing"), "volume_size");
  std::set<std::string> sourceIds;
  for (const auto &source : p.at("sizing").at("sources")) {
    const auto id = source.at("id").get<std::string>();
    if (id.empty() || !sourceIds.insert(id).second)
      throw std::runtime_error("duplicate sizing source id");
    positive(source, "target_size");
    if (source.at("target_size").get<double>() <
        positive(p.at("sizing"), "minimum_size"))
      throw std::runtime_error(
          "sizing source target_size is below minimum_size");
    const auto apply = source.value("apply_to", "all");
    if (apply != "all" && apply != "background" && apply != "surface" &&
        apply != "transition")
      throw std::runtime_error("sizing source apply_to must be background, "
                               "surface, transition or all");
    const auto kind = source.at("kind").get<std::string>();
    if (kind != "box" && kind != "sphere" && apply != "background")
      throw std::runtime_error("hybrid surface/transition size fields support "
                               "box or sphere regions");
    if (source.contains("transition_width") &&
        (!source.at("transition_width").is_number() ||
         !std::isfinite(source.at("transition_width").get<double>()) ||
         source.at("transition_width").get<double>() < 0))
      throw std::runtime_error(
          "transition_width must be finite and nonnegative");
    if (kind == "box") {
      const auto &b = source.at("region").at("bounds");
      for (int a = 0; a < 3; ++a)
        if (!std::isfinite(b.at("minimum")[a].get<double>()) ||
            !std::isfinite(b.at("maximum")[a].get<double>()) ||
            b.at("maximum")[a].get<double>() < b.at("minimum")[a].get<double>())
          throw std::runtime_error("invalid sizing box bounds");
    }
    if (kind == "sphere") {
      positive(source.at("region"), "radius");
      if (source.at("region").at("center").size() != 3)
        throw std::runtime_error("sizing sphere needs three-component center");
      for (const auto &c : source.at("region").at("center"))
        if (!c.is_number() || !std::isfinite(c.get<double>()))
          throw std::runtime_error("sphere centre must be finite");
    }
  }
  if (!p.at("background").at("uniform_cells").get<bool>())
    throw std::runtime_error(
        "hybrid background must be Cartesian with uniform_cells=true");
  if (p.at("quality").at("profile") != "strict-cfd-v1")
    throw std::runtime_error("unsupported hybrid quality profile");
  if (p.at("quality").at("mode") != "strict" &&
      p.at("quality").at("mode") != "diagnostic")
    throw std::runtime_error("quality.mode must be strict or diagnostic");
  const auto &counts = p.at("background").at("cell_counts");
  const auto &bounds = p.at("background").at("bounds");
  if (!counts.is_array() || counts.size() != 3 ||
      bounds.at("minimum").size() != 3 || bounds.at("maximum").size() != 3)
    throw std::runtime_error("background requires three coordinate components");
  double candidate = 1;
  for (std::size_t a = 0; a < 3; ++a) {
    if (!counts[a].is_number_integer() || counts[a].get<int>() <= 0)
      throw std::runtime_error(
          "background.cell_counts must be positive integers");
    candidate *= counts[a].get<double>();
    const double lo = bounds.at("minimum")[a].get<double>(),
                 hi = bounds.at("maximum")[a].get<double>();
    if (!std::isfinite(lo) || !std::isfinite(hi) || hi <= lo)
      throw std::runtime_error("invalid background bounds");
  }
  const double maxCells = positive(p.at("resources"), "max_cells");
  positive(p.at("resources"), "max_memory_mb");
  if (candidate > maxCells)
    throw std::runtime_error(
        "Cartesian background exceeds resources.max_cells before meshing");
  if (!p.at("resources").at("threads").is_number_integer() ||
      p.at("resources").at("threads").get<int>() < 1)
    throw std::runtime_error("resources.threads must be a positive integer");
  if (m.contains("core_buffer") &&
      (!m.at("core_buffer").is_number() ||
       !std::isfinite(m.at("core_buffer").get<double>()) ||
       m.at("core_buffer").get<double>() < 0))
    throw std::runtime_error(
        "meshing.core_buffer must be finite and nonnegative");
  if (m.at("dimension") == 2) {
    const auto axis = m.at("section_axis").get<std::string>();
    if (axis != "x" && axis != "y" && axis != "z")
      throw std::runtime_error("section_axis must be x, y or z");
    const std::size_t index = axis == "x" ? 0 : (axis == "y" ? 1 : 2);
    if (counts[index] != 1)
      throw std::runtime_error("2D generation requires exactly one background "
                               "cell along the extrusion axis");
    positive(m, "thickness");
    if (!m.at("section_coordinate").is_number() ||
        !std::isfinite(m.at("section_coordinate").get<double>()))
      throw std::runtime_error("section_coordinate must be finite");
    const double coordinate = m.at("section_coordinate").get<double>();
    if (coordinate <= bounds.at("minimum")[index].get<double>() ||
        coordinate >= bounds.at("maximum")[index].get<double>())
      throw std::runtime_error(
          "section plane must lie strictly inside background bounds");
    const double half = positive(m, "thickness") * .5;
    if (coordinate - half < bounds.at("minimum")[index].get<double>() ||
        coordinate + half > bounds.at("maximum")[index].get<double>())
      throw std::runtime_error(
          "2D extrusion thickness exceeds fixed background bounds");
  } else if (m.contains("section_axis") || m.contains("section_coordinate") ||
             m.contains("thickness"))
    throw std::runtime_error(
        "3D generation does not use section/extrusion parameters");
  for (const auto &t : p.at("boundary_layers").at("targets")) {
    positive(t, "first_layer");
    positive(t, "growth_ratio");
    if (!t.at("count").is_number_integer() || t.at("count").get<int>() < 1 ||
        t.at("count").get<int>() > 256)
      throw std::runtime_error(
          "boundary-layer count must be an integer in [1,256]");
    if (t.at("growth_ratio").get<double>() < 1 ||
        !std::isfinite(totalHeight(t, scaleOf(p))))
      throw std::runtime_error(
          "boundary-layer growth ratio must be >= 1 with finite total height");
  }
  if (p.at("boundary_layers").at("enabled").get<bool>() &&
      p.at("boundary_layers").at("targets").empty())
    throw std::runtime_error(
        "enabled boundary layers require explicit targets");
  for (const auto &f : p.at("output").at("formats"))
    if (f != "babelsim-v3" && f != "vtu" && f != "fluent-ascii" &&
        f != "cgns" && f != "cgns-poly")
      throw std::runtime_error("unsupported hybrid output format");
}

MixedMesh generateHybridMesh(const Json &raw, const CadModel &cad) {
  const Json p = normalizeHybridProject(raw);
  validateHybridProject(p);
  auto resourceCheck = [&]() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0 &&
        static_cast<double>(usage.ru_maxrss) / 1024 >
            p.at("resources").at("max_memory_mb").get<double>())
      throw std::runtime_error(
          "peak resident memory exceeds resources.max_memory_mb");
  };
  resourceCheck();
  const std::lock_guard<std::mutex> guard(gmshMutex);
  Session session;
  MixedMesh m;
  m.generation_dimension = p.at("meshing").at("dimension");
  const double scale = scaleOf(p);
  double h = positive(p.at("sizing"), "surface_size") * scale;
  if (m.generation_dimension == 2) {
    const auto axis = p.at("meshing").at("section_axis").get<std::string>();
    m.extrusion_axis = axis == "x" ? 0 : (axis == "y" ? 1 : 2);
    m.extrusion_thickness_m = positive(p.at("meshing"), "thickness") * scale;
  }
  const auto &backgroundBounds = p.at("background").at("bounds");
  std::array<double, 3> low{}, high{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    low[axis] = backgroundBounds.at("minimum")[axis].get<double>() * scale;
    high[axis] = backgroundBounds.at("maximum")[axis].get<double>() * scale;
  }
  const auto &domain = p.at("fluid_domain");
  if (domain.at("mode") == "external")
    for (std::size_t axis = 0; axis < 3; ++axis) {
      low[axis] = std::max(
          low[axis],
          domain.at("outer_bounds").at("minimum")[axis].get<double>() * scale);
      high[axis] = std::min(
          high[axis],
          domain.at("outer_bounds").at("maximum")[axis].get<double>() * scale);
    }
  TopoDS_Shape envelope = cad.fluidDomain();
  if (m.generation_dimension == 2) {
    const auto axis = static_cast<std::size_t>(m.extrusion_axis);
    const double coordinate =
        p.at("meshing").at("section_coordinate").get<double>() * scale;
    low[axis] = std::max(low[axis], coordinate - m.extrusion_thickness_m * .5);
    high[axis] =
        std::min(high[axis], coordinate + m.extrusion_thickness_m * .5);
    if (coordinate - m.extrusion_thickness_m * .5 < low[axis] ||
        coordinate + m.extrusion_thickness_m * .5 > high[axis])
      throw std::runtime_error(
          "2D extrusion thickness exceeds fixed computational domain bounds");
    envelope = sectionShape(cad.fluidDomain(), m.extrusion_axis, coordinate, p);
  }
  Bnd_Box box;
  BRepBndLib::AddOptimal(envelope, box, false, false);
  box.SetGap(0);
  Standard_Real x0, y0, z0, x1, y1, z1;
  box.Get(x0, y0, z0, x1, y1, z1);
  const std::array<double, 3> cadLow{x0, y0, z0}, cadHigh{x1, y1, z1};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    if (low[axis] >= high[axis])
      throw std::runtime_error(
          "background and computational domain have no valid intersection");
    if (m.generation_dimension == 2 &&
        static_cast<int>(axis) == m.extrusion_axis)
      continue;
    const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
                            std::max({std::abs(low[axis]), std::abs(high[axis]),
                                      high[axis] - low[axis]});
    if (cadLow[axis] < low[axis] - roundoff ||
        cadHigh[axis] > high[axis] + roundoff)
      throw std::runtime_error(
          "CAD computational domain exceeds fixed background/domain bounds; "
          "adjust the declared domain instead of expanding the grid");
  }
  m.metadata["domain_bounds_constraint"] = {
      {"minimum_m", low},
      {"maximum_m", high},
      {"background_bounds", backgroundBounds},
      {"source", "declared background and computational domain"},
      {"allow_expansion", false}};
  DomainQuery query(cad.fluidDomain());
  std::cout << "[1/7] Complete Cartesian background and local refinement\n"
            << std::flush;
  const auto grid = generateBackgroundGrid(p, cad);
  m.points_m = grid.points;
  m.metadata["background"] = grid.report;
  m.metadata["generation_pipeline"] = {
      "cartesian_background", "cad_and_layer_reservation",
      "background_carving",   "layer_and_cavity_fill",
      "topology_audit",       "quality_optimization",
      "cad_fit_audit"};
  resourceCheck();
  std::cout << "[2/7] CAD boundary and wall layer reservation geometry\n"
            << std::flush;
  auto skin = cadSkin(m, cad, p, h, query);
  std::size_t adjustment = 0;
  if (m.generation_dimension == 3) {
    const double tolerance =
        p.at("quality").value("max_surface_deviation_m", h * .15);
    for (;;) {
      double maximum = 0;
      for (const auto &face : skin) {
        const Vec centre = faceCentre(m, face.vertices);
        maximum = std::max(
            maximum, projectPointToFluidBoundary(cad, arr(centre)).distance_m);
        for (std::size_t j = 0; j < 3; ++j) {
          const Vec midpoint =
              (vec(m.points_m.at(static_cast<std::size_t>(face.vertices[j]))) +
               vec(m.points_m.at(
                   static_cast<std::size_t>(face.vertices[(j + 1) % 3])))) *
              .5;
          maximum = std::max(
              maximum,
              projectPointToFluidBoundary(cad, arr(midpoint)).distance_m);
        }
      }
      if (maximum <= tolerance * .8)
        break;
      if (adjustment >= 4 ||
          h * .7 < positive(p.at("sizing"), "minimum_size") * scale)
        throw std::runtime_error(
            "surface CAD fit needs finer discretization than minimum_size or "
            "the four-adjustment limit permits");
      h *= .7;
      ++adjustment;
      std::cout << "CAD fit adjustment " << adjustment << ": surface size " << h
                << " m (transition size stays fixed)\n"
                << std::flush;
      gmsh::clear();
      m.points_m = grid.points;
      m.patches.clear();
      skin = cadSkin(m, cad, p, h, query);
      resourceCheck();
    }
  }
  m.metadata["requested_surface_size_m"] =
      positive(p.at("sizing"), "surface_size") * scale;
  m.metadata["surface_fit_adjustments"] = adjustment;
  auto physicalSkin = skin;
  resourceCheck();
  std::vector<PlanarCell> planar;
  // The offset shell reserves the requested wall thickness before carving.
  // Commit its cells only after the complete background has been carved.
  if (m.generation_dimension == 2)
    layers2(m, skin, planar, p, query);
  else
    layers3(m, skin, p, query);
  auto wallCells = std::move(m.cells);
  m.cells.clear();
  auto wallPlanar = std::move(planar);
  planar.clear();
  std::cout
      << "[3/7] Carve solids and reserved wall/fit regions from background\n"
      << std::flush;
  bool retainOuter = true;
  std::array<Index, 6> outerPatches;
  outerPatches.fill(-1);
  if (retainOuter) {
    // Keep a box face only when its physical assignment is uniform. A CAD
    // group partition or layer on that face needs the full constrained skin.
    for (const auto &face : physicalSkin) {
      const int side = backgroundSide(m, face.vertices, p);
      if (side >= 0) {
        auto &patch = outerPatches[static_cast<std::size_t>(side)];
        if (patch >= 0 && patch != face.patch)
          retainOuter = false;
        patch = face.patch;
      }
    }
    const auto targetPatches = layerTargets(m, p);
    for (std::size_t side = 0; side < outerPatches.size(); ++side) {
      if (m.generation_dimension == 2 &&
          static_cast<int>(side / 2) == m.extrusion_axis)
        continue;
      const auto patch = outerPatches[side];
      if (patch < 0 || targetPatches.count(patch))
        retainOuter = false;
    }
  }
  auto core = coreCells(m, planar, cad, p, skin, query, grid, retainOuter);
  if (retainOuter) {
    std::array<double, 6> coverage{};
    for (const auto &face : core) {
      const int side = backgroundSide(m, face.vertices, p);
      if (side >= 0)
        coverage[static_cast<std::size_t>(side)] +=
            m.generation_dimension == 3
                ? faceArea(m, face.vertices)
                : (vec(m.points_m.at(
                       static_cast<std::size_t>(face.vertices[1]))) -
                   vec(m.points_m.at(
                       static_cast<std::size_t>(face.vertices[0]))))
                      .norm();
    }
    const auto &bounds = p.at("background").at("bounds");
    Vec span;
    for (int a = 0; a < 3; ++a)
      span[a] = (bounds.at("maximum")[a].get<double>() -
                 bounds.at("minimum")[a].get<double>()) *
                scale;
    for (int a = 0; a < 3; ++a)
      for (int side = 0; side < 2; ++side) {
        if (m.generation_dimension == 2 && a == m.extrusion_axis)
          continue;
        const double expected = m.generation_dimension == 3
                                    ? span[(a + 1) % 3] * span[(a + 2) % 3]
                                    : span[3 - a - m.extrusion_axis];
        if (std::abs(coverage[static_cast<std::size_t>(2 * a + side)] -
                     expected) > expected * 1e-8)
          retainOuter = false;
      }
    if (!retainOuter) {
      m.cells.clear();
      planar.clear();
      core = coreCells(m, planar, cad, p, skin, query, grid);
      m.metadata["outer_boundary_retention_reason"] =
          "partial background coverage; kept full CAD boundary and cavity fill";
    }
  }
  m.metadata["retained_cartesian_outer_boundary"] = retainOuter;
  if (retainOuter) {
    auto removeBox = [&](std::vector<SkinFace> &shell) {
      shell.erase(std::remove_if(shell.begin(), shell.end(),
                                 [&](const auto &f) {
                                   return backgroundSide(m, f.vertices, p) >= 0;
                                 }),
                  shell.end());
    };
    removeBox(skin);
    removeBox(physicalSkin);
    for (const auto &face : core) {
      const int side = backgroundSide(m, face.vertices, p);
      if (side >= 0)
        physicalSkin.push_back(
            {face.vertices, outerPatches[static_cast<std::size_t>(side)]});
    }
    removeBox(core);
    m.metadata["outer_boundary_retention_reason"] =
        "complete Cartesian coverage of uniformly assigned outer box "
        "boundaries";
    std::cout << "Retained Cartesian outer boundary; fill only carved internal "
                 "cavities\n"
              << std::flush;
  }
  m.metadata["cartesian_core_is_axis_aligned"] = true;
  m.metadata["surface_size_m"] = h;
  m.metadata["surface_boundary_elements"] = physicalSkin.size();
  std::map<Index, std::array<double, 3>> fixedCore;
  if (m.generation_dimension == 3) {
    for (const auto &cell : m.cells)
      if (cell.region == CellRegion::cartesian)
        for (auto node : cell.vertices)
          fixedCore[node] = m.points_m.at(static_cast<std::size_t>(node));
  } else
    for (const auto &cell : planar)
      if (cell.region == CellRegion::cartesian)
        for (auto node : cell.vertices)
          fixedCore[node] = m.points_m.at(static_cast<std::size_t>(node));
  const double volumeSize =
      p.at("sizing").value("volume_size",
                           p.at("sizing").at("maximum_size").get<double>()) *
      scale;
  std::cout << "[4/7] Commit wall layers and fill the carved cavity\n"
            << std::flush;
  m.cells.insert(m.cells.end(), std::make_move_iterator(wallCells.begin()),
                 std::make_move_iterator(wallCells.end()));
  planar.insert(planar.end(), std::make_move_iterator(wallPlanar.begin()),
                std::make_move_iterator(wallPlanar.end()));
  if (m.generation_dimension == 2)
    fill2(m, planar, skin, core, p, volumeSize);
  else
    fill3(m, skin, core, p, volumeSize);
  std::size_t moved = 0;
  for (const auto &node : fixedCore)
    if (node.second != m.points_m.at(static_cast<std::size_t>(node.first)))
      ++moved;
  if (moved)
    throw std::runtime_error("cavity backend moved fixed Cartesian nodes");
  m.metadata["cartesian_core_nodes_moved"] = moved;
  resourceCheck();
  BoundaryMap boundary;
  if (m.generation_dimension == 2)
    extrude2(m, planar, physicalSkin, boundary);
  else
    for (const auto &f : physicalSkin)
      boundaryFace(boundary, f.vertices, f.patch);
  const auto maximum = p.at("resources").at("max_cells").get<std::size_t>();
  if (m.cells.size() > maximum)
    throw std::runtime_error(
        "generated hybrid mesh exceeds resources.max_cells");
  std::cout << "[5/7] Shared faces, cell shells and volume checks\n"
            << std::flush;
  compactNodes(m, boundary);
  buildMixedFaces(m, boundary);
  resourceCheck();
  GProp_GProps props;
  if (m.generation_dimension == 3)
    BRepGProp::VolumeProperties(cad.fluidDomain(), props);
  else {
    const auto shape = sectionShape(
        cad.fluidDomain(), m.extrusion_axis,
        p.at("meshing").at("section_coordinate").get<double>() * scale, p);
    BRepGProp::SurfaceProperties(shape, props);
  }
  const double expected =
      std::abs(props.Mass()) *
      (m.generation_dimension == 2 ? m.extrusion_thickness_m : 1);
  const double actual =
      std::accumulate(m.volumes_m3.begin(), m.volumes_m3.end(), 0.0);
  m.metadata["cad_measure_volume_m3"] = expected;
  m.metadata["relative_volume_error"] = std::abs(actual - expected) / expected;
  std::cout << "[6/7] Constrained automatic quality optimization\n"
            << std::flush;
  optimizeMixedMesh(m, p.value("optimization", Json::object()));
  std::cout << "[7/7] Bidirectional CAD fit audit\n" << std::flush;
  HexMesh proxy;
  proxy.points_m = m.points_m;
  for (const auto &f : m.faces)
    proxy.faces.push_back({f.vertices, f.owner, f.neighbour, f.patch});
  const double tolerance =
      p.at("quality").value("max_surface_deviation_m", h * .15);
  const double spacing = std::max(tolerance, std::min(h * .5, tolerance * 2));
  CadModel reference = cad;
  if (m.generation_dimension == 2) {
    const auto section = sectionShape(
        cad.fluidDomain(), m.extrusion_axis,
        p.at("meshing").at("section_coordinate").get<double>() * scale, p);
    gp_Trsf transform;
    Vec shift = Vec::Zero();
    shift[m.extrusion_axis] = -m.extrusion_thickness_m * .5;
    transform.SetTranslation(gp_Vec(shift.x(), shift.y(), shift.z()));
    const auto bottom =
        BRepBuilderAPI_Transform(section, transform, true).Shape();
    shift[m.extrusion_axis] = m.extrusion_thickness_m;
    reference = makeFluidCadModel(
        BRepPrimAPI_MakePrism(bottom, gp_Vec(shift.x(), shift.y(), shift.z()))
            .Shape());
  }
  const auto fit = auditSurfaceDistance(proxy, reference, spacing);
  m.metadata["geometry_distance_audit"] = {
      {"method", "bidirectional sampled exact B-Rep distance"},
      {"continuous_bound_certified", false},
      {"mesh_to_cad_max_m", fit.mesh_to_cad_max_m},
      {"cad_to_mesh_max_m", fit.cad_to_mesh_max_m},
      {"tolerance_m", tolerance},
      {"sample_spacing_m", fit.requested_sample_spacing_m},
      {"cad_deflection_m", fit.cad_tessellation_deflection_m},
      {"mesh_to_cad_samples", fit.mesh_to_cad_samples},
      {"cad_to_mesh_samples", fit.cad_to_mesh_samples},
      {"passed",
       std::max(fit.mesh_to_cad_max_m, fit.cad_to_mesh_max_m) <= tolerance}};
  if (std::max(fit.mesh_to_cad_max_m, fit.cad_to_mesh_max_m) > tolerance)
    throw std::runtime_error("CAD fit exceeds quality.max_surface_deviation_m; "
                             "refine surface_size/local surface sources: " +
                             m.metadata.at("geometry_distance_audit").dump());
  m.metadata["collision_audit"] =
      "PLC recovery, shared-face consistency, convex shells and sampled layer "
      "clearance; global cell overlap audit not evaluated";
  resourceCheck();
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  m.metadata["resources"] = {
      {"peak_resident_memory_mb", static_cast<double>(usage.ru_maxrss) / 1024},
      {"memory_budget_mb", p.at("resources").at("max_memory_mb")},
      {"max_cells", maximum},
      {"checks",
       "preflight grid estimate and resident-memory checks between stages"}};
  const auto report = mixedMeshReport(m);
  if (p.at("quality").at("mode") == "strict" &&
      !report.at("quality").at("face_quality_thresholds_passed").get<bool>())
    throw std::runtime_error("hybrid face quality thresholds failed: " +
                             report.at("quality").dump());
  return m;
}
} // namespace babelmesh
