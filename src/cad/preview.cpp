#include "babelmesh/cad.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <cmath>

namespace babelmesh {
nlohmann::json cadPreview(const CadModel &model) {
  using Json = nlohmann::json;
  Json result = {{"preview_version", 1}, {"entities", Json::array()}};
  for (const auto &entry : model.entities()) {
    const auto &entity = entry.second;
    Json item = {{"id", entity.id},
                 {"operation", entity.operation},
                 {"inputs", entity.inputs},
                 {"points_m", entity.surface_mesh.vertices},
                 {"triangles", Json::array()},
                 {"edges_m", Json::array()},
                 {"face_centres_m", Json::array()}};
    for (const auto &triangle : entity.surface_mesh.triangles)
      item["triangles"].push_back(
          {{"vertices", triangle.vertices}, {"face", triangle.face_index}});
    for (TopExp_Explorer ex(entity.shape, TopAbs_EDGE); ex.More(); ex.Next()) {
      BRepAdaptor_Curve curve(TopoDS::Edge(ex.Current()));
      const double first = curve.FirstParameter(), last = curve.LastParameter();
      if (!std::isfinite(first) || !std::isfinite(last))
        continue;
      Json points = Json::array();
      for (int i = 0; i <= 32; ++i) {
        const gp_Pnt p = curve.Value(first + (last - first) * i / 32);
        points.push_back({p.X(), p.Y(), p.Z()});
      }
      item["edges_m"].push_back(points);
    }
    for (TopExp_Explorer ex(entity.shape, TopAbs_FACE); ex.More(); ex.Next()) {
      GProp_GProps props;
      BRepGProp::SurfaceProperties(ex.Current(), props);
      const auto c = props.CentreOfMass();
      item["face_centres_m"].push_back({c.X(), c.Y(), c.Z()});
    }
    result["entities"].push_back(std::move(item));
  }
  return result;
}
} // namespace babelmesh
