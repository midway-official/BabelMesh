#include "mesh_view.hpp"
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace {
QVector3D point(const nlohmann::json &p) {
  return {p[0].get<float>(), p[1].get<float>(), p[2].get<float>()};
}
QVector3D typeColour(int type) {
  switch (type) {
  case 12:
    return {.44f, .67f, .74f};
  case 13:
    return {.85f, .68f, .40f};
  case 14:
    return {.77f, .49f, .49f};
  case 42:
    return {.49f, .70f, .59f};
  default:
    return {.63f, .59f, .77f};
  }
}
QVector3D wallColour() { return {.85f, .87f, .84f}; }
// A regular Qt child keeps labels independent of the OpenGL vertex/state cache.
class ViewOverlay final : public QWidget {
public:
  explicit ViewOverlay(MeshView *view) : QWidget(view), view_(view) {
    setAttribute(Qt::WA_TransparentForMouseEvents);
    setAttribute(Qt::WA_NoSystemBackground);
    setAttribute(Qt::WA_TranslucentBackground);
  }

protected:
  void paintEvent(QPaintEvent *) override;

private:
  MeshView *view_;
};
} // namespace
MeshView::MeshView(QWidget *parent) : QOpenGLWidget(parent) {
  setMinimumSize(500, 400);
  setFocusPolicy(Qt::StrongFocus);
  rotation_ = QQuaternion::fromEulerAngles(35.264f, -45, 0);
  overlay_ = new ViewOverlay(this);
}
MeshView::~MeshView() {
  makeCurrent();
  buffer_.destroy();
  doneCurrent();
}
void MeshView::setCad(const nlohmann::json &preview) {
  source_ = preview;
  meshMode_ = false;
  typeCounts_.clear();
  points_.clear();
  faces_.clear();
  centres_.clear();
  cadEdges_.clear();
  clip_ = 1;
  std::set<std::string> consumed;
  for (const auto &e : preview.at("entities"))
    for (const auto &id : e.at("inputs"))
      consumed.insert(id.get<std::string>());
  int entityIndex = 0;
  for (const auto &e : preview.at("entities")) {
    const auto id = e.at("id").get<std::string>();
    if (consumed.count(id) && id != selected_)
      continue;
    const auto offset = static_cast<int>(points_.size());
    for (const auto &p : e.at("points_m"))
      points_.push_back(point(p));
    const QVector3D colour = entityIndex++ % 2 ? QVector3D(.45f, .72f, .72f)
                                               : QVector3D(.30f, .62f, .85f);
    for (const auto &t : e.at("triangles")) {
      auto ids = t.at("vertices").get<std::vector<int>>();
      for (auto &v : ids)
        v += offset;
      faces_.push_back({ids, colour, id, t.at("face"), -1, -1});
    }
    for (const auto &edge : e.at("edges_m"))
      for (std::size_t i = 1; i < edge.size(); ++i)
        cadEdges_.push_back({point(edge[i - 1]), point(edge[i])});
  }
  fit();
  rebuild();
}
void MeshView::setMesh(const nlohmann::json &checkpoint) {
  source_ = checkpoint;
  meshMode_ = true;
  typeCounts_.clear();
  points_.clear();
  faces_.clear();
  centres_.clear();
  cadEdges_.clear();
  selected_.clear();
  for (const auto &p : checkpoint.at("points_m"))
    points_.push_back(point(p));
  wallVertices_.assign(points_.size(), false);
  wallNormals_.assign(points_.size(), {});
  std::vector<int> types;
  for (const auto &cell : checkpoint.at("cells")) {
    QVector3D centre;
    const auto ids = cell.at("vertices").get<std::vector<int>>();
    for (auto i : ids)
      centre += points_.at(static_cast<std::size_t>(i));
    centres_.push_back(centre / static_cast<float>(ids.size()));
    types.push_back(cell.at("type"));
    ++typeCounts_[types.back()];
  }
  for (const auto &f : checkpoint.at("faces")) {
    const int owner = f.at("owner");
    faces_.push_back({f.at("vertices").get<std::vector<int>>(),
                      typeColour(types.at(static_cast<std::size_t>(owner))), "",
                      -1, owner, f.at("neighbour")});
    const int patch = f.value("patch", -1);
    if (patch >= 0 &&
        checkpoint.at("patches")[static_cast<std::size_t>(patch)].at("kind") ==
            "wall") {
      faces_.back().wall = true;
      const auto &ring = faces_.back().vertices;
      QVector3D area;
      for (std::size_t i = 1; i + 1 < ring.size(); ++i)
        area += QVector3D::crossProduct(
            points_[static_cast<std::size_t>(ring[i])] -
                points_[static_cast<std::size_t>(ring[0])],
            points_[static_cast<std::size_t>(ring[i + 1])] -
                points_[static_cast<std::size_t>(ring[0])]);
      for (int id : faces_.back().vertices)
        wallVertices_.at(static_cast<std::size_t>(id)) = true;
      for (int id : ring)
        wallNormals_.at(static_cast<std::size_t>(id)) += area;
    }
  }
  for (auto &normal : wallNormals_)
    normal.normalize();
  clipAxis_ = checkpoint.at("generation_dimension") == 2
                  ? checkpoint.at("extrusion_axis").get<int>()
                  : 2;
  clip_ = checkpoint.at("generation_dimension") == 2 ? 1 : .5;
  perspective_ = false;
  orient(checkpoint.at("generation_dimension") == 2 ? clipAxis_ : -1);
  fit();
  rebuild();
}
void MeshView::fit() {
  QVector3D lo(1e30f, 1e30f, 1e30f), hi(-1e30f, -1e30f, -1e30f);
  auto include = [&](const QVector3D &p) {
    for (int a = 0; a < 3; ++a) {
      lo[a] = std::min(lo[a], p[a]);
      hi[a] = std::max(hi[a], p[a]);
    }
  };
  for (const auto &p : points_)
    include(p);
  for (const auto &e : cadEdges_) {
    include(e.first);
    include(e.second);
  }
  if (points_.empty() && cadEdges_.empty()) {
    centre_ = {};
    extent_ = 1;
    boundsHalf_ = QVector3D(.5f, .5f, .5f);
  } else {
    centre_ = (lo + hi) * .5f;
    extent_ =
        std::max({hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z(), 1e-9f});
    boundsHalf_ = (hi - lo) / (2 * extent_);
  }
  pan_ = {};
  zoom_ = 1;
  dirty_ = true;
  update();
}
void MeshView::orient(int axis) {
  rotation_ = axis == 2   ? QQuaternion()
              : axis == 1 ? QQuaternion::fromEulerAngles(90, 0, 0)
              : axis == 0 ? QQuaternion::fromEulerAngles(0, -90, 0)
                          : QQuaternion::fromEulerAngles(35.264f, -45, 0);
  dirty_ = true;
  update();
}
void MeshView::setClip(double fraction) {
  clip_ = std::clamp(fraction, 0.0, 1.0);
  rebuild();
}
void MeshView::setClipAxis(int axis) {
  clipAxis_ = std::clamp(axis, 0, 2);
  rebuild();
}
void MeshView::setPlanarClip(bool enabled) {
  planarClip_ = enabled;
  rebuild();
}
void MeshView::setWireframe(bool enabled) {
  wireframe_ = enabled;
  update();
}
void MeshView::setOpacity(double value) {
  opacity_ = static_cast<float>(std::clamp(value, .12, 1.0));
  update();
}
void MeshView::setPerspective(bool enabled) {
  perspective_ = enabled;
  update();
}
void MeshView::setCellInset(double fraction) {
  cellInset_ = static_cast<float>(std::clamp(fraction, 0.0, .12));
  rebuild();
}
void MeshView::setFrosted(bool enabled) {
  frosted_ = enabled;
  rebuild();
}
nlohmann::json MeshView::displayState() const {
  auto types = nlohmann::json::object();
  for (const auto &entry : typeCounts_)
    types[std::to_string(entry.first)] = entry.second;
  auto caps = nlohmann::json::object();
  for (const auto &entry : capCounts_)
    caps[std::to_string(entry.first)] = entry.second;
  return {{"projection", perspective_ ? "perspective" : "orthographic"},
          {"material", frosted_ ? "frosted" : "solid"},
          {"opacity", opacity_},
          {"cell_inset_fraction", cellInset_},
          {"effective_cell_inset_fraction",
           meshMode_ && planarClip_ && clip_ < .999 ? 0.f : cellInset_},
          {"view_quaternion",
           {rotation_.scalar(), rotation_.x(), rotation_.y(), rotation_.z()}},
          {"rendered_triangles", triangleCount()},
          {"cells", centres_.size()},
          {"model", meshMode_ ? "mesh" : "cad"},
          {"cell_type_counts", types},
          {"edges_enabled", wireframe_},
          {"edge_colour", "#000000"},
          {"physical_wall_colour", "#d9ded6"},
          {"edge_style", "solid"},
          {"hidden_edges_suppressed", true},
          {"clip_fraction", clip_},
          {"clip_mode", planarClip_ ? "plane" : "whole_cells"},
          {"clip_axis", clipAxis_},
          {"clip_coordinate_m", cutCoordinate_},
          {"cap_cell_type_counts", caps},
          {"cap_vertices", capVertices_},
          {"max_cap_plane_error_m", maxCapPlaneError_},
          {"depth_sorted_transparency", true},
          {"generation_dimension",
           meshMode_ ? source_.value("generation_dimension", 3) : 3}};
}
void MeshView::selectEntity(const std::string &id) {
  selected_ = id;
  if (!meshMode_ && source_.contains("entities"))
    setCad(source_);
  else
    rebuild();
}
void MeshView::rebuild() {
  triangles_.clear();
  lines_.clear();
  capCounts_.clear();
  capVertices_ = 0;
  maxCapPlaneError_ = 0;
  const float epsilon = std::max(extent_ * 1e-6f, 1e-12f);
  const bool plane = meshMode_ && planarClip_ && clip_ < .999;
  float lo = 1e30f, hi = -1e30f;
  if (planarClip_ || centres_.empty())
    for (const auto &p : points_) {
      lo = std::min(lo, p[clipAxis_]);
      hi = std::max(hi, p[clipAxis_]);
    }
  else
    for (const auto &c : centres_) {
      lo = std::min(lo, c[clipAxis_]);
      hi = std::max(hi, c[clipAxis_]);
    }
  const float cut = lo + (hi - lo) * static_cast<float>(clip_);
  cutCoordinate_ = points_.empty() ? 0 : cut;
  std::vector<bool> visible(centres_.size()), cutCells(centres_.size());
  if (meshMode_)
    for (std::size_t i = 0; i < centres_.size(); ++i) {
      float minimum = 1e30f, maximum = -1e30f;
      for (int id : source_.at("cells")[i].at("vertices")) {
        minimum = std::min(minimum,
                           points_.at(static_cast<std::size_t>(id))[clipAxis_]);
        maximum = std::max(maximum,
                           points_.at(static_cast<std::size_t>(id))[clipAxis_]);
      }
      visible[i] = clip_ >= .999 || (plane ? minimum < cut - epsilon
                                           : centres_[i][clipAxis_] <= cut);
      cutCells[i] = plane && visible[i] && maximum >= cut - epsilon;
    }
  auto vertex = [&](const QVector3D &p, const QVector3D &n,
                    const QVector3D &c) {
    const QVector3D q = (p - centre_) / extent_;
    return Vertex{q.x(), q.y(), q.z(), n.x(), n.y(),
                  n.z(), c.x(), c.y(), c.z()};
  };
  struct ClipPoint {
    QVector3D position;
    bool wall{};
    QVector3D normal;
  };
  std::vector<std::vector<ClipPoint>> caps(centres_.size());
  auto appendPolygon = [&](const std::vector<ClipPoint> &polygon, int cell,
                           const QVector3D &colour, bool drawEdges,
                           bool wallSurface = false) {
    if (polygon.size() < 3)
      return;
    std::vector<QVector3D> ring;
    for (const auto &v : polygon) {
      auto p = v.position;
      if (meshMode_ && frosted_ && cellInset_ > 0 && !v.wall && !plane) {
        const auto c = centres_.at(static_cast<std::size_t>(cell));
        p = c + (p - c) * (1 - cellInset_);
      }
      ring.push_back(p);
    }
    QVector3D area;
    for (std::size_t j = 1; j + 1 < ring.size(); ++j)
      area += QVector3D::crossProduct(ring[j] - ring[0], ring[j + 1] - ring[0]);
    if (area.lengthSquared() <= epsilon * epsilon * epsilon * epsilon)
      return;
    const auto normal = area.normalized();
    auto shadingNormal = [&](std::size_t i) {
      const auto n = polygon[i].normal.normalized();
      // Smooth curvature shading only; preserve facets at sharp CAD corners.
      return wallSurface && QVector3D::dotProduct(n, normal) > .82f ? n
                                                                    : normal;
    };
    for (std::size_t j = 1; j + 1 < ring.size(); ++j) {
      triangles_.push_back(vertex(ring[0], shadingNormal(0), colour));
      triangles_.push_back(vertex(ring[j], shadingNormal(j), colour));
      triangles_.push_back(vertex(ring[j + 1], shadingNormal(j + 1), colour));
    }
    if (drawEdges)
      for (std::size_t j = 0; j < ring.size(); ++j) {
        lines_.push_back(vertex(ring[j], {}, QVector3D(0, 0, 0)));
        lines_.push_back(
            vertex(ring[(j + 1) % ring.size()], {}, QVector3D(0, 0, 0)));
      }
  };
  auto clipped = [&](const Face &face) {
    std::vector<ClipPoint> ring;
    for (int id : face.vertices)
      ring.push_back(
          {points_.at(static_cast<std::size_t>(id)),
           meshMode_ && wallVertices_.at(static_cast<std::size_t>(id)),
           {}});
    if (meshMode_)
      for (std::size_t i = 0; i < ring.size(); ++i)
        ring[i].normal =
            wallNormals_.at(static_cast<std::size_t>(face.vertices[i]));
    if (!plane)
      return ring;
    std::vector<ClipPoint> result;
    for (std::size_t j = 0; j < ring.size(); ++j) {
      const auto &a = ring[j], &b = ring[(j + 1) % ring.size()];
      const bool insideA = a.position[clipAxis_] <= cut;
      const bool insideB = b.position[clipAxis_] <= cut;
      if (insideA)
        result.push_back(a);
      if (insideA != insideB) {
        const float t = (cut - a.position[clipAxis_]) /
                        (b.position[clipAxis_] - a.position[clipAxis_]);
        auto p = a.position + (b.position - a.position) * t;
        p[clipAxis_] = cut;
        result.push_back(
            {p, a.wall && b.wall, a.normal + (b.normal - a.normal) * t});
      }
    }
    return result;
  };
  for (const auto &face : faces_) {
    auto ring = clipped(face);
    if (!meshMode_) {
      const auto colour = face.entity == selected_ && !selected_.empty()
                              ? QVector3D(.88f, .64f, .30f)
                              : face.colour;
      appendPolygon(ring, -1, colour, false);
      continue;
    }
    if (plane)
      for (int cell : {face.owner, face.neighbour})
        if (cell >= 0 && cutCells[static_cast<std::size_t>(cell)])
          for (const auto &v : ring)
            if (std::abs(v.position[clipAxis_] - cut) <= epsilon) {
              auto &cap = caps[static_cast<std::size_t>(cell)];
              const auto duplicate =
                  std::find_if(cap.begin(), cap.end(), [&](const auto &p) {
                    return (p.position - v.position).lengthSquared() <=
                           epsilon * epsilon;
                  });
              if (duplicate == cap.end())
                cap.push_back(v);
              else
                duplicate->wall = duplicate->wall || v.wall;
            }
    const bool owner = visible[static_cast<std::size_t>(face.owner)];
    const bool neighbour = face.neighbour >= 0 &&
                           visible[static_cast<std::size_t>(face.neighbour)];
    if (owner == neighbour || ring.size() < 3)
      continue;
    // A face coincident with the cut is emitted once through the cap path.
    if (plane && std::all_of(ring.begin(), ring.end(), [&](const auto &v) {
          return std::abs(v.position[clipAxis_] - cut) <= epsilon;
        }))
      continue;
    const int cell = owner ? face.owner : face.neighbour;
    if (!owner)
      std::reverse(ring.begin(), ring.end());
    const int type =
        source_.at("cells")[static_cast<std::size_t>(cell)].at("type");
    // Render the exposed clipped shell only. Buried cell faces would add ghost
    // triangles to translucent cavities despite having no visible boundary.
    const bool wallSurface =
        face.wall && source_.at("generation_dimension") == 3;
    appendPolygon(ring, cell, wallSurface ? wallColour() : typeColour(type),
                  true, wallSurface);
  }
  for (std::size_t cell = 0; cell < caps.size(); ++cell) {
    auto &cap = caps[cell];
    if (cap.size() < 3)
      continue;
    QVector3D centre;
    for (const auto &p : cap)
      centre += p.position;
    centre /= static_cast<float>(cap.size());
    const int a = (clipAxis_ + 1) % 3, b = (clipAxis_ + 2) % 3;
    std::sort(cap.begin(), cap.end(), [&](const auto &p, const auto &q) {
      return std::atan2(p.position[b] - centre[b], p.position[a] - centre[a]) <
             std::atan2(q.position[b] - centre[b], q.position[a] - centre[a]);
    });
    const int type = source_.at("cells")[cell].at("type");
    ++capCounts_[type];
    capVertices_ += cap.size();
    for (const auto &p : cap)
      maxCapPlaneError_ =
          std::max(maxCapPlaneError_, std::abs(p.position[clipAxis_] - cut));
    appendPolygon(cap, static_cast<int>(cell), typeColour(type), true);
  }
  for (const auto &edge : cadEdges_) {
    lines_.push_back(vertex(edge.first, {}, QVector3D(.09f, .15f, .21f)));
    lines_.push_back(vertex(edge.second, {}, QVector3D(.09f, .15f, .21f)));
  }
  dirty_ = true;
  update();
}
void MeshView::initializeGL() {
  initializeOpenGLFunctions();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);
  program_.addShaderFromSourceCode(
      QOpenGLShader::Vertex,
      "#version 120\nattribute vec3 position; attribute vec3 normal; attribute "
      "vec3 colour; "
      "uniform mat4 mvp; uniform mat4 rotation; varying vec3 n; varying vec3 "
      "c; void main(){gl_Position=mvp*vec4(position,1.0); "
      "n=mat3(rotation)*normal;c=colour;}");
  program_.addShaderFromSourceCode(
      QOpenGLShader::Fragment,
      "#version 120\nvarying vec3 n; varying vec3 c; uniform float shade; "
      "uniform float opacity; uniform float frosted; void main(){"
      "vec3 rgb=c; if(shade>0.5){vec3 nn=normalize(n);"
      "float light=0.76+0.20*abs(dot(nn,normalize(vec3(-0.4,0.6,1.0))))"
      "+0.07*abs(dot(nn,normalize(vec3(0.8,0.1,0.5))));"
      "float rim=pow(1.0-abs(nn.z),2.0);"
      "rgb=mix(c*light,vec3(0.93,0.96,0.98),0.12+0.16*rim);"
      "float "
      "grain=fract(sin(dot(gl_FragCoord.xy,vec2(12.9898,78.233)))*43758.5453);"
      "rgb+=(grain-0.5)*0.024*frosted;}"
      "gl_FragColor=vec4(rgb,opacity);}");
  if (!program_.link())
    qWarning("Mesh viewport shader failed: %s", qPrintable(program_.log()));
  buffer_.create();
  dirty_ = true;
}
void MeshView::resizeGL(int, int) { overlay_->setGeometry(rect()); }
float MeshView::frameHalf() const {
  const float aspect = static_cast<float>(width()) / std::max(height(), 1);
  float half = .1f;
  for (int i = 0; i < 8; ++i) {
    const QVector3D p = rotation_.rotatedVector(QVector3D(
        (i & 1 ? 1 : -1) * boundsHalf_.x(), (i & 2 ? 1 : -1) * boundsHalf_.y(),
        (i & 4 ? 1 : -1) * boundsHalf_.z()));
    half = std::max({half, std::abs(p.x()) / aspect, std::abs(p.y())});
  }
  return half * 1.32f;
}
QMatrix4x4 MeshView::matrix() const {
  QMatrix4x4 projection, model;
  const float aspect = static_cast<float>(width()) / std::max(height(), 1);
  const float half = frameHalf() / zoom_;
  if (perspective_) {
    projection.perspective(35, aspect, .02f, 100.f);
    model.translate(0, 0, -(half / std::tan(17.5f * .017453293f) + .8f));
  } else
    projection.ortho(-half * aspect, half * aspect, -half, half, -10, 10);
  model.translate(pan_);
  model.rotate(rotation_);
  return projection * model;
}
void MeshView::paintGL() {
  glDisable(GL_STENCIL_TEST);
  glDisable(GL_SCISSOR_TEST);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDepthMask(GL_TRUE);
  glDepthFunc(GL_LESS);
  glBlendEquation(GL_FUNC_ADD);
  glClearColor(.94f, .96f, .97f, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (!program_.isLinked()) {
    return;
  }
  glEnable(GL_DEPTH_TEST);
  glClear(GL_DEPTH_BUFFER_BIT);
  glEnable(GL_BLEND);
  // Qt composites this framebuffer as premultiplied RGBA. Preserve the opaque
  // clear alpha while blending translucent surfaces, or captures unpremultiply
  // already blended RGB again and produce saturated, incorrect colours.
  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE,
                      GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_MULTISAMPLE);
  if (meshMode_) {
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
  } else
    glDisable(GL_CULL_FACE);
  program_.bind();
  buffer_.bind();
  if (dirty_) {
    const auto bytes = (triangles_.size() + lines_.size()) * sizeof(Vertex);
    if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      buffer_.release();
      program_.release();
      return;
    }
    // Back-to-front triangle order is independent of insertion/owner order.
    std::vector<std::size_t> order(triangles_.size() / 3);
    std::iota(order.begin(), order.end(), 0);
    auto depth = [&](std::size_t triangle) {
      const auto &a = triangles_[3 * triangle],
                 &b = triangles_[3 * triangle + 1],
                 &c = triangles_[3 * triangle + 2];
      return rotation_
          .rotatedVector(
              QVector3D(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z))
          .z();
    };
    std::stable_sort(order.begin(), order.end(),
                     [&](auto a, auto b) { return depth(a) < depth(b); });
    std::vector<Vertex> sorted;
    sorted.reserve(triangles_.size());
    for (auto i : order)
      sorted.insert(
          sorted.end(), triangles_.begin() + static_cast<std::ptrdiff_t>(3 * i),
          triangles_.begin() + static_cast<std::ptrdiff_t>(3 * i + 3));
    buffer_.allocate(static_cast<int>(bytes));
    buffer_.write(0, sorted.data(),
                  static_cast<int>(triangles_.size() * sizeof(Vertex)));
    buffer_.write(static_cast<int>(triangles_.size() * sizeof(Vertex)),
                  lines_.data(),
                  static_cast<int>(lines_.size() * sizeof(Vertex)));
    dirty_ = false;
  }
  QMatrix4x4 rotation;
  rotation.rotate(rotation_);
  program_.setUniformValue("mvp", matrix());
  program_.setUniformValue("rotation", rotation);
  for (const auto &a : std::vector<std::pair<const char *, int>>{
           {"position", 0}, {"normal", 3}, {"colour", 6}}) {
    const int loc = program_.attributeLocation(a.first);
    program_.enableAttributeArray(loc);
    program_.setAttributeBuffer(loc, GL_FLOAT,
                                a.second * static_cast<int>(sizeof(float)), 3,
                                sizeof(Vertex));
  }
  program_.setUniformValue("shade", 1.0f);
  program_.setUniformValue("frosted", meshMode_ && frosted_ ? 1.f : 0.f);
  program_.setUniformValue("opacity", meshMode_ && frosted_ ? opacity_ : 1.f);
  glDepthMask(GL_FALSE);
  glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(triangles_.size()));
  // Depth-only surface pass suppresses hidden back and internal edge clutter,
  // while leaving the already-composited translucent colour untouched.
  glDepthMask(GL_TRUE);
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);
  glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(triangles_.size()));
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDisable(GL_POLYGON_OFFSET_FILL);
  if (wireframe_ || triangles_.empty()) {
    glDepthFunc(GL_LEQUAL);
    program_.setUniformValue("shade", 0.0f);
    program_.setUniformValue("opacity", 1.f);
    glLineWidth(meshMode_ ? 1.2f : 1.f);
    glDrawArrays(GL_LINES, static_cast<int>(triangles_.size()),
                 static_cast<int>(lines_.size()));
  }
  glDepthFunc(GL_LESS);
  glDisable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  for (const auto *attribute : {"position", "normal", "colour"})
    program_.disableAttributeArray(program_.attributeLocation(attribute));
  buffer_.release();
  program_.release();
  overlay_->update();
}
void ViewOverlay::paintEvent(QPaintEvent *) {
  const auto state = view_->displayState();
  const bool meshMode = state.at("model") == "mesh";
  const bool planar = meshMode && state.at("generation_dimension") == 2;
  const bool perspective = state.at("projection") == "perspective";
  const bool frosted = state.at("material") == "frosted";
  const bool cut =
      meshMode && !planar && state.at("clip_fraction").get<double>() < .999;
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);
  painter.setPen(Qt::NoPen);
  painter.setBrush(QColor(255, 255, 255, 222));
  painter.drawRoundedRect(QRectF(20, 18, 268, cut ? 84 : 64), 12, 12);
  painter.setPen(QColor("#243b47"));
  QFont title = font();
  title.setPointSize(11);
  title.setWeight(QFont::DemiBold);
  painter.setFont(title);
  painter.drawText(36, 44,
                   meshMode ? (planar ? QStringLiteral("二维截面 · 单层体网格")
                                      : QStringLiteral("三维混合体网格"))
                            : QStringLiteral("CAD 工作视图"));
  QFont caption = font();
  caption.setPointSize(9);
  painter.setFont(caption);
  painter.setPen(QColor("#708591"));
  const auto projection =
      perspective ? QStringLiteral("透视") : QStringLiteral("正交");
  painter.drawText(
      36, 65,
      meshMode ? QStringLiteral("%1 单元  ·  %2  ·  %3")
                     .arg(state.at("cells").get<qulonglong>())
                     .arg(projection, frosted ? QStringLiteral("磨砂")
                                              : QStringLiteral("实体"))
               : QStringLiteral("选择面或实体  ·  %1视图").arg(projection));
  if (cut)
    painter.drawText(
        36, 85,
        state.at("clip_mode") == "plane"
            ? QStringLiteral("平面剖切 · %1 = %2 m")
                  .arg(QChar("XYZ"[state.at("clip_axis").get<int>()]))
                  .arg(state.at("clip_coordinate_m").get<double>(), 0, 'g', 5)
            : QStringLiteral("按整单元筛选 · %1%")
                  .arg(state.at("clip_fraction").get<double>() * 100, 0, 'f',
                       0));
  if (meshMode) {
    const std::vector<std::pair<int, QString>> legend = {
        {12, QStringLiteral("六面体")},
        {13, QStringLiteral("三棱柱")},
        {10, QStringLiteral("四面体")},
        {14, QStringLiteral("四棱锥")},
        {42, QStringLiteral("多面体")}};
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(255, 255, 255, 218));
    const bool wallLegend = state.at("generation_dimension") == 3;
    painter.drawRoundedRect(
        QRectF(20, height() - 66,
               32 + 85 * (state.at("cell_type_counts").size() +
                          (wallLegend ? 1 : 0)),
               40),
        10, 10);
    int x = 36;
    for (const auto &entry : legend) {
      if (!state.at("cell_type_counts").contains(std::to_string(entry.first)))
        continue;
      const auto c = typeColour(entry.first);
      painter.setPen(Qt::NoPen);
      painter.setBrush(QColor::fromRgbF(c.x(), c.y(), c.z()));
      painter.drawEllipse(QPointF(x + 5, height() - 46), 4, 4);
      painter.setPen(QColor("#556e7b"));
      painter.drawText(x + 17, height() - 42, entry.second);
      x += 85;
    }
    if (wallLegend) {
      const auto c = wallColour();
      painter.setPen(Qt::NoPen);
      painter.setBrush(QColor::fromRgbF(c.x(), c.y(), c.z()));
      painter.drawEllipse(QPointF(x + 5, height() - 46), 4, 4);
      painter.setPen(QColor("#556e7b"));
      painter.drawText(x + 17, height() - 42, QStringLiteral("实体壁面"));
    }
  }
  // Camera-relative world-axis triad, with labels rather than a decorative
  // cube.
  const QPointF origin(width() - 56, height() - 60);
  const auto &q = state.at("view_quaternion");
  const QQuaternion rotation(q[0].get<float>(), q[1].get<float>(),
                             q[2].get<float>(), q[3].get<float>());
  for (int a = 0; a < 3; ++a) {
    QVector3D world;
    world[a] = 1;
    const auto d = rotation.rotatedVector(world);
    const QPointF endpoint = origin + QPointF(d.x() * 30, -d.y() * 30);
    const QColor colour(a == 0 ? "#bd7375" : a == 1 ? "#659686" : "#728eb7");
    painter.setPen(QPen(colour, 2, Qt::SolidLine, Qt::RoundCap));
    painter.drawLine(origin, endpoint);
    painter.drawText(endpoint + QPointF(4, -3), QString("XYZ"[a]));
  }
}
void MeshView::mousePressEvent(QMouseEvent *e) {
  last_ = press_ = e->position().toPoint();
  setFocus();
}
void MeshView::mouseMoveEvent(QMouseEvent *e) {
  const QPoint next = e->position().toPoint(), d = next - last_;
  last_ = next;
  if (e->buttons() & Qt::LeftButton)
    rotation_ =
        QQuaternion::fromEulerAngles(d.y() * .5f, d.x() * .5f, 0) * rotation_;
  if (e->buttons() & Qt::LeftButton)
    dirty_ = true;
  if (e->buttons() & Qt::RightButton)
    pan_ += QVector3D(d.x() * .0018f / zoom_, -d.y() * .0018f / zoom_, 0);
  update();
}
void MeshView::mouseReleaseEvent(QMouseEvent *e) {
  if ((e->position().toPoint() - press_).manhattanLength() < 4 &&
      e->button() == Qt::LeftButton)
    pick(e->position().toPoint());
}
void MeshView::wheelEvent(QWheelEvent *e) {
  zoom_ = std::clamp(
      zoom_ * std::pow(1.0015f, static_cast<float>(e->angleDelta().y())), .08f,
      40.0f);
  update();
}
void MeshView::pick(const QPoint &mouse) {
  if (meshMode_ || !picked)
    return;
  double best = 1e30;
  std::string entity;
  int face = -1;
  const auto transform = matrix();
  auto project = [&](int id) {
    const QVector3D p = transform.map(
        (points_.at(static_cast<std::size_t>(id)) - centre_) / extent_);
    return QVector3D((p.x() + 1) * width() / 2, (1 - p.y()) * height() / 2,
                     p.z());
  };
  for (const auto &f : faces_)
    for (std::size_t j = 1; j + 1 < f.vertices.size(); ++j) {
      const auto a = project(f.vertices[0]), b = project(f.vertices[j]),
                 c = project(f.vertices[j + 1]);
      const double det =
          (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
      if (std::abs(det) < 1e-8)
        continue;
      const double u = ((b.y() - c.y()) * (mouse.x() - c.x()) +
                        (c.x() - b.x()) * (mouse.y() - c.y())) /
                       det;
      const double v = ((c.y() - a.y()) * (mouse.x() - c.x()) +
                        (a.x() - c.x()) * (mouse.y() - c.y())) /
                       det;
      const double w = 1 - u - v, z = u * a.z() + v * b.z() + w * c.z();
      if (u >= 0 && v >= 0 && w >= 0 && z < best) {
        best = z;
        entity = f.entity;
        face = f.cadFace;
      }
    }
  if (!entity.empty())
    picked(entity, face);
}
