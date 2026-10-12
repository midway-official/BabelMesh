#pragma once
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <QQuaternion>
#include <QVector3D>
#include <functional>
#include <map>
#include <nlohmann/json.hpp>
#include <vector>

class MeshView final : public QOpenGLWidget, protected QOpenGLFunctions {
public:
  explicit MeshView(QWidget *parent = nullptr);
  ~MeshView() override;
  void setCad(const nlohmann::json &preview);
  void setMesh(const nlohmann::json &checkpoint);
  void setClip(double fraction);
  void setClipAxis(int axis);
  void setPlanarClip(bool enabled);
  void setWireframe(bool enabled);
  void setOpacity(double value);
  void setPerspective(bool enabled);
  void setCellInset(double fraction);
  void setFrosted(bool enabled);
  void selectEntity(const std::string &id);
  void fit();
  void orient(int axis);
  std::size_t triangleCount() const { return triangles_.size() / 3; }
  std::function<void(std::string, int)> picked;

  bool renderingReady() const { return program_.isLinked(); }
  nlohmann::json displayState() const;

protected:
  void initializeGL() override;
  void resizeGL(int, int) override;
  void paintGL() override;
  void mousePressEvent(QMouseEvent *) override;
  void mouseMoveEvent(QMouseEvent *) override;
  void mouseReleaseEvent(QMouseEvent *) override;
  void wheelEvent(QWheelEvent *) override;

private:
  struct Face {
    std::vector<int> vertices;
    QVector3D colour;
    std::string entity;
    int cadFace{-1};
    int owner{-1}, neighbour{-1};
    bool wall{};
  };
  struct Vertex {
    float x, y, z, nx, ny, nz, r, g, b;
  };
  std::vector<QVector3D> points_, centres_;
  std::vector<Face> faces_;
  std::vector<std::pair<QVector3D, QVector3D>> cadEdges_;
  std::vector<Vertex> triangles_, lines_;
  std::map<int, std::size_t> typeCounts_;
  std::vector<bool> wallVertices_;
  std::vector<QVector3D> wallNormals_;
  std::map<int, std::size_t> capCounts_;
  std::size_t capVertices_{};
  float cutCoordinate_{}, maxCapPlaneError_{};
  QWidget *overlay_{};
  nlohmann::json source_;
  QOpenGLShaderProgram program_;
  QOpenGLBuffer buffer_{QOpenGLBuffer::VertexBuffer};
  QQuaternion rotation_;
  QVector3D centre_, pan_, boundsHalf_;
  float extent_{1}, zoom_{1};
  bool wireframe_{true}, dirty_{true}, meshMode_{false};
  bool perspective_{false}, frosted_{true}, planarClip_{true};
  float opacity_{.62f}, cellInset_{.035f};
  double clip_{1};
  int clipAxis_{2};
  std::string selected_;
  QPoint last_, press_;
  QMatrix4x4 matrix() const;
  float frameHalf() const;
  void rebuild();
  void pick(const QPoint &);
};
