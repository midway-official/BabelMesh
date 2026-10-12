#pragma once
#include <QMainWindow>
#include <QProcess>
#include <QTemporaryDir>
#include <array>
#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

class MeshView;
class QTreeWidget;
class QPlainTextEdit;
class QComboBox;
class QDoubleSpinBox;
class QSpinBox;
class QCheckBox;
class QLabel;
class QProgressBar;
class QSlider;

class Studio final : public QMainWindow {
public:
  explicit Studio(const QString &cli, QWidget *parent = nullptr);
  ~Studio() override;
  bool loadProject(const QString &path);
  void smoke(const QString &project, const QString &mesh,
             const QString &screenshot);
  void previewDefault() { inspect(); }
  void setSmokePerspective(bool enabled) { smokePerspective_ = enabled; }
  void setSmokeViewAxis(int axis) { smokeViewAxis_ = axis; }

protected:
  void closeEvent(QCloseEvent *) override;

private:
  using Json = nlohmann::json;
  Json project_, preview_;
  std::vector<Json> undo_, redo_;
  QString cli_, projectPath_, jobOutput_, smokeMesh_, smokeScreenshot_,
      currentMesh_;
  std::unique_ptr<QTemporaryDir> jobDirectory_;
  QProcess *process_{};
  MeshView *view_{};
  QTreeWidget *tree_{};
  QPlainTextEdit *log_{}, *jsonEditor_{};
  QComboBox *division_{};
  QComboBox *viewProjection_{};
  QSlider *clipControl_{};
  QComboBox *clipModeControl_{}, *clipAxisControl_{};
  QComboBox *dimension_{}, *axis_{}, *quality_{}, *layerGroup_{};
  std::array<QSpinBox *, 3> counts_{};
  std::array<QDoubleSpinBox *, 3> minimum_{}, maximum_{}, spacing_{};
  QDoubleSpinBox *surface_{}, *coordinate_{}, *thickness_{}, *buffer_{},
      *firstLayer_{}, *growth_{}, *volumeSize_{}, *minimumSize_{},
      *cadTolerance_{}, *maxMemory_{};
  QSpinBox *layers_{}, *maxLevel_{}, *optIterations_{}, *maxCells_{},
      *threads_{};
  QTreeWidget *sources_{};
  QCheckBox *layerEnabled_{}, *overwrite_{}, *optEnabled_{}, *mergeOpt_{},
      *topologyOpt_{};
  QLabel *summary_{}, *layerSummary_{};
  QProgressBar *progress_{};
  bool modified_{false}, smokeMode_{false}, cancelled_{false},
      workflowMode_{false}, cancellationTest_{false}, smokePerspective_{false};
  QString jobKind_;
  int serial_{}, revision_{}, jobRevision_{}, smokeViewAxis_{-2};
  std::string selected_;
  int selectedFace_{-1};
  void buildUi();
  void showMesh(const Json &mesh);
  bool captureSmoke();
  void fresh(int dimension, bool internal = false);
  void refresh();
  bool changed(Json next, bool preview = true);
  bool applyParameters();
  void sizeSource(bool edit);
  void updateLayerSummary();
  void loadLayerParameters();
  bool save(bool saveAs = false);
  void primitive(const std::string &type);
  void sketch();
  void operation(const std::string &type);
  void editOperation();
  void importStep();
  void physicalGroup();
  void removeSelected();
  void inspect();
  void generate();
  void convertMesh();
  void optimizeMesh();
  void exportCad();
  void launchProcess(const QString &kind, const QStringList &arguments);
  void startJob(const QString &kind);
  void finishJob(int exitCode, QProcess::ExitStatus status);
  void cancel();
  void error(const QString &message);
  Json jobProject() const;
  QString uniqueId(const QString &prefix) const;
  bool readJson(const QString &path, Json &target);

public:
  void workflowSmoke(const QString &project, const QString &directory,
                     const QString &screenshot);
  void cancellationSmoke(const QString &project, const QString &directory);
};
