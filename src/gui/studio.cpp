#include "studio.hpp"
#include "babelmesh/hybrid_mesher.hpp"
#include "mesh_view.hpp"
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMenuBar>
#include <QMessageBox>
#include <QPainter>
#include <QPainterPath>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QSaveFile>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <set>

namespace {
using Json = nlohmann::json;
QIcon toolIcon(const QString &kind) {
  QPixmap pixmap(48, 48);
  pixmap.setDevicePixelRatio(2);
  pixmap.fill(Qt::transparent);
  QPainter p(&pixmap);
  p.setRenderHint(QPainter::Antialiasing);
  p.setPen(
      QPen(QColor("#587681"), 1.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  if (kind == "play") {
    p.setPen(Qt::NoPen);
    p.setBrush(Qt::white);
    p.drawPolygon(QPolygonF{QPointF(8, 5), QPointF(18, 12), QPointF(8, 19)});
  } else if (kind == "stop") {
    p.drawRoundedRect(QRectF(6, 6, 12, 12), 2, 2);
  } else if (kind == "import") {
    p.drawRoundedRect(QRectF(4, 4, 16, 16), 3, 3);
    p.drawLine(12, 2, 12, 13);
    p.drawLine(8, 9, 12, 13);
    p.drawLine(16, 9, 12, 13);
  } else if (kind == "sphere") {
    p.drawEllipse(QRectF(4, 4, 16, 16));
    p.drawEllipse(QRectF(8, 4, 8, 16));
    p.drawLine(4, 12, 20, 12);
  } else if (kind == "cylinder") {
    p.drawEllipse(QRectF(5, 4, 14, 6));
    p.drawLine(5, 7, 5, 17);
    p.drawLine(19, 7, 19, 17);
    p.drawArc(QRectF(5, 14, 14, 6), 180 * 16, 180 * 16);
  } else if (kind == "fit") {
    for (int x : {4, 20})
      for (int y : {4, 20}) {
        p.drawLine(x, y, x + (x == 4 ? 5 : -5), y);
        p.drawLine(x, y, x, y + (y == 4 ? 5 : -5));
      }
  } else if (kind == "refresh") {
    p.drawArc(QRectF(5, 5, 14, 14), 35 * 16, 290 * 16);
    p.drawLine(20, 4, 20, 10);
    p.drawLine(20, 10, 14, 10);
  } else {
    p.drawPolygon(QPolygonF{QPointF(12, 3), QPointF(21, 8), QPointF(21, 17),
                            QPointF(12, 22), QPointF(3, 17), QPointF(3, 8)});
    p.drawLine(3, 8, 12, 13);
    p.drawLine(21, 8, 12, 13);
    p.drawLine(12, 13, 12, 22);
    p.drawLine(12, 3, 12, 13);
  }
  return QIcon(pixmap);
}
QDoubleSpinBox *real(QWidget *parent, double minimum = -1e9,
                     double maximum = 1e9) {
  auto *s = new QDoubleSpinBox(parent);
  s->setDecimals(7);
  s->setRange(minimum, maximum);
  s->setSingleStep(.01);
  s->setKeyboardTracking(false);
  return s;
}
QWidget *vectorRow(QWidget *parent, std::array<QDoubleSpinBox *, 3> &widgets) {
  auto *box = new QWidget(parent);
  auto *layout = new QVBoxLayout(box);
  layout->setContentsMargins(0, 0, 0, 0);
  for (int i = 0; i < 3; ++i) {
    widgets[static_cast<std::size_t>(i)] = real(box);
    widgets[static_cast<std::size_t>(i)]->setPrefix(
        QString("%1 ").arg("xyz"[i]));
    layout->addWidget(widgets[static_cast<std::size_t>(i)]);
  }
  return box;
}
double scale(const Json &p) {
  const auto u = p.at("units").at("length").get<std::string>();
  return u == "mm" ? .001 : u == "cm" ? .01 : u == "inch" ? .0254 : 1;
}
bool writeJson(const QString &path, const Json &json) {
  QSaveFile file(path);
  if (!file.open(QIODevice::WriteOnly))
    return false;
  const auto text = json.dump(2);
  file.write(text.data(), static_cast<qint64>(text.size()));
  file.write("\n");
  return file.commit();
}
} // namespace
Studio::Studio(const QString &cli, QWidget *parent)
    : QMainWindow(parent), cli_(cli) {
  buildUi();
  fresh(3);
}
Studio::~Studio() {
  if (process_ && process_->state() != QProcess::NotRunning) {
    process_->kill();
    process_->waitForFinished(1500);
  }
}
void Studio::buildUi() {
  setWindowTitle(QStringLiteral("BabelMesh · CAD 与笛卡尔混合网格"));
  resize(1560, 980);
  setStyleSheet(QStringLiteral(R"(
    QMainWindow { background:#edf2f5; color:#304956; }
    QWidget { color:#304956; font-size:12px; }
    QMenuBar { background:#fff; padding:4px 12px; border-bottom:1px solid #e3eaee; }
    QMenuBar::item { padding:5px 10px; border-radius:5px; }
    QMenuBar::item:selected, QMenu::item:selected { background:#e4f1f1; color:#167d83; }
    QMenu { background:#fff; border:1px solid #dce5ea; padding:6px; }
    QMenu::item { padding:8px 22px; }
    QToolBar { background:#fff; spacing:7px; padding:9px 14px; border:0; border-bottom:1px solid #dce5ea; }
    QToolBar::separator { background:#e4ebef; width:1px; margin:6px 9px; }
    QToolButton { border:1px solid transparent; border-radius:7px; padding:7px 10px; color:#526c79; }
    QToolButton:hover { background:#edf4f6; border-color:#dbe6ea; }
    QToolButton[primary="true"] { color:white; background:#218b90; padding:8px 15px; font-weight:600; }
    QToolButton[primary="true"]:hover { background:#187b80; }
    QGroupBox { background:#fff; font-weight:600; border:1px solid #dce5ea; border-radius:10px; margin-top:22px; padding:18px 12px 12px; }
    QGroupBox::title { subcontrol-origin:margin; left:15px; padding:0 4px; color:#496572; }
    QLabel { background:transparent; }
    QAbstractSpinBox, QComboBox, QLineEdit { background:#f8fafb; border:1px solid #dce5ea; border-radius:6px; padding:6px 8px; min-height:17px; font-weight:400; }
    QAbstractSpinBox:focus, QComboBox:focus, QLineEdit:focus { border-color:#58a5aa; background:#fff; }
    QAbstractSpinBox:disabled, QComboBox:disabled { color:#a5b3bb; background:#f2f5f7; border-color:#e7edf0; }
    QComboBox::drop-down { border:0; width:23px; }
    QComboBox::down-arrow { image:url(:/studio/caret-down.svg); width:10px; height:6px; }
    QAbstractSpinBox::up-button { subcontrol-origin:border; subcontrol-position:top right; width:19px; border-left:1px solid #dce5ea; border-bottom:1px solid #e8eef1; border-top-right-radius:6px; background:#f1f6f8; }
    QAbstractSpinBox::down-button { subcontrol-origin:border; subcontrol-position:bottom right; width:19px; border-left:1px solid #dce5ea; border-bottom-right-radius:6px; background:#f1f6f8; }
    QAbstractSpinBox::up-button:hover, QAbstractSpinBox::down-button:hover { background:#dceeee; }
    QAbstractSpinBox::up-arrow { image:url(:/studio/caret-up.svg); width:8px; height:5px; }
    QAbstractSpinBox::down-arrow { image:url(:/studio/caret-down.svg); width:8px; height:5px; }
    QPushButton { background:#fff; border:1px solid #d8e3e8; border-radius:7px; padding:8px 12px; }
    QPushButton:hover { background:#edf6f6; border-color:#8dbbbd; }
    QPushButton[primary="true"] { color:#fff; background:#218b90; border-color:#218b90; font-weight:600; }
    QTreeWidget { background:#f8fafb; border:0; padding:9px; outline:0; }
    QTreeWidget::item { padding:7px 3px; border-radius:5px; }
    QTreeWidget::item:selected { background:#dceeee; color:#216c73; }
    QHeaderView::section { background:#f8fafb; color:#79919e; border:0; padding:12px 7px; font-weight:600; }
    QTabWidget::pane { border:0; background:#f0f4f6; }
    QTabBar::tab { background:#edf2f5; padding:12px 14px; border:0; color:#78909c; }
    QTabBar::tab:selected { color:#1e8288; background:#f9fbfc; border-bottom:2px solid #218b90; }
    QScrollArea { background:#f0f4f6; border:0; }
    QScrollBar:vertical { background:transparent; width:8px; margin:3px; }
    QScrollBar::handle:vertical { background:#c6d4db; border-radius:3px; min-height:28px; }
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height:0; }
    QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background:transparent; }
    QSplitter::handle { background:#dce5ea; }
    QSplitter::handle:horizontal { width:1px; }
    QSplitter::handle:vertical { height:1px; }
    QPlainTextEdit { background:#fff; border:1px solid #dce5ea; border-radius:7px; padding:9px; }
    QPlainTextEdit#taskLog { background:#263c47; color:#c4d6df; border:0; border-radius:0; padding:7px 14px; font-size:11px; }
    QWidget#viewTools { background:#f8fbfc; border-bottom:1px solid #dce5ea; }
    QStatusBar { background:#f8fbfc; border-top:1px solid #dce5ea; color:#68818e; padding:4px 12px; }
    QProgressBar { background:#e6eff2; border:0; border-radius:4px; max-height:8px; }
    QProgressBar::chunk { background:#69aaaf; border-radius:4px; }
    QSlider::groove:horizontal { height:4px; background:#dce6eb; border-radius:2px; }
    QSlider::handle:horizontal { background:#479b9f; width:12px; height:12px; margin:-4px 0; border-radius:6px; }
  )"));
  auto *file = menuBar()->addMenu(QStringLiteral("文件"));
  file->addAction(QStringLiteral("新建三维工程"), [this] { fresh(3); });
  file->addAction(QStringLiteral("新建二维挤出工程"), [this] { fresh(2); });
  file->addAction(QStringLiteral("新建二维内流工程"),
                  [this] { fresh(2, true); });
  file->addAction(QStringLiteral("新建三维内流工程"),
                  [this] { fresh(3, true); });
  file->addAction(
      QStringLiteral("打开工程…"),
      [this] {
        const auto path = QFileDialog::getOpenFileName(
            this, QStringLiteral("打开 BabelMesh 工程"), {},
            "BabelMesh (*.json)");
        if (!path.isEmpty())
          loadProject(path);
      },
      QKeySequence::Open);
  file->addAction(
      QStringLiteral("保存"), [this] { save(); }, QKeySequence::Save);
  file->addAction(
      QStringLiteral("另存为…"), [this] { save(true); }, QKeySequence::SaveAs);
  file->addSeparator();
  file->addAction(QStringLiteral("导入 STEP…"), [this] { importStep(); });
  file->addAction(QStringLiteral("导出流体域 STEP…"), [this] { exportCad(); });
  auto *meshMenu = menuBar()->addMenu(QStringLiteral("网格"));
  meshMenu->addAction(QStringLiteral("生成混合网格"), [this] { generate(); });
  meshMenu->addAction(QStringLiteral("导入 / 转换网格…"),
                      [this] { convertMesh(); });
  meshMenu->addAction(QStringLiteral("优化已有网格…"),
                      [this] { optimizeMesh(); });
  auto *edit = menuBar()->addMenu(QStringLiteral("编辑"));
  edit->addAction(
      QStringLiteral("撤销"),
      [this] {
        if (!undo_.empty()) {
          redo_.push_back(project_);
          project_ = undo_.back();
          ++revision_;
          undo_.pop_back();
          modified_ = true;
          refresh();
          inspect();
        }
      },
      QKeySequence::Undo);
  edit->addAction(
      QStringLiteral("重做"),
      [this] {
        if (!redo_.empty()) {
          undo_.push_back(project_);
          project_ = redo_.back();
          ++revision_;
          redo_.pop_back();
          modified_ = true;
          refresh();
          inspect();
        }
      },
      QKeySequence::Redo);
  edit->addAction(QStringLiteral("编辑选中操作…"), [this] { editOperation(); });
  edit->addAction(QStringLiteral("删除选中操作"), [this] { removeSelected(); });
  auto *cad = menuBar()->addMenu(QStringLiteral("CAD"));
  for (const auto &item : std::vector<std::pair<QString, std::string>>{
           {QStringLiteral("长方体 / 二维矩形…"), "box"},
           {QStringLiteral("球体…"), "sphere"},
           {QStringLiteral("圆柱 / 二维圆形…"), "cylinder"},
           {QStringLiteral("圆锥…"), "cone"}})
    cad->addAction(item.first, [this, type = item.second] { primitive(type); });
  cad->addAction(QStringLiteral("多边形草图与拉伸…"), [this] { sketch(); });
  cad->addSeparator();
  for (const auto &item : std::vector<std::pair<QString, std::string>>{
           {QStringLiteral("并集…"), "fuse"},
           {QStringLiteral("差集…"), "cut"},
           {QStringLiteral("交集…"), "common"},
           {QStringLiteral("变换…"), "transform"},
           {QStringLiteral("拉伸…"), "extrude"},
           {QStringLiteral("旋转成形…"), "revolve"},
           {QStringLiteral("放样…"), "loft"}})
    cad->addAction(item.first, [this, type = item.second] { operation(type); });
  cad->addAction(QStringLiteral("设置物理边界组…"),
                 [this] { physicalGroup(); });
  auto *toolbar = addToolBar(QStringLiteral("主工具"));
  toolbar->setMovable(false);
  toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  toolbar->setIconSize(QSize(20, 20));
  auto *brand = new QLabel(QStringLiteral("BABELMESH"));
  brand->setStyleSheet(
      "color:#2e5965;font-size:17px;font-weight:700;padding:0 24px 0 2px;");
  toolbar->addWidget(brand);
  toolbar->addAction(toolIcon("import"), QStringLiteral("导入 STEP"),
                     [this] { importStep(); });
  toolbar->addAction(toolIcon("box"), QStringLiteral("长方体"),
                     [this] { primitive("box"); });
  toolbar->addAction(toolIcon("cylinder"), QStringLiteral("圆柱"),
                     [this] { primitive("cylinder"); });
  toolbar->addAction(toolIcon("sphere"), QStringLiteral("球体"),
                     [this] { primitive("sphere"); });
  toolbar->addSeparator();
  toolbar->addAction(toolIcon("refresh"), QStringLiteral("刷新 CAD"),
                     [this] { inspect(); });
  auto *generateAction = toolbar->addAction(
      toolIcon("play"), QStringLiteral("生成混合网格"), [this] { generate(); });
  toolbar->widgetForAction(generateAction)->setProperty("primary", true);
  toolbar->addAction(toolIcon("stop"), QStringLiteral("取消任务"),
                     [this] { cancel(); });
  auto *vertical = new QSplitter(Qt::Vertical, this);
  auto *horizontal = new QSplitter(Qt::Horizontal, vertical);
  tree_ = new QTreeWidget(horizontal);
  tree_->setHeaderLabel(QStringLiteral("工程结构"));
  tree_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  tree_->setMinimumWidth(215);
  auto *viewPanel = new QWidget(horizontal);
  auto *viewLayout = new QVBoxLayout(viewPanel);
  viewLayout->setContentsMargins(0, 0, 0, 0);
  viewLayout->setSpacing(0);
  auto *viewTools = new QWidget(viewPanel);
  viewTools->setObjectName("viewTools");
  auto *vt = new QHBoxLayout(viewTools);
  vt->setContentsMargins(14, 8, 14, 8);
  vt->setSpacing(9);
  view_ = new MeshView(viewPanel);
  viewProjection_ = new QComboBox(viewTools);
  viewProjection_->addItems({QStringLiteral("正交"), QStringLiteral("透视")});
  vt->addWidget(viewProjection_);
  connect(viewProjection_, &QComboBox::currentIndexChanged, this,
          [this](int i) { view_->setPerspective(i == 1); });
  auto *material = new QComboBox(viewTools);
  material->addItems({QStringLiteral("磨砂单元"), QStringLiteral("实体表面")});
  vt->addWidget(material);
  connect(material, &QComboBox::currentIndexChanged, this,
          [this](int i) { view_->setFrosted(i == 0); });
  auto *opacityLabel = new QLabel(QStringLiteral("不透明度"), viewTools);
  vt->addWidget(opacityLabel);
  auto *opacity = new QSlider(Qt::Horizontal, viewTools);
  opacity->setRange(20, 95);
  opacity->setValue(62);
  opacity->setFixedWidth(74);
  opacity->setToolTip(QStringLiteral("磨砂单元不透明度；半透明绘制按深度排序"));
  vt->addWidget(opacity);
  connect(opacity, &QSlider::valueChanged, this,
          [this](int value) { view_->setOpacity(value / 100.0); });
  auto *gaps = new QCheckBox(QStringLiteral("单元间距"), viewTools);
  gaps->setChecked(true);
  gaps->setToolTip(QStringLiteral("仅缩小完整单元的显示外壳；平面剖切保留真实截"
                                  "面。不会改变网格节点与拓扑。"));
  vt->addWidget(gaps);
  connect(gaps, &QCheckBox::toggled, this,
          [this](bool enabled) { view_->setCellInset(enabled ? .035 : 0); });
  vt->addStretch();
  for (const auto &item : std::vector<std::pair<QString, int>>{
           {QStringLiteral("等轴"), -1}, {"XY", 2}, {"XZ", 1}, {"YZ", 0}}) {
    auto *button = new QToolButton(viewTools);
    button->setText(item.first);
    button->setToolTip(QStringLiteral("切换至 %1 视图").arg(item.first));
    vt->addWidget(button);
    connect(button, &QToolButton::clicked, this,
            [this, axis = item.second] { view_->orient(axis); });
  }
  auto *fitView = new QToolButton(viewTools);
  fitView->setIcon(toolIcon("fit"));
  fitView->setToolTip(QStringLiteral("适合窗口"));
  vt->addWidget(fitView);
  connect(fitView, &QToolButton::clicked, view_, &MeshView::fit);
  viewLayout->addWidget(viewTools);
  viewLayout->addWidget(view_, 1);
  auto *tabs = new QTabWidget(horizontal);
  tabs->setMinimumWidth(355);
  auto *scroll = new QScrollArea(tabs);
  scroll->setWidgetResizable(true);
  auto *params = new QWidget;
  auto *layout = new QVBoxLayout(params);
  auto *rules = new QGroupBox(QStringLiteral("生成规则"), params);
  auto *form = new QFormLayout(rules);
  dimension_ = new QComboBox(rules);
  dimension_->addItem(QStringLiteral("三维 · 六面体 + 四棱锥 + 四面体"), 3);
  dimension_->addItem(QStringLiteral("二维 · 一层六面体 + 三棱柱"), 2);
  dimension_->setSizeAdjustPolicy(
      QComboBox::AdjustToMinimumContentsLengthWithIcon);
  dimension_->setMinimumContentsLength(13);
  form->addRow(QStringLiteral("维数规则"), dimension_);
  axis_ = new QComboBox(rules);
  axis_->addItems({"x", "y", "z"});
  form->addRow(QStringLiteral("二维厚度轴"), axis_);
  coordinate_ = real(rules);
  form->addRow(QStringLiteral("截面位置"), coordinate_);
  thickness_ = real(rules, 1e-7);
  form->addRow(QStringLiteral("挤出厚度"), thickness_);
  layout->addWidget(rules);
  auto *grid = new QGroupBox(QStringLiteral("正交笛卡尔背景"), params);
  auto *gf = new QFormLayout(grid);
  division_ = new QComboBox(grid);
  division_->addItem(QStringLiteral("指定各方向分段数"), "counts");
  division_->addItem(QStringLiteral("指定最大背景间距"), "spacing");
  gf->addRow(QStringLiteral("划分方式"), division_);
  gf->addRow(QStringLiteral("最小角"), vectorRow(grid, minimum_));
  gf->addRow(QStringLiteral("最大角"), vectorRow(grid, maximum_));
  auto *countRow = new QWidget(grid);
  auto *cr = new QVBoxLayout(countRow);
  cr->setContentsMargins(0, 0, 0, 0);
  for (int i = 0; i < 3; ++i) {
    counts_[static_cast<std::size_t>(i)] = new QSpinBox(countRow);
    counts_[static_cast<std::size_t>(i)]->setRange(1, 100000);
    counts_[static_cast<std::size_t>(i)]->setPrefix(
        QString("%1 ").arg("xyz"[i]));
    cr->addWidget(counts_[static_cast<std::size_t>(i)]);
  }
  gf->addRow(QStringLiteral("分段数"), countRow);
  gf->addRow(QStringLiteral("最大间距"), vectorRow(grid, spacing_));
  connect(division_, &QComboBox::currentIndexChanged, this, [this](int) {
    const bool space = division_->currentData().toString() == "spacing";
    for (std::size_t a = 0; a < 3; ++a) {
      counts_[a]->setEnabled(!space);
      spacing_[a]->setEnabled(space);
    }
  });

  surface_ = real(grid, 1e-7);
  gf->addRow(QStringLiteral("贴体表面尺寸"), surface_);
  volumeSize_ = real(grid, 1e-7);
  gf->addRow(QStringLiteral("过渡区尺寸"), volumeSize_);
  minimumSize_ = real(grid, 1e-7);
  gf->addRow(QStringLiteral("加密尺寸下限"), minimumSize_);
  maxLevel_ = new QSpinBox(grid);
  maxLevel_->setRange(0, 20);
  gf->addRow(QStringLiteral("背景加密最大层级"), maxLevel_);

  buffer_ = real(grid, 0, 10);
  gf->addRow(QStringLiteral("核心间距系数"), buffer_);
  layout->addWidget(grid);
  auto *refinement = new QGroupBox(QStringLiteral("局部加密区域"), params);
  auto *rl = new QVBoxLayout(refinement);
  sources_ = new QTreeWidget(refinement);
  sources_->setHeaderLabels({"区域 ID", "作用范围", "目标尺寸"});
  sources_->setMaximumHeight(145);
  rl->addWidget(sources_);
  auto *buttons = new QHBoxLayout;
  auto *add = new QPushButton(QStringLiteral("新增区域"), refinement);
  auto *editSource = new QPushButton(QStringLiteral("编辑"), refinement);
  auto *remove = new QPushButton(QStringLiteral("删除"), refinement);
  buttons->addWidget(add);
  buttons->addWidget(editSource);
  buttons->addWidget(remove);
  rl->addLayout(buttons);
  layout->addWidget(refinement);
  connect(add, &QPushButton::clicked, this, [this] { sizeSource(false); });
  connect(editSource, &QPushButton::clicked, this,
          [this] { sizeSource(true); });
  connect(sources_, &QTreeWidget::itemDoubleClicked, this,
          [this](auto *, int) { sizeSource(true); });
  connect(remove, &QPushButton::clicked, this, [this] {
    auto *item = sources_->currentItem();
    if (!item)
      return;
    auto next = project_;
    next["sizing"]["sources"].erase(next["sizing"]["sources"].begin() +
                                    item->data(0, Qt::UserRole).toInt());
    changed(next, false);
  });
  auto *layer = new QGroupBox(QStringLiteral("壁面边界层"), params);
  auto *lf = new QFormLayout(layer);
  layerEnabled_ = new QCheckBox(QStringLiteral("生成指定高度的边界层"), layer);
  lf->addRow(layerEnabled_);
  layerGroup_ = new QComboBox(layer);
  lf->addRow(QStringLiteral("目标壁面组"), layerGroup_);
  connect(layerGroup_, &QComboBox::currentTextChanged, this,
          [this] { loadLayerParameters(); });
  firstLayer_ = real(layer, 1e-7);
  lf->addRow(QStringLiteral("首层高度"), firstLayer_);
  layers_ = new QSpinBox(layer);
  layers_->setRange(1, 256);
  lf->addRow(QStringLiteral("层数"), layers_);
  growth_ = real(layer, 1, 3);
  lf->addRow(QStringLiteral("增长率"), growth_);
  layerSummary_ = new QLabel(layer);
  layerSummary_->setWordWrap(true);
  lf->addRow(layerSummary_);
  auto *applyLayer = new QPushButton(QStringLiteral("保存当前组层参数"), layer);
  auto *removeLayer =
      new QPushButton(QStringLiteral("移除当前组层设置"), layer);
  lf->addRow(applyLayer);
  lf->addRow(removeLayer);
  connect(applyLayer, &QPushButton::clicked, this,
          [this] { applyParameters(); });
  connect(removeLayer, &QPushButton::clicked, this, [this] {
    auto next = project_;
    auto &targets = next["boundary_layers"]["targets"];
    const auto group = layerGroup_->currentText().toStdString();
    targets.erase(
        std::remove_if(targets.begin(), targets.end(),
                       [&](const auto &t) { return t.at("group") == group; }),
        targets.end());
    if (targets.empty())
      next["boundary_layers"]["enabled"] = false;
    changed(next, false);
  });
  connect(firstLayer_, &QDoubleSpinBox::valueChanged, this,
          [this](double) { updateLayerSummary(); });
  connect(growth_, &QDoubleSpinBox::valueChanged, this,
          [this](double) { updateLayerSummary(); });
  connect(layers_, &QSpinBox::valueChanged, this,
          [this](int) { updateLayerSummary(); });

  layout->addWidget(layer);
  auto *output = new QGroupBox(QStringLiteral("质量与输出"), params);
  auto *of = new QFormLayout(output);
  quality_ = new QComboBox(output);
  quality_->addItem(QStringLiteral("诊断 · 保留拓扑合法的网格和质量报告"),
                    "diagnostic");
  quality_->addItem(QStringLiteral("严格 · 面质量失败时停止导出"), "strict");
  quality_->setSizeAdjustPolicy(
      QComboBox::AdjustToMinimumContentsLengthWithIcon);
  quality_->setMinimumContentsLength(12);
  of->addRow(QStringLiteral("质量模式"), quality_);
  cadTolerance_ = real(output, 1e-9);
  cadTolerance_->setDecimals(9);
  of->addRow(QStringLiteral("CAD 距离容差（m）"), cadTolerance_);
  optEnabled_ =
      new QCheckBox(QStringLiteral("生成后自动优化（不增加单元）"), output);
  of->addRow(optEnabled_);
  topologyOpt_ = new QCheckBox(
      QStringLiteral("三维过渡区局部换边（不增加节点或单元）"), output);
  of->addRow(topologyOpt_);
  connect(dimension_, &QComboBox::currentIndexChanged, this, [this](int) {
    topologyOpt_->setEnabled(dimension_->currentData().toInt() == 3);
  });
  mergeOpt_ =
      new QCheckBox(QStringLiteral("尝试合并过渡四面体为凸多面体"), output);
  of->addRow(mergeOpt_);
  optIterations_ = new QSpinBox(output);
  optIterations_->setRange(1, 30);
  of->addRow(QStringLiteral("优化最大迭代次数"), optIterations_);
  overwrite_ = new QCheckBox(QStringLiteral("允许覆盖本工程已有结果"), output);
  of->addRow(overwrite_);
  auto *choose = new QPushButton(QStringLiteral("选择输出目录…"), output);
  of->addRow(choose);
  connect(choose, &QPushButton::clicked, this, [this] {
    const auto path =
        QFileDialog::getExistingDirectory(this, QStringLiteral("输出目录"));
    if (!path.isEmpty()) {
      auto next = project_;
      next["output"]["directory"] = path.toStdString();
      changed(next, false);
    }
  });
  layout->addWidget(output);
  auto *resources = new QGroupBox(QStringLiteral("生成预算"), params);
  auto *rf = new QFormLayout(resources);
  maxCells_ = new QSpinBox(resources);
  maxCells_->setRange(1, 2000000000);
  rf->addRow(QStringLiteral("最大单元数"), maxCells_);
  maxMemory_ = real(resources, 1, 1048576);
  rf->addRow(QStringLiteral("峰值内存预算（MB）"), maxMemory_);
  threads_ = new QSpinBox(resources);
  threads_->setRange(1, 256);
  rf->addRow(QStringLiteral("工作线程数"), threads_);
  layout->addWidget(resources);
  auto *apply = new QPushButton(QStringLiteral("应用参数"), params);
  apply->setProperty("primary", true);
  layout->addWidget(apply);
  connect(apply, &QPushButton::clicked, this, [this] { applyParameters(); });
  layout->addStretch();
  scroll->setWidget(params);
  tabs->addTab(scroll, QStringLiteral("网格参数"));
  auto *advanced = new QWidget(tabs);
  auto *av = new QVBoxLayout(advanced);
  jsonEditor_ = new QPlainTextEdit(advanced);
  jsonEditor_->setFont(QFont("monospace", 10));
  av->addWidget(jsonEditor_);
  auto *applyJson = new QPushButton(QStringLiteral("应用工程 JSON"), advanced);
  av->addWidget(applyJson);
  connect(applyJson, &QPushButton::clicked, this, [this] {
    try {
      changed(Json::parse(jsonEditor_->toPlainText().toStdString()));
    } catch (const std::exception &e) {
      error(QString::fromUtf8(e.what()));
    }
  });
  tabs->addTab(advanced, QStringLiteral("工程 JSON"));
  auto *result = new QWidget(tabs);
  auto *rv = new QVBoxLayout(result);
  summary_ = new QLabel(result);
  summary_->setWordWrap(true);
  summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  rv->addWidget(summary_);
  auto *wire = new QCheckBox(QStringLiteral("显示单元边线"), result);
  wire->setChecked(true);
  rv->addWidget(wire);
  connect(wire, &QCheckBox::toggled, view_, &MeshView::setWireframe);
  rv->addWidget(new QLabel(QStringLiteral("三维剖切"), result));
  clipModeControl_ = new QComboBox(result);
  clipModeControl_->addItems({QStringLiteral("平面剖切 · 平整截面"),
                              QStringLiteral("按整单元筛选 · 原始单元边界")});
  rv->addWidget(clipModeControl_);
  connect(clipModeControl_, &QComboBox::currentIndexChanged, this,
          [this](int index) { view_->setPlanarClip(index == 0); });
  clipAxisControl_ = new QComboBox(result);
  clipAxisControl_->addItems({QStringLiteral("X 法向"),
                              QStringLiteral("Y 法向"),
                              QStringLiteral("Z 法向")});
  clipAxisControl_->setCurrentIndex(2);
  rv->addWidget(clipAxisControl_);
  connect(clipAxisControl_, &QComboBox::currentIndexChanged, view_,
          &MeshView::setClipAxis);
  auto *clip = new QSlider(Qt::Horizontal, result);
  clipControl_ = clip;
  clip->setRange(0, 100);
  clip->setValue(50);
  rv->addWidget(clip);
  connect(clip, &QSlider::valueChanged, this,
          [this](int value) { view_->setClip(value / 100.0); });
  auto *loadMesh = new QPushButton(QStringLiteral("打开网格检查点…"), result);
  rv->addWidget(loadMesh);
  connect(loadMesh, &QPushButton::clicked, this, [this] {
    const auto path = QFileDialog::getOpenFileName(
        this, QStringLiteral("打开 mesh.json"), {}, "Mesh checkpoint (*.json)");
    if (path.isEmpty())
      return;
    Json j;
    if (readJson(path, j)) {
      try {
        showMesh(j);
        currentMesh_ = path;
      } catch (const std::exception &e) {
        error(QString::fromUtf8(e.what()));
      }
    }
  });
  rv->addStretch();
  tabs->addTab(result, QStringLiteral("网格结果"));
  horizontal->setSizes({225, 945, 390});
  log_ = new QPlainTextEdit(vertical);
  log_->setObjectName("taskLog");
  log_->setReadOnly(true);
  log_->setFont(QFont("monospace", 9));
  log_->setMaximumBlockCount(3000);
  vertical->setSizes({830, 105});
  setCentralWidget(vertical);
  progress_ = new QProgressBar;
  progress_->setRange(0, 1);
  progress_->setValue(0);
  progress_->setTextVisible(false);
  progress_->setMaximumWidth(180);
  statusBar()->addPermanentWidget(progress_);
  statusBar()->showMessage(
      QStringLiteral("就绪 · 工程长度单位随配置保存，输出坐标统一为米"));
  connect(tree_, &QTreeWidget::itemSelectionChanged, this, [this] {
    const auto items = tree_->selectedItems();
    if (items.isEmpty())
      return;
    selected_ = items.front()->data(0, Qt::UserRole).toString().toStdString();
    selectedFace_ = items.front()->data(0, Qt::UserRole + 1).toInt();
    view_->selectEntity(selected_);
  });
  connect(tree_, &QTreeWidget::itemDoubleClicked, this,
          [this](QTreeWidgetItem *, int) { editOperation(); });
  view_->picked = [this](const std::string &id, int face) {
    selected_ = id;
    selectedFace_ = face;
    view_->selectEntity(id);
    for (auto *item : tree_->findItems(QString::fromStdString(id),
                                       Qt::MatchContains | Qt::MatchRecursive))
      if (item->data(0, Qt::UserRole).toString().toStdString() == id) {
        tree_->setCurrentItem(item);
        break;
      }
    selectedFace_ = face;
    statusBar()->showMessage(QStringLiteral("已选择 %1 / 面 %2")
                                 .arg(QString::fromStdString(id))
                                 .arg(face));
  };
}
void Studio::showMesh(const Json &mesh) {
  view_->setMesh(mesh);
  if (smokeViewAxis_ >= -1)
    view_->orient(smokeViewAxis_);
  viewProjection_->setCurrentIndex(smokePerspective_ ? 1 : 0);
  view_->setPerspective(smokePerspective_);
  const QSignalBlocker blocker(clipControl_);
  const bool planar = mesh.value("generation_dimension", 3) == 2;
  clipControl_->setValue(planar ? 100 : 50);
  clipControl_->setEnabled(!planar);
  clipModeControl_->setEnabled(!planar);
  clipAxisControl_->setEnabled(!planar);
  clipAxisControl_->setCurrentIndex(planar ? mesh.value("extrusion_axis", 2)
                                           : 2);
}
bool Studio::captureSmoke() {
  const auto framebuffer = view_->grabFramebuffer();
  const auto pixels = framebuffer.convertToFormat(QImage::Format_RGBA8888);
  std::size_t dark = 0, tinted = 0, saturated = 0, nonOpaque = 0;
  for (int y = 0; y < pixels.height(); ++y) {
    const auto *row = pixels.constScanLine(y);
    for (int x = 0; x < pixels.width(); ++x) {
      const auto *pixel = row + 4 * x;
      const int lo = std::min({pixel[0], pixel[1], pixel[2]});
      const int hi = std::max({pixel[0], pixel[1], pixel[2]});
      dark += hi < 64;
      tinted += lo >= 50 && hi - lo > 12;
      saturated += lo < 40 && hi - lo > 180;
      nonOpaque += pixel[3] != 255;
    }
  }
  auto state = view_->displayState();
  const bool mesh = state.at("model") == "mesh";
  const bool coloursPassed = tinted > 1000 && saturated == 0 && nonOpaque == 0;
  const bool edgesPassed =
      !mesh || !state.at("edges_enabled").get<bool>() || dark > 100;
  const bool capsPassed =
      !mesh || state.at("generation_dimension") == 2 ||
      state.at("clip_mode") != "plane" ||
      state.at("clip_fraction").get<double>() >= .999 ||
      (state.at("cap_vertices").get<std::size_t>() > 0 &&
       state.at("max_cap_plane_error_m").get<double>() <= 1e-5);
  state["capture_checks"] = {{"dark_edge_pixels", dark},
                             {"tinted_surface_pixels", tinted},
                             {"saturated_primary_pixels", saturated},
                             {"non_opaque_composite_pixels", nonOpaque},
                             {"colours_passed", coloursPassed},
                             {"edges_passed", edgesPassed},
                             {"planar_caps_passed", capsPassed}};
  const bool ok = view_->renderingReady() && view_->triangleCount() > 0 &&
                  coloursPassed && edgesPassed && capsPassed;
  state["capture_checks"]["passed"] = ok;
  // The direct GL image checks colour/alpha independently of Qt window capture.
  return framebuffer.save(smokeScreenshot_ + ".framebuffer.png") &&
         writeJson(smokeScreenshot_ + ".viewport.json", state) &&
         grab().save(smokeScreenshot_) && ok;
}
bool Studio::readJson(const QString &path, Json &target) {
  QFile file(path);
  if (!file.open(QIODevice::ReadOnly)) {
    error(QStringLiteral("无法读取 %1").arg(path));
    return false;
  }
  try {
    target = Json::parse(file.readAll().toStdString());
    return true;
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
    return false;
  }
}
void Studio::fresh(int dimension, bool internal) {
  if (process_ && process_->state() != QProcess::NotRunning) {
    error(QStringLiteral("请先完成或取消当前任务"));
    return;
  }
  Json next;
  const QString example = internal
                              ? (dimension == 2 ? "channel2d" : "box3d")
                              : (dimension == 2 ? "cylinder2d" : "sphere3d");
  const QString file =
      QString(BABELMESH_SOURCE_DIR) + "/examples/" + example + ".json";
  if (!readJson(file, next))
    return;
  project_ = next;
  project_["output"]["directory"] =
      (QString(BABELMESH_SOURCE_DIR) + "/outputs/" + example).toStdString();
  projectPath_.clear();
  undo_.clear();
  redo_.clear();
  modified_ = false;
  refresh();
  if (isVisible())
    inspect();
}
bool Studio::loadProject(const QString &path) {
  if (process_ && process_->state() != QProcess::NotRunning) {
    error(QStringLiteral("请先完成或取消当前任务"));
    return false;
  }
  Json next;
  if (!readJson(path, next))
    return false;
  try {
    const auto base = QFileInfo(path).absolutePath();
    for (auto &imported : next.at("cad").at("imports")) {
      const auto source =
          QString::fromStdString(imported.at("file").get<std::string>());
      if (QFileInfo(source).isRelative())
        imported["file"] = QDir(base).absoluteFilePath(source).toStdString();
    }
    const auto output = QString::fromStdString(
        next.at("output").at("directory").get<std::string>());
    if (QFileInfo(output).isRelative())
      next["output"]["directory"] =
          QDir(base).absoluteFilePath(output).toStdString();
    if (next.at("schema_version") == 1) {
      next["schema_version"] = 2;
      next["meshing"] = {{"dimension", 3}, {"core_buffer", .25}};
      if (!next.at("sizing").contains("surface_size"))
        next["sizing"]["surface_size"] = next.at("sizing").at("maximum_size");
    }
    project_ = std::move(next);
    projectPath_ = QFileInfo(path).absoluteFilePath();
    undo_.clear();
    redo_.clear();
    modified_ = false;
    refresh();
    inspect();
    return true;
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
    return false;
  }
}
void Studio::refresh() {
  if (project_.is_null())
    return;
  try {
    tree_->clear();
    auto *imports = new QTreeWidgetItem(tree_, {QStringLiteral("STEP 模型")});
    auto *operations =
        new QTreeWidgetItem(tree_, {QStringLiteral("参数化操作")});
    auto *groups = new QTreeWidgetItem(tree_, {QStringLiteral("物理边界组")});
    for (const auto &item : project_.at("cad").at("imports")) {
      auto *node = new QTreeWidgetItem(
          imports, {QString::fromStdString(item.at("id")) + " · STEP"});
      node->setData(0, Qt::UserRole, QString::fromStdString(item.at("id")));
      node->setData(0, Qt::UserRole + 1, -1);
    }
    for (const auto &item : project_.at("cad").at("operations")) {
      auto *node = new QTreeWidgetItem(
          operations,
          {QString::fromStdString(item.at("id").get<std::string>() + " · " +
                                  item.at("operation").get<std::string>())});
      node->setData(0, Qt::UserRole, QString::fromStdString(item.at("id")));
      node->setData(0, Qt::UserRole + 1, -1);
    }
    const auto selectedLayer = layerGroup_->currentText();
    const QSignalBlocker layerSignalBlocker(layerGroup_);
    layerGroup_->clear();
    for (const auto &group : project_.at("cad").at("groups")) {
      new QTreeWidgetItem(
          groups,
          {QString::fromStdString(group.at("name").get<std::string>() + " · " +
                                  group.at("kind").get<std::string>())});
      if (group.at("kind") == "wall")
        layerGroup_->addItem(QString::fromStdString(group.at("name")));
    }
    tree_->expandAll();
    const auto &m = project_.at("meshing");
    const int dim = m.at("dimension");
    dimension_->setCurrentIndex(dimension_->findData(dim));
    const bool planar = dim == 2;
    axis_->setEnabled(planar);
    coordinate_->setEnabled(planar);
    thickness_->setEnabled(planar);
    axis_->setCurrentText(QString::fromStdString(m.value("section_axis", "z")));
    coordinate_->setValue(m.value("section_coordinate", 0.0));
    thickness_->setValue(m.value("thickness", .01));
    const auto normalized = babelmesh::normalizeHybridProject(project_);
    const bool space = project_.at("background").contains("cell_size");
    division_->setCurrentIndex(space ? 1 : 0);
    for (std::size_t a = 0; a < 3; ++a) {
      counts_[a]->setValue(normalized.at("background").at("cell_counts")[a]);
      minimum_[a]->setValue(
          project_.at("background").at("bounds").at("minimum")[a]);
      maximum_[a]->setValue(
          project_.at("background").at("bounds").at("maximum")[a]);
      spacing_[a]->setValue(
          space ? project_.at("background").at("cell_size")[a].get<double>()
                : (maximum_[a]->value() - minimum_[a]->value()) /
                      counts_[a]->value());
      spacing_[a]->setEnabled(space);
      counts_[a]->setEnabled(!space);
    }
    volumeSize_->setValue(project_.at("sizing").value(
        "volume_size", project_.at("sizing").at("maximum_size").get<double>()));
    minimumSize_->setValue(project_.at("sizing").at("minimum_size"));
    maxLevel_->setValue(project_.at("sizing").at("max_level"));
    const auto optimization = project_.value("optimization", Json::object());
    optEnabled_->setChecked(optimization.value("enabled", true));
    topologyOpt_->setChecked(optimization.value("topology_swaps", true));
    topologyOpt_->setEnabled(dim == 3);
    mergeOpt_->setChecked(optimization.value("merge_transition_cells", false));
    mergeOpt_->setEnabled(!planar);
    optIterations_->setValue(optimization.value("iterations", 6));
    sources_->clear();
    for (std::size_t i = 0; i < project_.at("sizing").at("sources").size();
         ++i) {
      const auto &source = project_.at("sizing").at("sources")[i];
      auto *item = new QTreeWidgetItem(
          sources_, {QString::fromStdString(source.at("id")),
                     QString::fromStdString(source.value("apply_to", "all")),
                     QString::number(source.at("target_size").get<double>())});
      item->setData(0, Qt::UserRole, static_cast<int>(i));
    }
    surface_->setValue(project_.at("sizing").at("surface_size"));
    buffer_->setValue(m.value("core_buffer", .35));
    quality_->setCurrentIndex(quality_->findData(
        QString::fromStdString(project_.at("quality").at("mode"))));
    cadTolerance_->setValue(project_.at("quality").value(
        "max_surface_deviation_m", surface_->value() * scale(project_) * .15));
    maxCells_->setValue(project_.at("resources").at("max_cells"));
    maxMemory_->setValue(project_.at("resources").at("max_memory_mb"));
    threads_->setValue(project_.at("resources").at("threads"));
    overwrite_->setChecked(project_.at("output").at("overwrite"));
    layerEnabled_->setChecked(project_.at("boundary_layers").at("enabled"));
    if (layerGroup_->findText(selectedLayer) >= 0) {
      layerGroup_->setCurrentText(selectedLayer);
    } else if (!project_.at("boundary_layers").at("targets").empty()) {
      const auto &t = project_.at("boundary_layers").at("targets")[0];
      layerGroup_->setCurrentText(QString::fromStdString(t.at("group")));
    }
    loadLayerParameters();
    jsonEditor_->setPlainText(QString::fromStdString(project_.dump(2)));
    setWindowTitle(QStringLiteral("BabelMesh · %1%2")
                       .arg(projectPath_.isEmpty()
                                ? QStringLiteral("未保存工程")
                                : QFileInfo(projectPath_).fileName(),
                            modified_ ? " *" : ""));
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
  }
}
bool Studio::changed(Json next, bool preview) {
  if (next == project_)
    return true;
  try {
    babelmesh::validateHybridProject(babelmesh::normalizeHybridProject(next));
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
    return false;
  }
  ++revision_;
  undo_.push_back(project_);
  if (undo_.size() > 100)
    undo_.erase(undo_.begin());
  redo_.clear();
  project_ = std::move(next);
  modified_ = true;
  refresh();
  if (preview)
    inspect();
  return true;
}
bool Studio::applyParameters() {
  try {
    auto next = project_;
    const int dimension = dimension_->currentData().toInt();
    next["meshing"] = {{"dimension", dimension},
                       {"core_buffer", buffer_->value()}};
    if (dimension == 2) {
      const auto axis = axis_->currentText().toStdString();
      const std::size_t a = axis == "x" ? 0 : axis == "y" ? 1 : 2;
      counts_[a]->setValue(1);
      next["meshing"].update({{"section_axis", axis},
                              {"section_coordinate", coordinate_->value()},
                              {"thickness", thickness_->value()}});
    }
    for (std::size_t a = 0; a < 3; ++a) {
      if (division_->currentData().toString() == "spacing")
        next["background"]["cell_size"][a] = spacing_[a]->value();
      else
        next["background"]["cell_counts"][a] = counts_[a]->value();
      next["background"]["bounds"]["minimum"][a] = minimum_[a]->value();
      next["background"]["bounds"]["maximum"][a] = maximum_[a]->value();
    }
    next["background"].erase(division_->currentData().toString() == "spacing"
                                 ? "cell_counts"
                                 : "cell_size");
    next["sizing"]["volume_size"] = volumeSize_->value();
    next["sizing"]["maximum_size"] =
        std::max(next.at("sizing").at("maximum_size").get<double>(),
                 volumeSize_->value());
    next["sizing"]["minimum_size"] = minimumSize_->value();
    next["sizing"]["max_level"] = maxLevel_->value();
    next["optimization"] = {
        {"topology_swaps", topologyOpt_->isChecked()},
        {"enabled", optEnabled_->isChecked()},
        {"iterations", optIterations_->value()},
        {"merge_transition_cells", mergeOpt_->isChecked()},
        {"max_merge_fraction", project_.value("optimization", Json::object())
                                   .value("max_merge_fraction", .1)},
        {"max_cell_growth_fraction", 0.0}};
    if (next.at("fluid_domain").at("mode") == "external")
      next["fluid_domain"]["outer_bounds"] = next.at("background").at("bounds");
    next["sizing"]["surface_size"] = surface_->value();
    next["quality"]["mode"] = quality_->currentData().toString().toStdString();
    next["quality"]["max_surface_deviation_m"] = cadTolerance_->value();
    next["resources"]["max_cells"] = maxCells_->value();
    next["resources"]["max_memory_mb"] = maxMemory_->value();
    next["resources"]["threads"] = threads_->value();
    next["output"]["overwrite"] = overwrite_->isChecked();
    next["boundary_layers"]["enabled"] = layerEnabled_->isChecked();
    if (layerEnabled_->isChecked()) {
      if (layerGroup_->currentText().isEmpty())
        throw std::runtime_error("请先设置壁面物理组");
      const Json target = {{"group", layerGroup_->currentText().toStdString()},
                           {"first_layer", firstLayer_->value()},
                           {"count", layers_->value()},
                           {"growth_ratio", growth_->value()}};
      auto &targets = next["boundary_layers"]["targets"];
      auto existing =
          std::find_if(targets.begin(), targets.end(), [&](const auto &t) {
            return t.at("group") == target.at("group");
          });
      if (existing == targets.end())
        targets.push_back(target);
      else
        *existing = target;
    }
    if (next != project_ && !changed(next, false))
      return false;
    return true;
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
    return false;
  }
}
bool Studio::save(bool saveAs) {
  if (!applyParameters())
    return false;
  QString path = projectPath_;
  if (saveAs || path.isEmpty())
    path = QFileDialog::getSaveFileName(this, QStringLiteral("保存工程"), path,
                                        "BabelMesh (*.json)");
  if (path.isEmpty())
    return false;
  if (!writeJson(path, project_)) {
    error(QStringLiteral("工程保存失败"));
    return false;
  }
  projectPath_ = path;
  modified_ = false;
  refresh();
  return true;
}
void Studio::loadLayerParameters() {
  if (!project_.contains("boundary_layers") || !firstLayer_)
    return;
  const auto name = layerGroup_->currentText().toStdString();
  const auto &targets = project_.at("boundary_layers").at("targets");
  const auto found =
      std::find_if(targets.begin(), targets.end(),
                   [&](const auto &t) { return t.at("group") == name; });
  if (found == targets.end()) {
    firstLayer_->setValue(surface_->value() * .2);
    growth_->setValue(1.2);
    layers_->setValue(3);
  } else {
    firstLayer_->setValue(found->at("first_layer"));
    growth_->setValue(found->at("growth_ratio"));
    layers_->setValue(found->at("count"));
  }
  updateLayerSummary();
}
void Studio::updateLayerSummary() {
  if (!layerSummary_)
    return;
  const double first = firstLayer_->value(), growth = growth_->value();
  double h = first, total = 0;
  for (int i = 0; i < layers_->value(); ++i) {
    total += h;
    h *= growth;
  }
  layerSummary_->setText(
      QStringLiteral(
          "总层厚：%1；末层高度：%2（工程单位）\n层高 = 首层 × 增长率^(层号−1)")
          .arg(total, 0, 'g', 6)
          .arg(first * std::pow(growth, layers_->value() - 1), 0, 'g', 6));
}
void Studio::sizeSource(bool edit) {
  try {
    const int selectedIndex =
        edit && sources_->currentItem()
            ? sources_->currentItem()->data(0, Qt::UserRole).toInt()
            : -1;
    if (edit && selectedIndex < 0)
      return;
    if (!applyParameters())
      return;
    const int index = selectedIndex;
    Json source =
        edit ? project_.at("sizing").at("sources")[index]
             : Json{{"id",
                     "region_" +
                         std::to_string(
                             project_.at("sizing").at("sources").size() + 1)},
                    {"kind", "box"},
                    {"apply_to", "all"},
                    {"target_size", surface_->value() * .5},
                    {"transition_width", surface_->value()},
                    {"region",
                     {{"bounds", project_.at("background").at("bounds")}}}};
    QDialog dialog(this);
    dialog.setWindowTitle(QStringLiteral("局部网格加密"));
    auto *form = new QFormLayout(&dialog);
    auto *id = new QLineEdit(QString::fromStdString(source.at("id")), &dialog);
    form->addRow("ID", id);
    auto *kind = new QComboBox(&dialog);
    kind->addItem(QStringLiteral("长方体区域"), "box");
    kind->addItem(QStringLiteral("球形区域 / 二维圆形区域"), "sphere");
    kind->setCurrentIndex(
        kind->findData(QString::fromStdString(source.at("kind"))));
    if (kind->currentIndex() < 0)
      throw std::runtime_error("这个 CAD 距离区域请在工程 JSON 中编辑");
    form->addRow(QStringLiteral("区域形状"), kind);
    auto *apply = new QComboBox(&dialog);
    apply->addItem(QStringLiteral("背景 + 表面 + 过渡区"), "all");
    apply->addItem(QStringLiteral("笛卡尔背景"), "background");
    apply->addItem(QStringLiteral("CAD 表面"), "surface");
    apply->addItem(QStringLiteral("三角形 / 四面体过渡区"), "transition");
    apply->setCurrentIndex(apply->findData(
        QString::fromStdString(source.value("apply_to", "all"))));
    form->addRow(QStringLiteral("作用范围"), apply);
    auto *target = real(&dialog, 1e-7);
    target->setValue(source.at("target_size"));
    form->addRow(QStringLiteral("目标尺寸"), target);
    auto *width = real(&dialog, 0);
    width->setValue(source.value("transition_width", 0.0));
    form->addRow(QStringLiteral("过渡宽度"), width);
    std::array<QDoubleSpinBox *, 3> low{}, high{}, centre{};
    form->addRow(QStringLiteral("最小角（长方体）"), vectorRow(&dialog, low));
    form->addRow(QStringLiteral("最大角（长方体）"), vectorRow(&dialog, high));
    form->addRow(QStringLiteral("球心（球形）"), vectorRow(&dialog, centre));
    auto *radius = real(&dialog, 1e-7);
    radius->setValue(surface_->value() * 3);
    form->addRow(QStringLiteral("半径（球形）"), radius);
    for (std::size_t a = 0; a < 3; ++a) {
      low[a]->setValue(
          source.at("kind") == "box"
              ? source.at("region").at("bounds").at("minimum")[a].get<double>()
              : minimum_[a]->value());
      high[a]->setValue(
          source.at("kind") == "box"
              ? source.at("region").at("bounds").at("maximum")[a].get<double>()
              : maximum_[a]->value());
      centre[a]->setValue(
          source.at("kind") == "sphere"
              ? source.at("region").at("center")[a].get<double>()
              : (minimum_[a]->value() + maximum_[a]->value()) * .5);
    }
    if (source.at("kind") == "sphere")
      radius->setValue(source.at("region").at("radius"));
    auto *ok = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
    form->addRow(ok);
    connect(ok, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(ok, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted)
      return;
    source["id"] = id->text().toStdString();
    source["kind"] = kind->currentData().toString().toStdString();
    source["apply_to"] = apply->currentData().toString().toStdString();
    source["target_size"] = target->value();
    source["transition_width"] = width->value();
    source["region"] = Json::object();
    if (source.at("kind") == "box") {
      source["region"]["bounds"]["minimum"] = {low[0]->value(), low[1]->value(),
                                               low[2]->value()};
      source["region"]["bounds"]["maximum"] = {
          high[0]->value(), high[1]->value(), high[2]->value()};
    } else {
      source["region"]["center"] = {centre[0]->value(), centre[1]->value(),
                                    centre[2]->value()};
      source["region"]["radius"] = radius->value();
    }
    auto next = project_;
    if (edit)
      next["sizing"]["sources"][index] = source;
    else
      next["sizing"]["sources"].push_back(source);
    changed(next, false);
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
  }
}
QString Studio::uniqueId(const QString &prefix) const {
  std::set<std::string> ids;
  for (const auto &i : project_.at("cad").at("imports"))
    ids.insert(i.at("id").get<std::string>());
  for (const auto &i : project_.at("cad").at("operations"))
    ids.insert(i.at("id").get<std::string>());
  for (int n = 1;; ++n) {
    const auto id = prefix + QString::number(n);
    if (!ids.count(id.toStdString()))
      return id;
  }
}
void Studio::primitive(const std::string &type) {
  QDialog dialog(this);
  dialog.setWindowTitle(
      QStringLiteral("创建参数化几何 · %1").arg(QString::fromStdString(type)));
  auto *layout = new QFormLayout(&dialog);
  std::map<std::string, std::array<QDoubleSpinBox *, 3>> vectors;
  std::map<std::string, QDoubleSpinBox *> scalars;
  const auto id = uniqueId(QString::fromStdString(type) + "_");
  layout->addRow(QStringLiteral("操作 ID"), new QLabel(id, &dialog));
  auto addVector = [&](const std::string &key, const QString &name,
                       const std::array<double, 3> &value) {
    layout->addRow(name, vectorRow(&dialog, vectors[key]));
    for (std::size_t i = 0; i < 3; ++i)
      vectors[key][i]->setValue(value[i]);
  };
  auto addReal = [&](const std::string &key, const QString &name, double value,
                     double minimum) {
    auto *input = real(&dialog, minimum);
    input->setValue(value);
    scalars[key] = input;
    layout->addRow(name, input);
  };
  if (type == "box") {
    addVector("minimum", QStringLiteral("最小角"), {-.3, -.3, -.3});
    addVector("maximum", QStringLiteral("最大角"), {.3, .3, .3});
  } else if (type == "sphere") {
    addVector("center", QStringLiteral("球心"), {0, 0, 0});
    addReal("radius", QStringLiteral("半径"), .3, 1e-7);
  } else {
    addVector("axis_origin", QStringLiteral("轴起点"), {0, 0, -.5});
    addVector("direction", QStringLiteral("轴方向"), {0, 0, 1});
    addReal("height", QStringLiteral("高度"), 1, 1e-7);
    if (type == "cylinder")
      addReal("radius", QStringLiteral("半径"), .3, 1e-7);
    else {
      addReal("bottom_radius", QStringLiteral("底部半径"), .3, 0);
      addReal("top_radius", QStringLiteral("顶部半径"), .1, 0);
    }
  }
  auto *buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  layout->addRow(buttons);
  connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if (dialog.exec() != QDialog::Accepted)
    return;
  Json params = Json::object();
  for (const auto &item : vectors)
    params[item.first] = {item.second[0]->value(), item.second[1]->value(),
                          item.second[2]->value()};
  for (const auto &item : scalars)
    params[item.first] = item.second->value();
  auto next = project_;
  next["cad"]["operations"].push_back({{"id", id.toStdString()},
                                       {"operation", type},
                                       {"inputs", Json::array()},
                                       {"parameters", params}});
  const auto key = next.at("fluid_domain").at("mode") == "external"
                       ? "solids"
                       : "fluid_solids";
  next["fluid_domain"][key].push_back(id.toStdString());
  next["cad"]["groups"].push_back(
      {{"name", id.toStdString() + "_wall"},
       {"kind", "wall"},
       {"selector",
        {{"operation_id", id.toStdString()}, {"entity_type", "face"}}}});
  changed(next);
}
void Studio::sketch() {
  QDialog dialog(this);
  dialog.setWindowTitle(QStringLiteral("多边形草图与拉伸"));
  auto *layout = new QVBoxLayout(&dialog);
  layout->addWidget(new QLabel(
      QStringLiteral("每行输入平面坐标 u,v。闭合多边形会沿厚度轴拉伸成实体。"),
      &dialog));
  auto *points =
      new QPlainTextEdit("-0.4,-0.3\n0.4,-0.3\n0.2,0.4\n-0.2,0.4", &dialog);
  layout->addWidget(points);
  auto *buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  layout->addWidget(buttons);
  connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if (dialog.exec() != QDialog::Accepted)
    return;
  try {
    const auto axis = project_.at("meshing").value("section_axis", "z");
    const std::size_t a = axis == "x"   ? 0
                          : axis == "y" ? 1
                                        : 2,
                      u = (a + 1) % 3, v = (a + 2) % 3;
    std::vector<std::array<double, 3>> coords;
    for (const auto &line :
         points->toPlainText().split('\n', Qt::SkipEmptyParts)) {
      const auto values = line.split(',');
      if (values.size() != 2)
        throw std::runtime_error("每行需要两个以逗号分隔的坐标");
      bool first = false, second = false;
      std::array<double, 3> p{};
      p[a] = project_.at("background").at("bounds").at("minimum")[a];
      p[u] = values[0].trimmed().toDouble(&first);
      p[v] = values[1].trimmed().toDouble(&second);
      if (!first || !second)
        throw std::runtime_error("草图坐标必须是数字");
      coords.push_back(p);
    }
    if (coords.size() < 3)
      throw std::runtime_error("多边形至少需要三个点");
    const auto base = uniqueId("sketch_").toStdString();
    auto next = project_;
    Json edges = Json::array();
    for (std::size_t i = 0; i < coords.size(); ++i) {
      const auto id = base + "_e" + std::to_string(i);
      edges.push_back(id);
      next["cad"]["operations"].push_back(
          {{"id", id},
           {"operation", "line"},
           {"inputs", Json::array()},
           {"parameters",
            {{"start", coords[i]}, {"end", coords[(i + 1) % coords.size()]}}}});
    }
    next["cad"]["operations"].push_back({{"id", base + "_wire"},
                                         {"operation", "wire"},
                                         {"inputs", edges},
                                         {"parameters", {{"closed", true}}}});
    next["cad"]["operations"].push_back({{"id", base + "_face"},
                                         {"operation", "plane"},
                                         {"inputs", {base + "_wire"}},
                                         {"parameters", Json::object()}});
    std::array<double, 3> direction{};
    direction[a] =
        project_.at("background").at("bounds").at("maximum")[a].get<double>() -
        coords[0][a];
    next["cad"]["operations"].push_back(
        {{"id", base},
         {"operation", "extrude"},
         {"inputs", {base + "_face"}},
         {"parameters", {{"direction", direction}}}});
    const auto key = next.at("fluid_domain").at("mode") == "external"
                         ? "solids"
                         : "fluid_solids";
    next["fluid_domain"][key].push_back(base);
    next["cad"]["groups"].push_back(
        {{"name", base + "_wall"},
         {"kind", "wall"},
         {"selector", {{"operation_id", base}, {"entity_type", "face"}}}});
    changed(next);
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
  }
}
void Studio::operation(const std::string &type) {
  QStringList ids;
  for (const auto &item : project_.at("cad").at("imports"))
    ids << QString::fromStdString(item.at("id"));
  for (const auto &item : project_.at("cad").at("operations"))
    ids << QString::fromStdString(item.at("id"));
  if (ids.isEmpty())
    return;
  QDialog dialog(this);
  dialog.setWindowTitle(
      QStringLiteral("CAD 操作 · %1").arg(QString::fromStdString(type)));
  auto *layout = new QFormLayout(&dialog);
  auto *first = new QComboBox(&dialog);
  first->addItems(ids);
  first->setCurrentText(QString::fromStdString(selected_));
  layout->addRow(QStringLiteral("第一输入"), first);
  const bool boolean = type == "cut" || type == "fuse" || type == "common";
  auto *second = new QComboBox(&dialog);
  second->addItems(ids);
  second->setCurrentIndex(ids.size() > 1 ? 1 : 0);
  if (boolean)
    layout->addRow(QStringLiteral("第二输入"), second);
  auto *params = new QPlainTextEdit(&dialog);
  Json initial = Json::object();
  if (type == "transform")
    initial = {
        {"translation", {0.1, 0, 0}}, {"rotation", {0, 0, 0}}, {"scale", 1}};
  else if (type == "extrude")
    initial = {{"direction", {0, 0, 1}}};
  else if (type == "revolve")
    initial = {{"axis_origin", {0, 0, 0}},
               {"direction", {0, 0, 1}},
               {"angle_radians", 6.283185307179586}};
  else if (type == "loft")
    initial = {{"solid", true}};
  params->setPlainText(QString::fromStdString(initial.dump(2)));
  layout->addRow(QStringLiteral("参数（工程单位；转角为弧度）"), params);
  auto *inputs = new QPlainTextEdit(&dialog);
  if (type == "loft") {
    inputs->setPlainText(ids.join(','));
    layout->addRow(QStringLiteral("放样线框 ID，逗号分隔"), inputs);
  }
  auto *buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  layout->addRow(buttons);
  connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if (dialog.exec() != QDialog::Accepted)
    return;
  try {
    Json input = {first->currentText().toStdString()};
    if (boolean)
      input.push_back(second->currentText().toStdString());
    if (type == "loft") {
      input = Json::array();
      for (const auto &id :
           inputs->toPlainText().split(',', Qt::SkipEmptyParts))
        input.push_back(id.trimmed().toStdString());
    }
    const auto id = uniqueId(QString::fromStdString(type) + "_").toStdString();
    auto next = project_;
    next["cad"]["operations"].push_back(
        {{"id", id},
         {"operation", type},
         {"inputs", input},
         {"parameters", Json::parse(params->toPlainText().toStdString())}});
    const auto key = next.at("fluid_domain").at("mode") == "external"
                         ? "solids"
                         : "fluid_solids";
    Json bodies = Json::array();
    for (const auto &body : next["fluid_domain"][key])
      if (std::find(input.begin(), input.end(), body) == input.end())
        bodies.push_back(body);
    bodies.push_back(id);
    next["fluid_domain"][key] = bodies;
    Json groups = Json::array();
    for (const auto &group : next["cad"]["groups"]) {
      const auto old = group.at("selector").value("operation_id", "");
      if (std::find(input.begin(), input.end(), Json(old)) == input.end())
        groups.push_back(group);
    }
    groups.push_back(
        {{"name", id + "_wall"},
         {"kind", "wall"},
         {"selector", {{"operation_id", id}, {"entity_type", "face"}}}});
    next["cad"]["groups"] = groups;
    changed(next);
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
  }
}
void Studio::editOperation() {
  if (selected_.empty())
    return;
  for (std::size_t i = 0; i < project_.at("cad").at("operations").size(); ++i)
    if (project_["cad"]["operations"][i].at("id") == selected_) {
      bool accepted = false;
      const auto text = QInputDialog::getMultiLineText(
          this,
          QStringLiteral("编辑操作 · %1")
              .arg(QString::fromStdString(selected_)),
          QStringLiteral("修改参数与输入，ID 保持不变"),
          QString::fromStdString(project_["cad"]["operations"][i].dump(2)),
          &accepted);
      if (!accepted)
        return;
      try {
        auto next = project_;
        const auto item = Json::parse(text.toStdString());
        if (item.at("id") != selected_)
          throw std::runtime_error("修改 ID 会影响依赖，请保留原 ID");
        next["cad"]["operations"][i] = item;
        changed(next);
      } catch (const std::exception &e) {
        error(QString::fromUtf8(e.what()));
      }
      return;
    }
}
void Studio::importStep() {
  const auto file =
      QFileDialog::getOpenFileName(this, QStringLiteral("导入 STEP"), {},
                                   "STEP (*.step *.stp *.STEP *.STP)");
  if (file.isEmpty())
    return;
  bool ok = false;
  const auto unit =
      QInputDialog::getItem(this, QStringLiteral("STEP 单位"),
                            QStringLiteral("选择文件声明的长度单位"),
                            {"mm", "m", "cm", "inch"}, 0, false, &ok);
  if (!ok)
    return;
  const auto id = uniqueId("step_").toStdString();
  auto next = project_;
  next["cad"]["imports"].push_back({{"id", id},
                                    {"file", file.toStdString()},
                                    {"length_unit", unit.toStdString()}});
  const auto key = next.at("fluid_domain").at("mode") == "external"
                       ? "solids"
                       : "fluid_solids";
  next["fluid_domain"][key].push_back(id);
  next["cad"]["groups"].push_back(
      {{"name", id + "_wall"},
       {"kind", "wall"},
       {"selector", {{"operation_id", id}, {"entity_type", "face"}}}});
  changed(next);
}
void Studio::physicalGroup() {
  if (selected_.empty())
    return;
  bool ok = false;
  const auto name = QInputDialog::getText(
      this, QStringLiteral("物理边界组"), QStringLiteral("组名（无空格）"),
      QLineEdit::Normal,
      QString::fromStdString(selected_) +
          (selectedFace_ >= 0 ? QString("_face_%1").arg(selectedFace_)
                              : "_wall"),
      &ok);
  if (!ok || name.isEmpty())
    return;
  const auto kind = QInputDialog::getItem(
      this, QStringLiteral("边界类型"), QStringLiteral("类型"),
      {"wall", "inlet", "outlet", "symmetry", "generic"}, 0, false, &ok);
  if (!ok)
    return;
  Json selector = {{"operation_id", selected_}, {"entity_type", "face"}};
  if (selectedFace_ >= 0) {
    for (const auto &entity : preview_.at("entities"))
      if (entity.at("id") == selected_) {
        const auto &centre = entity.at("face_centres_m")
                                 .at(static_cast<std::size_t>(selectedFace_));
        selector["center"] = {centre[0].get<double>() / scale(project_),
                              centre[1].get<double>() / scale(project_),
                              centre[2].get<double>() / scale(project_)};
        selector["tolerance"] = 1e-6 / scale(project_);
      }
  }
  auto next = project_;
  Json groups = Json::array();
  for (auto g : next["cad"]["groups"]) {
    if (g.at("name") == name.toStdString())
      continue;
    if (g.at("selector").value("operation_id", "") == selected_) {
      if (selectedFace_ < 0)
        continue;
      if (g.at("selector").contains("center") &&
          g.at("selector").at("center") == selector.at("center"))
        continue;
      g["selector"]["exclude_centers"].push_back(selector.at("center"));
    }
    groups.push_back(std::move(g));
  }
  groups.push_back({{"name", name.toStdString()},
                    {"kind", kind.toStdString()},
                    {"selector", selector}});
  next["cad"]["groups"] = groups;
  changed(next);
}
void Studio::removeSelected() {
  if (selected_.empty())
    return;
  for (const auto &item : project_.at("cad").at("operations"))
    for (const auto &input : item.at("inputs"))
      if (input == selected_) {
        error(QStringLiteral("该操作被 %1 引用，请先删除依赖操作")
                  .arg(QString::fromStdString(item.at("id"))));
        return;
      }
  auto next = project_;
  for (const char *key : {"operations", "imports"}) {
    Json keep = Json::array();
    for (const auto &i : next["cad"][key])
      if (i.at("id") != selected_)
        keep.push_back(i);
    next["cad"][key] = keep;
  }
  Json groups = Json::array();
  for (const auto &g : next["cad"]["groups"])
    if (g.at("selector").value("operation_id", "") != selected_)
      groups.push_back(g);
  next["cad"]["groups"] = groups;
  const auto key = next.at("fluid_domain").at("mode") == "external"
                       ? "solids"
                       : "fluid_solids";
  Json bodies = Json::array();
  for (const auto &b : next["fluid_domain"][key])
    if (b != selected_)
      bodies.push_back(b);
  next["fluid_domain"][key] = bodies;
  selected_.clear();
  changed(next);
}
Studio::Json Studio::jobProject() const {
  auto p = project_;
  for (auto &item : p["cad"]["imports"]) {
    const auto path = QString::fromStdString(item.at("file"));
    if (QFileInfo(path).isRelative())
      item["file"] =
          QDir(projectPath_.isEmpty() ? QString(BABELMESH_SOURCE_DIR)
                                      : QFileInfo(projectPath_).absolutePath())
              .absoluteFilePath(path)
              .toStdString();
  }
  return p;
}
void Studio::inspect() { startJob("inspect-cad"); }
void Studio::generate() {
  if (applyParameters())
    startJob("generate");
}
void Studio::convertMesh() {
  if (process_ && process_->state() != QProcess::NotRunning) {
    error(QStringLiteral("请先完成或取消当前任务"));
    return;
  }
  const auto input = QFileDialog::getOpenFileName(
      this, QStringLiteral("导入 / 转换体网格"), currentMesh_,
      "Volume meshes (*.bsm *.msh *.cas *.cgns *.json);;All files (*)");
  if (input.isEmpty())
    return;
  bool ok = false;
  const QStringList labels = {
      QStringLiteral("BabelSim v3 (.bsm)"),
      QStringLiteral("Fluent ASCII (.msh)"), QStringLiteral("CGNS (.cgns)"),
      QStringLiteral("CGNS NGON/NFACE (.cgns)"), QStringLiteral("VTU (.vtu)")};
  const QStringList formats = {"babelsim-v3", "fluent-ascii", "cgns",
                               "cgns-poly", "vtu"};
  const auto choice =
      QInputDialog::getItem(this, QStringLiteral("输出格式"),
                            QStringLiteral("转换目标"), labels, 0, false, &ok);
  if (!ok)
    return;
  const int index = labels.indexOf(choice);
  const QStringList suffixes = {"bsm", "msh", "cgns", "cgns", "vtu"};
  const auto output = QFileDialog::getSaveFileName(
      this, QStringLiteral("保存转换网格"),
      QFileInfo(input).absolutePath() + "/converted." + suffixes[index]);
  if (output.isEmpty())
    return;
  jobDirectory_ = std::make_unique<QTemporaryDir>();
  if (!jobDirectory_->isValid()) {
    error(QStringLiteral("无法创建转换任务目录"));
    return;
  }
  jobOutput_ = output;
  QStringList arguments = {
      "convert", input,  "--to",         formats[index],
      "--out",   output, "--checkpoint", jobDirectory_->path() + "/mesh.json"};
  if (overwrite_->isChecked())
    arguments.append("--overwrite");
  launchProcess("convert", arguments);
}
void Studio::optimizeMesh() {
  if (process_ && process_->state() != QProcess::NotRunning) {
    error(QStringLiteral("请先完成或取消当前任务"));
    return;
  }
  const auto input = QFileDialog::getOpenFileName(
      this, QStringLiteral("选择要优化的体网格"), currentMesh_,
      "Volume meshes (*.bsm *.msh *.cgns *.json)");
  if (input.isEmpty())
    return;
  const auto output = QFileDialog::getExistingDirectory(
      this, QStringLiteral("选择优化结果目录"));
  if (output.isEmpty())
    return;
  jobOutput_ = output;
  QStringList arguments = {
      "optimize", input,          "--out",
      output,     "--iterations", QString::number(optIterations_->value())};
  if (overwrite_->isChecked())
    arguments.append("--overwrite");
  launchProcess("optimize", arguments);
}
void Studio::exportCad() {
  if (process_ && process_->state() != QProcess::NotRunning) {
    error(QStringLiteral("请先完成或取消当前任务"));
    return;
  }
  if (!applyParameters())
    return;
  const auto output = QFileDialog::getSaveFileName(
      this, QStringLiteral("导出流体域 STEP"), {}, "STEP (*.step *.stp)");
  if (output.isEmpty())
    return;
  jobDirectory_ = std::make_unique<QTemporaryDir>();
  if (!jobDirectory_->isValid()) {
    error(QStringLiteral("无法创建 CAD 导出目录"));
    return;
  }
  const auto file = jobDirectory_->path() + "/project.json";
  if (!writeJson(file, jobProject())) {
    error(QStringLiteral("无法保存导出任务的输入"));
    return;
  }
  jobOutput_ = output;
  QStringList arguments = {"export-cad", file, "--out", output};
  if (overwrite_->isChecked())
    arguments.append("--overwrite");
  launchProcess("export-cad", arguments);
}
void Studio::startJob(const QString &kind) {
  if (process_ && process_->state() != QProcess::NotRunning) {
    statusBar()->showMessage(
        QStringLiteral("当前任务仍在运行；修改已保留，请完成后刷新"));
    return;
  }
  if (!QFileInfo::exists(cli_)) {
    error(QStringLiteral("找不到核心程序 %1").arg(cli_));
    return;
  }
  jobDirectory_ = std::make_unique<QTemporaryDir>(QDir::tempPath() +
                                                  "/babelmesh-gui-XXXXXX");
  if (!jobDirectory_->isValid()) {
    error(QStringLiteral("无法创建任务目录"));
    return;
  }
  auto snapshot = jobProject();
  jobRevision_ = revision_;
  const auto file = jobDirectory_->path() + "/project.json";
  if (kind == "inspect-cad") {
    snapshot["output"]["overwrite"] = true;
    jobOutput_ = jobDirectory_->path() + "/cad";
  } else {
    auto output = QString::fromStdString(snapshot.at("output").at("directory"));
    if (QFileInfo(output).isRelative())
      output =
          QDir(projectPath_.isEmpty() ? QString(BABELMESH_SOURCE_DIR)
                                      : QFileInfo(projectPath_).absolutePath())
              .absoluteFilePath(output);
    jobOutput_ = output;
  }
  if (!writeJson(file, snapshot)) {
    error(QStringLiteral("无法保存任务输入"));
    return;
  }
  launchProcess(kind, {kind, file, "--out", jobOutput_});
}
void Studio::launchProcess(const QString &kind, const QStringList &arguments) {
  if (process_ && process_->state() != QProcess::NotRunning) {
    statusBar()->showMessage(QStringLiteral("请等待当前任务完成或取消任务"));
    return;
  }
  if (!QFileInfo::exists(cli_)) {
    error(QStringLiteral("找不到核心程序 %1").arg(cli_));
    return;
  }
  jobRevision_ = revision_;
  if (process_)
    process_->deleteLater();
  process_ = new QProcess(this);
  jobKind_ = kind;
  cancelled_ = false;
  ++serial_;
  process_->setProcessChannelMode(QProcess::MergedChannels);
  process_->setWorkingDirectory(QDir(BABELMESH_SOURCE_DIR).absolutePath());
  connect(process_, &QProcess::readyReadStandardOutput, this, [this] {
    const auto text = QString::fromUtf8(process_->readAllStandardOutput());
    log_->moveCursor(QTextCursor::End);
    log_->insertPlainText(text);
    if (cancellationTest_ && jobKind_ == "generate" && text.contains("[1/7]"))
      cancel();
  });
  connect(process_, qOverload<int, QProcess::ExitStatus>(&QProcess::finished),
          this, [this](int code, QProcess::ExitStatus status) {
            finishJob(code, status);
          });
  connect(process_, &QProcess::errorOccurred, this,
          [this](QProcess::ProcessError e) {
            if (e == QProcess::FailedToStart) {
              progress_->setRange(0, 1);
              error(process_->errorString());
            }
          });
  log_->appendPlainText(QStringLiteral("\n%1 → %2").arg(kind, jobOutput_));
  progress_->setRange(0, 0);
  statusBar()->showMessage(
      QStringLiteral("正在执行 %1；界面可继续操作").arg(kind));
  process_->start(cli_, arguments);
}
void Studio::finishJob(int exitCode, QProcess::ExitStatus status) {
  progress_->setRange(0, 1);
  progress_->setValue(1);
  if (cancelled_) {
    statusBar()->showMessage(QStringLiteral("任务已取消"));
    if (cancellationTest_) {
      const bool unpublished = !QFileInfo::exists(jobOutput_ + "/mesh.json");
      const bool stopped = process_->state() == QProcess::NotRunning;
      const bool saved = writeJson(QFileInfo(jobOutput_).absolutePath() +
                                       "/cancellation-report.json",
                                   {{"cancelled", true},
                                    {"process_stopped", stopped},
                                    {"result_unpublished", unpublished},
                                    {"exit_code", exitCode}});
      QApplication::exit(unpublished && stopped && saved ? 0 : 1);
    }
    return;
  }
  if (exitCode || status != QProcess::NormalExit) {
    error(QStringLiteral(
              "%1 失败。详细原因见任务日志。已有结果未通过这次任务重新验收。")
              .arg(jobKind_));
    return;
  }
  try {
    if (jobKind_ == "export-cad") {
      statusBar()->showMessage(
          QStringLiteral("STEP 已导出：%1").arg(jobOutput_));
      return;
    }
    if (jobKind_ == "inspect-cad") {
      if (!readJson(jobOutput_ + "/cad_preview.json", preview_))
        return;
      view_->setCad(preview_);
      statusBar()->showMessage(
          jobRevision_ == revision_
              ? QStringLiteral("CAD 预览已更新")
              : QStringLiteral(
                    "此视图来自任务启动时的工程；请刷新 CAD 查看修改"));
      if (workflowMode_ || cancellationTest_) {
        surface_->setValue(surface_->value() * .98);
        if (!save())
          return;
        generate();
      } else if (smokeMode_) {
        if (!smokeMesh_.isEmpty()) {
          Json mesh;
          if (!readJson(smokeMesh_, mesh))
            return;
          showMesh(mesh);
        }
        QTimer::singleShot(800, this, [this] {
          const bool ok = captureSmoke();
          QApplication::exit(ok ? 0 : 1);
        });
      }
    } else {
      Json mesh, report;
      const auto meshFile = jobKind_ == "convert"
                                ? jobDirectory_->path() + "/mesh.json"
                                : jobOutput_ + "/mesh.json";
      const auto reportFile = jobKind_ == "convert"
                                  ? jobOutput_ + ".report.json"
                                  : jobOutput_ + "/mesh_report.json";
      if (!readJson(meshFile, mesh) || !readJson(reportFile, report))
        return;
      showMesh(mesh);
      currentMesh_ = jobKind_ == "convert" ? jobOutput_ : meshFile;
      const auto &q = report.at("quality");
      summary_->setText(
          QStringLiteral(
              "生成规则：%1D；输出：3D "
              "体网格\n\n单元数：%2\n六面体：%3\n三棱柱：%4\n四面体：%"
              "5\n四棱锥：%6\n多面体：%12\n\n最大非正交角：%7°\nP99 "
              "非正交角：%8°\n最大内部偏斜：%9\nP99 "
              "内部偏斜：%13\n面质量门槛：%10\n完整 CFD "
              "profile：尚未全部评估\n\n输出目录：\n%11")
              .arg(report.at("generation_dimension").get<int>())
              .arg(report.at("cells").get<qulonglong>())
              .arg(report.at("cell_types").at("hexahedron").get<qulonglong>())
              .arg(report.at("cell_types").at("prism").get<qulonglong>())
              .arg(report.at("cell_types").at("tetrahedron").get<qulonglong>())
              .arg(report.at("cell_types").at("pyramid").get<qulonglong>())
              .arg(q.at("max_non_orthogonality_degrees").get<double>(), 0, 'f',
                   2)
              .arg(q.at("p99_non_orthogonality_degrees").get<double>(), 0, 'f',
                   2)
              .arg(q.at("max_internal_skewness").get<double>(), 0, 'f', 3)
              .arg(q.at("face_quality_thresholds_passed").get<bool>()
                       ? QStringLiteral("通过")
                       : QStringLiteral("未通过；诊断结果"))
              .arg(jobOutput_)
              .arg(report.at("cell_types").at("polyhedron").get<qulonglong>())
              .arg(q.at("p99_internal_skewness").get<double>(), 0, 'f', 3));
      statusBar()->showMessage(
          jobRevision_ == revision_
              ? QStringLiteral("%1 完成 · %2")
                    .arg(jobKind_, QString::fromStdString(report.value(
                                       "status", "topology-validated")))
              : QStringLiteral(
                    "任务完成；结果来自启动时的参数，工程修改尚未重新生成"));
      if (workflowMode_)
        QTimer::singleShot(800, this, [this] {
          const bool ok = QFileInfo::exists(projectPath_) && captureSmoke();
          QApplication::exit(ok ? 0 : 1);
        });
    }
  } catch (const std::exception &e) {
    error(QString::fromUtf8(e.what()));
  }
}
void Studio::cancel() {
  if (!process_ || process_->state() == QProcess::NotRunning)
    return;
  cancelled_ = true;
  process_->terminate();
  const auto serial = serial_;
  QTimer::singleShot(1500, this, [this, serial] {
    if (serial_ == serial && process_ &&
        process_->state() != QProcess::NotRunning)
      process_->kill();
  });
}
void Studio::error(const QString &message) {
  log_->appendPlainText(QStringLiteral("错误：") + message);
  statusBar()->showMessage(message);
  if (smokeMode_) {
    qWarning("%s", qPrintable(message));
    QApplication::exit(1);
  } else
    QMessageBox::warning(this, QStringLiteral("BabelMesh"), message);
}
void Studio::smoke(const QString &project, const QString &mesh,
                   const QString &screenshot) {
  smokeMode_ = true;
  smokeMesh_ = mesh;
  smokeScreenshot_ = screenshot;
  loadProject(project);
}
void Studio::workflowSmoke(const QString &project, const QString &directory,
                           const QString &screenshot) {
  smokeMode_ = true;
  workflowMode_ = true;
  smokeScreenshot_ = screenshot;
  if (!QDir().mkpath(directory)) {
    error(QStringLiteral("无法创建界面验证目录"));
    return;
  }
  if (!loadProject(project))
    return;
  // The real save/parameter/job path is exercised after CAD preview completes.
  projectPath_ = QDir(directory).absoluteFilePath("edited-project.json");
  project_["output"]["directory"] =
      QDir(directory).absoluteFilePath("generated").toStdString();
  project_["output"]["overwrite"] = true;
  refresh();
}
void Studio::cancellationSmoke(const QString &project,
                               const QString &directory) {
  smokeMode_ = true;
  cancellationTest_ = true;
  if (!QDir().mkpath(directory) || !loadProject(project))
    return;
  projectPath_ = QDir(directory).absoluteFilePath("edited-project.json");
  project_["output"]["directory"] =
      QDir(directory).absoluteFilePath("generated").toStdString();
  project_["output"]["overwrite"] = false;
  refresh();
}
void Studio::closeEvent(QCloseEvent *e) {
  if (modified_ && !smokeMode_) {
    const auto answer = QMessageBox::question(
        this, QStringLiteral("未保存的工程"),
        QStringLiteral("保存当前工程修改？"),
        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    if (answer == QMessageBox::Cancel ||
        (answer == QMessageBox::Save && !save())) {
      e->ignore();
      return;
    }
  }
  if (process_ && process_->state() != QProcess::NotRunning) {
    cancel();
    process_->waitForFinished(1600);
  }
  e->accept();
}
