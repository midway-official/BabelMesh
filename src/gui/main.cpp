#include "studio.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFont>
#include <QSurfaceFormat>
#include <QTimer>

int main(int argc, char **argv) {
  QSurfaceFormat format;
  format.setDepthBufferSize(24);
  format.setStencilBufferSize(8);
  format.setSamples(4);
  format.setVersion(2, 1);
  QSurfaceFormat::setDefaultFormat(format);
  QApplication app(argc, argv);
  app.setStyle("Fusion");
  app.setFont(QFont("Noto Sans CJK SC", 10));
  app.setApplicationName("BabelMesh");
  app.setOrganizationName("BabelMesh");
  QCommandLineParser parser;
  parser.setApplicationDescription("Qt CAD and Cartesian hybrid meshing");
  parser.addHelpOption();
  parser.addOption({"cli", "Path to the BabelMesh command-line core", "path"});
  parser.addOption(
      {"perspective", "Use perspective for viewport smoke capture"});
  parser.addOption({"view-axis",
                    "View axis for smoke capture: auto, x, y, z, iso", "axis",
                    "auto"});
  parser.addOption({"smoke-test",
                    "Open a project and verify the rendered viewport",
                    "project"});
  parser.addOption(
      {"mesh", "Optional mesh checkpoint for the viewport smoke test", "path"});
  parser.addOption({"screenshot", "Save the smoke-test window image", "path"});
  parser.addOption({"workflow-test",
                    "Exercise project edit/save and GUI meshing", "project"});
  parser.addOption(
      {"output", "Output directory for the GUI workflow test", "directory"});
  parser.addOption(
      {"cancel-test", "Verify cancellation of a real meshing job", "project"});
  parser.addPositionalArgument("project", "BabelMesh project JSON",
                               "[project]");
  parser.process(app);
  const QString cli =
      parser.isSet("cli")
          ? parser.value("cli")
          : QDir(QCoreApplication::applicationDirPath()).filePath("babelmesh");
  Studio window(cli);
  window.setSmokePerspective(parser.isSet("perspective"));
  const auto axis = parser.value("view-axis");
  if (axis != "auto" && axis != "iso" && axis != "x" && axis != "y" &&
      axis != "z")
    parser.showHelp(2);
  window.setSmokeViewAxis(axis == "iso" ? -1
                          : axis == "x" ? 0
                          : axis == "y" ? 1
                          : axis == "z" ? 2
                                        : -2);
  window.show();
  if (parser.isSet("cancel-test")) {
    if (!parser.isSet("output"))
      return 2;
    QTimer::singleShot(0, &window, [&] {
      window.cancellationSmoke(parser.value("cancel-test"),
                               parser.value("output"));
    });
    QTimer::singleShot(60000, &app, [] { QApplication::exit(3); });
  } else if (parser.isSet("workflow-test")) {
    if (!parser.isSet("screenshot") || !parser.isSet("output"))
      return 2;
    QTimer::singleShot(0, &window, [&] {
      window.workflowSmoke(parser.value("workflow-test"),
                           parser.value("output"), parser.value("screenshot"));
    });
    QTimer::singleShot(240000, &app, [] { QApplication::exit(3); });
  } else if (parser.isSet("smoke-test")) {
    if (!parser.isSet("screenshot"))
      return 2;
    QTimer::singleShot(0, &window, [&] {
      window.smoke(parser.value("smoke-test"), parser.value("mesh"),
                   parser.value("screenshot"));
    });
    QTimer::singleShot(60000, &app, [] { QApplication::exit(3); });
  } else if (!parser.positionalArguments().isEmpty())
    QTimer::singleShot(0, &window, [&] {
      window.loadProject(parser.positionalArguments().front());
    });
  else
    QTimer::singleShot(0, &window, [&] { window.previewDefault(); });
  return app.exec();
}
