#include "babelmesh/cad.hpp"
#include "babelmesh/cartesian_mesh.hpp"
#include "babelmesh/hybrid_mesher.hpp"
#include "babelmesh/mesh_conversion.hpp"
#include "babelmesh/version.h"

#include <fstream>
#include <iostream>
#include <string_view>

namespace {

void printHelp(std::ostream &output) {
  output
      << "BabelMesh " << babelmesh::version() << "\n"
      << "Cartesian hybrid volume mesh generator (2D extrusion / 3D)\n\n"
      << "Usage:\n"
      << "  babelmesh --help\n"
      << "  babelmesh --version\n"
      << "  babelmesh <command> [arguments]\n\n"
      << "Commands:\n"
      << "  inspect-cad   inspect and report a CAD project (M1)\n"
      << "  generate      v2: Cartesian hybrid mesh; v1: legacy Cartesian/M4 "
         "path\n"
      << "  check         validate a mixed-mesh JSON checkpoint\n"
      << "  convert       convert Fluent ASCII, CGNS, BabelSim v3 and "
         "checkpoints\n"
      << "  export        alias for convert\n"
      << "  optimize      optimize transition vertices without adding cells\n"
      << "  export-cad    export an entity or the fluid domain to STEP\n"
      << "  verify-suite  reserved M8 command; use make test-hybrid today\n";
}

int runCommand(std::string_view command, int argc, char **argv) {
  if (command == "inspect-cad") {
    return babelmesh::runInspectCad(argc, argv);
  } else if (command == "generate") {
    if (argc >= 2) {
      try {
        std::ifstream input(argv[1]);
        nlohmann::json project;
        input >> project;
        if (project.value("schema_version", 1) == 2)
          return babelmesh::runGenerateHybrid(argc, argv);
      } catch (const std::exception &) { /* Legacy entry supplies its
                                            established input diagnostics. */
      }
    }
    return babelmesh::runGenerate(argc, argv);
  } else if (command == "check") {
    return babelmesh::runCheckMixed(argc, argv);
  } else if (command == "export") {
    return babelmesh::runConvertMesh(argc, argv);
  } else if (command == "convert") {
    return babelmesh::runConvertMesh(argc, argv);
  } else if (command == "optimize") {
    return babelmesh::runOptimizeMesh(argc, argv);
  } else if (command == "export-cad") {
    return babelmesh::runExportCad(argc, argv);
  } else if (command == "verify-suite") {
    std::cerr << "verify-suite is not implemented yet (planned for M8)\n";
  } else {
    std::cerr << "unknown command: " << command << "\n";
    return 2;
  }
  return 3;
}

} // namespace

int main(int argc, char **argv) {
  if (argc == 1 || (argc == 2 && std::string_view(argv[1]) == "--help")) {
    printHelp(std::cout);
    return 0;
  }
  if (argc == 2 && std::string_view(argv[1]) == "--version") {
    std::cout << babelmesh::version() << '\n';
    return 0;
  }
  if (argc >= 2 && argv[1][0] == '-') {
    std::cerr << "unknown option: " << argv[1] << "\n";
    return 2;
  }
  return runCommand(argv[1], argc - 1, argv + 1);
}
