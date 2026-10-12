#pragma once
#include "babelmesh/mixed_mesh.hpp"
namespace babelmesh {
void reconstructFaceCells(MixedMesh &mesh);
MixedMesh readMixedBabelSimV3(const std::filesystem::path &file);
MixedMesh readFluentAscii(const std::filesystem::path &file);
void writeFluentAscii(const MixedMesh &mesh, const std::filesystem::path &file);
MixedMesh readCgns(const std::filesystem::path &file);
void writeCgns(const MixedMesh &mesh, const std::filesystem::path &file,
               bool forcePolyhedral = false);
MixedMesh readExchangeMesh(const std::filesystem::path &file,
                           const std::string &format = "auto");
void writeExchangeMesh(const MixedMesh &mesh, const std::filesystem::path &file,
                       const std::string &format);
int runConvertMesh(int argc, char **argv);
int runOptimizeMesh(int argc, char **argv);
int runExportCad(int argc, char **argv);
} // namespace babelmesh
