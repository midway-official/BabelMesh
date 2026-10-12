#pragma once

#include "babelmesh/cad.hpp"
#include "babelmesh/mixed_mesh.hpp"

namespace babelmesh {

void validateHybridProject(const nlohmann::json &project);
nlohmann::json normalizeHybridProject(nlohmann::json project);
MixedMesh generateHybridMesh(const nlohmann::json &project,
                             const CadModel &cad);
int runGenerateHybrid(int argc, char **argv);
int runCheckMixed(int argc, char **argv);

} // namespace babelmesh
