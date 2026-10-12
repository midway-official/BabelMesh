#include "babelmesh/version.h"

#ifndef BABELMESH_VERSION
#define BABELMESH_VERSION "0.1.0"
#endif

namespace babelmesh {

std::string_view version() noexcept {
    return BABELMESH_VERSION;
}

}  // namespace babelmesh
