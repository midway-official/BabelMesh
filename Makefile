# BabelMesh build and verification entry points.
# Usage: make [all|debug|release|test|deps|clean] [BUILD_ROOT=build] [-j4]

SHELL := /bin/bash
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

CXX ?= c++
AR ?= ar
PYTHON ?= python3
BUILD_ROOT ?= build
MODE ?= debug
BABELSIM_ROOT ?= ../BabelSim

BUILD_DIR := $(BUILD_ROOT)/$(MODE)
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
INCLUDE_DIRS := -Iinclude
CPPFLAGS += $(INCLUDE_DIRS)
COMMON_CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -Wshadow -Wconversion
ifeq ($(MODE),release)
CXXFLAGS += $(COMMON_CXXFLAGS) -O3 -DNDEBUG
else ifeq ($(MODE),debug)
CXXFLAGS += $(COMMON_CXXFLAGS) -O0 -g3 -fno-omit-frame-pointer
else
$(error MODE must be debug or release, got '$(MODE)')
endif

# Ubuntu 22.04 / OCCT 7.5.1 library set. Override OCCT_INCLUDE_DIR,
# OCCT_LIBRARY_DIR, or OCCT_LIBS when using another supported installation.
OCCT_INCLUDE_DIR ?= /usr/include/opencascade
OCCT_LIBRARY_DIR ?= /usr/lib/x86_64-linux-gnu
OCCT_LIBS ?= -lTKXDESTEP -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase \
             -lTKXSBase -lTKXCAF -lTKCAF -lTKLCAF -lTKCDF -lTKV3d \
             -lTKService -lTKOffset -lTKBO -lTKBool -lTKShHealing -lTKTopAlgo \
             -lTKGeomAlgo -lTKGeomBase -lTKMesh -lTKPrim -lTKBRep \
             -lTKG3d -lTKG2d -lTKMath -lTKernel

CPPFLAGS += -I$(OCCT_INCLUDE_DIR) -I/usr/include/eigen3
LDLIBS += -L$(OCCT_LIBRARY_DIR) $(OCCT_LIBS) -pthread
GMSH_INCLUDE_DIR ?= $(firstword $(wildcard /usr/include/gmsh.h .cache/dev-packages/gmsh/usr/include/gmsh.h))
GMSH_INCLUDE_DIR := $(patsubst %/gmsh.h,%,$(GMSH_INCLUDE_DIR))
GMSH_LIBRARY ?= $(firstword $(wildcard /usr/lib/x86_64-linux-gnu/libgmsh.so /usr/lib/x86_64-linux-gnu/libgmsh.so.4.8))
CPPFLAGS += -I$(GMSH_INCLUDE_DIR)
LDLIBS += $(GMSH_LIBRARY)

CORE_SOURCES := $(sort $(wildcard src/core/*.cpp) $(wildcard src/cad/*.cpp))
CLI_SOURCES := $(sort $(wildcard src/cli/*.cpp))
CORE_OBJECTS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(CORE_SOURCES))
CLI_OBJECTS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(CLI_SOURCES))
CORE_LIBRARY := $(LIB_DIR)/libbabelmesh.a
CLI_BINARY := $(BIN_DIR)/babelmesh
DEPENDENCY_PROBE := $(BIN_DIR)/babelmesh-dependency-probe
M1_CAD_TEST := $(BIN_DIR)/m1_cad_test
M2_CARTESIAN_TEST := $(BIN_DIR)/m2_cartesian_test
M3_ADAPTIVE_TEST := $(BIN_DIR)/m3_adaptive_test
M4_DOMAIN_TEST := $(BIN_DIR)/m4_domain_test
POLYHEDRAL_AUDIT_TEST := $(BIN_DIR)/polyhedral_audit_test
M4_FEATURE_TEST := $(BIN_DIR)/m4_feature_test
BABELSIM_READER_SOURCES := $(BABELSIM_ROOT)/src/core/mesh.cpp $(BABELSIM_ROOT)/src/io/mesh_reader.cpp
BABELSIM_READER_HEADERS := $(wildcard $(BABELSIM_ROOT)/include/babelsim/*.h)

.PHONY: all debug release gui test-hybrid deps test test-polyhedral-audit test-m4-features test-surface-distance clean help

all: $(CLI_BINARY)

debug:
	$(MAKE) MODE=debug all

release:
	$(MAKE) MODE=release all

gui:
	cmake -S . -B $(BUILD_ROOT)/qt -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD_ROOT)/qt --target babelmesh-gui -j4

test-hybrid:
	cmake -S . -B $(BUILD_ROOT)/qt -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD_ROOT)/qt -j4
	ctest --test-dir $(BUILD_ROOT)/qt --output-on-failure
	$(PYTHON) tests/hybrid_integration_test.py $(BUILD_ROOT)/qt/bin/babelmesh $(BUILD_ROOT)/qt/bin/babelsim-mesh-check .

$(CORE_LIBRARY): $(CORE_OBJECTS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

$(CLI_BINARY): $(CORE_LIBRARY) $(CLI_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $(CLI_OBJECTS) $(CORE_LIBRARY) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(DEPENDENCY_PROBE): tools/dependency_probe.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

$(M1_CAD_TEST): tests/m1_cad_test.cpp $(CORE_LIBRARY)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

$(M2_CARTESIAN_TEST): tests/m2_cartesian_test.cpp $(CORE_LIBRARY) $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

$(M3_ADAPTIVE_TEST): tests/m3_adaptive_test.cpp $(CORE_LIBRARY) $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

$(M4_DOMAIN_TEST): tests/m4_domain_test.cpp $(CORE_LIBRARY) $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

$(POLYHEDRAL_AUDIT_TEST): tests/polyhedral_audit_test.cpp $(CORE_LIBRARY)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

test-polyhedral-audit: $(POLYHEDRAL_AUDIT_TEST)
	$(POLYHEDRAL_AUDIT_TEST)

$(M4_FEATURE_TEST): tests/m4_feature_test.cpp $(CORE_LIBRARY) $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

test-m4-features: $(M4_FEATURE_TEST)
	$(M4_FEATURE_TEST) acceptance/evidence/m4-rotated/project.json

$(BIN_DIR)/m4_solver_quality_probe: tools/m4_solver_quality_probe.cpp $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) -o $@

$(BIN_DIR)/m4_surface_distance_probe: tools/m4_surface_distance_probe.cpp $(CORE_LIBRARY) $(BABELSIM_READER_SOURCES) $(BABELSIM_READER_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -I$(BABELSIM_ROOT)/include $(CXXFLAGS) $< $(BABELSIM_READER_SOURCES) $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

$(BIN_DIR)/surface_distance_test: tests/surface_distance_test.cpp $(CORE_LIBRARY)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(CORE_LIBRARY) $(LDFLAGS) $(LDLIBS) -o $@

test-surface-distance: $(BIN_DIR)/surface_distance_test
	$(BIN_DIR)/surface_distance_test .

deps: $(DEPENDENCY_PROBE)
	$(PYTHON) tools/check_dependencies.py dependencies/ubuntu-jammy.lock.json
	$(DEPENDENCY_PROBE)

test: all deps $(M1_CAD_TEST) $(M2_CARTESIAN_TEST) $(M3_ADAPTIVE_TEST) $(M4_DOMAIN_TEST) test-polyhedral-audit test-m4-features test-surface-distance
	$(PYTHON) -c 'import jsonschema'
	$(PYTHON) tests/m0_contract_test.py $(CLI_BINARY) .
	$(M1_CAD_TEST) .
	$(M2_CARTESIAN_TEST)
	$(M3_ADAPTIVE_TEST)
	$(M4_DOMAIN_TEST)

clean:
	rm -rf $(BUILD_ROOT)

help:
	@printf '%s\n' \
	  'all      Build the default Debug CLI and core library' \
	  'debug    Build Debug artifacts' \
	  'release  Build optimized Release artifacts' \
	  'deps     Compile and run an OCCT/JSON/Eigen dependency probe' \
	  'test     Run M0-M4 CAD, Cartesian, adaptive topology and fluid-domain acceptance checks' \
	  'clean    Remove only $(BUILD_ROOT)'

-include $(CORE_OBJECTS:.o=.d) $(CLI_OBJECTS:.o=.d)

CGNS_INCLUDE_DIR := $(firstword $(wildcard /usr/include/cgnslib.h .cache/dev-packages/cgns/usr/include/cgnslib.h))
CGNS_LIBRARY := $(firstword $(wildcard /usr/lib/x86_64-linux-gnu/libcgns.so .cache/dev-packages/cgns/usr/lib/x86_64-linux-gnu/libcgns.so))
CPPFLAGS += -I$(dir $(CGNS_INCLUDE_DIR))
LDLIBS += $(CGNS_LIBRARY) -Wl,-rpath,$(abspath $(dir $(CGNS_LIBRARY)))
