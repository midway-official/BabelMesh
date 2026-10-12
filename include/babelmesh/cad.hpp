#pragma once

#include <nlohmann/json.hpp>
#include <TopoDS_Shape.hxx>

#include <array>
#include <cstddef>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace babelmesh {

struct SurfaceTriangle {
    std::array<std::size_t, 3> vertices{};
    std::size_t face_index{};
};

struct SurfaceBvhNode {
    std::array<double, 6> bounds{};
    std::size_t left{};
    std::size_t right{};
    std::size_t begin{};
    std::size_t count{};
};

struct SurfaceMesh {
    double linear_deflection_m{};
    std::size_t degenerate_triangles_skipped{};
    std::vector<std::array<double, 3>> vertices;
    std::vector<SurfaceTriangle> triangles;
    std::vector<SurfaceBvhNode> bvh_nodes;
    std::vector<std::size_t> bvh_triangle_order;
};
void buildSurfaceIndex(SurfaceMesh& surface);
bool pointInsideSurface(const SurfaceMesh& surface,const std::array<double,3>& point);

struct CadEntity {
    std::string id;
    std::string operation;
    TopoDS_Shape shape;
    std::vector<std::string> inputs;
    std::map<std::string, std::vector<TopoDS_Shape>> named_subshapes;
    std::string source_file;
    std::vector<std::string> source_length_units;
    SurfaceMesh surface_mesh;
};

struct PhysicalGroup {
    std::string name;
    std::string kind;
    std::vector<TopoDS_Shape> selected_subshapes;
};

class CadModel {
public:
    const std::map<std::string, CadEntity>& entities() const noexcept { return entities_; }
    const std::vector<PhysicalGroup>& groups() const noexcept { return groups_; }
    const TopoDS_Shape& fluidDomain() const noexcept { return fluid_domain_; }
    const SurfaceMesh& fluidSurfaceMesh() const noexcept { return fluid_surface_mesh_; }
    const std::string& domainMode() const noexcept { return domain_mode_; }

private:
    friend CadModel buildCadModel(const std::filesystem::path& project_file);
    friend CadModel makeFluidCadModel(const TopoDS_Shape& shape);
    std::map<std::string, CadEntity> entities_;
    std::vector<PhysicalGroup> groups_;
    TopoDS_Shape fluid_domain_;
    SurfaceMesh fluid_surface_mesh_;
    std::string domain_mode_;
};

CadModel buildCadModel(const std::filesystem::path& project_file);
CadModel makeFluidCadModel(const TopoDS_Shape& shape);
CadModel inspectCadProject(const std::filesystem::path& project_file);
nlohmann::json geometryReport(
    const CadModel& model,
    const std::filesystem::path& project_file);
nlohmann::json cadPreview(const CadModel& model);

int runInspectCad(int argc, char** argv);

}  // namespace babelmesh
