#pragma once

#include <nlohmann/json.hpp>

#include <array>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace babelmesh {

class CadModel;

struct MeshFace {
    std::vector<std::int32_t> vertices;
    std::int32_t owner{};
    std::int32_t neighbour{-1};
    std::int32_t patch{-1};
};

struct HexMesh {
    bool adaptive{};
    bool cad_clipping_applied{};
    bool cad_surface_refinement_applied{};
    bool cad_snapping_applied{};
    bool cad_surface_deviation_overridden{};
    bool cad_bidirectional_surface_audit{};
    double cad_surface_target_size_m{};
    double cad_mesh_to_surface_max_m{};
    double cad_surface_to_mesh_max_m{};
    double cad_surface_audit_sample_spacing_m{};
    int cad_surface_buffer_layers{};
    double cad_snap_fraction{};
    double cad_surface_deviation_m{};
    double cad_surface_deviation_tolerance_m{};
    std::size_t cad_removed_cells{};
    std::size_t cad_intersect_cells{};
    std::size_t cad_uncertain_cells{};
    std::size_t cad_seed_forced_cells{};
    std::size_t cad_selected_components{};
    std::size_t cad_snapping_vertices{};
    std::size_t cad_projection_faces{};
    std::size_t cad_surface_deviation_samples{};
    std::size_t cad_mesh_to_surface_samples{};
    std::size_t cad_surface_to_mesh_samples{};
    std::size_t cad_collapsed_faces_removed{};
    std::size_t cad_collapsed_face_edge_repairs{};
    std::size_t cad_consolidated_faces{};
    std::size_t cad_quality_rollback_count{};
    std::size_t cad_deformed_vertices{};
    std::size_t cad_deformation_iterations{};
    std::size_t cad_tangential_smoothing_iterations{};
    std::size_t cad_quality_optimization_attempts{};
    std::size_t cad_quality_optimization_iterations{};
    std::size_t cad_quality_optimization_proposed_vertices{};
    double cad_tangential_smoothing_max_move_m{};
    double cad_quality_optimization_objective_before{};
    double cad_quality_optimization_objective_after{};
    double cad_max_deformation_m{};
    double cad_core_deformation_fraction{};
    double cad_fluid_domain_volume_m3{};
    double cad_mesh_volume_relative_error{};
    std::array<std::int32_t, 3> cell_counts{};
    std::array<double, 6> bounds_m{};
    std::vector<std::array<double, 3>> points_m;
    std::vector<std::array<std::int32_t, 8>> cells;
    std::vector<std::int32_t> cell_levels;
    std::vector<bool> strict_hex;
    std::vector<double> cell_volumes_m3;
    std::vector<double> cell_aspect_ratios;
    std::map<std::string, std::size_t> source_intersection_counts;
    std::map<std::string, std::vector<std::uint8_t>> sizing_source_markers;
    std::map<std::string, nlohmann::json> sizing_source_details;
    std::vector<MeshFace> faces;
    std::vector<std::pair<std::string, std::string>> patches;
};

HexMesh makeUniformCartesianMesh(const nlohmann::json& project);
HexMesh makeAdaptiveCartesianMesh(const nlohmann::json& project, const CadModel* cadModel = nullptr,
                                  bool refineFluidBoundary = false);
nlohmann::json meshReport(const HexMesh& mesh, const std::filesystem::path& projectFile);
void writeBabelSimV3(const HexMesh& mesh, const std::filesystem::path& file);
void writeLegacyVtu(const HexMesh& mesh, const std::filesystem::path& file);
int runGenerate(int argc, char** argv);

}  // namespace babelmesh
