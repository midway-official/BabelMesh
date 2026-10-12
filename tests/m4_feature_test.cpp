#include "babelmesh/cad.hpp"
#include "babelmesh/cartesian_mesh.hpp"
#include "babelmesh/polyhedral_audit.hpp"
#include "babelsim/mesh_io.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: m4_feature_test PROJECT.json");
        nlohmann::json base;
        std::ifstream(argv[1]) >> base;
        base["sizing"]["surface_size"] = 0.0625;
        base["sizing"]["max_level"] = 2;
        for (const double degrees : {15.0, 30.0, 45.0}) {
            auto project = base;
            const double angle = degrees * std::acos(-1.0) / 180.0;
            auto& transform = project["cad"]["operations"][1]["parameters"];
            transform["rotation"] = {0.0, 0.0, angle};
            transform["translation"] = {1.0-std::cos(angle)+std::sin(angle), 1.0-std::sin(angle)-std::cos(angle), 0.0};
            const auto prefix = std::filesystem::temp_directory_path() / ("babelmesh-feature-"+std::to_string(static_cast<int>(degrees)));
            const auto jsonFile = prefix.string()+".json";
            const auto meshFile = prefix.string()+".bsm";
            std::ofstream(jsonFile) << project.dump(2) << '\n';
            const auto cad = babelmesh::buildCadModel(jsonFile);
            const auto mesh = babelmesh::makeAdaptiveCartesianMesh(project, &cad, true);
            babelmesh::requirePolyhedralShells(mesh);
            if (!mesh.cad_snapping_applied || mesh.cad_surface_deviation_m > mesh.cad_surface_deviation_tolerance_m)
                throw std::runtime_error("feature capture failed locked CAD deviation gate");
            const auto report = babelmesh::meshReport(mesh,jsonFile);
            if (report["quality"]["min_scaled_jacobian"].is_null() ||
                report["quality"]["max_non_orthogonality_degrees"].is_null() ||
                report["quality"]["max_internal_skewness"].is_null() ||
                report["quality"]["full_cfd_profile_evaluated"].get<bool>())
                throw std::runtime_error("mesh quality metrics must be computed while unsupported profile gates stay explicit");
            if (report["quality"]["min_scaled_jacobian"].get<double>() < 0.1)
                throw std::runtime_error("strict hex scaled Jacobian is below the locked 0.1 threshold: " +
                    report["quality"]["min_scaled_jacobian"].dump());
            for(const auto& [key,limit]:std::vector<std::pair<std::string,double>>{
                {"max_non_orthogonality_degrees",75},
                {"p99_non_orthogonality_degrees_face_count",60},
                {"max_internal_skewness",.5},
                {"p99_internal_skewness_face_count",.25}})
                if(report["quality"][key].get<double>()>limit)
                    throw std::runtime_error("angle="+std::to_string(degrees)+" failed "+key+": "+report["quality"][key].dump());
            babelmesh::writeBabelSimV3(mesh,meshFile);
            babelsim::readMeshFile(meshFile).validate();
            std::filesystem::remove(jsonFile);
            std::filesystem::remove(meshFile);
            std::cout << "Planar feature regression passed: angle=" << degrees << " cells=" << mesh.cells.size()
                      << " deviation=" << mesh.cad_surface_deviation_m << " volume_error=" << mesh.cad_mesh_volume_relative_error << '\n';
            std::cout << "quality max_nonorth=" << report["quality"]["max_non_orthogonality_degrees"]
                      << " p99_nonorth=" << report["quality"]["p99_non_orthogonality_degrees_face_count"]
                      << " max_skew=" << report["quality"]["max_internal_skewness"]
                      << " p99_skew=" << report["quality"]["p99_internal_skewness_face_count"]
                      << " min_jacobian=" << report["quality"]["min_scaled_jacobian"] << '\n';
        }
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
