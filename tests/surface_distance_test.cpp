#include "babelmesh/surface_distance.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char** argv){
    try {
        if(argc!=2)throw std::runtime_error("usage: surface_distance_test repository_root");
        std::ifstream input(std::filesystem::path(argv[1])/"acceptance/projects/cartesian-box.json");
        nlohmann::json project;input>>project;
        project["cad"]["operations"]=nlohmann::json::array({
            {{"id","box"},{"operation","box"},{"inputs",nlohmann::json::array()},
             {"parameters",{{"minimum",{0,0,0}},{"maximum",{1,1,1}}}}}});
        project["cad"]["groups"]=nlohmann::json::array();
        project["fluid_domain"]={{"mode","internal"},{"fluid_solids",{"box"}},
            {"seeds",{{0.5,0.5,0.5}}},{"caps",nlohmann::json::array()}};
        project["background"]["bounds"]={{"minimum",{0,0,0}},{"maximum",{1,1,1}}};
        project["background"]["cell_counts"]={1,1,1};
        project["sizing"]["maximum_size"]=1.0;
        const auto path=std::filesystem::temp_directory_path()/"babelmesh-surface-distance-test.json";
        {std::ofstream output(path);output<<project.dump();}
        const auto cad=babelmesh::buildCadModel(path);
        std::filesystem::remove(path);
        const auto mesh=babelmesh::makeUniformCartesianMesh(project);
        auto report=babelmesh::auditSurfaceDistance(mesh,cad,0.25);
        require(report.mesh_to_cad_max_m<1e-12&&report.cad_to_mesh_max_m<1e-12,
            "coincident box surfaces must have zero sampled distance in both directions");
        auto shifted=mesh;for(auto& p:shifted.points_m)p[0]+=0.002;
        report=babelmesh::auditSurfaceDistance(shifted,cad,0.25);
        require(std::abs(report.mesh_to_cad_max_m-0.002)<1e-12&&std::abs(report.cad_to_mesh_max_m-0.002)<1e-12,
            "translated box must recover the analytic 2 mm distance in both directions");
        auto missing=mesh;missing.faces.erase(missing.faces.begin());
        report=babelmesh::auditSurfaceDistance(missing,cad,0.25);
        require(report.mesh_to_cad_max_m<1e-12&&report.cad_to_mesh_max_m>0.2,
            "reverse distance must detect a missing boundary that one-way samples miss");
        bool rejected=false;
        try{(void)babelmesh::auditSurfaceDistance(mesh,cad,0.25,1);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"sample budget exhaustion must fail rather than report partial success");
        auto degenerate=mesh;degenerate.faces[0].vertices[1]=degenerate.faces[0].vertices[0];
        rejected=false;
        try{(void)babelmesh::auditSurfaceDistance(degenerate,cad,0.25);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"degenerate boundary triangles must be rejected");
        std::cout<<"Bidirectional surface distance: exact box, 2 mm translation, missing boundary and invalid input passed\n";
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
