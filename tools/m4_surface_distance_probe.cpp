#include "babelmesh/surface_distance.hpp"
#include "babelsim/mesh_io.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    try {
        if(argc!=5)throw std::runtime_error("usage: m4_surface_distance_probe project.json mesh.bsm spacing_m tolerance_m");
        const double spacing=std::stod(argv[3]),tolerance=std::stod(argv[4]);
        if(!(tolerance>0)||!std::isfinite(tolerance))throw std::runtime_error("tolerance must be finite and positive");
        const auto source=babelsim::readMeshFile(argv[2]);source.validate();
        babelmesh::HexMesh mesh;
        for(int v=0;v<source.vertexCount();++v){const auto p=source.vertex(v);mesh.points_m.push_back({p.x,p.y,p.z});}
        for(int f=0;f<source.faceCount();++f)if(source.neighbour(f)<0){
            const auto ring=source.faceVertices(f);
            mesh.faces.push_back({{ring.begin(),ring.end()},source.owner(f),-1,source.boundaryPatch(f)});
        }
        const auto cad=babelmesh::buildCadModel(argv[1]);
        const auto a=babelmesh::auditSurfaceDistance(mesh,cad,spacing);
        const bool passed=std::max(a.mesh_to_cad_max_m,a.cad_to_mesh_max_m)<=tolerance;
        std::cout<<nlohmann::json{{"project",argv[1]},{"mesh",argv[2]},
            {"mesh_to_cad_max_m",a.mesh_to_cad_max_m},{"cad_to_mesh_max_m",a.cad_to_mesh_max_m},
            {"mesh_to_cad_samples",a.mesh_to_cad_samples},{"cad_to_mesh_samples",a.cad_to_mesh_samples},
            {"sample_spacing_on_source_triangles_m",a.requested_sample_spacing_m},
            {"cad_reference_tessellation_deflection_m",a.cad_tessellation_deflection_m},
            {"cad_samples_projected_to_exact_brep",true},{"continuous_hausdorff_bound",false},
            {"tolerance_m",tolerance},{"sampled_bidirectional_passed",passed}}.dump(2)<<'\n';
        return passed?0:1;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
