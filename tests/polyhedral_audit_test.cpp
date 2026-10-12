#include "babelmesh/polyhedral_audit.hpp"
#include "babelmesh/polyhedral_reconstruction.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
void require(bool ok,const std::string& message) { if(!ok) throw std::runtime_error(message); }

babelmesh::HexMesh cube() {
    babelmesh::HexMesh m;
    m.points_m={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    m.cells.push_back({0,1,2,3,4,5,6,7});
    m.strict_hex={true};
    for(auto ring:std::vector<std::vector<std::int32_t>>{{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}})
        m.faces.push_back({ring,0,-1,0});
    m.patches={{"wall","wall"}};
    return m;
}

int main(int argc,char** argv) {
    try {
        if(argc==3 || argc==4) {
            // Replay a generated v3 file; this diagnostic parser is not a
            // replacement for the authoritative BabelSim reader.
            std::ifstream input(argv[1]);
            std::string label; int version=0; std::size_t n=0;
            input>>label>>version;
            require(label=="BABELSIM_MESH" && version==3,"expected v3 file");
            babelmesh::HexMesh m;
            input>>label>>n; require(label=="vertices","missing vertices");
            m.points_m.resize(n);
            for(auto& p:m.points_m) input>>p[0]>>p[1]>>p[2];
            input>>label>>n; require(label=="faces","missing faces");
            m.faces.resize(n);
            std::int32_t maxCell=-1;
            for(auto& f:m.faces) {
                input>>label>>n>>f.owner>>f.neighbour>>f.patch;
                require(label=="face","missing face");
                f.vertices.resize(n); for(auto& v:f.vertices) input>>v;
                maxCell=std::max(maxCell,std::max(f.owner,f.neighbour));
            }
            input>>label>>n; require(label=="patches","missing patches");
            m.patches.resize(n);
            for(auto& p:m.patches) { input>>label>>p.first>>p.second; require(label=="patch","missing patch"); }
            input>>label; require(input.good() && label=="end","truncated mesh");
            m.cells.resize(static_cast<std::size_t>(maxCell+1));
            const auto before=babelmesh::auditPolyhedralShells(m);
            std::size_t repairs=0;
            for(int pass=0;pass<4;++pass) {
                const auto changed=babelmesh::conformCollinearCellEdges(m,3.464101615e-11);
                repairs+=changed; if(!changed) break;
            }
            const auto after=babelmesh::auditPolyhedralShells(m);
            const auto merged=babelmesh::consolidateCoplanarFaces(m,1e-11,argc==4?std::stod(argv[3]):0.05);
            babelmesh::requirePolyhedralShells(m);
            std::cout<<"before_open="<<before.open_cells<<" after_open="<<after.open_cells
                     <<" before_halfspace="<<before.point_mean_halfspace_failures
                     <<" after_halfspace="<<after.point_mean_halfspace_failures<<" inserted="<<repairs
                     <<" merged="<<merged<<'\n';
            babelmesh::writeBabelSimV3(m,argv[2]);
            return 0;
        }
        auto m=cube();
        babelmesh::requirePolyhedralShells(m);
        auto tiled=m;
        auto bottom=tiled.faces.front();
        tiled.faces.erase(tiled.faces.begin());
        tiled.points_m.push_back({.5,.5,0});
        for(std::size_t j=0;j<4;++j)tiled.faces.push_back({{8,bottom.vertices[j],bottom.vertices[(j+1)%4]},0,-1,0});
        require(babelmesh::consolidateCoplanarFaces(tiled,1e-12)==3,
            "four coplanar fan triangles must become one physical face");
        babelmesh::requirePolyhedralShells(tiled);
        require(std::abs(babelmesh::auditMeshQuality(tiled).minimum_cell_volume-1)<1e-12,
            "coplanar consolidation must preserve cell volume");
        auto hangingCorner=cube();
        hangingCorner.points_m[1]={.5,.5,0};
        babelmesh::consolidateCoplanarFaces(hangingCorner,1e-12);
        require(!hangingCorner.strict_hex[0],
            "a collinear face corner is a hanging point, not a regular hex corner");
        babelmesh::requirePolyhedralShells(hangingCorner);
        require(babelmesh::auditMeshQuality(hangingCorner).minimum_cell_volume>0,
            "reclassified polyhedron must retain positive volume");
        auto poorJacobianHex=cube();
        for(auto& p:poorJacobianHex.points_m)p[0]+=20.0*p[2];
        babelmesh::consolidateCoplanarFaces(poorJacobianHex,1e-12);
        require(poorJacobianHex.strict_hex[0],
            "poor mapped quality must not alter the cell's topological classification");
        babelmesh::requirePolyhedralShells(poorJacobianHex);
        const auto poorJacobianQuality=babelmesh::auditMeshQuality(poorJacobianHex);
        require(poorJacobianQuality.regular_hex_cells==1 &&
                poorJacobianQuality.minimum_regular_hex_scaled_jacobian<0.1,
            "a poor mapped hex must remain visible to the strict quality gate");
        auto reconstructed=poorJacobianHex;
        const auto reconstructedParents=babelmesh::splitLowQualityHexesIntoPyramids(reconstructed);
        require(reconstructedParents==1 && reconstructed.cells.size()==6 &&
                std::count(reconstructed.strict_hex.begin(),reconstructed.strict_hex.end(),false)==6,
            "a low-quality hex must be replaced by six actual polyhedral pyramid cells");
        babelmesh::requirePolyhedralShells(reconstructed);
        const auto reconstructedQuality=babelmesh::auditMeshQuality(reconstructed);
        require(reconstructedQuality.negative_or_zero_volume_cells==0 &&
                reconstructedQuality.minimum_cell_volume>0 && reconstructedQuality.degenerate_faces==0,
            "pyramid reconstruction must retain closed positive-volume cells");
        require(std::abs(reconstructedQuality.minimum_cell_volume-1.0/6.0)<1e-12,
            "pyramid reconstruction must preserve the parent volume");
        auto highNonOrthogonalityHex=cube();
        for(auto& p:highNonOrthogonalityHex.points_m)p[0]+=4.0*p[2];
        const auto highNonOrthogonalityQuality=babelmesh::auditMeshQuality(highNonOrthogonalityHex);
        require(highNonOrthogonalityQuality.minimum_regular_hex_scaled_jacobian>0.1 &&
                highNonOrthogonalityQuality.maximum_non_orthogonality_degrees>75.0,
            "non-orthogonality fixture must be independent of the low-Jacobian reconstruction trigger");
        auto faceReconstructed=highNonOrthogonalityHex;
        require(babelmesh::splitLowQualityHexesIntoPyramids(faceReconstructed,0.1,
                std::numeric_limits<std::size_t>::max(),75.0)==1,
            "a cell incident to an over-threshold non-orthogonal face must be reconstructed");
        babelmesh::requirePolyhedralShells(faceReconstructed);
        const auto faceReconstructedQuality=babelmesh::auditMeshQuality(faceReconstructed);
        require(faceReconstructedQuality.negative_or_zero_volume_cells==0 &&
                faceReconstructedQuality.degenerate_faces==0 &&
                faceReconstructed.cells.size()==6,
            "face-targeted reconstruction must produce six valid positive-volume polyhedra");
        const auto cubeQuality=babelmesh::auditMeshQuality(m);
        require(cubeQuality.maximum_non_orthogonality_face_index<m.faces.size(),
            "quality audit must identify the face carrying maximum non-orthogonality");
        require(cubeQuality.minimum_regular_hex_scaled_jacobian>0.99,
            "regular hexahedron corner Jacobian orientation and scaling must be positive and normalized: "+
            std::to_string(cubeQuality.minimum_regular_hex_scaled_jacobian));
        auto reflected=m;
        for(auto& p:reflected.points_m) p[0]=-p[0];
        require(babelmesh::auditMeshQuality(reflected).minimum_regular_hex_scaled_jacobian<-.99,
            "reflection must retain negative Jacobian instead of masking it with an absolute value");
        auto sheared=m;
        for(auto& p:sheared.points_m) p[0]+=.5*p[1];
        require(std::abs(babelmesh::auditMeshQuality(sheared).minimum_regular_hex_scaled_jacobian-
            1.0/std::sqrt(1.25))<1e-12,"affine shear must match the analytic scaled Jacobian");
        m.points_m.push_back({.5,0,0});
        m.faces[2].vertices={0,8,1,5,4};
        auto a=babelmesh::auditPolyhedralShells(m);
        require(a.open_cells==1 && a.area_closure_failures==0,"area closure must not mask an unmatched subedge");
        require(babelmesh::conformCollinearCellEdges(m,1e-12)==1,"must insert midpoint into the other incident face");
        babelmesh::requirePolyhedralShells(m);
        require(babelmesh::conformCollinearCellEdges(m,1e-12)==0,"edge repair must be idempotent");
        m=cube(); m.points_m.push_back({.5,1e-5,0}); m.faces[2].vertices={0,8,1,5,4};
        require(babelmesh::conformCollinearCellEdges(m,1e-12)==0,"off-edge vertices must not be moved or inserted");
        m=cube(); std::reverse(m.faces[0].vertices.begin(),m.faces[0].vertices.end());
        a=babelmesh::auditPolyhedralShells(m);
        require(a.point_mean_halfspace_failures==1 && a.area_closure_failures==1,"inverted wall orientation must be rejected");
        bool rejected=false;
        try { babelmesh::requirePolyhedralShells(m); } catch(const std::runtime_error&) { rejected=true; }
        require(rejected,"strict validation must reject invalid shells");
        std::cout<<"Polyhedral shell audit and collinear conformance tests passed\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
