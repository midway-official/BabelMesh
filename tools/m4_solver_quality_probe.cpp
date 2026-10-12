#include "babelsim/mesh_io.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
using Json=nlohmann::json;
using babelsim::Vec3;
struct Samples {
    std::vector<double> angles, skew;
    int worst_angle{-1},worst_skew{-1};
    double max_angle{-1},max_skew{-1};
    void add(double a,double s,int face) {
        angles.push_back(a); if(a>max_angle){max_angle=a;worst_angle=face;}
        if(s>=0){skew.push_back(s);if(s>max_skew){max_skew=s;worst_skew=face;}}
    }
    Json report() const {
        auto stats=[](std::vector<double> v,double limit) -> Json {
            if(v.empty()) return nullptr;
            std::sort(v.begin(),v.end());
            return {{"count",v.size()},{"max",v.back()},
                {"p99",v[static_cast<std::size_t>(std::ceil(.99*static_cast<double>(v.size())))-1]},
                {"above_limit",std::count_if(v.begin(),v.end(),[&](double x){return x>limit;})}};
        };
        return {{"nonorth",stats(angles,75)},{"skew",stats(skew,.5)},
            {"worst_angle_face",worst_angle},{"worst_skew_face",worst_skew}};
    }
};
struct Group {Vec3 area{},center{};double weight{};int representative{-1};};
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("usage: m4_solver_quality_probe mesh.bsm");
        const auto m=babelsim::readMeshFile(argv[1]);m.validate();
        Samples all,internal,boundary,triangles,other,grouped;
        std::map<std::pair<int,int>,Group> groups;
        auto measure=[&](Vec3 sf,Vec3 cf,int o,int n){
            auto co=m.cellCentre(o);
            auto d=n>=0 ? m.cellCentre(n)-co : cf-co;
            const double angle=std::acos(std::clamp(babelsim::dot(sf,d)/(babelsim::norm(sf)*babelsim::norm(d)),-1.0,1.0))*180/std::acos(-1.0);
            double skew=-1;
            if(n>=0){auto q=co+d*(babelsim::dot(sf,cf-co)/babelsim::dot(sf,d));skew=babelsim::norm(cf-q)/std::sqrt(babelsim::norm(sf));}
            return std::pair<double,double>{angle,skew};
        };
        for(int f=0;f<m.faceCount();++f){
            auto sf=m.faceAreaVector(f),cf=m.faceCentre(f);int o=m.owner(f),n=m.neighbour(f);
            auto [a,s]=measure(sf,cf,o,n);all.add(a,s,f);
            (n>=0?internal:boundary).add(a,s,f);
            (m.faceVertices(f).size()==3?triangles:other).add(a,s,f);
            if(n>=0){auto& g=groups[{o,n}];g.area+=sf;g.center+=cf*m.faceArea(f);g.weight+=m.faceArea(f);g.representative=f;}
        }
        for(const auto& [ids,g]:groups){auto [a,s]=measure(g.area,g.center/g.weight,ids.first,ids.second);grouped.add(a,s,g.representative);}
        Json examples=Json::array();
        for(int f:{all.worst_angle,all.worst_skew}) if(f>=0){
            auto cf=m.faceCentre(f),co=m.cellCentre(m.owner(f));
            examples.push_back({{"face",f},{"owner",m.owner(f)},{"neighbour",m.neighbour(f)},
                {"vertex_count",m.faceVertices(f).size()},{"center",{cf.x,cf.y,cf.z}},
                {"owner_center",{co.x,co.y,co.z}}});
        }
        std::cout<<Json{{"mesh",argv[1]},{"cells",m.cellCount()},
            {"geometry","actual BabelSim face area centroids and signed-volume cell centroids"},
            {"skew_definition","distance(face centroid, cell-centre line intersection with face plane) / sqrt(norm(Sf))"},
            {"all",all.report()},{"internal",internal.report()},{"boundary",boundary.report()},
            {"triangular_faces",triangles.report()},{"other_faces",other.report()},
            {"combined_owner_neighbour_interfaces_diagnostic_only",grouped.report()},{"worst_examples",examples}}.dump(2)<<'\n';
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
