#include "babelmesh/surface_distance.hpp"
#include "babelmesh/fluid_domain.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace babelmesh {
namespace {
using Point = std::array<double,3>;
using Triangle = std::array<Point,3>;
Point sub(Point a,const Point& b) { for(std::size_t k=0;k<3;++k)a[k]-=b[k];return a; }
double dot(const Point& a,const Point& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double distance2(const Point& a,const Point& b) { const auto d=sub(a,b);return dot(d,d); }
Point along(Point a,const Point& b,double t) { for(std::size_t k=0;k<3;++k)a[k]+=t*b[k];return a; }
double triangleDistance2(const Point& p,const Triangle& t) {
    const auto a=t[0],b=t[1],c=t[2],ab=sub(b,a),ac=sub(c,a),ap=sub(p,a);
    const double d1=dot(ab,ap),d2=dot(ac,ap);
    if(d1<=0&&d2<=0)return distance2(p,a);
    const auto bp=sub(p,b);const double d3=dot(ab,bp),d4=dot(ac,bp);
    if(d3>=0&&d4<=d3)return distance2(p,b);
    const double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0)return distance2(p,along(a,ab,d1/(d1-d3)));
    const auto cp=sub(p,c);const double d5=dot(ab,cp),d6=dot(ac,cp);
    if(d6>=0&&d5<=d6)return distance2(p,c);
    const double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0)return distance2(p,along(a,ac,d2/(d2-d6)));
    const double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0)
        return distance2(p,along(b,sub(c,b),(d4-d3)/((d4-d3)+(d5-d6))));
    const double sum=va+vb+vc;
    if(!(sum>0))throw std::runtime_error("surface distance: degenerate boundary triangle");
    return distance2(p,along(along(a,ab,vb/sum),ac,vc/sum));
}
class TriangleTree {
    struct Node { Point low{},high{};std::size_t begin{},end{},left{},right{}; };
    const std::vector<Triangle>& triangles;
    std::vector<std::size_t> order;
    std::vector<Node> nodes;
    std::size_t build(std::size_t begin,std::size_t end) {
        const auto id=nodes.size();nodes.emplace_back();
        Node node;node.begin=begin;node.end=end;
        node.low.fill(std::numeric_limits<double>::infinity());
        node.high.fill(-std::numeric_limits<double>::infinity());
        for(auto i=begin;i<end;++i)for(const auto& p:triangles[order[i]])for(std::size_t k=0;k<3;++k){
            node.low[k]=std::min(node.low[k],p[k]);node.high[k]=std::max(node.high[k],p[k]);
        }
        if(end-begin>8){
            std::size_t axis=0;
            for(std::size_t k=1;k<3;++k)if(node.high[k]-node.low[k]>node.high[axis]-node.low[axis])axis=k;
            const auto middle=begin+(end-begin)/2;
            auto centroid=[&](std::size_t i){const auto& t=triangles[i];return t[0][axis]+t[1][axis]+t[2][axis];};
            std::nth_element(order.begin()+static_cast<std::ptrdiff_t>(begin),order.begin()+static_cast<std::ptrdiff_t>(middle),
                order.begin()+static_cast<std::ptrdiff_t>(end),[&](auto a,auto b){return centroid(a)<centroid(b);});
            node.left=build(begin,middle);node.right=build(middle,end);
        }
        nodes[id]=node;return id;
    }
    double bound(const Point& p,std::size_t i) const {
        double sum=0;const auto& node=nodes[i];
        for(std::size_t k=0;k<3;++k){const double d=std::max({node.low[k]-p[k],p[k]-node.high[k],0.0});sum+=d*d;}
        return sum;
    }
    void visit(const Point& p,std::size_t id,double& best) const {
        if(bound(p,id)>best)return;
        const auto& node=nodes[id];
        if(node.end-node.begin<=8){for(auto i=node.begin;i<node.end;++i)best=std::min(best,triangleDistance2(p,triangles[order[i]]));return;}
        auto a=node.left,b=node.right;if(bound(p,a)>bound(p,b))std::swap(a,b);
        visit(p,a,best);visit(p,b,best);
    }
public:
    explicit TriangleTree(const std::vector<Triangle>& input):triangles(input),order(input.size()) {
        if(input.empty())throw std::runtime_error("surface distance: empty boundary");
        std::iota(order.begin(),order.end(),0U);build(0,order.size());
    }
    double distance(const Point& p) const {double best=std::numeric_limits<double>::infinity();visit(p,0,best);return std::sqrt(best);}
};

template<class Evaluate>
void sample(const std::vector<Triangle>& triangles,double spacing,std::size_t budget,Evaluate evaluate,
            std::size_t& count,double& maximum) {
    std::set<Point> seen;
    auto measure=[&](const Point& p){
        if(!seen.insert(p).second)return;
        if(++count>budget)throw std::runtime_error("surface distance: sample budget exhausted");
        const double d=evaluate(p);
        if(!std::isfinite(d))throw std::runtime_error("surface distance: nonfinite sample distance");
        maximum=std::max(maximum,d);
    };
    for(const auto& t:triangles){
        const double edge=std::sqrt(std::max({distance2(t[0],t[1]),distance2(t[1],t[2]),distance2(t[2],t[0])}));
        const double divisions=std::max(1.0,std::ceil(edge/spacing));
        if(!std::isfinite(divisions)||divisions>static_cast<double>(budget))
            throw std::runtime_error("surface distance: sampling resolution exceeds budget");
        const auto n=static_cast<std::size_t>(divisions);
        if(static_cast<long double>(n+1)*(n+2)/2>static_cast<long double>(budget))
            throw std::runtime_error("surface distance: triangle sampling exceeds budget");
        const auto ab=sub(t[1],t[0]),ac=sub(t[2],t[0]);
        for(std::size_t i=0;i<=n;++i)for(std::size_t j=0;j<=n-i;++j)
            measure(along(along(t[0],ab,static_cast<double>(i)/static_cast<double>(n)),ac,static_cast<double>(j)/static_cast<double>(n)));
        measure(along(along(t[0],ab,1.0/3.0),ac,1.0/3.0));
    }
}
} // namespace

SurfaceDistanceAudit auditSurfaceDistance(const HexMesh& mesh,const CadModel& cad,double spacing,std::size_t budget) {
    if(!(spacing>0)||!std::isfinite(spacing)||budget==0)throw std::runtime_error("surface distance: invalid sampling controls");
    std::vector<Triangle> boundary;
    for(const auto& f:mesh.faces)if(f.neighbour<0){
        if(f.vertices.size()<3 || std::set<std::int32_t>(f.vertices.begin(),f.vertices.end()).size()!=f.vertices.size())
            throw std::runtime_error("surface distance: invalid or degenerate boundary face connectivity");
        std::size_t emitted=0;
        for(std::size_t j=1;j+1<f.vertices.size();++j){
            const Triangle t{mesh.points_m.at(static_cast<std::size_t>(f.vertices[0])),mesh.points_m.at(static_cast<std::size_t>(f.vertices[j])),mesh.points_m.at(static_cast<std::size_t>(f.vertices[j+1]))};
            const auto u=sub(t[1],t[0]),v=sub(t[2],t[0]);const Point area{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
            if(dot(area,area)<=1e-24*dot(u,u)*dot(v,v))continue;
            boundary.push_back(t);++emitted;
        }
        if(emitted==0)throw std::runtime_error("surface distance: degenerate boundary face");
    }
    for(const auto& t:boundary){
        for(const auto& p:t)for(double coordinate:p)if(!std::isfinite(coordinate))
            throw std::runtime_error("surface distance: nonfinite boundary coordinate");
        const auto u=sub(t[1],t[0]),v=sub(t[2],t[0]);
        const Point area{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
        if(!(dot(area,area)>0))throw std::runtime_error("surface distance: degenerate boundary triangle");
    }
    const TriangleTree tree(boundary);
    SurfaceDistanceAudit result;result.requested_sample_spacing_m=spacing;
    result.cad_tessellation_deflection_m=cad.fluidSurfaceMesh().linear_deflection_m;
    sample(boundary,spacing,budget,[&](const Point& p){return projectPointToFluidBoundary(cad,p).distance_m;},
        result.mesh_to_cad_samples,result.mesh_to_cad_max_m);
    std::vector<Triangle> reference;
    const auto& surface=cad.fluidSurfaceMesh();
    for(const auto& t:surface.triangles)reference.push_back({surface.vertices.at(t.vertices[0]),surface.vertices.at(t.vertices[1]),surface.vertices.at(t.vertices[2])});
    if(reference.empty())throw std::runtime_error("surface distance: empty CAD tessellation");
    sample(reference,spacing,budget,[&](const Point& p){return tree.distance(projectPointToFluidBoundary(cad,p).point_m);},
        result.cad_to_mesh_samples,result.cad_to_mesh_max_m);
    return result;
}
} // namespace babelmesh
