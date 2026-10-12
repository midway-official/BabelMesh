#include "babelmesh/polyhedral_audit.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <stdexcept>

namespace babelmesh {
namespace {
using P=std::array<double,3>;
P sub(P a,P b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
P cross(P a,P b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(P a,P b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
P normal(const HexMesh& m,const MeshFace& f){
    P n{},p=m.points_m[static_cast<std::size_t>(f.vertices[0])];
    for(std::size_t i=1;i+1<f.vertices.size();++i){
        auto c=cross(sub(m.points_m[static_cast<std::size_t>(f.vertices[i])],p),sub(m.points_m[static_cast<std::size_t>(f.vertices[i+1])],p));
        for(std::size_t k=0;k<3;++k)n[k]+=c[k];
    }
    return n;
}
bool rotateToReaderSafeFan(const HexMesh& m,MeshFace& face){
    if(face.vertices.size()<3)return false;
    auto n=normal(m,face);const double nl=std::sqrt(dot(n,n));if(!(nl>0))return false;
    for(auto& x:n)x/=nl;
    double scale=1.0;std::array<double,3> lo{INFINITY,INFINITY,INFINITY},hi{-INFINITY,-INFINITY,-INFINITY};
    for(auto id:face.vertices){const auto& p=m.points_m[static_cast<std::size_t>(id)];for(std::size_t k=0;k<3;++k){lo[k]=std::min(lo[k],p[k]);hi[k]=std::max(hi[k],p[k]);scale=std::max(scale,hi[k]-lo[k]);}}
    const double tolerance=1e-13*scale*scale;
    double best=0;std::size_t root=0;
    for(std::size_t r=0;r<face.vertices.size();++r){
        double minimum=INFINITY;const auto p=m.points_m[static_cast<std::size_t>(face.vertices[r])];
        for(std::size_t j=1;j+1<face.vertices.size();++j){
            const auto u=sub(m.points_m[static_cast<std::size_t>(face.vertices[(r+j)%face.vertices.size()])],p);
            const auto v=sub(m.points_m[static_cast<std::size_t>(face.vertices[(r+j+1)%face.vertices.size()])],p);
            minimum=std::min(minimum,dot(cross(u,v),n));
        }
        if(minimum>best){best=minimum;root=r;}
    }
    if(best<=tolerance)return false;
    std::rotate(face.vertices.begin(),face.vertices.begin()+static_cast<std::ptrdiff_t>(root),face.vertices.end());
    return true;
}
bool mergeMany(const HexMesh& m,const std::vector<const MeshFace*>& faces,double tol,double relativeWarp,MeshFace& result){
    const auto& a=*faces.front();
    auto n=normal(m,a); const auto reference=n;
    const double warp=a.neighbour>=0 ? relativeWarp : 0.0;
    const double referenceLength=std::sqrt(dot(n,n));
    if(!(referenceLength>0))return false;
    for(std::size_t i=1;i<faces.size();++i){
        const auto nb=normal(m,*faces[i]);
        if(dot(reference,nb)<(warp>0 ? 0.0 : (1-1e-10)*referenceLength*std::sqrt(dot(nb,nb))))return false;
        if(warp>0)for(std::size_t k=0;k<3;++k)n[k]+=nb[k];
    }
    const double length=std::sqrt(dot(n,n));
    if(!(length>0))return false;
    tol=std::max(tol,warp*std::sqrt(length/2));
    for(auto& x:n)x/=length;
    auto origin=m.points_m[static_cast<std::size_t>(a.vertices[0])];
    for(const auto* f:faces)for(auto v:f->vertices)
        if(std::abs(dot(n,sub(m.points_m[static_cast<std::size_t>(v)],origin)))>tol)return false;
    std::map<std::pair<int,int>,int> edges;
    int cancelled=0;
    for(const auto* f:faces)for(std::size_t j=0;j<f->vertices.size();++j){
        int u=f->vertices[j],v=f->vertices[(j+1)%f->vertices.size()];
        auto reverse=edges.find({v,u});
        if(reverse!=edges.end()){edges.erase(reverse);++cancelled;}else if(!edges.emplace(std::make_pair(u,v),1).second)return false;
    }
    if(!cancelled || edges.size()<3)return false;
    std::map<int,int> next;
    for(auto& [edge,count]:edges){(void)count;if(!next.emplace(edge.first,edge.second).second)return false;}
    std::vector<int> ring;int start=next.begin()->first,v=start;
    do{ring.push_back(v);auto it=next.find(v);if(it==next.end())return false;v=it->second;
    }while(v!=start && ring.size()<=next.size());
    if(v!=start || ring.size()!=next.size())return false;
    result=a;result.vertices=std::move(ring);return rotateToReaderSafeFan(m,result);
}
}
std::size_t consolidateCoplanarFaces(HexMesh& m,double tolerance,double relativeWarp){
    if(!std::isfinite(tolerance) || tolerance<0 || !std::isfinite(relativeWarp) || relativeWarp<0)
        throw std::invalid_argument("face consolidation tolerances must be finite and nonnegative");
    std::map<std::tuple<int,int,int>,std::vector<MeshFace>> groups;
    for(auto& f:m.faces)groups[{f.owner,f.neighbour,f.patch}].push_back(std::move(f));
    m.faces.clear();std::size_t merged=0;
    for(auto& [key,faces]:groups){
        (void)key;
        // A complete interface can have a valid polygon even when an
        // intermediate pair leaves an unsuitable fan root at a hanging node.
        if(faces.size()>2 && faces.front().neighbour>=0){
            std::vector<const MeshFace*> candidates;
            for(const auto& f:faces)candidates.push_back(&f);
            MeshFace combined;
            if(mergeMany(m,candidates,tolerance,relativeWarp,combined)){
                merged+=faces.size()-1;faces={std::move(combined)};
            }
        }
        bool changed=true;
        while(changed){changed=false;
            for(std::size_t a=0;a<faces.size()&&!changed;++a)for(std::size_t b=a+1;b<faces.size();++b){
                MeshFace combined;
                if(mergeMany(m,{&faces[a],&faces[b]},tolerance,relativeWarp,combined)){
                    faces[a]=std::move(combined);faces.erase(faces.begin()+static_cast<std::ptrdiff_t>(b));
                    ++merged;changed=true;break;
                }
            }
        }
        for(auto& f:faces){
            // Existing boundary rings can already contain collinear hanging
            // vertices. Rotate them too: BabelSim's convex-face reader uses a
            // fan from ring[0], so an otherwise valid face can fail if that
            // root lies on a subdivided straight edge.
            (void)rotateToReaderSafeFan(m,f);
            m.faces.push_back(std::move(f));
        }
    }
    std::vector<std::size_t> counts(m.cells.size());
    for(const auto& f:m.faces){
        ++counts[static_cast<std::size_t>(f.owner)];
        if(f.neighbour>=0)++counts[static_cast<std::size_t>(f.neighbour)];
        bool quadrilateral=f.vertices.size()==4;
        // A triangular face with a collinear hanging point has four stored
        // indices, but is not a quadrilateral of a trilinear hexahedron.
        // Retain the hanging point for shell conformity and export the cell
        // through its actual polyhedral faces. Never mask a reversed turn.
        if(quadrilateral)for(std::size_t i=0;i<4;++i){
            const auto p=m.points_m[static_cast<std::size_t>(f.vertices[i])];
            const auto u=sub(m.points_m[static_cast<std::size_t>(f.vertices[(i+3)%4])],p);
            const auto v=sub(m.points_m[static_cast<std::size_t>(f.vertices[(i+1)%4])],p);
            const auto cr=cross(u,v);
            // A face with a nearly vanished corner turn makes a formally
            // regular hex numerically singular (scaled Jacobian ~ 0). Keep
            // that cell as a general polyhedron so the polygonal face, not a
            // degenerate trilinear corner, defines its finite-volume geometry.
            if(dot(u,v)<0 && dot(cr,cr)<=1e-12*dot(u,u)*dot(v,v))quadrilateral=false;
        }
        if(!quadrilateral){
            if(static_cast<std::size_t>(f.owner)<m.strict_hex.size())m.strict_hex[static_cast<std::size_t>(f.owner)]=false;
            if(f.neighbour>=0 && static_cast<std::size_t>(f.neighbour)<m.strict_hex.size())m.strict_hex[static_cast<std::size_t>(f.neighbour)]=false;
        }
    }
    for(std::size_t c=0;c<m.strict_hex.size();++c)if(counts[c]!=6)m.strict_hex[c]=false;
    // Quality must not change topological classification. A six-face,
    // eight-corner cell remains in the regular-hex audit when its mapped
    // Jacobian is poor; the caller's locked quality gate must reject it.
    // Reclassify above only when its face topology is no longer hexahedral.
    return merged;
}
}
