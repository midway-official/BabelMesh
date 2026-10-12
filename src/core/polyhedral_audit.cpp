#include "babelmesh/polyhedral_audit.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace babelmesh {
namespace {
using Point = std::array<double, 3>;
using Edge = std::pair<std::int32_t, std::int32_t>;
Point subtract(const Point& a, const Point& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
double dot(const Point& a, const Point& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Edge key(std::int32_t a, std::int32_t b) { return {std::min(a,b), std::max(a,b)}; }
using Incidence = std::vector<std::vector<std::pair<std::size_t,int>>>;
Incidence incidence(const HexMesh& mesh) {
    Incidence result(mesh.cells.size());
    for (std::size_t i=0; i<mesh.faces.size(); ++i) {
        const auto& f=mesh.faces[i];
        if (f.vertices.size()<3 || f.owner<0 || static_cast<std::size_t>(f.owner)>=result.size() ||
            f.neighbour==f.owner || f.neighbour < -1 ||
            (f.neighbour>=0 && static_cast<std::size_t>(f.neighbour)>=result.size()))
            throw std::runtime_error("polyhedral audit: invalid face incidence at face "+std::to_string(i));
        std::set<std::int32_t> unique;
        for (auto v:f.vertices) {
            if (v<0 || static_cast<std::size_t>(v)>=mesh.points_m.size() || !unique.insert(v).second)
                throw std::runtime_error("polyhedral audit: invalid or repeated face vertex at face "+std::to_string(i));
            for (auto x:mesh.points_m[static_cast<std::size_t>(v)]) if (!std::isfinite(x))
                throw std::runtime_error("polyhedral audit: nonfinite vertex");
        }
        result[static_cast<std::size_t>(f.owner)].emplace_back(i,1);
        if(f.neighbour>=0) result[static_cast<std::size_t>(f.neighbour)].emplace_back(i,-1);
    }
    return result;
}
bool open(const HexMesh& mesh, const std::vector<std::pair<std::size_t,int>>& fs) {
    std::map<Edge,std::pair<int,int>> counts;
    for(const auto& [i,s]:fs) {
        const auto& vs=mesh.faces[i].vertices;
        for(std::size_t j=0;j<vs.size();++j) {
            auto a=vs[j], b=vs[(j+1)%vs.size()];
            auto& count=counts[key(a,b)];
            ++count.first; count.second+=s*(a<b ? 1 : -1);
        }
    }
    return fs.empty() || std::any_of(counts.begin(),counts.end(),[](const auto& e){
        return e.second.first!=2 || e.second.second!=0;
    });
}
std::pair<Point,Point> geometry(const HexMesh& mesh,const MeshFace& f) {
    Point center{},area{};
    const auto& origin=mesh.points_m[static_cast<std::size_t>(f.vertices[0])];
    for(std::size_t j=0;j<f.vertices.size();++j) {
        const auto& p=mesh.points_m[static_cast<std::size_t>(f.vertices[j])];
        const auto a=subtract(p,origin);
        const auto b=subtract(mesh.points_m[static_cast<std::size_t>(f.vertices[(j+1)%f.vertices.size()])],origin);
        for(std::size_t k=0;k<3;++k) {
            center[k]+=p[k]/static_cast<double>(f.vertices.size());
            area[k]+=.5*(a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3]);
        }
    }
    return {center,area};
}
} // namespace

PolyhedralAudit auditPolyhedralShells(const HexMesh& mesh) {
    const auto incident=incidence(mesh);
    PolyhedralAudit result;
    for(std::size_t c=0;c<incident.size();++c) {
        const auto& fs=incident[c];
        if(open(mesh,fs)) { ++result.open_cells; result.first_open_cell=std::min(result.first_open_cell,c); }
        Point center{},closure{};
        std::size_t n=0;
        for(const auto& [i,s]:fs) { (void)s; for(auto v:mesh.faces[i].vertices) {
            for(std::size_t k=0;k<3;++k) center[k]+=mesh.points_m[static_cast<std::size_t>(v)][k];
            ++n;
        }}
        if(n==0) continue;
        for(auto& x:center) x/=static_cast<double>(n);
        bool bad=false;
        double totalArea=0;
        for(const auto& [i,s]:fs) {
            const auto [fc,area]=geometry(mesh,mesh.faces[i]);
            bad=bad || s*dot(area,subtract(fc,center))<=0;
            totalArea+=std::sqrt(dot(area,area));
            for(std::size_t k=0;k<3;++k) closure[k]+=s*area[k];
        }
        if(bad) { ++result.point_mean_halfspace_failures; result.first_halfspace_cell=std::min(result.first_halfspace_cell,c); }
        if(std::sqrt(dot(closure,closure))>1e-9*std::max(totalArea,1.0)) ++result.area_closure_failures;
    }
    return result;
}

MeshQualityAudit auditMeshQuality(const HexMesh& mesh) {
    const auto incident=incidence(mesh);
    MeshQualityAudit q;
    q.minimum_cell_volume=std::numeric_limits<double>::infinity();
    q.minimum_regular_hex_scaled_jacobian=std::numeric_limits<double>::infinity();
    q.non_orthogonality_by_face.assign(mesh.faces.size(),0.0);
    std::vector<Point> centers(mesh.cells.size());
    for(std::size_t c=0;c<incident.size();++c) {
        const auto& fs=incident[c];
        std::set<std::int32_t> ids;
        for(const auto& [i,s]:fs) { (void)s; ids.insert(mesh.faces[i].vertices.begin(),mesh.faces[i].vertices.end()); }
        for(auto id:ids) for(std::size_t k=0;k<3;++k) centers[c][k]+=mesh.points_m[static_cast<std::size_t>(id)][k]/static_cast<double>(ids.size());
        long double volume=0;
        const auto reference=centers[c];
        Point moment{};
        for(const auto& [i,s]:fs) {
            const auto& f=mesh.faces[i];
            const auto& p0=mesh.points_m[static_cast<std::size_t>(f.vertices[0])];
            for(std::size_t j=1;j+1<f.vertices.size();++j) {
                const auto& a=mesh.points_m[static_cast<std::size_t>(f.vertices[j])];
                const auto& b=mesh.points_m[static_cast<std::size_t>(f.vertices[j+1])];
                const auto ra=subtract(a,reference),rb=subtract(b,reference),rp=subtract(p0,reference);
                const Point cr{ra[1]*rb[2]-ra[2]*rb[1],ra[2]*rb[0]-ra[0]*rb[2],ra[0]*rb[1]-ra[1]*rb[0]};
                const double tetra=static_cast<double>(s)*dot(rp,cr)/6.0;
                volume+=tetra;
                for(std::size_t k=0;k<3;++k) moment[k]+=tetra*(rp[k]+ra[k]+rb[k])/4.0;
            }
        }
        const double v=static_cast<double>(volume);
        q.cell_volumes.push_back(v); q.minimum_cell_volume=std::min(q.minimum_cell_volume,v);
        if(!(v>0) || !std::isfinite(v)) ++q.negative_or_zero_volume_cells;
        else for(std::size_t k=0;k<3;++k)centers[c][k]=reference[k]+moment[k]/v;
    }
    for(std::size_t faceIndex=0;faceIndex<mesh.faces.size();++faceIndex) {
        const auto& f=mesh.faces[faceIndex];
        auto [fc,sf]=geometry(mesh,f);
        Point faceMoment{};double triangleArea=0;
        const auto p0=mesh.points_m[static_cast<std::size_t>(f.vertices[0])];
        for(std::size_t j=1;j+1<f.vertices.size();++j){
            const auto a=mesh.points_m[static_cast<std::size_t>(f.vertices[j])],b=mesh.points_m[static_cast<std::size_t>(f.vertices[j+1])];
            const auto u=subtract(a,p0),v=subtract(b,p0);
            const Point cr{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
            const double w=.5*std::sqrt(dot(cr,cr));triangleArea+=w;
            for(std::size_t k=0;k<3;++k)faceMoment[k]+=w*(p0[k]+a[k]+b[k])/3.0;
        }
        if(triangleArea>0)for(std::size_t k=0;k<3;++k)fc[k]=faceMoment[k]/triangleArea;
        const double area=std::sqrt(dot(sf,sf));
        if(!(area>1e-30) || !std::isfinite(area)) { ++q.degenerate_faces; continue; }
        const auto co=centers[static_cast<std::size_t>(f.owner)];
        const Point d=f.neighbour>=0 ? subtract(centers[static_cast<std::size_t>(f.neighbour)],co) : subtract(fc,co);
        const double dlen=std::sqrt(dot(d,d));
        if(!(dlen>1e-30)) { ++q.degenerate_faces; continue; }
        const double cosine=std::clamp(dot(sf,d)/(area*dlen),-1.0,1.0);
        const double nonOrthogonality=std::acos(cosine)*180.0/std::acos(-1.0);
        q.non_orthogonality.push_back(nonOrthogonality);
        q.non_orthogonality_by_face[faceIndex]=nonOrthogonality;
        if(q.maximum_non_orthogonality_face_index==static_cast<std::size_t>(-1) ||
           nonOrthogonality>q.maximum_non_orthogonality_degrees) {
            q.maximum_non_orthogonality_degrees=nonOrthogonality;
            q.maximum_non_orthogonality_face_index=faceIndex;
        }
        const double weight=f.neighbour>=0 ? .5*(std::abs(q.cell_volumes[static_cast<std::size_t>(f.owner)])+
            std::abs(q.cell_volumes[static_cast<std::size_t>(f.neighbour)])) : std::abs(q.cell_volumes[static_cast<std::size_t>(f.owner)]);
        q.non_orthogonality_weights.push_back(weight);
        if(f.neighbour>=0) {
            q.internal_skewness_weights.push_back(weight);
            const double denom=dot(sf,d);
            if(std::abs(denom)<=1e-30) q.internal_skewness.push_back(std::numeric_limits<double>::infinity());
            else {
                const double t=dot(sf,subtract(fc,co))/denom;
                Point error=subtract(fc,{co[0]+t*d[0],co[1]+t*d[1],co[2]+t*d[2]});
                q.internal_skewness.push_back(std::sqrt(dot(error,error))/std::sqrt(area));
            }
        }
    }
    // Sort a copy: the original samples must stay aligned with face weights.
    auto summarize=[](std::vector<double> values,double& maximum,double& p99) {
        if(values.empty()) { maximum=0;p99=0;return; }
        std::sort(values.begin(),values.end());
        maximum=values.back();
        p99=values[static_cast<std::size_t>(std::ceil(.99*static_cast<double>(values.size())))-1U];
    };
    summarize(q.non_orthogonality,q.maximum_non_orthogonality_degrees,q.p99_non_orthogonality_degrees);
    summarize(q.internal_skewness,q.maximum_internal_skewness,q.p99_internal_skewness);
    for(std::size_t c=0;c<mesh.strict_hex.size() && c<mesh.cells.size();++c) if(mesh.strict_hex[c]) {
        const auto& ids=mesh.cells[c];
        double minimum=std::numeric_limits<double>::infinity();
        for(std::size_t corner=0;corner<8;++corner) {
            const auto p=mesh.points_m[static_cast<std::size_t>(ids[corner])];
            auto u=subtract(mesh.points_m[static_cast<std::size_t>(ids[corner^1U])],p);
            auto v=subtract(mesh.points_m[static_cast<std::size_t>(ids[corner^3U])],p);
            auto w=subtract(mesh.points_m[static_cast<std::size_t>(ids[corner^4U])],p);
            // VTK hex order: 000,100,110,010,001,101,111,011.
            if((corner&1U)^((corner>>1U)&1U)) for(auto& x:u) x=-x;
            if(corner&2U) for(auto& x:v) x=-x;
            if(corner&4U) for(auto& x:w) x=-x;
            const auto cr=Point{v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]};
            const double denom=std::sqrt(dot(u,u)*dot(v,v)*dot(w,w));
            minimum=std::min(minimum,denom>0 ? dot(u,cr)/denom : -1.0);
        }
        ++q.regular_hex_cells;
        q.minimum_regular_hex_scaled_jacobian=std::min(q.minimum_regular_hex_scaled_jacobian,minimum);
    }
    if(!std::isfinite(q.minimum_cell_volume)) q.minimum_cell_volume=0;
    if(!std::isfinite(q.minimum_regular_hex_scaled_jacobian)) q.minimum_regular_hex_scaled_jacobian=0;
    return q;
}

std::size_t conformCollinearCellEdges(HexMesh& mesh,double tolerance) {
    if(!std::isfinite(tolerance) || tolerance<=0) throw std::runtime_error("edge conformance requires positive tolerance");
    const auto incident=incidence(mesh);
    std::map<Edge,std::map<double,std::int32_t>> splits;
    for(const auto& fs:incident) {
        if(!open(mesh,fs)) continue;
        std::set<std::int32_t> vertices;
        for(const auto& [i,s]:fs) { (void)s; vertices.insert(mesh.faces[i].vertices.begin(),mesh.faces[i].vertices.end()); }
        for(const auto& [i,s]:fs) { (void)s;
            const auto& ring=mesh.faces[i].vertices;
            for(std::size_t j=0;j<ring.size();++j) {
                const auto e=key(ring[j],ring[(j+1)%ring.size()]);
                const auto& a=mesh.points_m[static_cast<std::size_t>(e.first)];
                const auto d=subtract(mesh.points_m[static_cast<std::size_t>(e.second)],a);
                const auto len2=dot(d,d);
                if(len2<=tolerance*tolerance) continue;
                for(auto v:vertices) {
                    if(std::find(ring.begin(),ring.end(),v)!=ring.end()) continue;
                    auto r=subtract(mesh.points_m[static_cast<std::size_t>(v)],a);
                    const double t=dot(r,d)/len2;
                    if(t<=tolerance/std::sqrt(len2) || t>=1-tolerance/std::sqrt(len2)) continue;
                    for(std::size_t k=0;k<3;++k) r[k]-=t*d[k];
                    if(dot(r,r)<=tolerance*tolerance) splits[e].emplace(t,v);
                }
            }
        }
    }
    std::size_t inserted=0;
    for(auto& f:mesh.faces) {
        std::vector<std::int32_t> ring;
        for(std::size_t j=0;j<f.vertices.size();++j) {
            auto a=f.vertices[j],b=f.vertices[(j+1)%f.vertices.size()];
            ring.push_back(a);
            const auto found=splits.find(key(a,b));
            if(found==splits.end()) continue;
            auto append=[&](const auto& item){
                if(std::find(f.vertices.begin(),f.vertices.end(),item.second)==f.vertices.end()) { ring.push_back(item.second); ++inserted; }
            };
            if(a<b) for(const auto& item:found->second) append(item);
            else for(auto it=found->second.rbegin();it!=found->second.rend();++it) append(*it);
        }
        if(ring.size()!=f.vertices.size()) {
            if(static_cast<std::size_t>(f.owner)<mesh.strict_hex.size()) mesh.strict_hex[static_cast<std::size_t>(f.owner)]=false;
            if(f.neighbour>=0 && static_cast<std::size_t>(f.neighbour)<mesh.strict_hex.size()) mesh.strict_hex[static_cast<std::size_t>(f.neighbour)]=false;
            f.vertices=std::move(ring);
        }
    }
    return inserted;
}

void requirePolyhedralShells(const HexMesh& mesh) {
    const auto a=auditPolyhedralShells(mesh);
    if(a.open_cells || a.area_closure_failures || a.point_mean_halfspace_failures)
        throw std::runtime_error("CAD post-snap shell validation failed: open_cells="+std::to_string(a.open_cells)+
            ", area_closure_failures="+std::to_string(a.area_closure_failures)+
            ", point_mean_halfspace_failures="+std::to_string(a.point_mean_halfspace_failures)+
            ", first_open_cell="+(a.open_cells ? std::to_string(a.first_open_cell) : "none")+
            ", first_halfspace_cell="+(a.point_mean_halfspace_failures ? std::to_string(a.first_halfspace_cell) : "none")+
            "; local topology reconstruction or constrained untangling required");
}
} // namespace babelmesh
