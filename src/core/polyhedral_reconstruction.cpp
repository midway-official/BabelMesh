#include "babelmesh/polyhedral_reconstruction.hpp"

#include "babelmesh/polyhedral_audit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace babelmesh {
namespace {
using Point = std::array<double, 3>;
using Edge = std::pair<std::int32_t, std::int32_t>;

Point sub(const Point& a, const Point& b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
double dot(const Point& a, const Point& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Point cross(const Point& a, const Point& b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
Edge edgeKey(std::int32_t a, std::int32_t b) { return {std::min(a,b),std::max(a,b)}; }

Point faceArea(const HexMesh& mesh, const MeshFace& face) {
    Point result{};
    const auto& origin=mesh.points_m.at(static_cast<std::size_t>(face.vertices.front()));
    for(std::size_t i=0;i<face.vertices.size();++i) {
        const auto a=sub(mesh.points_m.at(static_cast<std::size_t>(face.vertices[i])),origin);
        const auto b=sub(mesh.points_m.at(static_cast<std::size_t>(face.vertices[(i+1U)%face.vertices.size()])),origin);
        const auto term=cross(a,b);
        for(std::size_t axis=0;axis<3;++axis) result[axis]+=.5*term[axis];
    }
    return result;
}

double scaledJacobian(const HexMesh& mesh,const std::array<std::int32_t,8>& ids) {
    double minimum=std::numeric_limits<double>::infinity();
    for(std::size_t corner=0;corner<8;++corner) {
        const auto p=mesh.points_m.at(static_cast<std::size_t>(ids[corner]));
        auto u=sub(mesh.points_m.at(static_cast<std::size_t>(ids[corner^1U])),p);
        auto v=sub(mesh.points_m.at(static_cast<std::size_t>(ids[corner^3U])),p);
        auto w=sub(mesh.points_m.at(static_cast<std::size_t>(ids[corner^4U])),p);
        if((corner&1U)^((corner>>1U)&1U)) for(auto& x:u)x=-x;
        if(corner&2U)for(auto& x:v)x=-x;
        if(corner&4U)for(auto& x:w)x=-x;
        const double denominator=std::sqrt(dot(u,u)*dot(v,v)*dot(w,w));
        minimum=std::min(minimum,denominator>0?dot(u,cross(v,w))/denominator:-1.0);
    }
    return minimum;
}

Point meanPoints(const HexMesh& mesh,const std::set<std::int32_t>& ids) {
    Point result{};
    if(ids.empty())return result;
    for(auto id:ids)for(std::size_t axis=0;axis<3;++axis)
        result[axis]+=mesh.points_m.at(static_cast<std::size_t>(id))[axis]/static_cast<double>(ids.size());
    return result;
}
} // namespace

std::size_t splitLowQualityHexesIntoPyramids(HexMesh& mesh,double minimumJacobian,
        std::size_t maximumCells,double maximumFaceNonOrthogonalityDegrees) {
    if(!std::isfinite(minimumJacobian))throw std::runtime_error("polyhedral reconstruction threshold must be finite");
    if(std::isnan(maximumFaceNonOrthogonalityDegrees))throw std::runtime_error("polyhedral reconstruction non-orthogonality threshold must not be NaN");
    if(mesh.cells.empty()||mesh.strict_hex.size()!=mesh.cells.size())return 0;
    std::vector<std::vector<std::size_t>> incident(mesh.cells.size());
    for(std::size_t fi=0;fi<mesh.faces.size();++fi) {
        const auto& f=mesh.faces[fi];
        if(f.owner<0||static_cast<std::size_t>(f.owner)>=mesh.cells.size()||f.neighbour < -1||
           (f.neighbour>=0&&static_cast<std::size_t>(f.neighbour)>=mesh.cells.size()))return 0;
        incident[static_cast<std::size_t>(f.owner)].push_back(fi);
        if(f.neighbour>=0)incident[static_cast<std::size_t>(f.neighbour)].push_back(fi);
    }
    const auto originalQuality=auditMeshQuality(mesh);
    std::vector<bool> split(mesh.cells.size(),false);
    std::vector<Point> centers(mesh.cells.size());
    std::vector<std::size_t> childrenPerCell(mesh.cells.size(),0U);
    std::size_t selected=0;
    for(std::size_t c=0;c<mesh.cells.size();++c) {
        const bool lowJacobian=mesh.strict_hex[c]&&scaledJacobian(mesh,mesh.cells[c])<minimumJacobian;
        const bool highFaceNonOrthogonality=std::any_of(incident[c].begin(),incident[c].end(),[&](std::size_t fi) {
            return fi<originalQuality.non_orthogonality_by_face.size() &&
                originalQuality.non_orthogonality_by_face[fi]>maximumFaceNonOrthogonalityDegrees;
        });
        if(!lowJacobian&&!highFaceNonOrthogonality)continue;
        if(incident[c].size()<4U||incident[c].size()>32U)continue;
        bool valid=true;
        std::set<std::int32_t> vertices;
        for(auto fi:incident[c]) {
            const auto& f=mesh.faces[fi];
            if(f.vertices.size()<3U||f.vertices.size()+1U>8U||
               std::set<std::int32_t>(f.vertices.begin(),f.vertices.end()).size()!=f.vertices.size()) {valid=false;break;}
            vertices.insert(f.vertices.begin(),f.vertices.end());
        }
        if(!valid||vertices.size()<4U||vertices.size()>64U)continue;
        centers[c]=meanPoints(mesh,vertices);
        double scale=0;
        for(const auto& p:vertices)scale=std::max(scale,std::sqrt(dot(sub(mesh.points_m[static_cast<std::size_t>(p)],centers[c]),sub(mesh.points_m[static_cast<std::size_t>(p)],centers[c]))));
        for(auto fi:incident[c]) {
            const auto& f=mesh.faces[fi];
            const auto area=faceArea(mesh,f);
            Point fc{};
            for(auto v:f.vertices)for(std::size_t axis=0;axis<3;++axis)
                fc[axis]+=mesh.points_m[static_cast<std::size_t>(v)][axis]/static_cast<double>(f.vertices.size());
            const double orientation=f.owner==static_cast<std::int32_t>(c)?1.0:-1.0;
            const double inward=orientation*dot(area,sub(fc,centers[c]));
            if(!(inward>1.0e-12*std::max(1.0,scale)*std::sqrt(dot(area,area)))){valid=false;break;}
        }
        if(valid){split[c]=true;childrenPerCell[c]=incident[c].size();++selected;}
    }
    if(selected==0U)return 0;
    std::size_t additionalCells=0;
    for(std::size_t c=0;c<split.size();++c)if(split[c])additionalCells+=childrenPerCell[c]-1U;
    if(mesh.cells.size()+additionalCells>maximumCells)return 0;
    if(mesh.cells.size()+additionalCells>static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))return 0;

    HexMesh candidate=mesh;
    std::vector<std::int32_t> centerId(mesh.cells.size(),-1);
    for(std::size_t c=0;c<split.size();++c)if(split[c]) {
        centerId[c]=static_cast<std::int32_t>(candidate.points_m.size());
        candidate.points_m.push_back(centers[c]);
    }
    std::vector<std::array<std::int32_t,8>> cells;
    std::vector<bool> strict;
    std::vector<std::int32_t> levels;
    std::vector<double> references,aspects;
    std::map<std::string,std::vector<std::uint8_t>> markers;
    cells.reserve(mesh.cells.size()+additionalCells);strict.reserve(cells.capacity());
    levels.reserve(cells.capacity());references.reserve(cells.capacity());aspects.reserve(cells.capacity());
    for(const auto& item:mesh.sizing_source_markers)markers[item.first].reserve(cells.capacity());
    std::vector<std::int32_t> unchangedId(mesh.cells.size(),-1);
    std::vector<std::map<std::size_t,std::int32_t>> childForFace(mesh.cells.size());
    std::vector<Point> childCenters;
    childCenters.reserve(cells.capacity());
    for(std::size_t c=0;c<mesh.cells.size();++c) {
        if(!split[c]) {
            unchangedId[c]=static_cast<std::int32_t>(cells.size());
            cells.push_back(mesh.cells[c]);strict.push_back(mesh.strict_hex[c]);
            Point childCenter{};
            std::set<std::int32_t> ids(mesh.cells[c].begin(),mesh.cells[c].end());
            for(auto id:ids)for(std::size_t axis=0;axis<3;++axis)
                childCenter[axis]+=mesh.points_m[static_cast<std::size_t>(id)][axis]/static_cast<double>(ids.size());
            childCenters.push_back(childCenter);
            if(c<mesh.cell_levels.size())levels.push_back(mesh.cell_levels[c]);
            if(c<mesh.cell_volumes_m3.size())references.push_back(mesh.cell_volumes_m3[c]);
            if(c<mesh.cell_aspect_ratios.size())aspects.push_back(mesh.cell_aspect_ratios[c]);
            for(const auto& item:mesh.sizing_source_markers)markers[item.first].push_back(c<item.second.size()?item.second[c]:0U);
            continue;
        }
        const double parentVolume=c<mesh.cell_volumes_m3.size()?mesh.cell_volumes_m3[c]:0.0;
        for(auto fi:incident[c]) {
            const auto& ring=mesh.faces[fi].vertices;
            std::array<std::int32_t,8> pyramid{};
            std::copy(ring.begin(),ring.end(),pyramid.begin());
            std::fill(pyramid.begin()+static_cast<std::ptrdiff_t>(ring.size()),pyramid.end(),centerId[c]);
            const auto id=static_cast<std::int32_t>(cells.size());
            childForFace[c][fi]=id;
            cells.push_back(pyramid);strict.push_back(false);
            Point childCenter{};
            const double divisor=static_cast<double>(ring.size()+1U);
            for(auto vertex:ring)for(std::size_t axis=0;axis<3;++axis)
                childCenter[axis]+=mesh.points_m[static_cast<std::size_t>(vertex)][axis]/divisor;
            for(std::size_t axis=0;axis<3;++axis)
                childCenter[axis]+=centers[c][axis]/divisor;
            childCenters.push_back(childCenter);
            if(c<mesh.cell_levels.size())levels.push_back(mesh.cell_levels[c]);
            if(c<mesh.cell_volumes_m3.size())references.push_back(parentVolume/static_cast<double>(childrenPerCell[c]));
            if(c<mesh.cell_aspect_ratios.size())aspects.push_back(mesh.cell_aspect_ratios[c]);
            for(const auto& item:mesh.sizing_source_markers)markers[item.first].push_back(c<item.second.size()?item.second[c]:0U);
        }
    }

    std::vector<MeshFace> faces;
    std::size_t internalFaceCount=0;
    for(std::size_t c=0;c<split.size();++c)if(split[c]) {
        std::set<Edge> edges;
        for(auto fi:incident[c]) {
            const auto& ring=mesh.faces[fi].vertices;
            for(std::size_t j=0;j<ring.size();++j)edges.insert(edgeKey(ring[j],ring[(j+1U)%ring.size()]));
        }
        internalFaceCount+=edges.size();
    }
    faces.reserve(mesh.faces.size()+internalFaceCount);
    for(std::size_t fi=0;fi<mesh.faces.size();++fi) {
        auto face=mesh.faces[fi];
        const auto oldOwner=static_cast<std::size_t>(face.owner);
        face.owner=split[oldOwner]?childForFace[oldOwner].at(fi):unchangedId[oldOwner];
        if(face.neighbour>=0) {
            const auto oldNeighbour=static_cast<std::size_t>(face.neighbour);
            face.neighbour=split[oldNeighbour]?childForFace[oldNeighbour].at(fi):unchangedId[oldNeighbour];
        }
        faces.push_back(std::move(face));
    }
    for(std::size_t c=0;c<split.size();++c)if(split[c]) {
        std::map<Edge,std::vector<std::int32_t>> edgeChildren;
        for(auto fi:incident[c]) {
            const auto id=childForFace[c].at(fi);
            const auto& ring=mesh.faces[fi].vertices;
            for(std::size_t j=0;j<ring.size();++j)edgeChildren[edgeKey(ring[j],ring[(j+1U)%ring.size()])].push_back(id);
        }
        if(edgeChildren.empty())return 0;
        for(const auto& item:edgeChildren) {
            if(item.second.size()!=2U||item.second[0]==item.second[1])return 0;
            MeshFace internal{{centerId[c],item.first.first,item.first.second},item.second[0],item.second[1],-1};
            const auto ownerCenter=childCenters.at(static_cast<std::size_t>(internal.owner));
            const auto neighbourCenter=childCenters.at(static_cast<std::size_t>(internal.neighbour));
            const auto normal=faceArea(candidate,internal);
            if(dot(normal,sub(neighbourCenter,ownerCenter))<0)
                std::reverse(internal.vertices.begin(),internal.vertices.end());
            faces.push_back(std::move(internal));
        }
    }
    candidate.cells=std::move(cells);candidate.strict_hex=std::move(strict);
    if(!mesh.cell_levels.empty())candidate.cell_levels=std::move(levels);
    if(!mesh.cell_volumes_m3.empty())candidate.cell_volumes_m3=std::move(references);
    if(!mesh.cell_aspect_ratios.empty())candidate.cell_aspect_ratios=std::move(aspects);
    candidate.sizing_source_markers=std::move(markers);
    candidate.faces=std::move(faces);
    const auto shells=auditPolyhedralShells(candidate);
    const auto quality=auditMeshQuality(candidate);
    if(shells.open_cells||shells.area_closure_failures||shells.point_mean_halfspace_failures||
       quality.negative_or_zero_volume_cells||quality.degenerate_faces)return 0;
    if(!candidate.cell_volumes_m3.empty())candidate.cell_volumes_m3=quality.cell_volumes;
    mesh=std::move(candidate);
    return selected;
}
} // namespace babelmesh
