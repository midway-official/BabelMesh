#include "babelmesh/cad.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Builder.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_Curve.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColStd_SequenceOfAsciiString.hxx>
#include <TDataStd_Name.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDocStd_Document.hxx>
#include <Bnd_Box.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace babelmesh {
namespace {

using Json = nlohmann::json;

constexpr double kPi = 3.141592653589793238462643383279502884;

double unitScale(const std::string& unit) {
    if (unit == "m") return 1.0;
    if (unit == "mm") return 1.0e-3;
    if (unit == "cm") return 1.0e-2;
    if (unit == "inch") return 0.0254;
    throw std::runtime_error("unsupported length unit '" + unit + "'");
}

void checkKeys(const Json& object, const std::set<std::string>& allowed, const std::string& where) {
    if (!object.is_object()) throw std::runtime_error(where + " must be an object");
    for (auto it = object.begin(); it != object.end(); ++it) {
        if (allowed.find(it.key()) == allowed.end()) {
            throw std::runtime_error(where + " has unknown key '" + it.key() + "'");
        }
    }
}

void require(const Json& object, const char* key, const std::string& where) {
    if (!object.contains(key)) throw std::runtime_error(where + " is missing required key '" + key + "'");
}

std::array<double, 3> vector3(const Json& value, double scale, const std::string& where) {
    if (!value.is_array() || value.size() != 3) {
        throw std::runtime_error(where + " must be a three-component array");
    }
    std::array<double, 3> v{};
    for (std::size_t i = 0; i < 3; ++i) {
        if (!value[i].is_number()) throw std::runtime_error(where + " contains a non-number");
        v[i] = value[i].get<double>() * scale;
        if (!std::isfinite(v[i])) throw std::runtime_error(where + " contains a non-finite value");
    }
    return v;
}

gp_Pnt point(const Json& value, double scale, const std::string& where) {
    const auto v = vector3(value, scale, where);
    return gp_Pnt(v[0], v[1], v[2]);
}

gp_Dir direction(const Json& value, const std::string& where) {
    const auto v = vector3(value, 1.0, where);
    if (std::hypot(v[0], v[1], v[2]) < 1.0e-14) throw std::runtime_error(where + " must not be zero");
    return gp_Dir(v[0], v[1], v[2]);
}

const Json& paramsOf(const Json& operation) {
    static const Json empty = Json::object();
    return operation.contains("parameters") ? operation.at("parameters") : empty;
}

double number(const Json& params, const char* key, double scale, const std::string& where) {
    require(params, key, where);
    if (!params.at(key).is_number()) throw std::runtime_error(where + "." + key + " must be a number");
    const double result = params.at(key).get<double>() * scale;
    if (!std::isfinite(result)) throw std::runtime_error(where + "." + key + " must be finite");
    return result;
}

double positiveNumber(const Json& params, const char* key, double scale, const std::string& where) {
    const double result = number(params, key, scale, where);
    if (result <= 0.0) throw std::runtime_error(where + "." + key + " must be positive");
    return result;
}

TopoDS_Shape checkedShape(TopoDS_Shape shape, const std::string& id) {
    if (shape.IsNull()) throw std::runtime_error("CAD entity '" + id + "' produced a null shape");
    if (!BRepCheck_Analyzer(shape, Standard_True).IsValid()) {
        throw std::runtime_error("CAD entity '" + id + "' produced an invalid B-Rep");
    }
    return shape;
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string normalizeStepUnit(const std::string& raw) {
    const auto value = lower(raw);
    if (value.find("milli") != std::string::npos || value == "mm") return "mm";
    if (value.find("centi") != std::string::npos || value == "cm") return "cm";
    if (value.find("inch") != std::string::npos) return "inch";
    if (value == "m" || value.find("metre") != std::string::npos || value.find("meter") != std::string::npos) return "m";
    return "unknown";
}

std::mutex& stepMutex() {
    static std::mutex mutex;
    return mutex;
}

struct StepShape {
    TopoDS_Shape shape;
    std::map<std::string, std::vector<TopoDS_Shape>> names;
    std::vector<std::string> units;
};

StepShape readStep(const std::filesystem::path& file, const std::string& declaredUnit,
                   const std::string& wantedName) {
    if (!std::filesystem::is_regular_file(file)) {
        throw std::runtime_error("STEP file does not exist: " + file.string());
    }
    std::lock_guard<std::mutex> lock(stepMutex());
    STEPCAFControl_Reader reader;
    reader.SetNameMode(Standard_True);
    reader.SetColorMode(Standard_True);
    const Standard_CString previous = Interface_Static::CVal("xstep.cascade.unit");
    const std::string previousUnit = previous == nullptr ? "MM" : previous;
    if (!Interface_Static::SetCVal("xstep.cascade.unit", "M")) {
        throw std::runtime_error("OCCT cannot set STEP target unit to meters");
    }
    struct RestoreUnit {
        std::string value;
        ~RestoreUnit() { Interface_Static::SetCVal("xstep.cascade.unit", value.c_str()); }
    } restore{previousUnit};

    if (reader.ReadFile(file.string().c_str()) != IFSelect_RetDone) {
        throw std::runtime_error("OCCT failed to read STEP file: " + file.string());
    }
    TColStd_SequenceOfAsciiString lengthUnits, angleUnits, solidAngleUnits;
    reader.ChangeReader().FileUnits(lengthUnits, angleUnits, solidAngleUnits);
    std::vector<std::string> rawUnits;
    bool matchedDeclared = false;
    for (Standard_Integer i = 1; i <= lengthUnits.Length(); ++i) {
        const std::string raw = lengthUnits.Value(i).ToCString();
        rawUnits.push_back(raw);
        const std::string normalized = normalizeStepUnit(raw);
        if (normalized == "unknown") {
            throw std::runtime_error("STEP length unit is ambiguous or unsupported ('" + raw + "') in " + file.string());
        }
        if (normalized != declaredUnit) {
            throw std::runtime_error("STEP contains mixed length units ('" + raw + "') but import declares '" + declaredUnit + "' in " + file.string());
        }
        matchedDeclared = true;
    }
    if (rawUnits.empty()) throw std::runtime_error("STEP file has no declared length unit: " + file.string());
    if (!matchedDeclared) {
        throw std::runtime_error("STEP length unit does not match declared unit '" + declaredUnit + "' in " + file.string());
    }

    Handle(TDocStd_Document) document;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", document);
    if (!reader.Transfer(document)) throw std::runtime_error("OCCT could not transfer STEP product structure: " + file.string());
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(document->Main());
    TDF_LabelSequence freeShapes;
    shapeTool->GetFreeShapes(freeShapes);
    if (freeShapes.IsEmpty()) throw std::runtime_error("STEP file contains no free product shapes: " + file.string());

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    StepShape result;
    result.units = rawUnits;
    std::vector<TopoDS_Shape> selectedInstances;
    std::function<void(const TDF_Label&, const std::string&)> collectNames;
    collectNames = [&](const TDF_Label& label, const std::string& parentPath) {
        std::string name;
        Handle(TDataStd_Name) nameAttribute;
        if (label.FindAttribute(TDataStd_Name::GetID(), nameAttribute)) {
            TCollection_AsciiString utf8(nameAttribute->Get(), '?');
            name = utf8.ToCString();
        }
        if (name.empty()) name = "label-" + std::to_string(label.Tag());
        const std::string path = parentPath.empty() ? name : parentPath + "/" + name;
        const TopoDS_Shape namedShape = shapeTool->GetShape(label);
        if (!namedShape.IsNull()) {
            result.names[path].push_back(namedShape);
            if (path != name) result.names[name].push_back(namedShape);
            if (!wantedName.empty() && (wantedName == path || wantedName == name)) selectedInstances.push_back(namedShape);
        }
        TDF_LabelSequence components;
        if (XCAFDoc_ShapeTool::GetComponents(label, components, Standard_False)) {
            for (Standard_Integer j = 1; j <= components.Length(); ++j) collectNames(components.Value(j), path);
        }
    };
    for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
        const TDF_Label& label = freeShapes.Value(i);
        TopoDS_Shape root = shapeTool->GetShape(label);
        if (root.IsNull()) continue;
        if (wantedName.empty()) builder.Add(compound, root);
        collectNames(label, "");
    }
    if (wantedName.empty()) {
        if (freeShapes.IsEmpty()) throw std::runtime_error("STEP file has no transferable root shapes: " + file.string());
        result.shape = checkedShape(compound, file.filename().string());
    } else {
        if (selectedInstances.empty()) throw std::runtime_error("STEP instance/name '" + wantedName + "' was not found in " + file.string());
        if (selectedInstances.size() != 1) throw std::runtime_error("STEP instance/name '" + wantedName + "' is ambiguous in " + file.string());
        result.shape = checkedShape(selectedInstances.front(), file.filename().string());
    }
    return result;
}

TopAbs_ShapeEnum shapeType(const std::string& type) {
    if (type == "vertex") return TopAbs_VERTEX;
    if (type == "edge") return TopAbs_EDGE;
    if (type == "face") return TopAbs_FACE;
    if (type == "solid") return TopAbs_SOLID;
    throw std::runtime_error("unsupported physical-group entity_type '" + type + "'");
}

std::vector<TopoDS_Shape> enumerate(const TopoDS_Shape& shape, TopAbs_ShapeEnum type) {
    std::vector<TopoDS_Shape> items;
    for (TopExp_Explorer it(shape, type); it.More(); it.Next()) items.push_back(it.Current());
    return items;
}

gp_Pnt centerOfMass(const TopoDS_Shape& shape) {
    GProp_GProps props;
    if (shape.ShapeType() == TopAbs_VERTEX) return BRep_Tool::Pnt(TopoDS::Vertex(shape));
    if (shape.ShapeType() == TopAbs_EDGE) BRepGProp::LinearProperties(shape, props);
    else if (shape.ShapeType() == TopAbs_FACE) BRepGProp::SurfaceProperties(shape, props);
    else BRepGProp::VolumeProperties(shape, props);
    return props.CentreOfMass();
}

std::array<double, 6> boundsOf(const TopoDS_Shape& shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box, Standard_False);
    if (box.IsVoid()) throw std::runtime_error("shape has empty bounds");
    std::array<double, 6> b{};
    box.Get(b[0], b[1], b[2], b[3], b[4], b[5]);
    return b;
}

std::array<double, 3> wireAreaVector(const TopoDS_Wire& wire) {
    std::vector<std::array<double, 3>> points;
    constexpr int samplesPerEdge = 16;
    for (BRepTools_WireExplorer it(wire); it.More(); it.Next()) {
        const TopoDS_Edge edge = TopoDS::Edge(it.Current());
        TopLoc_Location location;
        Standard_Real first = 0.0, last = 0.0;
        const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, location, first, last);
        if (curve.IsNull() || !std::isfinite(first) || !std::isfinite(last) || last <= first) {
            throw std::runtime_error("planar face wire contains an edge without a finite 3D curve");
        }
        for (int i = 0; i < samplesPerEdge; ++i) {
            const double fraction = static_cast<double>(i) / static_cast<double>(samplesPerEdge);
            const double orientedFraction = edge.Orientation() == TopAbs_REVERSED ? 1.0 - fraction : fraction;
            const double parameter = first + orientedFraction * (last - first);
            gp_Pnt point = curve->Value(parameter).Transformed(location.Transformation());
            points.push_back({point.X(), point.Y(), point.Z()});
        }
    }
    if (points.size() < 3) throw std::runtime_error("planar face wire has fewer than three sampled points");
    const auto origin = points.front();
    std::array<double, 3> area{};
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& a = points[i];
        const auto& b = points[(i + 1U) % points.size()];
        const std::array<double, 3> u{a[0] - origin[0], a[1] - origin[1], a[2] - origin[2]};
        const std::array<double, 3> v{b[0] - origin[0], b[1] - origin[1], b[2] - origin[2]};
        area[0] += u[1] * v[2] - u[2] * v[1];
        area[1] += u[2] * v[0] - u[0] * v[2];
        area[2] += u[0] * v[1] - u[1] * v[0];
    }
    return area;
}

SurfaceMesh makeSurfaceMesh(const TopoDS_Shape& shape, const std::string& entityId) {
    SurfaceMesh result;
    const auto bounds = boundsOf(shape);
    const double diagonal = std::hypot(bounds[3] - bounds[0], bounds[4] - bounds[1], bounds[5] - bounds[2]);
    result.linear_deflection_m = std::max(diagonal * 1.0e-3, 1.0e-9);
    BRepMesh_IncrementalMesh mesher(shape, result.linear_deflection_m, Standard_False, 0.35, Standard_False);
    if (!mesher.IsDone()) throw std::runtime_error("OCCT surface triangulation failed for '" + entityId + "'");

    const double weldTolerance = std::max(result.linear_deflection_m * 1.0e-9, 1.0e-13);
    std::map<std::array<long long, 3>, std::size_t> vertexIndex;
    std::size_t faceIndex = 0;
    for (TopExp_Explorer faceIt(shape, TopAbs_FACE); faceIt.More(); faceIt.Next(), ++faceIndex) {
        const TopoDS_Face face = TopoDS::Face(faceIt.Current());
        TopLoc_Location location;
        const Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        if (triangulation.IsNull()) continue;
        std::vector<std::size_t> localToGlobal(static_cast<std::size_t>(triangulation->NbNodes()) + 1U);
        for (Standard_Integer node = 1; node <= triangulation->NbNodes(); ++node) {
            gp_Pnt p = triangulation->Node(node).Transformed(location.Transformation());
            const std::array<double, 3> xyz{p.X(), p.Y(), p.Z()};
            std::array<long long, 3> key{};
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const long double scaled = std::nearbyint(static_cast<long double>(xyz[axis]) / weldTolerance);
                if (std::abs(scaled) > static_cast<long double>(std::numeric_limits<long long>::max())) {
                    throw std::runtime_error("surface mesh coordinate exceeds vertex-welding key range");
                }
                key[axis] = static_cast<long long>(scaled);
            }
            auto found = vertexIndex.find(key);
            if (found == vertexIndex.end()) {
                const std::size_t index = result.vertices.size();
                result.vertices.push_back(xyz);
                found = vertexIndex.emplace(key, index).first;
            }
            localToGlobal[static_cast<std::size_t>(node)] = found->second;
        }
        for (Standard_Integer triangleIndex = 1; triangleIndex <= triangulation->NbTriangles(); ++triangleIndex) {
            Standard_Integer a = 0, b = 0, c = 0;
            triangulation->Triangle(triangleIndex).Get(a, b, c);
            if (face.Orientation() == TopAbs_REVERSED) std::swap(b, c);
            SurfaceTriangle triangle{{localToGlobal[static_cast<std::size_t>(a)],
                                      localToGlobal[static_cast<std::size_t>(b)],
                                      localToGlobal[static_cast<std::size_t>(c)]}, faceIndex};
            const auto& p0 = result.vertices[triangle.vertices[0]];
            const auto& p1 = result.vertices[triangle.vertices[1]];
            const auto& p2 = result.vertices[triangle.vertices[2]];
            const std::array<double, 3> e0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const std::array<double, 3> e1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const std::array<double, 3> cross{e0[1] * e1[2] - e0[2] * e1[1], e0[2] * e1[0] - e0[0] * e1[2], e0[0] * e1[1] - e0[1] * e1[0]};
            if (std::hypot(cross[0], cross[1], cross[2]) <= 1.0e-24) {
                ++result.degenerate_triangles_skipped;
                continue;
            }
            result.triangles.push_back(triangle);
        }
    }

    buildSurfaceIndex(result);

    return result;
}

bool insideBounds(const std::array<double, 6>& a, const Json& selector, double scale) {
    if (!selector.contains("bounds")) return true;
    const auto lo = vector3(selector.at("bounds").at("minimum"), scale, "selector.bounds.minimum");
    const auto hi = vector3(selector.at("bounds").at("maximum"), scale, "selector.bounds.maximum");
    return a[0] >= lo[0] && a[1] >= lo[1] && a[2] >= lo[2] &&
           a[3] <= hi[0] && a[4] <= hi[1] && a[5] <= hi[2];
}

std::vector<TopoDS_Shape> selectGroup(const CadModel& model, const Json& spec, double scale) {
    const Json& selector = spec.at("selector");
    checkKeys(selector, {"operation_id", "entity_type", "names", "center", "exclude_centers", "normal", "tolerance", "bounds", "adjacent_group"}, "group selector");
    if (selector.contains("adjacent_group")) {
        throw std::runtime_error("group selector adjacent_group is not implemented in M1; use explicit CAD selectors");
    }
    const TopAbs_ShapeEnum type = selector.contains("entity_type")
        ? shapeType(selector.at("entity_type").get<std::string>()) : TopAbs_SHAPE;
    std::vector<TopoDS_Shape> candidates;
    if (selector.contains("names")) {
        if (!selector.contains("operation_id")) throw std::runtime_error("named group selector requires operation_id");
        const auto entityIt = model.entities().find(selector.at("operation_id").get<std::string>());
        if (entityIt == model.entities().end()) throw std::runtime_error("group selector refers to unknown operation_id");
        for (const auto& name : selector.at("names")) {
            const auto found = entityIt->second.named_subshapes.find(name.get<std::string>());
            if (found == entityIt->second.named_subshapes.end()) {
                throw std::runtime_error("physical group name selector is ambiguous or absent: '" + name.get<std::string>() + "'");
            }
            for (const auto& namedShape : found->second) {
                if (type == TopAbs_SHAPE || namedShape.ShapeType() == type) candidates.push_back(namedShape);
                else {
                    const auto subshapes = enumerate(namedShape, type);
                    candidates.insert(candidates.end(), subshapes.begin(), subshapes.end());
                }
            }
        }
    } else {
        if (!selector.contains("operation_id")) throw std::runtime_error("physical group selector requires operation_id or names");
        const auto entityIt = model.entities().find(selector.at("operation_id").get<std::string>());
        if (entityIt == model.entities().end()) throw std::runtime_error("group selector refers to unknown operation_id");
        if (type == TopAbs_SHAPE) candidates.push_back(entityIt->second.shape);
        else candidates = enumerate(entityIt->second.shape, type);
    }

    const double tolerance = selector.value("tolerance", 1.0e-7 / scale) * scale;
    std::vector<TopoDS_Shape> selected;
    for (const auto& candidate : candidates) {
        if (candidate.IsNull() || (type != TopAbs_SHAPE && candidate.ShapeType() != type)) continue;
        if (!insideBounds(boundsOf(candidate), selector, scale)) continue;
        if (selector.contains("center")) {
            const gp_Pnt target = point(selector.at("center"), scale, "selector.center");
            if (centerOfMass(candidate).Distance(target) > tolerance) continue;
        }
        if (selector.contains("exclude_centers")) {
            bool excluded = false;
            for (const auto& value : selector.at("exclude_centers"))
                if (centerOfMass(candidate).Distance(point(value, scale, "selector.exclude_centers")) <= tolerance) excluded = true;
            if (excluded) continue;
        }
        if (selector.contains("normal")) {
            if (candidate.ShapeType() != TopAbs_FACE) throw std::runtime_error("normal selector requires face entities");
            const auto desired = direction(selector.at("normal"), "selector.normal");
            BRepAdaptor_Surface surface(TopoDS::Face(candidate), Standard_True);
            if (surface.GetType() != GeomAbs_Plane) continue;
            gp_Dir actual = surface.Plane().Axis().Direction();
            if (candidate.Orientation() == TopAbs_REVERSED) actual.Reverse();
            if (actual.Dot(desired) < 1.0 - std::max(1.0e-6, tolerance)) continue;
        }
        selected.push_back(candidate);
    }
    if (selected.empty()) throw std::runtime_error("physical group '" + spec.at("name").get<std::string>() + "' selected no subshapes");
    return selected;
}

TopoDS_Shape makeOperation(const Json& operation, CadModel& model, double scale, CadEntity& entity) {
    const std::string type = operation.at("operation").get<std::string>();
    const Json& p = paramsOf(operation);
    std::vector<TopoDS_Shape> inputs;
    for (const auto& input : operation.at("inputs")) {
        const std::string inputId = input.get<std::string>();
        const auto found = model.entities().find(inputId);
        if (found == model.entities().end()) throw std::runtime_error("CAD operation '" + entity.id + "' refers to unavailable input '" + inputId + "'");
        inputs.push_back(found->second.shape);
        entity.inputs.push_back(inputId);
    }
    const std::string where = "CAD operation '" + entity.id + "' (" + type + ")";
    if (type == "point") {
        require(p, "center", where + ".parameters");
        return BRepBuilderAPI_MakeVertex(point(p.at("center"), scale, where + ".center")).Vertex();
    }
    if (type == "line") {
        require(p, "start", where + ".parameters"); require(p, "end", where + ".parameters");
        return BRepBuilderAPI_MakeEdge(point(p.at("start"), scale, where + ".start"), point(p.at("end"), scale, where + ".end")).Edge();
    }
    if (type == "arc") {
        if (inputs.size() != 3 || inputs[0].ShapeType() != TopAbs_VERTEX || inputs[1].ShapeType() != TopAbs_VERTEX || inputs[2].ShapeType() != TopAbs_VERTEX) {
            throw std::runtime_error(where + " requires exactly three point-vertex inputs (start, on-arc, end)");
        }
        const auto arc = GC_MakeArcOfCircle(BRep_Tool::Pnt(TopoDS::Vertex(inputs[0])), BRep_Tool::Pnt(TopoDS::Vertex(inputs[1])), BRep_Tool::Pnt(TopoDS::Vertex(inputs[2])));
        if (!arc.IsDone()) throw std::runtime_error(where + " could not construct circular arc");
        return BRepBuilderAPI_MakeEdge(arc.Value()).Edge();
    }
    if (type == "circle") {
        require(p, "center", where + ".parameters"); require(p, "normal", where + ".parameters");
        const gp_Ax2 axis(point(p.at("center"), scale, where + ".center"), direction(p.at("normal"), where + ".normal"));
        return BRepBuilderAPI_MakeEdge(gp_Circ(axis, positiveNumber(p, "radius", scale, where + ".parameters"))).Edge();
    }
    if (type == "spline") {
        require(p, "points", where + ".parameters");
        if (p.at("points").size() < 2) throw std::runtime_error(where + " needs at least two spline points");
        Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, static_cast<Standard_Integer>(p.at("points").size()));
        for (std::size_t i = 0; i < p.at("points").size(); ++i) points->SetValue(static_cast<Standard_Integer>(i + 1), point(p.at("points")[i], scale, where + ".points"));
        GeomAPI_Interpolate interpolator(points, Standard_False, 1.0e-9);
        interpolator.Perform();
        if (!interpolator.IsDone()) throw std::runtime_error(where + " spline interpolation failed");
        return BRepBuilderAPI_MakeEdge(interpolator.Curve()).Edge();
    }
    if (type == "wire") {
        BRepBuilderAPI_MakeWire wire;
        if (inputs.empty()) throw std::runtime_error(where + " requires edge or wire inputs");
        for (const auto& input : inputs) {
            if (input.ShapeType() == TopAbs_EDGE) wire.Add(TopoDS::Edge(input));
            else if (input.ShapeType() == TopAbs_WIRE) wire.Add(TopoDS::Wire(input));
            else throw std::runtime_error(where + " accepts only edge and wire inputs");
        }
        if (!wire.IsDone()) throw std::runtime_error(where + " inputs do not form a connected wire");
        if (p.value("closed", false) && !BRep_Tool::IsClosed(wire.Wire())) throw std::runtime_error(where + " requested a closed wire but inputs are open");
        return wire.Wire();
    }
    if (type == "plane") {
        if (inputs.empty() || inputs[0].ShapeType() != TopAbs_WIRE) throw std::runtime_error(where + " requires an outer wire as its first input");
        const TopoDS_Wire outerWire = TopoDS::Wire(inputs[0]);
        const auto outerArea = wireAreaVector(outerWire);
        BRepBuilderAPI_MakeFace face(outerWire, Standard_True);
        for (std::size_t i = 1; i < inputs.size(); ++i) {
            if (inputs[i].ShapeType() != TopAbs_WIRE) throw std::runtime_error(where + " hole inputs must be wires");
            TopoDS_Wire hole = TopoDS::Wire(inputs[i]);
            const auto holeArea = wireAreaVector(hole);
            const double orientation = outerArea[0] * holeArea[0] + outerArea[1] * holeArea[1] + outerArea[2] * holeArea[2];
            if (std::abs(orientation) <= 1.0e-24) throw std::runtime_error(where + " has a degenerate or non-coplanar hole wire");
            if (orientation > 0.0) hole.Reverse();
            face.Add(hole);
        }
        if (!face.IsDone()) throw std::runtime_error(where + " failed to construct planar face with holes");
        return face.Face();
    }
    if (type == "shell") {
        if (inputs.empty()) throw std::runtime_error(where + " requires face inputs");
        const double sewingTolerance = p.contains("tolerance")
            ? positiveNumber(p, "tolerance", scale, where + ".parameters") : 1.0e-7;
        BRepBuilderAPI_Sewing sewing(sewingTolerance, Standard_True, Standard_True, Standard_True, Standard_False);
        for (const auto& input : inputs) {
            if (input.ShapeType() == TopAbs_FACE) sewing.Add(input);
            else for (TopExp_Explorer ex(input, TopAbs_FACE); ex.More(); ex.Next()) sewing.Add(ex.Current());
        }
        sewing.Perform();
        const TopoDS_Shape sewn = sewing.SewedShape();
        if (sewn.IsNull() || sewn.ShapeType() != TopAbs_SHELL) throw std::runtime_error(where + " did not produce one connected shell");
        return sewn;
    }
    if (type == "solid") {
        if (inputs.size() != 1 || inputs[0].ShapeType() != TopAbs_SHELL) throw std::runtime_error(where + " requires exactly one sewn shell input");
        const TopoDS_Shell shell = TopoDS::Shell(inputs[0]);
        if (!BRep_Tool::IsClosed(shell)) throw std::runtime_error(where + " cannot solidify an open shell");
        BRepBuilderAPI_MakeSolid solidMaker(shell);
        if (!solidMaker.IsDone()) throw std::runtime_error(where + " failed to build a solid from the closed shell");
        TopoDS_Solid solid = solidMaker.Solid();
        if (!BRepLib::OrientClosedSolid(solid)) throw std::runtime_error(where + " could not orient the closed shell as a solid");
        return solid;
    }
    if (type == "box") {
        require(p, "minimum", where + ".parameters"); require(p, "maximum", where + ".parameters");
        const gp_Pnt lo = point(p.at("minimum"), scale, where + ".minimum");
        const gp_Pnt hi = point(p.at("maximum"), scale, where + ".maximum");
        if (hi.X() <= lo.X() || hi.Y() <= lo.Y() || hi.Z() <= lo.Z()) throw std::runtime_error(where + " has non-positive box dimensions");
        return BRepPrimAPI_MakeBox(lo, hi).Solid();
    }
    if (type == "sphere") {
        require(p, "center", where + ".parameters");
        return BRepPrimAPI_MakeSphere(point(p.at("center"), scale, where + ".center"), positiveNumber(p, "radius", scale, where + ".parameters")).Solid();
    }
    if (type == "cylinder") {
        require(p, "axis_origin", where + ".parameters"); require(p, "direction", where + ".parameters");
        const gp_Ax2 axis(point(p.at("axis_origin"), scale, where + ".axis_origin"), direction(p.at("direction"), where + ".direction"));
        return BRepPrimAPI_MakeCylinder(axis, positiveNumber(p, "radius", scale, where + ".parameters"), positiveNumber(p, "height", scale, where + ".parameters")).Solid();
    }
    if (type == "cone") {
        require(p, "axis_origin", where + ".parameters"); require(p, "direction", where + ".parameters");
        const gp_Ax2 axis(point(p.at("axis_origin"), scale, where + ".axis_origin"), direction(p.at("direction"), where + ".direction"));
        const double bottom = number(p, "bottom_radius", scale, where + ".parameters");
        const double top = number(p, "top_radius", scale, where + ".parameters");
        if (bottom < 0.0 || top < 0.0 || (bottom == 0.0 && top == 0.0)) throw std::runtime_error(where + " cone radii must be nonnegative and not both zero");
        return BRepPrimAPI_MakeCone(axis, bottom, top, positiveNumber(p, "height", scale, where + ".parameters")).Solid();
    }
    if (type == "extrude") {
        if (inputs.size() != 1) throw std::runtime_error(where + " requires exactly one face or wire input");
        require(p, "direction", where + ".parameters");
        const auto v = vector3(p.at("direction"), scale, where + ".direction");
        if (std::hypot(v[0], v[1], v[2]) < 1.0e-14) throw std::runtime_error(where + ".direction must not be zero");
        return BRepPrimAPI_MakePrism(inputs[0], gp_Vec(v[0], v[1], v[2]), Standard_True).Shape();
    }
    if (type == "revolve") {
        if (inputs.size() != 1) throw std::runtime_error(where + " requires exactly one profile input");
        require(p, "axis_origin", where + ".parameters"); require(p, "direction", where + ".parameters"); require(p, "angle_radians", where + ".parameters");
        const gp_Ax1 axis(point(p.at("axis_origin"), scale, where + ".axis_origin"), direction(p.at("direction"), where + ".direction"));
        const double angle = p.at("angle_radians").get<double>();
        if (!std::isfinite(angle) || std::abs(angle) < 1.0e-12 || std::abs(angle) > 2.0 * kPi + 1.0e-12) throw std::runtime_error(where + " angle_radians must be in (0, 2*pi]");
        return BRepPrimAPI_MakeRevol(inputs[0], axis, angle, Standard_True).Shape();
    }
    if (type == "loft") {
        if (inputs.size() < 2) throw std::runtime_error(where + " requires at least two wire inputs");
        BRepOffsetAPI_ThruSections loft(p.value("solid", true), Standard_False, 1.0e-7);
        for (const auto& input : inputs) {
            if (input.ShapeType() != TopAbs_WIRE) throw std::runtime_error(where + " accepts only wire inputs");
            loft.AddWire(TopoDS::Wire(input));
        }
        loft.Build();
        if (!loft.IsDone()) throw std::runtime_error(where + " loft construction failed");
        return loft.Shape();
    }
    if (type == "transform" || type == "mirror") {
        if (inputs.size() != 1) throw std::runtime_error(where + " requires exactly one input");
        gp_Trsf trsf;
        if (type == "mirror") {
            require(p, "center", where + ".parameters"); require(p, "normal", where + ".parameters");
            trsf.SetMirror(gp_Ax2(point(p.at("center"), scale, where + ".center"), direction(p.at("normal"), where + ".normal")));
        } else {
            if (p.contains("scale")) trsf.SetScale(gp_Pnt(0, 0, 0), p.at("scale").get<double>());
            if (p.contains("rotation")) {
                const auto r = vector3(p.at("rotation"), 1.0, where + ".rotation");
                gp_Trsf rx, ry, rz;
                rx.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), r[0]);
                ry.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), r[1]);
                rz.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), r[2]);
                gp_Trsf rotation = rz * ry * rx;
                trsf = rotation * trsf;
            }
            if (p.contains("translation")) {
                const auto t = vector3(p.at("translation"), scale, where + ".translation");
                gp_Trsf translate;
                translate.SetTranslation(gp_Vec(t[0], t[1], t[2]));
                trsf = translate * trsf;
            }
        }
        return BRepBuilderAPI_Transform(inputs[0], trsf, Standard_True).Shape();
    }
    if (type == "fuse" || type == "cut" || type == "common") {
        if (inputs.size() < 2) throw std::runtime_error(where + " requires at least two inputs");
        TopoDS_Shape accumulated = inputs[0];
        std::map<std::string, std::vector<TopoDS_Shape>> nameMap;
        for (const auto& named : model.entities().at(entity.inputs.front()).named_subshapes) nameMap[named.first] = named.second;
        for (std::size_t i = 1; i < inputs.size(); ++i) {
            std::unique_ptr<BRepAlgoAPI_BuilderAlgo> operationBuilder;
            if (type == "fuse") operationBuilder = std::make_unique<BRepAlgoAPI_Fuse>(accumulated, inputs[i]);
            else if (type == "cut") operationBuilder = std::make_unique<BRepAlgoAPI_Cut>(accumulated, inputs[i]);
            else operationBuilder = std::make_unique<BRepAlgoAPI_Common>(accumulated, inputs[i]);
            operationBuilder->Build();
            if (!operationBuilder->IsDone()) throw std::runtime_error(where + " OCCT boolean operation failed");
            const TopoDS_Shape result = checkedShape(operationBuilder->Shape(), entity.id);
            std::map<std::string, std::vector<TopoDS_Shape>> nextNames;
            auto mapHistory = [&](const std::map<std::string, std::vector<TopoDS_Shape>>& sourceNames) {
                for (const auto& named : sourceNames) {
                    for (const auto& source : named.second) {
                        if (operationBuilder->IsDeleted(source)) continue;
                        const auto& modified = operationBuilder->Modified(source);
                        if (!modified.IsEmpty()) {
                            for (TopTools_ListIteratorOfListOfShape it(modified); it.More(); it.Next()) nextNames[named.first].push_back(it.Value());
                        } else {
                            bool stillPresent = false;
                            for (TopExp_Explorer ex(result, source.ShapeType()); ex.More(); ex.Next()) {
                                if (ex.Current().IsSame(source)) { stillPresent = true; break; }
                            }
                            if (stillPresent) nextNames[named.first].push_back(source);
                        }
                        const auto& generated = operationBuilder->Generated(source);
                        for (TopTools_ListIteratorOfListOfShape it(generated); it.More(); it.Next()) nextNames[named.first].push_back(it.Value());
                    }
                }
            };
            mapHistory(nameMap);
            const auto& toolNames = model.entities().at(entity.inputs[i]).named_subshapes;
            mapHistory(toolNames);
            nameMap = std::move(nextNames);
            accumulated = result;
        }
        entity.named_subshapes = std::move(nameMap);
        return accumulated;
    }
    if (type == "split") {
        if (inputs.size() < 2) throw std::runtime_error(where + " requires at least two inputs (argument followed by tools)");
        BRepAlgoAPI_Splitter splitter;
        TopTools_ListOfShape arguments, tools;
        arguments.Append(inputs[0]);
        for (std::size_t i = 1; i < inputs.size(); ++i) tools.Append(inputs[i]);
        splitter.SetArguments(arguments);
        splitter.SetTools(tools);
        splitter.Build();
        if (!splitter.IsDone()) throw std::runtime_error(where + " OCCT splitter failed");
        return splitter.Shape();
    }
    if (type == "import") throw std::runtime_error(where + " STEP imports must be declared in cad.imports");
    throw std::runtime_error(where + " is not implemented");
}

void validateProjectKeys(const Json& project) {
    checkKeys(project, {"schema_version", "units", "cad", "fluid_domain", "background", "sizing", "boundary_layers", "quality", "resources", "output", "meshing", "optimization"}, "project");
    for (const char* key : {"schema_version", "units", "cad", "fluid_domain", "background", "sizing", "boundary_layers", "quality", "resources", "output"}) require(project, key, "project");
    if (project.at("schema_version") != 1 && project.at("schema_version") != 2)
        throw std::runtime_error("unsupported project schema_version (expected 1 or 2)");
    if(project.at("schema_version")==1 && project.contains("meshing"))
        throw std::runtime_error("meshing requires schema_version 2");
    if(project.at("schema_version")==2) require(project,"meshing","project");
    checkKeys(project.at("units"), {"length", "output_length"}, "units");
    checkKeys(project.at("cad"), {"imports", "operations", "groups", "repair_tolerance"}, "cad");
    if (project.at("cad").contains("repair_tolerance")) {
        throw std::runtime_error("cad.repair_tolerance is reserved but automatic shape repair is not enabled in M1");
    }
    for (const auto& item : project.at("cad").at("imports")) checkKeys(item, {"id", "file", "length_unit", "instance"}, "cad import");
    for (const auto& item : project.at("cad").at("operations")) {
        checkKeys(item, {"id", "operation", "inputs", "parameters"}, "cad operation");
        if (item.contains("parameters")) {
            checkKeys(item.at("parameters"), {"minimum", "maximum", "center", "normal", "start", "end", "axis_origin", "direction", "translation", "rotation", "radius", "major_radius", "minor_radius", "bottom_radius", "top_radius", "height", "angle_radians", "closed", "solid", "scale", "tolerance", "wire", "profile", "shape", "points"}, "CAD operation parameters");
        }
    }
    for (const auto& group : project.at("cad").at("groups")) {
        checkKeys(group, {"name", "kind", "selector"}, "physical group");
        if (group.contains("selector")) checkKeys(group.at("selector"), {"operation_id", "entity_type", "names", "center", "exclude_centers", "normal", "tolerance", "bounds", "adjacent_group"}, "group selector");
    }
    const Json& domain = project.at("fluid_domain");
    if (domain.value("mode", "") == "external") checkKeys(domain, {"mode", "outer_bounds", "solids", "seeds"}, "fluid_domain");
    else if (domain.value("mode", "") == "internal") checkKeys(domain, {"mode", "fluid_solids", "seeds", "caps"}, "fluid_domain");
    checkKeys(project.at("background"), {"bounds", "cell_counts", "cell_size", "uniform_cells"}, "background");
    checkKeys(project.at("background").at("bounds"), {"minimum", "maximum"}, "background.bounds");
    checkKeys(project.at("sizing"), {"minimum_size", "maximum_size", "surface_size", "volume_size", "surface_buffer_layers", "max_level", "balance_face_neighbors", "sources"}, "sizing");
    if(project.at("sizing").contains("surface_buffer_layers")) {
        const auto& layers=project.at("sizing").at("surface_buffer_layers");
        if(!layers.is_number_integer() || layers.get<int>()<0 || layers.get<int>()>8)
            throw std::runtime_error("sizing.surface_buffer_layers must be an integer in [0,8]");
    }
    for (const auto& source : project.at("sizing").at("sources")) {
        checkKeys(source, {"id", "kind", "target_size", "region", "transition_width", "apply_to"}, "sizing source");
        if (source.contains("region")) checkKeys(source.at("region"), {"bounds", "center", "axis", "radius", "length", "cad_group", "other_cad_group", "distance", "max_radius"}, "sizing source region");
        if (source.value("kind", std::string{}) == "curvature") {
            if (!source.contains("region") || !source.at("region").contains("max_radius") ||
                !source.at("region").at("max_radius").is_number() ||
                !std::isfinite(source.at("region").at("max_radius").get<double>()) ||
                source.at("region").at("max_radius").get<double>() <= 0.0) {
                throw std::runtime_error("curvature sizing source requires a positive region.max_radius; target_size remains the requested cell size");
            }
        }
    }
    checkKeys(project.at("boundary_layers"), {"enabled", "targets", "strict_coverage"}, "boundary_layers");
    for (const auto& target : project.at("boundary_layers").at("targets")) checkKeys(target, {"group", "first_layer", "count", "growth_ratio"}, "boundary-layer target");
    checkKeys(project.at("quality"), {"profile", "mode", "max_surface_deviation_m"}, "quality");
    if (project.at("quality").contains("max_surface_deviation_m")) {
        const auto& value = project.at("quality").at("max_surface_deviation_m");
        if (!value.is_number() || !std::isfinite(value.get<double>()) || value.get<double>() <= 0.0)
            throw std::runtime_error("quality.max_surface_deviation_m must be finite and positive (meters)");
    }
    checkKeys(project.at("resources"), {"max_cells", "max_memory_mb", "threads"}, "resources");
    checkKeys(project.at("output"), {"directory", "formats", "overwrite"}, "output");
}

Json shapeMetrics(const TopoDS_Shape& shape) {
    GProp_GProps volumeProps, surfaceProps;
    GProp_GProps linearProps;
    BRepGProp::VolumeProperties(shape, volumeProps);
    BRepGProp::SurfaceProperties(shape, surfaceProps);
    BRepGProp::LinearProperties(shape, linearProps);
    const auto b = boundsOf(shape);
    const gp_Pnt c = volumeProps.Mass() > 1.0e-30 ? volumeProps.CentreOfMass() : centerOfMass(shape);
    Json result;
    result["valid"] = BRepCheck_Analyzer(shape, Standard_True).IsValid();
    const auto shells = enumerate(shape, TopAbs_SHELL);
    bool closed = shape.ShapeType() == TopAbs_SHELL && BRep_Tool::IsClosed(shape);
    if (!shells.empty()) {
        closed = std::all_of(shells.begin(), shells.end(), [](const TopoDS_Shape& shell) {
            return BRep_Tool::IsClosed(TopoDS::Shell(shell));
        });
    }
    result["closed"] = closed;
    result["shape_type"] = static_cast<int>(shape.ShapeType());
    result["solids"] = enumerate(shape, TopAbs_SOLID).size();
    result["faces"] = enumerate(shape, TopAbs_FACE).size();
    result["edges"] = enumerate(shape, TopAbs_EDGE).size();
    result["vertices"] = enumerate(shape, TopAbs_VERTEX).size();
    result["volume_m3"] = enumerate(shape, TopAbs_SOLID).empty()
        ? Json(nullptr) : Json(std::abs(volumeProps.Mass()));
    result["surface_area_m2"] = surfaceProps.Mass();
    result["curve_length_m"] = enumerate(shape, TopAbs_EDGE).empty()
        ? Json(nullptr) : Json(linearProps.Mass());
    result["center_of_mass_m"] = {c.X(), c.Y(), c.Z()};
    result["bounds_m"] = {{b[0], b[1], b[2]}, {b[3], b[4], b[5]}};
    return result;
}

bool pointInSolid(const TopoDS_Shape& shape, const gp_Pnt& p) {
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) {
        BRepClass3d_SolidClassifier classifier(TopoDS::Solid(it.Current()), p, 1.0e-9);
        if (classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON) return true;
    }
    return false;
}

void validateSeeds(const TopoDS_Shape& domain, const Json& seeds, double scale, const std::string& label) {
    if (seeds.empty()) throw std::runtime_error(label + " needs at least one seed point");
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        if (!pointInSolid(domain, point(seeds[i], scale, label + ".seeds"))) {
            throw std::runtime_error(label + " seed " + std::to_string(i) + " is outside the selected fluid volume");
        }
    }
}

}  // namespace

CadModel buildCadModel(const std::filesystem::path& project_file) {
    Json project;
    try {
        std::ifstream input(project_file);
        if (!input) throw std::runtime_error("cannot open project file: " + project_file.string());
        input >> project;
    } catch (const nlohmann::json::exception& e) {
        throw std::runtime_error("invalid project JSON: " + std::string(e.what()));
    }
    validateProjectKeys(project);
    const std::string projectUnit = project.at("units").at("length").get<std::string>();
    const double scale = unitScale(projectUnit);
    if (project.at("units").at("output_length") != "m") throw std::runtime_error("output_length must be 'm' for BabelSim compatibility");

    CadModel model;
    std::set<std::string> ids;
    for (const auto& item : project.at("cad").at("imports")) {
        checkKeys(item, {"id", "file", "length_unit", "instance"}, "cad import");
        for (const char* key : {"id", "file", "length_unit"}) require(item, key, "cad import");
        const std::string id = item.at("id").get<std::string>();
        if (!ids.insert(id).second) throw std::runtime_error("duplicate CAD entity id '" + id + "'");
        const auto declaredUnit = item.at("length_unit").get<std::string>();
        (void)unitScale(declaredUnit);
        std::filesystem::path path(item.at("file").get<std::string>());
        if (path.is_relative()) path = project_file.parent_path() / path;
        const auto step = readStep(path, declaredUnit, item.value("instance", ""));
        CadEntity entity;
        entity.id = id;
        entity.operation = "step-import";
        entity.shape = step.shape;
        entity.surface_mesh = makeSurfaceMesh(entity.shape, id);
        entity.source_file = std::filesystem::weakly_canonical(path).string();
        entity.source_length_units = step.units;
        entity.named_subshapes = step.names;
        model.entities_.emplace(id, std::move(entity));
    }

    for (const auto& operation : project.at("cad").at("operations")) {
        require(operation, "id", "cad operation"); require(operation, "operation", "cad operation"); require(operation, "inputs", "cad operation");
        const std::string id = operation.at("id").get<std::string>();
        if (!ids.insert(id).second) throw std::runtime_error("duplicate CAD entity id '" + id + "'");
        CadEntity entity;
        entity.id = id;
        entity.operation = operation.at("operation").get<std::string>();
        entity.shape = checkedShape(makeOperation(operation, model, scale, entity), id);
        entity.surface_mesh = makeSurfaceMesh(entity.shape, id);
        model.entities_.emplace(id, std::move(entity));
    }

    std::set<std::string> groupNames;
    for (const auto& group : project.at("cad").at("groups")) {
        require(group, "name", "physical group"); require(group, "kind", "physical group"); require(group, "selector", "physical group");
        const std::string name = group.at("name").get<std::string>();
        if (!groupNames.insert(name).second) throw std::runtime_error("duplicate physical group name '" + name + "'");
        PhysicalGroup physical;
        physical.name = name;
        physical.kind = group.at("kind").get<std::string>();
        physical.selected_subshapes = selectGroup(model, group, scale);
        model.groups_.push_back(std::move(physical));
    }

    const Json& domain = project.at("fluid_domain");
    model.domain_mode_ = domain.value("mode", "");
    if (model.domain_mode_ == "external") {
        checkKeys(domain, {"mode", "outer_bounds", "solids", "seeds"}, "fluid_domain");
        const gp_Pnt lo = point(domain.at("outer_bounds").at("minimum"), scale, "fluid_domain.outer_bounds.minimum");
        const gp_Pnt hi = point(domain.at("outer_bounds").at("maximum"), scale, "fluid_domain.outer_bounds.maximum");
        if (hi.X() <= lo.X() || hi.Y() <= lo.Y() || hi.Z() <= lo.Z()) throw std::runtime_error("external fluid_domain has invalid outer_bounds");
        model.fluid_domain_ = BRepPrimAPI_MakeBox(lo, hi).Solid();
        for (const auto& solidId : domain.at("solids")) {
            const auto solidEntityId = solidId.get<std::string>();
            const auto found = model.entities_.find(solidEntityId);
            if (found == model.entities_.end()) throw std::runtime_error("external fluid_domain refers to unknown solid '" + solidEntityId + "'");
            if (enumerate(found->second.shape, TopAbs_SOLID).empty()) throw std::runtime_error("external fluid_domain entity '" + solidEntityId + "' contains no solid");
            BRepAlgoAPI_Cut cut(model.fluid_domain_, found->second.shape);
            cut.Build();
            if (!cut.IsDone()) throw std::runtime_error("failed to subtract CAD solid '" + solidEntityId + "' from external domain");
            model.fluid_domain_ = checkedShape(cut.Shape(), "fluid_domain");
        }
        validateSeeds(model.fluid_domain_, domain.at("seeds"), scale, "external fluid_domain");
    } else if (model.domain_mode_ == "internal") {
        checkKeys(domain, {"mode", "fluid_solids", "seeds", "caps"}, "fluid_domain");
        BRep_Builder builder;
        TopoDS_Compound compound;
        builder.MakeCompound(compound);
        std::vector<TopoDS_Shape> sourceShells;
        std::vector<TopoDS_Shape> sourceSolids;
        for (const auto& solidId : domain.at("fluid_solids")) {
            const auto solidEntityId = solidId.get<std::string>();
            const auto found = model.entities_.find(solidEntityId);
            if (found == model.entities_.end()) throw std::runtime_error("internal fluid_domain refers to unknown solid '" + solidEntityId + "'");
            const auto solids = enumerate(found->second.shape, TopAbs_SOLID);
            const auto shells = enumerate(found->second.shape, TopAbs_SHELL);
            if (!solids.empty()) sourceSolids.insert(sourceSolids.end(), solids.begin(), solids.end());
            else if (!shells.empty()) sourceShells.insert(sourceShells.end(), shells.begin(), shells.end());
            else throw std::runtime_error("internal fluid_domain entity '" + solidEntityId + "' contains neither solids nor closed surfaces");
        }
        if (!sourceSolids.empty() && (!sourceShells.empty() || !domain.at("caps").empty())) {
            throw std::runtime_error("internal fluid_domain cannot mix prebuilt solids with shell cap construction");
        }
        if (!sourceSolids.empty()) {
            for (const auto& solid : sourceSolids) builder.Add(compound, solid);
        } else {
            BRepBuilderAPI_Sewing sewing(1.0e-7, Standard_True, Standard_True, Standard_True, Standard_False);
            for (const auto& shell : sourceShells) {
                for (TopExp_Explorer face(shell, TopAbs_FACE); face.More(); face.Next()) sewing.Add(face.Current());
            }
            std::set<std::string> existingGroupNames;
            for (const auto& current : model.groups_) existingGroupNames.insert(current.name);
            for (const auto& cap : domain.at("caps")) {
                checkKeys(cap, {"name", "kind", "selector"}, "fluid_domain cap");
                require(cap, "name", "fluid_domain cap"); require(cap, "kind", "fluid_domain cap"); require(cap, "selector", "fluid_domain cap");
                const std::string capName = cap.at("name").get<std::string>();
                if (!existingGroupNames.insert(capName).second) throw std::runtime_error("duplicate physical group name '" + capName + "'");
                PhysicalGroup capGroup;
                capGroup.name = capName;
                capGroup.kind = cap.at("kind").get<std::string>();
                capGroup.selected_subshapes = selectGroup(model, cap, scale);
                for (const auto& face : capGroup.selected_subshapes) {
                    if (face.ShapeType() != TopAbs_FACE) throw std::runtime_error("internal-domain cap '" + capName + "' must select faces");
                    sewing.Add(face);
                }
                model.groups_.push_back(std::move(capGroup));
            }
            sewing.Perform();
            const TopoDS_Shape sewn = sewing.SewedShape();
            if (sewn.IsNull()) throw std::runtime_error("internal fluid_domain cap sewing produced no shell");
            const auto closedShells = enumerate(sewn, TopAbs_SHELL);
            if (closedShells.size() != 1 || !BRep_Tool::IsClosed(TopoDS::Shell(closedShells.front()))) {
                throw std::runtime_error("internal fluid_domain shell remains open after applying caps");
            }
            BRepBuilderAPI_MakeSolid solidMaker(TopoDS::Shell(closedShells.front()));
            if (!solidMaker.IsDone()) throw std::runtime_error("internal fluid_domain could not build a solid from capped shell");
            TopoDS_Solid closedSolid = solidMaker.Solid();
            if (!BRepLib::OrientClosedSolid(closedSolid)) throw std::runtime_error("internal fluid_domain could not orient capped shell");
            builder.Add(compound, checkedShape(closedSolid, "internal fluid_domain"));
        }
        model.fluid_domain_ = checkedShape(compound, "fluid_domain");
        validateSeeds(model.fluid_domain_, domain.at("seeds"), scale, "internal fluid_domain");
    } else {
        throw std::runtime_error("fluid_domain.mode must be 'external' or 'internal'");
    }
    model.fluid_surface_mesh_ = makeSurfaceMesh(model.fluid_domain_, "fluid_domain");
    return model;
}

CadModel makeFluidCadModel(const TopoDS_Shape& shape) {
    CadModel result;
    result.fluid_domain_ = checkedShape(shape, "derived-fluid-domain");
    result.fluid_surface_mesh_ = makeSurfaceMesh(shape, "derived-fluid-domain");
    result.domain_mode_ = "derived";
    return result;
}

CadModel inspectCadProject(const std::filesystem::path& project_file) {
    return buildCadModel(project_file);
}

Json geometryReport(const CadModel& model, const std::filesystem::path& project_file) {
    Json report;
    report["schema_version"] = 1;
    report["report_type"] = "cad-geometry";
    report["project_file"] = std::filesystem::absolute(project_file).string();
    report["internal_length_unit"] = "m";
    report["status"] = "passed";
    report["entities"] = Json::array();
    report["surface_entity_order"] = Json::array();
    std::size_t surfaceEntityIndex = 0;
    for (const auto& pair : model.entities()) {
        Json entity;
        entity["id"] = pair.first;
        entity["operation"] = pair.second.operation;
        entity["inputs"] = pair.second.inputs;
        entity["metrics"] = shapeMetrics(pair.second.shape);
        entity["surface_mesh"] = {
            {"vertices", pair.second.surface_mesh.vertices.size()},
            {"triangles", pair.second.surface_mesh.triangles.size()},
            {"degenerate_triangles_skipped", pair.second.surface_mesh.degenerate_triangles_skipped},
            {"bvh_nodes", pair.second.surface_mesh.bvh_nodes.size()},
            {"linear_deflection_m", pair.second.surface_mesh.linear_deflection_m}
        };
        entity["source_file"] = pair.second.source_file.empty() ? Json(nullptr) : Json(pair.second.source_file);
        entity["source_length_units"] = pair.second.source_length_units;
        entity["named_subshapes"] = Json::array();
        for (const auto& named : pair.second.named_subshapes) entity["named_subshapes"].push_back({{"name", named.first}, {"shape_count", named.second.size()}});
        report["entities"].push_back(std::move(entity));
        report["surface_entity_order"].push_back({{"cad_entity_id", surfaceEntityIndex++}, {"entity_id", pair.first}});
    }
    report["physical_groups"] = Json::array();
    for (const auto& group : model.groups()) {
        report["physical_groups"].push_back({{"name", group.name}, {"kind", group.kind}, {"selected_subshape_count", group.selected_subshapes.size()}});
    }
    report["fluid_domain"] = {{"mode", model.domainMode()}, {"metrics", shapeMetrics(model.fluidDomain())},
        {"surface_mesh", {{"vertices", model.fluidSurfaceMesh().vertices.size()},
                           {"triangles", model.fluidSurfaceMesh().triangles.size()},
                           {"degenerate_triangles_skipped", model.fluidSurfaceMesh().degenerate_triangles_skipped},
                           {"bvh_nodes", model.fluidSurfaceMesh().bvh_nodes.size()},
                           {"linear_deflection_m", model.fluidSurfaceMesh().linear_deflection_m}}}};
    return report;
}

}  // namespace babelmesh
