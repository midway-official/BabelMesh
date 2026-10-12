#include "babelmesh/cad.hpp"

#include <BRepPrimAPI_MakeBox.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Writer.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Document.hxx>
#include <TopLoc_Location.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <gp_Trsf.hxx>

#include <TCollection_ExtendedString.hxx>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using Json = nlohmann::json;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& label) {
    if (std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << std::scientific << std::setprecision(12) << label << ": expected " << expected << ", got " << actual
                << " (tolerance " << tolerance << ')';
        throw std::runtime_error(message.str());
    }
}

Json baseProject(const std::string& unit = "m") {
    return Json{
        {"schema_version", 1},
        {"units", {{"length", unit}, {"output_length", "m"}}},
        {"cad", {{"imports", Json::array()}, {"operations", Json::array()}, {"groups", Json::array()}}},
        {"fluid_domain", {{"mode", "external"},
            {"outer_bounds", {{"minimum", {-10, -10, -10}}, {"maximum", {10, 10, 10}}}},
            {"solids", {"box24"}}, {"seeds", {{-5, -5, -5}}}}},
        {"background", {{"bounds", {{"minimum", {-10, -10, -10}}, {"maximum", {10, 10, 10}}}}, {"cell_counts", {1, 1, 1}}}},
        {"sizing", {{"minimum_size", 0.01}, {"maximum_size", 1.0}, {"max_level", 0}, {"sources", Json::array()}}},
        {"boundary_layers", {{"enabled", false}, {"targets", Json::array()}}},
        {"quality", {{"profile", "strict-cfd-v1"}, {"mode", "strict"}}},
        {"resources", {{"max_cells", 1000}, {"max_memory_mb", 256}, {"threads", 1}}},
        {"output", {{"directory", "out"}, {"formats", {"vtu"}}, {"overwrite", true}}}
    };
}

void addBox(Json& project, const std::string& id, Json lo, Json hi) {
    project["cad"]["operations"].push_back({
        {"id", id}, {"operation", "box"}, {"inputs", Json::array()},
        {"parameters", {{"minimum", std::move(lo)}, {"maximum", std::move(hi)}}}
    });
}

std::string addPlaneWire(Json& project, const std::string& prefix, const std::array<std::array<double, 3>, 4>& points) {
    Json edgeIds = Json::array();
    for (std::size_t i = 0; i < points.size(); ++i) {
        const std::string edgeId = prefix + "_edge_" + std::to_string(i);
        project["cad"]["operations"].push_back({{"id", edgeId}, {"operation", "line"}, {"inputs", Json::array()},
            {"parameters", {{"start", points[i]}, {"end", points[(i + 1U) % points.size()]}}}});
        edgeIds.push_back(edgeId);
    }
    const std::string wireId = prefix + "_wire";
    project["cad"]["operations"].push_back({{"id", wireId}, {"operation", "wire"}, {"inputs", edgeIds}, {"parameters", {{"closed", true}}}});
    return wireId;
}

void addPlaneFace(Json& project, const std::string& prefix, const std::array<std::array<double, 3>, 4>& points) {
    const std::string wireId = addPlaneWire(project, prefix, points);
    project["cad"]["operations"].push_back({{"id", prefix}, {"operation", "plane"}, {"inputs", {wireId}}});
}

std::filesystem::path writeProject(const std::filesystem::path& root, const std::string& name, const Json& project) {
    const auto path = root / name;
    std::ofstream out(path);
    out << project.dump(2) << '\n';
    if (!out) throw std::runtime_error("failed to write project fixture: " + path.string());
    return path;
}

void writeStepBox(const std::filesystem::path& path, const std::string& unit, double dimension) {
    STEPControl_Writer writer;
    const Standard_CString previous = Interface_Static::CVal("write.step.unit");
    const std::string previousUnit = previous == nullptr ? "MM" : previous;
    if (!Interface_Static::SetCVal("write.step.unit", unit.c_str())) {
        throw std::runtime_error("cannot configure STEP fixture writer unit '" + unit + "' (current='" + previousUnit + "')");
    }
    const auto box = BRepPrimAPI_MakeBox(dimension, dimension, dimension).Solid();
    if (writer.Transfer(box, STEPControl_AsIs) != IFSelect_RetDone || writer.Write(path.string().c_str()) != IFSelect_RetDone) {
        Interface_Static::SetCVal("write.step.unit", previousUnit.c_str());
        throw std::runtime_error("independent OCCT STEP writer failed to create fixture");
    }
    Interface_Static::SetCVal("write.step.unit", previousUnit.c_str());
}

void writeStepAssembly(const std::filesystem::path& path) {
    Handle(TDocStd_Document) document;
    XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", document);
    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(document->Main());
    const auto partShape = BRepPrimAPI_MakeBox(1000.0, 1000.0, 1000.0).Solid();
    const TDF_Label part = shapeTool->AddShape(partShape, Standard_False, Standard_False);
    TDataStd_Name::Set(part, TCollection_ExtendedString("CubePart", Standard_True));
    const TDF_Label assembly = shapeTool->NewShape();
    TDataStd_Name::Set(assembly, TCollection_ExtendedString("TwoCubes", Standard_True));
    const TDF_Label first = shapeTool->AddComponent(assembly, part, TopLoc_Location());
    TDataStd_Name::Set(first, TCollection_ExtendedString("CubeA", Standard_True));
    gp_Trsf translation;
    translation.SetTranslation(gp_Vec(2000.0, 0.0, 0.0));
    const TDF_Label second = shapeTool->AddComponent(assembly, part, TopLoc_Location(translation));
    TDataStd_Name::Set(second, TCollection_ExtendedString("CubeB", Standard_True));
    shapeTool->UpdateAssemblies();

    STEPCAFControl_Writer writer;
    writer.SetNameMode(Standard_True);
    if (!writer.Transfer(document, STEPControl_AsIs) || writer.Write(path.string().c_str()) != IFSelect_RetDone) {
        throw std::runtime_error("independent OCCT XCAF writer failed to create assembly fixture");
    }
}

Json entityMetrics(const babelmesh::CadModel& model, const std::string& id) {
    const auto report = babelmesh::geometryReport(model, "test-project.json");
    for (const auto& entity : report["entities"]) if (entity["id"] == id) return entity["metrics"];
    throw std::runtime_error("geometry report omits entity " + id);
}

double volume(const babelmesh::CadModel& model, const std::string& id) {
    return entityMetrics(model, id)["volume_m3"].get<double>();
}

void testPrimitivesBooleansAndGroups(const std::filesystem::path& root) {
    Json project = baseProject();
    addBox(project, "unit_a", {0, 0, 0}, {1, 1, 1});
    addBox(project, "unit_b", {0.5, 0, 0}, {1.5, 1, 1});
    addBox(project, "box24", {4, 5, 6}, {6, 8, 10});
    project["cad"]["operations"].push_back({{"id", "sewn_shell"}, {"operation", "shell"}, {"inputs", {"box24"}}});
    project["cad"]["operations"].push_back({{"id", "closed_solid"}, {"operation", "solid"}, {"inputs", {"sewn_shell"}}});
    project["cad"]["operations"].push_back({{"id", "common"}, {"operation", "common"}, {"inputs", {"unit_a", "unit_b"}}});
    project["cad"]["operations"].push_back({{"id", "union"}, {"operation", "fuse"}, {"inputs", {"unit_a", "unit_b"}}});
    project["cad"]["operations"].push_back({{"id", "difference"}, {"operation", "cut"}, {"inputs", {"unit_a", "unit_b"}}});
    project["cad"]["operations"].push_back({{"id", "cylinder"}, {"operation", "cylinder"}, {"inputs", Json::array()},
        {"parameters", {{"axis_origin", {0, 0, 0}}, {"direction", {0, 0, 1}}, {"radius", 1}, {"height", 2}}}});
    project["cad"]["operations"].push_back({{"id", "sphere"}, {"operation", "sphere"}, {"inputs", Json::array()},
        {"parameters", {{"center", {0, 0, 0}}, {"radius", 1}}}});
    project["cad"]["groups"].push_back({{"name", "box24_bottom"}, {"kind", "wall"},
        {"selector", {{"operation_id", "box24"}, {"entity_type", "face"}, {"center", {5, 6.5, 6}}, {"normal", {0, 0, -1}}, {"tolerance", 1.0e-6}}}});
    const auto projectPath = writeProject(root, "primitive-booleans.json", project);
    const auto model = babelmesh::buildCadModel(projectPath);
    near(volume(model, "box24"), 24.0, 1.0e-10, "2x3x4 box volume");
    near(volume(model, "closed_solid"), 24.0, 1.0e-10, "closed sewn shell solidification volume");
    near(volume(model, "common"), 0.5, 1.0e-9, "cube common volume");
    near(volume(model, "union"), 1.5, 1.0e-9, "cube union volume");
    near(volume(model, "difference"), 0.5, 1.0e-9, "cube cut volume");
    near(volume(model, "cylinder"), 2.0 * 3.14159265358979323846, 1.0e-8, "cylinder volume");
    near(volume(model, "sphere"), 4.0 * 3.14159265358979323846 / 3.0, 1.0e-8, "sphere volume");
    require(model.groups().size() == 1 && model.groups()[0].selected_subshapes.size() == 1,
            "normal+center selector should resolve exactly the bottom face");
    require(!model.fluidDomain().IsNull(), "external fluid domain must be built");
}

void testExtrudeAndLoft(const std::filesystem::path& root) {
    Json project = baseProject();
    addBox(project, "box24", {4, 5, 6}, {6, 8, 10});
    const std::array<std::array<double, 3>, 4> base{{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}}}};
    const std::array<std::array<double, 3>, 4> top{{{{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}}}};
    for (int ring = 0; ring < 2; ++ring) {
        const auto& points = ring == 0 ? base : top;
        for (int i = 0; i < 4; ++i) {
            project["cad"]["operations"].push_back({{"id", "p" + std::to_string(ring) + std::to_string(i)}, {"operation", "point"}, {"inputs", Json::array()},
                {"parameters", {{"center", points[static_cast<std::size_t>(i)]}}}});
        }
        for (int i = 0; i < 4; ++i) {
            project["cad"]["operations"].push_back({{"id", "e" + std::to_string(ring) + std::to_string(i)}, {"operation", "line"},
                {"inputs", Json::array()}, {"parameters", {{"start", points[static_cast<std::size_t>(i)]}, {"end", points[static_cast<std::size_t>((i + 1) % 4)]}}}});
        }
        Json edgeInputs = Json::array();
        for (int i = 0; i < 4; ++i) edgeInputs.push_back("e" + std::to_string(ring) + std::to_string(i));
        project["cad"]["operations"].push_back({{"id", "wire" + std::to_string(ring)}, {"operation", "wire"}, {"inputs", edgeInputs}, {"parameters", {{"closed", true}}}});
    }
    project["cad"]["operations"].push_back({{"id", "face"}, {"operation", "plane"}, {"inputs", {"wire0"}}});
    project["cad"]["operations"].push_back({{"id", "extrusion"}, {"operation", "extrude"}, {"inputs", {"face"}}, {"parameters", {{"direction", {0, 0, 1}}}}});
    project["cad"]["operations"].push_back({{"id", "loft"}, {"operation", "loft"}, {"inputs", {"wire0", "wire1"}}, {"parameters", {{"solid", true}}}});
    const auto model = babelmesh::buildCadModel(writeProject(root, "extrude-loft.json", project));
    near(volume(model, "extrusion"), 1.0, 1.0e-9, "unit-square extrusion volume");
    near(volume(model, "loft"), 1.0, 1.0e-8, "unit-square loft volume");
}

void testAdvancedCadOperations(const std::filesystem::path& root) {
    Json project = baseProject();
    addBox(project, "box24", {4, 5, 6}, {6, 8, 10});
    project["cad"]["operations"].push_back({{"id", "arc_start"}, {"operation", "point"}, {"inputs", Json::array()}, {"parameters", {{"center", {1, 0, 0}}}}});
    project["cad"]["operations"].push_back({{"id", "arc_mid"}, {"operation", "point"}, {"inputs", Json::array()}, {"parameters", {{"center", {std::sqrt(0.5), std::sqrt(0.5), 0}}}}});
    project["cad"]["operations"].push_back({{"id", "arc_end"}, {"operation", "point"}, {"inputs", Json::array()}, {"parameters", {{"center", {0, 1, 0}}}}});
    project["cad"]["operations"].push_back({{"id", "arc"}, {"operation", "arc"}, {"inputs", {"arc_start", "arc_mid", "arc_end"}}});
    project["cad"]["operations"].push_back({{"id", "circle"}, {"operation", "circle"}, {"inputs", Json::array()},
        {"parameters", {{"center", {0, 0, 0}}, {"normal", {0, 0, 1}}, {"radius", 1}}}});
    project["cad"]["operations"].push_back({{"id", "spline"}, {"operation", "spline"}, {"inputs", Json::array()},
        {"parameters", {{"points", {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}}}}}});
    project["cad"]["operations"].push_back({{"id", "cone"}, {"operation", "cone"}, {"inputs", Json::array()},
        {"parameters", {{"axis_origin", {0, 0, 0}}, {"direction", {0, 0, 1}}, {"bottom_radius", 1}, {"top_radius", 0}, {"height", 3}}}});
    addBox(project, "unit_box", {0, 0, 0}, {1, 1, 1});
    project["cad"]["operations"].push_back({{"id", "transformed"}, {"operation", "transform"}, {"inputs", {"unit_box"}},
        {"parameters", {{"translation", {3, 0, 0}}, {"scale", 2}}}});
    project["cad"]["operations"].push_back({{"id", "mirrored"}, {"operation", "mirror"}, {"inputs", {"unit_box"}},
        {"parameters", {{"center", {0, 0, 0}}, {"normal", {1, 0, 0}}}}});

    addPlaneFace(project, "revolve_profile", {{{{1, 0, 0}}, {{2, 0, 0}}, {{2, 0, 1}}, {{1, 0, 1}}}});
    project["cad"]["operations"].push_back({{"id", "revolved"}, {"operation", "revolve"}, {"inputs", {"revolve_profile"}},
        {"parameters", {{"axis_origin", {0, 0, 0}}, {"direction", {0, 0, 1}}, {"angle_radians", 2.0 * 3.14159265358979323846}}}});

    addPlaneFace(project, "split_box", {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}}}});
    project["cad"]["operations"].push_back({{"id", "split_box_solid"}, {"operation", "extrude"}, {"inputs", {"split_box"}}, {"parameters", {{"direction", {0, 0, 1}}}}});
    addPlaneFace(project, "split_plane", {{{{0.5, 0, 0}}, {{0.5, 1, 0}}, {{0.5, 1, 1}}, {{0.5, 0, 1}}}});
    project["cad"]["operations"].push_back({{"id", "split_result"}, {"operation", "split"}, {"inputs", {"split_box_solid", "split_plane"}}});

    const auto outerWire = addPlaneWire(project, "hole_outer", {{{{0, 0, 0}}, {{4, 0, 0}}, {{4, 4, 0}}, {{0, 4, 0}}}});
    const auto innerWire = addPlaneWire(project, "hole_inner", {{{{1, 1, 0}}, {{3, 1, 0}}, {{3, 3, 0}}, {{1, 3, 0}}}});
    project["cad"]["operations"].push_back({{"id", "plane_with_hole"}, {"operation", "plane"}, {"inputs", {outerWire, innerWire}}});

    const auto model = babelmesh::buildCadModel(writeProject(root, "advanced-operations.json", project));
    near(entityMetrics(model, "arc")["curve_length_m"].get<double>(), 0.5 * 3.14159265358979323846, 1.0e-7, "quarter circle arc length");
    near(entityMetrics(model, "circle")["curve_length_m"].get<double>(), 2.0 * 3.14159265358979323846, 1.0e-7, "circle length");
    require(entityMetrics(model, "spline")["curve_length_m"].get<double>() > 2.0, "spline must produce a measurable curve");
    near(volume(model, "cone"), 3.14159265358979323846, 1.0e-8, "cone volume");
    near(volume(model, "transformed"), 8.0, 1.0e-8, "scaled/transformed solid volume");
    near(volume(model, "mirrored"), 1.0, 1.0e-8, "mirrored solid volume");
    near(volume(model, "revolved"), 3.0 * 3.14159265358979323846, 1.0e-7, "full revolve analytic volume");
    near(volume(model, "split_result"), 1.0, 1.0e-8, "split total volume");
    require(entityMetrics(model, "split_result")["solids"] == 2, "split should expose two solid pieces");
    near(entityMetrics(model, "plane_with_hole")["surface_area_m2"].get<double>(), 12.0, 1.0e-9, "planar face with hole area");
}

void testMillimeterProjectUnits(const std::filesystem::path& root) {
    Json project = baseProject("mm");
    project["fluid_domain"]["outer_bounds"] = {{"minimum", {-10000, -10000, -10000}}, {"maximum", {10000, 10000, 10000}}};
    project["fluid_domain"]["seeds"] = {{-5000, -5000, -5000}};
    project["background"]["bounds"] = project["fluid_domain"]["outer_bounds"];
    addBox(project, "box24", {0, 0, 0}, {2000, 3000, 4000});
    const auto model = babelmesh::buildCadModel(writeProject(root, "millimeters.json", project));
    near(volume(model, "box24"), 24.0, 1.0e-9, "mm project unit normalization");
}

void testStepImportUnitsAndAssemblies(const std::filesystem::path& root, const std::filesystem::path& repositoryRoot) {
    const auto mmFile = repositoryRoot / "tests/fixtures/gmsh-box-1000mm.step";
    const auto mFile = root / "fixture-m.step";
    const auto assemblyFile = root / "fixture-assembly.step";
    writeStepBox(mFile, "M", 1000.0);
    writeStepAssembly(assemblyFile);
    for (const auto& test : {std::pair<std::filesystem::path, std::string>{mmFile, "mm"}, {mFile, "m"}}) {
        Json project = baseProject();
        project["cad"]["imports"].push_back({{"id", "imported"}, {"file", test.first.string()}, {"length_unit", test.second}});
        project["fluid_domain"]["solids"] = {"imported"};
        const auto model = babelmesh::buildCadModel(writeProject(root, test.second + "-import.json", project));
        near(volume(model, "imported"), 1.0, 1.0e-8, "STEP import normalized to meters");
        require(!model.entities().at("imported").source_length_units.empty(), "STEP unit metadata should be recorded");
    }

    Json assemblyProject = baseProject();
    assemblyProject["cad"]["imports"].push_back({{"id", "assembly"}, {"file", assemblyFile.string()}, {"length_unit", "mm"}});
    assemblyProject["fluid_domain"]["solids"] = {"assembly"};
    const auto fullAssembly = babelmesh::buildCadModel(writeProject(root, "assembly-import.json", assemblyProject));
    near(volume(fullAssembly, "assembly"), 2.0, 1.0e-8, "two-instance STEP assembly volume");
    const auto fullBounds = babelmesh::geometryReport(fullAssembly, "assembly-import.json")["entities"][0]["metrics"]["bounds_m"];
    near(fullBounds[0][0].get<double>(), 0.0, 1.0e-7, "assembly first instance lower x");
    near(fullBounds[1][0].get<double>(), 3.0, 1.0e-7, "assembly second instance upper x");
    require(fullAssembly.entities().at("assembly").named_subshapes.count("TwoCubes/CubeB") == 1,
            "STEP component name path should be retained");

    assemblyProject["cad"]["imports"][0]["instance"] = "TwoCubes/CubeB";
    const auto selectedAssemblyInstance = babelmesh::buildCadModel(writeProject(root, "assembly-instance-import.json", assemblyProject));
    near(volume(selectedAssemblyInstance, "assembly"), 1.0, 1.0e-8, "selected STEP assembly instance volume");
    const auto instanceBounds = babelmesh::geometryReport(selectedAssemblyInstance, "assembly-instance-import.json")["entities"][0]["metrics"]["bounds_m"];
    near(instanceBounds[0][0].get<double>(), 2.0, 2.0e-7, "selected assembly instance lower x");
    near(instanceBounds[1][0].get<double>(), 3.0, 2.0e-7, "selected assembly instance upper x");

    Json booleanTags = baseProject();
    booleanTags["cad"]["imports"].push_back({{"id", "assembly"}, {"file", assemblyFile.string()}, {"length_unit", "mm"}});
    addBox(booleanTags, "overlap", {0.5, 0, 0}, {1.5, 1, 1});
    booleanTags["cad"]["operations"].push_back({{"id", "fused"}, {"operation", "fuse"}, {"inputs", {"assembly", "overlap"}}});
    booleanTags["cad"]["groups"].push_back({{"name", "cube_b_tag"}, {"kind", "wall"},
        {"selector", {{"operation_id", "fused"}, {"entity_type", "solid"}, {"names", {"TwoCubes/CubeB"}}}}});
    booleanTags["fluid_domain"]["solids"] = {"fused"};
    const auto taggedBoolean = babelmesh::buildCadModel(writeProject(root, "assembly-boolean-tags.json", booleanTags));
    near(volume(taggedBoolean, "fused"), 2.5, 1.0e-8, "Boolean history tagged assembly volume");
    require(taggedBoolean.groups().size() == 1 && taggedBoolean.groups()[0].selected_subshapes.size() == 1,
            "boolean operation must preserve the unmodified assembly instance tag");
}

void testInvalidSeedAndDuplicateId(const std::filesystem::path& root) {
    Json invalidSeed = baseProject();
    addBox(invalidSeed, "box24", {4, 5, 6}, {6, 8, 10});
    invalidSeed["fluid_domain"]["seeds"] = {{20, 0, 0}};
    bool seedRejected = false;
    try { (void)babelmesh::buildCadModel(writeProject(root, "bad-seed.json", invalidSeed)); }
    catch (const std::exception& e) { seedRejected = std::string(e.what()).find("seed 0") != std::string::npos; }
    require(seedRejected, "outside-fluid seed must fail with indexed diagnostic");

    Json duplicate = baseProject();
    addBox(duplicate, "box24", {4, 5, 6}, {6, 8, 10});
    duplicate["cad"]["operations"].push_back(duplicate["cad"]["operations"][0]);
    bool duplicateRejected = false;
    try { (void)babelmesh::buildCadModel(writeProject(root, "duplicate-id.json", duplicate)); }
    catch (const std::exception& e) { duplicateRejected = std::string(e.what()).find("duplicate CAD entity id") != std::string::npos; }
    require(duplicateRejected, "duplicate CAD IDs must fail explicitly");

    Json internal = baseProject();
    addBox(internal, "box24", {4, 5, 6}, {6, 8, 10});
    internal["fluid_domain"] = {{"mode", "internal"}, {"fluid_solids", {"box24"}}, {"seeds", {{5, 6, 7}}}, {"caps", Json::array()}};
    const auto internalDomain = babelmesh::buildCadModel(writeProject(root, "internal-domain.json", internal));
    near(std::abs(babelmesh::geometryReport(internalDomain, "internal-domain.json")["fluid_domain"]["metrics"]["volume_m3"].get<double>()),
         24.0, 1.0e-9, "internal domain volume and seed");

    Json capped = baseProject();
    addBox(capped, "box24", {4, 5, 6}, {6, 8, 10});
    addPlaneFace(capped, "bottom", {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}}}});
    addPlaneFace(capped, "side_x0", {{{{0, 0, 0}}, {{0, 1, 0}}, {{0, 1, 1}}, {{0, 0, 1}}}});
    addPlaneFace(capped, "side_x1", {{{{1, 0, 0}}, {{1, 0, 1}}, {{1, 1, 1}}, {{1, 1, 0}}}});
    addPlaneFace(capped, "side_y0", {{{{0, 0, 0}}, {{0, 0, 1}}, {{1, 0, 1}}, {{1, 0, 0}}}});
    addPlaneFace(capped, "side_y1", {{{{0, 1, 0}}, {{1, 1, 0}}, {{1, 1, 1}}, {{0, 1, 1}}}});
    addPlaneFace(capped, "top_cap", {{{{0, 0, 1}}, {{0, 1, 1}}, {{1, 1, 1}}, {{1, 0, 1}}}});
    Json shellFaces = Json::array({"bottom", "side_x0", "side_x1", "side_y0", "side_y1"});
    capped["cad"]["operations"].push_back({{"id", "open_shell"}, {"operation", "shell"}, {"inputs", shellFaces}});
    capped["fluid_domain"] = {{"mode", "internal"}, {"fluid_solids", {"open_shell"}}, {"seeds", {{0.5, 0.5, 0.5}}},
        {"caps", {{{"name", "outlet_cap"}, {"kind", "outlet"}, {"selector", {{"operation_id", "top_cap"}, {"entity_type", "face"}}}}}}};
    const auto cappedModel = babelmesh::buildCadModel(writeProject(root, "internal-capped-domain.json", capped));
    const auto shellVolume = entityMetrics(cappedModel, "open_shell")["volume_m3"];
    require(shellVolume.is_null(), "open shell entity must not report a volume");
    near(babelmesh::geometryReport(cappedModel, "internal-capped-domain.json")["fluid_domain"]["metrics"]["volume_m3"].get<double>(),
         1.0, 1.0e-9, "internal cap should close unit-cube fluid domain");
    require(cappedModel.groups().size() == 1 && cappedModel.groups()[0].name == "outlet_cap",
            "internal cap group should become a named boundary group");

    Json openShell = capped;
    openShell["fluid_domain"]["caps"] = Json::array();
    bool openRejected = false;
    try { (void)babelmesh::buildCadModel(writeProject(root, "internal-open-domain.json", openShell)); }
    catch (const std::exception& e) { openRejected = std::string(e.what()).find("remains open") != std::string::npos; }
    require(openRejected, "internal open shell without caps must be rejected");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: m1_cad_test REPOSITORY_ROOT");
        const auto repositoryRoot = std::filesystem::absolute(argv[1]);
        const auto root = std::filesystem::temp_directory_path() / "babelmesh-m1-cad-test";
        std::filesystem::remove_all(root);
        std::filesystem::create_directories(root);
        testPrimitivesBooleansAndGroups(root);
        testExtrudeAndLoft(root);
        testAdvancedCadOperations(root);
        testMillimeterProjectUnits(root);
        testStepImportUnitsAndAssemblies(root, repositoryRoot);
        testInvalidSeedAndDuplicateId(root);
        std::cout << "M1 CAD acceptance passed: primitives, booleans, curves, faces with holes, transforms, extrude/revolve/loft/split, m-mm normalization, STEP units and assemblies, internal domains and invalid-input rejection\n";
        std::filesystem::remove_all(root);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "M1 CAD acceptance failed: " << e.what() << '\n';
        return 1;
    }
}
