#include <BRepAlgoAPI_Cut.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <GProp_GProps.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Version.hxx>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    const TopoDS_Shape outer = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape();
    const TopoDS_Shape cutter = BRepPrimAPI_MakeBox(1.0, 2.0, 2.0).Shape();
    BRepAlgoAPI_Cut cut(outer, cutter);
    if (!cut.IsDone() || !BRepCheck_Analyzer(cut.Shape()).IsValid()) {
        throw std::runtime_error("OCCT box/boolean/validity probe failed");
    }

    GProp_GProps volume;
    BRepGProp::VolumeProperties(cut.Shape(), volume);
    if (std::abs(volume.Mass() - 4.0) > 1e-12) {
        throw std::runtime_error("OCCT boolean volume probe failed");
    }

    Eigen::Vector3d vector(1.0, 2.0, 3.0);
    const nlohmann::json document = {{"ok", true}, {"value", vector.sum()}};
    if (!document.at("ok").get<bool>() || document.at("value").get<double>() != 6.0) {
        throw std::runtime_error("Eigen/nlohmann dependency probe failed");
    }

    // Force compile/link use of STEP reader entry points as well as CAD modules.
    STEPControl_Reader step_reader;
    STEPCAFControl_Reader assembly_reader;
    (void)step_reader;
    (void)assembly_reader;

    std::cout << "dependency_probe: OCCT=" << OCC_VERSION_COMPLETE
              << " Eigen=" << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.'
              << EIGEN_MINOR_VERSION << " nlohmann_json=" << NLOHMANN_JSON_VERSION_MAJOR << '.'
              << NLOHMANN_JSON_VERSION_MINOR << '.' << NLOHMANN_JSON_VERSION_PATCH << '\n';
}
