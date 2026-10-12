#!/usr/bin/env python3
"""Exercise real CAD, generation, solver reader and mesh exchange paths.

The suite distinguishes topology, sampled CAD fit, face quality and commercial
application import. It never promotes diagnostic meshes to CFD acceptance.
"""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Suite:
    def __init__(self, cli, reader, root, output):
        self.cli, self.reader, self.root, self.output = cli, reader, root, output
        self.number, self.checks, self.results = 0, [], {}
        self.start = time.monotonic()

    def run(self, name, arguments, expected=True, environment=None):
        self.number += 1
        log = self.output / "logs" / f"{self.number:03d}-{name}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run([str(a) for a in arguments], cwd=self.root,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=180, env=environment)
        log.write_text(result.stdout)
        self.checks.append({"name": name, "return_code": result.returncode,
                            "expected_success": expected,
                            "log": str(log.relative_to(self.output))})
        require((result.returncode == 0) == expected,
                f"{name}: unexpected status {result.returncode}:\n{result.stdout[-3500:]}")
        print(f"{name}: {'passed' if expected else 'rejected as expected'}", flush=True)
        return result.stdout

    def generate(self, name, project):
        project = copy.deepcopy(project)
        project["output"]["directory"] = str(self.output / name)
        project["output"]["overwrite"] = False
        path = self.output / "projects" / f"{name}.json"
        save(path, project)
        log = self.run(name, [self.cli, "generate", path, "--out", self.output / name,
                              "--overwrite"])
        require(log.index("[1/7] Complete Cartesian") < log.index("[3/7] Carve") <
                log.index("[4/7] Commit"), "background/carve/fill order is wrong")
        report, mesh = load(self.output / name / "mesh_report.json"), load(self.output / name / "mesh.json")
        require(report["quality"]["negative_or_zero_volume_cells"] == 0, "invalid volume")
        require(report["metadata"]["cartesian_core_nodes_moved"] == 0, "background nodes moved")
        locked = mesh["metadata"]["domain_bounds_constraint"]
        require(report["metadata"]["bounds_audit"]["passed"] and
                report["metadata"]["bounds_audit"]["maximum_excess_m"] == 0,
                "fixed bounds audit failed")
        require(all(locked["minimum_m"][a] <= point[a] <= locked["maximum_m"][a]
                    for point in mesh["points_m"] for a in range(3)),
                "generation/optimization expanded the declared domain")
        require(report["metadata"]["geometry_distance_audit"]["passed"], "CAD fit did not pass")
        require(report["metadata"]["optimization"]["added_cells"] == 0, "optimizer added cells")
        opt = report["metadata"]["optimization"]
        for metric in ["max_non_orthogonality_degrees", "p99_non_orthogonality_degrees",
                       "max_internal_skewness", "p99_internal_skewness"]:
            require(opt["after"][metric] <= opt["before"][metric] + 1e-8,
                    "optimization worsened " + metric)
        carve = report["metadata"]["background_carving"]
        require(sum(carve[k] for k in ["removed_nonfluid_cells", "reserved_fit_and_layer_cells",
                                      "reserved_refinement_transition_cells", "removed_nonmanifold_core_cells",
                                      "retained_cartesian_cells"]) == carve["candidate_cells"],
                "carving accounting is incomplete")
        if project["meshing"]["dimension"] == 2:
            require(report["cell_types"]["tetrahedron"] == report["cell_types"]["pyramid"] == 0,
                    "2D mode produced tetrahedra or pyramids")
            a = mesh["extrusion_axis"]
            levels = sorted(set(p[a] for p in mesh["points_m"]))
            require(len(levels) == 2 and math.isclose(levels[1] - levels[0],
                    mesh["extrusion_thickness_m"], rel_tol=1e-10), "not a single extrusion layer")
            kinds = {p["name"]: p["kind"] for p in mesh["patches"]}
            require(kinds.get("front") == kinds.get("back") == "symmetry", "wrong planar caps")
        self.probe(name, self.output / name / "mesh.bsm", report)
        self.results[name] = {"cells": report["cells"], "cell_types": report["cell_types"],
                              "regions": report["regions"], "quality": report["quality"],
                              "cad_fit": report["metadata"]["geometry_distance_audit"],
                              "boundary_layers": report["metadata"].get("boundary_layers", []),
                              "optimization": report["metadata"]["optimization"]}
        return mesh, report, path

    def probe(self, name, path, report):
        value = json.loads(self.run(name + "-BabelSim-reader", [self.reader, path]))
        save(self.output / "reader-reports" / (name + ".json"), value)
        require(value["cells"] == report["cells"], "solver cell count differs")
        metrics = [("max_non_orthogonality_degrees", value["all"]["nonorth"]["max"]),
                   ("p99_non_orthogonality_degrees", value["all"]["nonorth"]["p99"])]
        if value["internal"]["skew"]:
            metrics += [("max_internal_skewness", value["internal"]["skew"]["max"]),
                        ("p99_internal_skewness", value["internal"]["skew"]["p99"])]
        for key, expected in metrics:
            require(math.isclose(report["quality"][key], expected, abs_tol=1e-8),
                    f"producer/solver geometry differs for {name}: {key}")

    def exchange(self, name, mesh, report, formats):
        extensions = {"fluent-ascii": "msh", "cgns": "cgns", "cgns-poly": "cgns",
                       "babelsim-v3": "bsm"}
        for fmt in formats:
            stem = name + "-" + fmt
            file = self.output / "exchange" / (stem + "." + extensions[fmt])
            checkpoint = self.output / "exchange" / (stem + "-roundtrip.json")
            self.run(stem + "-write", [self.cli, "convert", self.output / name / "mesh.json",
                                       "--to", fmt, "--out", file, "--overwrite"])
            self.run(stem + "-read", [self.cli, "convert", file, "--to", "checkpoint",
                                      "--out", checkpoint, "--overwrite"])
            other, converted = load(checkpoint), load(str(checkpoint) + ".report.json")
            require(len(other["cells"]) == len(mesh["cells"]) and
                    len(other["faces"]) == len(mesh["faces"]), "conversion changed connectivity counts")
            require(converted["cell_types"] == report["cell_types"], "conversion changed actual cell types")
            require(math.isclose(converted["volume_m3"], report["volume_m3"],
                                 rel_tol=1e-11), "conversion changed volume/orientation")
            require({p["name"]: p["kind"] for p in other["patches"]} ==
                    {p["name"]: p["kind"] for p in mesh["patches"]}, "conversion changed boundary groups")
            require(other["generation_dimension"] == mesh["generation_dimension"],
                    "exchange metadata/inference lost generation rule")
            if fmt != "babelsim-v3":
                require([c["region"] for c in other["cells"]] ==
                        [c["region"] for c in mesh["cells"]], "exchange lost cell region labels")
                require(other["metadata"]["domain_bounds_constraint"] ==
                        mesh["metadata"]["domain_bounds_constraint"],
                        "exchange lost the fixed generation envelope")
            bsm = self.output / "exchange" / (stem + "-roundtrip.bsm")
            self.run(stem + "-solver-export", [self.cli, "convert", checkpoint, "--to", "babelsim-v3",
                                               "--out", bsm, "--overwrite"])
            self.probe(stem, bsm, converted)
            checker = self.root / ".cache/dev-packages/cgns-tools/usr/bin/cgnscheck"
            if fmt.startswith("cgns") and checker.exists():
                env = os.environ.copy()
                env["LD_LIBRARY_PATH"] = str(self.root / ".cache/dev-packages/cgns/usr/lib/x86_64-linux-gnu") + ":" + env.get("LD_LIBRARY_PATH", "")
                checked = self.run(stem + "-CGNS-independent-check", [checker, file], environment=env)
                require("ERROR:" not in checked and "WARNING:" not in checked,
                        "CGNS independent validator reported errors or warnings")

    def negative(self, name, project, reason):
        path, directory = self.output / "projects" / (name + ".json"), self.output / name
        save(path, project)
        text = self.run(name, [self.cli, "generate", path, "--out", directory], expected=False)
        require(reason in text, f"{name}: expected diagnostic {reason!r}")
        require(not (directory / "mesh.bsm").exists(), "failed generation published a mesh")

    def complete(self):
        report = {"status": "passed", "elapsed_seconds": time.monotonic() - self.start,
                  "commands": len(self.checks), "checks": self.checks, "cases": self.results,
                  "binaries": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [self.cli, self.reader]},
                  "commercial_application_import_tested": False,
                  "full_cfd_profile_evaluated": False,
                  "global_cell_overlap_audit_evaluated": False}
        save(self.output / "integration_report.json", report)
        print(f"{len(self.checks)} real command checks passed; evidence: {self.output}", flush=True)


def internal(project, axis):
    p = copy.deepcopy(project)
    p["cad"]["operations"] = [{"id": "channel", "operation": "box", "inputs": [],
                                 "parameters": {"minimum": [0, 0, 0], "maximum": [1, 1, 1]}}]
    p["cad"]["groups"] = [{"name": "channel_wall", "kind": "wall",
                            "selector": {"operation_id": "channel", "entity_type": "face"}}]
    p["fluid_domain"] = {"mode": "internal", "fluid_solids": ["channel"], "caps": [], "seeds": [[.5, .5, .5]]}
    p["background"]["bounds"] = {"minimum": [0, 0, 0], "maximum": [1, 1, 1]}
    p["background"]["cell_counts"] = [8, 8, 8]
    p["background"]["cell_counts"]["xyz".index(axis)] = 1
    p["meshing"].update({"section_axis": axis, "section_coordinate": .5, "thickness": .1})
    p["sizing"].update({"surface_size": .125, "volume_size": .125, "maximum_size": .125})
    p["boundary_layers"] = {"enabled": False, "targets": [], "strict_coverage": True}
    return p


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cli", type=Path)
    parser.add_argument("reader", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    temporary = tempfile.TemporaryDirectory(prefix="babelmesh-hybrid-suite-") if args.out is None else None
    output = (args.out or Path(temporary.name)).resolve()
    output.mkdir(parents=True, exist_ok=True)
    suite = Suite(args.cli.resolve(), args.reader.resolve(), args.root.resolve(), output)
    circle, sphere = load(suite.root / "examples/cylinder2d.json"), load(suite.root / "examples/sphere3d.json")
    circle["quality"]["mode"] = "strict"
    cm, cr, path = suite.generate("cylinder2d", circle)
    require(cr["quality"]["face_quality_thresholds_passed"], "default 2D face quality failed")
    for layer in cr["metadata"]["boundary_layers"][0]["height_measurements"]:
        require(layer["p95_relative_height_error"] < 1e-8, "2D actual layer height differs")
    digest = hashlib.sha256((output / "cylinder2d/mesh.bsm").read_bytes()).hexdigest()
    suite.run("overwrite-disabled", [suite.cli, "generate", path, "--out", output / "cylinder2d"], expected=False)
    require(digest == hashlib.sha256((output / "cylinder2d/mesh.bsm").read_bytes()).hexdigest(), "existing output changed")

    spacing = copy.deepcopy(circle)
    spacing["background"].pop("cell_counts")
    spacing["background"]["cell_size"] = [.2, .2, 1]
    _, spacing_report, _ = suite.generate("cylinder2d-spacing", spacing)
    require(spacing_report["cells"] == cr["cells"] and
            spacing_report["cell_types"] == cr["cell_types"], "spacing/count background mismatch")

    adaptive = copy.deepcopy(circle)
    adaptive["sizing"]["sources"] = [{"id": "wake_refinement", "kind": "box", "apply_to": "background",
                                      "target_size": .1, "transition_width": .1,
                                      "region": {"bounds": {"minimum": [-1.8, -1.8, -.5], "maximum": [-.9, -.9, .5]}}}]
    am, ar, _ = suite.generate("cylinder2d-adaptive", adaptive)
    require(ar["cell_types"]["polyhedron"] > 0, "2D hanging-node polygonal prisms missing")
    for axis in ["x", "y"]:
        m, r, _ = suite.generate("channel2d-" + axis, internal(circle, axis))
        require(math.isclose(r["volume_m3"], .1, rel_tol=1e-10), "internal section measure is wrong")

    boolean = internal(circle, "z")
    boolean["cad"]["operations"] += [
        {"id": "hole", "operation": "cylinder", "inputs": [], "parameters": {
            "axis_origin": [.5, .5, 0], "direction": [0, 0, 1], "radius": .15, "height": 1}},
        {"id": "fluid", "operation": "cut", "inputs": ["channel", "hole"], "parameters": {}}]
    boolean["fluid_domain"]["fluid_solids"] = ["fluid"]
    boolean["fluid_domain"]["seeds"] = [[.1, .1, .5]]
    boolean["cad"]["groups"] = []
    boolean["quality"]["mode"] = "diagnostic"
    _, br, _ = suite.generate("boolean-channel2d", boolean)
    require(abs(br["volume_m3"] / (.1 * (1 - math.pi * .15**2)) - 1) < .02, "boolean volume mismatch")

    sm, sr, _ = suite.generate("sphere3d", sphere)
    require(sr["metadata"]["retained_cartesian_outer_boundary"],
            "aligned outer Cartesian faces were unnecessarily carved")
    require(sr["cells"] < 8000, "unnecessary outer tetrahedral band returned")
    require(sr["cell_types"]["hexahedron"] > 0 and sr["cell_types"]["pyramid"] > 0 and
            sr["cell_types"]["tetrahedron"] > 0 and sr["cell_types"]["prism"] == 0, "3D hybrid types missing")
    layers = copy.deepcopy(sphere)
    layers["boundary_layers"] = {"enabled": True, "strict_coverage": True,
                                  "targets": [{"group": "body_wall", "first_layer": .035,
                                               "count": 2, "growth_ratio": 1.2}]}
    lm, lr, _ = suite.generate("sphere3d-layers", layers)
    require(all(c["type"] == 10 for c in lm["cells"] if c["region"] == 1), "3D wall layers must be tetrahedra")
    require(lr["regions"]["boundary_layer"] > 0, "actual wall layers missing")
    for layer in lr["metadata"]["boundary_layers"][0]["height_measurements"]:
        require(layer["p95_relative_height_error"] < .05, "3D normal height deviates too far")

    adaptive3 = copy.deepcopy(sphere)
    adaptive3["sizing"]["sources"] = [{"id": "local", "kind": "box", "apply_to": "background",
                                       "target_size": .15, "transition_width": .2,
                                       "region": {"bounds": {"minimum": [.3, .3, .3], "maximum": [1, 1, 1]}}}]
    a3m, a3r, _ = suite.generate("sphere3d-adaptive", adaptive3)
    require(a3r["metadata"]["background_carving"]["reserved_refinement_transition_cells"] > 0,
            "3D refinement transition cavity was not reserved")

    cylinder = copy.deepcopy(sphere)
    cylinder["cad"]["operations"][0].update({"operation": "cylinder", "parameters": {
        "axis_origin": [1.5, 1.5, 1], "direction": [0, 0, 1], "radius": .35, "height": 1}})
    cylinder["sizing"]["surface_size"] = .18
    _, c3r, _ = suite.generate("cylinder3d-sharp-caps", cylinder)
    require(c3r["metadata"]["relative_volume_error"] < .01, "sharp cylinder volume error")

    pierce = copy.deepcopy(cylinder)
    pierce["cad"]["operations"][0]["parameters"].update({"axis_origin": [1.5, 1.5, -.5], "height": 4})
    _, p3r, _ = suite.generate("cylinder3d-pierced-box", pierce)
    require(not p3r["metadata"]["retained_cartesian_outer_boundary"],
            "partial outer coverage must fall back to the full CAD skin")

    box = load(suite.root / "examples/box3d.json")
    _, bxr, _ = suite.generate("internal-box3d", box)
    require(math.isclose(bxr["volume_m3"], 1, rel_tol=1e-10), "internal box volume")
    require(bxr["cell_types"]["hexahedron"] == 216 and bxr["cells"] == 216,
            "aligned internal CAD box should retain the complete Cartesian grid")
    require(bxr["quality"]["face_quality_thresholds_passed"], "aligned box quality failed")

    step = output / "cad" / "sphere.step"
    suite.run("STEP-export", [suite.cli, "export-cad", suite.root / "examples/sphere3d.json",
                               "--entity", "body", "--out", step, "--overwrite"])
    imported = copy.deepcopy(sphere)
    imported["cad"]["operations"] = []
    imported["cad"]["imports"] = [{"id": "body", "file": str(step), "length_unit": "m"}]
    im, ir, _ = suite.generate("step-sphere3d", imported)
    require(abs(ir["metadata"]["cad_measure_volume_m3"] / sr["metadata"]["cad_measure_volume_m3"] - 1) < 1e-10,
            "STEP round-trip changed CAD units or volume")

    suite.exchange("cylinder2d", cm, cr, ["fluent-ascii", "cgns", "cgns-poly", "babelsim-v3"])
    suite.exchange("cylinder2d-adaptive", am, ar, ["fluent-ascii", "cgns"])
    suite.exchange("sphere3d", sm, sr, ["fluent-ascii", "cgns", "cgns-poly"])
    suite.exchange("sphere3d-adaptive", a3m, a3r, ["cgns"])

    scaled = output / "exchange" / "cylinder2d-scaled.json"
    suite.run("explicit-unit-scaling-preserves-bounds", [suite.cli, "convert",
              output / "cylinder2d/mesh.json", "--to", "checkpoint", "--out", scaled,
              "--length-scale", "0.001", "--overwrite"])
    changed = load(scaled)
    locked = changed["metadata"]["domain_bounds_constraint"]
    original_bounds = cm["metadata"]["domain_bounds_constraint"]
    for side in ["minimum_m", "maximum_m"]:
        require(locked[side] == [v * .001 for v in original_bounds[side]],
                "explicit unit scaling failed to rescale locked bounds")
    scaled_report = load(str(scaled) + ".report.json")
    require(math.isclose(scaled_report["volume_m3"], cr["volume_m3"] * 1e-9,
                         rel_tol=1e-11), "unit scaling changed normalized volume")

    # A hand-authored Fluent fixture is independent of this converter's writer.
    fluent = output / "exchange" / "external-tet.msh"
    fluent.write_text('(2 3)\n(10 (0 1 4 0 3))\n(10 (1 1 4 1 3) (0 0 0 1 0 0 0 1 0 0 0 1))\n'
                      '(12 (0 1 1 0))\n(12 (2 1 1 1 2))\n(13 (0 1 4 0))\n'
                      '(13 (3 1 4 3 3) (2 3 1 1 0 4 2 1 1 0 4 3 2 1 0 4 1 3 1 0))\n'
                      '(39 (2 fluid "fluid") ())\n(39 (3 wall "wall") ())\n')
    tet = output / "exchange" / "external-tet.bsm"
    suite.run("Fluent-homogeneous-external-tet", [suite.cli, "convert", fluent, "--to", "babelsim-v3", "--out", tet, "--overwrite"])
    tr = load(str(tet) + ".report.json")
    require(tr["cell_types"]["tetrahedron"] == 1 and math.isclose(tr["volume_m3"], 1/6, rel_tol=1e-12),
            "Fluent homogeneous cell or c0 orientation is wrong")
    suite.probe("external-tet", tet, tr)

    plain = copy.deepcopy(circle)
    plain["optimization"]["enabled"] = False
    pm, pr, _ = suite.generate("cylinder2d-unoptimized", plain)
    suite.run("standalone-optimization", [suite.cli, "optimize", output / "cylinder2d-unoptimized/mesh.json",
                                          "--out", output / "optimized", "--iterations", 8, "--overwrite"])
    om, oreport = load(output / "optimized/mesh.json"), load(output / "optimized/mesh_report.json")
    require(len(om["cells"]) == len(pm["cells"]), "standalone optimizer increased cell count")
    opt = oreport["metadata"]["optimization"]
    require(opt["objective_after"] <= opt["objective_before"] and opt["moved_vertices"] > 0,
            "optimizer did not improve this known improvable case")
    fixed = set(v for c in pm["cells"] if c["region"] != 2 for v in c["vertices"])
    require(all(pm["points_m"][v] == om["points_m"][v] for v in fixed), "optimizer moved background/layer vertices")
    suite.probe("optimized", output / "optimized/mesh.bsm", oreport)

    # Native face topology has no region tags. A reimport must preserve exact
    # Cartesian geometry, including hanging-face boxes, and pair planar nodes.
    for name, original, expected in [("cylinder2d-unoptimized", pm, pr),
                                     ("cylinder2d-adaptive", am, ar),
                                     ("sphere3d", sm, sr)]:
        destination = output / (name + "-native-optimized")
        suite.run(name + "-native-optimization", [suite.cli, "optimize",
                  output / name / "mesh.bsm", "--out", destination,
                  "--iterations", 8, "--overwrite"])
        changed = load(destination / "mesh.json")
        optimized_report = load(destination / "mesh_report.json")
        fixed = {node for c in original["cells"] if c["region"] == 0 for node in c["vertices"]}
        require(all(original["points_m"][v] == changed["points_m"][v] for v in fixed),
                "native reimport optimizer moved Cartesian/hanging-face vertices")
        require(len(original["cells"]) == len(changed["cells"]), "native optimizer added cells")
        require(original["points_m"] == changed["points_m"] and
                optimized_report["metadata"]["optimization"]["moved_vertices"] == 0,
                "native optimizer changed geometry without layer-region protection tags")
        require(changed["generation_dimension"] == original["generation_dimension"],
                "native reimport lost the generation rule")
        require(optimized_report["metadata"]["optimization"]["geometrically_protected_cartesian_cells"]
                >= expected["regions"]["cartesian"], "Cartesian geometric recognition was incomplete")
        if original["generation_dimension"] == 2:
            axis = original["extrusion_axis"]
            pairs = {}
            for i, point in enumerate(original["points_m"]):
                key = tuple(v for a, v in enumerate(point) if a != axis)
                pairs.setdefault(key, []).append(i)
            for pair in pairs.values():
                require(len(pair) == 2 and
                        all(changed["points_m"][pair[0]][a] == changed["points_m"][pair[1]][a]
                            for a in range(3) if a != axis), "native planar optimizer unpaired slab nodes")
        suite.probe(name + "-native-optimized", destination / "mesh.bsm", optimized_report)

    for axis in ["x", "y"]:
        checkpoint = output / "exchange" / ("channel2d-" + axis + "-native.json")
        suite.run("native-planar-inference-" + axis, [suite.cli, "convert",
                  output / ("channel2d-" + axis) / "mesh.bsm", "--to", "checkpoint",
                  "--out", checkpoint, "--overwrite"])
        recovered = load(checkpoint)
        require(recovered["generation_dimension"] == 2 and
                recovered["extrusion_axis"] == "xyz".index(axis), "native axis inference failed")

    unknown = output / "exchange" / "native-tags-unknown.msh"
    suite.run("native-tags-unknown-Fluent-write", [suite.cli, "convert",
              output / "cylinder2d/mesh.bsm", "--to", "fluent-ascii", "--out", unknown])
    unknown_checkpoint = output / "exchange" / "native-tags-unknown.json"
    suite.run("native-tags-unknown-Fluent-read", [suite.cli, "convert", unknown,
              "--to", "checkpoint", "--out", unknown_checkpoint])
    require(load(unknown_checkpoint)["metadata"]["generation_region_tags_available"] is False,
            "exchange relabelled unknown generation regions as available")

    bad = copy.deepcopy(circle); bad["background"]["cell_counts"][2] = 2
    suite.negative("reject-two-extrusion-layers", bad, "exactly one")
    bad = copy.deepcopy(circle); bad["meshing"]["thickness"] = 1.1
    suite.negative("reject-extrusion-outside-background", bad, "fixed background bounds")
    bad = copy.deepcopy(sphere); bad["background"]["bounds"]["maximum"][0] = 2.9
    suite.negative("reject-CAD-domain-outside-background", bad, "fixed background/domain bounds")
    bad = copy.deepcopy(cm); bad["points_m"][0][0] = 2.5
    save(output / "outside-bounds.json", bad)
    suite.run("reject-checkpoint-outside-fixed-bounds", [suite.cli, "check", output / "outside-bounds.json"], expected=False)
    bad_output = output / "outside-bounds-optimized"
    suite.run("reject-optimization-outside-fixed-bounds", [suite.cli, "optimize", output / "outside-bounds.json", "--out", bad_output], expected=False)
    require(not (bad_output / "mesh.json").exists(), "rejected out-of-bounds optimizer published a mesh")
    bad = copy.deepcopy(circle); bad["resources"]["max_cells"] = 10
    suite.negative("reject-cell-budget", bad, "max_cells")
    bad = copy.deepcopy(circle); bad["optimization"]["topology_swaps"] = "yes"
    suite.negative("reject-topology-option-type", bad, "topology_swaps must be boolean")
    bad = copy.deepcopy(circle); bad["resources"]["max_memory_mb"] = 1
    suite.negative("reject-memory-budget", bad, "memory")
    bad = copy.deepcopy(circle); bad["boundary_layers"]["targets"][0]["first_layer"] = 3
    suite.negative("reject-layer-clearance", bad, "clearance")
    bad = copy.deepcopy(circle); bad["meshing"]["thickness"] = .005
    suite.negative("reject-strict-face-quality", bad, "face quality thresholds failed")
    bad = copy.deepcopy(cm); bad["cells"][0]["type"] = 10
    save(output / "invalid-planar.json", bad)
    suite.run("reject-2D-tetrahedron-checkpoint", [suite.cli, "check", output / "invalid-planar.json"], expected=False)
    binary = output / "exchange" / "binary.msh"; binary.write_bytes(b'\x89HDF\r\n\x1a\n')
    suite.run("reject-Fluent-HDF5", [suite.cli, "convert", binary, "--to", "babelsim-v3",
                                     "--out", output / "bad-binary.bsm"], expected=False)
    suite.run("reject-output-checkpoint-alias", [suite.cli, "convert", output / "cylinder2d/mesh.json",
              "--to", "fluent-ascii", "--out", output / "alias.msh",
              "--checkpoint", output / "alias.msh", "--overwrite"], expected=False)
    suite.run("reject-checkpoint-input-alias", [suite.cli, "convert", output / "cylinder2d/mesh.json",
              "--to", "fluent-ascii", "--out", output / "alias2.msh",
              "--checkpoint", output / "cylinder2d/mesh.json", "--overwrite"], expected=False)
    suite.complete()


if __name__ == "__main__":
    main()
