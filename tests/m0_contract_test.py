#!/usr/bin/env python3
"""CLI/schema contract checks plus smoke coverage for the accepted M1/M2 gates."""

import json
import pathlib
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import jsonschema

BINARY = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
del sys.argv[1:3]


class M0ContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = BINARY
        cls.root = ROOT
        cls.schema = json.loads((ROOT / "schema/project.schema.json").read_text())
        cls.example = json.loads(
            (ROOT / "acceptance/projects/cartesian-box.json").read_text()
        )
        jsonschema.Draft7Validator.check_schema(cls.schema)

    def invoke(self, *args):
        return subprocess.run(
            [str(self.binary), *args], text=True, capture_output=True, check=False
        )

    def test_version_and_help(self):
        version = self.invoke("--version")
        self.assertEqual(version.returncode, 0, version.stderr)
        self.assertEqual(version.stdout.strip(), "0.1.0")
        help_result = self.invoke("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("generate", help_result.stdout)

    def test_later_commands_fail_explicitly(self):
        result = self.invoke("verify-suite")
        self.assertEqual(result.returncode, 3)
        self.assertIn("not implemented", result.stderr)

    def test_unknown_command_and_option_fail(self):
        self.assertEqual(self.invoke("not-a-command").returncode, 2)
        self.assertEqual(self.invoke("--not-a-real-option").returncode, 2)

    def test_example_matches_strict_schema(self):
        jsonschema.validate(self.example, self.schema)

    def test_surface_size_is_a_scalar_boundary_cell_cap(self):
        project = json.loads(json.dumps(self.example))
        project["sizing"]["surface_size"] = 0.02
        jsonschema.validate(project, self.schema)
        for invalid in (0, -0.01, "0.02", [0.02, 0.02, 0.02]):
            project["sizing"]["surface_size"] = invalid
            self.assertTrue(list(jsonschema.Draft7Validator(self.schema).iter_errors(project)))
        project["sizing"].pop("surface_size")
        project["sizing"]["surface_size_xyz"] = [0.02, 0.02, 0.02]
        self.assertTrue(list(jsonschema.Draft7Validator(self.schema).iter_errors(project)))

    def test_curvature_radius_is_separate_from_cell_size(self):
        project = json.loads(json.dumps(self.example))
        project["sizing"]["sources"] = [{
            "id": "curved_wall", "kind": "curvature", "target_size": 0.02,
            "region": {"cad_group": "walls", "max_radius": 0.1},
        }]
        jsonschema.validate(project, self.schema)
        project["sizing"]["sources"][0]["region"].pop("max_radius")
        self.assertTrue(list(jsonschema.Draft7Validator(self.schema).iter_errors(project)))

    def test_explicit_surface_tolerance_contract(self):
        project = json.loads(json.dumps(self.example))
        project["quality"]["max_surface_deviation_m"] = 0.001
        jsonschema.validate(project, self.schema)
        with tempfile.TemporaryDirectory(prefix="babelmesh-tolerance-") as tmp:
            path = pathlib.Path(tmp) / "project.json"
            for invalid in [0, -0.001, "0.001", None]:
                project["quality"]["max_surface_deviation_m"] = invalid
                self.assertTrue(list(jsonschema.Draft7Validator(self.schema).iter_errors(project)))
                path.write_text(json.dumps(project))
                result = self.invoke("generate", str(path), "--out", str(pathlib.Path(tmp) / "mesh"))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("max_surface_deviation_m", result.stderr)

    def test_unknown_keys_and_invalid_units_are_rejected(self):
        validator = jsonschema.Draft7Validator(self.schema)
        extra_key = {**self.example, "silently_ignored": True}
        self.assertGreater(len(list(validator.iter_errors(extra_key))), 0)
        bad_operation = json.loads(json.dumps(self.example))
        bad_operation["cad"]["operations"][0]["parameters"] = {"fudge_factor": 3}
        self.assertGreater(len(list(validator.iter_errors(bad_operation))), 0)
        bad_selector = json.loads(json.dumps(self.example))
        bad_selector["cad"]["groups"] = [
            {"name": "walls", "kind": "wall", "selector": {"future_magic": True}}
        ]
        self.assertGreater(len(list(validator.iter_errors(bad_selector))), 0)
        bad_units = json.loads(json.dumps(self.example))
        bad_units["units"]["length"] = "ft"
        self.assertGreater(len(list(validator.iter_errors(bad_units))), 0)

    def test_acceptance_gate_statuses_have_evidence(self):
        suite = json.loads((ROOT / "acceptance/suite.json").read_text())
        self.assertEqual(suite["status"], "in_progress")
        self.assertGreaterEqual(len(suite["gates"]), 10)
        self.assertTrue(all(gate["required"] for gate in suite["gates"]))
        statuses = {gate["id"]: gate["status"] for gate in suite["gates"]}
        self.assertEqual(statuses["M0-contract"], "passed")
        self.assertEqual(statuses["M1-cad-primitives-booleans-step"], "passed")
        self.assertEqual(statuses["M2-uniform-cartesian-and-v3-reader"], "passed")
        self.assertEqual(statuses["M3-adaptive-conforming-topology"], "passed")
        self.assertEqual(statuses["M4-surface-snapping-and-fluid-extraction"], "in_progress")
        self.assertTrue(all(statuses[gate_id] == "planned" for gate_id in statuses if gate_id.startswith(("M5", "M6", "M7", "M8"))))

    def test_inspect_cad_writes_geometry_report(self):
        with tempfile.TemporaryDirectory(prefix="babelmesh-cli-test-") as tmp:
            output_dir = pathlib.Path(tmp) / "cad-report"
            project = self.root / "acceptance/projects/cartesian-box.json"
            result = self.invoke("inspect-cad", str(project), "--out", str(output_dir))
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output_dir / "geometry_report.json").read_text())
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["entities"][0]["id"], "body")
            self.assertAlmostEqual(report["entities"][0]["metrics"]["volume_m3"], 1.0)
            surface = (output_dir / "cad_surface.vtu").read_text()
            self.assertIn('Name="cad_entity_id"', surface)
            self.assertIn('NumberOfCells="12"', surface)
            second = self.invoke("inspect-cad", str(project), "--out", str(output_dir))
            self.assertEqual(second.returncode, 1)
            self.assertIn("overwrite is false", second.stderr)

    def test_generate_writes_uniform_hex_outputs(self):
        with tempfile.TemporaryDirectory(prefix="babelmesh-generate-test-") as tmp:
            output_dir = pathlib.Path(tmp) / "mesh-output"
            project_file = self.root / "acceptance/projects/cartesian-box.json"
            result = self.invoke("generate", str(project_file), "--out", str(output_dir))
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output_dir / "mesh_report.json").read_text())
            self.assertEqual(report["generator"], "uniform-cartesian-hex")
            self.assertFalse(report["cad_clipping_applied"])
            self.assertIn("M4", report["limitations"][0])
            self.assertEqual(report["cells"], 12 ** 3)
            self.assertEqual(report["boundary_faces"], 6 * (12 ** 2))
            self.assertAlmostEqual(report["quality"]["min_scaled_jacobian"], 1.0)
            self.assertAlmostEqual(report["quality"]["max_cell_aspect_ratio"], 1.0)
            self.assertAlmostEqual(report["cell_volume_m3"], (0.25 ** 3))
            self.assertTrue((output_dir / "mesh.bsm").read_text().startswith("BABELSIM_MESH 3\n"))
            vtu = (output_dir / "mesh.vtu").read_text()
            self.assertIn('NumberOfCells="1728"', vtu)
            second = self.invoke("generate", str(project_file), "--out", str(output_dir))
            self.assertEqual(second.returncode, 1)
            self.assertIn("output already exists", second.stderr)

    def test_generate_applies_local_region_refinement(self):
        with tempfile.TemporaryDirectory(prefix="babelmesh-adaptive-cli-test-") as tmp:
            output_dir = pathlib.Path(tmp) / "adaptive"
            project = json.loads(json.dumps(self.example))
            project["sizing"]["sources"] = [
                {
                    "id": "small_box_region",
                    "kind": "box",
                    "target_size": 0.125,
                    "region": {
                        "bounds": {
                            "minimum": [0.01, 0.01, 0.01],
                            "maximum": [0.12, 0.12, 0.12],
                        }
                    },
                }
            ]
            project_file = pathlib.Path(tmp) / "adaptive-project.json"
            project_file.write_text(json.dumps(project))
            result = self.invoke("generate", str(project_file), "--out", str(output_dir))
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output_dir / "mesh_report.json").read_text())
            self.assertEqual(report["generator"], "adaptive-cartesian-octree")
            self.assertGreater(report["cells"], 12 ** 3)
            self.assertGreater(report["maximum_level"], 0)
            self.assertGreater(report["sizing_source_intersection_counts"]["small_box_region"], 0)
            self.assertEqual(
                report["sizing_sources"]["small_box_region"]["marked_cell_count"],
                report["sizing_source_intersection_counts"]["small_box_region"],
            )
            self.assertLessEqual(report["total_volume_m3"], 27.0 + 1.0e-12)
            vtk = ET.parse(output_dir / "mesh.vtu").getroot()
            arrays = {item.attrib.get("Name"): item.text or "" for item in vtk.findall(".//Cells/DataArray")}
            types = [int(value) for value in arrays["types"].split()]
            self.assertIn(42, types)
            self.assertEqual(len(arrays["faceoffsets"].split()), report["cells"])
            cell_data = {
                item.attrib.get("Name"): item.text or ""
                for item in vtk.findall(".//CellData/DataArray")
            }
            self.assertEqual(len(cell_data["level"].split()), report["cells"])
            self.assertEqual(len(cell_data["source_small_box_region"].split()), report["cells"])

    def test_generate_cad_staircase_extraction_is_reported(self):
        with tempfile.TemporaryDirectory(prefix="babelmesh-cad-staircase-cli-test-") as tmp:
            output_dir = pathlib.Path(tmp) / "fluid-mesh"
            project = json.loads(json.dumps(self.example))
            project["sizing"]["surface_size"] = 0.25
            project_file = pathlib.Path(tmp) / "cad-staircase-project.json"
            project_file.write_text(json.dumps(project))
            result = self.invoke("generate", str(project_file), "--cad-staircase", "--out", str(output_dir))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("sampled bidirectional surface-distance", result.stdout)
            report = json.loads((output_dir / "mesh_report.json").read_text())
            self.assertTrue(report["cad_clipping_applied"])
            self.assertTrue(report["cad_surface_refinement_applied"])
            self.assertTrue(report["cad_snapping_applied"])
            self.assertEqual(report["cad_snap_fraction"], 1.0)
            self.assertAlmostEqual(report["cad_surface_deviation_m"], 0.0)
            self.assertEqual(report["cad_removed_cells"], 64)
            self.assertEqual(report["cells"], 12 ** 3 - 64)
            self.assertAlmostEqual(report["cad_fluid_domain_volume_m3"], 26.0)
            self.assertAlmostEqual(report["total_volume_m3"], 26.0)
            self.assertAlmostEqual(report["cad_mesh_volume_relative_error"], 0.0)
            self.assertIn(["cad_boundary", "wall"], report["patches"])
            self.assertFalse(report["cad_surface_deviation_overridden"])
            project["quality"]["max_surface_deviation_m"] = 0.001
            project_file.write_text(json.dumps(project))
            explicit_dir = pathlib.Path(tmp) / "explicit-tolerance"
            result = self.invoke("generate", str(project_file), "--cad-staircase", "--out", str(explicit_dir))
            self.assertEqual(result.returncode, 0, result.stderr)
            explicit = json.loads((explicit_dir / "mesh_report.json").read_text())
            self.assertTrue(explicit["cad_surface_deviation_overridden"])
            self.assertEqual(explicit["cad_surface_deviation_tolerance_m"], 0.001)
            self.assertEqual(explicit["quality"], report["quality"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
