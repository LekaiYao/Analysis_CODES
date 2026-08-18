import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "fitER/psi2s_fit_scan_workflow.py"
SPEC = importlib.util.spec_from_file_location("psi2s_workflow", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Psi2SWorkflowTest(unittest.TestCase):
    def test_aggregate_marks_boundaries_as_warnings(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            points = []
            for index, value in enumerate((10, 15, 20, 25, 30, 35, 40)):
                key = f"psi2seff{value}"
                point = {
                    "key": key, "target_weighted_efficiency": value / 100,
                    "achieved_weighted_efficiency": value / 100,
                    "threshold": 1.0 - value / 100,
                }
                points.append(point)
                point_dir = output / key
                point_dir.mkdir()
                result = {
                    "fit_status": 0, "cov_qual": 3, "edm": 1e-5,
                    "mc_fit_status": 0, "mc_cov_qual": 3, "mc_edm": 1e-5,
                    "mc_parameter_boundary": False,
                    "parameter_boundary_flags": ["width_scale"] if index == 0 else [],
                    "chebyshev_a0": 0.1, "chebyshev_a1": -0.1,
                }
                payload = {
                    "point": key, "target_weighted_efficiency": value / 100,
                    "achieved_weighted_efficiency": value / 100,
                    "threshold": 1.0 - value / 100, "data_entries": 100 + index,
                    "mc_entries": 50, "mc_sumw": 40, "mc_sumw2": 35,
                    "mc_effective_entries": 45, "yield": 10, "yield_error": 3,
                    "fit_status": 0, "covQual": 3, "EDM": 1e-5, "mean": 3.686,
                    "width_scale": 1.0, "signal_over_background": 0.1,
                    "signal_over_sqrt_signal_plus_background": 1.0,
                    "parameter_boundary_flags": result["parameter_boundary_flags"],
                    "artifact_paths": {}, "fit_result": result,
                }
                (point_dir / "point_manifest.json").write_text(json.dumps(payload))
            manifest = {
                "train_tag": "Psi2S_test", "working_points": points,
                "nominal_fit_contract": {},
            }
            manifest_path = output / "input.json"
            manifest_path.write_text(json.dumps(manifest))
            with mock.patch.object(MODULE.subprocess, "check_output", return_value="test\n"):
                MODULE.aggregate(REPO, manifest_path, manifest, output)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "passed_with_warnings")
            self.assertIn("width_scale", validation["warnings"][0])

    def test_resolve_relative_input_against_manifest(self):
        manifest = Path("/tmp/example/manifest.json")
        self.assertEqual(
            MODULE.resolve(manifest, "DATA.root"), Path("/tmp/example/DATA.root")
        )


if __name__ == "__main__":
    unittest.main()
