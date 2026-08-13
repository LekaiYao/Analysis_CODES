import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "fitER" / "x_fit_scan_workflow.py"
ACTUAL_MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "XGBoost/output/selected"
    / "X_pb24_v18_fid19_6v5_rwr6range7dr025v1_xgb_v1"
    / "fit_scan_manifest.json"
)
NOMINAL_MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "XGBoost/output/selected"
    / "X_pb24_v13_fid14_6v5_rwr6range4v1_xgb_v1"
    / "fit_scan_manifest.data_only_nominal_v2.json"
)
SPEC = importlib.util.spec_from_file_location("x_fit_scan_workflow", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ManifestLoadTest(unittest.TestCase):
    def test_minimal_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps({
                "contract": "pbpb24_x_weighted_efficiency_fit_scan",
                "train_tag": "X_test", "inputs": {}, "fit_contract": {},
                "working_points": [{"key": "xeff10", "selection": "Prediction > 0.9"}],
            }))
            self.assertEqual(MODULE.load_manifest(path)["train_tag"], "X_test")

    def test_rejects_other_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps({
                "contract": "other", "train_tag": "X_test", "inputs": {},
                "fit_contract": {}, "working_points": [],
            }))
            with self.assertRaisesRegex(RuntimeError, "unsupported"):
                MODULE.load_manifest(path)

    def test_data_only_schema_v2_requires_narrow_sigma_contract(self):
        base = {
            "contract": MODULE.DATA_ONLY_CONTRACT,
            "schema_version": 2,
            "train_tag": "X_test", "inputs": {}, "working_points": [],
            "nominal_fit_contract": {
                "signal": {"sigma_gev": {"range": [0.002, 0.008]}}
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(base))
            self.assertEqual(MODULE.load_manifest(path)["train_tag"], "X_test")
            base["nominal_fit_contract"]["signal"]["sigma_gev"]["range"] = [0.0001, 0.01]
            path.write_text(json.dumps(base))
            with self.assertRaisesRegex(RuntimeError, "schema v2 requires sigma range"):
                MODULE.load_manifest(path)

    @unittest.skipUnless(NOMINAL_MANIFEST.is_file(), "nominal ML manifest is unavailable")
    def test_actual_data_only_nominal_manifest_loads_without_opening_root_inputs(self):
        manifest = MODULE.load_manifest(NOMINAL_MANIFEST)
        self.assertEqual(manifest["contract"], MODULE.DATA_ONLY_CONTRACT)
        self.assertEqual(manifest["schema_version"], 2)
        self.assertEqual(manifest["fit_contract"]["signal"]["model"], "single_gaussian")
        self.assertEqual(manifest["fit_contract"]["signal"]["sigma_gev"]["range"],
                         MODULE.DATA_ONLY_SIGMA_RANGE_V2)
        self.assertEqual([MODULE.point_threshold(point) for point in manifest["working_points"]],
                         [point["threshold"] for point in manifest["working_points"]])

    @unittest.skipUnless(ACTUAL_MANIFEST.is_file(), "actual ML manifest is unavailable")
    def test_actual_pbpb24_manifest_loads_without_opening_root_inputs(self):
        manifest = MODULE.load_manifest(ACTUAL_MANIFEST)
        self.assertEqual(manifest["schema_version"], 2)
        self.assertEqual(sorted(manifest["inputs"]), ["data", "signal_mc"])
        self.assertEqual([point["key"] for point in manifest["working_points"]], [
            "xeff10", "xeff15", "xeff20", "xeff25", "xeff30", "xeff35", "xeff40",
        ])


if __name__ == "__main__":
    unittest.main()
