import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SOURCE = Path(__file__).resolve().parents[1] / "scripts/submit_x_nominal_fit_manifest.py"
SPEC = importlib.util.spec_from_file_location("submit_x_nominal_fit_manifest", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ManifestTaskTest(unittest.TestCase):
    def manifest(self):
        return {
            "contract": MODULE.SUPPORTED_CONTRACT,
            "schema_version": 2,
            "train_tag": "X_pb24_test_xgb_v1",
            "inputs": {"signal_mc": {"path": "/tmp/mc.root", "tree": "ntmix_X3872"}},
            "fit_contract": {
                "data_mc_width_scale_range": [0.9, 1.5],
                "signal_model": "double Gaussian fitted to the scored signal MC",
            },
            "working_points": [{"key": key} for key in MODULE.POINTS],
        }

    def test_loads_supported_task_without_root_io(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(self.manifest()))
            _, tag = MODULE.load_task(path)
            self.assertEqual(tag, "X_pb24_test_xgb_v1")

    def test_rejects_changed_schema_or_width_scale(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            manifest = self.manifest()
            manifest["schema_version"] = 3
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "schema_version"):
                MODULE.load_task(path)
            manifest["schema_version"] = 2
            manifest["fit_contract"]["data_mc_width_scale_range"] = [0.9, 1.25]
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "width-scale range"):
                MODULE.load_task(path)

    def test_data_only_contract_requires_explicit_compatibility_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            manifest = {
                "contract": MODULE.DATA_ONLY_COMPAT_CONTRACT,
                "schema_version": 2,
                "train_tag": "X_pb24_test_xgb_v1",
                "nominal_fit_contract": {
                    "signal": {"sigma_gev": {"range": [0.002, 0.008]}}
                },
                "working_points": [{"key": key} for key in MODULE.POINTS],
            }
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "compatibility-only"):
                MODULE.load_task(path)
            _, tag = MODULE.load_task(path, allow_data_only_compat=True)
            self.assertEqual(tag, "X_pb24_test_xgb_v1")

    def test_layout_keeps_afs_submission_and_eos_results_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "manifest.json"
            path.write_text(json.dumps(self.manifest()))
            with mock.patch.object(MODULE, "AFS_ROOT", root / "afs"), \
                 mock.patch.object(MODULE, "EOS_RESULTS_ROOT", root / "eos"):
                task = MODULE.create_submission(path)
            self.assertTrue(Path(task["submission_dir"], "fit_scan.dag").is_file())
            self.assertFalse(Path(task["output_dir"]).exists())
            self.assertIn("/afs/", task["submission_dir"])
            self.assertIn("/eos/", task["output_dir"])
            self.assertTrue(task["output_dir"].endswith(
                "mc_shape_nominal_widthscale0p9_1p5_v1"
            ))


if __name__ == "__main__":
    unittest.main()
