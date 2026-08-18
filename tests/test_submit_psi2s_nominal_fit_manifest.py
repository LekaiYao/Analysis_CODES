import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "scripts/submit_psi2s_nominal_fit_manifest.py"
SPEC = importlib.util.spec_from_file_location("psi2s_submit", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def valid_manifest():
    fiducial = "Bpt > 10 && Btrk2dR <= 0.25"
    points = []
    for value, threshold in zip((10, 15, 20, 25, 30, 35, 40),
                                (0.97, 0.95, 0.93, 0.91, 0.89, 0.87, 0.85)):
        points.append({
            "key": f"psi2seff{value}", "target_weighted_efficiency": value / 100,
            "achieved_weighted_efficiency": value / 100 - 0.0001,
            "threshold": threshold,
            "selection": f"({fiducial}) && (Prediction > {threshold})",
            "fiducial_score_selected_data_entries": value * 10,
        })
    return {
        "contract": MODULE.SUPPORTED_CONTRACT, "schema_version": 1,
        "channel": "Psi2S", "system": "PbPb", "dataset": "pb24",
        "train_tag": "Psi2S_pb24_v1_fid1_6v1_rwr6range4v1_xgb_v1",
        "fiducial_selection": {"expression": fiducial, "profile": "pb24_fid1"},
        "inputs": {
            "data": {"path": "DATA.root", "tree": "ntmix", "mass_branch": "Bmass",
                     "score_branch": "Prediction", "event_weight": "unit"},
            "signal_mc": {"path": "MC.root", "tree": "ntmix_PSI2S",
                          "mass_branch": "Bmass", "score_branch": "Prediction",
                          "event_weight_branch": "Reweight",
                          "weight_usage": "signal_shape_and_efficiency"},
        },
        "score": {"branch": "Prediction", "comparison_operator": ">",
                  "equality_passes": False, "threshold_boundary": "exclusive"},
        "threshold_provenance": {"definition": "weighted signal efficiency",
                                 "event_weight_branch": "Reweight"},
        "nominal_fit_contract": {
            "version": 1, "fit_type": "extended_unbinned", "mass_range_gev": [3.6, 3.8],
            "data_event_weight": "unit",
            "signal": {"model": "double_gaussian_mc_shape", "shared_mean": True,
                       "shape_source": "weighted_signal_mc",
                       "event_weight_branch": "Reweight",
                       "fixed_from_mc": ["sigma1", "sigma2", "fraction"],
                       "data_mean_gev": {"range": [3.6811, 3.6911]},
                       "data_mc_width_scale": {"range": [0.9, 1.15]}},
            "background": {"model": "chebyshev", "order": 2,
                           "coefficient_ranges": {"a0": [-0.8, 0.8], "a1": [-0.8, 0.8]},
                           "additional_stability_models_required": False},
        },
        "working_points": points,
    }


class SubmitPsi2STest(unittest.TestCase):
    def write_manifest(self, directory, payload):
        path = Path(directory) / "manifest.json"
        path.write_text(json.dumps(payload))
        return path

    def test_load_accepts_exact_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, tag = MODULE.load_task(self.write_manifest(tmp, valid_manifest()))
            self.assertTrue(tag.startswith("Psi2S_pb24"))

    def test_rejects_wrong_weight_and_point_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = valid_manifest()
            payload["inputs"]["signal_mc"]["event_weight_branch"] = "unit"
            with self.assertRaisesRegex(RuntimeError, "Reweight"):
                MODULE.load_task(self.write_manifest(tmp, payload))
            payload = valid_manifest()
            payload["working_points"].reverse()
            with self.assertRaisesRegex(RuntimeError, "working points"):
                MODULE.load_task(self.write_manifest(tmp, payload))

    def test_creates_independent_dag(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = self.write_manifest(root, valid_manifest())
            with mock.patch.object(MODULE, "AFS_ROOT", root / "afs"), \
                 mock.patch.object(MODULE, "EOS_RESULTS_ROOT", root / "eos"):
                task = MODULE.create_submission(manifest, "test")
            dag = Path(task["submission_dir"], "fit_scan.dag").read_text()
            self.assertIn("JOB PREPARE", dag)
            self.assertEqual(dag.count("CATEGORY FIT"), 7)
            self.assertIn("FINAL AGGREGATE", dag)
            self.assertIn("source ROOT once", task["io_plan"])


if __name__ == "__main__":
    unittest.main()
