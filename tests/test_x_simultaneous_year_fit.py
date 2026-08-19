import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
ACTUAL_MANIFEST = (
    REPO.parent / "XGBoost/output/selected"
    / "X_pb23_v3_fid3_6v5_rwr6range5v1_xgb_v1"
    / "fit_scan_manifest.pb23_pb24_simultaneous_v1.json"
)
sys.path.insert(0, str(REPO / "fitER"))
WORKFLOW_SPEC = importlib.util.spec_from_file_location(
    "x_simultaneous_year_fit_workflow", REPO / "fitER/x_simultaneous_year_fit_workflow.py"
)
WORKFLOW = importlib.util.module_from_spec(WORKFLOW_SPEC)
WORKFLOW_SPEC.loader.exec_module(WORKFLOW)
SUBMIT_SPEC = importlib.util.spec_from_file_location(
    "submit_x_simultaneous_year_fit_manifest", REPO / "scripts/submit_x_simultaneous_year_fit_manifest.py"
)
SUBMIT = importlib.util.module_from_spec(SUBMIT_SPEC)
SUBMIT_SPEC.loader.exec_module(SUBMIT)


@unittest.skipUnless(ACTUAL_MANIFEST.is_file(), "ML simultaneous-year manifest unavailable")
class SimultaneousYearContractTest(unittest.TestCase):
    def test_actual_manifest_loads_without_root_io(self):
        manifest = WORKFLOW.load_manifest(ACTUAL_MANIFEST)
        self.assertEqual(tuple(p["key"] for p in manifest["working_points"]), WORKFLOW.POINTS)
        self.assertEqual(manifest["nominal_fit_contract"]["signal"]["sigma_gev"]["range"], [0.002, 0.008])

    def test_rejects_changed_fit_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            manifest = json.loads(ACTUAL_MANIFEST.read_text())
            manifest["nominal_fit_contract"]["fit_type"] = "merged_unbinned"
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "fit_type"):
                WORKFLOW.load_manifest(path)

    def test_dag_has_seven_fits_and_no_toys(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with mock.patch.object(SUBMIT, "AFS_ROOT", root / "afs"), mock.patch.object(SUBMIT, "EOS_RESULTS_ROOT", root / "eos"):
                task = SUBMIT.create_submission(ACTUAL_MANIFEST, "phase1_test")
            dag = Path(task["submission_dir"], "fit_scan.dag").read_text()
            self.assertEqual(dag.count("JOB FIT_"), 7)
            self.assertIn("PARENT PREPARE CHILD", dag)
            self.assertIn("FINAL AGGREGATE", dag)
            self.assertNotIn("TOY", dag.upper())
            self.assertEqual(task["toy_count"], 0)
            self.assertFalse(Path(task["output_dir"]).exists())

    def test_phase1_semantics_are_explicit(self):
        source = (REPO / "fitER/x_simultaneous_year_fit_workflow.py").read_text()
        self.assertIn('CALIBRATION = "none_sqrt_q0_heuristic"', source)
        self.assertIn('"p0": None', source)
        self.assertIn('"merged_mass_distribution"', source)
        self.assertIn('"display_only"', source)


if __name__ == "__main__":
    unittest.main()
