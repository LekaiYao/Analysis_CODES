import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/submit_psi2s_simultaneous_year_fit_manifest.py"
REAL_MANIFEST = Path(
    "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/"
    "Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/"
    "fit_scan_manifest.pb23_pb24_psi2s_simultaneous_v1.json"
)


def load_module():
    spec = importlib.util.spec_from_file_location("submit_psi2s_simultaneous", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@unittest.skipUnless(REAL_MANIFEST.is_file(), "official ML manifest unavailable")
class Psi2SSubmissionTest(unittest.TestCase):
    def test_official_manifest_validates(self):
        module = load_module()
        manifest, tag = module.load_task(REAL_MANIFEST)
        self.assertEqual(manifest["contract"], module.CONTRACT)
        self.assertTrue(tag.startswith("Psi2S_pb23_"))
        self.assertEqual(len(manifest["working_points"]), 7)

    def test_width_scale_contract_is_strict(self):
        module = load_module()
        manifest = copy.deepcopy(json.loads(REAL_MANIFEST.read_text()))
        manifest["nominal_fit_contract"]["signal"]["data_fit"][
            "category_width_scale"
        ]["range"] = [0.9, 1.5]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "width-scale range"):
                module.load_task(path)

    def test_submission_dag_has_prepare_seven_fits_and_final(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with mock.patch.object(module, "AFS_ROOT", root / "afs"), \
                    mock.patch.object(module, "EOS_RESULTS_ROOT", root / "results"):
                task = module.create_submission(REAL_MANIFEST, "test_v1")
            dag = (Path(task["submission_dir"]) / "fit_scan.dag").read_text()
            self.assertEqual(dag.count("JOB FIT_"), 7)
            self.assertIn("JOB PREPARE prepare.sub", dag)
            self.assertIn("FINAL AGGREGATE aggregate.sub", dag)
            self.assertNotIn("TOY", dag)
            self.assertTrue(task["significance_scope"].startswith("fit-only"))


if __name__ == "__main__":
    unittest.main()
