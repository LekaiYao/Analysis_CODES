import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
MANIFEST = Path(
    "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/"
    "Psi2S_pb24_v1_fid1_6v1_rwr6range4v1_xgb_v1/"
    "fit_scan_manifest.psi2s_nominal_v1.json"
)
CONFIG = REPO / "fitER/configs/pbpb24_psi2s_xeff45_70_exploratory_v1.json"


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


workflow = load(
    "psi2s_high_eff_workflow_test",
    "fitER/psi2s_high_efficiency_exploratory_workflow.py",
)
submitter = load(
    "psi2s_high_eff_submitter_test",
    "scripts/submit_psi2s_high_efficiency_exploratory.py",
)


class Psi2SHighEfficiencyExploratoryTest(unittest.TestCase):
    def test_real_config_is_one_off_and_exact(self):
        manifest, config = workflow.load_inputs(REPO, MANIFEST, CONFIG)
        self.assertEqual(manifest["train_tag"], config["train_tag"])
        self.assertTrue(config["excluded_from_regular_workflow"])
        self.assertEqual(tuple(p["key"] for p in config["working_points"]), workflow.POINTS)
        self.assertGreater(config["working_points"][0]["threshold"],
                           config["working_points"][-1]["threshold"])

    def test_regular_seven_point_workflow_is_unchanged(self):
        source = (REPO / "scripts/submit_psi2s_data_gaussian_manifest.py").read_text()
        self.assertIn("(10, 15, 20, 25, 30, 35, 40)", source)
        self.assertNotIn("45, 50, 55, 60, 65, 70", source)

    def test_dag_has_six_exploratory_fit_nodes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with mock.patch.object(submitter, "AFS_ROOT", root / "afs"), \
                    mock.patch.object(submitter, "EOS_ROOT", root / "eos"):
                task = submitter.create_submission(MANIFEST, CONFIG, "test")
            dag = Path(task["submission_dir"], "fit_scan.dag").read_text()
            self.assertEqual(dag.count("CATEGORY FIT"), 6)
            self.assertIn("PARENT PREPARE CHILD", dag)
            self.assertTrue(task["excluded_from_regular_workflow"])


if __name__ == "__main__":
    unittest.main()
