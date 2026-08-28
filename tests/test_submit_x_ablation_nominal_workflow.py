import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
SOURCE = REPO / "scripts/submit_x_ablation_nominal_workflow.py"
SPEC = importlib.util.spec_from_file_location("submit_x_ablation_nominal_workflow", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FullNominalDagTest(unittest.TestCase):
    def test_comparison_runs_after_seven_point_aggregate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text("{}")
            simultaneous_output = root / "eos-sim"
            submission_dir = root / "afs-submit"
            fake_task = {
                "contract": "old", "schema_version": 2,
                "anchor_train_tag": "X_pb23_test",
                "input_manifest": str(manifest),
                "output_dir": str(simultaneous_output),
                "submission_dir": str(submission_dir),
            }

            def fake_create(*_args, **_kwargs):
                (submission_dir / "logs").mkdir(parents=True)
                for name in ("prepare.sub", "fit.sub", "aggregate.sub"):
                    (submission_dir / name).write_text(
                        'should_transfer_files = YES\ntransfer_output_files = ""\n'
                    )
                return dict(fake_task)

            with mock.patch.object(
                MODULE.simultaneous_submit, "create_submission", side_effect=fake_create
            ), mock.patch.object(MODULE, "COMPARISON_RESULTS_ROOT", root / "eos-compare"):
                task = MODULE.create_submission(manifest)
            dag = (submission_dir / "fit_scan.dag").read_text()
            self.assertEqual(dag.count("JOB FIT_"), 7)
            self.assertIn("JOB AGGREGATE aggregate.sub", dag)
            self.assertIn("PARENT AGGREGATE CHILD COMPARISON", dag)
            self.assertNotIn("FINAL AGGREGATE", dag)
            comparison = (submission_dir / "comparison.sub").read_text()
            self.assertIn("--point best", comparison)
            self.assertIn("should_transfer_files = YES", comparison)
            self.assertIn('transfer_output_files = ""', comparison)
            persisted = json.loads((submission_dir / "task.json").read_text())
            self.assertEqual(persisted["contract"], task["contract"])


if __name__ == "__main__":
    unittest.main()
