import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


workflow = load_module(
    "psi2s_splot_workflow_test",
    "plotER/Validation/psi2s_pbpb_splot_validation_workflow.py",
)
submitter = load_module(
    "psi2s_splot_submitter_test",
    "scripts/submit_psi2s_pbpb_splot_validation.py",
)


class Psi2SSPlotValidationTest(unittest.TestCase):
    def test_contract_has_twelve_variables_and_no_bootstrap(self):
        self.assertEqual(len(workflow.VARIABLES), 12)
        self.assertIn("Btrk2dR", workflow.VARIABLES)
        source = (REPO / "plotER/Validation/psi2s_pbpb_splot_validation_workflow.py").read_text()
        self.assertNotIn("bootstrap_samples", source)

    def test_aggregate_ranks_within_category(self):
        manifest = {"train_tag": "tag"}
        config_path = REPO / "fitER/configs/pbpb24_psi2s_nominal_v2.json"
        result_path = Path("/tmp/result_manifest.json")
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            for index, variable in enumerate(workflow.VARIABLES):
                directory = output / "variables" / variable
                directory.mkdir(parents=True)
                category = "in_model_R6" if index < 6 else "held_out_transfer"
                payload = {
                    "variable": variable, "expression": variable, "category": category,
                    "agreement_score": float(12 - index),
                    "unit_agreement_score": float(13 - index),
                    "delta_discrepancy_unit_minus_weighted": 1.0,
                    "splot_sensitive": variable == "Btrk2dR",
                    "weighted_mass_pearson": 0.0, "weighted_mass_spearman": 0.0,
                    "weighted_mass_slice_max_l1": 0.0,
                    "weighted_mass_slice_max_cdf": 0.0,
                }
                (directory / "metrics.json").write_text(json.dumps(payload))
            manifest_file = output / "input.json"
            manifest_file.write_text("{}")
            result_path = output / "source_result.json"
            result_path.write_text("{}")
            fake = (manifest, {}, config_path, {}, result_path, {})
            with mock.patch.object(workflow, "load_inputs", return_value=fake):
                workflow.aggregate(REPO, manifest_file, output)
            rows = json.loads((output / "agreement_summary.json").read_text())
            self.assertEqual(len(rows), 12)
            self.assertEqual(rows[0]["rank_within_category"], 1)
            self.assertEqual(rows[6]["rank_within_category"], 1)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "passed_with_warnings")
            self.assertFalse(validation["bootstrap"])

    def test_dag_topology_has_one_prepare_splot_and_twelve_analyzers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "input.json"
            manifest.write_text("{}")
            fake = ({"train_tag": "safe_tag"}, {}, Path("config"), {},
                    Path("result"), {"workspace": "workspace"})
            with mock.patch.object(submitter, "AFS_ROOT", root / "afs"), \
                    mock.patch.object(submitter, "EOS_RESULTS_ROOT", root / "eos"), \
                    mock.patch.object(submitter.workflow, "load_inputs", return_value=fake):
                task = submitter.create_submission(manifest, "test")
            dag = Path(task["submission_dir"], "validation.dag").read_text()
            self.assertEqual(dag.count("JOB ANALYZE_"), 12)
            self.assertIn("PARENT PREPARE CHILD SPLOT", dag)
            self.assertIn("PARENT SPLOT CHILD ANALYZE_", dag)
            self.assertNotIn("BOOTSTRAP", dag)


if __name__ == "__main__":
    unittest.main()
