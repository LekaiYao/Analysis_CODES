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

    def test_aggregate_reports_cdf_without_ranking(self):
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
                    "weighted_10bin": {"cdf": 0.1, "l1": 0.2, "chi2": 9.0, "ndf": 9},
                    "weighted_5bin": {"cdf": 0.2, "l1": 0.3, "chi2": 4.0, "ndf": 4},
                    "unit_10bin": {"cdf": 0.3, "l1": 0.4, "chi2": 9.0, "ndf": 9},
                    "unit_5bin": {"cdf": 0.4, "l1": 0.5, "chi2": 4.0, "ndf": 4},
                }
                (directory / "metrics.json").write_text(json.dumps(payload))
            manifest_file = output / "input.json"
            manifest_file.write_text("{}")
            result_path = output / "source_result.json"
            result_path.write_text("{}")
            splot = output / "splot_corrected_v2"
            splot.mkdir()
            (splot / "sweight_quality.json").write_text(json.dumps({
                "yield_fit_status": 0, "yield_fit_cov_qual": 3,
                "yield_fit_edm": 1e-8, "relative_yield_closure": 1e-14,
                "effective_entries": 26.7,
            }))
            (splot / "sweighted_data.root").touch()
            (output / "cache").mkdir()
            (output / "cache/MC_xeff30.root").touch()
            fake = (manifest, {}, config_path, {}, result_path, {})
            with mock.patch.object(workflow, "load_inputs", return_value=fake):
                workflow.aggregate(
                    REPO, manifest_file, output, output,
                    "splot_corrected_v2", True,
                )
            rows = json.loads((output / "variable_results.json").read_text())
            self.assertEqual(len(rows), 12)
            self.assertNotIn("rank_within_category", rows[0])
            csv_text = (output / "variable_results.csv").read_text()
            self.assertIn("weighted_cdf_10bin", csv_text)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "passed_with_warnings")
            self.assertFalse(validation["bootstrap"])
            self.assertEqual(validation["ranking"], "none")

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

    def test_low_neff_resume_has_no_prepare_splot_or_ranking(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "input.json"
            manifest.write_text("{}")
            source = root / "source"
            source.mkdir()
            fake = ({"train_tag": "safe_tag"}, {}, Path("config"), {},
                    Path("result"), {"workspace": "workspace"})
            quality = ({"effective_entries": 26.7}, Path("sw.root"), Path("mc.root"))
            with mock.patch.object(submitter, "AFS_ROOT", root / "afs"), \
                    mock.patch.object(submitter, "EOS_RESULTS_ROOT", root / "eos"), \
                    mock.patch.object(submitter.workflow, "load_inputs", return_value=fake), \
                    mock.patch.object(submitter.workflow, "validate_sweight_input", return_value=quality):
                task = submitter.create_low_neff_resume_submission(
                    manifest, "exploratory", source
                )
            dag = Path(task["submission_dir"], "validation.dag").read_text()
            self.assertEqual(dag.count("JOB ANALYZE_"), 12)
            self.assertNotIn("PREPARE", dag)
            self.assertNotIn("SPLOT", dag)
            self.assertEqual(task["ranking"], "none")


if __name__ == "__main__":
    unittest.main()
