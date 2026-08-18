import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.test_submit_psi2s_nominal_fit_manifest import valid_manifest


REPO = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, REPO / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


WORKFLOW = load_module(
    "psi2s_data_gaussian_workflow", "fitER/psi2s_data_gaussian_workflow.py"
)
SUBMITTER = load_module(
    "submit_psi2s_data_gaussian", "scripts/submit_psi2s_data_gaussian_manifest.py"
)


def make_inputs(root):
    manifest = valid_manifest()
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    cache_dir = root / "source_cache"
    cache_dir.mkdir()
    data_cache = cache_dir / "DATA_fit_cache.root"
    data_cache.write_bytes(b"compact-root-cache")
    broadest = manifest["working_points"][-1]
    metadata = {
        "source_data": str((root / "DATA.root").resolve()),
        "source_data_tree": "ntmix",
        "selection": (
            f"({broadest['selection']}) && Bmass>3.600000 && Bmass<3.800000"
        ),
        "data_cache": str(data_cache.resolve()),
        "data_entries": broadest["fiducial_score_selected_data_entries"],
        "data_branches": ["Bmass", "Prediction", "source_entry"],
    }
    (cache_dir / "cache_metadata.json").write_text(json.dumps(metadata))
    return manifest_path, cache_dir, manifest


class Psi2SDataGaussianTest(unittest.TestCase):
    def test_candidate_contract_is_exact(self):
        contract = WORKFLOW.candidate_contract(valid_manifest())
        self.assertEqual(contract["signal"]["mean_initial_gev"], 3.686097)
        self.assertEqual(contract["signal"]["mean_range_gev"], [3.681097, 3.691097])
        self.assertEqual(contract["signal"]["sigma_range_gev"], [0.001, 0.008])
        self.assertEqual(contract["background"]["order"], 2)

    def test_cache_validation_is_strict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path, cache_dir, manifest = make_inputs(root)
            WORKFLOW.validate_cache(manifest_path, manifest, cache_dir)
            metadata_path = cache_dir / "cache_metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["source_data"] = "/wrong/DATA.root"
            metadata_path.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(RuntimeError, "source DATA"):
                WORKFLOW.validate_cache(manifest_path, manifest, cache_dir)

    def test_dag_reuses_cache_without_source_prepare(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path, cache_dir, _ = make_inputs(root)
            with mock.patch.object(SUBMITTER, "AFS_ROOT", root / "afs"), \
                 mock.patch.object(SUBMITTER, "EOS_RESULTS_ROOT", root / "eos"):
                task = SUBMITTER.create_submission(manifest_path, cache_dir, "test")
            dag = Path(task["submission_dir"], "fit_scan.dag").read_text()
            self.assertIn("JOB CACHE", dag)
            self.assertEqual(dag.count("CATEGORY FIT"), 7)
            self.assertNotIn("source ROOT", dag)
            self.assertIn("source DATA and MC ROOT are not read", task["io_plan"])

    def test_direct_cli_validate_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path, cache_dir, _ = make_inputs(root)
            result = subprocess.run(
                ["python3", str(REPO / "scripts/submit_psi2s_data_gaussian_manifest.py"),
                 str(manifest_path), "--cache-dir", str(cache_dir), "--validate-only"],
                cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn('"status": "valid"', result.stdout)

    def test_aggregate_marks_sigma_boundary_as_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            manifest = valid_manifest()
            manifest_path = output / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            cache_dir = output / "cache"
            cache_dir.mkdir()
            (output / "run_context.json").write_text(json.dumps({
                "reused_cache": str(cache_dir / "DATA_fit_cache.root"),
                "reused_cache_sha256": "abc",
                "reused_cache_metadata": str(cache_dir / "cache_metadata.json"),
                "root_version": "6.32.02",
            }))
            for index, point in enumerate(manifest["working_points"]):
                point_dir = output / point["key"]
                point_dir.mkdir()
                flags = ["sigma"] if index == 0 else []
                result = {
                    "fit_status": 0, "cov_qual": 3, "edm": 1.e-5,
                    "null_fit_status": 0, "null_cov_qual": 3, "null_edm": 1.e-5,
                    "signal_yield": 10.0,
                    "signal_yield_error": 3.0, "local_significance": 2.0,
                    "mean": 3.686, "sigma": 0.004, "chi2_ndf": 1.0,
                    "chebyshev_a0": 0.1, "chebyshev_a1": -0.1,
                    "signal_over_background": 0.1,
                    "signal_over_sqrt_signal_plus_background": 1.0,
                    "parameter_boundary_flags": flags,
                }
                row = {
                    "point": point["key"],
                    "target_weighted_efficiency": point["target_weighted_efficiency"],
                    "achieved_weighted_efficiency": point["achieved_weighted_efficiency"],
                    "threshold": point["threshold"], "data_entries": 100 + index,
                    "yield": 10.0, "yield_error": 3.0,
                    "local_significance": 2.0, "fit_status": 0,
                    "covQual": 3, "EDM": 1.e-5, "mean": 3.686,
                    "sigma": 0.004, "chi2_ndf": 1.0,
                    "parameter_boundary_flags": flags,
                    "artifact_paths": {}, "fit_result": result,
                }
                (point_dir / "point_manifest.json").write_text(json.dumps(row))
            with mock.patch.object(WORKFLOW.subprocess, "check_output", return_value="test\n"):
                WORKFLOW.aggregate(REPO, manifest_path, manifest, cache_dir, output)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "passed_with_warnings")
            self.assertIn("sigma", validation["warnings"][0])


if __name__ == "__main__":
    unittest.main()
