import importlib.util
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
    "xeff45_workflow_test",
    "plotER/Validation/psi2s_xeff45_splot_closure_workflow.py",
)
submitter = load_module(
    "xeff45_submitter_test", "scripts/submit_psi2s_xeff45_splot_closure.py"
)


class Xeff45SPlotClosureTest(unittest.TestCase):
    def test_contract_is_one_off_and_strict_gate(self):
        self.assertEqual(workflow.POINT, "psi2seff45")
        self.assertEqual(workflow.NEFF_GATE, 30.0)
        self.assertEqual(len(workflow.VARIABLES), 12)

    def test_dag_stops_children_when_splot_gate_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.json"
            manifest.write_text("{}")
            fake = ({"train_tag": "tag"}, {}, Path("config"), {}, Path("result"),
                    {"workspace": "workspace"})
            with mock.patch.object(submitter, "AFS_ROOT", root / "afs"), \
                    mock.patch.object(submitter, "EOS_ROOT", root / "eos"), \
                    mock.patch.object(submitter.workflow, "load_inputs", return_value=fake):
                task = submitter.create_submission(manifest)
            dag = Path(task["submission_dir"], "validation.dag").read_text()
            self.assertEqual(dag.count("JOB ANALYZE_"), 12)
            self.assertIn("PARENT SPLOT CHILD ANALYZE_", dag)
            self.assertIn("JOB AGGREGATE aggregate.sub", dag)
            self.assertNotIn("FINAL", dag)
            self.assertNotIn("BOOTSTRAP", dag)

    def test_regular_xeff30_defaults_unchanged(self):
        source = (REPO / "scripts/submit_psi2s_pbpb_splot_validation.py").read_text()
        self.assertIn('default="xeff30_nominal_v2"', source)


if __name__ == "__main__":
    unittest.main()
