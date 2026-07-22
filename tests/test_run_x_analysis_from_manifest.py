import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("adapter", ROOT / "scripts/run_x_analysis_from_manifest.py")
adapter = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = adapter
spec.loader.exec_module(adapter)


def model():
    return {"profile": "X_test", "channel": "X",
        "data": {"path": "/tmp/DATA.root", "tree": "ntmix"},
        "mc": {"path": "/tmp/MC.root", "tree": "ntmix_X3872"},
        "score_branch": "xgb_score", "pre_cut": "BQvalue < 0.15 && Bpt > 7.5"}


class AdapterTest(unittest.TestCase):
    def test_native_contract(self):
        adapter.validate_native_x_contract(model(), False)
        bad = model(); bad["mc"] = {"path": "/tmp/MC.root", "tree": "wrong"}
        with self.assertRaisesRegex(adapter.AdapterError, "ntmix_X3872"):
            adapter.validate_native_x_contract(bad, False)

    def test_optimal_exact_replacements(self):
        source = adapter.OPTIMAL_TEMPLATE.read_text()
        patched, _ = adapter.patch_optimal_template(source, model())
        self.assertIn('TString dataPath = "/tmp/DATA.root";', patched)
        self.assertIn('const TString baseCut = "BQvalue < 0.15 && Bpt > 7.5";', patched)
        self.assertEqual(patched.count("xgb_score"), 4)
        self.assertNotIn("Prediction", patched)

    def test_fit_leaves_cuts_unchanged(self):
        source = adapter.FIT_TEMPLATE.read_text()
        patched, _ = adapter.patch_fit_template(source, model())
        self.assertIn('MC="/tmp/MC.root"', patched)
        self.assertEqual([x for x in source.splitlines() if x.startswith("CUTs")],
                         [x for x in patched.splitlines() if x.startswith("CUTs")])

    def test_template_drift_fails_closed(self):
        source = adapter.OPTIMAL_TEMPLATE.read_text().replace("TString dataPath", "TString renamed")
        with self.assertRaisesRegex(adapter.AdapterError, "expected exactly once"):
            adapter.patch_optimal_template(source, model())

    def test_unsafe_value_rejected(self):
        bad = model(); bad["score_branch"] = 'score";bad'
        with self.assertRaisesRegex(adapter.AdapterError, "unsafe"):
            adapter.patch_optimal_template(adapter.OPTIMAL_TEMPLATE.read_text(), bad)

    def test_parse_native_optimal_results(self):
        output = "\n".join([
            "ppRef 7.5 < p_{T} [GeV/c] < 50.0 (incl.): best threshold = 0.58, FOM = 0.002580",
            "ppRef 7.5 < p_{T} [GeV/c] < 12.5: best threshold = 0.24, FOM = 0.002899",
            "ppRef 12.5 < p_{T} [GeV/c] < 17.5: best threshold = 0.38, FOM = 0.004206",
            "ppRef 17.5 < p_{T} [GeV/c] < 22.5: best threshold = 0.44, FOM = 0.006895",
            "ppRef 22.5 < p_{T} [GeV/c] < 50.0: best threshold = 0.00, FOM = 0.007713",
        ])
        results = adapter.parse_optimal_results(output)
        self.assertEqual([row["cut"] for row in results], [0.58, 0.24, 0.38, 0.44, 0.0])

    def test_parse_rejects_incomplete_output(self):
        with self.assertRaisesRegex(adapter.AdapterError, "expected 5"):
            adapter.parse_optimal_results("ppRef only: best threshold = 0.58, FOM = 0.1")

    def test_fit_cuts_come_from_optimal_results(self):
        results = [{"cut": value} for value in (0.58, 0.24, 0.38, 0.44, 0.0)]
        inclusive, binned = adapter.fit_cut_lines(results, "Prediction", "BQvalue < 0.15")
        self.assertEqual(inclusive, 'CUTs_INC="(Prediction > 0.58) && BQvalue < 0.15"')
        self.assertIn("Bpt > 22.5 && Bpt < 50 && Prediction > 0.00", binned)
        source = adapter.FIT_TEMPLATE.read_text()
        patched = adapter.patch_fit_cuts(source, inclusive, binned)
        self.assertIn(inclusive, patched)
        self.assertIn(binned, patched)


if __name__ == "__main__":
    unittest.main()
