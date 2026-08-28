import importlib.util
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "fitER"))
sys.path.insert(0, str(REPO / "scripts"))
SOURCE = REPO / "fitER/run_x_two_year_fit_strategy_comparison.py"
SPEC = importlib.util.spec_from_file_location("x_two_year_fit_strategy_comparison", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PointSelectionTest(unittest.TestCase):
    def summary(self):
        return [
            {"point": "xeff20", "fit_result": {"Z_approx": 4.7}},
            {"point": "xeff35", "fit_result": {"Z_approx": 9.0}},
            {"point": "xeff40", "fit_result": {"Z_approx": 8.0}},
        ]

    def test_best_excludes_only_failed_points(self):
        selected, excluded = MODULE.select_point_key(
            "best", self.summary(),
            {"failures": ["xeff35: DATA status", "xeff40: DATA EDM"]},
            ("xeff20", "xeff35", "xeff40"),
        )
        self.assertEqual(selected, "xeff20")
        self.assertEqual(set(excluded), {"xeff35", "xeff40"})

    def test_explicit_failed_point_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "requested working point"):
            MODULE.select_point_key(
                "xeff35", self.summary(), {"failures": ["xeff35: DATA status"]},
                ("xeff20", "xeff35", "xeff40"),
            )

    def test_global_failure_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "global failures"):
            MODULE.select_point_key(
                "best", self.summary(), {"failures": ["pb23: thresholds not monotonic"]},
                ("xeff20", "xeff35", "xeff40"),
            )


if __name__ == "__main__":
    unittest.main()
