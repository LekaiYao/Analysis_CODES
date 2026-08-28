import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "scripts/submit_x_two_year_fit_strategy_comparison.py"
sys.path.insert(0, str(REPO / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "submit_x_two_year_fit_strategy_comparison", SOURCE
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ComparisonSubmissionTest(unittest.TestCase):
    def test_explicitly_disables_scratch_output_transfer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            task = {
                "input_manifest": "/eos/input.json",
                "cache_root": "/eos/cache",
                "simultaneous_root": "/eos/simultaneous",
                "output_dir": str(root / "eos-output"),
                "submission_dir": str(root / "afs-submit"),
                "point": "xeff20",
            }
            MODULE.create_submission(task)
            submit = (Path(task["submission_dir"]) / "comparison.sub").read_text()
            self.assertIn('transfer_output_files = ""', submit)
            self.assertIn("should_transfer_files = YES", submit)
            self.assertFalse(Path(task["output_dir"]).exists())


if __name__ == "__main__":
    unittest.main()
