import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "scripts" / "import_ml_manifest.py"
SPEC = importlib.util.spec_from_file_location("import_ml_manifest", MODULE_PATH)
IMPORTER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(IMPORTER)
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "analysis_manifest.example.json"


class ImportManifestTest(unittest.TestCase):
    def load_raw(self):
        return json.loads(FIXTURE.read_text(encoding="utf-8"))

    def test_validation_resolves_paths_and_converts_boolean_words(self):
        model = IMPORTER.validate_manifest(self.load_raw(), FIXTURE.resolve())
        self.assertEqual(model["profile"], "Example_Bu_pp")
        self.assertEqual(model["pre_cut"], "(abs(By) < 2.4) && !(Bpt <= 7.5)")
        self.assertEqual(model["data"]["path"], str((FIXTURE.parent / "DATA_with_score.root").resolve()))

    def test_duplicate_target_sections_collapse_and_preserve_cuts(self):
        model = IMPORTER.validate_manifest(self.load_raw(), FIXTURE.resolve())
        existing = """# local config
[other]
value=keep

[Example_Bu_pp]
optimalCUT_punzi=0.91
optimalCUT_precut=0.60

[Example_Bu_pp]
stale=yes
optimalCUT_fom=0.82
"""
        updated = IMPORTER.update_config(existing, model)
        self.assertEqual(updated.count("[Example_Bu_pp]"), 1)
        self.assertIn("[other]\nvalue=keep", updated)
        self.assertIn("optimalCUT_punzi=0.91", updated)
        self.assertIn("optimalCUT_fom=0.82", updated)
        self.assertIn("optimalCUT_precut=0.60", updated)
        self.assertNotIn("stale=yes", updated)

    def test_suggestions_are_comments_not_active_physics_keys(self):
        model = IMPORTER.validate_manifest(self.load_raw(), FIXTURE.resolve())
        rendered = IMPORTER.render_profile(model, {})
        self.assertIn("# suggested_mass_range=", rendered)
        self.assertIn("# suggested_signal_region=", rendered)
        self.assertIn("# suggested_bin_width=", rendered)
        active = [line for line in rendered.splitlines() if line and not line.startswith("#")]
        self.assertFalse(any(line.startswith(("mass_range=", "fsRegion=", "bin_width=")) for line in active))

    def test_rejects_missing_required_field(self):
        raw = self.load_raw()
        del raw["profiles"]["selection"]["signal_expression"]
        with self.assertRaises(IMPORTER.ManifestError):
            IMPORTER.validate_manifest(raw, FIXTURE.resolve())

    def test_rejects_bad_schema_and_path_base(self):
        raw = self.load_raw()
        for key, value in (("schema_version", 2), ("path_base", "cwd")):
            bad = copy.deepcopy(raw)
            bad[key] = value
            with self.subTest(key=key), self.assertRaises(IMPORTER.ManifestError):
                IMPORTER.validate_manifest(bad, FIXTURE.resolve())

    def test_end_to_end_update_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "optimalCUT.local.conf"
            IMPORTER.import_manifest(FIXTURE, output)
            first = output.read_text(encoding="utf-8")
            IMPORTER.import_manifest(FIXTURE, output)
            self.assertEqual(output.read_text(encoding="utf-8"), first)


if __name__ == "__main__":
    unittest.main()
