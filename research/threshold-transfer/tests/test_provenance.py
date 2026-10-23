"""Reject malformed saved-model inputs before any executable deserialization."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import warnings
import zipfile

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("threshold_hosted_runner", HERE / "run.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def encoded(value):
    return json.dumps(value, sort_keys=True).encode()


class ModelProvenanceTests(unittest.TestCase):
    def fixture(self, root, *, duplicate=False, bad_completion=False):
        model = b"deliberately not an executable pickle"
        members = {
            "frozen-models.pkl": model,
            "frozen-models.json": encoded({"model_sha256": hashlib.sha256(model).hexdigest()}),
            "identity.json": encoded({"git_revision": "trusted", "github_run_id": "123", "dependencies": {}, "source_sha256": {}}),
            "protocol.json": b"{}",
            "report.json": encoded({"gate": {"status": "inconclusive_validation_screen"}}),
            "input-recheck.json": b"{}",
        }
        def inventory():
            return {name: {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()} for name, raw in members.items()}
        members["run-manifest.json"] = encoded({"files": inventory()})
        members["completion.json"] = encoded({"status": "completed", "input_identity_rechecked": True,
            "resource_checks_passed": True,
            "report_sha256": "bad" if bad_completion else hashlib.sha256(members["report.json"]).hexdigest(),
            "run_manifest_sha256": hashlib.sha256(members["run-manifest.json"]).hexdigest()})
        archive = root / "source.zip"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with zipfile.ZipFile(archive, "w") as bundle:
                for name, raw in members.items():
                    bundle.writestr("calibration/" + name, raw)
                if duplicate:
                    bundle.writestr("calibration/frozen-models.pkl", model)
        source = {"artifact_sha256": runner.digest(archive), "prefix": "calibration/",
                  "files": inventory(), "commit": "trusted", "run_id": 123}
        return archive, source

    def test_wrong_archive_hash_never_deserializes_or_creates_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive, source = self.fixture(root)
            source["artifact_sha256"] = "0" * 64
            with patch.object(runner.pickle, "load", side_effect=AssertionError("must not unpickle")):
                with self.assertRaisesRegex(ValueError, "ZIP digest"):
                    runner.verify_model_archive(archive, root / "output", source)
            self.assertFalse((root / "output").exists())

    def test_duplicate_pinned_member_never_deserializes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive, source = self.fixture(root, duplicate=True)
            with patch.object(runner.pickle, "load", side_effect=AssertionError("must not unpickle")):
                with self.assertRaisesRegex(ValueError, "duplicate"):
                    runner.verify_model_archive(archive, root / "output", source)

    def test_invalid_completion_never_deserializes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive, source = self.fixture(root, bad_completion=True)
            with patch.object(runner.pickle, "load", side_effect=AssertionError("must not unpickle")):
                with self.assertRaisesRegex(ValueError, "completion"):
                    runner.verify_model_archive(archive, root / "output", source)

    def test_existing_member_mutation_never_deserializes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive, source = self.fixture(root)
            output = root / "output"
            output.mkdir()
            (output / "frozen-models.pkl").write_bytes(b"mutated")
            with patch.object(runner.pickle, "load", side_effect=AssertionError("must not unpickle")):
                with self.assertRaisesRegex(ValueError, "Existing trusted"):
                    runner.verify_model_archive(archive, output, source)


if __name__ == "__main__":
    unittest.main()
