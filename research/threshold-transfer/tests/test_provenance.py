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


class ExecutionWindowTests(unittest.TestCase):
    def fixture(self, root):
        protocol = root / "protocol.json"
        protocol.write_bytes(encoded({"deadline_utc": "2020-01-01T00:00:00Z"}))
        window = runner.read(HERE / "execution-window.json")
        window.update(started_utc="2020-01-02T00:00:00Z", deadline_utc="2020-01-02T03:00:00Z",
                      schema="rosalind-threshold-resumption-v1", previous_hosted_attempts=2,
                      additional_hosted_attempts=1, previous_runner_seconds=275,
                      maximum_additional_runner_seconds=10800,
                      original_deadline_utc="2020-01-01T00:00:00Z",
                      frozen_protocol_sha256=runner.digest(protocol))
        path = root / "execution-window.json"
        path.write_bytes(encoded(window))
        return protocol, path, window

    def test_resumption_leaves_original_scientific_bytes_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            protocol, path, window = self.fixture(Path(directory))
            original = protocol.read_bytes()
            self.assertEqual(runner.execution_deadline(protocol, path).isoformat(), "2020-01-02T03:00:00+00:00")
            self.assertEqual(protocol.read_bytes(), original)
            path.unlink()
            self.assertEqual(runner.execution_deadline(protocol, path).isoformat(), "2020-01-01T00:00:00+00:00")

    def test_changed_protocol_and_excess_time_attempts_or_compute_fail_closed(self):
        changes = {"frozen_protocol_sha256": "0" * 64, "deadline_utc": "2020-01-02T03:00:01Z",
                   "additional_hosted_attempts": 2, "previous_hosted_attempts": 1,
                   "maximum_additional_runner_seconds": 10801, "previous_runner_seconds": 21600,
                   "maximum_total_runner_seconds": 21601}
        for key, value in changes.items():
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                protocol, path, window = self.fixture(Path(directory))
                window[key] = value
                path.write_bytes(encoded(window))
                with self.assertRaisesRegex(ValueError, "Operational resumption differs"):
                    runner.execution_deadline(protocol, path)

    def test_new_input_continuation_retains_overall_compute_ceiling(self):
        with tempfile.TemporaryDirectory() as directory:
            protocol, path, window = self.fixture(Path(directory))
            window.update(schema="rosalind-threshold-preflight-continuation-v1",
                          previous_hosted_attempts=3, additional_hosted_attempts=2,
                          previous_runner_seconds=618, maximum_additional_runner_seconds=18000,
                          deadline_utc="2020-01-02T04:00:00Z")
            path.write_bytes(encoded(window))
            self.assertEqual(runner.execution_deadline(protocol, path).isoformat(), "2020-01-02T04:00:00+00:00")
            window["previous_runner_seconds"] = 3601
            path.write_bytes(encoded(window))
            with self.assertRaisesRegex(ValueError, "remaining budget"):
                runner.execution_deadline(protocol, path)


if __name__ == "__main__":
    unittest.main()
