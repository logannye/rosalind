"""No scientific fitting: fail-closed input provenance and safe extraction."""
import hashlib
import importlib.util
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("calibration_runner", Path(__file__).parents[1] / "run.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class ProvenanceTests(unittest.TestCase):
    def fixture(self, directory):
        source = {"repository": "owner/repo", "run_id": 4, "commit": "abc", "artifact_id": 5}
        files = {
            "records.jsonl": b'{"untrusted":true}\n',
            "study/frozen-models.pkl": b"not a pickle; extraction must never deserialize",
            "identity.json": json.dumps({"git_revision": "abc", "github_run_id": "4"}).encode(),
            "report.json": json.dumps({"status": "completed", "gate": {"status": "stop_no_signal"},
                                       "p1_authorized_by_gate": False}).encode(),
        }
        entry = lambda b: {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
        files["run-manifest.json"] = json.dumps({"files": {k: entry(v) for k, v in files.items()}}).encode()
        # The source pilot explicitly wrote final telemetry after its inventory.
        files["telemetry.json"] = b'{"final":true}\n'
        source["files"] = {k: entry(v) for k, v in files.items()}
        archive = directory / "input.zip"
        with zipfile.ZipFile(archive, "w") as z:
            for k, v in files.items():
                z.writestr("hg002-p0/" + k, v)
            z.writestr("../must-not-extract", b"arbitrary member")
        source["artifact_sha256"] = runner.digest(archive)
        return archive, source

    def test_pinned_extraction_ignores_unconsumed_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive, source = self.fixture(root)
            runner.extract_verified(archive, root / "saved", source)
            self.assertEqual(set(p.relative_to(root / "saved").as_posix()
                                 for p in (root / "saved").rglob("*") if p.is_file()), set(source["files"]))
            self.assertFalse((root / "must-not-extract").exists())

    def test_zip_corruption_rejected_before_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive, source = self.fixture(root)
            with archive.open("ab") as h:
                h.write(b"changed")
            with self.assertRaisesRegex(ValueError, "ZIP digest"):
                runner.extract_verified(archive, root / "saved", source)
            self.assertFalse((root / "saved").exists())

    def test_member_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive, source = self.fixture(root)
            source["files"]["records.jsonl"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                runner.extract_verified(archive, root / "saved", source)

    def test_provider_identity_must_match_successful_run(self):
        source = {"repository": "owner/repo", "run_id": 4, "commit": "abc", "artifact_id": 5,
                  "artifact_sha256": "123"}
        run = {"id": 4, "head_sha": "abc", "conclusion": "success", "repository": {"full_name": "owner/repo"}}
        artifact = {"id": 5, "expired": False, "digest": "sha256:123", "workflow_run": {"id": 4, "head_sha": "abc"}}
        runner.verify_provider(run, artifact, source)
        for key, value in (("conclusion", "failure"), ("head_sha", "changed"), ("id", 6)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                runner.verify_provider(run | {key: value}, artifact, source)
        with self.assertRaises(ValueError):
            runner.verify_provider(run, artifact | {"expired": True}, source)


if __name__ == "__main__":
    unittest.main()
