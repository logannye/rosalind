"""Pure stdlib artifact, mutation, and admission checks; no native execution."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from contract import Candidate
from planner import MAX_FILES, MIB, admit, windows
from run_manifest import (SCHEMA, InputLedger, digest_json, file_digest, load_json,
                          output_inventory, publish_new, verify_run, write_json)


class TemporaryFiles(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()

    def make_file(self, name, data=b"evidence"):
        path = self.root / name
        path.write_bytes(data)
        return path


class InputLedgerTests(TemporaryFiles):
    def test_nonregular_inputs_are_refused_before_open(self):
        paths = [self.root]
        if Path("/dev/null").exists():
            paths.append(Path("/dev/null"))
        for path in paths:
            with self.subTest(path=path), \
                    patch.object(Path, "open", side_effect=AssertionError("must not open")) as opened:
                with self.assertRaisesRegex(ValueError, "must be a regular file"):
                    file_digest(path)
                opened.assert_not_called()

    @unittest.skipUnless(hasattr(os, "mkfifo"), "named pipes require POSIX")
    def test_fifo_is_refused_without_waiting_for_a_writer(self):
        fifo = self.root / "candidate-input.fifo"
        os.mkfifo(fifo)
        ledger = InputLedger()
        # An attempted open fails immediately even if the regular-file guard
        # regresses; this test can never hang waiting for a FIFO writer.
        with patch.object(Path, "open", side_effect=AssertionError("must not open FIFO")) as opened:
            with self.assertRaisesRegex(ValueError, "must be a regular file"):
                ledger.pin(fifo)
            opened.assert_not_called()
        self.assertEqual(ledger.as_list(), [])
        self.assertEqual(ledger.bytes_hashed, 0)

    def test_final_union_detects_mutation_of_an_earlier_consumed_child(self):
        first = self.make_file("earlier.arrow")
        last = self.make_file("last.arrow", b"last-child")
        ledger = InputLedger()
        ledger.pin(first)
        ledger.pin(last)
        first.write_bytes(b"modified")
        with self.assertRaisesRegex(ValueError, "changed before publication"):
            ledger.recheck()

    def test_repeated_consumption_refuses_changed_identity(self):
        path = self.make_file("partition.arrow")
        ledger = InputLedger()
        ledger.pin(path)
        path.write_bytes(b"new-data")
        with self.assertRaisesRegex(ValueError, "repeated input identity changed"):
            ledger.pin(path)

    def test_same_bytes_at_retargeted_symlink_are_still_refused(self):
        first = self.make_file("first.arrow")
        second = self.make_file("second.arrow")
        access = self.root / "access.arrow"
        access.symlink_to(first)
        ledger = InputLedger()
        ledger.pin(access)
        access.unlink()
        access.symlink_to(second)
        with self.assertRaisesRegex(ValueError, "path retargeted"):
            ledger.recheck()

    def test_expected_digest_and_inventory_cap(self):
        first = self.make_file("first.arrow")
        second = self.make_file("second.arrow")
        ledger = InputLedger(max_files=1)
        with self.assertRaisesRegex(ValueError, "verified input identity changed"):
            ledger.pin(first, "0" * 64)
        self.assertEqual(ledger.as_list(), [])
        expected, size = file_digest(first)
        ledger.pin(first, expected)
        ledger.pin(first, expected)
        self.assertEqual(len(ledger.as_list()), 1)
        with self.assertRaisesRegex(ValueError, "inventory limit"):
            ledger.pin(second)
        ledger.recheck()
        self.assertEqual(ledger.bytes_hashed, 4 * size)


class PublicationAndVerificationTests(TemporaryFiles):
    def create_run(self):
        directory = self.root / "run"
        directory.mkdir()
        write_json(directory / "scientific.json", {"schema": "rosalind-cohort-scientific-v1", "test": 1})
        write_json(directory / "plan.json", {"verification_level": "metadata_only"})
        (directory / "summary.jsonl").write_text('{"n_observed":0}\n')
        (directory / "coverage.jsonl").write_text('{"missing_loci":1}\n')
        outputs, size = output_inventory(directory, 10, 1024)
        manifest = {"schema": SCHEMA, "status": "completed", "outputs": outputs,
                    "scientific_sha256": digest_json(load_json(directory / "scientific.json")),
                    "retained_bytes_excluding_manifest": size}
        self.write_manifest(directory, manifest)
        return directory

    def write_manifest(self, directory, manifest):
        manifest = dict(manifest)
        manifest.pop("manifest_sha256", None)
        manifest["manifest_sha256"] = digest_json(manifest)
        path = directory / "run-manifest.json"
        path.unlink(missing_ok=True)
        write_json(path, manifest)

    @unittest.skipUnless(sys.platform in ("linux", "darwin"), "atomic publication supports Linux/macOS")
    def test_atomic_publication_preserves_existing_destination(self):
        staging, destination = self.root / ".staging", self.root / "published"
        staging.mkdir()
        destination.mkdir()
        (staging / "new").write_bytes(b"new")
        (destination / "existing").write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            publish_new(staging, destination)
        self.assertEqual((destination / "existing").read_bytes(), b"keep")
        self.assertEqual((staging / "new").read_bytes(), b"new")
        empty_destination = self.root / "existing-empty"
        empty_destination.mkdir()
        with self.assertRaises(FileExistsError):
            publish_new(staging, empty_destination)
        fresh = self.root / "fresh"
        publish_new(staging, fresh)
        self.assertFalse(staging.exists())
        self.assertEqual((fresh / "new").read_bytes(), b"new")

    def test_verify_requires_completed_manifest_and_unchanged_outputs(self):
        run = self.create_run()
        self.assertEqual(verify_run(run)["status"], "completed")
        (run / "summary.jsonl").write_text('{"n_observed":1}\n')
        with self.assertRaisesRegex(ValueError, "output identity mismatch"):
            verify_run(run)
        manifest = load_json(run / "run-manifest.json")
        manifest["status"] = "failed"
        self.write_manifest(run, manifest)
        with self.assertRaisesRegex(ValueError, "not a completed"):
            verify_run(run)

    def test_verify_rejects_unbound_manifest_edit_and_missing_required_output(self):
        run = self.create_run()
        manifest = load_json(run / "run-manifest.json")
        manifest["status"] = "planned"
        (run / "run-manifest.json").write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            verify_run(run)
        manifest["status"] = "completed"
        manifest["outputs"] = [entry for entry in manifest["outputs"] if entry["path"] != "summary.jsonl"]
        self.write_manifest(run, manifest)
        with self.assertRaisesRegex(ValueError, "lacks required outputs"):
            verify_run(run)

    def test_verify_checks_scientific_identity_against_retained_scientific_document(self):
        run = self.create_run()
        manifest = load_json(run / "run-manifest.json")
        manifest["scientific_sha256"] = "0" * 64
        self.write_manifest(run, manifest)
        with self.assertRaisesRegex(ValueError, "scientific.*identity|science.*identity"):
            verify_run(run)

    def test_verify_rejects_traversal_and_duplicate_inventory_paths(self):
        run = self.create_run()
        original = load_json(run / "run-manifest.json")
        for path in ("../outside.json", "/absolute.json", original["outputs"][1]["path"]):
            manifest = json.loads(json.dumps(original))
            manifest["outputs"][0]["path"] = path
            self.write_manifest(run, manifest)
            with self.subTest(path=path), self.assertRaises(ValueError):
                verify_run(run)

    def test_inventory_limits_and_symlinks_fail(self):
        directory = self.root / "files"
        directory.mkdir()
        (directory / "first").write_bytes(b"1234")
        with self.assertRaisesRegex(ValueError, "output envelope"):
            output_inventory(directory, max_files=1, max_bytes=3)
        (directory / "second").write_bytes(b"5678")
        with self.assertRaisesRegex(ValueError, "output envelope"):
            output_inventory(directory, max_files=1, max_bytes=100)
        (directory / "link").symlink_to(directory / "first")
        with self.assertRaisesRegex(ValueError, "symlink"):
            output_inventory(directory, max_files=10, max_bytes=100)
        (directory / "link").unlink()
        subdirectory = self.root / "outside-directory"
        subdirectory.mkdir()
        (directory / "directory-link").symlink_to(subdirectory, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink"):
            output_inventory(directory, max_files=10, max_bytes=100)

    def test_metadata_parser_is_byte_bounded_and_write_is_create_new(self):
        path = self.root / "metadata.json"
        write_json(path, {"test": 1})
        with self.assertRaisesRegex(ValueError, "metadata exceeds"):
            load_json(path, limit=3)
        with self.assertRaises(FileExistsError):
            write_json(path, {"test": 2})
        self.assertEqual(load_json(path), {"test": 1})


class AdmissionTests(unittest.TestCase):
    def admission(self, **changes):
        parameters = dict(parent_mb=128, metadata_bytes=4096, candidate_rows=4,
                          thresholds=[5, 10, 20], member_count=3, window_size=2,
                          unique_loci=4, inventory_files=12)
        parameters.update(changes)
        return admit(**parameters)

    def test_admitted_outputs_and_launches_have_explicit_bounds(self):
        result = self.admission()
        self.assertEqual(result["native_child_invocations_upper_bound"], 6)
        self.assertEqual(result["output_rows"], 12)
        self.assertLessEqual(result["modeled_parent_bytes"], 128 * MIB)
        self.assertIn("not an OS RSS cap", result["assurance"])

    def test_oversized_metadata_and_child_inventory_are_refused(self):
        for changes in (dict(parent_mb=64), dict(metadata_bytes=10 * MIB),
                        dict(candidate_rows=100001), dict(member_count=129),
                        dict(inventory_files=MAX_FILES + 1), dict(window_size=0),
                        dict(member_count=128, unique_loci=33, window_size=1),
                        dict(thresholds=list(range(1, 10)))):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.admission(**changes)

    def test_locus_windows_never_split_alternate_rows(self):
        candidates = [Candidate(0, 0, "C", "chr", "A"), Candidate(0, 0, "T", "chr", "A"),
                      Candidate(0, 2, "G", "chr", "A")]
        self.assertEqual(list(windows(candidates, 1)), [candidates[:2], candidates[2:]])
        self.assertEqual(list(windows([], 1)), [])


if __name__ == "__main__":
    unittest.main()
