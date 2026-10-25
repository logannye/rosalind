"""Executable/source substitutions and silent test skips must fail closed."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "case_build_identity", Path(__file__).resolve().parents[1] / "build_identity.py")
IDENTITY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IDENTITY)


class BuildIdentityTests(unittest.TestCase):
    def fixture(self, root):
        binary = root / "explicit-binary"
        binary.write_bytes(b"authored executable bytes, never executed")
        identity = {"schema_version": IDENTITY.SCHEMA,
                    "source_revision": "a" * 40, "source_tree": "b" * 40,
                    "binary_sha256": IDENTITY.sha256(binary),
                    "binary_size_bytes": binary.stat().st_size,
                    "binary_version_assertion": "rosalind 0.5.0",
                    "environment": {"dependencies": {"test": "1"}},
                    "native_build_metadata": {"code_git_sha": "a" * 40},
                    "rustc_verbose": "recorded compiler"}
        path = root / "identity.json"
        path.write_text(json.dumps(identity))
        return binary, path, identity

    def authenticate(self, binary, path, record, *, source=None, native=None):
        with patch.object(IDENTITY, "source_identity", return_value=source or {
                "source_revision": "a" * 40, "source_tree": "b" * 40}), \
             patch.object(IDENTITY, "dependencies", return_value={"test": "1"}), \
             patch.object(IDENTITY, "_command", return_value="rosalind 0.5.0"), \
             patch.object(IDENTITY, "native_metadata", return_value=native or record["native_build_metadata"]), \
             patch.object(IDENTITY, "_validate_native"):
            return IDENTITY.authenticate(binary, path, "a" * 40)

    def test_matching_build_record_accepts_selected_binary(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary, path, record = self.fixture(Path(temporary))
            self.assertEqual(self.authenticate(binary, path, record), record)

    def test_wrong_executable_cannot_hide_behind_matching_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary, path, record = self.fixture(Path(temporary))
            binary.write_bytes(b"another executable")
            with self.assertRaisesRegex(ValueError, "digest"):
                self.authenticate(binary, path, record)

    def test_wrong_source_or_compiled_metadata_refuses(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary, path, record = self.fixture(Path(temporary))
            with self.assertRaisesRegex(ValueError, "source_revision"):
                self.authenticate(binary, path, record, source={"source_revision": "c" * 40})
            with self.assertRaisesRegex(ValueError, "compiled metadata"):
                self.authenticate(binary, path, record, native={"code_git_sha": "c" * 40})

    def test_prefix_and_user_labels_are_not_exact_source_identities(self):
        for value in ("a" * 8, "HEAD", "source-preview", "A" * 40):
            with self.subTest(value=value), self.assertRaises(ValueError):
                IDENTITY._source_sha(value)

    def test_exact_checkout_must_be_clean(self):
        with patch.object(IDENTITY, "_command", side_effect=["a" * 40, " M src/main.rs"]):
            with self.assertRaisesRegex(ValueError, "clean"):
                IDENTITY.source_identity("a" * 40)

    def test_existing_report_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = Path(temporary) / "report.json"
            report.write_text("prior failure")
            with self.assertRaises(FileExistsError):
                IDENTITY._write_new(report, {"passed": True})
            self.assertEqual(report.read_text(), "prior failure")

    def test_skipped_native_case_cannot_make_green_gate(self):
        class SkippedNative(unittest.TestCase):
            @unittest.skip("required native binary absent")
            def test_native(self):
                pass

        with tempfile.TemporaryDirectory() as temporary:
            report = Path(temporary) / "tests.json"
            suite = unittest.defaultTestLoader.loadTestsFromTestCase(SkippedNative)
            with patch.object(unittest.defaultTestLoader, "discover", return_value=suite), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertFalse(IDENTITY.run_test_suite("unused", "test_*.py", report))
            result = json.loads(report.read_text())
            self.assertEqual(result["tests_run"], 1)
            self.assertEqual(len(result["skips"]), 1)
            self.assertFalse(result["passed"])

    def test_empty_suite_cannot_make_green_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(unittest.defaultTestLoader, "discover", return_value=unittest.TestSuite()), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertFalse(IDENTITY.run_test_suite("unused", "test_*.py", Path(temporary) / "tests.json"))


if __name__ == "__main__":
    unittest.main()
