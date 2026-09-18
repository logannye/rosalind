import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("reproduce_demo", Path(__file__).with_name("reproduce_demo.py"))
DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEMO)


class DemonstrationGuards(unittest.TestCase):
    def binary(self, root):
        executable = root / "selected-binary"
        executable.write_text(f"#!{sys.executable}\nimport sys\nassert sys.argv[1:] == ['--version']\nprint('rosalind 0.5.0')\n")
        executable.chmod(0o755)
        return executable

    def test_wrong_version_cannot_create_a_demo_or_fall_back_to_another_binary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as refusal:
                DEMO.main(["--binary", str(self.binary(root)), "--expected-version", "0.5.0-rc.2",
                           "--output", str(root / "demo")])
            self.assertEqual(refusal.exception.code, 2)
            self.assertFalse((root / "demo").exists())

    def test_existing_output_is_never_reused_or_deleted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "existing"
            output.mkdir()
            sentinel = output / "keep-me"
            sentinel.write_text("prior work")
            with patch.dict(sys.modules, {"pysam": SimpleNamespace(__version__="0.23.3")}), \
                    self.assertRaises(FileExistsError):
                DEMO.main(["--binary", str(self.binary(root)), "--expected-version", "0.5.0",
                           "--output", str(output)])
            self.assertEqual(sentinel.read_text(), "prior work")
            self.assertEqual(list(output.iterdir()), [sentinel])

    def test_transcript_hides_machine_paths_and_preserves_public_hashes(self):
        raw = "/machine/demo/portable/a /machine/bin/rosalind /machine/demo abcdef123"
        sanitized = DEMO.normalize(raw, {"/machine/demo": "$DEMO", "/machine/bin/rosalind": "$ROSALIND"})
        self.assertEqual(sanitized, "$DEMO/portable/a $ROSALIND $DEMO abcdef123")


if __name__ == "__main__":
    unittest.main()
