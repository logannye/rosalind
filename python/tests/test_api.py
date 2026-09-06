import importlib
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
from pathlib import Path

from rosalind.features import FeatureRun, _normalized_version


class FeatureRunLifecycleTest(unittest.TestCase):
    def test_rc_versions_match_without_accepting_other_candidates(self):
        self.assertEqual(_normalized_version("0.4.0-rc.1"), "0.4.0rc1")
        self.assertEqual(_normalized_version("0.4.0"), "0.4.0")
        self.assertNotEqual(_normalized_version("0.4.0-rc.1"), _normalized_version("0.4.0-rc.2"))
        self.assertNotEqual(_normalized_version("0.4.0-rc.1"), _normalized_version("0.4.0"))

    def test_successful_consumption_closes_reader_and_stdout_even_when_run_is_retained(self):
        class EmptyReader:
            closed = False
            def __iter__(self):
                return iter(())
            def close(self):
                self.closed = True
        reader = EmptyReader()
        arrow = SimpleNamespace(ipc=SimpleNamespace(open_stream=lambda _: reader))
        process = subprocess.Popen([sys.executable, "-c", "pass"], stdout=subprocess.PIPE)
        with patch.object(importlib.import_module("rosalind.features"), "_require_pyarrow", return_value=arrow):
            with FeatureRun(process, Path("unused.manifest.json")) as run:
                self.assertEqual(list(run), [])
                self.assertEqual(run.result.returncode, 0)
                self.assertTrue(process.stdout.closed)
                self.assertTrue(reader.closed)
        # Holding the completed run for receipt inspection must not retain FDs.
        self.assertIsNotNone(run.result)
        self.assertTrue(process.stdout.closed)

    def test_early_close_terminates_the_native_child(self):
        process = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            stdout=subprocess.PIPE,
        )
        run = FeatureRun(process, Path(tempfile.gettempdir()) / "unused.manifest.json")
        run.close()
        self.assertIsNotNone(process.poll())


if __name__ == "__main__":
    unittest.main()
