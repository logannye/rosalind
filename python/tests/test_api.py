import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rosalind.features import FeatureRun


class FeatureRunLifecycleTest(unittest.TestCase):
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
