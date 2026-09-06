import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import pyarrow as pa

from rosalind import EvidenceRun, EvidenceProcessError
from rosalind.evidence import _command


class EvidenceLifecycleTest(unittest.TestCase):
    def test_native_stream_is_lazy_bounded_and_finalizes(self):
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "started"
            manifest = Path(directory) / "receipt.json"
            code = (
                "import pathlib,sys,pyarrow as pa; "
                "pathlib.Path(sys.argv[1]).write_text('started'); "
                "schema=pa.schema([('depth',pa.uint64())]); "
                "writer=pa.ipc.new_stream(sys.stdout.buffer,schema); "
                "writer.write_batch(pa.record_batch([pa.array([1]*1024,type=pa.uint64())],schema=schema)); "
                "writer.write_batch(pa.record_batch([pa.array([0],type=pa.uint64())],schema=schema)); "
                "writer.close(); pathlib.Path(sys.argv[2]).write_text('{}')"
            )
            with EvidenceRun([sys.executable, "-c", code, str(marker), str(manifest)], manifest) as run:
                self.assertFalse(marker.exists())
                self.assertIsNone(run.result)
                batches = iter(run)
                self.assertEqual(next(batches).num_rows, 1024)
                self.assertIsNone(run.result)
                self.assertEqual(next(batches).num_rows, 1)
                with self.assertRaises(StopIteration):
                    next(batches)
                self.assertEqual(run.result.returncode, 0)
                self.assertIsNone(run.result.artifact_path)
                self.assertTrue(run.result.manifest_path.is_file())

    def test_native_errors_preserve_exit_code_and_diagnostics(self):
        for code in (2, 3, 4, 5):
            with self.subTest(code=code):
                run = EvidenceRun([sys.executable, "-c", f"import sys; print('precise child failure',file=sys.stderr); sys.exit({code})"], Path("unused.json"))
                with self.assertRaises(EvidenceProcessError) as caught:
                    list(run)
                self.assertEqual(caught.exception.returncode, code)
                self.assertIn("precise child failure", str(caught.exception))
                self.assertIsNone(run.result)
                self.assertIsNotNone(run._process.poll())

    def test_partial_consumption_cancels_child_and_has_no_success_result(self):
        code = (
            "import sys,time,pyarrow as pa; schema=pa.schema([('depth',pa.uint64())]); "
            "writer=pa.ipc.new_stream(sys.stdout.buffer,schema); "
            "writer.write_batch(pa.record_batch([pa.array([1],type=pa.uint64())],schema=schema)); "
            "sys.stdout.buffer.flush(); time.sleep(60)"
        )
        with EvidenceRun([sys.executable, "-c", code], Path("unused.json")) as run:
            iterator = iter(run)
            self.assertEqual(next(iterator).num_rows, 1)
        self.assertIsNotNone(run._process.poll())
        self.assertIsNone(run.result)
        iterator.close()

    def test_malformed_producer_cannot_deadlock_on_its_output_pipe(self):
        code = (
            "import sys,time; "
            "sys.stdout.buffer.write(b'\\xff'*4+(8).to_bytes(4,'little')+b'\\x00'*8); "
            "sys.stdout.buffer.flush(); time.sleep(60)"
        )
        run = EvidenceRun([sys.executable, "-c", code], Path("unused.json"))
        started = time.monotonic()
        with self.assertRaises((pa.ArrowInvalid, pa.ArrowIOError, OSError)):
            list(run)
        self.assertLess(time.monotonic() - started, 10)
        self.assertIsNotNone(run._process.poll())
        self.assertIsNone(run.result)

    def test_exclusive_selection_is_checked_before_starting_native_code(self):
        with self.assertRaises(ValueError):
            _command("reference.fa", "reads.bam")
        with self.assertRaises(ValueError):
            _command("reference.fa", "reads.bam", sites="sites.vcf", regions="panel.bed")
        with self.assertRaises(ValueError):
            _command(None, "reads.bam", regions="panel.bed")


if __name__ == "__main__":
    unittest.main()
