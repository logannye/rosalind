"""Scientific contract parity through the real native CLI, not a mock producer."""

import csv
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import pysam

from rosalind import EvidenceProcessError, materialize_evidence, panel_qc


class EvidenceContractTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.binary = Path(os.environ.get(
            "ROSALIND_BINARY", Path(__file__).resolve().parents[2] / "target/debug/rosalind"
        )).resolve()
        if not self.binary.is_file():
            self.fail("build rosalind and set ROSALIND_BINARY to test the native contract")
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">chr1\n" + "A" * 64 + "\n")
        pysam.faidx(str(self.reference))
        self.bam = self.root / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 64}],
                  "RG": [{"ID": "a", "SM": "alpha"}, {"ID": "b", "SM": "beta"}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as writer:
            for position, depth in enumerate((9, 10, 19, 20)):
                for i in range(depth + 1):
                    read = pysam.AlignedSegment(writer.header)
                    read.query_name = f"r{position}-{i}"
                    read.query_sequence = "A"
                    read.query_qualities = [40]
                    read.reference_id = 0
                    read.reference_start = position
                    read.mapping_quality = 60
                    read.cigarstring = "1M"
                    read.set_tag("RG", "a" if i < depth else "b")
                    writer.write(read)
        pysam.index(str(self.bam))
        self.bed = self.root / "targets.bed"
        self.bed.write_text("".join(f"chr1\t{p}\t{p+1}\tdepth-{d}\n"
                                    for p, d in enumerate((9, 10, 19, 20))))

    def panel(self, name, **options):
        result = panel_qc(self.bam, self.bed, self.root / name,
                          reference=self.reference, binary=self.binary, **options)
        with result.artifact_path.open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        return result, rows

    def test_python_and_cli_defaults_match_at_callability_boundaries(self):
        result, rows = self.panel("python.tsv", sample="alpha")
        self.assertEqual([r["callable_positions"] for r in rows], ["0", "1", "1", "1"])
        self.assertEqual({r["callable_threshold"] for r in rows}, {"10"})
        native = self.root / "native.tsv"
        subprocess.run([str(self.binary), "analyze", "panel-qc", "--alignments", str(self.bam),
                        "--reference", str(self.reference), "--regions", str(self.bed),
                        "--sample", "alpha", "--output", str(native)], check=True,
                       capture_output=True)
        self.assertEqual(result.artifact_path.read_bytes(), native.read_bytes())
        _, explicit = self.panel("explicit.tsv", sample="alpha", min_callable_depth=20)
        self.assertEqual([r["callable_positions"] for r in explicit], ["0", "0", "0", "1"])

    def test_multiple_samples_require_explicit_scope_and_cache_keeps_it(self):
        with self.assertRaises(EvidenceProcessError) as caught:
            self.panel("implicit.tsv")
        self.assertEqual(caught.exception.returncode, 2)
        self.assertFalse((self.root / "implicit.tsv").exists())
        cache = self.root / "cache"
        alpha, alpha_rows = self.panel("alpha.tsv", sample="alpha", cache_dir=cache)
        beta, beta_rows = self.panel("beta.tsv", sample="beta", cache_dir=cache, resume=True)
        pooled, pooled_rows = self.panel("pooled.tsv", pool_samples=True)
        self.assertEqual([int(r["callable_depth_sum"]) for r in alpha_rows], [9, 10, 19, 20])
        self.assertEqual([int(r["callable_depth_sum"]) for r in beta_rows], [1, 1, 1, 1])
        self.assertEqual([int(r["callable_depth_sum"]) for r in pooled_rows], [10, 11, 20, 21])
        receipts = [json.loads(r.manifest_path.read_text()) for r in (alpha, beta, pooled)]
        self.assertEqual(len({r["params"]["evidence.science_blake3"] for r in receipts}), 3)
        for receipt in receipts:
            self.assertIn("evidence.sample_scope", receipt["params"])
        subprocess.run([str(self.binary), "reproduce", "--manifest", str(alpha.manifest_path),
                        "--inputs", str(self.root), "--binary", str(self.binary), "--no-attest"],
                       check=True, capture_output=True)

    def test_conflicting_sample_options_fail_before_execution(self):
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            materialize_evidence(self.reference, self.bam, self.root / "conflict.arrow",
                                 regions=self.bed, binary=self.binary,
                                 sample="alpha", pool_samples=True)
        self.assertFalse((self.root / "conflict.arrow").exists())


if __name__ == "__main__":
    unittest.main()
