"""Physical projection through native extraction, persisted cache, and Python."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pysam

from rosalind import iter_evidence, materialize_evidence, panel_qc


class EvidenceProjectionTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.binary = Path(os.environ['ROSALIND_BINARY']).resolve()
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>chr1\n' + 'A' * 131073 + '\n')
        pysam.faidx(str(self.reference))
        self.bam = self.root / 'reads.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 131073}]}
        starts = [0, 16370, 32750, 49000, 65000, 81910, 98290, 114670, 131040]
        with pysam.AlignmentFile(str(self.bam), 'wb', header=header) as out:
            for position in starts:
                for count in range(20):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = f'r{position}-{count}'
                    read.query_sequence = ('A' if count < 12 else 'C') * 30
                    read.query_qualities = [30 + count % 10] * 30
                    read.flag = 16 if count % 2 else 0
                    read.reference_id = 0
                    read.reference_start = position
                    read.mapping_quality = 60
                    read.cigarstring = '30M'
                    out.write(read)
        pysam.index(str(self.bam))
        self.bed = self.root / 'sites.bed'
        selected = sorted({p for start in starts for p in (start, start + 2, start + 29, start + 31)})
        self.bed.write_text(''.join(f'chr1\t{p}\t{p+1}\tlocus-{p}\n' for p in selected))
        self.options = dict(regions=self.bed, binary=self.binary)

    def table(self, result):
        with pa.ipc.open_stream(result.artifact_path) as stream:
            return stream.read_all()

    def test_projection_survives_workers_budgets_cache_replay_and_diff(self):
        full = materialize_evidence(self.reference, self.bam, self.root / 'full.arrow', **self.options)
        baseline = None
        for workers, budget, tile in [(1, 512, 1), (2, 768, 137), (8, 1024, 16384)]:
            cache = self.root / f'cache-{workers}'
            result = materialize_evidence(self.reference, self.bam, self.root / f'projected-{workers}.arrow',
                fields=['depths', 'alleles'], workers=workers, memory_budget_mb=budget, tile_bases=tile,
                cache_dir=cache, **self.options)
            table = self.table(result)
            self.assertNotIn('base_quality_histogram', table.column_names)
            self.assertTrue(table.equals(self.table(full).select(table.column_names), check_metadata=False))
            data = result.artifact_path.read_bytes()
            if baseline is None:
                baseline = data
            self.assertEqual(data, baseline)
            self.assertLess(len(data), full.artifact_path.stat().st_size)
            resumed = materialize_evidence(self.reference, self.bam, self.root / f'resumed-{workers}.arrow',
                fields='depths,alleles', cache_dir=cache, resume=True, workers=workers, **self.options)
            self.assertEqual(resumed.artifact_path.read_bytes(), baseline)
            manifest = json.loads(resumed.manifest_path.read_text())
            self.assertEqual(manifest['params']['evidence.fields'], '3')
            self.assertEqual(manifest['params']['evidence.schema'], '2')
            self.assertEqual(manifest['measurements']['execution.computed_partitions'], '0')
            subprocess.run([str(self.binary), 'reproduce', '--manifest', str(resumed.manifest_path),
                '--inputs', str(self.root), '--binary', str(self.binary), '--no-attest'], check=True, capture_output=True)
        delta = self.root / 'same-loci.tsv'
        compared = subprocess.run([str(self.binary), 'diff', str(result.manifest_path),
            str(resumed.manifest_path), '--loci-output', str(delta)], capture_output=True)
        self.assertIn(compared.returncode, (0, 1), compared.stderr.decode())
        self.assertEqual(len(delta.read_text().splitlines()), 1)
        mismatch = self.root / 'mismatched-loci.tsv'
        refused = subprocess.run([str(self.binary), 'diff', str(full.manifest_path),
            str(result.manifest_path), '--loci-output', str(mismatch)], capture_output=True)
        self.assertEqual(refused.returncode, 3, refused.stderr.decode())
        self.assertFalse(mismatch.exists())
        with iter_evidence(self.reference, self.bam, fields=['depths'], **self.options) as run:
            self.assertGreater(sum(batch.num_rows for batch in run), 0)
            self.assertIsNotNone(run.result)

    def test_panel_projection_is_default_and_matches_explicit_full(self):
        results = []
        plans = []
        for name, fields in [('default', None), ('full', 'all')]:
            result = panel_qc(self.bam, self.bed, self.root / f'panel-{name}.tsv',
                reference=self.reference, binary=self.binary, fields=fields)
            results.append(result.artifact_path.read_bytes())
            command = [str(self.binary), 'analyze', 'panel-qc', '--alignments', str(self.bam),
                       '--reference', str(self.reference), '--regions', str(self.bed), '--plan']
            if fields:
                command += ['--fields', fields]
            plans.append(json.loads(subprocess.check_output(command)))
        self.assertEqual(results[0], results[1])
        self.assertEqual(plans[0]['fields'], 9)
        self.assertEqual(plans[1]['fields'], 63)
        self.assertLess(plans[0]['bytes_per_locus'], plans[1]['bytes_per_locus'])
        for fields in [['alleles'], ['unknown']]:
            with self.assertRaises(ValueError):
                panel_qc(self.bam, self.bed, self.root / 'invalid.tsv', binary=self.binary, fields=fields)


if __name__ == '__main__':
    unittest.main()
