"""Portable native evidence reuse, physical projection, and language interoperability."""
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from rosalind import EvidenceProcessError, materialize_evidence, open_dataset, panel_qc


class PersistedDatasetTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.binary = Path(os.environ['ROSALIND_BINARY']).resolve()
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>chr1\n' + 'A' * 17020 + '\n')
        pysam.faidx(str(self.reference))
        self.bam = self.root / 'reads.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 17020}]}
        with pysam.AlignmentFile(str(self.bam), 'wb', header=header) as out:
            for start in [0, 1018, 16370, 17000]:
                for count in range(12):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = f'{start}-{count}'
                    read.query_sequence = ('A' if count < 8 else 'C') * 20
                    read.query_qualities = [30 + count] * 20
                    read.flag = 16 if count % 2 else 0
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigarstring = '20M'
                    out.write(read)
        pysam.index(str(self.bam))
        self.bed = self.root / 'panel.bed'
        self.bed.write_text('chr1\t0\t17020\tfull\nchr1\t1018\t1030\tinner\n')
        self.subset = self.root / 'subset.bed'
        self.subset.write_text('chr1\t5\t30\tfirst\nchr1\t16370\t16400\tboundary\n')
        self.cache = self.root / 'cache'
        result = materialize_evidence(self.reference, self.bam, self.root / 'full.arrow',
            regions=self.bed, fields=['depths', 'alleles', 'quality-sums', 'allele-quality'],
            cache_dir=self.cache, binary=self.binary)
        receipt = json.loads(result.manifest_path.read_text())
        self.manifest = Path(receipt['measurements']['execution.evidence_dataset_manifest'])
        self.dataset = open_dataset(self.manifest, binary=self.binary)

    def test_portable_subset_matches_fresh_across_budgets_and_replays(self):
        fresh = materialize_evidence(self.reference, self.bam, self.root / 'fresh.arrow',
            regions=self.subset, fields=['depths', 'alleles'], binary=self.binary)
        live_panel = panel_qc(self.bam, self.subset, self.root / 'live-panel.tsv',
            reference=self.reference, binary=self.binary)
        portable = self.root / 'portable'
        shutil.move(str(self.manifest.parent), portable)
        self.dataset = open_dataset(portable / self.manifest.name, binary=self.binary)
        self.reference.unlink()
        self.bam.unlink()
        info = self.dataset.inspect()
        self.assertEqual(info['fields'], 75)
        verified = self.dataset.verify()
        self.assertEqual(verified['verified_loci'], 17020)
        self.assertFalse(verified['original_sources_rehashed'])
        science = []
        for budget in [128, 256, 512]:
            output = self.dataset.materialize(self.root / f'query-{budget}.arrow',
                regions=self.subset, fields=['depths', 'alleles'], memory_budget_mb=budget)
            self.assertEqual(output.artifact_path.read_bytes(), fresh.artifact_path.read_bytes())
            receipt = json.loads(output.manifest_path.read_text())
            science.append(receipt['params']['science.blake3'])
            self.assertEqual(receipt['measurements']['execution.alignment_record_visits'], '0')
            replay = subprocess.run([str(self.binary), 'reproduce', '--manifest', str(output.manifest_path),
                '--inputs', str(self.root), '--binary', str(self.binary), '--no-attest'], capture_output=True)
            self.assertEqual(replay.returncode, 0, replay.stderr.decode() + replay.stdout.decode())
        self.assertEqual(len(set(science)), 1)
        saved_panel = self.dataset.panel_qc(self.subset, self.root / 'saved-panel.tsv')
        self.assertEqual(saved_panel.artifact_path.read_bytes(), live_panel.artifact_path.read_bytes())
        plan = self.dataset.plan(fields=['depths'])
        self.assertEqual(plan['source_fields'], 75)
        self.assertEqual(plan['fields'], 1)
        self.assertGreater(plan['source_decoder_bytes'], plan['projection_bytes'])

    def test_iterator_lifecycle_empty_and_missing_evidence(self):
        workdir = self.root / 'iterator'
        with self.dataset.batches(fields=['depths'], workdir=workdir) as run:
            self.assertIsNone(run._process)
            sizes = [batch.num_rows for batch in run]
            self.assertEqual(sum(sizes), 17020)
            self.assertLessEqual(max(sizes), 1024)
            self.assertTrue(run.result.manifest_path.is_file())
        with self.dataset.batches(workdir=self.root / 'cancel') as run:
            iterator = iter(run)
            next(iterator)
        iterator.close()
        self.assertIsNone(run.result)
        empty = self.root / 'empty.bed'
        empty.write_text('')
        exported = self.dataset.materialize(self.root / 'empty.arrow', regions=empty)
        with pa.ipc.open_stream(exported.artifact_path) as stream:
            self.assertEqual(stream.read_all().num_rows, 0)
        with self.assertRaises(EvidenceProcessError):
            self.dataset.materialize(self.root / 'absent.arrow', fields=['quality-histograms'])
        self.assertFalse((self.root / 'absent.arrow').exists())
        with self.assertRaises(EvidenceProcessError) as refused:
            self.dataset.materialize(self.root / 'low.arrow', memory_budget_mb=1)
        self.assertIn(refused.exception.returncode, (3, 4))
        self.assertFalse((self.root / 'low.arrow').exists())

    def test_parquet_parts_preserve_uint64_lists_and_are_deterministic(self):
        outputs = [self.dataset.export_parquet(self.root / f'parquet-{n}', memory_budget_mb=budget)
                   for n, budget in enumerate([128, 256, 512])]
        reference = None
        for exported in outputs:
            self.assertTrue(exported.manifest_path.is_file())
            files = sorted(exported.directory.glob('*.parquet'))
            self.assertEqual(len(files), 2)
            for path in files:
                parquet = pq.ParquetFile(path)
                self.assertLessEqual(parquet.metadata.num_rows, 16384)
                self.assertTrue(all(parquet.metadata.row_group(i).num_rows <= 1024
                                    for i in range(parquet.metadata.num_row_groups)))
            table = pa.concat_tables([pq.read_table(path) for path in files])
            self.assertEqual(table.num_rows, 17020)
            self.assertEqual(table.schema.field('callable_depth').type, pa.uint64())
            self.assertEqual(table.schema.field('allele_base_quality_sum').type.value_type, pa.uint64())
            data = [path.read_bytes() for path in files]
            if reference is None:
                reference = data
            self.assertEqual(data, reference)
            receipt = json.loads(exported.manifest_path.read_text())
            self.assertEqual(len(receipt['outputs']), 2)
        with self.assertRaises(EvidenceProcessError):
            self.dataset.export_parquet(outputs[0].directory)
        self.assertFalse(any(p.name.endswith('.partial') for p in self.root.iterdir()))
        plan = subprocess.run([str(self.binary), 'dataset', 'export', '--dataset', str(self.manifest),
            '-o', str(self.root / 'plan-only'), '--plan'], capture_output=True)
        self.assertEqual(plan.returncode, 0, plan.stderr.decode())
        self.assertFalse((self.root / 'plan-only').exists())

    def test_partial_reuse_computes_holes_and_keeps_scientific_identity(self):
        sparse_cache = self.root / 'sparse-cache'
        source = materialize_evidence(self.reference, self.bam, self.root / 'sparse.arrow',
            regions=self.subset, fields=['depths', 'alleles', 'quality-sums', 'allele-quality'],
            cache_dir=sparse_cache, binary=self.binary)
        source_manifest = Path(json.loads(source.manifest_path.read_text())['measurements']['execution.evidence_dataset_manifest'])
        fresh = materialize_evidence(self.reference, self.bam, self.root / 'fresh-all.arrow',
            regions=self.bed, fields=['depths', 'alleles'], binary=self.binary)
        for budget, tile in [(256, 1), (384, 137), (512, 16384)]:
            reused = materialize_evidence(self.reference, self.bam, self.root / f'reuse-{budget}.arrow',
                regions=self.bed, fields=['depths', 'alleles'], reuse_dataset=source_manifest,
                memory_budget_mb=budget, tile_bases=tile, binary=self.binary)
            self.assertEqual(reused.artifact_path.read_bytes(), fresh.artifact_path.read_bytes())
            receipt = json.loads(reused.manifest_path.read_text())
            self.assertEqual(receipt['measurements']['execution.reused_loci'], '55')
            self.assertEqual(receipt['measurements']['execution.computed_loci'], str(17020 - 55))
            fresh_receipt = json.loads(fresh.manifest_path.read_text())
            self.assertEqual(receipt['params']['science.blake3'], fresh_receipt['params']['science.blake3'])
        replay = subprocess.run([str(self.binary), 'reproduce', '--manifest', str(reused.manifest_path),
            '--inputs', str(self.root), '--binary', str(self.binary), '--no-attest'], capture_output=True)
        self.assertEqual(replay.returncode, 0, replay.stderr.decode() + replay.stdout.decode())
        for options in [dict(mapq=21), dict(fields=['quality-histograms'])]:
            with self.assertRaises(EvidenceProcessError):
                materialize_evidence(self.reference, self.bam, self.root / 'incompatible.arrow',
                    regions=self.bed, reuse_dataset=source_manifest, binary=self.binary, **options)
            self.assertFalse((self.root / 'incompatible.arrow').exists())
        all_reused = materialize_evidence(self.reference, self.bam, self.root / 'all-reused.arrow',
            regions=self.subset, fields=['depths'], reuse_dataset=self.manifest, binary=self.binary)
        receipt = json.loads(all_reused.manifest_path.read_text())
        self.assertEqual(receipt['measurements']['execution.reused_loci'], '55')
        self.assertEqual(receipt['measurements']['execution.computed_loci'], '0')

    def test_mutation_refuses_without_overwriting_previous_artifacts(self):
        saved = self.dataset.materialize(self.root / 'saved.arrow')
        before = saved.artifact_path.read_bytes()
        before_receipt = saved.manifest_path.read_bytes()
        part = next(self.manifest.parent.glob('*/evidence.arrow'))
        data = bytearray(part.read_bytes())
        data[len(data) // 2] ^= 1
        part.write_bytes(data)
        with self.assertRaises(EvidenceProcessError):
            self.dataset.materialize(saved.artifact_path, force=True)
        self.assertEqual(saved.artifact_path.read_bytes(), before)
        self.assertEqual(saved.manifest_path.read_bytes(), before_receipt)
        with self.assertRaises(EvidenceProcessError):
            self.dataset.export_parquet(self.root / 'corrupt-export')
        self.assertFalse((self.root / 'corrupt-export').exists())


if __name__ == '__main__':
    unittest.main()
