"""Independent typed-variant checks for evidence annotation and atomic replay."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pysam

from rosalind import materialize_evidence, EvidenceProcessError


class EvidenceAnnotationTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.binary = Path(os.environ['ROSALIND_BINARY']).resolve()
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>chr1\n' + 'A' * 40000 + '\n')
        pysam.faidx(str(self.reference))
        self.bam = self.root / 'reads.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': 40000}],
                  'RG': [{'ID': 'rg1', 'SM': 'alignment-sample'}]}
        with pysam.AlignmentFile(str(self.bam), 'wb', header=header) as out:
            for start in [0, 16999]:
                for i, (base, bq, mq, flags) in enumerate([
                    ('A', 30, 60, 0), ('C', 40, 50, 16), ('T', 25, 30, 0),
                    ('C', 5, 60, 0), ('T', 40, 60, 1024), ('A', 255, 60, 0)]):
                    r = pysam.AlignedSegment(out.header)
                    r.query_name = f'{start}-{i}'
                    r.query_sequence = base * 4
                    r.query_qualities = [bq] * 4
                    r.flag = flags
                    r.reference_id = 0
                    r.reference_start = start
                    r.mapping_quality = mq
                    r.cigarstring = '4M'
                    r.set_tag('RG', 'rg1')
                    out.write(r)
        pysam.index(str(self.bam))
        self.sites = self.root / 'sites.vcf'
        self.sites.write_text(
            '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=40000>\n'
            '##INFO=<ID=OLD,Number=1,Type=String,Description="Keep">\n'
            '##FILTER=<ID=q10,Description="Original filter">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Original depth">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
            'chr1\t35001\tzero\tA\tG\t.\t.\tOLD=zero\tGT:DP\t./.:.\t0/0:0\n'
            'chr1\t2\tmulti\tA\tT,C\t31\tq10\tOLD=first\tGT:DP\t2|1:7\t0/2:3\n'
            'chr1\t17001\tlater\tA\tC\t17\tPASS\tOLD=middle\tGT:DP\t0/1:2\t1|1:8\n'
            'chr1\t2\tduplicate\tA\tC,T\t9\tPASS\tOLD=last\tGT:DP\t1/2:4\t0|1:1\n')

    def run_cli(self, *args):
        result = subprocess.run([str(self.binary), *map(str, args)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        return result

    def records(self, path):
        with pysam.VariantFile(str(path)) as f:
            return [r.copy() for r in f]

    def test_record_order_allele_order_samples_and_exact_sums_survive_formats_and_replay(self):
        original = self.records(self.sites)
        baseline_evidence = None
        for index, (extension, mode, budget, workers, tile) in enumerate([
            ('vcf', 'w', 512, 1, 1), ('vcf.gz', 'wz', 768, 2, 137), ('bcf', 'wb', 1024, 8, 16384)]):
            source = self.root / f'source-{index}.{extension}'
            with pysam.VariantFile(str(self.sites)) as inp:
                with pysam.VariantFile(str(source), mode, header=inp.header) as out:
                    for r in inp:
                        out.write(r)
            annotated = self.root / f'annotated-{index}.{extension}'
            result = materialize_evidence(self.reference, self.bam, self.root / f'evidence-{index}.arrow',
                sites=source, annotated_variants=annotated, fields='all-supported',
                binary=self.binary, memory_budget_mb=budget, workers=workers, tile_bases=tile)
            self.assertEqual(result.annotated_variants_path, annotated)
            data = result.artifact_path.read_bytes()
            if baseline_evidence is None:
                baseline_evidence = data
            self.assertEqual(data, baseline_evidence)
            with pa.ipc.open_stream(result.artifact_path) as stream:
                table = stream.read_all()
            self.assertEqual(table['pos'].to_pylist(), [2, 17001, 35001])
            self.assertEqual(table['allele_base_quality_sum'].to_pylist(), [[30, 40, 0, 25], [30, 40, 0, 25], [0, 0, 0, 0]])
            actual = self.records(annotated)
            self.assertEqual(len(actual), len(original))
            for before, after in zip(original, actual):
                self.assertEqual((after.contig, after.pos, after.id, after.alleles, after.qual),
                                 (before.contig, before.pos, before.id, before.alleles, before.qual))
                self.assertEqual(list(after.filter), list(before.filter))
                self.assertEqual(after.info['OLD'], before.info['OLD'])
                self.assertEqual(list(after.samples), list(before.samples))
                for sample in before.samples:
                    self.assertEqual(dict(after.samples[sample]), dict(before.samples[sample]))
                    self.assertEqual(after.samples[sample].phased, before.samples[sample].phased)
            self.assertEqual(actual[0].info['RSL_DP'], 0)
            self.assertEqual(actual[0].info['RSL_AD'], (0, 0))
            self.assertEqual(actual[1].info['RSL_AD'], (1, 1, 1))
            self.assertEqual(actual[1].info['RSL_ADF'], (1, 1, 0))
            self.assertEqual(actual[1].info['RSL_ADR'], (0, 0, 1))
            self.assertEqual(actual[1].info['RSL_BQS'], ('30', '25', '40'))
            self.assertEqual(actual[3].info['RSL_BQS'], ('30', '40', '25'))
            self.assertEqual(actual[1].info['RSL_MQS'], ('60', '30', '50'))
            self.assertEqual(actual[1].info['RSL_RPS'], ('1', '1', '2'))
            self.assertEqual(actual[1].info['RSL_RLS'], ('4', '4', '4'))
            self.assertEqual(actual[1].info['RSL_PF'], 6)
            self.assertEqual(actual[1].info['RSL_ED'], 5)
            manifest = json.loads(result.manifest_path.read_text())
            self.assertEqual(manifest['params']['annotation.semantics'], 'snv-read-evidence-info-v1')
            self.assertEqual(manifest['params']['evidence.fields_version'], '2')
            self.run_cli('verify', '--manifest', result.manifest_path)
            self.run_cli('reproduce', '--manifest', result.manifest_path, '--inputs', self.root, '--binary', self.binary, '--no-attest')

    def test_annotation_bytes_are_invariant_and_zero_is_distinct_from_absent(self):
        baseline = None
        for i, (budget, workers, tile) in enumerate([(512, 1, 1), (768, 2, 137), (1024, 8, 16384)]):
            output = self.root / f'counts-{i}.vcf.gz'
            result = materialize_evidence(self.reference, self.bam, self.root / f'counts-{i}.arrow',
                sites=self.sites, annotated_variants=output, fields='depths,alleles',
                memory_budget_mb=budget, workers=workers, tile_bases=tile, binary=self.binary)
            if baseline is None:
                baseline = output.read_bytes()
            self.assertEqual(output.read_bytes(), baseline)
            self.assertEqual(self.records(output)[0].info['RSL_DP'], 0)
            self.assertNotIn('RSL_BQS', self.records(output)[0].info)
            self.run_cli('verify', '--manifest', result.manifest_path)

    def test_invalid_annotation_preserves_all_existing_destinations(self):
        evidence = self.root / 'keep.arrow'
        annotated = self.root / 'keep.vcf'
        receipt = Path(str(evidence) + '.manifest.json')
        for path in [evidence, annotated, receipt]:
            path.write_text('original')
        # A header conflict is caught before the extraction or artifact staging.
        text = self.sites.read_text().replace('#CHROM', '##INFO=<ID=RSL_DP,Number=1,Type=Integer,Description="Existing">\n#CHROM')
        self.sites.write_text(text)
        with self.assertRaises(EvidenceProcessError) as error:
            materialize_evidence(self.reference, self.bam, evidence, sites=self.sites,
                annotated_variants=annotated, binary=self.binary, force=True)
        self.assertIn('already exists', str(error.exception))
        for path in [evidence, annotated, receipt]:
            self.assertEqual(path.read_text(), 'original')
        self.assertFalse(list(self.root.glob('*.partial')))
        self.assertFalse(list(self.root.glob('*.tmp*')))

    def test_empty_selection_and_late_resource_failure_have_explicit_artifacts(self):
        self.sites.write_text(''.join(line + '\n' for line in self.sites.read_text().splitlines() if line.startswith('#')))
        output = self.root / 'empty.arrow'
        annotated = self.root / 'empty.bcf'
        result = materialize_evidence(self.reference, self.bam, output, sites=self.sites,
            annotated_variants=annotated, binary=self.binary, fields='depths,alleles')
        self.assertEqual(self.records(annotated), [])
        self.run_cli('verify', '--manifest', result.manifest_path)
        self.run_cli('reproduce', '--manifest', result.manifest_path, '--inputs', self.root,
                     '--binary', self.binary, '--no-attest')

        primary = self.root / 'failed.arrow'
        annotation = self.root / 'failed.vcf.gz'
        command = [str(self.binary), 'analyze', 'evidence', '--reference', str(self.reference),
                   '--alignments', str(self.bam), '--sites', str(self.sites), '--fields', 'depths,alleles',
                   '--memory-budget-mb', '512', '--format', 'arrow-ipc', '-o', str(primary),
                   '--annotated-variants', str(annotation)]
        result = subprocess.run(command, capture_output=True, text=True,
            env=dict(os.environ, ROSALIND_FORCE_FINAL_RSS_BYTES=str(4 << 30)))
        self.assertEqual(result.returncode, 4, result.stderr)
        self.assertFalse(primary.exists())
        self.assertFalse(annotation.exists())
        self.assertTrue(Path(str(primary) + '.partial').is_file())
        self.assertTrue(Path(str(annotation) + '.partial').is_file())
        manifest = Path(str(primary) + '.manifest.json')
        recorded = json.loads(manifest.read_text())
        self.assertEqual(recorded['params']['run_status'], 'resource-failed')
        self.assertTrue(all(out['path'].endswith('.partial') for out in recorded['outputs']))
        verified = subprocess.run([str(self.binary), 'verify', '--manifest', str(manifest)],
                                  capture_output=True, text=True)
        self.assertEqual(verified.returncode, 5, verified.stderr)
        self.assertIn('recorded peak', verified.stderr)
        self.assertNotIn('hash mismatch', verified.stderr.lower())

    def test_oversized_declared_annotation_envelope_cannot_wrap_admission(self):
        result = subprocess.run([str(self.binary), 'analyze', 'evidence',
            '--reference', str(self.reference), '--alignments', str(self.bam),
            '--sites', str(self.sites), '--fields', 'depths,alleles',
            '--annotated-variants', str(self.root / 'too-large.vcf'),
            '-o', str(self.root / 'too-large.arrow'), '--plan',
            '--max-variant-header-bytes', str((1 << 64) - 1)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn('memory envelope is too large', result.stderr)
        self.assertEqual(result.stdout, '')
        self.assertFalse((self.root / 'too-large.arrow').exists())
        self.assertFalse((self.root / 'too-large.vcf').exists())


if __name__ == '__main__':
    unittest.main()
