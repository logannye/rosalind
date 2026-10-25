"""Bundle completeness, frozen requests, and failed-publication regressions."""
import csv
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import case_common as common
import checked_parquet
import controls
import verify_bundle as replay


class BundleFixture:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir()
        (self.root / 'native/part').mkdir(parents=True)
        (self.root / 'native/part/evidence.arrow').write_bytes(b'authored-fake-payload-for-control-test')
        (self.root / 'native/part/receipt.json').write_text('{}\n')
        (self.root / 'native/evidence-dataset.manifest.json').write_text('{}\n')
        self.contigs = [{'name': 'chr20', 'length': 63025520}]
        descriptor = {'sources': [{'role': 'alignments', 'path': str(self.root.parent / 'absent.bam')}],
                      'partitions': [{'contig': 0, 'start': 9994240,
                                      'intervals': [{'start': 10000000, 'end': 10000001}],
                                      'arrow': {'path': 'part/evidence.arrow'},
                                      'receipt': {'path': 'part/receipt.json'}}]}
        common.write_json(self.root / 'native/dataset.descriptor.json', descriptor)
        (self.root / 'build-identity.json').write_text('{}\n')
        shutil.copy2(HERE / 'protocol-v1.json', self.root / 'protocol-v1.json')
        source = self.root.parent / 'probes.vcf'
        rows = [dict(contig='chr20', pos=10000001+index*1000, ref='A', alt='C') for index in range(1000)]
        common.write_sites(source, rows, self.contigs)
        requests = common.freeze_queries(source, self.root / 'requests', self.contigs)
        header = ['#contig', 'pos', 'ref', 'requested_alts', *checked_parquet.SCALARS,
                  *checked_parquet.HISTOGRAMS]
        with (self.root / 'oracle.tsv').open('w') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            for index, row in enumerate(rows):
                depth = index % 21
                values = {name: 0 for name in checked_parquet.SCALARS}
                values.update(prefilter_depth=depth, aligned_depth=depth, callable_depth=depth,
                              a=depth, a_fwd=depth, base_quality_sum=depth*30, mapping_quality_sum=depth*60,
                              read_length_sum=depth)
                writer.writerow([row['contig'], row['pos'], 'A', 'C', *values.values(),
                                 f'30:{depth}' if depth else '.', f'60:{depth}' if depth else '.'])
        (self.root / 'native-results').mkdir()
        for request in requests:
            common.write_json(self.root / request['expected'], common.oracle_answers(self.root/'oracle.tsv', request))
            (self.root / f"native-results/{request['id']}.tsv").write_text('retained native query\n')
            common.write_json(self.root / f"native-results/{request['id']}.tsv.manifest.json", {'inputs': []})
        (self.root / 'parquet').mkdir()
        (self.root / 'parquet/inventory.json').write_text('{}\n')
        (self.root / 'synthetic-zero/native').mkdir(parents=True)
        (self.root / 'synthetic-zero/native/evidence-dataset.manifest.json').write_text('{}\n')
        common.write_json(self.root / 'synthetic-zero/native/dataset.descriptor.json', {'partitions': []})
        (self.root / 'synthetic-zero/zero.bed').write_text('synthetic\t39\t40\n')
        (self.root / 'synthetic-zero/absent.bed').write_text('synthetic\t40\t41\n')
        shutil.copy2(HERE.parents[1] / 'examples/cohort-reanalysis/expected.json',
                     self.root / 'synthetic-zero/authored-expected.json')
        self.state = {'protocol_sha256': common.sha256(self.root/'protocol-v1.json'), 'requests': requests,
            'native': 'native/evidence-dataset.manifest.json', 'parquet': 'parquet',
            'parquet_inventory_sha256': common.sha256(self.root/'parquet/inventory.json'),
            'expected_native': {q['id']: 'native-results/'+q['id']+'.tsv' for q in requests},
            'oracle': 'oracle.tsv', 'probe_sites': 'requests/probes.vcf', 'contigs': self.contigs,
            'build_identity': 'build-identity.json', 'source_revision': 'a'*40}
        common.write_json(self.root / 'state.json', self.state)
        common.write_json(self.root / 'report.json', {'passed': True, 'status': 'completed'})
        (self.root / 'REPORT.md').write_text('Synthetic test report.\n')
        (self.root / 'trials.csv').write_text('test-fixture\n')
        self.seal()

    def seal(self):
        common.write_json(self.root/'inventory.json', {'files': common.inventory(self.root, exclude=('inventory.json',)),
            'source_revision': self.state['source_revision'], 'protocol_sha256': self.state['protocol_sha256'],
            'build_identity_sha256': common.sha256(self.root/'build-identity.json')})


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = BundleFixture(self.root / 'bundle')
        self.bundle = self.fixture.root

    def tearDown(self):
        self.temporary.cleanup()

    def test_complete_inventory_and_regenerated_queries(self):
        state, _, _, digest = replay.inspect_bundle(self.bundle, common.sha256(self.bundle/'inventory.json'))
        requests, expected = replay.reconstruct_requests(self.bundle, state, self.root/'regenerated')
        self.assertEqual([len(q['rows']) for q in requests], [250, 500, 1000])
        self.assertEqual(expected['Q3'][0]['callable_depth'], 0)
        self.assertFalse(expected['Q3'][0]['depth_ge_5'])
        self.assertEqual(digest, common.sha256(self.bundle/'inventory.json'))

    def test_removed_partition_rejected_despite_producer_pass(self):
        (self.bundle/'native/part/evidence.arrow').unlink()
        with self.assertRaisesRegex(ValueError, 'inventory mismatch'):
            replay.inspect_bundle(self.bundle)
        self.assertTrue(common.read_json(self.bundle/'report.json')['passed'])

    def test_missing_required_receipt_rejected_with_rewritten_inventory(self):
        (self.bundle/'native-results/Q1.tsv.manifest.json').unlink()
        self.fixture.seal()
        with self.assertRaisesRegex(ValueError, 'required bundle artifact'):
            replay.inspect_bundle(self.bundle)

    def test_missing_partition_rejected_even_with_rewritten_inventory(self):
        (self.bundle/'native/part/evidence.arrow').unlink()
        self.fixture.seal()
        with self.assertRaisesRegex(ValueError, 'descriptor-required partition'):
            replay.inspect_bundle(self.bundle)

    def test_changed_request_rejected_even_after_self_consistent_reseal(self):
        state = common.read_json(self.bundle/'state.json')
        state['requests'][0]['rows'][0]['alt'] = 'T'
        common.write_json(self.bundle/'state.json', state)
        self.fixture.seal()
        replay.inspect_bundle(self.bundle)
        with self.assertRaisesRegex(ValueError, 'reconstructed frozen definitions'):
            replay.reconstruct_requests(self.bundle, state, self.root/'regenerated')

    def test_changed_expected_answer_rejected_by_independent_derivation(self):
        expected = common.read_json(self.bundle/'expected/Q1.json')
        expected[0]['alt_count'] = 123
        common.write_json(self.bundle/'expected/Q1.json', expected)
        self.fixture.seal()
        with self.assertRaisesRegex(ValueError, 'independent oracle derivation'):
            replay.reconstruct_requests(self.bundle, self.fixture.state, self.root/'regenerated')

    def test_external_anchor_rejects_resealed_payload(self):
        digest = common.sha256(self.bundle/'inventory.json')
        (self.bundle/'native/part/evidence.arrow').write_bytes(b'changed')
        self.fixture.seal()
        with self.assertRaisesRegex(ValueError, 'external digest'):
            replay.inspect_bundle(self.bundle, digest)

    def test_unlisted_file_and_path_escape_are_rejected(self):
        (self.bundle/'unlisted').write_text('new')
        with self.assertRaisesRegex(ValueError, 'complete inventory'):
            replay.inspect_bundle(self.bundle)
        (self.bundle/'unlisted').unlink()
        self.fixture.state['native'] = '../outside.json'
        common.write_json(self.bundle/'state.json', self.fixture.state)
        self.fixture.seal()
        with self.assertRaises(ValueError):
            replay.inspect_bundle(self.bundle)

    def test_existing_output_preserved_and_failure_report_retained(self):
        output = self.root/'existing'
        output.mkdir()
        (output/'user-data').write_text('preserve')
        with self.assertRaises(FileExistsError):
            replay.verify_bundle(self.bundle, self.root/'binary', self.bundle/'build-identity.json', output)
        self.assertEqual((output/'user-data').read_text(), 'preserve')
        wrong_identity = self.root/'wrong-build.json'
        wrong_identity.write_text('{"wrong":true}\n')
        with self.assertRaisesRegex(ValueError, 'external build identity'):
            replay.verify_bundle(self.bundle, self.root/'binary', wrong_identity, self.root/'failed')
        self.assertEqual(common.read_json(self.root/'failed/report.json')['status'], 'failed')

    def test_raw_source_present_invalidates_saved_only_claim(self):
        (self.root/'absent.bam').write_text('original input still accessible')
        with self.assertRaisesRegex(ValueError, 'still access an original'):
            replay.source_absence(self.bundle, self.fixture.state)

    def fake_command(self, argv, destination, expect_failure=False, publish=False):
        self.assertTrue(expect_failure)
        destination = Path(destination)
        destination.mkdir(parents=True)
        text = 'evidence payload hash mismatch'
        if '--regions' in argv:
            text = 'query includes 1 loci absent from the dataset'
        if '--sites' in argv and str(argv[argv.index('--sites')+1]).endswith('wrong-ref.vcf'):
            text = 'query SNV REF differs from stored reference base'
        (destination/'stderr.log').write_text(text)
        (destination/'stdout.log').write_text('')
        if publish and '--output' in argv:
            Path(argv[argv.index('--output')+1]).write_text('improper success publication')
        return {'exit_code': 2, 'wall_seconds': 0.01, 'argv': list(map(str, argv))}

    def test_consumed_corruption_preserves_original_manifests_and_reports_refusal(self):
        before = common.inventory(self.bundle)
        def parquet_refusal(directory, rows, **kwargs):
            raise ValueError('requested REF mismatch' if rows[0]['ref'] != 'A' else 'missing requested locus')
        with patch.object(controls, 'command', side_effect=self.fake_command), \
                patch.object(checked_parquet, 'query', side_effect=parquet_refusal), \
                patch.object(controls, 'synthetic_zero', return_value={'status': 'completed'}):
            result = controls.run_controls('explicit-binary', self.bundle, self.root/'controls')
        self.assertEqual(result['status'], 'expected_negative_passed')
        tamper = result['checks']['tamper']
        self.assertEqual(tamper['consumed_locus']['pos'], 10000001)
        self.assertNotEqual(tamper['before_sha256'], tamper['after_sha256'])
        self.assertTrue(tamper['manifests_unchanged'])
        self.assertEqual(before, common.inventory(self.bundle))

    def test_nonzero_exit_with_success_output_is_rejected(self):
        def bad_command(argv, destination, expect_failure=False):
            return self.fake_command(argv, destination, expect_failure, publish=True)
        with patch.object(controls, 'command', side_effect=bad_command):
            with self.assertRaisesRegex(ValueError, 'published a final primary'):
                controls.run_controls('explicit-binary', self.bundle, self.root/'controls')
        self.assertEqual(common.read_json(self.root/'controls/report.json')['status'], 'failed')


if __name__ == '__main__':
    unittest.main()
