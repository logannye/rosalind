from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('prepare_representative', Path(__file__).with_name('prepare_representative.py'))
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def lock_rows():
    return [[name, filename, 'a' * 64, 'https://example.invalid/' + filename]
            for name, filename in [('reads_bam', 'reads.bam'), ('reads_bai', 'reads.bai'), ('reference', 'reference.fa.gz')]]


class PreparationTests(unittest.TestCase):
    def test_lock_rejects_malformed_rows_hashes_and_aliased_paths_before_download(self):
        changes = [
            lambda rows: rows[0].pop(),
            lambda rows: rows[0].append('extra'),
            lambda rows: rows[0].__setitem__(2, 'g' * 64),
            lambda rows: rows[0].__setitem__(2, 'a' * 63),
            lambda rows: rows[0].__setitem__(1, '../escape'),
            lambda rows: rows[0].__setitem__(1, '/absolute'),
            lambda rows: rows[0].__setitem__(1, 'folder\\escape'),
            lambda rows: rows[0].__setitem__(1, 'READS.BAI'),
            lambda rows: rows[0].__setitem__(3, 'file:///local'),
            lambda rows: rows[0].__setitem__(3, 'https://@example.invalid/file'),
            lambda rows: rows.append(rows[0].copy()),
            lambda rows: rows.pop(),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = SimpleNamespace(contig='chr20', start=0, end=12000, sites=10,
                                   source_lock=root / 'lock.tsv', download_cache=root / 'cache', output=root / 'out')
            for mutate in changes:
                rows = lock_rows()
                mutate(rows)
                args.source_lock.write_text('\n'.join('\t'.join(row) for row in rows))
                with patch.object(prepare, 'download') as download:
                    with self.assertRaises(ValueError):
                        prepare.prepare(args)
                    download.assert_not_called()
                self.assertFalse(args.download_cache.exists())
                self.assertFalse(args.output.exists())

    def test_valid_lock_can_include_other_unique_locked_resources(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'lock.tsv'
            rows = lock_rows() + [['truth', 'truth.vcf.gz', 'b' * 64, 'https://example.invalid/truth']]
            raw = '# id\tname\tsha256\turl\n' + '\n'.join('\t'.join(row) for row in rows)
            path.write_text(raw)
            selected, captured = prepare.source_lock(path)
            self.assertEqual(selected, rows[:3])
            self.assertEqual(captured, raw.encode())

    def test_invalid_selection_is_rejected_without_reading_sources(self):
        for values in [dict(start=-1), dict(end=0), dict(end=11999), dict(sites=0),
                       dict(sites=12001), dict(contig='bad\tname'), dict(end=1 << 63)]:
            args = SimpleNamespace(contig='chr20', start=0, end=12000, sites=10)
            args.__dict__.update(values)
            with patch.object(prepare, 'source_lock') as lock:
                with self.assertRaises(ValueError):
                    prepare.prepare(args)
                lock.assert_not_called()
        targets = prepare.panel_targets(0, 12000)
        self.assertEqual(sum(end - start for start, end in targets), 10000)
        self.assertEqual(targets[0], (1000, 1500))
        self.assertEqual(targets[-1], (10500, 11000))

    def test_concurrent_downloaders_publish_one_verified_object(self):
        content = b'content-locked public fixture' * 100
        row = ['reads_bam', 'reads.bam', hashlib.sha256(content).hexdigest(), 'https://example.invalid/read']
        barrier = threading.Barrier(2)
        def opened(*args, **kwargs):
            barrier.wait(timeout=5)
            return io.BytesIO(content)
        with tempfile.TemporaryDirectory() as directory, patch.object(prepare.urllib.request, 'urlopen', opened):
            cache = Path(directory)
            with ThreadPoolExecutor(max_workers=2) as pool:
                values = list(pool.map(lambda _: prepare.download(row, cache), range(2)))
            self.assertEqual((cache / row[1]).read_bytes(), content)
            self.assertTrue(all(value['sha256'] == row[2] for value in values))
            self.assertEqual([path.name for path in cache.iterdir()], [row[1]])

    def test_concurrent_conflicting_object_is_never_overwritten(self):
        content = b'expected bytes'
        row = ['reads_bam', 'reads.bam', hashlib.sha256(content).hexdigest(), 'https://example.invalid/read']
        def conflicting_writer(source, destination):
            destination.write_bytes(b'concurrent conflicting object')
            raise FileExistsError(str(destination))
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(prepare.urllib.request, 'urlopen', return_value=io.BytesIO(content)), \
                patch.object(prepare.os, 'link', side_effect=conflicting_writer):
            cache = Path(directory)
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                prepare.download(row, cache)
            self.assertEqual((cache / row[1]).read_bytes(), b'concurrent conflicting object')
            self.assertEqual([path.name for path in cache.iterdir()], [row[1]])

    def test_bad_download_and_existing_symlink_do_not_publish(self):
        row = ['reference', 'reference.fa.gz', 'a' * 64, 'https://example.invalid/reference']
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            with patch.object(prepare.urllib.request, 'urlopen', return_value=io.BytesIO(b'wrong')):
                with self.assertRaisesRegex(ValueError, 'download SHA256 mismatch'):
                    prepare.download(row, cache)
            self.assertEqual(list(cache.iterdir()), [])
            (cache / 'target').write_bytes(b'data')
            (cache / row[1]).symlink_to(cache / 'target')
            with patch.object(prepare.urllib.request, 'urlopen') as download:
                with self.assertRaisesRegex(ValueError, 'not a symlink'):
                    prepare.download(row, cache)
                download.assert_not_called()

    def test_tiny_preparation_keeps_executed_source_and_lock_snapshots(self):
        try:
            import pysam
        except ImportError:
            self.skipTest('optional preparation integration requires pysam0.23.3')
        if pysam.__version__ != '0.23.3':
            self.skipTest('preparation pins pysam0.23.3')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / 'cache'
            cache.mkdir()
            reference = cache / 'reference.fa.gz'
            with gzip.open(reference, 'wb') as out:
                out.write(b'>chr20\n' + b'A' * 12000 + b'\n')
            bam = cache / 'reads.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': 'chr20', 'LN': 12000}], 'RG': [{'ID': 'rg1', 'SM': 'HG002'}]}
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                record = pysam.AlignedSegment(out.header)
                record.query_name = 'synthetic-preparation-test'
                record.query_sequence = 'A' * 35
                record.query_qualities = [30] * 35
                record.reference_id = 0
                record.reference_start = 1100
                record.mapping_quality = 60
                record.cigarstring = '35M'
                record.set_tag('RG', 'rg1')
                out.write(record)
            pysam.index(str(bam))
            rows = [[identifier, path.name, prepare.digest(path), 'https://example.invalid/' + path.name]
                    for identifier, path in [('reads_bam', bam), ('reads_bai', Path(str(bam) + '.bai')), ('reference', reference)]]
            lock = root / 'source-lock.tsv'
            lock.write_text('\n'.join('\t'.join(row) for row in rows))
            initial_lock = lock.read_bytes()
            script = root / 'source.py'
            script.write_text('# initial source snapshot\n')
            initial_script = script.read_bytes()
            original_download = prepare.download
            changed = threading.Event()
            mutation_lock = threading.Lock()
            def downloaded(row, cache):
                with mutation_lock:
                    if not changed.is_set():
                        script.write_text('# source changed while preparation runs\n')
                        lock.write_text('# source lock changed while preparation runs\n')
                        changed.set()
                return original_download(row, cache)
            args = SimpleNamespace(contig='chr20', start=0, end=12000, sites=10,
                                   source_lock=lock, download_cache=cache, output=root / 'prepared')
            with patch.object(prepare, '__file__', str(script)), patch.object(prepare, 'download', downloaded):
                prepare.prepare(args)
            provenance = json.loads((args.output / 'preparation.json').read_text())
            self.assertEqual((args.output / provenance['preparation_script_snapshot']).read_bytes(), initial_script)
            self.assertEqual((args.output / provenance['source_lock_snapshot']).read_bytes(), initial_lock)
            self.assertEqual(provenance['preparation_script_sha256'], hashlib.sha256(initial_script).hexdigest())
            self.assertEqual(provenance['source_lock_sha256'], hashlib.sha256(initial_lock).hexdigest())
            self.assertEqual(provenance['window_records'], 1)


if __name__ == '__main__':
    unittest.main()
