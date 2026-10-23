"""Pure Python tests: network and pysam are replaced by finite fixtures."""

import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    "sufficiency_prepare", Path(__file__).resolve().parents[1] / "prepare.py")
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def protocol():
    return {"sample": "HG002", "intervals": [
        {"contig": "chr20", "start": 10, "end": 20, "partition": "train", "block_id": "a"},
        {"contig": "chr20", "start": 30, "end": 40, "partition": "test", "block_id": "b"}],
        "downsample": {"seed": 42, "numerator": 10, "denominator": 35},
        "limits": {"download_bytes": 10000, "scratch_bytes": 100000,
                   "analysis_seconds": 60, "network_timeout_seconds": 2}}


class Response(io.BytesIO):
    def __init__(self, body, length=None):
        super().__init__(body)
        self.headers = {"Content-Length": str(len(body) if length is None else length)}


class Header:
    def to_dict(self):
        return {"RG": [{"ID": "one", "SM": "HG002"}]}


class Read:
    def __init__(self, name, start, end, rg="one"):
        self.query_name = name
        self.reference_start = start
        self.reference_end = end
        self.query_length = 5
        self.rg = rg

    def has_tag(self, key):
        return self.rg is not None

    def get_tag(self, key):
        return self.rg


class Source:
    references = ("chr20", "other")
    lengths = (100, 10)
    header = Header()

    def __init__(self, reads=()):
        self.reads = reads
        self.fetch_calls = []

    def fetch(self, *args):
        self.fetch_calls.append(args)
        return iter(self.reads)

    def has_index(self):
        return True

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class PreparationTests(unittest.TestCase):
    def prepared_state(self, root):
        def identity(name):
            path = root / name
            body = name.encode()
            path.write_bytes(body)
            return {"path": str(path), "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()}

        sources = {key: identity("source-" + key) for key in prepare.REQUIRED}
        identities = {"sources": sources, "original_bam": dict(sources["reads_bam"])}
        for key in ("reference", "reference_fai", "downsampled_bam", "downsampled_bai", "regions"):
            identities[key] = identity(key)
        return {"identities": identities}

    def test_final_input_verification_hashes_duplicate_path_only_once(self):
        with tempfile.TemporaryDirectory() as directory:
            state = self.prepared_state(Path(directory))
            original = prepare._digest
            with patch.object(prepare, "_digest", wraps=original) as hashed:
                report = prepare.verify_prepared_inputs(state)
            self.assertEqual(report["files_checked"], 10)
            self.assertEqual(hashed.call_count, 10)
            self.assertEqual(report["bytes_checked"],
                             sum(path.stat().st_size for path in Path(directory).iterdir()))

    def test_final_input_verification_refuses_mutation_missing_and_conflicting_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = self.prepared_state(root)
            path = Path(state["identities"]["sources"]["truth_bed"]["path"])
            original = path.read_bytes()
            path.write_bytes(b"X" * len(original))
            with self.assertRaisesRegex(ValueError, "SHA256 changed"):
                prepare.verify_prepared_inputs(state)
            path.write_bytes(original)
            state["identities"]["original_bam"]["sha256"] = "a" * 64
            with self.assertRaisesRegex(ValueError, "conflicting"):
                prepare.verify_prepared_inputs(state)
            state["identities"]["original_bam"] = dict(state["identities"]["sources"]["reads_bam"])
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                prepare.verify_prepared_inputs(state)

    def test_final_input_verification_refuses_missing_identity_and_symlink(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = self.prepared_state(root)
            identity = state["identities"].pop("reference_fai")
            with self.assertRaisesRegex(ValueError, "missing prepared identity"):
                prepare.verify_prepared_inputs(state)
            state["identities"]["reference_fai"] = identity
            path = Path(identity["path"])
            target = root / "elsewhere"
            path.rename(target)
            path.symlink_to(target)
            with self.assertRaisesRegex(ValueError, "regular"):
                prepare.verify_prepared_inputs(state)

    def test_import_does_not_import_pysam(self):
        with patch.object(prepare.importlib, "import_module") as native:
            spec.loader.exec_module(prepare)
        native.assert_not_called()

    def test_protocol_rejects_overlap_and_unsafe_coordinates(self):
        for mutation in (
                lambda p: p["intervals"][1].update(start=19),
                lambda p: p["intervals"][0].update(start=True),
                lambda p: p["intervals"][0].update(contig="chr21"),
                lambda p: p["intervals"][1].update(block_id="a"),
                lambda p: p["downsample"].update(numerator=35),
                lambda p: p["limits"].update(download_bytes=5 << 30)):
            value = protocol()
            mutation(value)
            with self.assertRaises(ValueError):
                prepare._validate_protocol(value)

    def test_half_open_overlap_keeps_spanning_record_once(self):
        self.assertFalse(prepare._overlaps(0, 10, [10, 30], [20, 40]))
        self.assertFalse(prepare._overlaps(20, 30, [10, 30], [20, 40]))
        self.assertTrue(prepare._overlaps(19, 31, [10, 30], [20, 40]))
        self.assertTrue(prepare._overlaps(0, 50, [10, 30], [20, 40]))

    def test_template_key_uses_source_rg_name_and_seed_exactly(self):
        key = json.dumps([42, "a" * 64, "one", "template"], separators=(",", ":")).encode()
        value = int.from_bytes(hashlib.sha256(key).digest(), "big")
        self.assertEqual(prepare._keep_template(42, "a" * 64, "one", "template", 10, 35),
                         value * 35 < 10 * (1 << 256))
        with self.assertRaises(ValueError):
            prepare._keep_template(42, "a" * 64, "one", None, 10, 35)

    def test_dictionary_order_sample_and_missing_index_refuse(self):
        intervals = protocol()["intervals"]
        source = Source()
        self.assertEqual(prepare._audit_inputs(source, Source(), "HG002", intervals), {"one"})
        wrong = SimpleNamespace(references=("other", "chr20"), lengths=(10, 100))
        with self.assertRaisesRegex(ValueError, "ordered"):
            prepare._audit_inputs(source, wrong, "HG002", intervals)
        with self.assertRaisesRegex(ValueError, "sole named"):
            prepare._audit_inputs(source, Source(), "HG001", intervals)
        source.has_index = lambda: False
        with self.assertRaisesRegex(ValueError, "index"):
            prepare._audit_inputs(source, Source(), "HG002", intervals)

    def test_lock_accepts_extra_pins_but_download_selects_five(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lock.tsv"
            rows = [f"{key}\t{key}.dat\t{'a' * 64}\thttps://example.invalid/{key}"
                    for key in (*prepare.REQUIRED, "stratum")]
            path.write_text("\n".join(rows))
            result, raw = prepare._source_lock(path)
            self.assertEqual([x["id"] for x in result], list(prepare.REQUIRED))
            self.assertEqual(raw, path.read_bytes())
            path.write_text("\n".join(rows + [rows[0]]))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                prepare._source_lock(path)

    def test_hash_rejects_symlink_and_fifo(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "file").write_bytes(b"fixture")
            (root / "link").symlink_to(root / "file")
            with self.assertRaisesRegex(ValueError, "regular"):
                prepare._digest(root / "link")
            prepare.os.mkfifo(root / "fifo")
            with self.assertRaisesRegex(ValueError, "regular"):
                prepare._digest(root / "fifo")

    def test_download_verified_no_overwrite(self):
        body = b"verified source"
        row = {"id": "reads_bam", "filename": "reads.bam", "url": "https://example.invalid/a",
               "sha256": hashlib.sha256(body).hexdigest()}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            budget = prepare._Budget(root, protocol()["limits"])
            result = prepare._download(row, root, budget, lambda *a, **kw: Response(body))
            self.assertEqual(result["bytes"], len(body))
            self.assertEqual(budget.download_bytes, len(body))
            with self.assertRaises(FileExistsError):
                prepare._download(row, root, budget, lambda *a, **kw: Response(b"replacement"))
            self.assertEqual((root / "reads.bam").read_bytes(), body)

    def test_download_retry_bytes_count_against_one_budget(self):
        body = b"12345678"
        row = {"id": "reads_bam", "filename": "reads.bam", "url": "https://example.invalid/a",
               "sha256": hashlib.sha256(body).hexdigest()}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            limits = {**protocol()["limits"], "download_bytes": 12}
            budget = prepare._Budget(root, limits)
            responses = iter([Response(body[:5], 8), Response(body)])
            with self.assertRaisesRegex(ValueError, "download byte budget"):
                prepare._download(row, root, budget, lambda *a, **kw: next(responses))
            self.assertEqual(budget.download_bytes, 5)
            self.assertEqual(list(root.iterdir()), [])

    def test_bad_hash_never_publishes_or_retries(self):
        row = {"id": "reference", "filename": "ref.gz", "url": "https://example.invalid/a",
               "sha256": "a" * 64}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            budget = prepare._Budget(root, protocol()["limits"])
            with self.assertRaisesRegex(ValueError, "SHA256"):
                prepare._download(row, root, budget, lambda *a, **kw: Response(b"wrong"))
            self.assertEqual(budget.attempts, 1)
            self.assertEqual(list(root.iterdir()), [])

    def test_deadline_and_scratch_exhaustion_refuse(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            budget = prepare._Budget(root, protocol()["limits"])
            with self.assertRaisesRegex(ValueError, "scratch"):
                budget.check_scratch(100001)
            budget.deadline = 0
            with self.assertRaises(TimeoutError):
                budget.check_time()

    def test_downsample_preserves_physical_duplicates_and_template_mates(self):
        source_sha = "a" * 64
        retained_name = next(str(i) for i in range(100) if prepare._keep_template(
            42, source_sha, "one", str(i), 10, 35))
        removed_name = next(str(i) for i in range(100) if not prepare._keep_template(
            42, source_sha, "one", str(i), 10, 35))
        duplicated = Read(retained_name, 15, 35)
        mate = Read(retained_name, 35, 39)
        source = Source([Read(retained_name, 0, 5), duplicated, duplicated, mate,
                         Read(removed_name, 11, 16)])
        written = []

        class Writer(Source):
            def write(self, record):
                written.append(record)

        native = SimpleNamespace(AlignmentFile=lambda *a, **kw: Writer())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metrics = prepare._downsample(source, root / "out.bam", native, protocol()["intervals"],
                {"one"}, source_sha, protocol()["downsample"], prepare._Budget(root, protocol()["limits"]))
        self.assertEqual(source.fetch_calls, [("chr20",)])
        self.assertEqual(written, [duplicated, duplicated, mate])
        self.assertEqual(metrics["selected_records"], 4)
        self.assertEqual(metrics["retained_records"], 3)

    def test_prepare_complete_mocked_native_keeps_original_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            resources = root / "resources.tsv"
            bodies = {"reads_bam": b"original bam", "reads_bai": b"source index",
                      "reference": gzip.compress(b">chr20\nACGT\n>other\nACGT\n"),
                      "truth_vcf": b"vcf fixture", "truth_bed": b"bed fixture"}
            resources.write_text("\n".join(
                f"{key}\t{key}.dat\t{hashlib.sha256(body).hexdigest()}\thttps://example.invalid/{key}"
                for key, body in bodies.items()))
            original = Source([Read("test", 11, 16)])

            class Writer(Source):
                def __init__(self, path):
                    self.path = Path(path)
                    self.path.write_bytes(b"new bam header")

                def write(self, record):
                    with self.path.open("ab") as output:
                        output.write(record.query_name.encode())

            native = SimpleNamespace(__version__="mock",
                faidx=lambda path: Path(path + ".fai").write_bytes(b"reference index"),
                index=lambda path: Path(path + ".bai").write_bytes(b"new index"),
                FastaFile=lambda *a: Source(),
                AlignmentFile=lambda path, mode, **kw: original if mode == "rb" else Writer(path))

            actual_download = prepare._download
            def bounded_download(row, path, budget):
                return actual_download(row, path, budget,
                    lambda *a, **kw: Response(bodies[row["id"]]))

            with patch.object(prepare, "_download", bounded_download), \
                    patch.object(prepare.importlib, "import_module", return_value=native):
                result = prepare.prepare_inputs(resources, protocol(), root / "out")
            self.assertEqual(Path(result["original_bam"]).read_bytes(), bodies["reads_bam"])
            self.assertEqual(result["sample"], "HG002")
            self.assertEqual(result["metrics"]["interval_bases"], 20)
            self.assertEqual(result["metrics"]["download_bytes"], sum(map(len, bodies.values())))
            self.assertEqual(Path(result["reference"]).read_bytes(), gzip.decompress(bodies["reference"]))
            self.assertEqual(json.loads((root / "out/preparation.json").read_text()), result)
            self.assertEqual(result["identities"]["dictionary_contigs"], 2)
            with self.assertRaises(FileExistsError):
                prepare.prepare_inputs(resources, protocol(), root / "out")


if __name__ == "__main__":
    unittest.main()
