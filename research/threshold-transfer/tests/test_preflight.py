"""Hosted metadata tests with synthetic records; never fetch public data."""
from contextlib import contextmanager
import copy
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import preflight
import prepare_fresh as fresh

SAVED_HEADER = HERE / "reports" / "attempt-3" / "fresh" / "source-header.sam"


class Header:
    def __init__(self, raw=None):
        self.raw = SAVED_HEADER.read_text() if raw is None else raw
        self.dictionary = {"HD": {}, "SQ": []}
        for line in self.raw.splitlines():
            if line.startswith(("@HD\t", "@SQ\t")):
                values = dict(field.split(":", 1) for field in line.split("\t")[1:])
                if line.startswith("@HD"):
                    self.dictionary["HD"] = values
                else:
                    values["LN"] = int(values["LN"])
                    self.dictionary["SQ"].append(values)

    def __str__(self):
        return self.raw

    def to_dict(self):
        return copy.deepcopy(self.dictionary)


class Record:
    def __init__(self, qname="HWI-D00360:5:H814YADXX:1:1101:1319:2041", *, rg=None, start=55000000,
                 paired=True, read1=True, read2=False, secondary=False, supplementary=False):
        self.query_name, self.rg, self.reference_start = qname, rg, start
        self.is_paired, self.is_read1, self.is_read2 = paired, read1, read2
        self.is_secondary, self.is_supplementary = secondary, supplementary
        self.flag = 1 * paired + 64 * read1 + 128 * read2 + 256 * secondary + 2048 * supplementary

    def has_tag(self, tag):
        assert tag == "RG"
        return self.rg is not None

    def get_tag(self, tag):
        assert tag == "RG"
        return self.rg


def protocol():
    return {"fresh": {"intervals": [{"contig": "chr20", "start": 55000000, "end": 59000000,
                                     "partition": "test", "block_id": "fresh-hg001-55m"}],
                      "conditions": [{"name": "original", "numerator": 35, "denominator": 300}],
                      "sources": copy.deepcopy(fresh.DEFAULT_SOURCES)},
            "limits": {"download_bytes": 4 << 30, "scratch_bytes": 30 << 30,
                       "analysis_seconds": 7200, "network_timeout_seconds": 60}}


class MetadataIdentityTests(unittest.TestCase):
    def test_exact_retained_header_source_and_dictionary_bind_technical_namespace(self):
        result = preflight.validate_header_identity(Header(), preflight.EXPECTED_BAM)
        self.assertEqual(result["source_header_sha256"], preflight.HEADER_SHA256)
        self.assertEqual(result["dictionary_sha256"], preflight.DICTIONARY_SHA256)
        self.assertEqual(result["contigs"], 2580)
        self.assertEqual(result["chr20_length"], 64444167)
        self.assertEqual(result["declared_sort_order"], "unsorted")
        self.assertFalse(result["reference_sequence_verified"])
        self.assertEqual(result["technical_read_group"], preflight.technical_read_group_id(preflight.EXPECTED_BAM))
        self.assertIn("original library/read-group identity unavailable", result["technical_group_scope"])

    def test_different_object_header_dictionary_or_invented_rg_is_rejected(self):
        for key in preflight.EXPECTED_BAM:
            source = dict(preflight.EXPECTED_BAM, **{key: "changed"})
            with self.subTest(source=key), self.assertRaisesRegex(ValueError, "exact reviewed"):
                preflight.validate_header_identity(Header(), source)
        with self.assertRaisesRegex(ValueError, "header differs"):
            preflight.validate_header_identity(Header(SAVED_HEADER.read_text() + "@CO\tchanged\n"), preflight.EXPECTED_BAM)
        for kind in ("dictionary", "read_group"):
            header = Header()
            if kind == "dictionary":
                header.dictionary["SQ"][0]["LN"] += 1
            else:
                header.dictionary["RG"] = [{"ID": "invented", "SM": "HG001"}]
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                preflight.validate_header_identity(header, preflight.EXPECTED_BAM)

    def test_full_qname_is_unchanged_including_namespace_and_signed_coordinates(self):
        names = ("HWI-D00360:5:H814YADXX:1:1101:1319:2041",
                 "D_00360.1:005:FLOW-CELL_1:2:1101:-1319:-2041")
        for name in names:
            self.assertEqual(preflight.validate_cluster_qname(name), name)
            self.assertEqual(preflight.validate_no_rg_record_identity(Record(name)), name)
        self.assertNotEqual(preflight.validate_cluster_qname(names[0]),
                            preflight.validate_cluster_qname(names[0].replace("H814YADXX", "H947YADXX")))
        for name in ("1101:1319:2041", names[0] + "/1", names[0] + " 1:N:0:ACGT",
                     names[0].replace(":1:1101:", ":0:1101:"), "", None, "@" + names[0]):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "seven-field"):
                preflight.validate_cluster_qname(name)

    def test_record_rg_and_invalid_primary_pair_flags_fail_without_overrestricting_supplementary(self):
        for record in (Record(rg="unexplained"), Record(paired=False),
                       Record(read1=False), Record(read2=True)):
            with self.subTest(record=record.__dict__), self.assertRaises(ValueError):
                preflight.validate_no_rg_record_identity(record)
        for record in (Record(paired=False, read1=False, secondary=True),
                       Record(paired=False, read1=False, supplementary=True)):
            self.assertEqual(preflight.validate_no_rg_record_identity(record), record.query_name)

    def test_metadata_budget_caps_bytes_and_charges_enclosing_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = fresh._Budget(Path(directory), protocol()["limits"])
            with patch.object(preflight, "MAX_DOWNLOAD_BYTES", 5):
                budget = preflight._MetadataBudget(Path(directory), parent)
            response = io.BytesIO(b"abcdefgh")
            self.assertEqual(budget.read(response, 100), b"abcde")
            self.assertEqual(parent.download_bytes, 5)
            self.assertEqual(budget.download_bytes, 5)
            self.assertLessEqual(budget.deadline, parent.deadline)
            with self.assertRaisesRegex(ValueError, "budget"):
                budget.read(response, 1)
            self.assertEqual(response.read(), b"fgh")


class ProbeWorkerTests(unittest.TestCase):
    def run_worker(self, directory, records, header=None, indexed=True):
        source = SimpleNamespace(header=header or Header(), has_index=lambda: indexed)
        fetched = []

        def fetch(contig, start, end):
            self.assertEqual((contig, start, end), ("chr20", 55000000, 55001000))
            for record in records:
                fetched.append(record)
                yield record

        source.fetch = fetch

        @contextmanager
        def alignment_file(*args, **kwargs):
            yield source

        request = {"url": "http://127.0.0.1/fixture.bam", "source": preflight.EXPECTED_BAM,
                   "bai": str(Path(directory) / "fixture.bai"), "out": str(directory),
                   "region": {"contig": "chr20", "start": 55000000, "end": 55001000}}
        request_path = Path(directory) / "request.json"
        request_path.write_text(json.dumps(request))
        with patch.dict(sys.modules, {"pysam": SimpleNamespace(AlignmentFile=alignment_file)}):
            preflight._probe_worker(request_path)
        return json.loads(Path(directory, "read-metadata.json").read_text()), fetched

    def test_capped_indexed_probe_keeps_physical_duplicates_and_never_loads_outcomes(self):
        repeated = Record()
        records = [repeated] * 70
        with tempfile.TemporaryDirectory() as directory:
            result, fetched = self.run_worker(directory, records)
            self.assertTrue(result["accepted"])
            self.assertEqual(len(fetched), 64)
            self.assertEqual(result["probe_records"], 64)
            self.assertEqual(result["repeated_primary_end_observations"], 63)
            self.assertEqual(result["primary_templates_missing_mate_in_capped_probe"], 1)
            self.assertFalse(result["physical_records_deduplicated"])
            self.assertFalse(result["outcome_data_loaded"])
            self.assertTrue(all("sequence" not in row and "quality" not in row for row in result["records"]))
            self.assertEqual(Path(directory, "source-header.sam").read_text(), SAVED_HEADER.read_text())

    def test_rejected_identity_preserves_header_and_all_capped_metadata(self):
        for records in ([Record(rg="unexpected")], [Record("short-name")],
                        [Record(start=55000002), Record(start=55000001)],
                        [Record(secondary=True)], []):
            with self.subTest(records=records), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    self.run_worker(directory, records)
                result = json.loads(Path(directory, "read-metadata.json").read_text())
                self.assertFalse(result["accepted"])
                self.assertEqual(len(result["records"]), len(records))
                self.assertTrue(Path(directory, "source-header.sam").exists())

    def test_header_mismatch_and_unavailable_index_stop_before_any_record(self):
        for header, indexed in ((Header(SAVED_HEADER.read_text() + "@CO\tchanged\n"), True), (Header(), False)):
            with self.subTest(indexed=indexed), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    self.run_worker(directory, [Record()], header, indexed)
                result = json.loads(Path(directory, "read-metadata.json").read_text())
                self.assertEqual(result["records"], [])
                self.assertFalse(result["accepted"])
                self.assertEqual(Path(directory, "source-header.sam").read_text(), str(header))


class PreflightOrchestrationTests(unittest.TestCase):
    def test_not_hosted_fails_before_io(self):
        with patch.dict("os.environ", {"GITHUB_ACTIONS": "false"}), \
                patch.object(preflight, "_adapter") as adapter:
            with self.assertRaisesRegex(RuntimeError, "hosted"):
                preflight.run_preflight(protocol(), "/unused")
            adapter.assert_not_called()

    def test_failed_worker_retains_final_range_and_bai_receipts_without_large_downloads(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "preflight"
            downloaded = []
            receipt = {"received_bytes": 4, "received_sha256": hashlib.sha256(b"test").hexdigest()}

            def download(source, path, budget, receipt_path):
                downloaded.append(source)
                path.write_bytes(b"BAI fixture")
                return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "bytes": len(path.read_bytes()), "published_md5_verified": True}

            @contextmanager
            def relay(source, budget):
                ranges = []
                try:
                    yield "http://127.0.0.1/fixture.bam", ranges, []
                finally:
                    ranges.append(receipt)

            with patch.dict("os.environ", {"GITHUB_ACTIONS": "true"}), \
                    patch.object(fresh, "_download", side_effect=download), \
                    patch.object(fresh, "_range_relay", side_effect=relay), \
                    patch.object(preflight.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "probe")):
                with self.assertRaises(subprocess.CalledProcessError):
                    preflight.run_preflight(protocol(), out)
            self.assertEqual(downloaded, [preflight.EXPECTED_BAI])
            result = json.loads((out / "preflight.json").read_text())
            self.assertFalse(result["accepted"])
            self.assertEqual(result["received_bam_bytes"], 4)
            self.assertFalse(result["remote_bam_full_md5_verified"])
            self.assertTrue(Path(result["bai_download"]["path"]).is_file())
            self.assertEqual(json.loads((out / "bam-byte-ranges.json").read_text()), [receipt])


class SavedReceiptTests(unittest.TestCase):
    def write_fixture(self, directory, bai_source):
        root = Path(directory)
        (root / "source-header.sam").write_text(SAVED_HEADER.read_text())
        (root / "source.bam.bai").write_bytes(b"synthetic BAI metadata fixture")
        bai = preflight._file_identity(root / "source.bam.bai")
        identity = preflight.validate_header_identity(Header(), preflight.EXPECTED_BAM)
        region = {"contig": "chr20", "start": 55000000, "end": 55001000}
        worker = dict(identity, accepted=True, outcome_data_loaded=False, probe_region=region,
            read_identity_errors=[], probe_records=1,
            records=[{"qname": Record().query_name, "flag": 65, "has_rg": False, "rg_value": None,
                      "reference_start": 55000000}])
        report = dict(identity, accepted=True, source=preflight.EXPECTED_BAM,
            outcome_data_loaded=False, reference_sequence_verified=False, remote_bam_full_md5_verified=False,
            probe_records=1, wall_seconds=1, received_bam_bytes=4, download_bytes=bai["bytes"] + 4,
            bai_download=dict(bai_source, sha256=bai["sha256"], published_md5_verified=True))
        download = {"source": bai_source, "received_sha256": bai["sha256"], "received_md5": bai["md5"],
                    "received_bytes": bai["bytes"], "published_md5_verified": True, "failure": None}
        for name, value in (("read-metadata.json", worker), ("preflight.json", report),
                            ("download-reads_bai.json", download),
                            ("probe-request.json", {"source": preflight.EXPECTED_BAM, "region": region}),
                            ("bam-byte-ranges.json", [{"received_bytes": 4}]), ("range-errors.json", [])):
            preflight._write(root / name, value)
        self.rehash(root)

    def rehash(self, root):
        files = {p.name: preflight._file_identity(p) for p in root.iterdir() if p.name != "completion.json"}
        preflight._write(root / "completion.json", {"schema": "rosalind-metadata-preflight-v1", "files": files})

    @contextmanager
    def fixture_sources(self):
        content = b"synthetic BAI metadata fixture"
        bai = dict(preflight.EXPECTED_BAI, bytes=len(content), md5=hashlib.md5(content).hexdigest())
        sources = copy.deepcopy(fresh.DEFAULT_SOURCES)
        sources["reads_bai"] = bai
        with patch.object(preflight, "EXPECTED_BAI", bai), patch.object(fresh, "DEFAULT_SOURCES", sources):
            yield bai

    def test_relocated_receipt_revalidates_exact_bytes_and_reusable_bai(self):
        with self.fixture_sources() as bai, tempfile.TemporaryDirectory() as directory:
            self.write_fixture(directory, bai)
            receipt = preflight.validate_preflight_receipt(directory, protocol())
            self.assertTrue(receipt["accepted"])
            self.assertEqual(receipt["bai_path"], str(Path(directory) / "source.bam.bai"))
            self.assertEqual(len(receipt["verified_preflight_files"]), 8)
            Path(directory, "source.bam.bai").write_bytes(b"modified")
            with self.assertRaisesRegex(ValueError, "identity changed"):
                preflight.validate_preflight_receipt(directory, protocol())

    def test_receipt_hashes_do_not_replace_semantic_source_and_budget_admission(self):
        changes = (("preflight.json", "accepted", False),
                   ("preflight.json", "technical_read_group", "invented-library"),
                   ("preflight.json", "wall_seconds", 121),
                   ("preflight.json", "download_bytes", (64 << 20) + 1),
                   ("download-reads_bai.json", "received_md5", "changed"),
                   ("read-metadata.json", "read_identity_errors", ["unexplained RG"]))
        for filename, key, value in changes:
            with self.subTest(file=filename, key=key), self.fixture_sources() as bai, tempfile.TemporaryDirectory() as directory:
                self.write_fixture(directory, bai)
                path = Path(directory) / filename
                contents = json.loads(path.read_text())
                contents[key] = value
                preflight._write(path, contents)
                self.rehash(Path(directory))
                with self.assertRaises(ValueError):
                    preflight.validate_preflight_receipt(directory, protocol())


if __name__ == "__main__":
    unittest.main()
