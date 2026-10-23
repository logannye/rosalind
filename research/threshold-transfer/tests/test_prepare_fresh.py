"""Hosted adapter checks: no real source downloads or model fits in fixtures."""
import copy
from contextlib import contextmanager
import hashlib
import http.client
import importlib.util
import io
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import urllib.request

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import prepare_fresh as fresh
from features import feature_vectors, FILTER_NAMES, SUM_NAMES


def protocol():
    return {"fresh": {"intervals": [{"contig": "chr20", "start": 55000000,
                                    "end": 59000000, "partition": "test", "block_id": "fresh"}],
                      "conditions": [{"name": "original", "numerator": 35, "denominator": 300},
                                     {"name": "downsampled", "numerator": 10, "denominator": 300}],
                      "seed": 42, "sources": copy.deepcopy(fresh.DEFAULT_SOURCES)},
            "limits": {"download_bytes": 4 << 30, "scratch_bytes": 30 << 30,
                       "analysis_seconds": 60, "network_timeout_seconds": 5}}


def raw_row():
    row = {"ref": "A", "observed": True, "prefilter_depth": 10,
           "aligned_depth": 9, "callable_depth": 8,
           "a": 5, "c": 3, "g": 0, "t": 0,
           "a_fwd": 2, "a_rev": 3, "c_fwd": 1, "c_rev": 2,
           "g_fwd": 0, "g_rev": 0, "t_fwd": 0, "t_rev": 0}
    row.update({"filtered_" + name: 0 for name in FILTER_NAMES})
    row["filtered_secondary"] = 1
    row["filtered_low_base_quality"] = 1
    for name, mean in zip(SUM_NAMES, (30, 60, 50, 100)):
        row[name + "_sum"] = 8 * mean
        row["allele_" + name + "_sum"] = [5 * mean, 3 * mean, 0, 0]
    return row


class SourceBoundaryTests(unittest.TestCase):
    def test_only_frozen_test_intervals_and_source_metadata_are_accepted(self):
        valid = protocol()
        intervals, conditions = fresh._validate(valid)
        self.assertEqual(intervals, valid["fresh"]["intervals"])
        self.assertEqual(conditions, valid["fresh"]["conditions"])
        changes = (("partition", "train"), ("contig", "chr1"), ("end", 70000000))
        for key, value in changes:
            altered = copy.deepcopy(valid)
            altered["fresh"]["intervals"][0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                fresh._validate(altered)
        altered = copy.deepcopy(valid)
        altered["fresh"]["sources"]["reads_bam"]["last_modified"] = "changed"
        with self.assertRaises(ValueError):
            fresh._validate(altered)

    def test_duplicate_conditions_and_wrong_source_depth_are_rejected(self):
        for change in ("duplicate", "denominator", "overlap"):
            altered = protocol()
            if change == "duplicate":
                altered["fresh"]["conditions"][1]["name"] = "original"
            elif change == "denominator":
                altered["fresh"]["conditions"][0]["denominator"] = 35
            else:
                altered["fresh"]["intervals"].append(dict(altered["fresh"]["intervals"][0]))
            with self.subTest(change=change), self.assertRaises(ValueError):
                fresh._validate(altered)

    def test_bam_object_and_exact_range_metadata_must_match(self):
        source = {"bytes": 100, "last_modified": "fixed"}
        fresh._check_bam_headers({"Content-Length": "100", "Last-Modified": "fixed"}, source)
        fresh._check_bam_headers({"Content-Range": "bytes 10-19/100", "Last-Modified": "fixed"}, source, 10, 19)
        for headers in ({"Content-Range": "bytes 10-20/100", "Last-Modified": "fixed"},
                        {"Content-Range": "bytes 10-19/101", "Last-Modified": "fixed"},
                        {"Content-Range": "bytes 10-19/100", "Last-Modified": "changed"}):
            with self.subTest(headers=headers), self.assertRaises(ValueError):
                fresh._check_bam_headers(headers, source, 10, 19)

    def test_preconsumption_input_snapshot_detects_later_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"original")
            snapshot = {"input": fresh._identity(path)}
            self.assertTrue(fresh.verify_input_snapshot(snapshot)["passed"])
            path.write_bytes(b"modified")
            with self.assertRaisesRegex(ValueError, "changed"):
                fresh.verify_input_snapshot(snapshot)

    def test_counts_features_match_rich_f0_without_any_rich_fields(self):
        rich = raw_row()
        counts = {key: rich[key] for key in ("ref", "observed", "callable_depth", "a", "c", "g", "t")}
        for alt in "CGT":
            with self.subTest(alt=alt):
                self.assertEqual(fresh.counts_features(counts, alt), feature_vectors(rich, alt)["f0"])
        counts.pop("c")
        with self.assertRaisesRegex(ValueError, "absent"):
            fresh.counts_features(counts, "C")

    def test_counts_features_preserve_undefined_zero_depth_ratios(self):
        row = {"ref": "A", "callable_depth": 0, "a": 0, "c": 0, "g": 0, "t": 0}
        result = fresh.counts_features(row, "C")
        self.assertIsNone(result["alt_fraction"])
        self.assertEqual(result["alt_fraction_undefined"], 1.0)
        row["a"] = 1
        with self.assertRaisesRegex(ValueError, "reconcile"):
            fresh.counts_features(row, "C")


class FakeHeader:
    def __init__(self, groups):
        self.value = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr20", "LN": 100}], "RG": groups,
                      "PG": [{"ID": "novoalign", "VN": "3.02.07"}]}

    def to_dict(self):
        return copy.deepcopy(self.value)

    def __str__(self):
        return json.dumps(self.value, sort_keys=True) + "\n"


class SpecimenAdmissionTests(unittest.TestCase):
    def test_exact_library_and_canonical_labels_can_coexist_without_changing_source(self):
        labels = sorted(fresh.HG001_LIBRARY_ALIASES | {"HG001", "NA12878"})
        self.assertEqual(len(fresh.HG001_LIBRARY_ALIASES), 14)
        header = FakeHeader([{"ID": f"original-rg-{i}", "SM": label, "LB": label, "PL": "ILLUMINA"}
                             for i, label in enumerate(labels)])
        original = header.to_dict()
        with tempfile.TemporaryDirectory() as directory:
            normalized, admission = fresh._admit_source_header(header, fresh.DEFAULT_SOURCES["reads_bam"], directory)
            self.assertTrue(admission["accepted"])
            self.assertEqual(Path(directory, "source-header.sam").read_text(), str(header))
            self.assertEqual(admission["observed_sample_labels"], labels)
            self.assertEqual([r["original_read_group"] for r in admission["read_groups"]], original["RG"])
            self.assertEqual(json.loads(Path(directory, "source-sample-admission.json").read_text()), admission)
        self.assertEqual(header.to_dict(), original)
        expected = copy.deepcopy(original)
        for group in expected["RG"]:
            group["SM"] = "HG001"
        self.assertEqual(normalized, expected)

    def test_unknown_prefix_case_empty_and_mixed_donor_labels_fail_with_audit(self):
        for label in ("Sample_U0a", "U0A", "U6a", "HG002", "NA24385", "HG001+HG002", "", None):
            header = FakeHeader([{"ID": "known", "SM": "HG001"}, {"ID": "other", "SM": label}])
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, "unmapped SM"):
                    fresh._admit_source_header(header, fresh.DEFAULT_SOURCES["reads_bam"], directory)
                audit = json.loads(Path(directory, "source-sample-admission.json").read_text())
                self.assertFalse(audit["accepted"])
                self.assertIsNone(audit["read_groups"][1]["canonical_sample"])
                self.assertEqual(audit["read_groups"][1]["original_read_group"]["SM"], label)
                self.assertEqual(Path(directory, "source-header.sam").read_text(), str(header))

    def test_missing_and_duplicate_read_group_ids_fail_with_retained_originals(self):
        groups = ([], [{"SM": "HG001"}], [{"ID": "", "SM": "HG001"}],
                  [{"ID": " ", "SM": "HG001"}],
                  [{"ID": "same", "SM": "U0a"}, {"ID": "same", "SM": "NA12878"}])
        for values in groups:
            with self.subTest(groups=values), tempfile.TemporaryDirectory() as directory:
                header = FakeHeader(values)
                with self.assertRaises(ValueError):
                    fresh._admit_source_header(header, fresh.DEFAULT_SOURCES["reads_bam"], directory)
                audit = json.loads(Path(directory, "source-sample-admission.json").read_text())
                self.assertFalse(audit["accepted"])
                self.assertEqual([r["original_read_group"] for r in audit["read_groups"]], values)

    def test_alias_rule_cannot_be_reused_for_another_object(self):
        for key in ("url", "bytes", "md5", "last_modified"):
            source = dict(fresh.DEFAULT_SOURCES["reads_bam"])
            source[key] = "changed"
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, "exact reviewed"):
                    fresh._admit_source_header(FakeHeader([{"ID": "rg", "SM": "U0a"}]), source, directory)
                audit = json.loads(Path(directory, "source-sample-admission.json").read_text())
                self.assertFalse(audit["source_matches_reviewed_object"])

    def test_worker_changes_only_derived_sm_and_preserves_template_hash_inputs(self):
        import pysam
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            reference = root / "reference.fa"
            reference.write_text(">chr20\nACGTACGTACGT\n")
            pysam.faidx(str(reference))
            source_path = root / "fixture.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr20", "LN": 12}],
                      "RG": [{"ID": "rg1", "SM": "U0a", "LB": "library1"},
                             {"ID": "rg2", "SM": "NA12878", "LB": "library2"}]}
            with pysam.AlignmentFile(str(source_path), "wb", header=header) as writer:
                for i, group in enumerate(header["RG"]):
                    record = pysam.AlignedSegment()
                    record.query_name, record.query_sequence = f"original-qname-{i}", "ACGT"
                    record.flag, record.reference_id, record.reference_start = 0, 0, i
                    record.mapping_quality, record.cigar = 60, [(0, 4)]
                    record.query_qualities = pysam.qualitystring_to_array("IIII")
                    record.set_tag("RG", group["ID"])
                    writer.write(record)
            pysam.index(str(source_path))
            source_sha = fresh.pilot.digest(source_path)
            request = {"url": str(source_path), "source": fresh.DEFAULT_SOURCES["reads_bam"],
                       "bai": str(source_path) + ".bai", "reference": str(reference),
                       "out": str(inputs), "audit_out": str(root), "seed": 42,
                       "intervals": [{"contig": "chr20", "start": 0, "end": 12}],
                       "conditions": [{"name": "original", "numerator": 35, "denominator": 300}]}
            request_path = root / "request.json"
            request_path.write_text(json.dumps(request))
            with patch.object(fresh.prepare, "_keep_template", return_value=True) as keep:
                fresh._extract_bams_worker(request_path)
            self.assertEqual(fresh.pilot.digest(source_path), source_sha)
            identity = hashlib.sha256(fresh.pilot.canonical(request["source"]).encode()).hexdigest()
            self.assertEqual([call.args for call in keep.call_args_list], [
                (42, identity, "rg1", "original-qname-0", 35, 300),
                (42, identity, "rg2", "original-qname-1", 35, 300)])
            with pysam.AlignmentFile(str(inputs / "original.bam"), "rb") as derived:
                expected = copy.deepcopy(header)
                for group in expected["RG"]:
                    group["SM"] = "HG001"
                self.assertEqual(derived.header.to_dict(), expected)
                self.assertEqual([(r.query_name, r.get_tag("RG")) for r in derived.fetch("chr20")],
                                 [("original-qname-0", "rg1"), ("original-qname-1", "rg2")])
            result = json.loads((inputs / "bam-preparation.json").read_text())
            self.assertEqual(result["sample"], "HG001")
            self.assertTrue(result["full_ordered_bam_reference_dictionary_match"])


class FakeResponse(io.BytesIO):
    def __init__(self, status, headers, content=b""):
        super().__init__(content)
        self.status, self.headers = status, headers


class RangeRelayTests(unittest.TestCase):
    def test_worker_failure_retains_range_receipts_and_combined_download_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            budget = fresh._Budget(root, protocol()["limits"])
            budget.download_bytes = 123
            receipt = {"start": 0, "requested_end_inclusive": 7, "received_bytes": 8,
                       "received_sha256": hashlib.sha256(b"abcdefgh").hexdigest(), "complete": True}

            @contextmanager
            def relay(source, supplied_budget):
                self.assertIs(supplied_budget, budget)
                logs = []
                try:
                    yield "http://127.0.0.1/source.bam", logs, []
                finally:
                    logs.append(receipt)

            request = {"source": fresh.DEFAULT_SOURCES["reads_bam"], "out": str(inputs)}
            with patch.object(fresh, "_range_relay", side_effect=relay), \
                    patch.object(fresh.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "worker")):
                with self.assertRaises(subprocess.CalledProcessError):
                    fresh._run_bam_extraction(request, budget, root)
            self.assertEqual(json.loads((root / "bam-byte-ranges.json").read_text()), [receipt])
            audit = json.loads((root / "bam-acquisition.json").read_text())
            self.assertFalse(audit["completed_without_error"])
            self.assertFalse(audit["remote_bam_full_md5_verified"])
            self.assertEqual(audit["combined_download_bytes"], 123)
            self.assertEqual(audit["received_bam_bytes"], 8)
            self.assertEqual(audit["failure"]["type"], "CalledProcessError")
            self.assertEqual(json.loads((inputs / "extraction-request.json").read_text())["audit_out"], str(root))

    def test_relay_admission_failure_is_preserved_without_claiming_acquisition(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            budget = fresh._Budget(root, protocol()["limits"])
            with patch.object(fresh, "_range_relay", side_effect=ValueError("remote BAM size changed")):
                with self.assertRaisesRegex(ValueError, "size changed"):
                    fresh._run_bam_extraction({"source": fresh.DEFAULT_SOURCES["reads_bam"]}, budget, root)
            audit = json.loads((root / "bam-acquisition.json").read_text())
            self.assertFalse(audit["completed_without_error"])
            self.assertEqual(audit["received_bam_bytes"], 0)
            self.assertEqual(json.loads((root / "bam-byte-ranges.json").read_text()), [])

    def test_failed_full_resource_download_retains_received_hash_and_rejects_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            budget = fresh._Budget(root, protocol()["limits"])
            source = {"url": "https://example.invalid/resource", "bytes": 8, "md5": "wrong"}
            with patch.object(fresh, "_request", return_value=FakeResponse(200, {"Content-Length": "8"}, b"abcdefgh")):
                with self.assertRaisesRegex(ValueError, "MD5"):
                    fresh._download(source, root / "resource", budget, root / "receipt.json")
            receipt = json.loads((root / "receipt.json").read_text())
            self.assertEqual(receipt["received_bytes"], 8)
            self.assertEqual(receipt["received_sha256"], hashlib.sha256(b"abcdefgh").hexdigest())
            self.assertFalse(receipt["published_md5_verified"])
            self.assertEqual(receipt["failure"]["type"], "ValueError")

    def relay_case(self, directory, *, response_status=206, byte_limit=100):
        source = {"url": "https://example.invalid/source.bam", "bytes": 8, "last_modified": "fixed"}
        limits = protocol()["limits"] | {"download_bytes": byte_limit}
        budget = fresh._Budget(Path(directory), limits)
        calls = []

        def request(url, **kwargs):
            calls.append((url, kwargs))
            if kwargs.get("method") == "HEAD":
                return FakeResponse(200, {"Content-Length": "8", "Accept-Ranges": "bytes", "Last-Modified": "fixed"})
            self.assertEqual(kwargs["headers"]["Range"], "bytes=0-7")
            return FakeResponse(response_status,
                                {"Content-Range": "bytes 0-7/8", "Last-Modified": "fixed"}, b"abcdefgh")

        return source, budget, request, calls

    def test_received_range_bytes_are_accounted_and_hashed(self):
        with tempfile.TemporaryDirectory() as directory:
            source, budget, request, calls = self.relay_case(directory)
            with patch.object(fresh, "_request", side_effect=request):
                with fresh._range_relay(source, budget) as (url, logs, errors):
                    with urllib.request.urlopen(url, timeout=3) as response:
                        self.assertEqual(response.read(), b"abcdefgh")
            self.assertFalse(errors)
            self.assertEqual(budget.download_bytes, 8)
            self.assertEqual(logs[0]["received_sha256"], hashlib.sha256(b"abcdefgh").hexdigest())
            self.assertEqual(len(calls), 2)

    def test_http_200_upstream_cannot_fall_back_to_whole_genome(self):
        with tempfile.TemporaryDirectory() as directory:
            source, budget, request, _ = self.relay_case(directory, response_status=200)
            with patch.object(fresh, "_request", side_effect=request):
                with fresh._range_relay(source, budget) as (url, logs, errors):
                    with urllib.request.urlopen(url, timeout=3) as response:
                        with self.assertRaises(http.client.IncompleteRead):
                            response.read()
            self.assertEqual(budget.download_bytes, 0)
            self.assertTrue(any("whole genome" in error["reason"] for error in errors))

    def test_combined_budget_does_not_read_past_remaining_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            source, budget, request, _ = self.relay_case(directory, byte_limit=5)
            with patch.object(fresh, "_request", side_effect=request):
                with fresh._range_relay(source, budget) as (url, logs, errors):
                    with urllib.request.urlopen(url, timeout=3) as response:
                        with self.assertRaises(http.client.IncompleteRead):
                            response.read()
            self.assertEqual(budget.download_bytes, 5)
            self.assertEqual(logs[0]["received_bytes"], 5)
            self.assertTrue(any("budget" in error["reason"] for error in errors))


if __name__ == "__main__":
    unittest.main()
