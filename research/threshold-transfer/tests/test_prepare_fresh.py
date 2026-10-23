"""Hosted adapter checks: no real source downloads or model fits in fixtures."""
import copy
import hashlib
import http.client
import importlib.util
import io
from pathlib import Path
import sys
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


class FakeResponse(io.BytesIO):
    def __init__(self, status, headers, content=b""):
        super().__init__(content)
        self.status, self.headers = status, headers


class RangeRelayTests(unittest.TestCase):
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
