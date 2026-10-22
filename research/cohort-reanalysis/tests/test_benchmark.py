"""Independent pilot gates without native execution or optional binary packages."""
from __future__ import annotations

import contextlib
import json
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import benchmark


def read(position, base="A", *, flag=0, mapq=60, bq=30, cigar="1M"):
    return SimpleNamespace(reference_start=position - 1, query_sequence=base,
        query_qualities=[bq], cigarstring=cigar, flag=flag, mapping_quality=mapq)


class IndependentOracleTests(unittest.TestCase):
    def independent(self, reads, stored):
        fake = SimpleNamespace(AlignmentFile=lambda *_a, **_kw: contextlib.nullcontext(iter(reads)))
        with patch.dict(sys.modules, {"pysam": fake}):
            return list(benchmark.independent_rows({"bam": "same-generated.bam", "stored": stored}))

    def test_coverage_mask_stored_zero_and_exclusions(self):
        reads = [read(10, "A"), read(10, "C"), read(10, "C", flag=1024),
                 read(10, "G", flag=256), read(10, "G", flag=2048),
                 read(10, "G", flag=512), read(10, "T", mapq=19),
                 read(10, "G", mapq=255), read(10, "C", bq=19), read(10, "N"),
                 read(20, "T"), read(40, "G", bq=10)]
        actual = self.independent(reads, [10, 30, 40, 50])
        self.assertEqual([r["pos"] for r in actual], [10, 30, 40, 50])
        self.assertEqual(actual[0], {"pos": 10, "callable_depth": 2, "a": 1, "c": 1, "g": 0, "t": 0})
        self.assertEqual([r["callable_depth"] for r in actual], [2, 0, 0, 0])
        self.assertNotIn(20, [r["pos"] for r in actual], "unmeasured nonzero source data must remain missing")

    def test_empty_source_preserves_only_the_stored_coordinates(self):
        self.assertEqual(self.independent([], []), [])
        rows = self.independent([], [5, 9])
        self.assertEqual([(r["pos"], r["callable_depth"]) for r in rows], [(5, 0), (9, 0)])

    def test_oracle_rejects_unsupported_cigar(self):
        with self.assertRaisesRegex(ValueError, "only supports 1M"):
            self.independent([read(10, "AA", cigar="2M")], [10])


class BaselineSummaryTests(unittest.TestCase):
    def test_independent_reducer_matches_unchanged_authored_fixture(self):
        fixture = json.loads((HERE.parents[1] / "examples/cohort-reanalysis/expected.json").read_text())
        rows = []
        for text, alt in fixture["candidates"].items():
            position = int(text)
            request = {"candidates": [(position, alt)], "alt": alt, "thresholds": [10]}
            counters = benchmark.counter_rows(request, 3)
            for member in fixture["members"].values():
                if position not in member["stored_positions"]:
                    continue
                full = next(row for row in member["full_rows"] if row[0] == position)
                benchmark.add(counters, request, {"pos": position, "callable_depth": full[3],
                                                 **dict(zip("acgt", full[4:]))})
            rows.extend(benchmark.finish(counters, 3))
        columns = ["pos", "alt", "n_requested", "n_observed", "n_depth_eligible", "n_alt_supported",
                   "total_callable_depth", "total_alt_count", "eligible_callable_depth", "eligible_alt_count"]
        self.assertEqual([[row[column] for column in columns] for row in rows], fixture["partial_summary"])

    def test_threshold_screens_do_not_replace_observed_denominator(self):
        request = {"candidates": [(10, "G"), (20, "G")], "alt": "G", "thresholds": [5, 10, 20]}
        counters = benchmark.counter_rows(request, 3)
        for depth, alt in [(0, 0), (9, 3), (10, 0)]:
            benchmark.add(counters, request, {"pos": 10, "callable_depth": depth, "g": alt})
        rows = benchmark.finish(counters, 3)
        self.assertEqual([r["n_observed"] for r in rows[:3]], [3, 3, 3])
        self.assertEqual([r["n_depth_eligible"] for r in rows[:3]], [2, 1, 0])
        self.assertEqual([r["n_alt_supported"] for r in rows[:3]], [1, 0, 0])
        self.assertEqual([r["total_callable_depth"] for r in rows[:3]], [19, 19, 19])
        self.assertEqual(rows[3]["n_observed"], 0)
        self.assertEqual(rows[3]["reference_validation"], "unavailable")
        self.assertIsNone(rows[3]["alt_fraction_numerator"])
        self.assertIsNone(rows[3]["alt_fraction_denominator"])

    def test_counts_above_double_precision_and_overflow_are_exact(self):
        request = {"candidates": [(10, "C")], "alt": "C", "thresholds": [10]}
        counters = benchmark.counter_rows(request, 1)
        depth = (1 << 53) + 1
        benchmark.add(counters, request, {"pos": 10, "callable_depth": depth, "c": depth - 1})
        result = benchmark.finish(counters, 1)[0]
        self.assertEqual(result["total_callable_depth"], depth)
        self.assertEqual(result["total_alt_count"], depth - 1)
        with self.assertRaises(OverflowError):
            benchmark.add(counters, request, {"pos": 10, "callable_depth": (1 << 64) - 1, "c": 0})
        for value in [1.0, True, -1, 1 << 64]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                benchmark.add(benchmark.counter_rows(request, 1), request,
                              {"pos": 10, "callable_depth": value, "c": 0})
        with self.assertRaises(ValueError):
            benchmark.add(benchmark.counter_rows(request, 1), request,
                          {"pos": 10, "callable_depth": 10, "c": -1})

    def test_worker_only_reveals_current_request_and_honors_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            requests = [{"id": name} for name in ("first", "second", "unrevealed")]
            state = {"requests": requests, "members": [], "config": {}}
            path = root / "state.json"
            path.write_text(json.dumps(state))
            calls = []

            def query(visible, binary, request, destination):
                self.assertNotIn("requests", visible)
                calls.append(request["id"])
                return [{"request": request["id"]}], {}

            args = SimpleNamespace(state=path, worker="saved", output=root / "result", prefix=2, binary="binary")
            with patch.object(benchmark, "saved_request", side_effect=query):
                benchmark.worker(args)
            self.assertEqual(calls, ["first", "second"])
            self.assertEqual(json.loads((root / "result/answers.json").read_text()),
                             [[{"request": "first"}], [{"request": "second"}]])

    def test_saved_baseline_matches_masks_by_specimen_not_list_position(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "native").mkdir()
            (root / "native/members.json").write_text(json.dumps({"version": 1, "members": [
                {"specimen_id": "b", "dataset": "b-leaf"},
                {"specimen_id": "a", "dataset": "a-leaf"}]}))
            state = {"root": str(root), "config": {"reference_length": 64, "native_memory_budget_mb": 256},
                     "members": [{"specimen_id": "a", "stored": [10]}, {"specimen_id": "b", "stored": [20]}]}
            request = {"candidates": [(10, "C"), (20, "C")], "alt": "C", "thresholds": [10]}
            calls = []

            class Stream:
                result = True

                def __init__(self, leaf, selection):
                    positions = [int(line.split("\t")[1]) for line in selection.read_text().splitlines()
                                 if not line.startswith("#")]
                    calls.append((leaf, positions))
                    self.rows = [{"pos": pos, "callable_depth": 10, "c": 1} for pos in positions]

                def __enter__(self):
                    return self

                def __exit__(self, *_):
                    pass

                def __iter__(self):
                    yield SimpleNamespace(to_pylist=lambda: self.rows)

            def open_dataset(leaf, **_kwargs):
                return SimpleNamespace(batches=lambda **kwargs: Stream(leaf, kwargs["sites"]))

            fake = SimpleNamespace(open_dataset=open_dataset)
            with patch.dict(sys.modules, {"rosalind": fake}):
                answer, _metrics = benchmark.saved_request(state, "binary", request, root)
            self.assertEqual(calls, [("b-leaf", [20]), ("a-leaf", [10])])
            self.assertEqual([r["n_observed"] for r in answer], [1, 1])
            self.assertEqual([r["total_callable_depth"] for r in answer], [10, 10])

    def test_columnar_projects_only_requested_alt_and_performance_mode_does_not_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "columnar").mkdir()
            (root / "columnar/inventory.json").write_text(json.dumps({"files": [
                {"path": "selected.parquet", "member": "a", "min_pos": 10, "max_pos": 10}]}))
            state = {"root": str(root), "members": [{"specimen_id": "a"}]}
            request = {"candidates": [(10, "G")], "alt": "G", "thresholds": [10]}
            options = []

            def scanner(**kwargs):
                options.append(kwargs)
                batch = SimpleNamespace(to_pylist=lambda: [{"pos": 10, "callable_depth": 10, "g": 2}])
                return SimpleNamespace(to_batches=lambda: iter([batch]))

            dataset_module = ModuleType("pyarrow.dataset")
            dataset_module.dataset = lambda *_a, **_kw: SimpleNamespace(scanner=scanner)
            dataset_module.field = lambda _name: SimpleNamespace(isin=lambda values: tuple(values))
            arrow_module = ModuleType("pyarrow")
            arrow_module.dataset = dataset_module
            with patch.dict(sys.modules, {"pyarrow": arrow_module, "pyarrow.dataset": dataset_module}):
                with patch.object(benchmark, "sha256", side_effect=AssertionError("performance path must not hash")):
                    answer, metrics = benchmark.columnar_request(state, request, root, False)
            self.assertEqual(options[0]["columns"], ["pos", "callable_depth", "g"])
            self.assertEqual(options[0]["batch_size"], 1024)
            self.assertEqual(options[0]["fragment_readahead"], 1)
            self.assertEqual(answer[0]["total_alt_count"], 2)
            self.assertEqual(metrics["bytes_hashed"], 0)

    def test_source_identity_includes_exact_bam_index_reference_and_mask(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference, bam, bed = root / "ref.fa", root / "reads.bam", root / "stored.bed"
            paths = [reference, Path(str(reference) + ".fai"), bam, Path(str(bam) + ".bai"), bed]
            for index, path in enumerate(paths):
                path.write_bytes(bytes([index]))
            state = {"reference": str(reference), "members": [{"bam": str(bam), "regions": str(bed)}]}
            before = benchmark.source_identities(state)
            self.assertEqual(set(before), set(map(str, paths)))
            bed.write_text("changed original stored mask")
            after = benchmark.source_identities(state)
            self.assertNotEqual(before[str(bed)], after[str(bed)])
            self.assertEqual(before[str(bam)], after[str(bam)])


if __name__ == "__main__":
    unittest.main()
