import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from labels import TruthIndex


def record(pos0, ref="A", alts=("C",), gt=(0, 1), **extra):
    return {"contig": "chr20", "pos0": pos0, "ref": ref, "alts": alts, "gt": gt, **extra}


class TruthLabelTests(unittest.TestCase):
    def test_genotype_aware_multiallelic_labels(self):
        index = TruthIndex.from_records([record(50, alts=("C", "G"), gt=(0, 2))], [("chr20", 0, 100)])
        self.assertEqual(index.label("chr20", 50, "A", "G"), (1, "benchmark_supported_alt"))
        self.assertEqual(index.label("chr20", 50, "A", "C"), (0, "benchmark_region_unmatched_alt"))
        self.assertEqual(index.positives, {("chr20", 50, "A", "G")})

    def test_homozygous_reference_is_not_a_positive_alt(self):
        index = TruthIndex.from_records([record(50, gt=(0, 0))], [("chr20", 0, 100)])
        self.assertEqual(index.label("chr20", 50, "A", "C")[0], 0)
        self.assertEqual(index.positives, set())

    def test_benchmark_and_selection_half_open_boundaries(self):
        index = TruthIndex.from_records([], [("chr20", 10, 20)], [("chr20", 5, 30)])
        self.assertEqual(index.label("chr20", 10, "A", "C")[0], 0)
        self.assertEqual(index.label("chr20", 19, "A", "C")[0], 0)
        self.assertEqual(index.label("chr20", 20, "A", "C"), (None, "outside_benchmark_regions"))
        self.assertEqual(index.label("chr20", 4, "A", "C"), (None, "outside_selected_intervals"))

    def test_partial_missing_invalid_and_filtered_gt_exclusion(self):
        for rec in (
            record(50, gt=(0, None)), record(50, gt=None), record(50, gt=(0, 2)),
            record(50, filters=("LowQual",)),
        ):
            index = TruthIndex.from_records([rec], [("chr20", 0, 100)])
            self.assertEqual(index.label("chr20", 40, "A", "C")[1], "excluded_truth_event")
            self.assertEqual(index.label("chr20", 60, "A", "C")[1], "excluded_truth_event")
            self.assertEqual(index.label("chr20", 39, "A", "C")[0], 0)
            self.assertEqual(index.label("chr20", 61, "A", "C")[0], 0)
            self.assertEqual(index.positives, set())

    def test_complex_span_end_and_neighboring_positive_excluded(self):
        records = [record(120), record(100, alts=("<DEL>",), end=130)]
        index = TruthIndex.from_records(records, [("chr20", 0, 200)])
        self.assertEqual(index.excluded_intervals, [{"contig": "chr20", "start": 90, "end": 140}])
        self.assertEqual(index.label("chr20", 139, "A", "C")[0], None)
        self.assertEqual(index.label("chr20", 140, "A", "C")[0], 0)
        self.assertNotIn(("chr20", 120, "A", "C"), index.positives)
        mnp = TruthIndex.from_records([record(40, ref="AC", alts=("GT",))], [("chr20", 0, 100)])
        self.assertEqual(mnp.excluded_intervals[0]["end"], 52)

    def test_event_starting_outside_selection_still_excludes_its_overlap(self):
        index = TruthIndex.from_records(
            [record(70, alts=("<DEL>",), end=110), record(500)],
            [("chr20", 0, 1000)], [("chr20", 100, 105)],
        )
        self.assertEqual(index.records_considered, 1)
        self.assertEqual(index.label("chr20", 102, "A", "C")[0], None)

    def test_optional_reference_validation_and_ref_mismatch(self):
        def fetch(contig, start, end):
            if contig != "chr20":
                raise KeyError(contig)
            return "A" * (end - start)
        index = TruthIndex.from_records(
            [record(50, ref="G", alts=("T",))], [("chr20", 0, 100)], reference_fetch=fetch,
        )
        self.assertEqual(index.exclusion_counts["unresolved_truth_reference"], 1)
        self.assertEqual(index.label("chr20", 50, "A", "C")[0], None)
        self.assertEqual(index.label("chr20", 20, "G", "C")[1], "candidate_reference_mismatch")
        index = TruthIndex.from_records([record(50)], [("chr20", 0, 100)])
        self.assertEqual(index.label("chr20", 50, "G", "C")[1], "candidate_reference_mismatch")

    def test_positive_denominators_respect_masks_and_partitions(self):
        intervals = [
            {"contig": "chr20", "start": 0, "end": 50, "partition": "train", "block_id": "t"},
            {"contig": "chr20", "start": 50, "end": 100, "partition": "test", "block_id": "e"},
        ]
        index = TruthIndex.from_records([record(10), record(60), record(90)], [("chr20", 0, 80)], intervals)
        self.assertEqual(index.truth_keys("train"), {("chr20", 10, "A", "C")})
        self.assertEqual(index.truth_keys("test"), {("chr20", 60, "A", "C")})
        self.assertEqual([row["n_positive"] for row in index.positive_counts_by_interval()], [1, 1])

    def test_distant_records_streamed_without_retention(self):
        def records():
            for pos in range(1000, 2000):
                yield record(pos)
            yield record(20)
        index = TruthIndex.from_records(records(), [("chr20", 0, 3000)], [("chr20", 0, 100)])
        self.assertEqual(index.records_considered, 1)
        self.assertEqual(index.positives, {("chr20", 20, "A", "C")})
        self.assertEqual(len(index._truth_refs), 1)

    def test_domain_report_intersects_and_subtracts_interval_masks(self):
        intervals = [
            {"contig": "chr20", "start": 10, "end": 50, "partition": "train", "block_id": "a"},
            {"contig": "chr20", "start": 50, "end": 90, "partition": "test", "block_id": "b"},
            {"contig": "chr21", "start": 0, "end": 1_000_000_000, "block_id": "empty"},
        ]
        index = TruthIndex.from_records(
            [record(20), record(42), record(45, ref="AC", alts=("GT",)),
             record(75, alts=("C", "G"), gt=(1, 2))],
            [("chr20", 0, 30), ("chr20", 40, 60), ("chr20", 70, 100)], intervals,
        )
        report = index.domain_report()
        self.assertEqual(report["exclusion_padding_bp"], 10)
        self.assertEqual(report["records_considered"], 4)
        self.assertEqual(report["exclusion_counts"], {"complex_or_unsupported_record": 1})
        self.assertEqual(report["excluded_intervals"], [{"contig": "chr20", "start": 35, "end": 57}])
        train, test, empty = report["blocks"]
        self.assertEqual((train["selected_bases"], train["benchmark_confident_bases"],
                          train["excluded_confident_bases"], train["assessed_bases"],
                          train["assessed_truth_snvs"]), (40, 30, 10, 20, 1))
        self.assertEqual((test["selected_bases"], test["benchmark_confident_bases"],
                          test["excluded_confident_bases"], test["assessed_bases"],
                          test["assessed_truth_snvs"]), (40, 30, 7, 23, 2))
        self.assertEqual(empty["selected_bases"], 1_000_000_000)
        self.assertEqual(empty["assessed_bases"], 0)
        self.assertEqual(empty["assessed_truth_snvs"], 0)

    def test_unsupported_candidates_and_contradictory_truth(self):
        index = TruthIndex.from_records([record(50), record(50, ref="G")], [("chr20", 0, 100)])
        self.assertEqual(index.positives, set())
        for pos, ref, alt in ((-1, "A", "C"), (5, "N", "C"), (5, "A", "A"), (5, "A", "AC")):
            self.assertEqual(index.label("chr20", pos, ref, alt), (None, "unsupported_candidate"))


if __name__ == "__main__":
    unittest.main()
