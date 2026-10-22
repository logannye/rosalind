"""Independent synthetic and unchanged fixture-oracle checks of exact reduction."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from contract import Candidate, ContractError, UINT64_MAX
from reduce import WindowReducer, checked_add


def candidate(pos=10, alt="C"):
    return Candidate(0, pos - 1, alt, "synthetic", "A")


def observation(pos=10, depth=10, alt=1, **kwargs):
    return {"contig": "synthetic", "pos": pos, "ref": "A", "callable_depth": depth,
            "a": depth - alt, "c": alt, "g": 0, "t": 0, **kwargs}


class ReducerTests(unittest.TestCase):
    def test_unchanged_three_specimen_oracle(self):
        path = Path(__file__).resolve().parents[3] / "examples/cohort-reanalysis/expected.json"
        expected = json.loads(path.read_text())
        reducer = WindowReducer([candidate(int(pos), alt) for pos, alt in expected["candidates"].items()],
                                [10], 3)
        for specimen, member in sorted(expected["members"].items()):
            reducer.begin_member(specimen)
            for pos, prefilter, aligned, depth, a, c, g, t in member["full_rows"]:
                if pos in member["stored_positions"]:
                    reducer.observe({"contig": "synthetic", "pos": pos, "ref": "A",
                                     "callable_depth": depth, "a": a, "c": c, "g": g, "t": t})
            self.assertEqual(reducer.end_member(), len(member["stored_positions"]))
        result = reducer.finish()
        columns = ("pos", "alt", "n_requested", "n_observed", "n_depth_eligible", "n_alt_supported",
                   "total_callable_depth", "total_alt_count", "eligible_callable_depth", "eligible_alt_count")
        self.assertEqual([[row[key] for key in columns] for row in result], expected["partial_summary"])

    def test_missing_zero_and_positive_below_threshold_remain_distinct(self):
        reducer = WindowReducer([candidate(10), candidate(20), candidate(30)], [10], 2)
        reducer.begin_member("a")
        reducer.observe(observation(10, 0, 0))
        reducer.observe(observation(20, 5, 0))
        reducer.end_member()
        reducer.begin_member("b")
        reducer.end_member()
        zero, low, absent = reducer.finish()
        self.assertEqual((zero["n_observed"], low["n_observed"], absent["n_observed"]), (1, 1, 0))
        self.assertEqual((zero["total_callable_depth"], low["total_callable_depth"], absent["total_callable_depth"]), (0, 5, 0))
        self.assertIsNone(zero["alt_fraction_numerator"])
        self.assertEqual((low["alt_fraction_numerator"], low["alt_fraction_denominator"]), (0, 5))
        self.assertIsNone(absent["alt_fraction_denominator"])
        self.assertEqual(zero["reference_validation"], "checked_on_observed_rows")
        self.assertEqual(absent["reference_validation"], "unavailable")
        with self.assertRaisesRegex(ContractError, "strict missingness"):
            reducer.finish(missingness="strict")

    def test_counts_above_2_to_53_never_pass_through_float(self):
        depth = (1 << 53) + 1
        reducer = WindowReducer([candidate()], [10], 2)
        for specimen, observed in (("a", observation(depth=depth, alt=depth - 1)),
                                   ("b", observation(depth=2, alt=1))):
            reducer.begin_member(specimen)
            reducer.observe(observed)
            reducer.end_member()
        row = reducer.finish(missingness="strict")[0]
        self.assertEqual(row["total_callable_depth"], (1 << 53) + 3)
        self.assertEqual(row["total_alt_count"], (1 << 53) + 1)
        self.assertEqual(json.loads(json.dumps(row))["total_callable_depth"], (1 << 53) + 3)

    def test_sum_overflow_refuses_without_partial_update(self):
        reducer = WindowReducer([candidate()], [10], 2)
        reducer.begin_member("a")
        reducer.observe(observation(depth=UINT64_MAX, alt=UINT64_MAX))
        reducer.end_member()
        reducer.begin_member("b")
        with self.assertRaisesRegex(ContractError, "overflow"):
            reducer.observe(observation(depth=1, alt=1))
        # Catching the error cannot turn a failed run into a partial success.
        with self.assertRaisesRegex(ContractError, "prevents summary publication"):
            reducer.finish()
        self.assertEqual(checked_add(UINT64_MAX, 0), UINT64_MAX)

    def test_malformed_counts_ref_mismatch_and_duplicate_observations_fail(self):
        invalid = [observation(depth=1.0), observation(callable_depth=True),
                   observation(c=-1), observation(c=UINT64_MAX + 1),
                   observation(a=1), observation(ref="T"), observation(pos=30)]
        for value in invalid:
            with self.subTest(value=value):
                reducer = WindowReducer([candidate()], [10], 1)
                reducer.begin_member("a")
                with self.assertRaises(ContractError):
                    reducer.observe(value)
        reducer = WindowReducer([candidate()], [10], 1)
        reducer.begin_member("a")
        reducer.observe(observation())
        with self.assertRaisesRegex(ContractError, "duplicate observed"):
            reducer.observe(observation())

    def test_multiple_alts_and_thresholds_repeat_but_do_not_combine_denominators(self):
        reducer = WindowReducer([candidate(10, "T"), candidate(10, "C")], [20, 5, 10], 2)
        reducer.begin_member("a")
        reducer.observe(observation(depth=10, alt=2, a=5, t=3))
        reducer.end_member()
        reducer.begin_member("b")
        reducer.observe(observation(depth=5, alt=0))
        reducer.end_member()
        rows = reducer.finish(missingness="strict")
        self.assertEqual([(r["alt"], r["min_callable_depth"]) for r in rows],
                         [("C", 5), ("C", 10), ("C", 20), ("T", 5), ("T", 10), ("T", 20)])
        self.assertEqual([r["total_callable_depth"] for r in rows], [15] * 6)
        self.assertEqual([r["total_alt_count"] for r in rows], [2] * 3 + [3] * 3)
        self.assertEqual([r["n_depth_eligible"] for r in rows], [2, 1, 0] * 2)
        self.assertIsNone(rows[2]["eligible_alt_fraction_denominator"])

    def test_window_partitioning_has_identical_scientific_output(self):
        def run(candidates):
            reducer = WindowReducer(candidates, [5, 10], 2)
            wanted = {c.pos for c in candidates}
            for specimen, depths in (("a", [0, 5, 10]), ("b", [5, 9, 11])):
                reducer.begin_member(specimen)
                for pos, depth in zip((10, 20, 30), depths):
                    if pos in wanted:
                        reducer.observe(observation(pos, depth, min(depth, 1)))
                reducer.end_member()
            return reducer.finish(missingness="strict")
        candidates = [candidate(pos) for pos in (10, 20, 30)]
        self.assertEqual(run(candidates), [row for c in candidates for row in run([c])])

    def test_empty_request_and_completion_guard(self):
        reducer = WindowReducer([], [10], 1)
        with self.assertRaisesRegex(ContractError, "every requested member"):
            reducer.finish()
        reducer.begin_member("a")
        reducer.end_member()
        self.assertEqual(reducer.finish(missingness="strict"), [])
        with self.assertRaisesRegex(ContractError, "current reducer state"):
            reducer.begin_member("b")


if __name__ == "__main__":
    unittest.main()
