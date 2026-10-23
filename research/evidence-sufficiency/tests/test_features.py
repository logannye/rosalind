import copy
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FILTER_NAMES, SUM_NAMES, feature_vectors


def measured_row():
    row = {
        "contig": "chr20", "pos": 100, "ref": "A", "requested_alts": "C,G",
        "prefilter_depth": 13, "aligned_depth": 11, "callable_depth": 10,
        "a": 6, "c": 3, "g": 1, "t": 0,
        "a_fwd": 4, "a_rev": 2, "c_fwd": 1, "c_rev": 2,
        "g_fwd": 1, "g_rev": 0, "t_fwd": 0, "t_rev": 0,
        "allele_base_quality_sum": [180, 60, 25, 0],
        "allele_mapping_quality_sum": [360, 90, 40, 0],
        "allele_read_position_sum": [60, 90, 50, 0],
        "allele_read_length_sum": [600, 450, 100, 0],
    }
    for name in FILTER_NAMES:
        row["filtered_" + name] = 0
    row["filtered_duplicate"] = 2
    row["filtered_low_base_quality"] = 1
    for name in SUM_NAMES:
        row[name + "_sum"] = sum(row["allele_" + name + "_sum"])
    return row


class FeatureTests(unittest.TestCase):
    def test_independent_expected_ratios_and_alt_vector_index(self):
        result = feature_vectors(measured_row(), "C")
        f0, f1 = result["f0"], result["f1"]
        self.assertEqual(f0["callable_depth"], 10.0)
        self.assertEqual((f0["ref_count"], f0["alt_count"], f0["other_count"]), (6.0, 3.0, 1.0))
        self.assertEqual(f0["alt_fraction"], 0.3)
        self.assertEqual(f0["alt_fraction_undefined"], 0.0)
        self.assertAlmostEqual(f0["alt_count_log1p"], math.log(4))
        self.assertEqual(f1["alt_base_quality_mean"], 20.0)
        self.assertEqual(f1["ref_base_quality_mean"], 30.0)
        self.assertEqual(f1["alt_minus_ref_base_quality_mean"], -10.0)
        self.assertEqual(f1["alt_mapping_quality_mean"], 30.0)
        self.assertEqual(f1["alt_read_position_mean"], 30.0)
        self.assertEqual(f1["alt_read_length_mean"], 150.0)
        self.assertAlmostEqual(f1["alt_forward_fraction"], 1 / 3)
        self.assertAlmostEqual(f1["filtered_duplicate_fraction"], 2 / 13)
        self.assertAlmostEqual(f1["filtered_low_base_quality_fraction"], 1 / 13)
        self.assertEqual(f0["ref_A"], 1.0)
        self.assertEqual(f0["alt_C"], 1.0)
        self.assertIsNone(f1["t_forward_fraction"])
        self.assertEqual(f1["t_forward_fraction_undefined"], 1.0)
        self.assertEqual(feature_vectors(measured_row(), "G")["f1"]["alt_base_quality_mean"], 25.0)

    def test_reference_and_alt_roles_are_not_fixed_to_a_c(self):
        row = measured_row()
        row["ref"] = ord("G")
        result = feature_vectors(row, "A")
        self.assertEqual(result["f0"]["ref_count"], 1.0)
        self.assertEqual(result["f0"]["alt_count"], 6.0)
        self.assertEqual(result["f1"]["ref_mapping_quality_mean"], 40.0)
        self.assertEqual(result["f1"]["alt_mapping_quality_mean"], 60.0)

    def test_zero_observation_is_not_a_missing_row(self):
        row = measured_row()
        for key, value in list(row.items()):
            if isinstance(value, int):
                row[key] = 0
            elif isinstance(value, list):
                row[key] = [0] * 4
        result = feature_vectors(row, "C")
        self.assertEqual(result["f0"]["alt_count"], 0.0)
        self.assertIsNone(result["f0"]["alt_fraction"])
        self.assertEqual(result["f0"]["alt_fraction_undefined"], 1.0)
        self.assertIsNone(result["f1"]["base_quality_mean"])
        self.assertIsNone(result["f1"]["filtered_duplicate_fraction"])
        for missing in (None, {}, {**row, "observed": False}):
            with self.assertRaises(ValueError):
                feature_vectors(missing, "C")

    def test_feature_schema_stable_and_identity_independent(self):
        row = measured_row()
        first = feature_vectors(row, "C")
        row.update(contig="chr1", pos=5000, specimen_id="private-person", partition="test")
        self.assertEqual(first, feature_vectors(row, "C"))
        row["ref"] = "G"
        second = feature_vectors(row, "T")
        self.assertEqual(set(first["f0"]), set(second["f0"]))
        self.assertEqual(set(first["f1"]), set(second["f1"]))
        self.assertTrue(set(first["f0"]) <= set(first["f1"]))
        self.assertFalse(set(first["f1"]) & {"contig", "pos", "partition", "specimen_id"})

    def test_projection_missingness_and_invalid_integer_refused(self):
        for column in ("a", "a_fwd", "read_length_sum", "allele_base_quality_sum", "filtered_duplicate"):
            row = measured_row()
            del row[column]
            with self.assertRaises(ValueError):
                feature_vectors(row, "C")
        for invalid in (True, 1.0, -1, 1 << 64, None):
            row = measured_row()
            row["callable_depth"] = invalid
            with self.assertRaises(ValueError):
                feature_vectors(row, "C")

    def test_inconsistent_observation_groups_refused(self):
        for column in ("a_fwd", "a", "filtered_duplicate", "base_quality_sum"):
            row = measured_row()
            row[column] += 1
            with self.assertRaises(ValueError):
                feature_vectors(row, "C")
        row = measured_row()
        row["allele_base_quality_sum"] = [1, 2, 3]
        with self.assertRaises(ValueError):
            feature_vectors(row, "C")

    def test_input_not_mutated_and_large_raw_integers_retained(self):
        row = measured_row()
        before = copy.deepcopy(row)
        feature_vectors(row, "C")
        self.assertEqual(row, before)
        # Deliberate float model features do not overwrite exact integer input.
        large = (1 << 53) + 1
        row = measured_row()
        for key, value in list(row.items()):
            if isinstance(value, int):
                row[key] = 0
            elif isinstance(value, list):
                row[key] = [0] * 4
        row.update(prefilter_depth=large, aligned_depth=large, callable_depth=large,
                   a=large, a_fwd=large)
        row["read_length_sum"] = large
        row["allele_read_length_sum"] = [large, 0, 0, 0]
        result = feature_vectors(row, "C")
        self.assertEqual(row["a"], large)
        self.assertTrue(math.isfinite(result["f0"]["callable_depth_log1p"]))

    def test_non_snv_and_unknown_reference_refused(self):
        for ref, alt in (("N", "C"), ("A", "A"), ("A", "CT"), ("A", "<DEL>")):
            row = measured_row()
            row["ref"] = ref
            with self.assertRaises(ValueError):
                feature_vectors(row, alt)


if __name__ == "__main__":
    unittest.main()
