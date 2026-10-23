"""Scientific boundary tests; execute with binary dependencies on hosted CI only."""
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import pickle
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "evidence-sufficiency"))
import calibration
import study


def fixture_records():
    """Correlated depth derivatives, unknowns, and distinct genomic role blocks."""
    records = []
    roles = (("train", "train-a"), ("train", "train-b"),
             ("calibration", "cal-a"), ("threshold", "threshold-a"),
             ("test", "test-a"))
    for block_index, (role, block) in enumerate(roles):
        for index in range(42):
            label = None if index >= 40 else index % 2
            # Unknown candidates still have finite, plausible feature vectors.
            signal = index % 2
            for condition in ("original", "downsampled"):
                depth = (32 if condition == "original" else 10) + index % 3
                fraction = .24 + .1 * signal + .005 * (index % 5)
                # Keep a transferable ranking signal for the full-run boundary
                # test. A raw alt-count split can collapse every reduced-depth
                # score to one leaf, which the positive-slope guard rightly rejects.
                f0 = {"callable_depth": depth, "alt_fraction": fraction,
                      "constant": 1.0, "undefined": None}
                f1 = dict(f0, quality_difference=.2 + .6 * signal,
                          strand_fraction=.4 + .1 * (index % 2))
                row = {"candidate_id": f"{block}:{index}", "contig": "chr20",
                       "pos0": block_index * 1000000 + index * 20000,
                       "ref": "A", "alt": "C", "partition": role,
                       "block_id": block, "condition": condition,
                       "label": label, "label_reason": "fixture",
                       "f0": f0, "f1": f1}
                if condition == "original":
                    row["caller_score"] = 10.0 + 40.0 * signal
                records.append(row)
    return records


def experiment_protocol():
    return json.loads((HERE / "protocol.json").read_text())


def metadata_bundle(records):
    """Enough trusted-bundle metadata to exercise validation without fitting."""
    return {"schema": "rosalind-evidence-sufficiency-models-v1",
            "original_condition": "original", "feature_names": study._validate_records(records),
            "arms": {arm: {"feature_group": group, "learner": learner,
                           "operating_point": {"eligible": True, "threshold": .5}}
                     for arm, group, learner in study.ARMS}}


def gain_comparison():
    return {"matched_recall": {"absolute_fp_reduction_counts_minus_rich": 20,
                              "relative_fp_reduction_counts_minus_rich": .5,
                              "recall_rich_minus_counts": 0.0},
            "fixed_operating_points": {"recall_rich_minus_counts": 0.0},
            "bootstrap": {"eligible": True, "block_count": 20,
                          "requested_replicates": 2000, "valid_replicates": 2000,
                          "skipped_replicates": 0,
                          "absolute_fp_reduction_counts_minus_rich_95ci":
                              {"lower": 10.0, "median": 20.0, "upper": 30.0},
                          "fixed_threshold_recall_rich_minus_counts_95ci":
                              {"lower": -.004, "median": 0.0, "upper": .004}}}


class ThresholdBoundaryTests(unittest.TestCase):
    def test_highest_target_threshold_accepts_entire_tie(self):
        labels = [1, 1, 1, 0, 0]
        scores = [.9, .5, .5, .5, .1]
        point = study.select_score_threshold(labels, scores, target=.66)
        self.assertEqual(point["threshold"], .5)
        self.assertEqual(point["true_positives"], 3)
        self.assertEqual(point["false_positives"], 1)
        self.assertEqual(point["recall"], 1)
        self.assertGreater(point["recall"], point["target_recall"])

    def test_no_positive_class_is_ineligible_not_perfect(self):
        point = study.select_score_threshold([0, 0], [.8, .1])
        self.assertFalse(point["eligible"])
        self.assertNotIn("threshold", point)

    def test_unknowns_do_not_enter_supervised_metrics_or_become_negatives(self):
        predictions = [{"partition": "test", "condition": "downsampled", "arm": "C",
                        "label": label, "accepted": accepted, "probability": probability}
                       for label, accepted, probability in
                       [(1, True, .9), (0, True, .8), (None, True, .99)]]
        result = study.evaluate_predictions(predictions)["test"]["downsampled"]["C"]
        self.assertEqual(result["candidate_count"], 3)
        self.assertEqual(result["evaluated"], 2)
        self.assertEqual(result["unknown_labels"], 1)
        self.assertEqual(result["false_positives"], 1)
        self.assertAlmostEqual(result["brier_score"], (.1 ** 2 + .8 ** 2) / 2)

    def test_descriptive_threshold_is_explicitly_not_a_deployable_policy(self):
        result = calibration.descriptive_matched_recall([1, 1, 0], [.9, .5, .5], .99)
        self.assertTrue(result["threshold_selected_on_evaluation_labels"])
        self.assertFalse(result["deployable_threshold_validation"])
        self.assertEqual(result["false_positives"], 1)


class InputBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.records = fixture_records()
        self.bundle = metadata_bundle(self.records)

    def validate(self, records=None):
        return calibration._validate_inputs(self.records if records is None else records,
                                            self.bundle, "original", "downsampled")

    def test_counts_report_correlated_derivatives_and_unknowns_explicitly(self):
        result = self.validate()
        self.assertEqual(result["record_count_including_depth_derivatives"], 420)
        self.assertEqual(result["unique_candidate_count"], 210)
        self.assertEqual(result["unique_genomic_positions"], 210)
        self.assertFalse(result["paired_depth_derivatives_are_independent_samples"])
        for condition in ("original", "downsampled"):
            self.assertEqual(result["partition_condition_counts"]["test"][condition]["unknown"], 2)

    def test_missing_or_duplicate_derivative_is_rejected(self):
        for changed in (self.records[:-1], self.records + [copy.deepcopy(self.records[0])]):
            with self.subTest(length=len(changed)), self.assertRaises(ValueError):
                self.validate(changed)

    def test_depth_derivatives_cannot_cross_roles_or_change_labels(self):
        for field, value in (("partition", "test"), ("label", 1)):
            changed = copy.deepcopy(self.records)
            changed[1][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate(changed)

    def test_different_alts_at_one_locus_cannot_cross_roles(self):
        changed = copy.deepcopy(self.records)
        pair = copy.deepcopy(changed[:2])
        for row in pair:
            row.update(candidate_id="second-alt", alt="G", partition="test", block_id="test-extra")
        with self.assertRaisesRegex(ValueError, "crosses partition"):
            self.validate(changed + pair)

    def test_genomic_block_cannot_cross_roles_even_at_different_loci(self):
        changed = copy.deepcopy(self.records)
        for row in changed:
            if row["partition"] == "test":
                row["block_id"] = "train-a"
        with self.assertRaisesRegex(ValueError, "crosses partition"):
            self.validate(changed)

    def test_wrong_identity_or_unexpected_condition_is_rejected(self):
        for field, value in (("alt", "G"), ("condition", "new-depth")):
            changed = copy.deepcopy(self.records)
            changed[1][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate(changed)

    def test_feature_leakage_nonfinite_and_f0_disagreement_are_rejected(self):
        mutations = (("f1", "label", 1), ("f1", "quality_difference", float("nan")),
                     ("f1", "callable_depth", 999), ("f0", "alt_fraction", float("inf")))
        for group, key, value in mutations:
            changed = copy.deepcopy(self.records)
            changed[0][group][key] = value
            with self.subTest(group=group, key=key), self.assertRaises(ValueError):
                self.validate(changed)

    def test_missing_development_class_cannot_be_replaced_by_test_labels(self):
        changed = copy.deepcopy(self.records)
        for row in changed:
            if row["partition"] == "calibration" and row["label"] is not None:
                row["label"] = 0
        with self.assertRaisesRegex(ValueError, "both labelled classes"):
            self.validate(changed)


class GateBoundaryTests(unittest.TestCase):
    def gate(self, comparison=None, protocol=None):
        return calibration.evaluate_gain_gate(gain_comparison() if comparison is None else comparison,
                                             {"feature_storage": {"measured": True},
                                              "selected_inference": {"measured": True}},
                                             experiment_protocol() if protocol is None else protocol)

    def test_predictive_nomination_never_authorizes_p1_or_selective_escalation(self):
        result = self.gate()
        self.assertTrue(result["nominate_fresh_validation"])
        self.assertFalse(result["p1_authorized"])
        self.assertFalse(result["selective_escalation_authorized"])
        self.assertFalse(result["rich_evidence_total_cost_advantage_established"])
        self.assertTrue(result["exploratory_reused_test_only"])

    def test_positive_point_gain_with_missing_or_crossing_interval_cannot_pass(self):
        for interval in (None, {"lower": -1.0, "median": 20.0, "upper": 40.0},
                         {"lower": 0.0, "median": 20.0, "upper": 40.0}):
            comparison = gain_comparison()
            comparison["bootstrap"]["absolute_fp_reduction_counts_minus_rich_95ci"] = interval
            with self.subTest(interval=interval):
                self.assertFalse(self.gate(comparison)["nominate_fresh_validation"])

    def test_insufficient_bootstrap_support_cannot_pass_even_with_positive_interval(self):
        for key, value in (("eligible", False), ("block_count", 19), ("valid_replicates", 1899)):
            comparison = gain_comparison()
            comparison["bootstrap"][key] = value
            with self.subTest(key=key):
                self.assertFalse(self.gate(comparison)["nominate_fresh_validation"])

    def test_gate_uses_exact_protocol_keys_instead_of_silently_using_defaults(self):
        for key, value in (("minimum_relative_fp_reduction", .6),
                           ("minimum_absolute_fp_reduction", 21), ("maximum_recall_loss", .001)):
            protocol = experiment_protocol()
            protocol["gates"][key] = value
            with self.subTest(key=key):
                self.assertFalse(self.gate(protocol=protocol)["nominate_fresh_validation"])

    def test_fixed_threshold_recall_regression_or_uncertainty_blocks_nomination(self):
        for point, interval in ((-.006, {"lower": -.004, "upper": .001}),
                                (0.0, {"lower": -.006, "upper": .001}),
                                (0.0, {"lower": -.005, "upper": .001}), (0.0, None)):
            comparison = gain_comparison()
            comparison["fixed_operating_points"]["recall_rich_minus_counts"] = point
            comparison["bootstrap"]["fixed_threshold_recall_rich_minus_counts_95ci"] = interval
            with self.subTest(point=point, interval=interval):
                self.assertFalse(self.gate(comparison)["nominate_fresh_validation"])

    def test_tied_empirical_curve_cannot_win_by_materially_lower_attained_recall(self):
        comparison = gain_comparison()
        comparison["matched_recall"]["recall_rich_minus_counts"] = -.006
        self.assertFalse(self.gate(comparison)["nominate_fresh_validation"])

    def test_zero_denominator_and_small_absolute_gain_cannot_pass(self):
        for key, value in (("relative_fp_reduction_counts_minus_rich", None),
                           ("absolute_fp_reduction_counts_minus_rich", 9)):
            comparison = gain_comparison()
            comparison["matched_recall"][key] = value
            with self.subTest(key=key):
                self.assertFalse(self.gate(comparison)["nominate_fresh_validation"])


@unittest.skipUnless(importlib.util.find_spec("sklearn") and importlib.util.find_spec("scipy"),
                     "Binary dependencies required; run scientific checks on hosted validation")
class EstimatorBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dep = study._dependencies()

    def test_monotonic_recalibration_preserves_ranking_ties_and_target_decisions(self):
        np = self.dep["np"]
        raw = np.asarray([-3, -2, -2, -1, 0, 1, 2, 2, 3], dtype=float)
        labels = [0, 0, 1, 0, 0, 1, 0, 1, 1]
        calibrator = calibration._fit_sigmoid(raw, labels, self.dep)
        self.assertGreater(float(calibrator.coef_[0, 0]), 0)
        probabilities = calibrator.predict_proba(raw.reshape(-1, 1))[:, 1].tolist()
        self.assertEqual(probabilities[1], probabilities[2])
        self.assertEqual(probabilities[6], probabilities[7])
        self.assertEqual(study.average_precision(labels, raw.tolist()),
                         study.average_precision(labels, probabilities))
        for target in (.5, .95, .99):
            before = study.select_score_threshold(labels, raw.tolist(), target)
            after = study.select_score_threshold(labels, probabilities, target)
            self.assertEqual([s >= before["threshold"] for s in raw],
                             [p >= after["threshold"] for p in probabilities])

    def test_negative_or_zero_calibration_slope_is_not_silently_called_monotonic(self):
        np = self.dep["np"]
        for raw, labels in (([-2, -1, 1, 2], [1, 1, 0, 0]), ([0, 0, 0, 0], [0, 1, 0, 1])):
            with self.subTest(raw=raw), self.assertRaisesRegex(ValueError, "positive slope"):
                calibration._fit_sigmoid(np.asarray(raw, dtype=float), labels, self.dep)

    def bootstrap_fixture(self):
        rows, counts, rich, count_accepted, rich_accepted = [], [], [], [], []
        for block in range(12):
            for index, (label, c, r) in enumerate(zip([1, 1, 0, 0], [.9, .8, .85, .1], [.9, .8, .2, .1])):
                rows.append({"contig": "chr20", "pos0": block * 100000 + index, "label": label})
                counts.append(c)
                rich.append(r)
                count_accepted.append(c >= .8)
                rich_accepted.append(r >= .8)
        return rows, counts, rich, count_accepted, rich_accepted

    def bootstrap(self, values, **overrides):
        config = {"block_bp": 100000, "replicates": 200, "seed": 42} | overrides
        return calibration._paired_bootstrap(*values, .99, config, self.dep)

    def test_paired_window_bootstrap_is_deterministic_and_preserves_gain_direction(self):
        values = self.bootstrap_fixture()
        one = self.bootstrap(values)
        self.assertEqual(one, self.bootstrap(values))
        self.assertTrue(one["eligible"])
        self.assertEqual(one["block_count"], 12)
        self.assertEqual(one["valid_replicates"] + one["skipped_replicates"], 200)
        gain = one["absolute_fp_reduction_counts_minus_rich_95ci"]
        self.assertEqual((gain["lower"], gain["upper"]), (12.0, 12.0))
        rows, counts, rich, ca, ra = values
        reversed_result = self.bootstrap((rows, rich, counts, ra, ca))
        reverse_gain = reversed_result["absolute_fp_reduction_counts_minus_rich_95ci"]
        self.assertEqual((reverse_gain["lower"], reverse_gain["upper"]), (-12.0, -12.0))

    def test_zero_denominator_ratios_remain_unavailable_while_absolute_gain_survives(self):
        rows, _, rich, _, ra = self.bootstrap_fixture()
        result = self.bootstrap((rows, rich, rich, ra, ra))
        self.assertEqual(result["relative_fp_gain_valid_replicates"], 0)
        self.assertIsNone(result["relative_fp_reduction_counts_minus_rich_95ci"])
        self.assertEqual(result["absolute_fp_reduction_counts_minus_rich_95ci"]["lower"], 0)

    def test_empirical_target_is_reselected_inside_each_bootstrap_resample(self):
        rows = [{"contig": "chr20", "pos0": block * 100000 + index, "label": label}
                for block in range(2) for index, label in enumerate([1, 1, 0, 0])]
        counts = [.9, .8, .85, .1, .3, .2, .25, .1]
        rich = [.9, .8, .2, .1, .3, .2, .1, 0]
        result = self.bootstrap((rows, counts, rich, [s >= .2 for s in counts],
                                 [s >= .2 for s in rich]))
        gain = result["absolute_fp_reduction_counts_minus_rich_95ci"]
        # A mixed-window draw has one extra rich FP at the shared .2 threshold;
        # a repeated-window draw reselects that window's threshold and has none.
        self.assertEqual((gain["lower"], gain["upper"]), (1.0, 2.0))

    def test_bootstrap_reports_no_class_resamples_and_single_window_inadequacy(self):
        values = self.bootstrap_fixture()
        for row in values[0]:
            row["label"] = 0
        result = self.bootstrap(values)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["valid_replicates"], 0)
        self.assertEqual(result["skipped_replicates"], 200)
        self.assertIsNone(result["absolute_fp_reduction_counts_minus_rich_95ci"])
        one_window = tuple(v[:4] for v in self.bootstrap_fixture())
        self.assertFalse(self.bootstrap(one_window)["eligible"])

    def test_development_only_fits_freeze_before_test_and_preserve_unknown_predictions(self):
        records = fixture_records()
        for row in records:
            value = {"train": 1, "calibration": 1111, "threshold": 2222, "test": 9999}[row["partition"]]
            row["f0"]["constant"] = row["f1"]["constant"] = value
        protocol = experiment_protocol()
        # Scientific logic stays unchanged; lower replication count bounds this synthetic fixture.
        protocol["bootstrap"]["replicates"] = 200
        protocol["timing"]["rounds"] = 2
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            study.run_study(records, root / "original", {})
            original = study.load_bundle(root / "original")
            output = root / "one"
            score_function = calibration._scores
            observed_test_calls = []

            def score_after_freeze(entry, rows, *args, **kwargs):
                if any(row["partition"] == "test" for row in rows):
                    self.assertTrue((output / "frozen-models.pkl").exists())
                    self.assertTrue((output / "frozen-models.json").exists())
                    observed_test_calls.append(len(rows))
                return score_function(entry, rows, *args, **kwargs)

            with patch.object(calibration, "_scores", side_effect=score_after_freeze):
                first = calibration.run_experiment(records, original, output, protocol)
            self.assertTrue(observed_test_calls)
            frozen = json.loads((output / "frozen-models.json").read_text())
            model_bytes = (output / "frozen-models.pkl").read_bytes()
            self.assertEqual(hashlib.sha256(model_bytes).hexdigest(), frozen["model_sha256"])
            bundle = pickle.loads(model_bytes)  # Trusted fixture generated in this test process.
            for arm, group, _ in study.ARMS:
                names = bundle["feature_names"][group]
                fitted = bundle["arms"][arm + "_pooled"]["model"]
                self.assertEqual(fitted.named_steps["imputer"].statistics_[names.index("constant")], 1)
            self.assertFalse(first["selected_comparators"]["test_outcomes_used"])
            for group, key in (("f0", "counts"), ("f1", "rich")):
                candidates = [arm for arm, entry in bundle["arms"].items()
                              if arm.startswith("M") and entry["feature_group"] == group
                              and entry["variant"] != "frozen"
                              and entry["operating_points"]["downsampled"]["0.99"]["eligible"]]

                def selection_key(arm):
                    point = bundle["arms"][arm]["operating_points"]["downsampled"]["0.99"]
                    return point["false_positives"], -point["recall"], arm

                self.assertEqual(first["selected_comparators"][key], min(candidates, key=selection_key))
            self.assertNotIn("B1", first["metrics"]["test"]["downsampled"])
            self.assertIn("B1", first["metrics"]["test"]["original"])
            for condition in ("original", "downsampled"):
                for metric in first["metrics"]["test"][condition].values():
                    self.assertEqual(metric["candidate_count"], 42)
                    self.assertEqual(metric["unknown_labels"], 2)
                    self.assertEqual(metric["evaluated"], 40)
            predictions = [json.loads(line) for line in (output / "predictions.jsonl").read_text().splitlines()]
            self.assertEqual(len(predictions), len(records))
            self.assertEqual(sum(row["label"] is None for row in predictions), 20)
            for check in first["monotonic_calibration_checks"].values():
                self.assertTrue(check["raw_scores_equal"])
                self.assertTrue(check["all_target_decisions_equal"])
                self.assertGreater(check["recalibration_slope"], 0)
            changed = copy.deepcopy(records)
            for row in changed:
                if row["partition"] == "test" and row["label"] is not None:
                    row["label"] = 1 - row["label"]
            second = calibration.run_experiment(changed, original, root / "two", protocol)
            self.assertEqual(first["frozen_model"], second["frozen_model"])
            self.assertEqual(first["selected_comparators"], second["selected_comparators"])
            self.assertNotEqual(first["metrics"]["test"], second["metrics"]["test"])
            second_predictions = [json.loads(line) for line in
                                  (root / "two" / "predictions.jsonl").read_text().splitlines()]
            without_label = lambda rows: [{key: value for key, value in row.items() if key != "label"}
                                          for row in rows]
            self.assertEqual(without_label(predictions), without_label(second_predictions))
            for report in (first, second):
                self.assertFalse(report["gate"]["p1_authorized"])
                self.assertFalse(report["gate"]["selective_escalation_authorized"])
            with self.assertRaises(FileExistsError):
                calibration.run_experiment(records, original, output, protocol)


if __name__ == "__main__":
    unittest.main()
