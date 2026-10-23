"""Independent policy and spatial-uncertainty boundaries; run on hosted CI only.

These fixtures do not fit estimators or reuse HG001 outcomes. Synthetic bootstrap
sizes are intentionally small; production settings remain frozen in protocol.json.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "depth-calibration"))
sys.path.insert(0, str(HERE.parent / "evidence-sufficiency"))

import thresholds


def synthetic_protocol():
    return {
        "target_recall": 0.99,
        "confidence_level": 0.95,
        "bootstrap": {
            "block_bp": 100_000,
            "development_replicates": 80,
            "fresh_replicates": 80,
            "seed": 42,
            "minimum_nonempty_blocks": 3,
            "minimum_valid_fraction": 0.9,
            "minimum_error_blocks": 2,
        },
        "gates": {
            "confidence_level": 0.95,
            "minimum_candidate_recall": 0.99,
            "maximum_recall_loss": 0.005,
        },
    }


def spatial_rows(blocks=6, positives=20, negatives=3,
                 conditions=("original", "downsampled")):
    """Identical loci and outcomes across conditions, grouped by physical window."""
    rows = []
    for block in range(blocks):
        for index in range(positives + negatives):
            for condition in conditions:
                rows.append({
                    "candidate_id": f"candidate-{block}-{index}",
                    "contig": "chr20", "pos0": block * 100_000 + index,
                    "ref": "A", "alt": "C", "block_id": f"block-{block}",
                    "partition": "threshold", "condition": condition,
                    "label": int(index < positives),
                })
    return rows


def support(result, policy="policy", condition="original"):
    return result["policies"][policy]["conditions"][condition]


class FrozenScoreModel:
    """Deterministic already-fitted stand-in; fitting and test leakage are errors."""
    def __init__(self, column):
        self.column = column

    def fit(self, *_args, **_kwargs):
        raise AssertionError("The threshold experiment must never refit a base model")

    def predict_proba(self, matrix):
        import numpy as np
        score = matrix[:, self.column]
        if np.any(np.abs(score) > 100):
            raise AssertionError("Forbidden test-role features were scored during policy freezing")
        probability = 1 / (1 + np.exp(-score))
        return np.column_stack((1 - probability, probability))


def policy_fixture():
    records = []
    for role_index, role in enumerate(("calibration", "threshold", "test")):
        for block in range(6):
            for index in range(24):
                label = None if index == 23 else int(index < 20)
                score = 1.0 + index / 100 if label == 1 else -1.0 - index / 100
                for condition in ("original", "downsampled"):
                    depth = (30 if condition == "original" else 10) + index % 3
                    f0 = {"callable_depth": float(depth), "score": score}
                    f1 = dict(f0, rich_score=score + 0.1)
                    records.append({
                        "candidate_id": f"{role}-{block}-{index}", "contig": "chr20",
                        "pos0": role_index * 10_000_000 + block * 100_000 + index,
                        "ref": "A", "alt": "C", "block_id": f"{role}-{block}",
                        "partition": role, "condition": condition,
                        "label": label, "label_reason": "synthetic", "f0": f0, "f1": f1,
                    })
    names = {"f0": ["callable_depth", "score"],
             "f1": ["callable_depth", "rich_score", "score"]}
    bundle = {"schema": "rosalind-depth-calibration-models-v1", "feature_names": names, "arms": {}}
    for model, group in (("M10_pooled", "f0"), ("M11_pooled", "f1")):
        bundle["arms"][model] = {
            "feature_group": group, "variant": "pooled", "model": FrozenScoreModel(1),
            "operating_points": {condition: {"0.99": {"eligible": True, "threshold": 0.0}}
                                 for condition in ("original", "downsampled")},
        }
    protocol = synthetic_protocol()
    protocol.update({
        "development_conditions": ["original", "downsampled"],
        "required_fresh_conditions": ["original", "downsampled", "novel_low"],
        "primary_condition": "downsampled",
        "condition_sources": {"original": "original", "downsampled": "downsampled", "novel_low": "downsampled"},
        "policy_recall_targets": [0.99, 0.995, 0.999, 1.0],
        "depth_edges": [0, 5, 10, 20, 40, None],
        "minimum_calibration_positives_per_depth_bin": 300,
        "low_depth_passthrough_below": 5,
        "timing": {"rounds": 3, "repeated_queries": [1, 5, 20]},
    })
    return records, bundle, protocol


def fresh_fixture(records):
    fresh = [copy.deepcopy(r) for r in records if r["partition"] == "test"]
    for row in list(fresh):
        row["sample_id"] = "synthetic-fresh-specimen"
        if row["condition"] == "downsampled":
            novel = copy.deepcopy(row)
            novel["condition"] = "novel_low"
            novel["f0"]["callable_depth"] = novel["f1"]["callable_depth"] = 4.0
            fresh.append(novel)
    return fresh


class PolicyDecisionTests(unittest.TestCase):
    def test_retain_all_is_explicit_unresolved_work(self):
        policy = {"id": "retain", "kind": "retain_all"}
        result = thresholds.policy_decision(policy, -100.0, 12.0, "original")
        self.assertTrue(result["retained"])
        self.assertTrue(result["deferred"])
        self.assertTrue(result["reason"])

    def test_global_threshold_accepts_equality_without_deferral(self):
        policy = {"id": "global", "kind": "global", "threshold": 0.25,
                  "passthrough_below": None}
        below = thresholds.policy_decision(policy, 0.249, 10.0, "original")
        equal = thresholds.policy_decision(policy, 0.25, 10.0, "original")
        above = thresholds.policy_decision(policy, 0.251, 10.0, "original")
        self.assertFalse(below["retained"])
        self.assertTrue(equal["retained"])
        self.assertTrue(above["retained"])
        self.assertFalse(equal["deferred"])

    def test_low_depth_is_passed_through_and_charged_as_deferred(self):
        policy = {"id": "guarded", "kind": "global", "threshold": 5.0,
                  "passthrough_below": 5}
        low = thresholds.policy_decision(policy, -100.0, 4.0, "original")
        boundary = thresholds.policy_decision(policy, -100.0, 5.0, "original")
        self.assertTrue(low["retained"])
        self.assertTrue(low["deferred"])
        self.assertFalse(boundary["retained"])
        self.assertFalse(boundary["deferred"])

    def test_missing_calibration_threshold_is_unresolved_not_negative(self):
        policy = {"id": "sparse-global", "kind": "global", "threshold": None,
                  "passthrough_below": None}
        result = thresholds.policy_decision(policy, -100.0, 10.0, "original")
        self.assertTrue(result["retained"])
        self.assertTrue(result["deferred"])

    def test_depth_bins_have_half_open_boundaries(self):
        policy = {"id": "table", "kind": "depth_table", "passthrough_below": None,
                  "bins": [
                      {"lower": 0, "upper": 5, "threshold": None, "reason": "sparse"},
                      {"lower": 5, "upper": 10, "threshold": 1.0, "reason": "estimated"},
                      {"lower": 10, "upper": None, "threshold": -1.0, "reason": "estimated"},
                  ]}
        low = thresholds.policy_decision(policy, 0.0, 4.0, "original")
        middle = thresholds.policy_decision(policy, 0.0, 5.0, "original")
        upper = thresholds.policy_decision(policy, 0.0, 10.0, "original")
        self.assertEqual((low["retained"], low["deferred"]), (True, True))
        self.assertEqual((middle["retained"], middle["deferred"]), (False, False))
        self.assertEqual((upper["retained"], upper["deferred"]), (True, False))

    def test_new_global_policy_is_condition_blind(self):
        policy = {"id": "shared", "kind": "global", "threshold": 0.1,
                  "passthrough_below": 5}
        original = thresholds.policy_decision(policy, 0.2, 8, "original")
        reduced = thresholds.policy_decision(policy, 0.2, 8, "downsampled")
        novel = thresholds.policy_decision(policy, 0.2, 8, "untrained-depth")
        self.assertEqual(original, reduced)
        self.assertEqual(original, novel)

    def test_legacy_condition_requires_explicit_mapping(self):
        policy = {"id": "old", "kind": "legacy", "passthrough_below": None,
                  "legacy_thresholds": {"original": 0.5, "downsampled": -0.5}}
        self.assertFalse(thresholds.policy_decision(policy, 0, 10, "original")["retained"])
        self.assertTrue(thresholds.policy_decision(policy, 0, 10, "downsampled")["retained"])
        with self.assertRaises((ValueError, KeyError)):
            thresholds.policy_decision(policy, 0, 10, "unmapped-depth")

    def test_nonfinite_observations_fail_instead_of_becoming_rejections(self):
        policy = {"id": "valid", "kind": "global", "threshold": 0.0,
                  "passthrough_below": 5}
        for score, depth in ((math.nan, 10), (math.inf, 10), (0, math.nan),
                             (0, math.inf), (0, -1)):
            with self.subTest(score=score, depth=depth):
                with self.assertRaises(ValueError):
                    thresholds.policy_decision(policy, score, depth, "original")

    def test_uncovered_depth_bin_is_not_silently_discarded(self):
        policy = {"id": "gap", "kind": "depth_table", "passthrough_below": None,
                  "bins": [{"lower": 5, "upper": 10, "threshold": 0.0, "reason": "estimated"}]}
        with self.assertRaises(ValueError):
            thresholds.policy_decision(policy, 1.0, 15.0, "original")


class SpatialSupportTests(unittest.TestCase):
    def bounds(self, rows, decisions, *, structural=False, protocol=None, stage="development"):
        return thresholds.simultaneous_recall_bounds(
            rows, {"policy": {"decisions": decisions, "structural": structural}},
            synthetic_protocol() if protocol is None else protocol, stage=stage)

    def test_nonstructural_zero_errors_do_not_certify_perfect_recall(self):
        rows = spatial_rows()
        result = self.bounds(rows, [True] * len(rows))
        for condition in ("original", "downsampled"):
            value = support(result, condition=condition)
            self.assertEqual(value["recall"], 1.0)
            self.assertEqual(value["error_block_count"], 0)
            self.assertFalse(value["supported"])
            self.assertFalse(value["structural"])

    def test_structural_retention_is_distinct_from_bootstrap_certainty(self):
        rows = spatial_rows()
        result = self.bounds(rows, [True] * len(rows), structural=True)
        for condition in ("original", "downsampled"):
            value = support(result, condition=condition)
            self.assertTrue(value["structural"])
            self.assertTrue(value["supported"])
            self.assertEqual(value["recall"], 1.0)
            self.assertEqual(value["lower_bound"], 1.0)

    def test_structural_flag_cannot_hide_a_rejected_candidate(self):
        rows = spatial_rows()
        decisions = [True] * len(rows)
        decisions[0] = False
        with self.assertRaises(ValueError):
            self.bounds(rows, decisions, structural=True)

    def test_errors_confined_to_one_window_are_insufficient_support(self):
        rows = spatial_rows()
        decisions = [not (r["pos0"] < 100_000 and r["label"] == 1) for r in rows]
        result = self.bounds(rows, decisions)
        for condition in ("original", "downsampled"):
            self.assertEqual(support(result, condition=condition)["error_block_count"], 1)
            self.assertFalse(support(result, condition=condition)["supported"])

    def test_two_error_windows_can_support_an_empirical_bound(self):
        rows = spatial_rows()
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        result = self.bounds(rows, decisions)
        for condition in ("original", "downsampled"):
            value = support(result, condition=condition)
            self.assertEqual(value["error_block_count"], 2)
            self.assertTrue(value["supported"])
            self.assertAlmostEqual(value["recall"], 118 / 120)
            self.assertLessEqual(value["lower_bound"], value["recall"])

    def test_positive_counts_do_not_pool_correlated_depth_conditions(self):
        rows = spatial_rows(blocks=6, positives=20)
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        result = self.bounds(rows, decisions)
        self.assertEqual(support(result, condition="original")["positive_count"], 120)
        self.assertEqual(support(result, condition="downsampled")["positive_count"], 120)

    def test_unknown_truth_does_not_become_negative_or_enter_recall_support(self):
        rows = spatial_rows()
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        before = self.bounds(rows, decisions)
        unknown = dict(rows[0], candidate_id="unassessed", pos0=50_000, label=None)
        after = self.bounds(rows + [unknown], decisions + [False])
        for condition in ("original", "downsampled"):
            old = support(before, condition=condition)
            new = support(after, condition=condition)
            for key in ("positive_count", "error_block_count", "recall", "lower_bound"):
                self.assertEqual(old[key], new[key])

    def test_lower_bounds_are_identical_for_identical_paired_conditions(self):
        rows = spatial_rows()
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        result = self.bounds(rows, decisions)
        original = support(result, condition="original")
        reduced = support(result, condition="downsampled")
        self.assertEqual(original["recall"], reduced["recall"])
        self.assertEqual(original["lower_bound"], reduced["lower_bound"])

    def test_condition_specific_failure_cannot_be_hidden_by_pooling(self):
        rows = spatial_rows()
        decisions = [not (r["condition"] == "downsampled" and r["label"] == 1)
                     for r in rows]
        result = self.bounds(rows, decisions)
        self.assertEqual(support(result, condition="original")["recall"], 1.0)
        self.assertEqual(support(result, condition="downsampled")["recall"], 0.0)
        self.assertLess(support(result, condition="downsampled")["lower_bound"], 0.99)

    def test_minimum_nonempty_windows_is_enforced(self):
        rows = spatial_rows(blocks=2)
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        result = self.bounds(rows, decisions)
        self.assertFalse(support(result)["supported"])

    def test_deterministic_resampling_is_reproducible(self):
        rows = spatial_rows()
        decisions = [not (r["pos0"] in (0, 100_000)) for r in rows]
        first = self.bounds(rows, decisions)
        second = self.bounds(rows, decisions)
        self.assertEqual(first, second)

    def test_length_mismatch_cannot_drop_or_duplicate_predictions(self):
        rows = spatial_rows()
        with self.assertRaises(ValueError):
            self.bounds(rows, [True] * (len(rows) - 1))

    def test_extra_search_candidates_cannot_loosen_an_existing_bound(self):
        rows = spatial_rows()
        original = [not (r["pos0"] in (0, 100_000)) for r in rows]
        competitor = [not (r["pos0"] // 100_000 in (0, 1, 2) and r["label"] == 1)
                      for r in rows]
        one = self.bounds(rows, original)
        two = thresholds.simultaneous_recall_bounds(rows, {
            "policy": {"decisions": original, "structural": False},
            "additional-search-option": {"decisions": competitor, "structural": False},
        }, synthetic_protocol(), stage="development")
        self.assertLessEqual(support(two)["lower_bound"], support(one)["lower_bound"])


class FreezeAndEvaluationTests(unittest.TestCase):
    def test_freeze_never_scores_test_features_or_uses_test_labels(self):
        records, bundle, protocol = policy_fixture()
        guarded = copy.deepcopy(records)
        flipped = copy.deepcopy(records)
        for row in guarded:
            if row["partition"] == "test":
                row["f0"]["score"] = row["f1"]["score"] = 1_000_000.0
                row["f1"]["rich_score"] = 1_000_000.0
        for row in flipped:
            if row["partition"] == "test" and row["label"] is not None:
                row["label"] = 1 - row["label"]
        with tempfile.TemporaryDirectory() as tmp, patch.object(
                thresholds.study, "_new_model", side_effect=AssertionError("Base-model creation is forbidden")):
            root = Path(tmp)
            first = thresholds.freeze_policies(records, bundle, root / "first", protocol)
            hidden = thresholds.freeze_policies(guarded, bundle, root / "hidden", protocol)
            changed = thresholds.freeze_policies(flipped, bundle, root / "changed", protocol)
        self.assertEqual(first["selected"], hidden["selected"])
        self.assertEqual(first["policies"], hidden["policies"])
        self.assertEqual(first["policy_sha256"], hidden["policy_sha256"])
        self.assertEqual(first["policy_sha256"], changed["policy_sha256"])
        self.assertEqual(first["base_model_refits"], 0)
        self.assertFalse(first["test_labels_used_for_fitting_or_selection"])
        self.assertTrue(all(first["policies"][p]["eligible_for_selection"] for p in first["selected"].values()))

    def test_fresh_labels_change_metrics_but_never_frozen_actions(self):
        records, bundle, protocol = policy_fixture()
        fresh = fresh_fixture(records)
        changed = copy.deepcopy(fresh)
        for row in changed:
            if row["label"] is not None:
                row["label"] = 1 - row["label"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen = thresholds.freeze_policies(records, bundle, root / "freeze", protocol)
            original_freeze = copy.deepcopy(frozen)
            first = thresholds.evaluate_frozen(fresh, bundle, frozen, root / "first", protocol)
            second = thresholds.evaluate_frozen(changed, bundle, frozen, root / "changed", protocol)
            def predictions(directory):
                rows = [json.loads(line) for line in (directory / "predictions.jsonl").read_text().splitlines()]
                return [{k: v for k, v in row.items() if k not in ("label", "label_reason")} for row in rows]
            self.assertEqual(predictions(root / "first"), predictions(root / "changed"))
        self.assertEqual(frozen, original_freeze)
        self.assertEqual(first["selected"], second["selected"])
        self.assertFalse(first["fresh_label_threshold_fitting"])
        self.assertEqual(first["base_model_refits"], 0)
        self.assertFalse(first["gate"]["selective_escalation_authorized"])
        self.assertFalse(first["gate"]["proceed_to_p1"])

    def test_modified_policy_or_protocol_is_rejected_before_scoring(self):
        records, bundle, protocol = policy_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen = thresholds.freeze_policies(records, bundle, root / "freeze", protocol)
            changed = copy.deepcopy(frozen)
            chosen = changed["selected"]["counts"]
            changed["policies"][chosen]["passthrough_below"] = 50
            with self.assertRaisesRegex(ValueError, "policy"):
                thresholds.evaluate_frozen(fresh_fixture(records), bundle, changed, root / "bad-policy", protocol)
            with self.assertRaisesRegex(ValueError, "protocol|model"):
                thresholds.evaluate_frozen(fresh_fixture(records), bundle, frozen, root / "bad-protocol",
                                           dict(protocol, target_recall=0.98))

    def test_standalone_counts_prediction_does_not_require_rich_fields_or_labels(self):
        records, bundle, _ = policy_fixture()
        rows = [{"candidate_id": r["candidate_id"], "condition": r["condition"], "f0": r["f0"]}
                for r in records if r["partition"] == "test"]
        policy = {"id": "counts-global", "kind": "global", "model_id": "M10_pooled",
                  "threshold": 0.0, "passthrough_below": 5}
        manifest = {"selected": {"counts": "counts-global"}, "policies": {"counts-global": policy}}
        result = thresholds.predict_policy(rows, bundle, manifest, "counts")
        self.assertEqual(len(result), len(rows))
        self.assertTrue(any(a["retained"] for a in result))
        self.assertTrue(any(not a["retained"] for a in result))
        self.assertTrue(all(not a["deferred"] for a in result))

    def test_parquet_cost_path_preserves_standalone_predictions(self):
        import costs
        records, bundle, protocol = policy_fixture()
        fresh = fresh_fixture(records)
        policies = {
            label: {"id": label, "kind": "global", "model_id": model,
                    "threshold": 0.0, "passthrough_below": 5}
            for label, model in (("counts", "M10_pooled"), ("rich", "M11_pooled"))
        }
        manifest = {"selected": {"counts": "counts", "rich": "rich"}, "policies": policies}
        with tempfile.TemporaryDirectory() as tmp:
            result = costs.measure_costs(fresh, bundle, manifest, Path(tmp) / "costs", protocol)
        self.assertTrue(result["native_parquet_prediction_parity"])
        self.assertFalse(result["net_end_user_benefit_established"])
        self.assertIsNone(result["manual_review_seconds"])
        self.assertEqual(result["storage"]["counts"]["candidate_rows_including_unknown"],
                         result["storage"]["rich"]["candidate_rows_including_unknown"])


def passing_gate_fixture():
    _, _, protocol = policy_fixture()
    protocol["gates"].update({
        "minimum_relative_fp_reduction": 0.25,
        "minimum_absolute_fp_reduction": 10,
        "minimum_assessed_retained_reduction": 0.01,
        "minimum_total_retained_reduction": 0.01,
    })
    selected = {"counts": "counts-selected", "rich": "rich-selected"}
    conditions = protocol["required_fresh_conditions"]
    comparisons = {condition: {
        "absolute_fp_reduction": 30, "relative_fp_reduction": 0.3,
        "recall_rich_minus_counts": 0.0,
        "assessed_workload_reduction": 0.02, "all_workload_reduction": 0.02,
        "bootstrap": {
            "block_count": 6, "requested_replicates": 80,
            "valid_replicates": 80, "invalid_replicates": 0,
            "absolute_fp_reduction_95ci": {"lower": 10, "upper": 50},
            "recall_rich_minus_counts_95ci": {"lower": -0.001, "upper": 0.001},
        },
    } for condition in conditions}
    bounds = {"policies": {policy: {"conditions": {
        condition: {"supported": True, "lower_bound": 0.995}
        for condition in conditions}}
        for policy in selected.values()}}
    return comparisons, bounds, selected, protocol


class AdvancementGateTests(unittest.TestCase):
    def gate(self, comparisons=None, bounds=None, selected=None, protocol=None):
        c, b, s, p = passing_gate_fixture()
        return thresholds.evaluate_gate(c if comparisons is None else comparisons,
                                        b if bounds is None else bounds,
                                        s if selected is None else selected,
                                        p if protocol is None else protocol)

    def test_passing_predictive_screen_never_authorizes_a_router_or_p1(self):
        result = self.gate()
        self.assertTrue(result["predictive_and_workload_screen_passed"])
        self.assertFalse(result["selective_escalation_authorized"])
        self.assertFalse(result["proceed_to_p1"])
        self.assertFalse(result["total_workflow_cost_advantage_established"])

    def test_matching_poor_recall_cannot_pass_noninferiority_alone(self):
        _, bounds, selected, _ = passing_gate_fixture()
        for policy in selected.values():
            bounds["policies"][policy]["conditions"]["novel_low"]["lower_bound"] = 0.98
        self.assertFalse(self.gate(bounds=bounds)["predictive_and_workload_screen_passed"])

    def test_an_unsupported_absolute_bound_is_not_a_pass(self):
        _, bounds, selected, _ = passing_gate_fixture()
        entry = bounds["policies"][selected["rich"]]["conditions"]["original"]
        entry.update(supported=False, lower_bound=1.0)
        self.assertFalse(self.gate(bounds=bounds)["predictive_and_workload_screen_passed"])

    def test_inadequate_paired_comparison_support_prevents_advancement(self):
        for field, value in (("block_count", 2), ("valid_replicates", 71)):
            with self.subTest(field=field):
                comparisons, _, _, _ = passing_gate_fixture()
                comparisons["downsampled"]["bootstrap"][field] = value
                self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])

    def test_missing_novel_condition_prevents_advancement(self):
        comparisons, _, _, _ = passing_gate_fixture()
        del comparisons["novel_low"]
        self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])

    def test_recall_difference_bound_must_strictly_exceed_margin(self):
        comparisons, _, _, _ = passing_gate_fixture()
        comparisons["downsampled"]["bootstrap"]["recall_rich_minus_counts_95ci"]["lower"] = -0.005
        self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])

    def test_unknown_workload_cannot_hide_behind_assessed_only_savings(self):
        comparisons, _, _, _ = passing_gate_fixture()
        comparisons["downsampled"]["all_workload_reduction"] = -0.02
        self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])

    def test_each_workload_threshold_uses_its_declared_configuration_key(self):
        for key in ("minimum_assessed_retained_reduction", "minimum_total_retained_reduction"):
            with self.subTest(key=key):
                _, _, _, protocol = passing_gate_fixture()
                protocol["gates"][key] = 0.03
                self.assertFalse(self.gate(protocol=protocol)["predictive_and_workload_screen_passed"])

    def test_required_condition_configuration_cannot_be_silently_ignored(self):
        _, _, _, protocol = passing_gate_fixture()
        protocol["required_fresh_conditions"] += ["missing_required_condition"]
        self.assertFalse(self.gate(protocol=protocol)["predictive_and_workload_screen_passed"])

    def test_fp_reduction_requires_magnitude_and_uncertainty_support(self):
        for field, value in (("absolute_fp_reduction", 9), ("relative_fp_reduction", 0.24)):
            with self.subTest(field=field):
                comparisons, _, _, _ = passing_gate_fixture()
                comparisons["downsampled"][field] = value
                self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])
        comparisons, _, _, _ = passing_gate_fixture()
        comparisons["downsampled"]["bootstrap"]["absolute_fp_reduction_95ci"]["lower"] = 0
        self.assertFalse(self.gate(comparisons=comparisons)["predictive_and_workload_screen_passed"])


if __name__ == "__main__":
    unittest.main()
