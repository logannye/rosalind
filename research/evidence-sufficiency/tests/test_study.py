"""Small expected-value and fit-boundary tests; no native executable required."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import study


def fixture_records():
    records = []
    roles = [("train", "train-a"), ("train", "train-b"),
             ("calibration", "cal-a"), ("threshold", "threshold-a"), ("test", "test-a")]
    for block_index, (partition, block) in enumerate(roles):
        for index in range(40):
            label = index % 2
            identity = f"{block}:{index}"
            for condition in ("original", "downsampled"):
                depth = 30 if condition == "original" else 10
                f0 = {"callable_depth": depth, "alt_fraction": .1,
                      "constant": 1.0, "undefined": None}
                f1 = dict(f0, quality_difference=float(label), strand_fraction=.5)
                records.append({"candidate_id": identity, "contig": "chr20",
                                "pos0": block_index * 100000 + index, "ref": "A", "alt": "C",
                                "partition": partition, "block_id": block, "condition": condition,
                                "label": label, "label_reason": "fixture", "f0": f0, "f1": f1})
    return records


class StudyPureTests(unittest.TestCase):
    def test_protocol_cannot_silently_change_frozen_parameters(self):
        valid = {"gates": {"min_development_negatives": 50, "min_mean_ap_gain": .005,
                           "target_recall": .99}, "model_seed": 42}
        study._validate_protocol(valid)
        for key, value in (("min_development_negatives", 49), ("min_mean_ap_gain", .001),
                           ("target_recall", .95)):
            changed = copy.deepcopy(valid)
            changed["gates"][key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "frozen P0"):
                study.run_study([], Path("must-not-be-created"), changed)
        with self.assertRaisesRegex(ValueError, "model_seed"):
            study._validate_protocol({"model_seed": 43})

    def test_caller_predictions_exclude_unprovided_condition(self):
        records = fixture_records()[:2]
        records[0]["caller_score"] = None  # Original: explicitly not emitted.
        bundle = {"feature_names": study._validate_records(records), "arms": {
            "B1": {"kind": "caller", "available_conditions": ["original"],
                   "operating_point": {"eligible": True, "threshold": .5}}}}
        with patch.object(study, "_dependencies", return_value={}), \
                patch.object(study, "_calibrated", return_value=[.2]) as calibrated:
            predictions = study.predict_records(bundle, records)
        self.assertEqual(len(predictions), 1)
        self.assertEqual(predictions[0]["condition"], "original")
        self.assertEqual(calibrated.call_args.args[1], [records[0]])
        records[0].pop("caller_score")
        with patch.object(study, "_dependencies", return_value={}), \
                self.assertRaisesRegex(ValueError, "explicit caller_score"):
            study.predict_records(bundle, records)

    def test_average_precision_respects_ties(self):
        self.assertEqual(study.average_precision([1, 0], [.5, .5]), .5)
        self.assertEqual(study.average_precision([1, 0], [.9, .1]), 1)
        self.assertAlmostEqual(study.average_precision([0, 1, 1], [.9, .8, .7]), (0.5 + 2 / 3) / 2)
        self.assertIsNone(study.average_precision([0, 0], [.9, .1]))

    def test_highest_threshold_keeps_ties(self):
        result = study.select_score_threshold([1, 1, 0, 0], [.9, .5, .5, .1])
        self.assertEqual(result["threshold"], .5)
        self.assertEqual(result["false_positives"], 1)
        self.assertEqual(result["recall"], 1)
        self.assertFalse(study.select_score_threshold([0], [.9])["eligible"])

    def test_depth_grid_and_ineligible_zero_depth(self):
        rows = [{"label": y, "f0": {"callable_depth": d, "alt_fraction": f}}
                for y, d, f in [(1, 20, .2), (1, 10, .2), (0, 20, .1), (0, 4, .3)]]
        selected = study.select_b0(rows)
        self.assertEqual((selected["depth_threshold"], selected["fraction_threshold"]), (5, .2))
        self.assertEqual(selected["false_positives"], 0)
        self.assertFalse(study.select_b0([{"label": 1, "f0": {"callable_depth": 0, "alt_fraction": None}}])["eligible"])

    def test_development_gate_never_uses_test_or_depth_derivatives(self):
        records = fixture_records()
        original = study.development_adequacy(records)
        self.assertTrue(original["adequate"])
        for row in records:
            if row["partition"] == "test" or row["condition"] == "downsampled":
                row["label"] = None
        self.assertEqual(study.development_adequacy(records), original)

    def test_one_class_fold_fails_adequacy_even_with_large_pooled_counts(self):
        records = fixture_records()
        for row in records:
            if row["block_id"] == "train-a":
                row["label"] = 0
        self.assertFalse(study.development_adequacy(records)["adequate"])

    def test_refuse_split_leakage_and_f0_f1_disagreement(self):
        records = fixture_records()
        study._validate_records(records)
        records[1]["partition"] = "test"
        with self.assertRaisesRegex(ValueError, "crosses partition"):
            study._validate_records(records)
        records = fixture_records()
        records[0]["f1"]["callable_depth"] = 999
        with self.assertRaisesRegex(ValueError, "F0 values"):
            study._validate_records(records)

    def test_unknowns_excluded_from_metrics_without_removing_candidates(self):
        predictions = [{"partition": "test", "condition": "original", "arm": "M00",
                        "label": y, "accepted": a, "probability": p}
                       for y, a, p in [(1, True, .9), (0, True, .8), (None, False, .1)]]
        metric = study.evaluate_predictions(predictions)["test"]["original"]["M00"]
        self.assertEqual((metric["candidate_count"], metric["evaluated"], metric["unknown_labels"]), (3, 2, 1))
        self.assertEqual(metric["false_positives"], 1)
        self.assertAlmostEqual(metric["brier_score"], .325)

    def test_insufficient_data_skips_dependency_import_and_fits(self):
        records = fixture_records()[:4]
        with tempfile.TemporaryDirectory() as directory:
            report = study.run_study(records, Path(directory), {})
            self.assertEqual(report["gate"]["status"], "insufficient_data")
            self.assertEqual(report["costs"]["model_fits"], 0)
            self.assertFalse((Path(directory) / "frozen-models.pkl").exists())
            self.assertEqual((Path(directory) / "predictions.jsonl").read_text(), "")
            with self.assertRaises(FileExistsError):
                study.run_study(records, Path(directory), {})


@unittest.skipUnless(importlib.util.find_spec("sklearn") and importlib.util.find_spec("scipy"),
                     "sklearn/SciPy wheels required; run estimator checks on hosted validation")
class StudyEstimatorTests(unittest.TestCase):
    def test_fit_boundaries_replay_and_test_independent_gate(self):
        records = fixture_records()
        # Test covariates differ dramatically; the fitted imputer must not see them.
        for row in records:
            if row["partition"] == "test":
                row["f0"]["constant"] = row["f1"]["constant"] = 10000
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = study.run_study(records, root / "one", {})
            self.assertEqual(report["gate"]["status"], "proceed_to_p1")
            bundle = study.load_bundle(root / "one")
            names = bundle["feature_names"]["f0"]
            imputer = bundle["arms"]["M00"]["model"].named_steps["imputer"]
            self.assertEqual(imputer.statistics_[names.index("constant")], 1)
            original_predictions = [json.loads(line) for line in (root / "one" / "predictions.jsonl").read_text().splitlines()]
            self.assertEqual(study.predict_records(root / "one", records), original_predictions)
            changed = copy.deepcopy(records)
            for row in changed:
                if row["partition"] == "test":
                    row["label"] = 1 - row["label"]
            second = study.run_study(changed, root / "two", {})
            self.assertEqual(report["gate"], second["gate"])
            self.assertEqual(report["frozen_model"], second["frozen_model"])
            self.assertNotEqual(report["metrics"]["test"], second["metrics"]["test"])

    def test_model_corruption_refused_before_unpickling(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            study.run_study(fixture_records(), path, {})
            model = path / "frozen-models.pkl"
            model.write_bytes(model.read_bytes() + b"mutation")
            with self.assertRaisesRegex(ValueError, "identity"):
                study.load_bundle(path)


if __name__ == "__main__":
    unittest.main()
