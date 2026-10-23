"""Frozen P0 models and development-only go/stop decision.

This module consumes prepared candidate feature dictionaries. It does not read
alignments, generate labels, or establish native/Parquet equivalence. Native
correctness and backend equivalence remain independent prerequisites. The only
serialized executable objects are locally created sklearn models: never load an
untrusted pickle with ``load_bundle``.
"""
from __future__ import annotations

import hashlib
import json
import math
import pickle
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ARMS = (("M00", "f0", "LR"), ("M01", "f1", "LR"),
        ("M10", "f0", "HGB"), ("M11", "f1", "HGB"))
PARTITIONS = ("train", "calibration", "threshold", "test")
DEPTH_GRID = (0, 5, 10, 20)
FRACTION_GRID = (0.02, 0.05, 0.10, 0.20, 0.30)
FORBIDDEN_FEATURES = {
    "candidate_id", "contig", "chromosome", "chrom", "pos", "pos0", "position",
    "sample", "sample_id", "specimen_id", "subject", "family", "partition",
    "block_id", "condition", "label", "label_reason", "benchmark_membership",
    "caller_score",
}
OUTPUTS = ("frozen-models.pkl", "frozen-models.json", "predictions.jsonl",
           "metrics.json", "gate.json", "report.json")


def _validate_protocol(protocol: dict) -> None:
    """P0 parameters are frozen; a protocol edit must not silently be ignored."""
    gates = protocol.get("gates", {})
    if not isinstance(gates, dict):
        raise ValueError("protocol gates must be an object")
    for key, expected in (("min_development_negatives", 50), ("min_mean_ap_gain", .005),
                          ("target_recall", .99)):
        value = gates.get(key, expected)
        if type(value) not in (int, float) or value != expected:
            raise ValueError(f"frozen P0 gate {key} must equal {expected}")
    seed = protocol.get("model_seed", 42)
    if type(seed) is not int or seed != 42:
        raise ValueError("frozen P0 model_seed must equal 42")
    if "target_recall" in protocol:
        value = protocol["target_recall"]
        if type(value) not in (int, float) or value != .99:
            raise ValueError("frozen P0 target_recall must equal 0.99")


def _canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _write_json(path: Path, value: Any) -> None:
    with path.open("xb") as handle:
        handle.write(_canonical(value))


def _labelled(records: list[dict]) -> list[dict]:
    return [row for row in records if row["label"] is not None]


def class_counts(records: list[dict]) -> dict:
    counts = Counter(row["label"] for row in records)
    return {"negative": counts[0], "positive": counts[1], "unknown": counts[None],
            "total": len(records)}


def _validate_records(records: list[dict]) -> dict[str, list[str]]:
    """Refuse duplicated observations and locus/block role leakage, including derivatives."""
    seen, loci, blocks, names = set(), {}, {}, {}
    candidate_identity = {}
    for row in records:
        for key in ("candidate_id", "contig", "ref", "alt", "block_id", "condition"):
            if not isinstance(row.get(key), str) or not row[key]:
                raise ValueError(f"record requires nonempty {key}")
        if type(row.get("pos0")) is not int or row["pos0"] < 0:
            raise ValueError("pos0 must be a nonnegative integer")
        if row.get("partition") not in PARTITIONS:
            raise ValueError("unknown partition")
        if row.get("label") is not None and (type(row["label"]) is not int or row["label"] not in (0, 1)):
            raise ValueError("label must be 0, 1, or None")
        key = (row["candidate_id"], row["condition"])
        if key in seen:
            raise ValueError("duplicate candidate/condition observation")
        seen.add(key)
        locus = (row["contig"], row["pos0"])
        role = row["partition"]
        if loci.setdefault(locus, role) != role or blocks.setdefault(row["block_id"], role) != role:
            raise ValueError("locus or genomic block crosses partition boundaries")
        identity = (row["contig"], row["pos0"], row["ref"], row["alt"], role, row["block_id"])
        if candidate_identity.setdefault(row["candidate_id"], identity) != identity:
            raise ValueError("candidate identity differs between conditions")
        for group in ("f0", "f1"):
            features = row.get(group)
            if not isinstance(features, dict) or not features:
                raise ValueError(f"record requires nonempty {group} features")
            keys = sorted(features)
            if any(not isinstance(name, str) or name in FORBIDDEN_FEATURES for name in keys):
                raise ValueError("metadata or forbidden field in model features")
            if group in names and names[group] != keys:
                raise ValueError(f"inconsistent feature schema for {group}")
            names[group] = keys
            for value in features.values():
                if value is not None and (type(value) not in (int, float)
                                          or not math.isfinite(value)):
                    raise ValueError("features must be finite numeric values or None")
        if not set(row["f0"]).issubset(row["f1"]):
            raise ValueError("F1 must include every F0 feature")
        if any(row["f1"][key] != value for key, value in row["f0"].items()):
            raise ValueError("F0 values changed in F1")
    return names


def development_adequacy(records: list[dict], original_condition: str = "original") -> dict:
    """Only original-condition development labels can affect this decision."""
    development = [r for r in records if r["condition"] == original_condition
                   and r["partition"] != "test"]
    counts = {role: class_counts([r for r in development if r["partition"] == role])
              for role in PARTITIONS[:-1]}
    train_blocks = sorted({r["block_id"] for r in development if r["partition"] == "train"})
    block_counts = {block: class_counts([r for r in development
                                       if r["partition"] == "train" and r["block_id"] == block])
                    for block in train_blocks}
    reasons = []
    if len(train_blocks) != 2:
        reasons.append("P0 requires exactly two training blocks for two-fold development CV")
    for name, count in list(counts.items()) + list(block_counts.items()):
        if not count["negative"] or not count["positive"]:
            reasons.append(f"both labelled classes required in development role/block: {name}")
    negative_count = sum(count["negative"] for count in counts.values())
    if negative_count < 50:
        reasons.append("fewer than 50 labelled negative candidates in original-condition development")
    return {"adequate": not reasons, "reasons": reasons, "partition_counts": counts,
            "training_block_counts": block_counts, "development_negative_count": negative_count,
            "original_condition": original_condition,
            "test_outcomes_used": False}


def average_precision(labels: list[int], scores: list[float]) -> float | None:
    """Tie-aware AP, calculated at distinct descending score levels."""
    if len(labels) != len(scores):
        raise ValueError("label/score length mismatch")
    positives = sum(labels)
    if not positives:
        return None
    groups = defaultdict(list)
    for label, score in zip(labels, scores):
        if label not in (0, 1) or not math.isfinite(score):
            raise ValueError("invalid metric input")
        groups[float(score)].append(label)
    true_positives = total = 0
    result = 0.0
    for score in sorted(groups, reverse=True):
        added = sum(groups[score])
        true_positives += added
        total += len(groups[score])
        result += (added / positives) * (true_positives / total)
    return result


def select_score_threshold(labels: list[int], scores: list[float], target: float = 0.99) -> dict:
    """Highest observed score threshold reaching target recall; ties remain together."""
    if len(labels) != len(scores) or not 0 < target <= 1:
        raise ValueError("invalid threshold inputs")
    if any(label not in (0, 1) for label in labels) or any(not math.isfinite(s) for s in scores):
        raise ValueError("invalid label/score")
    positive_scores = sorted((float(s) for y, s in zip(labels, scores) if y), reverse=True)
    if not positive_scores:
        return {"eligible": False, "reason": "no positive threshold-selection candidates"}
    threshold = positive_scores[max(0, math.ceil(target * len(positive_scores)) - 1)]
    decisions = [s >= threshold for s in scores]
    tp = sum(y == 1 and d for y, d in zip(labels, decisions))
    fp = sum(y == 0 and d for y, d in zip(labels, decisions))
    return {"eligible": True, "threshold": threshold, "target_recall": target,
            "recall": tp / len(positive_scores), "true_positives": tp,
            "false_positives": fp, "positive_count": len(positive_scores)}


def _b0_values(row: dict, depth_key: str, fraction_key: str) -> tuple[Any, Any]:
    if depth_key not in row["f0"] or fraction_key not in row["f0"]:
        raise ValueError("B0 depth/fraction keys missing from F0")
    return row["f0"][depth_key], row["f0"][fraction_key]


def select_b0(records: list[dict], target: float = 0.99, *,
              depth_key: str = "callable_depth", fraction_key: str = "alt_fraction") -> dict:
    records = _labelled(records)
    positives = sum(row["label"] for row in records)
    if not positives:
        return {"eligible": False, "reason": "no positive threshold-selection candidates"}
    choices = []
    for depth in DEPTH_GRID:
        for fraction in FRACTION_GRID:
            accepted = []
            for row in records:
                d, f = _b0_values(row, depth_key, fraction_key)
                accepted.append(d is not None and f is not None and d >= depth and f >= fraction)
            tp = sum(row["label"] == 1 and value for row, value in zip(records, accepted))
            fp = sum(row["label"] == 0 and value for row, value in zip(records, accepted))
            recall = tp / positives
            if recall >= target:
                choices.append((fp, -recall, depth, fraction, tp))
    if not choices:
        return {"eligible": False, "reason": "fixed B0 grid cannot reach target recall",
                "target_recall": target}
    fp, negative_recall, depth, fraction, tp = min(choices)
    return {"eligible": True, "depth_threshold": depth, "fraction_threshold": fraction,
            "recall": -negative_recall, "true_positives": tp, "false_positives": fp,
            "positive_count": positives, "target_recall": target,
            "depth_key": depth_key, "fraction_key": fraction_key}


def decision_metrics(labels: list[int], accepted: list[bool],
                     probabilities: list[float] | None = None) -> dict:
    if len(labels) != len(accepted) or (probabilities is not None and len(labels) != len(probabilities)):
        raise ValueError("metric lengths differ")
    if any(type(y) is not int or y not in (0, 1) for y in labels):
        raise ValueError("metrics require binary labels")
    tp = sum(y == 1 and d for y, d in zip(labels, accepted))
    fp = sum(y == 0 and d for y, d in zip(labels, accepted))
    fn = sum(y == 1 and not d for y, d in zip(labels, accepted))
    tn = sum(y == 0 and not d for y, d in zip(labels, accepted))
    result = {"evaluated": len(labels), "true_positives": tp, "false_positives": fp,
              "false_negatives": fn, "true_negatives": tn,
              "recall": tp / (tp + fn) if tp + fn else None,
              "precision": tp / (tp + fp) if tp + fp else None,
              "accepted_fraction": sum(accepted) / len(labels) if labels else None,
              "average_precision": None, "brier_score": None, "log_loss": None,
              "reliability_bins": None}
    if probabilities is not None and labels:
        if any(not math.isfinite(p) or not 0 <= p <= 1 for p in probabilities):
            raise ValueError("invalid probabilities")
        result["average_precision"] = average_precision(labels, probabilities)
        result["brier_score"] = sum((p - y) ** 2 for y, p in zip(labels, probabilities)) / len(labels)
        result["log_loss"] = -sum(y * math.log(max(1e-15, p)) + (1 - y) * math.log(max(1e-15, 1 - p))
                                  for y, p in zip(labels, probabilities)) / len(labels)
        bins = []
        for index in range(10):
            selected = [(y, p) for y, p in zip(labels, probabilities)
                        if min(9, int(p * 10)) == index]
            bins.append({"lower": index / 10, "upper": (index + 1) / 10,
                         "count": len(selected),
                         "mean_probability": sum(p for _, p in selected) / len(selected) if selected else None,
                         "positive_fraction": sum(y for y, _ in selected) / len(selected) if selected else None})
        result["reliability_bins"] = bins
    return result


def _dependencies() -> dict:
    import numpy as np
    import sklearn
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    return locals()


def _matrix(records: list[dict], group: str, names: list[str], np):
    return np.array([[math.nan if row[group][name] is None else float(row[group][name])
                      for name in names] for row in records], dtype=np.float64)


def _new_model(learner: str, dep: dict):
    imputer = dep["SimpleImputer"](strategy="median", add_indicator=True, keep_empty_features=True)
    if learner == "LR":
        estimator = dep["LogisticRegression"](C=1.0, max_iter=2000, solver="lbfgs", random_state=42)
        return dep["Pipeline"]([("imputer", imputer), ("scaler", dep["StandardScaler"]()),
                                ("classifier", estimator)])
    estimator = dep["HistGradientBoostingClassifier"](max_iter=200, learning_rate=0.05,
                    max_leaf_nodes=15, l2_regularization=1.0, early_stopping=False, random_state=42)
    return dep["Pipeline"]([("imputer", imputer), ("classifier", estimator)])


def _raw_scores(model, matrix, np):
    p = np.clip(model.predict_proba(matrix)[:, 1], 1e-12, 1 - 1e-12)
    return np.log(p / (1 - p))


def _calibrated(entry: dict, records: list[dict], names: dict, dep: dict) -> list[float]:
    if not records:
        return []
    np = dep["np"]
    if entry.get("kind") == "caller":
        raw = np.asarray([0.0 if row.get("caller_score") is None else float(row["caller_score"])
                          for row in records])
    else:
        matrix = _matrix(records, entry["feature_group"], names[entry["feature_group"]], np)
        raw = _raw_scores(entry["model"], matrix, np)
    return [float(value) for value in entry["calibrator"].predict_proba(raw.reshape(-1, 1))[:, 1]]


def _cross_validate(training: list[dict], names: dict, dep: dict) -> dict:
    blocks = sorted({r["block_id"] for r in training})
    if len(blocks) != 2:
        raise ValueError("P0 grouped CV requires exactly two training blocks")
    result = {arm: {"feature_group": group, "learner": learner, "folds": []}
              for arm, group, learner in ARMS}
    for held_out in blocks:
        fit = [r for r in training if r["block_id"] != held_out]
        validation = [r for r in training if r["block_id"] == held_out]
        for arm, group, learner in ARMS:
            model = _new_model(learner, dep)
            model.fit(_matrix(fit, group, names[group], dep["np"]), [r["label"] for r in fit])
            probability = [float(p) for p in model.predict_proba(
                _matrix(validation, group, names[group], dep["np"]))[:, 1]]
            metrics = decision_metrics([r["label"] for r in validation], [p >= .5 for p in probability], probability)
            result[arm]["folds"].append({"held_out_block": held_out, "fit_block_ids": sorted({r["block_id"] for r in fit}),
                                         "validation_count": len(validation),
                                         "average_precision": metrics["average_precision"], "log_loss": metrics["log_loss"]})
    for entry in result.values():
        entry["mean_average_precision"] = sum(f["average_precision"] for f in entry["folds"]) / 2
        entry["mean_log_loss"] = sum(f["log_loss"] for f in entry["folds"]) / 2
    gains = {"LR": result["M01"]["mean_average_precision"] - result["M00"]["mean_average_precision"],
             "HGB": result["M11"]["mean_average_precision"] - result["M10"]["mean_average_precision"]}
    return {"arms": result, "mean_ap_gain_f1_minus_f0": gains, "fold_count": 2,
            "selection_inputs": "original-condition training blocks only", "test_outcomes_used": False}


def predict_records(bundle_or_path, records: list[dict]) -> list[dict]:
    """Apply a frozen bundle to either backend's same candidate rows.

    Labels are copied for later evaluation but never enter prediction. Only load
    bundles this study created locally; pickle is an executable trusted format.
    """
    bundle = load_bundle(bundle_or_path) if isinstance(bundle_or_path, (str, Path)) else bundle_or_path
    names = _validate_records(records)
    if records and names != bundle["feature_names"]:
        raise ValueError("feature schema differs from frozen model")
    dep = _dependencies()
    results = []
    for arm in sorted(bundle["arms"]):
        entry = bundle["arms"][arm]
        operating = entry["operating_point"]
        arm_records = records
        if arm == "B1":
            arm_records = [r for r in records if r["condition"] in entry["available_conditions"]]
            if any("caller_score" not in r for r in arm_records):
                raise ValueError("B1 requires explicit caller_score, with None for a nonemitted candidate")
        probabilities = _calibrated(entry, arm_records, names, dep) if arm != "B0" else [None] * len(arm_records)
        for row, probability in zip(arm_records, probabilities):
            accepted = None
            if operating["eligible"]:
                if arm == "B0":
                    depth, fraction = _b0_values(row, operating["depth_key"], operating["fraction_key"])
                    accepted = (depth is not None and fraction is not None
                                and depth >= operating["depth_threshold"]
                                and fraction >= operating["fraction_threshold"])
                else:
                    accepted = probability >= operating["threshold"]
            results.append({key: row.get(key) for key in ("candidate_id", "contig", "pos0", "ref", "alt",
                            "partition", "block_id", "condition", "label", "label_reason")}
                           | {"arm": arm, "probability": probability, "accepted": accepted,
                              "operating_point_eligible": operating["eligible"]})
    return sorted(results, key=lambda row: (row["candidate_id"], row["condition"], row["arm"]))


def evaluate_predictions(predictions: list[dict]) -> dict:
    groups = defaultdict(list)
    for row in predictions:
        groups[(row["partition"], row["condition"], row["arm"])].append(row)
    result = {}
    for (partition, condition, arm), rows in sorted(groups.items()):
        known = [r for r in rows if r["label"] is not None and r["accepted"] is not None]
        probabilities = [r["probability"] for r in known]
        metrics = decision_metrics([r["label"] for r in known], [r["accepted"] for r in known],
                                   probabilities if all(p is not None for p in probabilities) else None)
        metrics.update({"candidate_count": len(rows), "unknown_labels": sum(r["label"] is None for r in rows),
                        "ineligible_predictions": sum(r["accepted"] is None for r in rows)})
        result.setdefault(partition, {}).setdefault(condition, {})[arm] = metrics
    return result


def load_bundle(path: str | Path) -> dict:
    """Load a trusted local study bundle after verifying its adjacent digest.

    Hash integrity does not make an untrusted pickle safe to execute.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "frozen-models.pkl"
    manifest = json.loads((path.parent / "frozen-models.json").read_text())
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != manifest["model_sha256"]:
        raise ValueError("frozen model bytes differ from recorded identity")
    bundle = pickle.loads(data)
    if bundle.get("schema") != "rosalind-evidence-sufficiency-models-v1":
        raise ValueError("unsupported model bundle")
    return bundle


def _finish(output_dir: Path, report: dict, predictions: list[dict]) -> dict:
    with (output_dir / "predictions.jsonl").open("xb") as handle:
        for row in predictions:
            handle.write(_canonical(row))
    for name in ("metrics", "gate"):
        _write_json(output_dir / f"{name}.json", report[name])
    _write_json(output_dir / "report.json", report)
    return report


def run_study(records: list[dict], output_dir: Path, protocol: dict) -> dict:
    """Run the fixed P0 factorial and write a development-only exploratory gate.

    This writes a frozen model bundle and manifest before scoring test rows. The
    test is descriptive P0 evidence; its outcomes never decide the P0 gate.
    """
    started, cpu_started = time.perf_counter(), time.process_time()
    _validate_protocol(protocol)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any((output_dir / name).exists() for name in OUTPUTS):
        raise FileExistsError("study output already exists; use a new run directory")
    names = _validate_records(records)
    original = protocol.get("original_condition", "original")
    target = float(protocol.get("gates", {}).get("target_recall", 0.99))
    adequacy = development_adequacy(records, original)
    report = {"schema": "rosalind-evidence-sufficiency-p0-study-v1", "protocol_sha256": _digest(protocol),
              "adequacy": adequacy, "metrics": {}, "costs": {}, "frozen_model": None,
              "scope": "One-person P0 feasibility; no P1 predictive or population claim",
              "test_outcomes_used_for_gate": False}
    if not adequacy["adequate"]:
        report["gate"] = {"status": "insufficient_data", "proceed_to_p1": False,
                          "reasons": adequacy["reasons"], "test_outcomes_used": False,
                          "required_external_gates": ["native_correctness", "parquet_equivalence", "saved_evidence_replay"]}
        report["costs"] = {"elapsed_seconds": time.perf_counter() - started,
                           "cpu_seconds": time.process_time() - cpu_started, "model_fits": 0}
        return _finish(output_dir, report, [])

    dep = _dependencies()
    np = dep["np"]
    development = {role: _labelled([r for r in records if r["condition"] == original and r["partition"] == role])
                   for role in PARTITIONS[:-1]}
    cv_started = time.perf_counter()
    cv = _cross_validate(development["train"], names, dep)
    cv_seconds = time.perf_counter() - cv_started
    arms = {}
    fit_costs = {}
    for arm, group, learner in ARMS:
        fit_started = time.perf_counter()
        model = _new_model(learner, dep)
        model.fit(_matrix(development["train"], group, names[group], np),
                  [r["label"] for r in development["train"]])
        calibration_scores = _raw_scores(model, _matrix(development["calibration"], group, names[group], np), np)
        calibrator = dep["LogisticRegression"](C=1e6, max_iter=2000, solver="lbfgs", random_state=42)
        calibrator.fit(calibration_scores.reshape(-1, 1), [r["label"] for r in development["calibration"]])
        entry = {"feature_group": group, "learner": learner, "model": model, "calibrator": calibrator}
        probability = _calibrated(entry, development["threshold"], names, dep)
        entry["operating_point"] = select_score_threshold([r["label"] for r in development["threshold"]], probability, target)
        arms[arm] = entry
        fit_costs[arm] = time.perf_counter() - fit_started
    b0 = protocol.get("b0", {})
    arms["B0"] = {"kind": "threshold_grid", "operating_point": select_b0(development["threshold"], target,
                        depth_key=b0.get("depth_key", "callable_depth"), fraction_key=b0.get("fraction_key", "alt_fraction"))}
    caller_status = {"status": "not_supplied", "reason": "no caller_score supplied; B1 is not silently substituted"}
    if any("caller_score" in r for role in development.values() for r in role):
        for row in records:
            if row["condition"] != original:
                continue
            if "caller_score" not in row:
                raise ValueError("B1 requires explicit caller_score on every original-condition candidate")
            value = row.get("caller_score")
            if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
                raise ValueError("caller_score must be nonnegative finite QUAL or None for nonemitted candidates")
        raw = np.asarray([0.0 if r.get("caller_score") is None else float(r["caller_score"])
                          for r in development["calibration"]])
        calibrator = dep["LogisticRegression"](C=1e6, max_iter=2000, solver="lbfgs", random_state=42)
        calibrator.fit(raw.reshape(-1, 1), [r["label"] for r in development["calibration"]])
        entry = {"kind": "caller", "calibrator": calibrator, "available_conditions": [original]}
        entry["operating_point"] = select_score_threshold([r["label"] for r in development["threshold"]],
                                                _calibrated(entry, development["threshold"], names, dep), target)
        arms["B1"] = entry
        caller_status = {"status": "included", "missing_score": 0.0,
                         "missing_policy": "explicit None retained with lowest nonnegative QUAL score",
                         "available_conditions": [original],
                         "other_conditions": "not evaluated; caller was not run for these inputs"}

    rank = lambda arm: (-cv["arms"][arm]["mean_average_precision"], cv["arms"][arm]["mean_log_loss"], arm)
    rich_challenger = min(("M01", "M11"), key=rank)
    probability_comparator = min(("M00", "M10"), key=rank)
    simple_order = ("B0", "B1", "M00", "M10")
    eligible = [arm for arm in simple_order if arm in arms and arms[arm]["operating_point"]["eligible"]]
    strongest_simple = min(eligible, key=lambda arm: (arms[arm]["operating_point"]["false_positives"],
                            -arms[arm]["operating_point"]["recall"], simple_order.index(arm))) if eligible else None
    bundle = {"schema": "rosalind-evidence-sufficiency-models-v1", "arms": arms,
              "feature_names": names, "original_condition": original, "protocol_sha256": _digest(protocol)}
    model_bytes = pickle.dumps(bundle, protocol=5)
    frozen = {"schema": "rosalind-evidence-sufficiency-model-freeze-v1",
              "model_sha256": hashlib.sha256(model_bytes).hexdigest(), "protocol_sha256": _digest(protocol),
              "feature_names": names, "sklearn_version": dep["sklearn"].__version__,
              "model_settings": {
                  "LR": {"C": 1.0, "max_iter": 2000, "solver": "lbfgs", "random_state": 42},
                  "HGB": {"max_iter": 200, "learning_rate": .05, "max_leaf_nodes": 15,
                          "l2_regularization": 1.0, "early_stopping": False, "random_state": 42},
                  "imputer": {"strategy": "median", "add_indicator": True, "keep_empty_features": True},
                  "scaler": "StandardScaler for LR only; fitted on training rows",
                  "sigmoid_calibrator": {"C": 1e6, "max_iter": 2000, "solver": "lbfgs", "random_state": 42},
                  "B0": {"depth_grid": list(DEPTH_GRID), "alt_fraction_grid": list(FRACTION_GRID)},
              },
              "development_rows_sha256": _digest([r for r in records if r["condition"] == original and r["partition"] != "test"]),
              "stage": "frozen_before_test_predictions_and_metrics",
              "operating_points": {arm: entry["operating_point"] for arm, entry in arms.items()},
              "rich_challenger": rich_challenger, "probabilistic_comparator": probability_comparator,
              "strongest_simple_baseline": strongest_simple, "cv": cv, "caller_baseline": caller_status,
              "training_policy": "original condition only; imputation/scaling fit on training only",
              "calibration_policy": "sigmoid fitted on calibration role only",
              "threshold_policy": "threshold role only; tied scores accepted together"}
    # Publish the exact frozen identity before evaluating any test outcomes.
    with (output_dir / "frozen-models.pkl").open("xb") as handle:
        handle.write(model_bytes)
    _write_json(output_dir / "frozen-models.json", frozen)
    freeze_sha = hashlib.sha256((output_dir / "frozen-models.json").read_bytes()).hexdigest()
    signal = any(gain >= 0.005 for gain in cv["mean_ap_gain_f1_minus_f0"].values())
    status = "proceed_to_p1" if signal else "stop_no_signal"
    if strongest_simple is None:
        status = "inconclusive_operating_point"
    report["gate"] = {"status": status, "proceed_to_p1": status == "proceed_to_p1",
                      "ap_gain_threshold": 0.005, "mean_ap_gain_f1_minus_f0": cv["mean_ap_gain_f1_minus_f0"],
                      "test_outcomes_used": False,
                      "interpretation": "Exploratory development signal only; advancement also requires external correctness gates",
                      "required_external_gates": ["native_correctness", "parquet_equivalence", "saved_evidence_replay"]}
    report["frozen_model"] = {"model_sha256": frozen["model_sha256"], "manifest_sha256": freeze_sha}
    report["cross_validation"] = cv
    report["caller_baseline"] = caller_status
    report["selected_comparators"] = {"rich_challenger": rich_challenger,
                    "probabilistic_comparator": probability_comparator, "strongest_simple": strongest_simple}
    prediction_started = time.perf_counter()
    predictions = predict_records(bundle, records)
    report["metrics"] = evaluate_predictions(predictions)
    report["costs"] = {"elapsed_seconds": time.perf_counter() - started,
                       "cpu_seconds": time.process_time() - cpu_started,
                       "cv_seconds": cv_seconds, "fit_calibrate_threshold_seconds_by_arm": fit_costs,
                       "prediction_and_metrics_seconds": time.perf_counter() - prediction_started,
                       "base_model_fits": 12, "calibrator_fits": 4 + int("B1" in arms),
                       "timing_scope": "model study only; extraction/acquisition/native verification excluded"}
    return _finish(output_dir, report, predictions)
