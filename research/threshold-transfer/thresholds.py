"""Freeze finite HG002 retention policies and evaluate them without refitting.

Recall here means retaining assessed positive candidates, including unresolved
pass-throughs. The genomic-window bootstrap describes empirical spatial
uncertainty. It is neither an iid/conformal guarantee nor a guarantee for a new
person, library, depth, or platform. Only retain-all has structural recall one,
and only for the supplied candidate list.
"""
from __future__ import annotations

import hashlib
import json
import math
import pickle
import time
from collections import defaultdict
from pathlib import Path

import study


MODELS = {"counts": "M10_pooled", "rich": "M11_pooled"}
DEFAULT_EDGES = (0, 5, 10, 20, 40, None)
DEFAULT_QUANTILES = (.99, .995, .999, 1.0)


def _canonical(value) -> bytes:
    return study._canonical(value)


def _write(path: Path, value) -> None:
    with path.open("xb") as handle:
        handle.write(_canonical(value))


def _digest(value) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _model_digest(bundle: dict) -> str:
    contents = {"feature_names": bundle["feature_names"],
                "models": {name: bundle["arms"][name]["model"] for name in MODELS.values()}}
    return hashlib.sha256(pickle.dumps(contents, protocol=5)).hexdigest()


def _conditions(protocol: dict) -> tuple[str, ...]:
    return tuple(protocol.get("development_conditions", ("original", "downsampled")))


def _validate(records: list[dict], bundle: dict) -> dict:
    names = study._validate_records(records)
    if not records:
        raise ValueError("nonempty records required")
    if bundle.get("schema") != "rosalind-depth-calibration-models-v1":
        raise ValueError("trusted depth-calibration model bundle required")
    if names != bundle.get("feature_names"):
        raise ValueError("record feature names differ from frozen models")
    for category, model_id in MODELS.items():
        entry = bundle.get("arms", {}).get(model_id, {})
        expected_group = "f0" if category == "counts" else "f1"
        if entry.get("feature_group") != expected_group or entry.get("variant") != "pooled":
            raise ValueError(f"expected fixed pooled model: {model_id}")
    identity = {}
    conditions = defaultdict(set)
    for row in records:
        depth = row["f0"].get("callable_depth")
        if type(depth) not in (float, int) or not math.isfinite(depth) or depth < 0:
            raise ValueError("finite nonnegative callable_depth required")
        key = row["candidate_id"]
        value = (row["contig"], row["pos0"], row["ref"], row["alt"], row["partition"], row["label"])
        if identity.setdefault(key, value) != value:
            raise ValueError("candidate identities or labels differ across derivatives")
        conditions[key].add(row["condition"])
    return {"feature_names": names, "record_count": len(records),
            "unique_candidates": len(identity),
            "unique_loci": len({(r.get("sample_id", "unspecified"), r["contig"], r["pos0"]) for r in records}),
            "unknown_record_count": sum(r["label"] is None for r in records),
            "depth_derivatives_are_independent_biological_samples": False}


def _scores(records: list[dict], bundle: dict) -> dict[str, list[float]]:
    """All inputs to prediction are saved numeric features; labels are unused."""
    dep = study._dependencies()
    result = {}
    for model_id in MODELS.values():
        entry = bundle["arms"][model_id]
        matrix = study._matrix(records, entry["feature_group"], bundle["feature_names"][entry["feature_group"]], dep["np"])
        values = study._raw_scores(entry["model"], matrix, dep["np"])
        result[model_id] = [float(value) for value in values]
        if any(not math.isfinite(value) for value in result[model_id]):
            raise ValueError("nonfinite frozen-model score")
    return result


def policy_decision(policy: dict, score: float, depth: float,
                    condition_source: str) -> dict:
    """Retain exact ties; pass-through is explicitly unresolved work."""
    if not math.isfinite(score) or not math.isfinite(depth) or depth < 0:
        raise ValueError("finite score and finite nonnegative depth required")
    kind = policy["kind"]
    if kind == "retain_all":
        return {"retained": True, "deferred": True, "reason": "retain_all_unresolved"}
    cutoff = policy.get("passthrough_below")
    if cutoff is not None and depth < cutoff:
        return {"retained": True, "deferred": True, "reason": "low_depth_pass_through"}
    reason = "score_threshold"
    if kind == "legacy":
        if condition_source not in policy["legacy_thresholds"]:
            raise ValueError("legacy diagnostic condition lacks explicit source mapping")
        threshold = policy["legacy_thresholds"][condition_source]
    elif kind == "global":
        threshold = policy["threshold"]
    elif kind == "depth_table":
        bins = [b for b in policy["bins"] if depth >= b["lower"]
                and (b["upper"] is None or depth < b["upper"])]
        if len(bins) != 1:
            raise ValueError("depth table must cover each depth exactly once")
        threshold = bins[0]["threshold"]
        reason = bins[0].get("reason", "depth_bin_threshold")
    else:
        raise ValueError("unknown policy kind")
    if threshold is None:
        return {"retained": True, "deferred": True, "reason": reason + "_unresolved"}
    if not math.isfinite(threshold):
        raise ValueError("finite threshold or explicit None pass-through required")
    retained = score >= threshold
    return {"retained": retained, "deferred": False,
            "reason": "model_retained" if retained else "model_filtered"}


def _structural(policy: dict) -> bool:
    return (policy["kind"] == "retain_all"
            or (policy["kind"] == "global" and policy["threshold"] is None)
            or (policy["kind"] == "depth_table" and all(b["threshold"] is None for b in policy["bins"])))


def _actions(records: list[dict], scores: dict, policies: dict, protocol: dict) -> dict:
    sources = protocol.get("condition_sources", {"original": "original", "downsampled": "downsampled",
                                                  "novel_low": "downsampled"})
    result = {}
    for policy_id, policy in policies.items():
        result[policy_id] = [policy_decision(policy, score, row["f0"]["callable_depth"],
                               sources.get(row["condition"], row["condition"]))
                            for row, score in zip(records, scores[policy["model_id"]])]
    return result


def predict_policy(rows: list[dict], bundle: dict, manifest: dict, label: str) -> list[dict]:
    """Label-free selected-path scoring, including matrix construction.

    Count-only rows need F0 only; rich rows need F1 only. This makes the same
    entry point suitable for matched standalone-table cost measurements.
    """
    if label not in MODELS:
        raise ValueError("policy label must be counts or rich")
    policy = manifest["policies"][manifest["selected"][label]]
    entry = bundle["arms"][policy["model_id"]]
    group = entry["feature_group"]
    if _structural(policy):
        return [policy_decision(policy, 0.0, row[group]["callable_depth"], row.get("condition", "")) for row in rows]
    dep = study._dependencies()
    matrix = study._matrix(rows, group, bundle["feature_names"][group], dep["np"])
    scores = study._raw_scores(entry["model"], matrix, dep["np"])
    return [policy_decision(policy, float(score), row[group]["callable_depth"], row.get("condition", ""))
            for row, score in zip(rows, scores)]


def _workload_metrics(rows: list[dict], actions: list[dict]) -> dict:
    known = [(r["label"], a) for r, a in zip(rows, actions) if r["label"] is not None]
    metric = study.decision_metrics([y for y, _ in known], [a["retained"] for _, a in known])
    metric["candidate_retention_recall"] = metric.pop("recall")
    metric["retention_precision"] = metric.pop("precision")
    metric.update({"candidate_count": len(rows), "unknown_labels": sum(r["label"] is None for r in rows),
        "all_retained_including_unknown": sum(a["retained"] for a in actions),
        "assessed_retained": sum(a["retained"] for _, a in known),
        "unknown_retained": sum(r["label"] is None and a["retained"] for r, a in zip(rows, actions)),
        "deferred_unresolved": sum(a["deferred"] for a in actions),
        "known_positive_deferred": sum(y == 1 and a["deferred"] for y, a in known),
        "known_negative_deferred": sum(y == 0 and a["deferred"] for y, a in known),
        "unknown_deferred": sum(r["label"] is None and a["deferred"] for r, a in zip(rows, actions)),
        "autonomous_retained_positive_candidates": sum(y == 1 and a["retained"] and not a["deferred"] for y, a in known),
        "autonomous_retained_negative_candidates": sum(y == 0 and a["retained"] and not a["deferred"] for y, a in known),
        "autonomous_variant_call_accuracy_established": False,
        "interpretation": "All retained candidates, including deferred and unknown candidates, remain downstream workload"})
    return metric


def simultaneous_recall_bounds(rows: list[dict], decisions_by_policy: dict,
                               protocol: dict, stage: str = "development") -> dict:
    """One-sided simultaneous centered block-bootstrap lower recall bands.

    Each resample uses identical genomic-window multiplicities for every model,
    finite policy, and depth condition. Selecting among this finite family is
    therefore included in the empirical maximum-deviation correction. Sparse
    error support is marked unsupported rather than treating bootstrap zeros as
    certainty. The correction does not establish coverage under person shift.
    """
    import numpy as np
    if not rows or not decisions_by_policy:
        raise ValueError("nonempty rows and finite policy family required")
    settings = protocol.get("bootstrap", {})
    block_bp = int(settings.get("block_bp", 100000))
    repetitions = int(settings.get("development_replicates" if stage == "development" else "fresh_replicates",
                                   10000 if stage == "development" else 2000))
    seed = int(settings.get("seed", 42))
    minimum_blocks = int(settings.get("minimum_nonempty_blocks", 20))
    minimum_errors = int(settings.get("minimum_error_blocks", 2))
    valid_fraction = float(settings.get("minimum_valid_fraction", .95))
    confidence = float(protocol.get("confidence_level", .95))
    if block_bp < 1 or repetitions < 1 or not 0 < confidence < 1 or not 0 < valid_fraction <= 1:
        raise ValueError("invalid uncertainty configuration")
    conditions = sorted({r["condition"] for r in rows})
    policy_ids = sorted(decisions_by_policy)
    for policy_id, entry in decisions_by_policy.items():
        if len(entry["decisions"]) != len(rows) or any(type(x) is not bool for x in entry["decisions"]):
            raise ValueError(f"aligned boolean decisions required: {policy_id}")
        if entry.get("structural") and not all(entry["decisions"]):
            raise ValueError("structural retain-all policy contains a filtered candidate")
    block_keys = sorted({(r.get("sample_id", "unspecified"), r["contig"], r["pos0"] // block_bp) for r in rows})
    block_index = {key: i for i, key in enumerate(block_keys)}
    ci = {c: i for i, c in enumerate(conditions)}
    positives = np.zeros((len(conditions), len(block_keys)), dtype=float)
    retained = np.zeros((len(policy_ids), len(conditions), len(block_keys)), dtype=float)
    for i, row in enumerate(rows):
        if row["label"] not in (0, 1, None):
            raise ValueError("binary or unknown labels required")
        if row["label"] != 1:
            continue
        c = ci[row["condition"]]
        b = block_index[(row.get("sample_id", "unspecified"), row["contig"], row["pos0"] // block_bp)]
        positives[c, b] += 1
        for j, policy_id in enumerate(policy_ids):
            retained[j, c, b] += decisions_by_policy[policy_id]["decisions"][i]
    totals = positives.sum(axis=1)
    point = np.divide(retained.sum(axis=2), totals[None, :],
                      out=np.full((len(policy_ids), len(conditions)), np.nan), where=totals[None, :] > 0)
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(len(block_keys), np.full(len(block_keys), 1 / len(block_keys)), size=repetitions)
    bootstrap_totals = weights @ positives.T
    valid = np.all(bootstrap_totals > 0, axis=1)
    valid_count = int(valid.sum())
    max_deviations = []
    if valid_count:
        selected_weights = weights[valid]
        for offset in range(0, valid_count, 250):
            w = selected_weights[offset:offset + 250]
            den = w @ positives.T
            numer = np.einsum("rb,pcb->rpc", w, retained)
            recalled = numer / den[:, None, :]
            # Basic bootstrap: error of estimate is approximated by estimate* - estimate.
            max_deviations.extend(np.max(recalled - point[None, :, :], axis=(1, 2)).tolist())
    margin = max(0.0, float(np.quantile(max_deviations, confidence, method="higher"))) if max_deviations else None
    result = {"stage": stage, "method": "one-sided simultaneous centered genomic-window bootstrap band",
              "confidence_level": confidence, "block_bp": block_bp, "block_count": len(block_keys),
              "requested_replicates": repetitions, "valid_replicates": valid_count,
              "invalid_replicates": repetitions - valid_count, "seed": seed,
              "family_size": len(policy_ids) * len(conditions), "policy_ids": policy_ids,
              "conditions": conditions, "simultaneous_margin": margin,
              "threshold_selection_accounting": "Maximum deviation over the entire prespecified finite policy/model/condition family",
              "iid_or_distribution_free_guarantee": False, "person_shift_guarantee": False,
              "support_guard": {"minimum_nonempty_positive_blocks": minimum_blocks,
                                "minimum_error_blocks_for_nonstructural": minimum_errors,
                                "minimum_valid_fraction": valid_fraction}, "policies": {}}
    for j, policy_id in enumerate(policy_ids):
        result["policies"][policy_id] = {"conditions": {}}
        for c, condition in enumerate(conditions):
            structural = bool(decisions_by_policy[policy_id].get("structural", False))
            positive_blocks = int(np.sum(positives[c] > 0))
            errors = positives[c] - retained[j, c]
            error_blocks = int(np.sum(errors > 0))
            reasons = []
            if totals[c] == 0:
                reasons.append("no assessed positive candidates")
            if not structural:
                if positive_blocks < minimum_blocks:
                    reasons.append("too few nonempty positive-bearing genomic windows")
                if error_blocks < minimum_errors:
                    reasons.append("too few error-bearing windows for bootstrap tail support")
                if valid_count < math.ceil(valid_fraction * repetitions):
                    reasons.append("too few valid simultaneous resamples")
                if margin is None:
                    reasons.append("bootstrap band unavailable")
            recall = float(point[j, c]) if totals[c] else None
            lower = (1.0 if structural and totals[c] else
                     max(0.0, recall - margin) if recall is not None and margin is not None else None)
            result["policies"][policy_id]["conditions"][condition] = {
                "recall": recall, "lower_bound": lower, "supported": not reasons,
                "reasons": reasons, "positive_count": int(totals[c]),
                "positive_block_count": positive_blocks, "error_block_count": error_blocks,
                "false_negatives": int(errors.sum()), "structural": structural,
                "bound_kind": "structural candidate-list retention" if structural else "empirical simultaneous spatial lower bound"}
    return result


def _candidate_policies(rows: list[dict], scores: dict, bundle: dict, protocol: dict) -> dict:
    quantiles = tuple(float(q) for q in protocol.get("policy_recall_targets", DEFAULT_QUANTILES))
    edges = tuple(protocol.get("depth_edges", DEFAULT_EDGES))
    minimum = int(protocol.get("minimum_calibration_positives_per_depth_bin", 300))
    low = float(protocol.get("low_depth_passthrough_below", 5))
    if edges != DEFAULT_EDGES or any(not 0 < q <= 1 for q in quantiles):
        raise ValueError("fixed depth edges and valid recall targets required")
    policies = {}
    known_positives = [i for i, r in enumerate(rows) if r["label"] == 1]
    if not known_positives:
        raise ValueError("calibration role requires assessed positives")
    for model_id in MODELS.values():
        common = {"model_id": model_id, "fit_role": "calibration", "condition_agnostic": True}
        legacy_thresholds = {condition: bundle["arms"][model_id]["operating_points"][condition]["0.99"]["threshold"]
                             for condition in _conditions(protocol)}
        for legacy_id, cutoff in (("existing_global", None), ("existing_global_low_depth_passthrough", low)):
            policy_id = f"{model_id}::{legacy_id}"
            policies[policy_id] = common | {"id": policy_id, "kind": "legacy", "passthrough_below": cutoff,
                "legacy_thresholds": legacy_thresholds, "eligible_for_selection": False,
                "condition_agnostic": False,
                "selection_exclusion": "Legacy thresholds were already selected using the HG002 threshold role"}
        for q in quantiles:
            suffix = format(q, ".12g")
            nonlow = [i for i in known_positives if rows[i]["f0"]["callable_depth"] >= low]
            threshold = (study.select_score_threshold([1] * len(nonlow), [scores[model_id][i] for i in nonlow], q)["threshold"]
                         if nonlow else None)
            policy_id = f"{model_id}::global_q{suffix}"
            policies[policy_id] = common | {"id": policy_id, "kind": "global", "threshold": threshold,
                "passthrough_below": low, "quantile_target": q, "eligible_for_selection": True,
                "calibration_positive_rows_above_passthrough": len(nonlow)}
            bins = []
            for lower, upper in zip(edges[:-1], edges[1:]):
                indexes = [i for i in known_positives if rows[i]["f0"]["callable_depth"] >= lower
                           and (upper is None or rows[i]["f0"]["callable_depth"] < upper)]
                sparse = len(indexes) < minimum or lower < low
                threshold = (None if sparse else study.select_score_threshold(
                    [1] * len(indexes), [scores[model_id][i] for i in indexes], q)["threshold"])
                bins.append({"lower": lower, "upper": upper, "threshold": threshold,
                             "calibration_positive_rows": len(indexes),
                             "reason": "low_depth_pass_through" if lower < low else
                                       "sparse_calibration_bin_pass_through" if sparse else "depth_bin_threshold"})
            policy_id = f"{model_id}::depth_table_q{suffix}"
            policies[policy_id] = common | {"id": policy_id, "kind": "depth_table", "bins": bins,
                "passthrough_below": low, "quantile_target": q, "eligible_for_selection": True}
        policy_id = f"{model_id}::retain_all"
        policies[policy_id] = common | {"id": policy_id, "kind": "retain_all", "eligible_for_selection": True}
    return policies


def _metrics_by_condition(rows: list[dict], actions: dict) -> dict:
    result = {}
    for condition in sorted({r["condition"] for r in rows}):
        indexes = [i for i, r in enumerate(rows) if r["condition"] == condition]
        result[condition] = {policy_id: _workload_metrics([rows[i] for i in indexes], [a[i] for i in indexes])
                             for policy_id, a in actions.items()}
    return result


def freeze_policies(hg002_records: list[dict], trusted_depth_bundle: dict,
                    output_dir: Path, protocol: dict) -> dict:
    """Construct policies on calibration rows, select on threshold rows, freeze."""
    start = time.perf_counter()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("frozen-policies.json", "development-report.json"):
        if (output_dir / name).exists():
            raise FileExistsError("frozen policy outputs already exist")
    inputs = _validate(hg002_records, trusted_depth_bundle)
    conditions = set(_conditions(protocol))
    calibration_rows = [r for r in hg002_records if r["partition"] == "calibration" and r["condition"] in conditions]
    selection_rows = [r for r in hg002_records if r["partition"] == "threshold" and r["condition"] in conditions]
    if {r["condition"] for r in calibration_rows} != conditions or {r["condition"] for r in selection_rows} != conditions:
        raise ValueError("all development conditions required in both calibration and threshold roles")
    # No test-role scoring occurs while constructing or selecting a policy.
    calibration_scores = _scores(calibration_rows, trusted_depth_bundle)
    policies = _candidate_policies(calibration_rows, calibration_scores, trusted_depth_bundle, protocol)
    selection_scores = _scores(selection_rows, trusted_depth_bundle)
    actions = _actions(selection_rows, selection_scores, policies, protocol)
    # Correct over every policy eligible for selection. Legacy diagnostics never
    # enter this search family and cannot change the nomination's uncertainty.
    decisions = {p: {"decisions": [a["retained"] for a in values], "structural": _structural(policies[p])}
                 for p, values in actions.items() if policies[p]["eligible_for_selection"]}
    bounds = simultaneous_recall_bounds(selection_rows, decisions, protocol, "development")
    metrics = _metrics_by_condition(selection_rows, actions)
    target = float(protocol.get("target_recall", .99))
    eligible, selected = {}, {}
    for category, model_id in MODELS.items():
        candidates = []
        for policy_id, policy in policies.items():
            if not policy["eligible_for_selection"]:
                eligible[policy_id] = False
                continue
            checks = bounds["policies"][policy_id]["conditions"]
            passes = policy["eligible_for_selection"] and all(checks[c]["supported"]
                and checks[c]["lower_bound"] is not None and checks[c]["lower_bound"] >= target for c in conditions)
            eligible[policy_id] = passes
            if passes and policy["model_id"] == model_id:
                candidates.append(policy_id)
        if not candidates:
            raise ValueError("no policy meets selection constraints, including structural retain-all")
        selected[category] = min(candidates, key=lambda p: (
            sum(metrics[c][p]["all_retained_including_unknown"] for c in conditions), p))
    development_report = {"schema": "rosalind-threshold-transfer-development-v1", "inputs": inputs,
        "metrics": metrics, "bounds": bounds, "eligible": eligible, "selected": selected,
        "construction_role": "calibration", "selection_role": "threshold",
        "selection_objective": "Least total retained workload including unknowns across both development conditions, then policy ID",
        "base_model_refits": 0, "threshold_family_test_outcomes_used": False,
        "calibration_rows_sha256": _digest(calibration_rows), "selection_rows_sha256": _digest(selection_rows),
        "wall_seconds": time.perf_counter() - start}
    _write(output_dir / "development-report.json", development_report)
    frozen = {"schema": "rosalind-threshold-transfer-policies-v1", "policies": policies, "selected": selected,
        "selected_structural_retain_all": {category: _structural(policies[policy_id]) for category, policy_id in selected.items()},
        "protocol_sha256": _digest(protocol), "model_sha256": _model_digest(trusted_depth_bundle),
        "feature_names": trusted_depth_bundle["feature_names"], "model_ids": MODELS,
        "calibration_rows_sha256": development_report["calibration_rows_sha256"],
        "selection_rows_sha256": development_report["selection_rows_sha256"],
        "development_evidence_sha256": _digest({key: development_report[key] for key in
            ("metrics", "bounds", "eligible", "selected", "calibration_rows_sha256", "selection_rows_sha256")}),
        "stage": "frozen_before_fresh_scoring", "base_model_refits": 0,
        "test_labels_used_for_fitting_or_selection": False,
        "bound_interpretation": "Empirical simultaneous spatial uncertainty, not a distribution-free 99% guarantee",
        "recall_estimand": "Retention of assessed positive candidates including deferred unresolved work"}
    frozen["policy_sha256"] = _digest(frozen)
    _write(output_dir / "frozen-policies.json", frozen)
    return frozen


def _paired_comparison(rows: list[dict], counts: list[dict], rich: list[dict], protocol: dict) -> dict:
    """Frozen decisions only; no evaluation-label threshold optimization."""
    import numpy as np
    settings = protocol.get("bootstrap", {})
    block_bp = int(settings.get("block_bp", 100000))
    reps = int(settings.get("fresh_replicates", 2000))
    rng = np.random.default_rng(int(settings.get("seed", 42)))
    blocks = defaultdict(list)
    for i, row in enumerate(rows):
        blocks[(row.get("sample_id", "unspecified"), row["contig"], row["pos0"] // block_bp)].append(i)
    vectors = []
    for key in sorted(blocks):
        indexes = blocks[key]
        positives = sum(rows[i]["label"] == 1 for i in indexes)
        ctp = sum(rows[i]["label"] == 1 and counts[i]["retained"] for i in indexes)
        rtp = sum(rows[i]["label"] == 1 and rich[i]["retained"] for i in indexes)
        cfp = sum(rows[i]["label"] == 0 and counts[i]["retained"] for i in indexes)
        rfp = sum(rows[i]["label"] == 0 and rich[i]["retained"] for i in indexes)
        vectors.append((positives, ctp, rtp, cfp, rfp))
    table = np.asarray(vectors, dtype=float)
    weights = rng.multinomial(len(table), np.full(len(table), 1 / len(table)), size=reps)
    sampled = weights @ table
    valid = sampled[:, 0] > 0
    differences = sampled[valid, 2] / sampled[valid, 0] - sampled[valid, 1] / sampled[valid, 0]
    gains = sampled[valid, 3] - sampled[valid, 4]
    relative = [float((v[3] - v[4]) / v[3]) for v in sampled[valid] if v[3] > 0]

    def interval(values):
        if len(values) == 0:
            return None
        return {"lower": float(np.quantile(values, .025)), "median": float(np.quantile(values, .5)),
                "upper": float(np.quantile(values, .975))}

    c, r = _workload_metrics(rows, counts), _workload_metrics(rows, rich)
    gain = c["false_positives"] - r["false_positives"]
    return {"counts": c, "rich": r,
            "absolute_fp_reduction": gain,
            "relative_fp_reduction": gain / c["false_positives"] if c["false_positives"] else None,
            "recall_rich_minus_counts": r["candidate_retention_recall"] - c["candidate_retention_recall"]
                if c["candidate_retention_recall"] is not None and r["candidate_retention_recall"] is not None else None,
            "assessed_workload_reduction": (c["assessed_retained"] - r["assessed_retained"]) / c["assessed_retained"] if c["assessed_retained"] else None,
            "all_workload_reduction": (c["all_retained_including_unknown"] - r["all_retained_including_unknown"]) / c["all_retained_including_unknown"] if c["all_retained_including_unknown"] else None,
            "bootstrap": {"block_count": len(table), "requested_replicates": reps,
                "valid_replicates": int(valid.sum()), "invalid_replicates": int((~valid).sum()),
                "relative_fp_valid_replicates": len(relative),
                "absolute_fp_reduction_95ci": interval(gains), "relative_fp_reduction_95ci": interval(relative),
                "recall_rich_minus_counts_95ci": interval(differences),
                "thresholds_refitted_in_bootstrap": False}}


def evaluate_gate(comparisons: dict, bounds: dict, selected: dict, protocol: dict) -> dict:
    gates = protocol.get("gates", {})
    required = tuple(protocol.get("required_fresh_conditions", ("original", "downsampled", "novel_low")))
    primary = protocol.get("primary_condition", "downsampled")
    target = float(protocol.get("target_recall", .99))
    tolerance = float(gates.get("maximum_recall_loss", .005))
    checks = {"all_required_conditions_present": all(c in comparisons for c in required)}
    for condition in required:
        for category in ("counts", "rich"):
            entry = bounds.get("policies", {}).get(selected[category], {}).get("conditions", {}).get(condition, {})
            checks[f"{condition}_{category}_absolute_recall"] = (entry.get("supported", False)
                and entry.get("lower_bound") is not None and entry["lower_bound"] >= target)
        comparison = comparisons.get(condition, {})
        bootstrap = comparison.get("bootstrap", {})
        settings = protocol.get("bootstrap", {})
        checks[f"{condition}_comparison_uncertainty_supported"] = (
            bootstrap.get("block_count", 0) >= int(settings.get("minimum_nonempty_blocks", 20))
            and bootstrap.get("valid_replicates", 0) >= math.ceil(float(settings.get("minimum_valid_fraction", .95))
                * int(settings.get("fresh_replicates", 2000))))
        delta = comparison.get("recall_rich_minus_counts")
        interval = comparison.get("bootstrap", {}).get("recall_rich_minus_counts_95ci")
        checks[f"{condition}_recall_noninferiority"] = (delta is not None and delta >= -tolerance
            and interval is not None and interval["lower"] > -tolerance)
    main = comparisons.get(primary, {})
    absolute = main.get("absolute_fp_reduction")
    relative = main.get("relative_fp_reduction")
    interval = main.get("bootstrap", {}).get("absolute_fp_reduction_95ci")
    checks["primary_absolute_fp_gain"] = absolute is not None and absolute >= int(gates.get("minimum_absolute_fp_reduction", 10))
    checks["primary_relative_fp_gain"] = relative is not None and relative >= float(gates.get("minimum_relative_fp_reduction", .25))
    checks["primary_fp_gain_ci_positive"] = interval is not None and interval["lower"] > 0
    for key, setting in (("assessed_workload_reduction", "minimum_assessed_retained_reduction"),
                         ("all_workload_reduction", "minimum_total_retained_reduction")):
        checks["primary_" + key] = main.get(key) is not None and main[key] >= float(gates.get(setting, .01))
    positive = all(checks.values())
    return {"status": "predictive_and_workload_screen_passed" if positive else "prespecified_screen_not_met",
            "predictive_and_workload_screen_passed": positive, "checks": checks,
            "nominate_fresh_workflow_validation": positive,
            "selective_escalation_authorized": False, "proceed_to_p1": False,
            "total_workflow_cost_advantage_established": False,
            "reason": "Full acquisition, deferred review and fallback costs remain external requirements",
            "recall_estimand": "Candidate retention including unresolved pass-throughs",
            "empirical_spatial_bounds_are_not_person_shift_guarantees": True}


def _diagnostic_baselines(records: list[dict], bundle: dict, protocol: dict) -> dict:
    result = {}
    sources = protocol.get("condition_sources", {"original": "original", "downsampled": "downsampled",
                                                  "novel_low": "downsampled"})
    for condition in sorted({r["condition"] for r in records}):
        rows = [r for r in records if r["condition"] == condition]
        result[condition] = {}
        point = bundle.get("arms", {}).get("B0", {}).get("operating_points", {}).get(sources.get(condition, condition), {}).get("0.99")
        if point is None or not point.get("eligible"):
            result[condition]["B0"] = {"available": False, "eligible_for_nomination": False,
                "reason": "No eligible frozen HG002 depth/fraction operating point for mapped condition"}
        else:
            actions = []
            for row in rows:
                depth, fraction = study._b0_values(row, point["depth_key"], point["fraction_key"])
                retained = depth is not None and fraction is not None and depth >= point["depth_threshold"] and fraction >= point["fraction_threshold"]
                actions.append({"retained": retained, "deferred": False})
            metrics = _workload_metrics(rows, actions)
            result[condition]["B0"] = {"available": True, "metrics": metrics,
                "operating_point": point, "source_condition": sources.get(condition, condition),
                "eligible_for_nomination": False,
                "observed_recall_target_met": metrics["candidate_retention_recall"] is not None and metrics["candidate_retention_recall"] >= protocol.get("target_recall", .99)}
        if not all("caller_score" in row for row in rows):
            result[condition]["caller_emitted"] = {"available": False,
                "reason": "Caller emission not supplied on every candidate; no substitution", "eligible_for_nomination": False}
        else:
            if any(row["caller_score"] is not None and (not math.isfinite(float(row["caller_score"])) or float(row["caller_score"]) < 0) for row in rows):
                raise ValueError("caller QUAL must be finite nonnegative or explicit None")
            actions = [{"retained": row["caller_score"] is not None, "deferred": False} for row in rows]
            metrics = _workload_metrics(rows, actions)
            result[condition]["caller_emitted"] = {"available": True, "metrics": metrics,
                "decision": "Retain candidates emitted by the fresh caller, without a learned QUAL threshold",
                "eligible_for_nomination": False,
                "observed_recall_target_met": metrics["candidate_retention_recall"] is not None and metrics["candidate_retention_recall"] >= protocol.get("target_recall", .99),
                "scope": "Conditional on candidate membership; caller-only variants and de novo recovery reported by acquisition adapter"}
    return result


def evaluate_frozen(fresh_records: list[dict], trusted_depth_bundle: dict,
                    frozen_policy: dict, output_dir: Path, protocol: dict) -> dict:
    """Evaluate all frozen policies; never select thresholds from fresh labels."""
    start, cpu = time.perf_counter(), time.process_time()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("predictions.jsonl", "report.json"):
        if (output_dir / name).exists():
            raise FileExistsError("fresh evaluation outputs already exist")
    inputs = _validate(fresh_records, trusted_depth_bundle)
    if frozen_policy.get("schema") != "rosalind-threshold-transfer-policies-v1":
        raise ValueError("unsupported frozen policy schema")
    claimed = frozen_policy["policy_sha256"]
    if claimed != _digest({k: v for k, v in frozen_policy.items() if k != "policy_sha256"}):
        raise ValueError("frozen policy bytes or fields changed")
    if frozen_policy["protocol_sha256"] != _digest(protocol) or frozen_policy["model_sha256"] != _model_digest(trusted_depth_bundle):
        raise ValueError("protocol or base model differs from frozen policy identity")
    scores = _scores(fresh_records, trusted_depth_bundle)
    policies, selected = frozen_policy["policies"], frozen_policy["selected"]
    actions = _actions(fresh_records, scores, policies, protocol)
    # Fresh labels evaluate only the two already selected policies for primary
    # uncertainty; all other frozen policies remain descriptive, not new winners.
    decisions = {p: {"decisions": [a["retained"] for a in values], "structural": _structural(policies[p])}
                 for p, values in actions.items() if p in selected.values()}
    bounds = simultaneous_recall_bounds(fresh_records, decisions, protocol, "fresh")
    metrics = _metrics_by_condition(fresh_records, actions)
    comparisons, depth_metrics = {}, {}
    for condition in sorted(metrics):
        indexes = [i for i, r in enumerate(fresh_records) if r["condition"] == condition]
        rows = [fresh_records[i] for i in indexes]
        comparisons[condition] = _paired_comparison(rows, [actions[selected["counts"]][i] for i in indexes],
                                                    [actions[selected["rich"]][i] for i in indexes], protocol)
        depth_metrics[condition] = {}
        for category, policy_id in selected.items():
            depth_metrics[condition][category] = []
            for low, high in zip(DEFAULT_EDGES[:-1], DEFAULT_EDGES[1:]):
                indexes = [i for i, r in enumerate(fresh_records) if r["condition"] == condition
                           and r["f0"]["callable_depth"] >= low and (high is None or r["f0"]["callable_depth"] < high)]
                depth_metrics[condition][category].append({"lower": low, "upper": high,
                    "metrics": _workload_metrics([fresh_records[i] for i in indexes], [actions[policy_id][i] for i in indexes]),
                    "bound_support": "descriptive depth slice; aggregate empirical bound does not certify this stratum"})
    with (output_dir / "predictions.jsonl").open("xb") as handle:
        for i, row in enumerate(fresh_records):
            result = {key: row.get(key) for key in ("candidate_id", "sample_id", "contig", "pos0", "ref", "alt", "partition", "condition", "label", "label_reason")}
            result["scores"] = {model_id: values[i] for model_id, values in scores.items()}
            result["policies"] = {policy_id: values[i] for policy_id, values in actions.items()}
            handle.write(_canonical(result))
    report = {"schema": "rosalind-threshold-transfer-evaluation-v1", "inputs": inputs,
        "protocol_sha256": _digest(protocol), "policy_sha256": claimed,
        "selected": selected, "metrics": metrics, "bounds": bounds, "comparisons": comparisons,
        "diagnostic_baselines": _diagnostic_baselines(fresh_records, trusted_depth_bundle, protocol),
        "selected_measured_depth_slices": depth_metrics, "gate": evaluate_gate(comparisons, bounds, selected, protocol),
        "base_model_refits": 0, "fresh_label_threshold_fitting": False,
        "costs": {"wall_seconds": time.perf_counter() - start, "cpu_seconds": time.process_time() - cpu,
                  "scope": "Saved-feature scoring, policy application, metrics and bootstrap; acquisition/review/fallback excluded"},
        "limitations": ["Candidate-list retention is not end-to-end variant sensitivity",
            "Pass-throughs are unresolved workload, not autonomous correct variant calls",
            "One-sided bounds summarize spatial variability without guaranteeing person or platform shift",
            "Policies were fixed on HG002; no replacement winner is selected using fresh outcomes"]}
    _write(output_dir / "report.json", report)
    return report
