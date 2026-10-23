"""Exploratory depth calibration using the immutable HG002 P0 feature rows.

The runner verifies and supplies the trusted original pickle and makes the P0
``study`` module importable. This module never reads alignments or runs a caller.
The saved test regions have already been inspected: their results are
descriptive, even though all choices in this experiment precede test scoring.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import pickle
import statistics
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import study


VARIANTS = ("frozen", "retuned", "recalibrated", "pooled")
OUTPUTS = ("frozen-models.pkl", "frozen-models.json", "predictions.jsonl",
           "report.json", "features-f0.jsonl.gz", "features-f1.jsonl.gz")
META_KEYS = ("candidate_id", "contig", "pos0", "ref", "alt", "partition",
             "block_id", "condition")


def _clock() -> tuple[float, float]:
    return time.perf_counter(), time.process_time()


def _elapsed(start: tuple[float, float]) -> dict:
    return {"wall_seconds": time.perf_counter() - start[0],
            "cpu_seconds": time.process_time() - start[1]}


def _target_key(target: float) -> str:
    return format(target, ".12g")


def _canonical(value: Any) -> bytes:
    return study._canonical(value)


def _write(path: Path, value: Any) -> None:
    with path.open("xb") as handle:
        handle.write(_canonical(value))


def _validate_inputs(records: list[dict], original_bundle: dict,
                     original: str, reduced: str) -> dict:
    names = study._validate_records(records)
    if not records or original == reduced:
        raise ValueError("nonempty records and two distinct conditions required")
    if original_bundle.get("schema") != "rosalind-evidence-sufficiency-models-v1":
        raise ValueError("unsupported trusted original model bundle")
    if original_bundle.get("feature_names") != names:
        raise ValueError("saved feature names differ from original model")
    if original_bundle.get("original_condition") != original:
        raise ValueError("original model condition differs from experiment")
    conditions = {row["condition"] for row in records}
    if conditions != {original, reduced}:
        raise ValueError("exactly original and reduced-depth conditions required")
    candidates = defaultdict(dict)
    for row in records:
        depth = row["f0"].get("callable_depth")
        if type(depth) not in (int, float) or not math.isfinite(depth) or depth < 0:
            raise ValueError("finite nonnegative callable_depth required")
        candidates[row["candidate_id"]][row["condition"]] = row
    for pair in candidates.values():
        if set(pair) != conditions:
            raise ValueError("every candidate must have both paired depth conditions")
        if pair[original]["label"] != pair[reduced]["label"]:
            raise ValueError("paired candidates have different labels")
    counts = {}
    for partition in study.PARTITIONS:
        counts[partition] = {}
        for condition in (original, reduced):
            rows = [r for r in records if r["partition"] == partition
                    and r["condition"] == condition]
            counts[partition][condition] = study.class_counts(rows)
            if partition != "test" and (not counts[partition][condition]["positive"]
                                         or not counts[partition][condition]["negative"]):
                raise ValueError(f"both labelled classes required: {partition}/{condition}")
    for base, group, learner in study.ARMS:
        entry = original_bundle.get("arms", {}).get(base, {})
        if entry.get("feature_group") != group or entry.get("learner") != learner:
            raise ValueError(f"original bundle lacks expected arm {base}")
        if not entry.get("operating_point", {}).get("eligible"):
            raise ValueError(f"original arm {base} has no eligible frozen threshold")
    return {"feature_names": names, "partition_condition_counts": counts,
            "record_count_including_depth_derivatives": len(records),
            "unique_candidate_count": len(candidates),
            "unique_genomic_positions": len({(r["contig"], r["pos0"]) for r in records}),
            "paired_depth_derivatives_are_independent_samples": False}


def _raw(entry: dict, rows: list[dict], names: dict, dep: dict):
    if not rows:
        return dep["np"].empty(0, dtype=float)
    return study._raw_scores(entry["model"], study._matrix(
        rows, entry["feature_group"], names[entry["feature_group"]], dep["np"]), dep["np"])


def _fit_sigmoid(raw, labels: list[int], dep: dict):
    """Fit the existing sigmoid specification; reject any reversal of ranking."""
    model = dep["LogisticRegression"](C=1e6, max_iter=2000,
                                       solver="lbfgs", random_state=42)
    model.fit(raw.reshape(-1, 1), labels)
    slope = float(model.coef_[0, 0])
    if not math.isfinite(slope) or slope <= 0:
        raise ValueError("sigmoid calibration must have a strictly positive slope")
    return model


def _depth_matrix(raw, rows: list[dict], dep: dict):
    np = dep["np"]
    depth = np.log1p(np.asarray([r["f0"]["callable_depth"] for r in rows], dtype=float))
    return np.column_stack((raw, depth, raw * depth))


def _fit_depth_calibrator(raw, rows: list[dict], dep: dict):
    calibrator = dep["Pipeline"]([
        ("scaler", dep["StandardScaler"]()),
        ("classifier", dep["LogisticRegression"](
            C=1e6, max_iter=2000, solver="lbfgs", random_state=42)),
    ])
    calibrator.fit(_depth_matrix(raw, rows, dep), [r["label"] for r in rows])
    return calibrator


def _scores(entry: dict, rows: list[dict], condition: str, names: dict,
            dep: dict) -> tuple[list[float | None], list[float | None]]:
    """Return decision scores and probabilities without consulting labels."""
    if entry.get("kind") == "threshold_grid":
        return [None] * len(rows), [None] * len(rows)
    if entry.get("kind") == "caller":
        if any("caller_score" not in r for r in rows):
            raise ValueError("explicit original caller score required")
        raw = dep["np"].asarray([0.0 if r["caller_score"] is None
                                else float(r["caller_score"]) for r in rows])
        probabilities = entry["calibrators"][condition].predict_proba(raw.reshape(-1, 1))[:, 1]
        return [float(p) for p in probabilities], [float(p) for p in probabilities]
    raw = _raw(entry, rows, names, dep)
    calibrator = entry["calibrators"][condition]
    if entry["variant"] == "depth_calibrated":
        probabilities = calibrator.predict_proba(_depth_matrix(raw, rows, dep))[:, 1]
        # decision_function avoids artificial ties from saturated probabilities.
        score = calibrator.decision_function(_depth_matrix(raw, rows, dep))
    else:
        probabilities = calibrator.predict_proba(raw.reshape(-1, 1))[:, 1]
        score = probabilities if entry["variant"] == "frozen" else raw
    return [float(s) for s in score], [float(p) for p in probabilities]


def _decisions(entry: dict, rows: list[dict], scores: list,
               condition: str, target: float) -> list[bool | None]:
    operating = entry["operating_points"][condition][_target_key(target)]
    if not operating["eligible"]:
        return [None] * len(rows)
    if entry.get("kind") == "threshold_grid":
        result = []
        for row in rows:
            depth, fraction = study._b0_values(row, operating["depth_key"],
                                               operating["fraction_key"])
            result.append(depth is not None and fraction is not None
                          and depth >= operating["depth_threshold"]
                          and fraction >= operating["fraction_threshold"])
        return result
    return [s >= operating["threshold"] for s in scores]


def _metrics(rows: list[dict], scores: list, probabilities: list,
             accepted: list) -> dict:
    known = [i for i, row in enumerate(rows)
             if row["label"] is not None and accepted[i] is not None]
    labels = [rows[i]["label"] for i in known]
    prob = [probabilities[i] for i in known]
    result = study.decision_metrics(labels, [accepted[i] for i in known],
                                   prob if all(p is not None for p in prob) else None)
    rank_scores = [scores[i] for i in known]
    if all(s is not None for s in rank_scores):
        result["decision_score_average_precision"] = study.average_precision(labels, rank_scores)
    else:
        result["decision_score_average_precision"] = None
    result.update({"candidate_count": len(rows),
                   "unknown_labels": sum(r["label"] is None for r in rows),
                   "ineligible_predictions": sum(a is None for a in accepted),
                   "unknown_accepted": sum(r["label"] is None and a is True
                                           for r, a in zip(rows, accepted))})
    return result


def _depth_bins(rows: list[dict], values: dict, target: float) -> list[dict]:
    """Descriptive coverage slices; no fitting or operating-point selection."""
    boundaries = (0, 5, 10, 20, 40, None)
    result = []
    for low, high in zip(boundaries[:-1], boundaries[1:]):
        selected = [i for i, row in enumerate(rows)
                    if row["f0"]["callable_depth"] >= low
                    and (high is None or row["f0"]["callable_depth"] < high)]
        result.append({"callable_depth_lower_inclusive": low,
                       "callable_depth_upper_exclusive": high,
                       "upper_unbounded": high is None,
                       "metrics": _metrics([rows[i] for i in selected],
                                           [values["score"][i] for i in selected],
                                           [values["probability"][i] for i in selected],
                                           [values["decisions"][_target_key(target)][i] for i in selected])})
    return result


def descriptive_matched_recall(labels: list[int], scores: list[float],
                               target: float) -> dict:
    """Empirical curve point chosen using evaluation labels, never deployable."""
    operating = study.select_score_threshold(labels, scores, target)
    result = dict(operating)
    result["threshold_selected_on_evaluation_labels"] = True
    result["deployable_threshold_validation"] = False
    if operating["eligible"]:
        result.update(study.decision_metrics(labels,
                                            [s >= operating["threshold"] for s in scores]))
    return result


def _paired_bootstrap(rows: list[dict], count_score: list[float],
                      rich_score: list[float], count_accepted: list[bool],
                      rich_accepted: list[bool], target: float,
                      configuration: dict, dep: dict) -> dict:
    """Resample shared 100kb blocks, with all score ties retained per replicate."""
    np = dep["np"]
    block_bp = int(configuration.get("block_bp", 100000))
    replicates = int(configuration.get("replicates", 2000))
    seed = int(configuration.get("seed", 42))
    if block_bp <= 0 or replicates < 1:
        raise ValueError("invalid bootstrap block size or replicates")
    blocks = defaultdict(list)
    for i, row in enumerate(rows):
        blocks[(row["contig"], row["pos0"] // block_bp)].append(i)
    block_ids = sorted(blocks)
    arrays = [np.asarray(blocks[key], dtype=int) for key in block_ids]
    if len(arrays) < 2:
        return {"eligible": False, "reason": "at least two genomic blocks required",
                "block_count": len(arrays), "requested_replicates": replicates}
    labels = np.asarray([r["label"] for r in rows], dtype=int)
    count_score = np.asarray(count_score, dtype=float)
    rich_score = np.asarray(rich_score, dtype=float)
    count_accepted = np.asarray(count_accepted, dtype=bool)
    rich_accepted = np.asarray(rich_accepted, dtype=bool)
    random = np.random.default_rng(seed)
    absolute_gains, relative_gains, fixed_recall_deltas, matched_recall_deltas = [], [], [], []
    skipped = 0
    for _ in range(replicates):
        sample = np.concatenate([arrays[int(i)] for i in random.integers(0, len(arrays), len(arrays))])
        y = labels[sample]
        positives = int(y.sum())
        if positives == 0 or positives == len(y):
            skipped += 1
            continue
        mask = y == 1
        rank = max(0, math.ceil(target * positives) - 1)
        c = count_score[sample]
        r = rich_score[sample]
        ct = np.sort(c[mask])[::-1][rank]
        rt = np.sort(r[mask])[::-1][rank]
        ca, ra = c >= ct, r >= rt
        cf, rf = int(np.sum(ca & ~mask)), int(np.sum(ra & ~mask))
        absolute_gains.append(cf - rf)
        if cf:
            relative_gains.append((cf - rf) / cf)
        matched_recall_deltas.append(float(np.sum(ra & mask) - np.sum(ca & mask)) / positives)
        fixed_recall_deltas.append(float(np.sum(rich_accepted[sample] & mask)
                                         - np.sum(count_accepted[sample] & mask)) / positives)

    def interval(values):
        if not values:
            return None
        return {"lower": float(np.percentile(values, 2.5)),
                "median": float(np.percentile(values, 50)),
                "upper": float(np.percentile(values, 97.5))}

    return {"eligible": bool(absolute_gains), "method": "paired genomic-block percentile bootstrap",
            "block_bp": block_bp, "block_count": len(arrays),
            "block_ids": [f"{contig}:{index * block_bp}-{(index + 1) * block_bp}"
                          for contig, index in block_ids],
            "seed": seed, "requested_replicates": replicates,
            "valid_replicates": len(absolute_gains), "skipped_replicates": skipped,
            "relative_fp_gain_valid_replicates": len(relative_gains),
            "absolute_fp_reduction_counts_minus_rich_95ci": interval(absolute_gains),
            "relative_fp_reduction_counts_minus_rich_95ci": interval(relative_gains),
            "fixed_threshold_recall_rich_minus_counts_95ci": interval(fixed_recall_deltas),
            "matched_recall_rich_minus_counts_95ci": interval(matched_recall_deltas),
            "rethreshold_each_replicate": True, "preserves_score_ties": True,
            "limitation": "Conditional one-specimen spatial uncertainty; excludes donor, library, platform, selection and shift uncertainty"}


def _feature_storage(records: list[dict], names: dict, output: Path) -> dict:
    result = {}
    for group in ("f0", "f1"):
        path = output / f"features-{group}.jsonl.gz"
        uncompressed = 0
        started = _clock()
        # Identical metadata, row order, JSON format, compression level and mtime.
        # Labels and label reasons do not belong in a reusable predictor artifact.
        with path.open("xb") as destination:
            with gzip.GzipFile(filename="", fileobj=destination, mode="wb",
                               compresslevel=9, mtime=0) as handle:
                for row in records:
                    encoded = _canonical({key: row[key] for key in META_KEYS}
                                         | {"features": row[group]})
                    uncompressed += len(encoded)
                    handle.write(encoded)
        result[group] = {"feature_count": len(names[group]), "row_count": len(records),
                         "gzip_bytes": path.stat().st_size, "uncompressed_bytes": uncompressed,
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "creation_cost": _elapsed(started)}
    result["rich_to_counts_gzip_ratio"] = result["f1"]["gzip_bytes"] / result["f0"]["gzip_bytes"]
    result["scope"] = "Canonical standalone feature tables only; evidence extraction and raw fallback costs excluded"
    return result


def _score_costs(arms: dict, selected: dict, rows: list[dict], condition: str,
                  names: dict, target: float, rounds: int, dep: dict,
                  output_dir: Path) -> dict:
    if rounds < 1:
        raise ValueError("timing rounds must be positive")
    order = (selected["counts"], selected["rich"])
    # Force identical candidate rows, include Python matrix assembly and decisions.
    for arm in order:
        scores, _ = _scores(arms[arm], rows, condition, names, dep)
        _decisions(arms[arm], rows, scores, condition, target)
    timings = {arm: [] for arm in order}
    for round_index in range(rounds):
        for arm in order if round_index % 2 == 0 else tuple(reversed(order)):
            start = _clock()
            scores, _ = _scores(arms[arm], rows, condition, names, dep)
            accepted = _decisions(arms[arm], rows, scores, condition, target)
            duration = _elapsed(start)
            duration["accepted_count_including_unknown"] = sum(a is True for a in accepted)
            timings[arm].append(duration)
    result = {"condition": condition, "candidate_count_including_unknown": len(rows),
              "rounds": rounds, "warmup_rounds": 1,
              "order": "counts/rich then rich/counts, alternating",
              "scope": "model inference, matrix assembly, calibration and decision; in-memory rows; excludes I/O and evidence acquisition",
              "arms": {}}
    for arm, measurements in timings.items():
        result["arms"][arm] = {"rounds": measurements,
            "median_wall_seconds": statistics.median(v["wall_seconds"] for v in measurements),
            "median_cpu_seconds": statistics.median(v["cpu_seconds"] for v in measurements),
            "wall_seconds_range": [min(v["wall_seconds"] for v in measurements),
                                   max(v["wall_seconds"] for v in measurements)],
            "cpu_seconds_range": [min(v["cpu_seconds"] for v in measurements),
                                  max(v["cpu_seconds"] for v in measurements)]}
    for measure in ("wall", "cpu"):
        denominator = result["arms"][selected["counts"]][f"median_{measure}_seconds"]
        result[f"rich_to_counts_{measure}_ratio"] = (
            result["arms"][selected["rich"]][f"median_{measure}_seconds"] / denominator
            if denominator else None)
    # A second counterbalanced measurement starts from standalone compressed
    # features. Both paths scan the same full pool and select the same test rows.
    # This captures actual decompression/JSON costs separately from warm inference.
    identities = {(row["candidate_id"], row["condition"]) for row in rows}

    def parse_and_score(arm):
        entry = arms[arm]
        group = entry["feature_group"]
        start = _clock()
        parsed, scanned = [], 0
        with gzip.open(output_dir / f"features-{group}.jsonl.gz", "rt", encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                scanned += 1
                if (row["candidate_id"], row["condition"]) not in identities:
                    continue
                features = row.pop("features")
                row[group] = features
                if group == "f1":
                    row["f0"] = {key: features[key] for key in names["f0"]}
                parsed.append(row)
        parse_cost = _elapsed(start)
        if [(r["candidate_id"], r["condition"]) for r in parsed] != [
                (r["candidate_id"], r["condition"]) for r in rows]:
            raise AssertionError("parsed timing candidate set or order differs")
        scores, _ = _scores(entry, parsed, condition, names, dep)
        accepted = _decisions(entry, parsed, scores, condition, target)
        overall_cost = _elapsed(start)
        return overall_cost | {"parse_only": parse_cost, "input_rows_scanned": scanned,
                               "accepted_count_including_unknown": sum(a is True for a in accepted)}

    for arm in order:
        parse_and_score(arm)
    parsed_timings = {arm: [] for arm in order}
    for round_index in range(rounds):
        for arm in order if round_index % 2 == 0 else tuple(reversed(order)):
            parsed_timings[arm].append(parse_and_score(arm))
    result["standalone_parse_and_inference"] = {
        "scope": "Read and decompress entire same candidate pool, parse JSON, select the same test candidates, build matrices, calibrate and decide",
        "warmup_rounds": 1, "rounds": rounds, "arms": {}}
    for arm, measurements in parsed_timings.items():
        result["standalone_parse_and_inference"]["arms"][arm] = {
            "rounds": measurements,
            "median_wall_seconds": statistics.median(v["wall_seconds"] for v in measurements),
            "median_cpu_seconds": statistics.median(v["cpu_seconds"] for v in measurements),
            "wall_seconds_range": [min(v["wall_seconds"] for v in measurements),
                                   max(v["wall_seconds"] for v in measurements)],
            "cpu_seconds_range": [min(v["cpu_seconds"] for v in measurements),
                                  max(v["cpu_seconds"] for v in measurements)]}
    for measure in ("wall", "cpu"):
        values = result["standalone_parse_and_inference"]["arms"]
        denominator = values[selected["counts"]][f"median_{measure}_seconds"]
        result["standalone_parse_and_inference"][f"rich_to_counts_{measure}_ratio"] = (
            values[selected["rich"]][f"median_{measure}_seconds"] / denominator if denominator else None)
    return result


def evaluate_gain_gate(comparison: dict, costs: dict, protocol: dict) -> dict:
    """Nominate fresh validation only. Sunk saved-evidence costs cannot earn a router."""
    gates = protocol.get("gates", {})
    if gates.get("confidence_level", .95) != .95 or gates.get("requires_positive_fp_gain_ci_lower", True) is not True:
        raise ValueError("frozen screen requires 95% intervals with positive FP gain lower bound")
    relative_min = float(gates.get("minimum_relative_fp_reduction", .25))
    absolute_min = int(gates.get("minimum_absolute_fp_reduction", 10))
    recall_tolerance = float(gates.get("maximum_recall_loss", .005))
    matched = comparison["matched_recall"]
    fixed = comparison["fixed_operating_points"]
    bootstrap = comparison["bootstrap"]
    absolute = matched["absolute_fp_reduction_counts_minus_rich"]
    relative = matched["relative_fp_reduction_counts_minus_rich"]
    gain_ci = bootstrap.get("absolute_fp_reduction_counts_minus_rich_95ci")
    recall_ci = bootstrap.get("fixed_threshold_recall_rich_minus_counts_95ci")
    configuration = protocol.get("bootstrap", {})
    adequate = (bootstrap.get("eligible", False)
                and bootstrap.get("block_count", 0) >= configuration.get("minimum_nonempty_blocks", 20)
                and bootstrap.get("valid_replicates", 0) >= configuration.get("minimum_valid_replicates", 1900))
    costs_recorded = bool(costs.get("feature_storage") and costs.get("selected_inference"))
    checks = {
        "bootstrap_adequate": adequate,
        "saved_cost_measurements_present": costs_recorded,
        "relative_fp_reduction": relative is not None and relative >= relative_min,
        "absolute_fp_reduction": absolute >= absolute_min,
        "positive_absolute_fp_gain_ci": adequate and gain_ci is not None and gain_ci["lower"] > 0,
        "matched_attained_recall_not_materially_lower":
            matched["recall_rich_minus_counts"] >= -recall_tolerance,
        "frozen_operating_point_recall_not_materially_lower":
            fixed["recall_rich_minus_counts"] >= -recall_tolerance,
        "frozen_operating_point_recall_ci":
            adequate and recall_ci is not None and recall_ci["lower"] > -recall_tolerance,
    }
    positive = all(checks.values())
    point_margin_pass = checks["relative_fp_reduction"] and checks["absolute_fp_reduction"]
    uncertainty_pass = checks["positive_absolute_fp_gain_ci"] and checks["frozen_operating_point_recall_ci"]
    status = ("nominate_fresh_validation" if positive else
              "inconclusive_validation_screen" if not adequate or not costs_recorded or (point_margin_pass and not uncertainty_pass)
              else "stop_prespecified_screen_not_met")
    return {"status": status,
            "nominate_fresh_validation": positive, "checks": checks,
            "criteria": {"minimum_relative_fp_reduction": relative_min,
                         "minimum_absolute_fp_reduction": absolute_min,
                         "maximum_recall_loss": recall_tolerance,
                         "absolute_fp_gain_ci_lower_must_exceed": 0},
            "test_outcomes_used": True, "exploratory_reused_test_only": True,
            "selective_escalation_authorized": False, "p1_authorized": False,
            "proceed_to_selective_escalation": False, "proceed_to_p1": False,
            "rich_evidence_total_cost_advantage_established": False,
            "costs_recorded": costs_recorded,
            "limitations": ["Saved HG002 test evidence was already inspected in P0",
                            "Rich extraction costs are sunk; current timing is incremental replay only",
                            "No reduced-depth caller rerun or fallback workflow measured",
                            "No untouched independent specimen or library validation",
                            "Candidate list was generated at original depth; no de novo low-depth sensitivity claim"],
            "next_action": ("Pre-register fresh cost and independent-specimen validation before considering selective escalation"
                            if positive else "Retain the counts baseline; do not implement selective escalation; this is not evidence of equivalence")}


def run_experiment(records: list[dict], original_bundle: dict,
                   output_dir: Path, protocol: dict) -> dict:
    """Fit development-only corrections, freeze, then describe reused test results."""
    began = _clock()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any((output_dir / name).exists() for name in OUTPUTS):
        raise FileExistsError("experiment outputs exist; use a new run directory")
    original = protocol.get("original_condition", "original")
    reduced = protocol.get("reduced_condition", "downsampled")
    main_target = float(protocol.get("target_recall", .99))
    targets = tuple(dict.fromkeys([main_target] + [float(x) for x in
                          protocol.get("descriptive_recall_targets", [.95, .995])]))
    if main_target != .99 or any(not 0 < target <= 1 for target in targets):
        raise ValueError("primary target must be .99 and all targets in (0, 1]")
    records = sorted(records, key=lambda row: (row["candidate_id"], row["condition"]))
    inputs = _validate_inputs(records, original_bundle, original, reduced)
    names = inputs["feature_names"]
    dep = study._dependencies()
    conditions = (original, reduced)
    development = {(role, condition): study._labelled([r for r in records
                    if r["partition"] == role and r["condition"] == condition])
                   for role in study.PARTITIONS[:-1] for condition in conditions}
    pooled_train = development[("train", original)] + development[("train", reduced)]
    pooled_calibration = development[("calibration", original)] + development[("calibration", reduced)]
    arms, fits, rank_equivalence = {}, {}, {}
    base_fits = calibrator_fits = 0
    for base, group, learner in study.ARMS:
        old = original_bundle["arms"][base]
        for variant in VARIANTS + (("depth_calibrated",) if learner == "LR" else ()):
            arm = f"{base}_{variant}"
            start = _clock()
            entry = {"feature_group": group, "learner": learner, "variant": variant,
                     "model": old["model"], "calibrators": {}, "operating_points": {},
                     "available_conditions": list(conditions)}
            if variant == "pooled":
                entry["model"] = study._new_model(learner, dep)
                entry["model"].fit(study._matrix(pooled_train, group, names[group], dep["np"]),
                                   [r["label"] for r in pooled_train])
                base_fits += 1
            if variant == "depth_calibrated":
                calibration_raw = _raw(entry, pooled_calibration, names, dep)
                depth_calibrator = _fit_depth_calibrator(calibration_raw, pooled_calibration, dep)
                calibrator_fits += 1
            for condition in conditions:
                cal = development[("calibration", condition)]
                threshold_rows = development[("threshold", condition)]
                if variant in ("frozen", "retuned"):
                    entry["calibrators"][condition] = old["calibrator"]
                elif variant == "depth_calibrated":
                    entry["calibrators"][condition] = depth_calibrator
                else:
                    entry["calibrators"][condition] = _fit_sigmoid(
                        _raw(entry, cal, names, dep), [r["label"] for r in cal], dep)
                    calibrator_fits += 1
                # The frozen control transfers the old original operating point.
                selection_rows = development[("threshold", original)] if variant == "frozen" else threshold_rows
                scores, _ = _scores(entry, selection_rows, condition, names, dep)
                entry["operating_points"][condition] = {}
                for target in targets:
                    point = (dict(old["operating_point"]) if variant == "frozen" and target == main_target
                             else study.select_score_threshold([r["label"] for r in selection_rows], scores, target))
                    point["selection_partition"] = "threshold"
                    point["selection_condition"] = original if variant == "frozen" else condition
                    point["score_kind"] = ("probability" if variant == "frozen" else
                                           "depth_calibration_logit" if variant == "depth_calibrated" else "raw_model_logit")
                    entry["operating_points"][condition][_target_key(target)] = point
            fits[arm] = _elapsed(start) | {"base_model_fitted": variant == "pooled",
                "calibrator_fits": 2 if variant in ("pooled", "recalibrated") else int(variant == "depth_calibrated")}
            arms[arm] = entry
    arms["B0"] = {"kind": "threshold_grid", "feature_group": "f0", "variant": "grid",
                   "available_conditions": list(conditions), "operating_points": {
                       condition: {_target_key(target): study.select_b0(development[("threshold", condition)], target)
                                   for target in targets} for condition in conditions}}
    caller = {"status": "not_supplied", "reduced_depth_evaluated": False}
    if "B1" in original_bundle.get("arms", {}):
        old = original_bundle["arms"]["B1"]
        entry = {"kind": "caller", "variant": "frozen", "available_conditions": [original],
                 "calibrators": {original: old["calibrator"]}, "operating_points": {original: {}}}
        threshold_rows = development[("threshold", original)]
        score, _ = _scores(entry, threshold_rows, original, names, dep)
        for target in targets:
            entry["operating_points"][original][_target_key(target)] = (
                dict(old["operating_point"]) if target == main_target else
                study.select_score_threshold([r["label"] for r in threshold_rows], score, target))
        arms["B1"] = entry
        caller = {"status": "trusted_original_only", "available_conditions": [original],
                  "reduced_depth_evaluated": False,
                  "limitation": "No raw-read caller rerun occurs in this saved-feature experiment"}

    # Select both candidates before any test prediction; arm IDs break ties.
    stable_order = sorted(arms)
    selectable = [a for a in arms if a.startswith("M") and arms[a]["variant"] != "frozen"
                  and arms[a]["operating_points"][reduced][_target_key(main_target)]["eligible"]]
    selected = {}
    for group, key in (("f0", "counts"), ("f1", "rich")):
        candidates = [a for a in selectable if arms[a]["feature_group"] == group]
        if not candidates:
            raise ValueError(f"no eligible {group} comparator on reduced-depth threshold role")
        def rank(arm):
            point = arms[arm]["operating_points"][reduced][_target_key(main_target)]
            return point["false_positives"], -point["recall"], arm
        selected[key] = min(candidates, key=rank)
    selection = selected | {"selection_partition": "threshold", "selection_condition": reduced,
                            "criterion": "fewest false positives at target recall, highest attained recall, arm ID ascending",
                            "stable_arm_order": stable_order, "test_outcomes_used": False,
                            "frozen_controls_excluded": True}
    bundle = {"schema": "rosalind-depth-calibration-models-v1", "arms": arms,
              "feature_names": names, "protocol_sha256": study._digest(protocol),
              "selected_comparators": selection, "targets": targets}
    model_bytes = pickle.dumps(bundle, protocol=5)
    model_costs = {arm: {"bytes": len(pickle.dumps(entry, protocol=5)),
                          "feature_count": len(names[entry["feature_group"]]) if "feature_group" in entry else 1}
                   for arm, entry in arms.items()}
    frozen = {"schema": "rosalind-depth-calibration-freeze-v1",
              "model_sha256": hashlib.sha256(model_bytes).hexdigest(),
              "protocol_sha256": study._digest(protocol), "sklearn_version": dep["sklearn"].__version__,
              "feature_names": names, "selected_comparators": selection,
              "operating_points": {arm: entry["operating_points"] for arm, entry in arms.items()},
              "development_rows_sha256": study._digest([r for r in records if r["partition"] != "test"]),
              "stage": "frozen_before_test_scoring", "model_sizes": model_costs,
              "training_policy": "Existing original P0 fits reused; pooled models use train role only; paired depths never split",
              "calibration_policy": "Calibration role only; condition-specific positive-slope sigmoid; depth scaler fitted only on pooled calibration role",
              "threshold_policy": "Threshold role only; condition-specific unless frozen control; all tied scores accepted together",
              "monotonic_calibration_policy": "Use original raw model logits for decisions; calibration changes probabilities only",
              "depth_calibration_features": ["raw_score", "log1p(callable_depth)", "raw_score*log1p(callable_depth)"],
              "test_scope": "Previously inspected HG002 P0 test; descriptive, not untouched validation"}
    with (output_dir / "frozen-models.pkl").open("xb") as handle:
        handle.write(model_bytes)
    _write(output_dir / "frozen-models.json", frozen)
    freeze_sha = hashlib.sha256((output_dir / "frozen-models.json").read_bytes()).hexdigest()

    metrics, target_metrics, cache = {}, {}, {}
    prediction_start = _clock()
    # Per-condition arrays preserve candidate identity for paired comparisons.
    for partition in study.PARTITIONS:
        metrics[partition], target_metrics[partition] = {}, {}
        for condition in conditions:
            rows = [r for r in records if r["partition"] == partition and r["condition"] == condition]
            metrics[partition][condition], target_metrics[partition][condition] = {}, {}
            cache[(partition, condition)] = {"rows": rows, "arms": {}}
            for arm, entry in arms.items():
                if condition not in entry["available_conditions"]:
                    continue
                score, probability = _scores(entry, rows, condition, names, dep)
                decisions = {_target_key(target): _decisions(entry, rows, score, condition, target)
                             for target in targets}
                cache[(partition, condition)]["arms"][arm] = {
                    "score": score, "probability": probability, "decisions": decisions}
                arm_metrics = {_target_key(target): _metrics(rows, score, probability, decisions[_target_key(target)])
                               for target in targets}
                target_metrics[partition][condition][arm] = arm_metrics
                metrics[partition][condition][arm] = arm_metrics[_target_key(main_target)]
            for base, _, _ in study.ARMS:
                a = cache[(partition, condition)]["arms"][f"{base}_retuned"]
                b = cache[(partition, condition)]["arms"][f"{base}_recalibrated"]
                score_equal = a["score"] == b["score"]
                decisions_equal = a["decisions"] == b["decisions"]
                if not score_equal or not decisions_equal:
                    raise AssertionError("positive sigmoid changed ranking or threshold decisions")
                rank_equivalence[f"{partition}/{condition}/{base}"] = {
                    "raw_scores_equal": score_equal, "all_target_decisions_equal": decisions_equal,
                    "recalibration_slope": float(arms[f"{base}_recalibrated"]["calibrators"][condition].coef_[0, 0])}
    with (output_dir / "predictions.jsonl").open("xb") as handle:
        for partition in study.PARTITIONS:
            for condition in conditions:
                group = cache[(partition, condition)]
                for i, row in enumerate(group["rows"]):
                    prediction = {key: row[key] for key in META_KEYS}
                    prediction.update({"label": row["label"], "label_reason": row.get("label_reason"),
                        "arms": {arm: {"score": values["score"][i], "probability": values["probability"][i],
                                        "accepted": values["decisions"][_target_key(main_target)][i],
                                        "accepted_by_target": {key: values[i] for key, values in values["decisions"].items()}}
                                 for arm, values in group["arms"].items()}})
                    handle.write(_canonical(prediction))
    prediction_cost = _elapsed(prediction_start)
    curves, comparisons, depth_diagnostics = {}, {}, {}
    bootstrap_start = _clock()
    for condition in conditions:
        group = cache[("test", condition)]
        known = [i for i, row in enumerate(group["rows"]) if row["label"] is not None]
        rows = [group["rows"][i] for i in known]
        labels = [r["label"] for r in rows]
        curves[condition] = {}
        for arm, values in group["arms"].items():
            if arm == "B0":
                continue
            score = [values["score"][i] for i in known]
            curves[condition][arm] = {_target_key(target): descriptive_matched_recall(labels, score, target)
                                      for target in targets}
        c, r = selected["counts"], selected["rich"]
        depth_diagnostics[condition] = {arm: _depth_bins(group["rows"], group["arms"][arm], main_target)
                                       for arm in (c, r)}
        cm = curves[condition][c][_target_key(main_target)]
        rm = curves[condition][r][_target_key(main_target)]
        if not cm["eligible"] or not rm["eligible"]:
            raise ValueError("test candidates cannot support descriptive matched recall")
        fixed_c, fixed_r = metrics["test"][condition][c], metrics["test"][condition][r]
        bootstrap = _paired_bootstrap(rows,
            [group["arms"][c]["score"][i] for i in known],
            [group["arms"][r]["score"][i] for i in known],
            [group["arms"][c]["decisions"][_target_key(main_target)][i] for i in known],
            [group["arms"][r]["decisions"][_target_key(main_target)][i] for i in known],
            main_target, protocol.get("bootstrap", {}), dep)
        gain = cm["false_positives"] - rm["false_positives"]
        comparisons[condition] = {
            "counts_arm": c, "rich_arm": r, "target_recall": main_target,
            "matched_recall": {"counts": cm, "rich": rm,
                "absolute_fp_reduction_counts_minus_rich": gain,
                "relative_fp_reduction_counts_minus_rich": gain / cm["false_positives"] if cm["false_positives"] else None,
                "recall_rich_minus_counts": rm["recall"] - cm["recall"],
                "test_label_selected_thresholds": True, "deployable_validation": False},
            "fixed_operating_points": {"counts": fixed_c, "rich": fixed_r,
                "recall_rich_minus_counts": fixed_r["recall"] - fixed_c["recall"],
                "absolute_fp_reduction_counts_minus_rich": fixed_c["false_positives"] - fixed_r["false_positives"],
                "thresholds_selected_on": "separate threshold role"},
            "bootstrap": bootstrap}
    bootstrap_cost = _elapsed(bootstrap_start)
    feature_storage = _feature_storage(records, names, output_dir)
    score_timings = _score_costs(arms, selected, cache[("test", reduced)]["rows"], reduced,
        names, main_target, int(protocol.get("timing", {}).get("rounds", 7)), dep, output_dir)
    costs = {"fit_calibrate_threshold_by_arm": fits, "new_base_model_fits": base_fits,
             "new_calibrator_fits": calibrator_fits, "all_predictions_metrics_and_write": prediction_cost,
             "bootstrap": bootstrap_cost, "feature_storage": feature_storage,
             "models": model_costs, "frozen_bundle_bytes": len(model_bytes),
             "selected_inference": score_timings,
             "scope": "Saved evidence only; original extraction, download, native build and fallback caller excluded"}
    gate = evaluate_gain_gate(comparisons[reduced], costs, protocol)
    report = {"schema": "rosalind-depth-calibration-experiment-v1",
              "protocol_sha256": study._digest(protocol), "inputs": inputs,
              "scope": "Exploratory recalibration of paired saved HG002 candidates, no independent test or clinical claim",
              "primary_condition": reduced, "primary_recall_target": main_target,
              "descriptive_recall_targets": list(targets), "metrics": metrics,
              "metrics_by_target": target_metrics, "selected_comparators": selection,
              "descriptive_test_matched_recall_curves": curves,
              "selected_test_measured_depth_bins": depth_diagnostics,
              "comparisons": comparisons, "monotonic_calibration_checks": rank_equivalence,
              "caller_baseline": caller, "costs": costs, "gate": gate,
              "frozen_model": {"model_sha256": frozen["model_sha256"], "manifest_sha256": freeze_sha},
              "limitations": gate["limitations"]}
    report["costs"]["total_engine"] = _elapsed(began)
    _write(output_dir / "report.json", report)
    return report
