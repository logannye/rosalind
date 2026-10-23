"""Identity-bound hosted threshold freezing and untouched HG001 evaluation."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import pickle
import platform
import signal
import subprocess
import sys
import traceback
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "evidence-sufficiency"))
sys.path.insert(0, str(ROOT / "python"))


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def write(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def read(path):
    return json.loads(Path(path).read_text())


def execution_deadline(protocol_path, window_path):
    """Bind an explicit operational resumption to unchanged scientific bytes."""
    protocol = read(protocol_path)
    original = datetime.fromisoformat(protocol["deadline_utc"].replace("Z", "+00:00"))
    if not window_path.exists():
        return original
    window = read(window_path)
    start = datetime.fromisoformat(window["started_utc"].replace("Z", "+00:00"))
    end = datetime.fromisoformat(window["deadline_utc"].replace("Z", "+00:00"))
    if (window["schema"] != "rosalind-threshold-resumption-v1"
            or window["frozen_protocol_sha256"] != digest(protocol_path)
            or window["original_deadline_utc"] != protocol["deadline_utc"]
            or not original < start < end
            or (end - start).total_seconds() > 10800
            or window["previous_hosted_attempts"] != 2
            or window["additional_hosted_attempts"] != 1
            or not 0 < window["maximum_additional_runner_seconds"] <= 10800
            or window["previous_runner_seconds"] < 0
            or window["previous_runner_seconds"] + window["maximum_additional_runner_seconds"] > 21600
            or window["maximum_total_runner_seconds"] != 21600):
        raise ValueError("Operational resumption differs from the authorized remaining budget or frozen protocol")
    if datetime.now(timezone.utc) < start:
        raise ValueError("Operational resumption window has not started")
    return end


def provenance_module():
    spec = importlib.util.spec_from_file_location("saved_hg002_provenance", HERE.parent / "depth-calibration/run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_model_archive(archive, destination, source):
    """Extract only pinned model members; deserialize only after full validation."""
    if digest(archive) != source["artifact_sha256"]:
        raise ValueError("Trusted model ZIP digest differs")
    if destination.exists():
        for relative, expected in source["files"].items():
            path = destination / relative
            if path.stat().st_size != expected["bytes"] or digest(path) != expected["sha256"]:
                raise ValueError(f"Existing trusted model member differs: {relative}")
    else:
        destination.mkdir(parents=True)
        with zipfile.ZipFile(archive) as bundle:
            names = bundle.namelist()
            for relative, expected in source["files"].items():
                wanted = source["prefix"] + relative
                if names.count(wanted) != 1 or bundle.getinfo(wanted).file_size != expected["bytes"]:
                    raise ValueError(f"Missing, duplicate or oversized trusted model member: {relative}")
                path = destination / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                with bundle.open(wanted) as src, path.open("xb") as dst:
                    for chunk in iter(lambda: src.read(1 << 20), b""):
                        dst.write(chunk)
                if digest(path) != expected["sha256"]:
                    raise ValueError(f"Trusted model member checksum differs: {relative}")
    identity = read(destination / "identity.json")
    complete = read(destination / "completion.json")
    manifest = read(destination / "run-manifest.json")
    if (identity["git_revision"] != source["commit"]
            or int(identity["github_run_id"]) != source["run_id"]
            or complete.get("status") != "completed"
            or complete.get("input_identity_rechecked") is not True
            or complete.get("resource_checks_passed") is not True
            or complete["report_sha256"] != digest(destination / "report.json")
            or complete["run_manifest_sha256"] != digest(destination / "run-manifest.json")):
        raise ValueError("Saved model execution/completion identity mismatch")
    for relative, expected in source["files"].items():
        if relative not in ("completion.json", "run-manifest.json") and manifest["files"].get(relative) != expected:
            raise ValueError("Trusted model inventory differs")
    frozen = read(destination / "frozen-models.json")
    if frozen["model_sha256"] != digest(destination / "frozen-models.pkl"):
        raise ValueError("Trusted model freeze differs")
    if read(destination / "report.json")["gate"]["status"] != "inconclusive_validation_screen":
        raise ValueError("Historical calibration outcome differs")
    dependencies = {name: importlib.metadata.version(name) for name in identity["dependencies"]}
    if dependencies != identity["dependencies"]:
        raise ValueError("Installed wheels differ from trusted model environment")
    # Feature and estimator semantics must match the source that created the models.
    for relative in ("research/evidence-sufficiency/study.py", "research/evidence-sufficiency/features.py"):
        expected = identity["source_sha256"].get(relative)
        if expected is not None and digest(ROOT / relative) != expected:
            raise ValueError(f"Model semantics changed: {relative}")
    with (destination / "frozen-models.pkl").open("rb") as handle:
        bundle = pickle.load(handle)
    if bundle.get("schema") != "rosalind-depth-calibration-models-v1":
        raise ValueError("Unsupported trusted model schema")
    return bundle, identity


def source_snapshot(binary, protocol_path):
    paths = sorted(HERE.glob("*.py")) + sorted((HERE.parent / "evidence-sufficiency").glob("*.py"))
    paths += sorted((ROOT / "python/rosalind").glob("*.py"))
    paths += [protocol_path, ROOT / ".github/workflows/threshold-transfer.yml"]
    if (HERE / "execution-window.json").exists():
        paths.append(HERE / "execution-window.json")
    return {
        "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "github_run_id": os.environ["GITHUB_RUN_ID"],
        "github_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "platform": platform.platform(), "python": sys.version,
        "native_sha256": digest(binary),
        "native_version": subprocess.check_output([str(binary), "--version"], text=True).strip(),
        "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in paths},
    }


def recheck(identity, binary):
    if digest(binary) != identity["native_sha256"]:
        raise ValueError("Native binary changed")
    for relative, expected in identity["source_sha256"].items():
        if digest(ROOT / relative) != expected:
            raise ValueError(f"Analysis source changed: {relative}")


def load_sources(args, protocol):
    provenance = provenance_module()
    for kind, source in (("source", protocol["source"]), ("model", protocol["model_source"])):
        provenance.verify_provider(read(args.inputs / f"{kind}-run.json"),
                                   read(args.inputs / f"{kind}-artifact.json"), source)
    model_dir = args.inputs / "verified-models"
    bundle, identity = verify_model_archive(args.inputs / "model.zip", model_dir, protocol["model_source"])
    return provenance, bundle, identity


def freeze(args, protocol, binary):
    from thresholds import freeze_policies
    args.output.mkdir(parents=True, exist_ok=False)
    identity = source_snapshot(binary, args.protocol)
    write(args.output / "identity.json", identity)
    write(args.output / "protocol.json", protocol)
    if (HERE / "execution-window.json").exists():
        write(args.output / "execution-window.json", read(HERE / "execution-window.json"))
    provenance, bundle, model_identity = load_sources(args, protocol)
    saved = args.inputs / "verified-hg002"
    original_identity = provenance.extract_verified(args.inputs / "source.zip", saved, protocol["source"])
    for relative in ("research/evidence-sufficiency/features.py", "research/evidence-sufficiency/study.py"):
        if digest(ROOT / relative) != original_identity["source_sha256"][relative]:
            raise ValueError(f"Feature/model semantics differ from original training source: {relative}")
    with (saved / "records.jsonl").open() as handle:
        records = [json.loads(line) for line in handle if line.strip()]
    policy = freeze_policies(records, bundle, args.output / "policies", protocol)
    write(args.output / "frozen-policy.json", policy)
    recheck(identity, binary)
    write(args.output / "freeze-completion.json", {
        "status": "frozen_before_fresh_acquisition", "frozen_policy_sha256": digest(args.output / "frozen-policy.json"),
        "protocol_sha256": digest(args.output / "protocol.json"),
        "model_pickle_sha256": protocol["model_source"]["files"]["frozen-models.pkl"]["sha256"],
        "model_source_revision": model_identity["git_revision"],
        "fresh_test_scored": False, "base_model_fits": 0,
    })
    print(json.dumps({"phase": "freeze_completed", "selected": policy.get("selected", policy.get("selected_policies"))}))


def fresh(args, protocol, binary):
    from prepare_fresh import prepare_fresh, verify_input_snapshot
    from thresholds import evaluate_frozen
    from run_pilot import Telemetry
    from costs import measure_costs
    freeze_receipt = read(args.output / "freeze-completion.json")
    if (freeze_receipt["status"] != "frozen_before_fresh_acquisition"
            or freeze_receipt["frozen_policy_sha256"] != digest(args.output / "frozen-policy.json")
            or freeze_receipt["protocol_sha256"] != digest(args.protocol)
            or read(args.output / "protocol.json") != protocol):
        raise ValueError("Frozen threshold/protocol identity differs")
    identity = read(args.output / "identity.json")
    recheck(identity, binary)
    _, bundle, _ = load_sources(args, protocol)
    policy = read(args.output / "frozen-policy.json")
    out = args.output / "fresh"
    if out.exists():
        raise FileExistsError("Fresh run exists; never silently rescore or overwrite it")
    with Telemetry(args.output, protocol) as telemetry:
        with telemetry.phase("fresh_acquisition_extraction_and_oracles"):
            data = prepare_fresh(protocol, out, binary)
        write(args.output / "fresh-evaluation-started.json", {
            "fresh_test_consumed": True, "frozen_policy_sha256": freeze_receipt["frozen_policy_sha256"],
            "rule": "No threshold, model or region changes after this marker; execution-only repairs preserve frozen bytes",
        })
        with telemetry.phase("frozen_fresh_evaluation"):
            report = evaluate_frozen(data["records"], bundle, policy, args.output / "evaluation", protocol)
        with telemetry.phase("same_candidate_cost_measurements"):
            costs = measure_costs(data["records"], bundle, policy, args.output / "costs", protocol)
        report["fresh_data"] = {key: value for key, value in data.items() if key not in ("records", "state")}
        report["evaluation_costs"] = report.get("costs", {})
        report["costs"] = costs
        report["scope"] = protocol["scope"]
        report["cost_scope"] = protocol["practical_scope"]
        report["p1_authorized"] = False
        report["selective_escalation_authorized"] = False
        report["product_promotion_authorized"] = False
        report["full_user_workflow_advantage_established"] = False
        recovery = {}
        for condition, regions in data["candidate_recall"].items():
            denominator = sum(region["truth_snvs"] for region in regions.values())
            proposed = sum(region["proposed_truth_snvs"] for region in regions.values())
            arms = {}
            for label in ("counts", "rich"):
                values = report["comparisons"][condition][label]
                if values["true_positives"] + values["false_negatives"] != proposed:
                    raise ValueError("Candidate-generation and filtering truth denominators differ")
                arms[label] = {"retained_truth_snvs": values["true_positives"],
                               "missed_by_candidate_generation": denominator-proposed,
                               "filtered_true_candidates": values["false_negatives"],
                               "truth_retention_fraction": values["true_positives"]/denominator if denominator else None}
            caller = data["metrics"]["condition_costs"][condition]["caller_unthresholded_snv_calls"]
            recovery[condition] = {"assessed_truth_snvs": denominator, "proposed_truth_snvs": proposed,
                "candidate_generation_recall": proposed/denominator if denominator else None,
                "selected_policies": arms, "caller_emitted_truth_snvs": caller["true"],
                "caller_truth_recovery": caller["true"]/denominator if denominator else None,
                "scope": "Policy truth retention includes unresolved pass-throughs; it is not final called-variant recall"}
        report["end_to_end_truth_retention"] = recovery
        scientific_gate = report["gate"]
        correctness = data["metrics"].get("correctness", {})
        external_checks = {
            "fresh_independent_counts": correctness.get("independent_counts") is True,
            "fresh_native_parquet_features": correctness.get("native_parquet_features") is True,
            "standalone_counts_features_equal": correctness.get("standalone_counts_f0_equal") is True,
            "raw_unavailable_replay": correctness.get("saved_only_replay") is True,
            "standalone_columnar_policy_parity": costs.get("native_parquet_prediction_parity") is True,
            "all_declared_depth_costs_present": set(protocol["required_fresh_conditions"]).issubset(data["metrics"].get("condition_costs", {})),
            "saved_feature_costs_recorded": bool(costs.get("storage") and costs.get("measurements")),
            "fresh_inputs_unchanged": verify_input_snapshot(data["identities"]["inputs"]).get("passed") is True,
        }
        scientific_gate["external_checks"] = external_checks
        scientific_gate["nominate_fresh_workflow_validation"] = (
            scientific_gate.get("predictive_and_workload_screen_passed") is True and all(external_checks.values()))
        if not all(external_checks.values()):
            scientific_gate["status"] = "incomplete_correctness_or_cost_evidence"
        recheck(identity, binary)
        if digest(args.output / "frozen-policy.json") != freeze_receipt["frozen_policy_sha256"]:
            raise ValueError("Frozen policy changed during fresh evaluation")
        write(args.output / "report.json", report)
        write(args.output / "input-recheck.json", {"source_and_frozen_policy_unchanged": True})
    # Raw acquisition scratch is excluded; data adapter records its identities.
    files = {str(p.relative_to(args.output)): {"bytes": p.stat().st_size, "sha256": digest(p)}
             for p in sorted(args.output.rglob("*")) if p.is_file()
             and "inputs" not in p.relative_to(args.output).parts}
    total = sum(entry["bytes"] for entry in files.values())
    if total > protocol["limits"]["retained_evidence_bytes"]:
        raise ValueError("Retained evidence budget exceeded")
    write(args.output / "run-manifest.json", {"files": files, "retained_bytes": total})
    write(args.output / "completion.json", {
        "status": "completed", "report_sha256": digest(args.output / "report.json"),
        "run_manifest_sha256": digest(args.output / "run-manifest.json"),
        "frozen_policy_sha256": freeze_receipt["frozen_policy_sha256"],
        "fresh_test_consumed": True, "input_identity_rechecked": True,
        "p1_authorized": False, "selective_escalation_authorized": False,
    })
    print(json.dumps({"phase": "fresh_completed", "gate": report.get("gate")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("freeze", "fresh"), required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if os.environ.get("GITHUB_ACTIONS") != "true" or platform.system() != "Linux":
        raise RuntimeError("Scientific execution requires the authorized hosted Linux workflow")
    protocol = read(args.protocol)
    remaining = (execution_deadline(args.protocol, HERE / "execution-window.json")
                 - datetime.now(timezone.utc)).total_seconds()
    if remaining <= 0:
        raise TimeoutError("Overnight execution deadline passed; no new study phase")
    binary = args.binary.resolve(strict=True)
    def interrupted(signum, _frame):
        raise TimeoutError(f"Study interrupted by signal {signum}")
    signal.signal(signal.SIGALRM, interrupted)
    signal.signal(signal.SIGTERM, interrupted)
    signal.alarm(max(1, min(protocol["limits"]["analysis_seconds"], int(remaining))))
    try:
        (freeze if args.phase == "freeze" else fresh)(args, protocol, binary)
    except BaseException as exc:
        args.output.mkdir(parents=True, exist_ok=True)
        path = args.output / f"failure-{args.phase}.json"
        if not path.exists():
            write(path, {"status": "invalid", "phase": args.phase, "reason": repr(exc), "traceback": traceback.format_exc()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
