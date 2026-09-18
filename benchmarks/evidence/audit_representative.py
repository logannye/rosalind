#!/usr/bin/env python3
"""Audit a finished representative run without changing its frozen measurements.

This is an after-run identity/resource audit, not another timed extraction or a
replacement for `rosalind verify`. Original-source identities are SHA-256 checked
again; receipt verification results remain the separately measured native proof.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import statistics
import time


MAX_METADATA_BYTES = 32 << 20
MAX_ROWS = 100000
MAX_STORAGE_FILES = 100000
HASH = frozenset("0123456789abcdef")


def read_json(path, limit=MAX_METADATA_BYTES):
    with Path(path).open("rb") as source:
        data = source.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"metadata exceeds {limit} bytes: {path}")
    return json.loads(data)


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    after = path.stat()
    stamp = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    if stamp(before) != stamp(after):
        raise ValueError(f"file changed during audit hashing: {path}")
    return dict(sha256=digest.hexdigest(), bytes=after.st_size)


def integer(value, label):
    if isinstance(value, bool):
        raise ValueError(f"invalid integer for {label}")
    parsed = int(value)
    if str(parsed) != str(value) or not 0 <= parsed <= (1 << 64) - 1:
        raise ValueError(f"invalid unsigned integer for {label}")
    return parsed


def operand(argv, flag):
    positions = [index for index, value in enumerate(argv) if value == flag]
    if len(positions) != 1 or positions[0] + 1 >= len(argv):
        raise ValueError(f"expected one {flag} operand")
    return argv[positions[0] + 1]


def local_path(root, value):
    path = (root / value).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"artifact is outside retained run directory: {value}")
    return path


def storage_inventory(directory):
    """Count retained bytes, not physical device I/O; refuse symlink traversal."""
    result = dict(path=str(directory), files=0, logical_bytes=0, allocated_bytes=0)
    if not directory.exists():
        return result
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError(f"invalid retained-storage directory: {directory}")
    pending = [directory]
    entries = 0
    while pending:
        with os.scandir(pending.pop()) as children:
            for child in children:
                entries += 1
                if entries > MAX_STORAGE_FILES:
                    raise ValueError("retained-storage entry limit exceeded")
                if child.is_symlink():
                    raise ValueError(f"symlink in retained storage: {child.path}")
                if child.is_dir(follow_symlinks=False):
                    pending.append(Path(child.path))
                elif child.is_file(follow_symlinks=False):
                    stat = child.stat(follow_symlinks=False)
                    result["files"] += 1
                    result["logical_bytes"] += stat.st_size
                    result["allocated_bytes"] += getattr(stat, "st_blocks", 0) * 512
                else:
                    raise ValueError(f"nonregular retained storage: {child.path}")
    return result


def audit(report_path, manifest_path, binary=None, harness_dir=None):
    report_path = Path(report_path).resolve()
    root = report_path.parent
    result = dict(schema=1, kind="representative-post-run-audit", status="failed",
                  report=str(report_path), issues=[], identities=[], measurements=[],
                  prediction_underestimates=[], execution=[], verification_summaries=[],
                  limitations=[
                      "After-run hashes detect changed final bytes; they do not prove files never changed and reverted between invocations.",
                      "Receipt verification is the retained timed native verification, not independently reimplemented here.",
                      "Microtile counts and record visits corroborate work changes; admitted widths and worker counts do not prove actual worker utilization.",
                      "Retained storage sizes are current logical/allocated file bytes, not physical I/O traffic or peak temporary disk use.",
                      "Prediction underestimates are reported separately; only exceeding an admitted budget is a budget breach."])
    started = time.perf_counter()
    try:
        result["report_identity"] = fingerprint(report_path)
        report = read_json(report_path)
        if report.get("status") == "running":
            raise ValueError("matrix is still running; audit only after it finishes")
        rows, inventory = report["measurements"], report["job_inventory"]
        if not isinstance(rows, list) or not isinstance(inventory, list) or max(len(rows), len(inventory)) > MAX_ROWS:
            raise ValueError("invalid or excessive measurement inventory")
        result["base_status"] = report.get("status")
        if report.get("status") != "passed":
            result["issues"].append("base matrix did not pass; its failure evidence is preserved")
        expected, observed = Counter(job["label"] for job in inventory), Counter(row["label"] for row in rows)
        result["inventory"] = dict(expected=len(inventory), recorded=len(rows),
                                   missing=list((expected - observed).elements()),
                                   unexpected=list((observed - expected).elements()),
                                   duplicate_labels=[key for key, count in observed.items() if count > 1])
        if expected != observed or any(count != 1 for count in expected.values()):
            result["issues"].append("invocation inventory differs from retained measurements")
        environment = report["environment"]
        identities = [("source:" + row["role"], Path(row["path"]), row["sha256"]) for row in report["inputs"]]
        if len(identities) > MAX_ROWS:
            raise ValueError("excessive source identity inventory")
        identities.append(("manifest", Path(manifest_path), environment["manifest_sha256"]))
        native_paths = {Path(row["argv"][0]) for row in rows if row["kind"] != "oracle" and row.get("argv")}
        if binary is not None:
            native_paths = {Path(binary)}
        if not native_paths:
            raise ValueError("no binary path in retained invocations; supply --binary")
        identities.extend(("binary", path, environment["binary_sha256"]) for path in sorted(native_paths))
        oracle_paths = {Path(row["argv"][1]).parent for row in rows if row["kind"] == "oracle" and row.get("argv")}
        if harness_dir is not None:
            oracle_paths = {Path(harness_dir)}
        if not oracle_paths:
            raise ValueError("no oracle source path in retained invocations; supply --harness-dir")
        for directory in sorted(oracle_paths):
            for name, digest in environment["harness_sha256"].items():
                if Path(name).name != name:
                    raise ValueError("invalid harness filename")
                identities.append(("harness:" + name, directory / name, digest))
        seen = {}
        for role, path, expected_hash in identities:
            identity = dict(role=role, path=str(path), expected_sha256=expected_hash)
            try:
                if len(expected_hash) != 64 or not set(expected_hash) <= HASH:
                    raise ValueError("invalid startup SHA-256")
                key = path.resolve()
                if key not in seen:
                    seen[key] = fingerprint(path)
                identity.update(seen[key])
                identity["matches_startup"] = identity["sha256"] == expected_hash
                if not identity["matches_startup"]:
                    raise ValueError("startup identity mismatch")
            except (OSError, ValueError) as error:
                identity["error"] = str(error)
                result["issues"].append(f"{role}: {path}: {error}")
            result["identities"].append(identity)
        result["identity_hashing_seconds"] = time.perf_counter() - started
        workloads = {row["id"]: row for row in report["manifest"]["workloads"]}
        cases = {row["id"]: row for row in report["manifest"]["cases"]}
        execution, verification = defaultdict(list), defaultdict(list)
        for row in rows:
            checked = dict(label=row["label"], kind=row["kind"], status=row.get("status"), issues=[])
            result["measurements"].append(checked)
            if row.get("status") != "completed" or row.get("exit_code") != 0 or not row.get("valid"):
                checked["issues"].append("invocation failed, refused, or did not pass semantic verification")
                result["issues"].append(f"{row['label']}: incomplete/invalid invocation retained")
                continue
            try:
                output = local_path(root, row["artifact"])
                observed_output = fingerprint(output)
                if observed_output["sha256"] != row["output_sha256"] or observed_output["bytes"] != row["output_bytes"]:
                    raise ValueError("output identity differs from retained measurement")
                with output.open("rb") as source:
                    count = sum(block.count(b"\n") for block in iter(lambda: source.read(1 << 20), b"")) - 1
                checked["observed_loci"] = count
                expected_rows = workloads[row["workload"]].get("expected_rows")
                if expected_rows is not None and count != expected_rows:
                    raise ValueError("actual output denominator differs from workload")
                if row["kind"] == "oracle":
                    oracle = row["oracle_measurements"]
                    if integer(oracle["selected_loci"], "oracle rows") != count:
                        raise ValueError("oracle row count differs from actual output")
                    continue
                receipt_path = local_path(root, operand(row["argv"], "--manifest"))
                receipt = read_json(receipt_path)
                values, params = receipt["measurements"], receipt["params"]
                if values != row.get("receipt_measurements"):
                    raise ValueError("receipt measurements changed after timed run")
                budget = integer(row["declared_budget_mib"], "declared budget") * 1048576
                if budget != integer(cases[row["case"]]["budget_mib"], "case budget") * 1048576 or budget != integer(
                        operand(row["argv"], "--memory-budget-mb"), "argv budget") * 1048576:
                    checked["issues"].append("measurement budget differs from manifest case or executed argv")
                receipt_budget = (integer(params["memory_budget_bytes"], "receipt budget") if "memory_budget_bytes" in params
                                  else integer(params["memory_budget_mb"], "receipt budget") * 1048576)
                checked.update(declared_budget_bytes=budget, receipt_budget_bytes=receipt_budget,
                               observed_peak_rss_bytes=integer(row["peak_rss_bytes"], "process RSS"),
                               receipt_peak_rss_bytes=integer(values["peak_rss_bytes"], "receipt RSS"))
                if not budget or budget != receipt_budget:
                    checked["issues"].append("receipt budget differs from requested budget")
                if params.get("run_status") != "completed":
                    checked["issues"].append("receipt does not declare completed status")
                verdict = values.get("contract_verdict", params.get("contract_verdict"))
                checked["receipt_contract_verdict"] = verdict
                if verdict is not None and verdict != "within":
                    checked["issues"].append("receipt contract verdict is not within")
                checked["budget_breach"] = max(checked["observed_peak_rss_bytes"], checked["receipt_peak_rss_bytes"]) > budget
                if checked["budget_breach"]:
                    checked["issues"].append("observed process or receipt RSS exceeds declared budget")
                if integer(values["execution.emitted_loci"], "receipt rows") != count:
                    checked["issues"].append("receipt denominator differs from actual output")
                if "predicted_peak_rss_bytes" in values:
                    predicted = integer(values["predicted_peak_rss_bytes"], "prediction")
                    checked["predicted_peak_rss_bytes"] = predicted
                    if max(checked["observed_peak_rss_bytes"], checked["receipt_peak_rss_bytes"]) > predicted:
                        result["prediction_underestimates"].append(dict(label=row["label"], predicted_bytes=predicted,
                            observed_bytes=checked["observed_peak_rss_bytes"], receipt_bytes=checked["receipt_peak_rss_bytes"],
                            budget_bytes=budget, budget_breach=checked["budget_breach"]))
                verified = row["verification"]
                if verified["exit_code"] != 0:
                    checked["issues"].append("retained native verification did not pass")
                verification_path = local_path(root, f"raw/{row['label']}.verify.json")
                if fingerprint(verification_path)["sha256"] != verified["stdout_sha256"]:
                    checked["issues"].append("retained verification output changed")
                verification[(row["workload"], row["kind"], row.get("case"), row.get("phase"))].append(verified["wall_seconds"])
                if row["kind"] == "rosalind":
                    execution[row["workload"]].append(dict(label=row["label"], phase=row.get("phase"),
                        case=row.get("case"), repeat=row.get("repeat"), loci=count,
                        admitted_tile_bases=integer(values["execution.microtile_bases"], "tile width"),
                        admitted_workers=integer(values.get("execution.worker_count", 1), "worker count"),
                        actual_microtiles=integer(values["execution.microtiles"], "microtiles"),
                        alignment_record_visits=integer(values["execution.record_visits"], "record visits")))
            except (OSError, ValueError, KeyError, TypeError) as error:
                checked["issues"].append(str(error))
            result["issues"].extend(f"{row['label']}: {issue}" for issue in checked["issues"])
        for workload, observations in sorted(execution.items()):
            computed = [row for row in observations if row["phase"] != "resumed"]
            signatures = {(row["actual_microtiles"], row["alignment_record_visits"]) for row in computed}
            result["execution"].append(dict(workload=workload, observations=observations,
                distinct_observed_work_signatures=len(signatures), observed_execution_changed=len(signatures) > 1,
                denominator_consistent=len({row["loci"] for row in observations}) == 1))
            if workloads[workload].get("gates", {}).get("min_effective_tiles", 1) > 1 and len(signatures) < 2:
                result["issues"].append(f"{workload}: distinct admitted widths have no observed work difference")
        for key, times in sorted(verification.items()):
            result["verification_summaries"].append(dict(zip(("workload", "kind", "case", "phase"), key),
                runs=len(times), wall_seconds=dict(median=statistics.median(times), minimum=min(times), maximum=max(times))))
        result["retained_storage"] = [storage_inventory(root / name) for name in ("cache", "artifacts", "raw")]
        if fingerprint(report_path) != result["report_identity"]:
            result["issues"].append("base report changed during audit")
        result["status"] = "passed" if not result["issues"] else "failed"
    except (OSError, ValueError, KeyError, TypeError) as error:
        result["issues"].append(f"audit could not complete: {error}")
    result["elapsed_seconds"] = time.perf_counter() - started
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path, help="original workload manifest, checked against startup hash")
    parser.add_argument("--output", required=True, type=Path, help="new supplementary JSON; existing files are preserved")
    parser.add_argument("--binary", type=Path, help="content-equivalent relocated binary")
    parser.add_argument("--harness-dir", type=Path, help="content-equivalent relocated frozen harness sources")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("supplementary audit destination already exists")
    result = audit(args.report, args.manifest, args.binary, args.harness_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(f"{result['status']}: {args.output}")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
