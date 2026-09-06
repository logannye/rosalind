#!/usr/bin/env python3
"""Run an explicit, content-locked evidence matrix with a bounded pysam oracle."""
import argparse
import json
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys
import time

from run import compare, run_measured, sha256

MAX_MANIFEST_BYTES = 4 << 20
FILE_ROLES = ("alignments", "alignment_index", "reference", "reference_fai", "selection")
IDENTIFIER = re.compile(r"[a-z0-9][a-z0-9-]{0,63}\Z")


def file_stamp(path):
    stat = Path(path).stat()
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def file_identity(path):
    path = Path(path)
    before = file_stamp(path)
    digest = sha256(path)
    after = file_stamp(path)
    if before != after:
        raise ValueError(f"input changed during harness fingerprinting: {path}")
    return dict(path=str(path), sha256=digest, bytes=after[2], snapshot=after)


def captured_identities(report):
    return {row["path"]: row for row in report.get("inputs", []) + report.get("execution_inputs", [])}


def check_snapshots(report):
    for path, identity in captured_identities(report).items():
        if file_stamp(path) != identity["snapshot"]:
            raise ValueError(f"immutable benchmark input changed: {path}")


def final_identity_check(report):
    started, checks = time.perf_counter(), []
    for path, expected in captured_identities(report).items():
        check = dict(path=path, expected_sha256=expected["sha256"], valid=False)
        try:
            observed = file_identity(path)
            check.update(observed)
            check["valid"] = observed["sha256"] == expected["sha256"] and observed["snapshot"] == expected["snapshot"]
            if not check["valid"]:
                check["error"] = "input content or metadata changed after startup"
        except (OSError, ValueError) as error:
            check["error"] = str(error)
        checks.append(check)
    report["final_input_checks"] = checks
    report["harness_final_hashing_seconds"] = time.perf_counter() - started
    return bool(checks) and all(check["valid"] for check in checks)


def read_json(path, limit=MAX_MANIFEST_BYTES):
    with Path(path).open("rb") as source:
        data = source.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"metadata exceeds {limit} bytes: {path}")
    return json.loads(data)


def positive(value, label, maximum=1 << 31):
    if type(value) is not int or not 1 <= value <= maximum:
        raise ValueError(f"{label} must be a positive integer <= {maximum}")
    return value


def identifier(value):
    if not isinstance(value, str) or not IDENTIFIER.fullmatch(value):
        raise ValueError(f"invalid case/workload identifier: {value!r}")
    return value


def load_manifest(path, smoke=False):
    path = Path(path).resolve()
    manifest = read_json(path)
    if manifest.get("schema") != 1:
        raise ValueError("unsupported representative workload schema")
    repeats = positive(manifest.get("repeats", 3), "repeats", 100)
    if repeats < 3 and not smoke:
        raise ValueError("representative measurements require at least three repeats; --smoke is diagnostic only")
    manifest["repeats"] = repeats
    if type(manifest.get("seed", 17)) is not int:
        raise ValueError("seed must be an integer")
    oracle = manifest.setdefault("oracle", {})
    for name, default, maximum in [("tile_bases", 1024, 16384), ("max_read_len", 250, 1000000),
                                   ("max_record_bytes", 1048576, 1 << 30)]:
        oracle[name] = positive(oracle.get(name, default), "oracle." + name, maximum)
    cases = manifest.get("cases", [])
    workloads = manifest.get("workloads", [])
    if not 1 <= len(cases) <= 256 or not 1 <= len(workloads) <= 32:
        raise ValueError("manifest requires 1..256 cases and 1..32 workloads")
    case_ids = set()
    for case in cases:
        name = identifier(case["id"])
        if name in case_ids:
            raise ValueError("duplicate case id")
        case_ids.add(name)
        positive(case["budget_mib"], "budget_mib")
        positive(case["tile_bases"], "tile_bases", 16384)
        positive(case["workers"], "workers", 64)
        if case.get("cache", "none") not in ("none", "cold-resume", "cold-resume-dataset"):
            raise ValueError("case cache must be none, cold-resume or cold-resume-dataset")
    workload_ids = set()
    for workload in workloads:
        name = identifier(workload["id"])
        if name in workload_ids:
            raise ValueError("duplicate workload id")
        workload_ids.add(name)
        identifier(workload.get("equivalence_group", name))
        selected = workload.get("cases", sorted(case_ids))
        if not selected or len(set(selected)) != len(selected) or not set(selected) <= case_ids:
            raise ValueError("workload cases must name unique declared cases")
        workload["cases"] = selected
        if workload.get("sample") and workload.get("pool_samples"):
            raise ValueError("sample and pool_samples are mutually exclusive")
        for role in FILE_ROLES:
            file = workload[role]
            if not isinstance(file.get("sha256"), str) or not re.fullmatch("[0-9a-f]{64}", file["sha256"]):
                raise ValueError(f"{name}.{role} requires a lowercase SHA-256")
            file["path"] = str((path.parent / file["path"]).resolve())
        if workload["selection"].get("kind") not in ("sites", "regions"):
            raise ValueError("selection.kind must be sites or regions")
        if "expected_rows" in workload:
            positive(workload["expected_rows"], "expected_rows")
        gates = workload.setdefault("gates", {})
        gates.setdefault("min_admitted_budgets", 3)
        gates.setdefault("min_effective_tiles", 1)
        gates.setdefault("workers", [1])
        positive(gates["min_admitted_budgets"], "min_admitted_budgets", 256)
        positive(gates["min_effective_tiles"], "min_effective_tiles", 16384)
        for worker in gates["workers"]:
            positive(worker, "gate worker count", 64)
    return manifest


def job_groups(manifest):
    """Randomize independent groups; preserve cold -> resumed dependencies."""
    cases = {case["id"]: case for case in manifest["cases"]}
    result = []
    for repeat in range(1, manifest["repeats"] + 1):
        groups = []
        for workload in manifest["workloads"]:
            groups.append([dict(workload=workload["id"], kind="oracle", repeat=repeat,
                                label=f"{workload['id']}-oracle-r{repeat}")])
            for case_id in workload["cases"]:
                case = cases[case_id]
                base = f"{workload['id']}-{case_id}-r{repeat}"
                phases = {"cold-resume": ("cold", "resumed"),
                          "cold-resume-dataset": ("cold", "resumed", "persisted")}.get(case.get("cache"), ("fresh",))
                groups.append([dict(workload=workload["id"], case=case_id,
                                    kind="dataset" if phase == "persisted" else "rosalind",
                                    repeat=repeat, phase=phase, cache_key=base,
                                    label=base + ("-" + phase if phase != "fresh" else ""))
                               for phase in phases])
        random.Random(manifest.get("seed", 17) + repeat).shuffle(groups)
        result.extend(groups)
    return result


def verify_inputs(manifest):
    identities, seen = [], {}
    started = time.perf_counter()
    for workload in manifest["workloads"]:
        for role in FILE_ROLES:
            source = workload[role]
            path = Path(source["path"])
            if path not in seen:
                seen[path] = file_identity(path)
            if seen[path]["sha256"] != source["sha256"]:
                raise ValueError(f"content fingerprint mismatch: {path}")
            identities.append(dict(workload=workload["id"], role=role, **seen[path]))
    return identities, time.perf_counter() - started


def native_command(binary, workload, case, oracle, output, receipt, cache=None, resume=False):
    argv = [str(binary), "analyze", "evidence", "--fields", "all", "--format", "tsv",
            "--mapq-threshold", "20", "--base-quality-threshold", "20"]
    for role in FILE_ROLES[:-1]:
        argv += ["--" + role.replace("_", "-"), workload[role]["path"]]
    argv += ["--" + workload["selection"]["kind"], workload["selection"]["path"],
             "--memory-budget-mb", str(case["budget_mib"]), "--tile-bases", str(case["tile_bases"]),
             "--workers", str(case["workers"]), "--max-read-len", str(oracle["max_read_len"]),
             "--max-record-bytes", str(oracle["max_record_bytes"]),
             "--output", str(output), "--manifest", str(receipt)]
    if workload.get("sample"):
        argv += ["--sample", workload["sample"]]
    if workload.get("pool_samples"):
        argv += ["--pool-samples"]
    if cache is not None:
        argv += ["--cache-dir", str(cache)]
    if resume:
        argv += ["--resume"]
    return argv


def oracle_command(workload, oracle, stats):
    argv = [sys.executable, str(Path(__file__).with_name("streaming_pysam.py")),
            "--reference", workload["reference"]["path"], "--reference-fai", workload["reference_fai"]["path"],
            "--alignments", workload["alignments"]["path"],
            "--alignment-index", workload["alignment_index"]["path"],
            "--selection", workload["selection"]["path"], "--selection-kind", workload["selection"]["kind"],
            "--stats", str(stats)]
    for key, value in oracle.items():
        if key in ("tile_bases", "max_read_len", "max_record_bytes"):
            argv += ["--" + key.replace("_", "-"), str(value)]
    if workload.get("sample"):
        argv += ["--sample", workload["sample"]]
    if workload.get("pool_samples"):
        argv += ["--pool-samples"]
    return argv


def dataset_command(binary, workload, case, parent, output, receipt):
    return [str(binary), "dataset", "extract", "--dataset", str(parent),
            "--fields", "all", "--format", "tsv", "--" + workload["selection"]["kind"],
            workload["selection"]["path"], "--memory-budget-mb", str(case["budget_mib"]),
            "--output", str(output), "--manifest", str(receipt)]


def extra_time_metrics(path):
    text = Path(path).read_text()
    metrics = {}
    patterns = {
        "user_cpu_seconds": [r"User time \(seconds\):\s*([\d.]+)", r"([\d.]+)\s+user"],
        "system_cpu_seconds": [r"System time \(seconds\):\s*([\d.]+)", r"([\d.]+)\s+sys"],
        "filesystem_input_operations": [r"File system inputs:\s*(\d+)", r"(\d+)\s+block input operations"],
        "filesystem_output_operations": [r"File system outputs:\s*(\d+)", r"(\d+)\s+block output operations"],
    }
    for key, candidates in patterns.items():
        for pattern in candidates:
            match = re.search(pattern, text)
            if match:
                metrics[key] = float(match[1]) if key.endswith("seconds") else int(match[1])
                break
    return metrics


def assess_gates(manifest, measurements):
    results = {}
    for workload in manifest["workloads"]:
        rows = [row for row in measurements if row["workload"] == workload["id"] and row["kind"] == "rosalind"]
        successful = [row for row in rows if row.get("valid")]
        effective = [row["effective"] for row in successful]
        budgets = sorted({row["declared_budget_mib"] for row in successful})
        workers = sorted({value["workers"] for value in effective})
        tiles = sorted({value["tile_bases"] for value in effective})
        work = sorted({(row["effective"]["microtiles"], row["effective"]["record_visits"])
                       for row in successful if row.get("phase") != "resumed"})
        gates = workload["gates"]
        checks = {
            "all_requested_runs_completed_verified_and_equal": bool(rows) and len(rows) == len(successful),
            "three_or_more_repeats": manifest["repeats"] >= 3,
            "admitted_budgets": len(budgets) >= gates["min_admitted_budgets"],
            "effective_tiles": len(tiles) >= gates["min_effective_tiles"],
            "observed_work_changes": gates["min_effective_tiles"] < 2 or len(work) >= 2,
            "effective_workers": set(gates["workers"]) <= set(workers),
            "resume_performs_no_indexed_extraction": all(
                row.get("effective", {}).get("record_visits") == 0
                and row.get("effective", {}).get("computed_partitions") == 0
                and row.get("effective", {}).get("reused_partitions", 0) > 0
                for row in rows if row.get("phase") == "resumed"),
        }
        results[workload["id"]] = dict(passed=all(checks.values()), checks=checks,
                                       observed_budgets_mib=budgets, observed_workers=workers,
                                       observed_tile_bases=tiles, observed_work_signatures=work)
    return results


def summaries(measurements):
    grouped = {}
    for row in measurements:
        key = (row["workload"], row["kind"], row.get("case", "oracle"), row.get("phase", "oracle"))
        grouped.setdefault(key, []).append(row)
    result = []
    for key, rows in grouped.items():
        summary = dict(zip(("workload", "kind", "case", "phase"), key))
        summary.update(requested_repeats=len(rows), completed_repeats=sum(row.get("valid", False) for row in rows))
        for metric in ("wall_seconds", "peak_rss_bytes", "output_bytes", "user_cpu_seconds",
                       "system_cpu_seconds", "filesystem_input_operations", "filesystem_output_operations",
                       "indexed_alignment_record_visits", "decoder_bytes", "cram_validation_records",
                       "cram_validation_bases", "cram_validation_wall_seconds"):
            values = [row[metric] for row in rows if row.get("valid") and row.get(metric) is not None]
            if values:
                summary[metric] = dict(median=statistics.median(values), minimum=min(values), maximum=max(values))
        result.append(summary)
    return result


def write_report(root, report):
    temporary = root / "report.json.tmp"
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(root / "report.json")


def unsigned(value, label):
    parsed = int(value)
    if isinstance(value, bool) or str(parsed) != str(value) or not 0 <= parsed <= (1 << 64) - 1:
        raise ValueError(f"invalid unsigned integer: {label}")
    return parsed


def validate_resources(row, claim, case):
    values, params = claim["measurements"], claim["params"]
    budget = case["budget_mib"] * 1048576
    receipt_budget = (unsigned(params["memory_budget_bytes"], "receipt budget") if "memory_budget_bytes" in params
                      else unsigned(params["memory_budget_mb"], "receipt budget") * 1048576)
    measured_peak = unsigned(row["peak_rss_bytes"], "observed RSS")
    receipt_peak = unsigned(values["peak_rss_bytes"], "receipt RSS")
    verdict = values.get("contract_verdict", params.get("contract_verdict"))
    checks = dict(receipt_budget_matches=receipt_budget == budget,
                  receipt_completed=params.get("run_status") == "completed",
                  observed_rss_within_budget=measured_peak <= budget,
                  receipt_rss_within_budget=receipt_peak <= budget,
                  contract_verdict_within=(verdict == "within" if row["kind"] == "rosalind" else verdict in (None, "within")))
    row["resource_validation"] = dict(passed=all(checks.values()), checks=checks,
        declared_budget_bytes=budget, receipt_budget_bytes=receipt_budget,
        measured_peak_rss_bytes=measured_peak, receipt_peak_rss_bytes=receipt_peak,
        contract_verdict=verdict)
    if "predicted_peak_rss_bytes" in values:
        predicted = unsigned(values["predicted_peak_rss_bytes"], "predicted RSS")
        row["prediction_validation"] = dict(predicted_peak_rss_bytes=predicted,
            observed_peak_rss_bytes=max(measured_peak, receipt_peak),
            underestimated=max(measured_peak, receipt_peak) > predicted)


def measure_job(manifest, binary, root, job, workload, cases, row, measure):
    """Keep the process result even if metadata parsing or verification fails."""
    label = job["label"]
    output = root / "artifacts" / (label + ".tsv")
    receipt = root / "artifacts" / (label + ".tsv.manifest.json")
    row["artifact"] = str(output.relative_to(root))
    began = time.perf_counter()
    try:
        if job["kind"] == "oracle":
            stats_path = root / "artifacts" / (label + ".oracle.json")
            argv = oracle_command(workload, manifest["oracle"], stats_path)
            row["argv"] = argv
            row.update(measure(root, label, argv, output))
            if stats_path.is_file():
                row["oracle_measurements"] = read_json(stats_path)
                if row["exit_code"] == 0:
                    for key in ("selected_loci", "record_visits", "windows"):
                        unsigned(row["oracle_measurements"][key], "oracle " + key)
            elif row["exit_code"] == 0:
                raise ValueError("successful oracle did not produce statistics")
        else:
            case = cases[job["case"]]
            row["declared_budget_mib"] = case["budget_mib"]
            row["requested"] = (dict(workers=1, stored_dataset=True) if job["kind"] == "dataset" else
                                dict(tile_bases=case["tile_bases"], workers=case["workers"]))
            if job["kind"] == "dataset":
                cold_receipt = root / "artifacts" / (job["cache_key"] + "-cold.tsv.manifest.json")
                try:
                    cold = read_json(cold_receipt, 32 << 20)
                    parent = cold["measurements"]["execution.evidence_dataset_manifest"]
                    if cold.get("params", {}).get("run_status") != "completed":
                        raise ValueError("cold receipt is not completed")
                except (OSError, ValueError, KeyError) as error:
                    row.update(status="dependency-failed", error=f"cold dataset publication unavailable: {error}")
                    return
                argv = dataset_command(binary, workload, case, parent, output, receipt)
            else:
                cache = root / "cache" / job["cache_key"] if case.get("cache", "none") != "none" else None
                argv = native_command(binary, workload, case, manifest["oracle"], output, receipt,
                                      cache, job["phase"] == "resumed")
            row["argv"] = argv
            row.update(measure(root, label, argv, root / "raw" / (label + ".stdout.txt")))
            if receipt.is_file():
                claim = read_json(receipt, 32 << 20)
                values = claim.get("measurements", {})
                row["receipt_measurements"] = values
                if "execution.decoder_model" in values:
                    row["decoder_model"] = values["execution.decoder_model"]
                if "execution.decoder_bytes" in values:
                    row["decoder_bytes"] = unsigned(values["execution.decoder_bytes"], "decoder bytes")
                cram = {key.removeprefix("execution.cram."): value for key, value in values.items()
                        if key.startswith("execution.cram.")}
                if cram:
                    row["cram_decoder_measurements"] = cram
                if "validated_records" in cram:
                    row["cram_validation_records"] = unsigned(cram["validated_records"], "CRAM validation records")
                if "validated_bases" in cram:
                    row["cram_validation_bases"] = unsigned(cram["validated_bases"], "CRAM validation bases")
                if "validation_wall_micros" in cram:
                    row["cram_validation_wall_seconds"] = unsigned(cram["validation_wall_micros"], "CRAM validation time") / 1000000
                row["producer"] = {key: value for key, value in claim.get("params", {}).items()
                                   if key.startswith(("producer.", "code_", "deps_", "target_", "rustc_"))}
                if job["kind"] == "dataset":
                    row["effective"] = dict(workers=1, record_visits=unsigned(values["execution.alignment_record_visits"], "record visits"),
                        rows=unsigned(values["execution.emitted_loci"], "rows"), original_sources_rehashed=values.get("original_sources_rehashed"))
                else:
                    row["effective"] = {
                        "workers": unsigned(values.get("execution.worker_count", 1), "workers"),
                        "tile_bases": unsigned(values["execution.microtile_bases"], "tile width"),
                        "microtiles": unsigned(values["execution.microtiles"], "microtiles"),
                        "record_visits": unsigned(values["execution.record_visits"], "record visits"),
                        "rows": unsigned(values["execution.emitted_loci"], "rows"),
                        "computed_partitions": unsigned(values.get("execution.computed_partitions", 0), "computed partitions"),
                        "reused_partitions": unsigned(values.get("execution.reused_partitions", 0), "reused partitions"),
                    }
                row["indexed_alignment_record_visits"] = row["effective"]["record_visits"]
                if row["exit_code"] == 0:
                    validate_resources(row, claim, case)
            elif row["exit_code"] == 0:
                raise ValueError("successful native process did not produce a receipt")
            if row["exit_code"] == 0:
                row["verification"] = measure(root, label + "-verify",
                    [str(binary), "verify", "--manifest", str(receipt), "--json"],
                    root / "raw" / (label + ".verify.json"))
        row["status"] = {0: "completed", 3: "refused", 4: "resource-failed"}.get(row["exit_code"], "failed")
        row["valid"] = row["exit_code"] == 0 and row.get("verification", {"exit_code": 0})["exit_code"] == 0
        if job["kind"] != "oracle":
            row["valid"] &= row.get("resource_validation", {}).get("passed", False)
        if job["kind"] == "dataset":
            row["valid"] &= (row.get("effective", {}).get("record_visits") == 0
                             and row.get("effective", {}).get("original_sources_rehashed") == "false")
    except BaseException as error:
        row.update(valid=False, status="interrupted" if isinstance(error, KeyboardInterrupt) else "metadata-or-harness-failed",
                   error=f"{type(error).__name__}: {error}")
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
    finally:
        row["harness_job_elapsed_seconds"] = time.perf_counter() - began
        try:
            if (root / "raw" / (label + ".time.txt")).is_file():
                row.setdefault("raw_time", str(Path("raw") / (label + ".time.txt")))
            if output.is_file():
                row.update(output_bytes=output.stat().st_size, output_sha256=sha256(output))
            for suffix in (".partial", ".manifest.json.partial"):
                partial = Path(str(output) + suffix)
                if partial.is_file():
                    row.setdefault("partials", []).append(dict(path=str(partial.relative_to(root)),
                        bytes=partial.stat().st_size, sha256=sha256(partial)))
            if row.get("raw_time"):
                row.update(extra_time_metrics(root / row["raw_time"]))
        except (OSError, ValueError) as error:
            row.update(valid=False, error=f"failed to retain artifact metadata: {error}")


def execute(manifest, binary, root, manifest_path, smoke=False, measure=run_measured):
    import pysam
    root.mkdir(parents=True, exist_ok=False)
    (root / "raw").mkdir()
    (root / "artifacts").mkdir()
    groups = job_groups(manifest)
    report = dict(schema=1, label=manifest.get("label", "representative-evidence"), status="running",
                  scientific_profile="shortread-dna-readcount-v1", fields=63, mapq=20, base_quality=20,
                  mode="diagnostic-smoke" if smoke else "representative", manifest=manifest,
                  job_inventory=[job for group in groups for job in group], measurements=[],
                  limitations=[
                      "No OS-enforced memory ceiling is imposed; declared native budgets are cooperative.",
                      "Application cold/resumed cache is distinct from filesystem cache; filesystem caches are not flushed.",
                      "Pysam times cover extraction and TSV encoding; it does not hash sources, seal receipts or verify outputs.",
                      "Native whole-process times include input hashing, extraction, encoding and artifact publication; verification is measured separately.",
                      "Native receipt phase timings combine analysis and encoding; no isolated kernel-only claim is made.",
                      "Zero indexed extraction on resume does not imply zero alignment decoding: CRAM whole-file validation is reported separately when available.",
                      "Setup timings include hashing; they are not disjoint phases. Startup/final harness hashes are outside invocation timing.",
                      "Oracle arrays are bounded by tile width; native decoder/header/aux allocations are observed, with post-decode checks.",
                      "OS filesystem operation counters are retained as operations, not converted to physical byte traffic.",
                      "Measured RSS is process high-water RSS, not a universal hard-RAM guarantee."])
    write_report(root, report)
    try:
        if pysam.__version__ != "0.23.3":
            raise ValueError("oracle requires pysam==0.23.3")
        report["inputs"], report["harness_input_hashing_seconds"] = verify_inputs(manifest)
        harness_paths = [Path(__file__).with_name(name) for name in ("representative.py", "streaming_pysam.py", "run.py")]
        report["execution_inputs"] = [file_identity(path) for path in (binary, manifest_path, *harness_paths)]
        captured = {row["path"]: row for row in report["execution_inputs"]}
        stat = os.statvfs(root)
        report["environment"] = dict(platform=platform.platform(), machine=platform.machine(),
            processor=platform.processor(), logical_cpus=os.cpu_count(), python=sys.version,
            pysam=pysam.__version__, htslib=pysam.__samtools_version__,
            filesystem_block_bytes=stat.f_frsize, filesystem_available_bytes=stat.f_bavail * stat.f_frsize,
            binary_sha256=captured[str(binary)]["sha256"], manifest_sha256=captured[str(manifest_path)]["sha256"],
            binary_version=subprocess.check_output([str(binary), "--version"], text=True).strip(),
            harness_sha256={path.name: captured[str(path)]["sha256"] for path in harness_paths})
        workloads = {workload["id"]: workload for workload in manifest["workloads"]}
        cases = {case["id"]: case for case in manifest["cases"]}
        for group in groups:
            for job in group:
                row = dict(job, status="not-started", valid=False, exit_code=None)
                report["measurements"].append(row)
                print(f"[{len(report['measurements'])}/{len(report['job_inventory'])}] {job['label']}", flush=True)
                try:
                    check_snapshots(report)
                    measure_job(manifest, binary, root, job, workloads[job["workload"]], cases, row, measure)
                    check_snapshots(report)
                except BaseException as error:
                    row.update(valid=False, status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                               error=f"{type(error).__name__}: {error}")
                    raise
                finally:
                    write_report(root, report)
        immutable = final_identity_check(report)
        if not immutable:
            for row in report["measurements"]:
                row["valid"] = False
                row["input_identity_valid"] = False
        comparisons, anchors = [], {}
        for row in report["measurements"]:
            if row["kind"] == "oracle" and row["valid"]:
                anchors.setdefault(row["workload"], row["artifact"])
        for row in report["measurements"]:
            anchor = anchors.get(row["workload"])
            if not row["valid"] or not anchor:
                row["valid"] = False
                continue
            comparison = compare(root / anchor, root / row["artifact"])
            comparison.update(label=row["label"], reference=anchor, artifact=row["artifact"])
            if "expected_rows" in workloads[row["workload"]]:
                comparison["expected_rows_match"] = comparison["lines_including_header"] - 1 == workloads[row["workload"]]["expected_rows"]
            emitted = (row.get("oracle_measurements", {}).get("selected_loci") if row["kind"] == "oracle" else
                       row.get("effective", {}).get("rows"))
            comparison["reported_rows_match"] = emitted == comparison["lines_including_header"] - 1
            row["valid"] &= comparison["equal"] and comparison.get("expected_rows_match", True) and comparison["reported_rows_match"]
            comparisons.append(comparison)
        equivalences, group_anchors = [], {}
        for workload in manifest["workloads"]:
            anchor = anchors.get(workload["id"])
            if not anchor:
                continue
            group = workload.get("equivalence_group", workload["id"])
            previous = group_anchors.setdefault(group, anchor)
            equality = compare(root / previous, root / anchor)
            equivalences.append(dict(group=group, left=previous, right=anchor, **equality))
        report["semantic_comparisons"] = comparisons
        report["format_equivalences"] = equivalences
        report["gates"] = assess_gates(manifest, report["measurements"])
        report["summaries"] = summaries(report["measurements"])
        scientific = all(row["valid"] for row in report["measurements"]) and all(row["equal"] for row in equivalences)
        admitted = all(value["passed"] for value in report["gates"].values())
        report["status"] = ("smoke-passed" if smoke else "passed") if scientific and (admitted or smoke) else "failed-or-refused"
    except BaseException as error:
        report["status"] = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        recorded = {row["label"] for row in report["measurements"]}
        report["measurements"].extend(dict(job, status="not-run", valid=False, exit_code=None,
            error="matrix aborted before invocation") for job in report["job_inventory"] if job["label"] not in recorded)
        if captured_identities(report) and "final_input_checks" not in report:
            final_identity_check(report)
        write_report(root, report)
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
    write_report(root, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plan-only", action="store_true", help="print explicit ordered jobs without hashing or running")
    parser.add_argument("--smoke", action="store_true", help="diagnostic only: permit fewer repeats/incomplete scheduling gates")
    args = parser.parse_args()
    try:
        manifest = load_manifest(args.manifest, args.smoke)
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    if args.plan_only:
        print(json.dumps(dict(schema=1, jobs=[job for group in job_groups(manifest) for job in group]), indent=2))
        return 0
    report = execute(manifest, args.binary.resolve(), args.output.resolve(), args.manifest.resolve(), args.smoke)
    print(f"{report['status']}: {args.output.resolve() / 'report.json'}")
    return 0 if report["status"] in ("passed", "smoke-passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
