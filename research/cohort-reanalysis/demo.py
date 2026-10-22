#!/usr/bin/env python3
"""Retain a known-answer plan, summary, verification and source-free replay.

Requires the already prepared and independently validated three-specimen fixture.
This script invokes an explicit built binary; it never compiles or installs it.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORACLE = ROOT / "examples/cohort-reanalysis/expected.json"

from run_manifest import file_digest, load_json, source_identity, write_json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_summary(path, expected):
    """Compare the retained scientific table to the unchanged authored oracle."""
    aliases = {"callable_total": "total_callable_depth", "alt_total": "total_alt_count",
               "eligible_callable_total": "eligible_callable_depth",
               "eligible_alt_total": "eligible_alt_count"}
    rows = []
    with Path(path).open() as stream:
        for line in stream:
            require(len(rows) < len(expected["partial_summary"]), "unexpected extra summary rows")
            row = json.loads(line)
            require(row["contig"] == expected["contig"] and row["ref"] == expected["reference"],
                    "candidate contig/REF differs from the authored oracle")
            require(row["min_callable_depth"] == expected["thresholds"]["min_callable_depth"],
                    "depth threshold differs from the authored oracle")
            require(row["reference_validation"] == "checked_on_observed_rows",
                    "fixture candidates must have an observed REF check")
            rows.append([row[aliases.get(column, column)] for column in expected["summary_columns"]])
    require(rows == expected["partial_summary"], "summary differs from authored partial-summary oracle")
    return len(rows)


def run_demo(args):
    prepared = args.prepared.resolve(strict=True)
    binary = args.binary.resolve(strict=True)
    output = args.output.absolute()
    require(output.parent.is_dir(), "demo output parent must exist")
    output.mkdir(exist_ok=False)
    logs = output / "commands"
    logs.mkdir()
    report = {"schema": "rosalind-cohort-retained-demo-v1", "status": "running",
              "scope": "synthetic fixture correctness and saved-data replay; no performance claim",
              "prepared": str(prepared), "commands": [], "checks": [],
              "analyzer": source_identity(ROOT)}

    def run(arguments):
        command = list(map(str, arguments))
        index = len(report["commands"])
        stdout = logs / f"{index:02d}.stdout.txt"
        stderr = logs / f"{index:02d}.stderr.txt"
        item = {"argv": command, "stdout": str(stdout.relative_to(output)),
                "stderr": str(stderr.relative_to(output)), "exit_code": None}
        report["commands"].append(item)
        with stdout.open("x") as out, stderr.open("x") as err:
            result = subprocess.run(command, stdout=out, stderr=err, check=False)
        item["exit_code"] = result.returncode
        require(result.returncode == 0, f"demo command {index} failed; see {stderr}")
        return stdout

    raw = prepared / "raw"
    hidden = prepared / "raw.demo-unavailable"
    raw_hidden = False
    try:
        expected = load_json(ORACLE)
        require(len(expected["members"]) == 3 and len(expected["candidates"]) == 4,
                "demo accepts only the authored three-member/four-locus fixture")
        validation_path = prepared / "validation/report.json"
        validation = load_json(validation_path)
        oracle_hash, _ = file_digest(ORACLE)
        binary_hash, binary_size = file_digest(binary)
        require(validation.get("passed") is True, "independent fixture validation must pass first")
        require(validation.get("expectations_sha256") == oracle_hash,
                "independent validation used a different authored oracle")
        require(validation.get("binary_sha256") == binary_hash,
                "demo binary differs from the independently validated artifact")
        require(raw.is_dir() and not raw.is_symlink(), "prepared/raw must be a generated directory")
        require(not hidden.exists() and not hidden.is_symlink(), "temporary raw-hiding path already exists")
        report.update({"oracle": {"path": str(ORACLE), "sha256": oracle_hash},
                       "independent_validation": {"path": str(validation_path),
                                                  "sha256": file_digest(validation_path)[0]},
                       "native_binary": {"path": str(binary), "sha256": binary_hash,
                                         "bytes": binary_size, "version": validation.get("binary_version")}})
        members = []
        for specimen in sorted(expected["members"]):
            leaf = prepared / "validation/portable" / specimen / "evidence-dataset.manifest.json"
            require(leaf.is_file(), f"explicit portable fixture leaf missing: {leaf}")
            members.append({"specimen_id": specimen,
                            "dataset": os.path.relpath(leaf, output)})
        members_path = output / "members.json"
        write_json(members_path, {"version": 1, "members": members})
        report["members_manifest_sha256"] = file_digest(members_path)[0]
        report["leaf_manifest_sha256"] = {
            member["specimen_id"]: file_digest(output / member["dataset"])[0] for member in members}
        sites = ROOT / "examples/cohort-reanalysis/candidates.vcf"
        report["candidate_vcf_sha256"] = file_digest(sites)[0]
        base = [sys.executable, HERE / "analyze.py", "--binary", binary,
                "--members", members_path, "--sites", sites, "--missingness", "partial",
                "--window-size", "2", "--min-callable-depth", "10",
                "--native-memory-budget-mb", "256", "--parent-memory-budget-mb", "128"]
        plan_output, run_output, replay_output = output / "plan", output / "run", output / "replay"
        run([*base, "--plan", "--output", plan_output])
        require(not (plan_output / "run-manifest.json").exists(), "metadata plan must not claim completed aggregation")
        plan = load_json(plan_output / "plan.json")
        require(plan["verification_level"] == "metadata_only", "plan verification scope changed")
        report["checks"].append("retained metadata-only plan")
        run([*base, "--output", run_output])
        report["summary_rows"] = check_summary(run_output / "summary.jsonl", expected)
        report["checks"].append("partial summary equals unchanged authored oracle")
        run([sys.executable, HERE / "analyze.py", "verify", run_output])
        report["checks"].append("research output manifest and retained bytes verified")

        # Only the fixture generator's raw directory is renamed. Every portable
        # leaf stays in place; try/finally restores raw even on ordinary failures.
        raw.rename(hidden)
        raw_hidden = True
        try:
            require(not raw.exists(), "original generated raw paths must be unavailable")
            for member in members:
                leaf = (output / member["dataset"]).resolve(strict=True)
                verified_stdout = run([binary, "dataset", "verify", "--dataset", leaf,
                                       "--memory-budget-mb", "256"])
                verified = load_json(verified_stdout)
                require(verified["verified_loci"] == len(expected["members"][member["specimen_id"]]["stored_positions"]),
                        "portable leaf verified denominator differs from fixture")
                require(verified.get("original_sources_rehashed") is False,
                        "native saved-data verification scope changed")
            report["checks"].append("all portable leaf payloads verified with original raw paths unavailable")
            run([sys.executable, HERE / "analyze.py", "replay", run_output,
                 "--members", members_path, "--binary", binary, "--window-size", "1",
                 "--output", replay_output])
            check_summary(replay_output / "summary.jsonl", expected)
            require(file_digest(run_output / "summary.jsonl")[0] == file_digest(replay_output / "summary.jsonl")[0],
                    "window-one offline replay differs from window-two summary")
            run([sys.executable, HERE / "analyze.py", "verify", replay_output])
            report["checks"].append("source-free window-one replay reproduces canonical science")
        finally:
            if raw_hidden:
                hidden.rename(raw)
                raw_hidden = False
        report["raw_restored"] = raw.is_dir()
        report["outputs"] = {name: {"path": name, "sha256": file_digest(output / name)[0]}
                             for name in ("plan/plan.json", "run/summary.jsonl", "run/scientific.json",
                                          "run/run-manifest.json", "replay/summary.jsonl", "replay.replay.json")}
        report["verification_scope"] = {
            "native": "dataset verify checks every portable leaf row/partition; original raw sources not rehashed",
            "outer": "research manifest and retained-output byte checks; native receipts do not certify outer reducer",
            "replay": "canonical summary/scientific metadata agreement across windows 2 and 1"}
        report["status"] = "passed"
    except BaseException as error:
        report["status"] = "failed"
        report["error"] = {"type": type(error).__name__, "message": str(error)}
        report["raw_restored"] = raw.is_dir() and not raw_hidden
        raise
    finally:
        if raw_hidden:
            hidden.rename(raw)
            report["raw_restored"] = True
        write_json(output / "report.json", report)
    print(output / "report.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    def cancelled(signum, _frame):
        raise KeyboardInterrupt(f"demo cancelled by signal {signum}")

    signal.signal(signal.SIGTERM, cancelled)
    run_demo(args)


if __name__ == "__main__":
    main()
