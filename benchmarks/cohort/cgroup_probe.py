#!/usr/bin/env python3
"""Compare saved-only cohort reports inside a real 512 MiB Linux cgroup.

Consumes a passed synthetic cohort demonstration; no original alignments or
references are mounted. Reports record observations, not a general memory claim.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

EVIDENCE_HELPER = Path(__file__).resolve().parents[1] / "evidence" / "cgroup_probe.py"
_spec = importlib.util.spec_from_file_location("evidence_cgroup_helpers", EVIDENCE_HELPER)
evidence = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(evidence)

# Same immutable Linux/amd64 runtime used by the retained single-sample probes.
DEFAULT_IMAGE = "rust:1.95-bookworm@sha256:6258907abe69656e41cd992e0b705cdcfabcbbe3db374f92ed2d47121282d4a1"
HARD_LIMIT_BYTES = 512 << 20


def require(condition, message):
    if not condition:
        raise ValueError(message)


def command(binary, operation, cohort, snapshot, sites, output, *, budget=512, os_limit=False, pairs=None):
    require((operation == "compare-pairs") == (pairs is not None), "compare-pairs requires an explicit pair table; other operations do not")
    argv = [str(binary), "cohort", operation, "--cohort", str(cohort), "--snapshot", snapshot,
            "--sites", str(sites), "--missing", "partial", "--fields", "depths,alleles",
            "--format", "tsv", "--output", str(output), "--memory-budget-mb", str(budget), "--enforce"]
    if os_limit:
        argv.append("--require-os-limit")
    if pairs is not None:
        argv += ["--pairs", str(pairs)]
    return argv


def container_arguments(name, image, binary, cohort, sites, destination, argv, *, pairs=None):
    arguments = ["create", "--name", name, "--platform", "linux/amd64", "--network", "none",
            "--memory", "512m", "--memory-swap", "512m", "--cpus", "2", "--pids-limit", "64",
            "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--read-only",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--mount", f"type=bind,src={binary},dst=/inputs/rosalind,readonly",
            "--mount", f"type=bind,src={cohort},dst=/cohort,readonly",
            "--mount", f"type=bind,src={sites},dst=/inputs/candidates.vcf,readonly",
            "--mount", f"type=bind,src={destination},dst=/output"]
    if pairs is not None:
        arguments += ["--mount", f"type=bind,src={pairs},dst=/inputs/pairs.tsv,readonly"]
    return arguments + [image, "/bin/sh", "-c", evidence.CONTROLLER, "cohort-cgroup-controller", *argv]


def validate_case(item, expected_sha256):
    errors = list(item["errors"])
    if item["memory_peak_bytes"] is None:
        errors.append("cgroup memory.peak was not observed")
    if any(not isinstance(item[key], dict) or "oom_kill" not in item[key]
           for key in ("events_before", "events_after")):
        errors.append("before/after memory.events oom_kill counters were not observed")
    if item["analysis_exit_code"] != item["container_exit_code"]:
        errors.append("native and container exit codes disagree")
    if item["case"]["expected"] == "completed":
        if item.get("measurements", {}).get("execution.original_alignment_records_decoded") != "0":
            errors.append("saved-only zero original-alignment decoding measurement is absent")
        observed = item["files"].get("result.tsv", {}).get("sha256")
        if observed != expected_sha256:
            errors.append("cohort TSV differs from the outside-probe-cgroup baseline")
        if item.get("measurements", {}).get("resource.os_limit_bytes") != str(HARD_LIMIT_BYTES):
            errors.append("native receipt did not record the observed 512 MiB OS limit")
    return errors


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--demo-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pairs", type=Path, help="Also exercise explicit pairs, including tiny-budget refusal")
    parser.add_argument("--docker-context", default="default")
    parser.add_argument("--image", default=DEFAULT_IMAGE)
    parser.add_argument("--timeout", default=180, type=int)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    (output / "harness.py").write_bytes(Path(__file__).read_bytes())
    (output / "evidence_cgroup_helpers.py").write_bytes(EVIDENCE_HELPER.read_bytes())
    raw = output / "raw"
    raw.mkdir()
    report = {"schema": 1, "kind": "rosalind-cohort-cgroup-probe", "status": "running",
              "cases": [], "baselines": {}, "docker_context": args.docker_context,
              "runtime_image": args.image, "harness_sha256": evidence.sha256(Path(__file__)),
              "helper_sha256": evidence.sha256(EVIDENCE_HELPER),
              "limitations": ["Three authored synthetic members; not a scalability or performance claim.",
                              "Cgroup memory includes page cache and the controller; receipt peak is process RSS.",
                              "Only saved cohort objects, candidates and the binary are mounted; no original sources.",
                              "The baseline runs on the CI host outside these probe containers, not on an uncapped OS."]}
    report_path = output / "report.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    save()
    docker = evidence.Docker(args.docker_context, raw)
    created = []
    identity = uuid.uuid4().hex[:12]
    try:
        require(args.timeout > 0, "timeout must be positive")
        require(re.fullmatch(r"[^\s@]+@sha256:[0-9a-f]{64}", args.image), "runtime image must use an immutable digest")
        binary = args.binary.resolve(strict=True)
        with binary.open("rb") as stream:
            header = stream.read(20)
        require(header[:4] == b"\x7fELF" and header[18:20] == b">\x00", "binary must be Linux/amd64 ELF")
        demonstration = evidence.read_json(args.demo_report)
        require(demonstration.get("status") == "passed", "cohort demonstration must pass first")
        cohort = Path(demonstration["cohort"]).resolve(strict=True)
        snapshot = demonstration["parent_snapshot"]
        require(re.fullmatch(r"[0-9a-f]{64}", snapshot), "invalid snapshot identity")
        sites = (args.demo_report.resolve().parent / "first.vcf").resolve(strict=True)
        pairs = args.pairs.resolve(strict=True) if args.pairs is not None else None
        require(binary.is_file() and cohort.is_dir() and sites.is_file(), "required inputs are absent")
        for path in (binary, cohort, sites, output, *((pairs,) if pairs is not None else ())):
            require("," not in str(path), "Docker bind paths containing commas are unsupported")
        report["inputs"] = {"binary_sha256": evidence.sha256(binary), "snapshot_id": snapshot,
                            "candidate_sha256": evidence.sha256(sites),
                            "demo_report_sha256": evidence.sha256(args.demo_report)}
        if pairs is not None:
            require(pairs.is_file(), "pair table must be a file")
            report["inputs"]["pairs_sha256"] = evidence.sha256(pairs)
        report["binary_version"] = subprocess.check_output([str(binary), "--version"], text=True, timeout=args.timeout).strip()
        operations = ["extract", "summarize"] + (["compare-pairs"] if pairs is not None else [])
        for operation in operations:
            baseline = output / f"baseline-{operation}"
            baseline.mkdir()
            result_path = baseline / "result.tsv"
            invocation = command(binary, operation, cohort, snapshot, sites, result_path,
                                 pairs=pairs if operation == "compare-pairs" else None)
            (baseline / "argv.json").write_text(json.dumps(invocation, indent=2) + "\n")
            started = time.monotonic()
            with (baseline / "stdout").open("wb") as stdout, (baseline / "stderr").open("wb") as stderr:
                proc = subprocess.run(invocation, stdout=stdout, stderr=stderr, timeout=args.timeout, check=False)
            report["baselines"][operation] = {"exit_code": proc.returncode, "elapsed_seconds": time.monotonic() - started}
            save()
            require(proc.returncode == 0 and result_path.is_file(), f"{operation} baseline failed")
            report["baselines"][operation]["tsv_sha256"] = evidence.sha256(result_path)
        docker.call(["pull", "--platform", "linux/amd64", args.image], timeout=600)
        _, image_json = docker.call(["image", "inspect", args.image])
        image = json.loads(image_json)[0]
        require(image.get("Architecture") == "amd64" and image.get("Os") == "linux", "runtime image is not Linux/amd64")
        report["image"] = {key: image.get(key) for key in ("Id", "RepoDigests", "Architecture", "Os")}
        _, daemon = docker.call(["info", "--format",
            '{"kernel":{{json .KernelVersion}},"architecture":{{json .Architecture}},"os":{{json .OperatingSystem}},"cgroup_version":{{json .CgroupVersion}}}'])
        (raw / "daemon.json").write_text(daemon)
        daemon_info = json.loads(daemon)
        report["docker_daemon"] = daemon_info
        require(daemon_info.get("cgroup_version") == "2", "Linux cgroup v2 is required")
        cases = [
            {"name": "extract-os-limit", "operation": "extract", "expected": "completed", "budget_mib": 512, "require_os_limit": True},
            {"name": "summarize-os-limit", "operation": "summarize", "expected": "completed", "budget_mib": 512, "require_os_limit": True},
            {"name": "tiny-budget-refusal", "operation": "extract", "expected": "refused", "budget_mib": 1},
        ]
        if pairs is not None:
            cases += [
                {"name": "pairs-os-limit", "operation": "compare-pairs", "expected": "completed", "budget_mib": 512, "require_os_limit": True},
                {"name": "pairs-tiny-budget-refusal", "operation": "compare-pairs", "expected": "refused", "budget_mib": 1},
            ]
        for case in cases:
            case["hard_limit_bytes"] = HARD_LIMIT_BYTES
            destination = output / case["name"]
            destination.mkdir()
            name = f"rosalind-cohort-{identity}-{case['name']}"
            invocation = command("/inputs/rosalind", case["operation"], "/cohort", snapshot,
                                 "/inputs/candidates.vcf", "/output/result.tsv", budget=case["budget_mib"],
                                 os_limit=case.get("require_os_limit", False),
                                 pairs="/inputs/pairs.tsv" if case["operation"] == "compare-pairs" else None)
            case["argv"] = invocation
            docker.call(container_arguments(name, args.image, binary, cohort, sites, destination, invocation,
                                            pairs=pairs if case["operation"] == "compare-pairs" else None))
            created.append(name)
            timed_out = False
            started = time.monotonic()
            try:
                docker.call(["start", "--attach", name], check=False, timeout=args.timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                docker.call(["kill", name], check=False)
            _, inspection = docker.call(["inspect", name])
            (destination / "docker-inspect.json").write_text(inspection)
            item = evidence.collect(destination, case, json.loads(inspection)[0], timed_out)
            item["wall_seconds"] = time.monotonic() - started
            item["errors"] = validate_case(item, report["baselines"][case["operation"]]["tsv_sha256"])
            item["passed"] = not item["errors"]
            report["cases"].append(item)
            save()
            docker.call(["rm", name])
            created.remove(name)
        report["status"] = "passed" if all(case["passed"] for case in report["cases"]) else "failed"
    except Exception as error:
        report["status"] = "failed"
        report["failure"] = f"{type(error).__name__}: {error}"
    finally:
        for name in reversed(created):
            try:
                code, _ = docker.call(["rm", "--force", name], check=False)
                if code:
                    raise RuntimeError(f"cleanup exit {code}")
            except Exception as error:
                report.setdefault("cleanup_errors", []).append(f"{name}: {error}")
                report["status"] = "failed"
        save()
    print(json.dumps({"status": report["status"], "report": str(report_path)}))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
