#!/usr/bin/env python3
"""Bind an explicitly selected executable to this case study's exact source.

This checks consistency with a retained, trusted build record. It is not a
signature, software supply-chain attestation, or protection against replacing
both the executable and all its purported provenance.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "research/reproducibility-case-study/protocol-v1.json"
REQUIREMENTS = ROOT / "research/cohort-reanalysis/workloads/requirements.lock"
SCHEMA = "rosalind-case-study-build-v1"
NATIVE_KEYS = ("code_git_sha", "code_dirty", "rustc_version", "target_triple", "deps_lock_blake3")
SOURCE_PATHS = (
    "Cargo.toml", "Cargo.lock", "build.rs", "src", "crates", "tests", "xtask", "python/rosalind",
    "benchmarks/evidence", "benchmarks/giab/resources.tsv", "scripts/reproduce_demo.py",
    "scripts/reproduce_demo.sh", "examples/research-filter", "examples/cohort-reanalysis",
    "research/cohort-reanalysis/benchmark.py",
    "research/cohort-reanalysis/workloads/requirements.lock", "research/reproducibility-case-study",
    ".github/workflows/reproducibility-case-study.yml",
)


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _command(argv: list[str], *, cwd: Path | None = None) -> str:
    return subprocess.run(argv, cwd=cwd, check=True, capture_output=True,
                          text=True, timeout=120).stdout.strip()


def _source_sha(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{40}", value):
        raise ValueError("expected source must be a complete lowercase 40-character commit SHA")
    return value


def source_identity(expected_source: str) -> dict:
    expected_source = _source_sha(expected_source)
    revision = _command(["git", "rev-parse", "HEAD"], cwd=ROOT)
    if revision != expected_source:
        raise ValueError(f"checkout/source mismatch: {revision} != {expected_source}")
    status = _command(["git", "status", "--porcelain", "--untracked-files=all"], cwd=ROOT)
    if status:
        raise ValueError(f"source checkout must be clean: {status}")
    tracked = _command(["git", "ls-files", "--", *SOURCE_PATHS], cwd=ROOT).splitlines()
    hashes = {name: sha256(ROOT / name) for name in tracked}
    for required in ("Cargo.lock", "research/reproducibility-case-study/protocol-v1.json",
                     ".github/workflows/reproducibility-case-study.yml"):
        if required not in hashes:
            raise ValueError(f"required case-study source is not tracked: {required}")
    return {"source_revision": revision,
            "source_tree": _command(["git", "rev-parse", "HEAD^{tree}"], cwd=ROOT),
            "source_status": status, "source_sha256": hashes,
            "cargo_lock_sha256": sha256(ROOT / "Cargo.lock"),
            "protocol_sha256": sha256(PROTOCOL), "requirements_sha256": sha256(REQUIREMENTS)}


def dependencies() -> dict[str, str]:
    expected = dict(line.split("==") for line in REQUIREMENTS.read_text().splitlines()
                    if line and not line.startswith("#"))
    observed = {name: importlib.metadata.version(name) for name in expected}
    if observed != expected:
        raise ValueError(f"pinned wheel versions differ: {observed!r} != {expected!r}")
    if sys.version_info[:2] != (3, 11):
        raise ValueError("case-study environment requires Python 3.11")
    return observed


def environment_identity() -> dict:
    event = {}
    event_path = os.environ.get("GITHUB_EVENT_PATH")
    if event_path:
        payload = json.loads(Path(event_path).read_text())
        pr = payload.get("pull_request", {})
        event = {"pull_request_head_sha": pr.get("head", {}).get("sha"),
                 "pull_request_merge_sha": pr.get("merge_commit_sha")}
    cpu = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next((line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                    if line.startswith("model name")), cpu)
    compiler = shutil.which("rustc")
    return {"recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version, "dependencies": dependencies(), "platform": platform.platform(),
            "machine": platform.machine(), "cpu": cpu, "logical_cpu_count": os.cpu_count(),
            "available_rustc_verbose": _command([compiler, "--version", "--verbose"]) if compiler else None,
            "github": {key.lower(): os.environ.get(key) for key in (
                "GITHUB_SHA", "GITHUB_WORKFLOW_SHA", "GITHUB_WORKFLOW_REF", "GITHUB_RUN_ID",
                "GITHUB_RUN_ATTEMPT", "GITHUB_EVENT_NAME", "GITHUB_REPOSITORY", "RUNNER_OS",
                "RUNNER_ARCH", "ImageOS", "ImageVersion")}, **event}


def native_metadata(binary: Path | str) -> dict:
    """Observe compiled metadata from a fresh tiny receipt, outside measurements.

    There is no native build-info CLI. Reference-pack construction stamps the
    same five compiled fields used by scientific receipts and needs no raw data.
    """
    binary = Path(binary).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="rosalind-identity-") as directory:
        root = Path(directory)
        reference = root / "reference.fa"
        reference.write_text(">identity\nACGT\n")
        output = root / "reference.rref"
        _command([str(binary), "reference", "build", "--fasta", str(reference), "--output", str(output)])
        receipt = Path(str(output) + ".manifest.json")
        params = json.loads(receipt.read_text())["params"]
        observed = {key: params.get(key) for key in NATIVE_KEYS}
        if any(value in (None, "", "unknown") for value in observed.values()):
            raise ValueError(f"missing compiled native identity: {observed}")
        _command([str(binary), "verify", "--manifest", str(receipt),
                  "--expect-code", observed["code_git_sha"]])
        return observed


def _validate_native(native: dict, expected_source: str, rustc: str) -> None:
    import blake3

    expected = {"code_git_sha": expected_source, "code_dirty": "false",
                "rustc_version": rustc.splitlines()[0],
                "deps_lock_blake3": blake3.blake3((ROOT / "Cargo.lock").read_bytes()).hexdigest()}
    for name, value in expected.items():
        if native.get(name) != value:
            raise ValueError(f"compiled native {name} mismatch: {native.get(name)!r} != {value!r}")
    target = re.search(r"^host: (.+)$", rustc, re.MULTILINE)
    if not target or native.get("target_triple") != target.group(1):
        raise ValueError("compiled target differs from recorded native compiler host")


def create_record(binary: Path | str, expected_source: str) -> dict:
    binary = Path(binary).resolve(strict=True)
    source = source_identity(expected_source)
    digest = sha256(binary)
    rustc = _command(["rustc", "--version", "--verbose"])
    native = native_metadata(binary)
    _validate_native(native, expected_source, rustc)
    version = _command([str(binary), "--version"])
    if sha256(binary) != digest:
        raise ValueError("executable changed during build identity inspection")
    return {"schema_version": SCHEMA, **source, "binary_sha256": digest,
            "binary_version_assertion": version, "binary_size_bytes": binary.stat().st_size,
            "native_build_metadata": native, "rustc_verbose": rustc,
            "environment": environment_identity(),
            "trust_boundary": "consistency with a retained trusted build record; not a signature"}


def authenticate(binary: Path | str, identity_path: Path | str, expected_source: str) -> dict:
    """Authenticate explicit binary/source binding; no PATH fallback or build."""
    binary = Path(binary).resolve(strict=True)
    identity = json.loads(Path(identity_path).read_text())
    if identity.get("schema_version") != SCHEMA:
        raise ValueError("unsupported build identity schema")
    current_source = source_identity(expected_source)
    for key, value in current_source.items():
        if identity.get(key) != value:
            raise ValueError(f"build identity {key} differs from the exact checked-out source")
    if sha256(binary) != identity.get("binary_sha256"):
        raise ValueError("selected executable digest differs from the trusted build record")
    if binary.stat().st_size != identity.get("binary_size_bytes"):
        raise ValueError("selected executable size differs from the trusted build record")
    if dependencies() != identity["environment"]["dependencies"]:
        raise ValueError("replay dependencies differ from the producer")
    if _command([str(binary), "--version"]) != identity.get("binary_version_assertion"):
        raise ValueError("selected executable version assertion differs")
    native = native_metadata(binary)
    if native != identity.get("native_build_metadata"):
        raise ValueError("selected executable compiled metadata differs from the build record")
    _validate_native(native, expected_source, identity["rustc_verbose"])
    if sha256(binary) != identity["binary_sha256"]:
        raise ValueError("executable changed during authentication")
    return identity


def resource_preflight(scratch: Path) -> dict:
    import psutil

    scratch.mkdir(parents=True, exist_ok=True)
    report = {"scratch_free_bytes": shutil.disk_usage(scratch).free,
              "ram_available_bytes": psutil.virtual_memory().available,
              "ram_total_bytes": psutil.virtual_memory().total,
              "required_scratch_free_bytes": 20 * 1024**3,
              "required_ram_available_bytes": 8 * 1024**3,
              "expected_public_download_bytes": 1980479678,
              "planned_scratch_bytes": {"downloads": 2 * 1024**3,
                                         "uncompressed_reference_and_preparation": 7 * 1024**3,
                                         "native_target_reserve": 8 * 1024**3,
                                         "reports_and_headroom": 3 * 1024**3},
              "note": "planning reserves, not measured build/workload requirements"}
    report["passed"] = (report["scratch_free_bytes"] >= report["required_scratch_free_bytes"]
                        and report["ram_available_bytes"] >= report["required_ram_available_bytes"])
    return report


def run_test_suite(start: str, pattern: str, output: Path) -> bool:
    """A green process with skipped native tests is not a passing hosted gate."""
    suite = unittest.defaultTestLoader.discover(start, pattern=pattern)
    started = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = {"start": start, "pattern": pattern, "tests_run": result.testsRun,
              "elapsed_seconds": time.perf_counter() - started,
              "skips": [{"test": str(test), "reason": reason} for test, reason in result.skipped],
              "failures": [{"test": str(test), "detail": detail} for test, detail in result.failures],
              "errors": [{"test": str(test), "detail": detail} for test, detail in result.errors],
              "expected_failures": [str(test) for test, _ in result.expectedFailures],
              "unexpected_successes": [str(test) for test in result.unexpectedSuccesses]}
    report["passed"] = bool(result.testsRun and result.wasSuccessful()
                            and not result.skipped and not result.expectedFailures)
    _write_new(output, report)
    return report["passed"]


def _write_new(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    record = commands.add_parser("record")
    record.add_argument("--binary", required=True, type=Path)
    record.add_argument("--expected-source", required=True)
    record.add_argument("--output", required=True, type=Path)
    check = commands.add_parser("authenticate")
    check.add_argument("--binary", required=True, type=Path)
    check.add_argument("--build-identity", required=True, type=Path)
    check.add_argument("--expected-source", required=True)
    check.add_argument("--output", required=True, type=Path)
    preflight = commands.add_parser("preflight")
    preflight.add_argument("--scratch", required=True, type=Path)
    preflight.add_argument("--expected-source", required=True)
    preflight.add_argument("--output", required=True, type=Path)
    tests = commands.add_parser("test")
    tests.add_argument("--start", required=True)
    tests.add_argument("--pattern", default="test_*.py")
    tests.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "test":
        return 0 if run_test_suite(args.start, args.pattern, args.output) else 1
    if args.command == "record":
        value = create_record(args.binary, args.expected_source)
    elif args.command == "authenticate":
        value = {"build_identity": authenticate(args.binary, args.build_identity, args.expected_source),
                 "verifier_environment": environment_identity(), "passed": True}
    else:
        value = {**resource_preflight(args.scratch), **source_identity(args.expected_source),
                 "environment": environment_identity()}
    _write_new(args.output, value)
    return 0 if value.get("passed", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
