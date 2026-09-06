#!/usr/bin/env python3
"""Prepare only a wheel build checkout; tracked release versions stay unchanged."""

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


def prepare(root: Path, requested: str, candidate: str) -> dict:
    actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if not re.fullmatch(r"[0-9a-f]{40}", candidate) or candidate != actual:
        raise ValueError("wheel candidate must be the exact checked-out commit SHA")
    manifest_path, lock_path = root / "Cargo.toml", root / "Cargo.lock"
    manifest, lock = manifest_path.read_text(), lock_path.read_text()
    pattern = r'(\[package\]\s*\nname = "rosalind-bio"\s*\nversion = ")([^"]+)(")'
    found = re.search(pattern, manifest)
    if not found:
        raise ValueError("cannot locate root package version")
    base = found[2]
    version = requested or base
    if version != base and not re.fullmatch(re.escape(base) + r"-rc\.[1-9][0-9]*", version):
        raise ValueError("wheel version must match Cargo or be its numbered release candidate")
    lock_pattern = r'(\[\[package\]\]\s*\nname = "rosalind-bio"\s*\nversion = ")' + re.escape(base) + r'(")'
    if len(re.findall(lock_pattern, lock)) != 1:
        raise ValueError("Cargo.lock root package does not match Cargo.toml")
    if version != base:
        manifest_path.write_text(re.sub(pattern, lambda match: match[1] + version + match[3], manifest, count=1))
        lock_path.write_text(re.sub(lock_pattern, lambda match: match[1] + version + match[2], lock, count=1))
    return {
        "schema": 1,
        "source_commit": actual,
        "source_version": base,
        "rust_version": version,
        "python_version": version.replace("-rc.", "rc"),
        "derived_manifest": version != base,
        "cargo_toml_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "cargo_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="")
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--output", default="wheel-build.json")
    args = parser.parse_args()
    try:
        report = prepare(Path.cwd(), args.version, args.candidate_sha)
    except ValueError as error:
        parser.error(str(error))
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["python_version"])


if __name__ == "__main__":
    main()
