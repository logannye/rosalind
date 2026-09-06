#!/usr/bin/env python3
"""Bind tested wheel bytes to a candidate and stage only a complete release set."""

import argparse
import email.parser
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile

TARGETS = {
    "aarch64-apple-darwin": r"macosx_\d+_\d+_arm64",
    "x86_64-apple-darwin": r"macosx_\d+_\d+_x86_64",
    "x86_64-unknown-linux-gnu": r"manylinux[^.]*_x86_64",
}


def python_version(version: str) -> str:
    if not re.fullmatch(r"\d+\.\d+\.\d+(-rc\.[1-9]\d*)?", version):
        raise ValueError("expected a stable version or numbered release candidate")
    return version.replace("-rc.", "rc")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def wheel_identity(path: Path, target: str, version: str) -> dict:
    expected_version = python_version(version)
    if target not in TARGETS or path.is_symlink() or not path.is_file():
        raise ValueError("expected a regular wheel for a supported target")
    parts = path.name.removesuffix(".whl").split("-")
    if len(parts) != 5 or parts[:4] != ["rosalind_bio", expected_version, "py3", "none"]:
        raise ValueError(f"{path.name}: unexpected package, version, or wheel tags")
    if not all(re.fullmatch(TARGETS[target], tag) for tag in parts[4].split(".")):
        raise ValueError(f"{path.name}: wheel platform does not match {target}")
    with zipfile.ZipFile(path) as archive:
        metadata_paths = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata_paths) != 1:
            raise ValueError(f"{path.name}: expected one distribution metadata record")
        metadata = email.parser.Parser().parsestr(archive.read(metadata_paths[0]).decode("utf-8"))
    if metadata["Name"] != "rosalind-bio" or metadata["Version"] != expected_version:
        raise ValueError(f"{path.name}: distribution metadata does not match requested package/version")
    return {"filename": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def record(directory: Path, build_report: Path, target: str) -> dict:
    report = json.loads(build_report.read_text())
    if not re.fullmatch(r"[0-9a-f]{40}", report.get("source_commit", "")):
        raise ValueError("build report lacks an exact candidate commit")
    version = report["rust_version"]
    if report.get("python_version") != python_version(version):
        raise ValueError("native and Python build versions differ")
    wheels = list(directory.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError("each target build must produce exactly one wheel")
    report["target"] = target
    report["wheel"] = wheel_identity(wheels[0], target, version)
    return report


def stage(artifacts: Path, output: Path, candidate: str, version: str, index: str) -> list[dict]:
    expected_version = python_version(version)
    if not re.fullmatch(r"[0-9a-f]{40}", candidate):
        raise ValueError("candidate must be a full commit SHA")
    if index not in {"pypi", "testpypi"} or ("-rc." in version) != (index == "testpypi"):
        raise ValueError("numbered candidates publish to TestPyPI; stable versions publish to PyPI")
    reports = {path.name: path for path in artifacts.glob("wheel-build-*.json")}
    expected_reports = {f"wheel-build-{target}.json" for target in TARGETS}
    if set(reports) != expected_reports:
        raise ValueError("expected exactly one build report for each supported target")
    verified, filenames = [], set()
    for target in TARGETS:
        report_path = reports[f"wheel-build-{target}.json"]
        if report_path.is_symlink():
            raise ValueError("build reports must be regular files")
        report = json.loads(report_path.read_text())
        if (report.get("source_commit"), report.get("rust_version"), report.get("python_version"), report.get("target")) != (
            candidate, version, expected_version, target
        ):
            raise ValueError(f"{report_path.name}: candidate, version, or target differs")
        identity = report.get("wheel", {})
        filename = identity.get("filename", "")
        if not filename or Path(filename).name != filename or filename in filenames:
            raise ValueError("wheel identities must name distinct files inside the artifact directory")
        actual = wheel_identity(artifacts / filename, target, version)
        if identity != actual:
            raise ValueError(f"{filename}: bytes differ from the tested candidate wheel")
        filenames.add(filename)
        verified.append(actual)
    if {path.name for path in artifacts.glob("*.whl")} != filenames:
        raise ValueError("unrecorded wheel found among release artifacts")
    # Validate the complete set before creating any upload directory.
    output.mkdir(parents=True, exist_ok=False)
    for identity in verified:
        shutil.copyfile(artifacts / identity["filename"], output / identity["filename"])
    return verified


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    capture = commands.add_parser("record")
    capture.add_argument("--directory", type=Path, required=True)
    capture.add_argument("--build-report", type=Path, required=True)
    capture.add_argument("--target", choices=TARGETS, required=True)
    capture.add_argument("--output", type=Path, required=True)
    publish = commands.add_parser("stage")
    publish.add_argument("--artifacts", type=Path, required=True)
    publish.add_argument("--output", type=Path, required=True)
    publish.add_argument("--candidate-sha", required=True)
    publish.add_argument("--version", required=True)
    publish.add_argument("--index", choices=["pypi", "testpypi"], required=True)
    args = parser.parse_args()
    try:
        if args.command == "record":
            report = record(args.directory, args.build_report, args.target)
            with args.output.open("x") as stream:
                stream.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
        else:
            print(json.dumps(stage(args.artifacts, args.output, args.candidate_sha, args.version, args.index), indent=2))
    except (KeyError, ValueError, OSError, zipfile.BadZipFile) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
