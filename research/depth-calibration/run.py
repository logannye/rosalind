"""Hosted, identity-bound reanalysis of saved HG002 features; no native rebuild."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "evidence-sufficiency"))


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: dict) -> None:
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def verify_file(path: Path, expected: dict) -> None:
    if path.stat().st_size != expected["bytes"] or digest(path) != expected["sha256"]:
        raise ValueError(f"Saved evidence identity mismatch: {path.name}")


def verify_provider(run: dict, artifact: dict, source: dict) -> None:
    """Bind bytes to the authorized repository/run, not just a self-signed manifest."""
    if (run.get("id") != source["run_id"] or run.get("head_sha") != source["commit"]
            or run.get("conclusion") != "success"
            or run.get("repository", {}).get("full_name") != source["repository"]):
        raise ValueError("Unexpected source workflow identity or unsuccessful run")
    arun = artifact.get("workflow_run", {})
    if (artifact.get("id") != source["artifact_id"] or artifact.get("expired")
            or artifact.get("digest") != "sha256:" + source["artifact_sha256"]
            or arun.get("id") != source["run_id"] or arun.get("head_sha") != source["commit"]):
        raise ValueError("Unexpected source artifact provenance")


def extract_verified(archive: Path, destination: Path, source: dict) -> dict:
    if digest(archive) != source["artifact_sha256"]:
        raise ValueError("Source ZIP digest mismatch")
    destination.mkdir(parents=True, exist_ok=False)
    # Never extract arbitrary ZIP paths or deserialize an unverified pickle.
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        for relative, expected in source["files"].items():
            wanted = "hg002-p0/" + relative
            if names.count(wanted) != 1:
                raise ValueError(f"Missing or duplicate saved artifact member: {wanted}")
            info = bundle.getinfo(wanted)
            if info.file_size != expected["bytes"]:
                raise ValueError(f"Saved member size mismatch: {wanted}")
            path = destination / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            with bundle.open(info) as src, path.open("xb") as dst:
                for block in iter(lambda: src.read(1 << 20), b""):
                    dst.write(block)
            verify_file(path, expected)
    manifest = json.loads((destination / "run-manifest.json").read_text())
    for relative, expected in source["files"].items():
        # P0 writes final telemetry after its inventory and explicitly excludes
        # that file (run_pilot.py). Its ZIP and member hashes are still pinned
        # and verified above; no model/data file receives this exception.
        if relative == "run-manifest.json" or (relative == "telemetry.json" and relative not in manifest["files"]):
            continue
        if manifest["files"].get(relative) != expected:
            raise ValueError(f"Pinned identity differs from original inventory: {relative}")
    identity = json.loads((destination / "identity.json").read_text())
    if (identity["git_revision"] != source["commit"]
            or int(identity["github_run_id"]) != source["run_id"]):
        raise ValueError("Saved execution identity mismatch")
    original = json.loads((destination / "report.json").read_text())
    if (original["status"] != "completed" or original["gate"]["status"] != "stop_no_signal"
            or original["p1_authorized_by_gate"]):
        raise ValueError("Original stop decision differs from the frozen premise")
    return identity


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--run-metadata", type=Path, required=True)
    parser.add_argument("--artifact-metadata", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if os.environ.get("GITHUB_ACTIONS") != "true" or platform.system() != "Linux":
        raise RuntimeError("Scientific execution requires the authorized hosted Linux workflow")
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    try:
        protocol = json.loads(args.protocol.read_text())
        run = json.loads(args.run_metadata.read_text())
        artifact = json.loads(args.artifact_metadata.read_text())
        verify_provider(run, artifact, protocol["source"])
        saved = output.parent / "verified-source"
        source_identity = extract_verified(args.archive, saved, protocol["source"])
        # Reuse exactly the installed binary wheels used to serialize the trusted models.
        dependencies = {name: importlib.metadata.version(name)
                        for name in source_identity["dependencies"]}
        if dependencies != source_identity["dependencies"]:
            raise ValueError("Dependency versions differ from saved model environment")
        import psutil
        import study
        from calibration import run_experiment

        if digest(HERE.parent / "evidence-sufficiency" / "study.py") != source_identity["source_sha256"]["research/evidence-sufficiency/study.py"]:
            raise ValueError("Original model helper source changed")
        stop = threading.Event()
        sampled = {"peak_rss_bytes": 0}

        def monitor() -> None:
            process = psutil.Process()
            while not stop.wait(.1):
                sampled["peak_rss_bytes"] = max(sampled["peak_rss_bytes"], process.memory_info().rss)
                if (sampled["peak_rss_bytes"] > protocol["limits"]["rss_bytes"]
                        or time.perf_counter() - started > protocol["limits"]["analysis_seconds"]):
                    write_json(output / "failure.json", {"reason": "resource_limit", **sampled})
                    os._exit(2)

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()
        write_json(output / "protocol.json", protocol)
        source_paths = [HERE / "run.py", HERE / "calibration.py", HERE / "protocol.json",
                        HERE / "PROTOCOL.md", HERE.parent / "evidence-sufficiency" / "study.py",
                        ROOT / ".github/workflows/depth-calibration.yml"]
        write_json(output / "identity.json", {
            "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
            "github_run_id": os.environ["GITHUB_RUN_ID"], "github_run_attempt": os.environ["GITHUB_RUN_ATTEMPT"],
            "platform": platform.platform(), "python": sys.version, "dependencies": dependencies,
            "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in source_paths},
            "input_source": protocol["source"], "original_p0_status": "stop_no_signal",
        })
        with (saved / "records.jsonl").open() as handle:
            records = [json.loads(line) for line in handle if line.strip()]
        bundle = study.load_bundle(saved / "study")
        write_json(output / "source-truth-recovery.json", {
            "candidate_recovery": json.loads((saved / "candidate-recall.json").read_text()),
            "original_p0_truth_recovery": json.loads((saved / "truth-recovery.json").read_text()),
            "scope": "Historical P0 recovery only. New calibration metrics remain conditional on the same original candidate list; no de novo reduced-depth recovery is measured.",
        })
        report = run_experiment(records, bundle, output, protocol)
        for relative, expected in protocol["source"]["files"].items():
            verify_file(saved / relative, expected)
        write_json(output / "input-recheck.json", {"passed": True, "files_verified": len(protocol["source"]["files"])})
        stop.set()
        thread.join(timeout=1)
        retained = sum(p.stat().st_size for p in output.rglob("*") if p.is_file())
        if retained > protocol["limits"]["output_bytes"]:
            raise ValueError("Output storage ceiling exceeded")
        write_json(output / "telemetry.json", {
            "elapsed_seconds": time.perf_counter() - started, **sampled,
            "retained_bytes_before_inventory": retained,
            "scope": "Saved artifact verification, feature loading and model reanalysis; excludes source ZIP download, dependency installation, native compilation and raw read acquisition",
            "input_archive_bytes": args.archive.stat().st_size,
        })
        write_json(output / "run-manifest.json", {"files": {
            str(p.relative_to(output)): {"bytes": p.stat().st_size, "sha256": digest(p)}
            for p in sorted(output.rglob("*")) if p.is_file()}})
        write_json(output / "completion.json", {
            "status": "completed", "input_identity_rechecked": True,
            "resource_checks_passed": True, "report_sha256": digest(output / "report.json"),
            "run_manifest_sha256": digest(output / "run-manifest.json"),
            "p1_authorized": False, "selective_escalation_authorized": False,
        })
        print(json.dumps({"gate": report["gate"], "selected_comparators": report["selected_comparators"]}, indent=2))
    except BaseException as exc:
        if not (output / "failure.json").exists():
            write_json(output / "failure.json", {"type": type(exc).__name__, "reason": str(exc)})
        raise


if __name__ == "__main__":
    main()
