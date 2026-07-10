#!/usr/bin/env python3
"""Compare attested GIAB evidence without ever mutating the accepted baseline."""

import argparse
import json
from pathlib import Path
import shutil

EXACT_FIELDS = (
    "calls_filter_all",
    "calls_filter_pass",
    "external_happy_vcfeval",
    "memory",
    "receipt_claim",
    "producer_identity",
    "command_argv",
    "data_manifest",
)


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".new")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def compare(baseline_path: Path, latest_path: Path, results: Path) -> str:
    baseline = json.loads(baseline_path.read_text())
    latest = json.loads(latest_path.read_text())
    missing = sorted(set(EXACT_FIELDS) - latest.keys())
    if missing:
        raise ValueError(f"latest report lacks exact-comparison fields: {missing}")
    results.mkdir(parents=True, exist_ok=True)
    candidate = results / "baseline-candidate.json"
    if candidate.exists():
        candidate.unlink()

    if baseline.get("status") == "not-yet-established":
        status, matching = "candidate-pending-review", None
        shutil.copyfile(latest_path, candidate)
    else:
        baseline_missing = sorted(set(EXACT_FIELDS) - baseline.keys())
        if baseline_missing:
            raise ValueError(
                f"accepted baseline lacks exact-comparison fields: {baseline_missing}"
            )
        matching = all(baseline[field] == latest[field] for field in EXACT_FIELDS)
        status = "reproduced" if matching else "diverged"
        if not matching:
            shutil.copyfile(latest_path, candidate)

    write_json(
        results / "comparison.json",
        {"schema": 1, "status": status, "matching": matching},
    )
    external = latest["external_happy_vcfeval"]
    memory = latest["memory"]
    lines = [
        "# HG002 v5.0q GRCh38 chr20 credibility report",
        "",
        f"- Status: **{status}**",
        f"- Receipt claim: `{latest['receipt_claim']}`",
        f"- Evaluator image: `{external['container']}`",
        (
            "- Predicted / realized RSS: "
            f"{memory['predicted_peak_rss_bytes']} / {memory['peak_rss_bytes']} bytes"
        ),
        "",
        "This is an honest SNV-focused baseline, not a competitive threshold or an indel-quality claim.",
    ]
    (results / "credibility.md").write_text("\n".join(lines) + "\n")
    return status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--latest", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    try:
        status = compare(args.baseline, args.latest, args.results)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"GIAB comparison failed: {error}")
        return 2
    print(status)
    return 1 if status == "diverged" else 0


if __name__ == "__main__":
    raise SystemExit(main())
