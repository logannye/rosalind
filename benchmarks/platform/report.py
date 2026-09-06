#!/usr/bin/env python3
"""Derive a platform report solely from retained command outputs and measurements."""

import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(os.environ.get("RESULTS", "/results"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tool_version(label: str, argv: list[str]) -> str:
    # Some packaged tools include arbitrary compiler/build bytes after their
    # version header. Retain them exactly; only the header is UTF-8 text.
    output = subprocess.check_output(argv)
    (ROOT / "raw" / f"version-{label}.stdout").write_bytes(output)
    (ROOT / "raw" / f"version-{label}.argv.json").write_text(json.dumps(argv) + "\n")
    return output.splitlines()[0].decode("utf-8").strip()


def measurement(label: str) -> dict:
    text = (ROOT / "raw" / f"{label}.time.txt").read_text()
    def value(pattern: str):
        found = re.search(pattern, text, re.MULTILINE)
        return float(found.group(1)) if found else None
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ((?:\d+:)?[0-9.]+)$", text, re.MULTILINE)
    wall_seconds = None
    if elapsed:
        parts = [float(part) for part in elapsed.group(1).split(":")]
        wall_seconds = sum(value * (60 ** index) for index, value in enumerate(reversed(parts)))
    output = ROOT / "outputs" / f"{label}.tsv"
    return {
        "argv": json.loads((ROOT / "raw" / f"{label}.argv.json").read_text()),
        "exit_code": int((ROOT / "raw" / f"{label}.exit").read_text()),
        "wall_seconds": wall_seconds,
        "peak_rss_kib": value(r"Maximum resident set size \(kbytes\): (\d+)$"),
        "output_bytes": output.stat().st_size if output.exists() else 0,
        "sha256": sha256(output) if output.exists() else None,
    }


def rows(path: Path, contig_order: dict):
    """Yield strict reference-order rows using constant memory."""
    with path.open(newline="", encoding="utf-8") as stream:
        previous = None
        for row in csv.reader(stream, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            key = (contig_order[row[0]], int(row[1]))
            if previous is not None and key <= previous:
                raise ValueError(f"{path}: duplicated or out-of-order position {row[:2]}")
            previous = key
            yield key, row


def compare(left: Path, right: Path, columns: list[int], contig_order: dict) -> dict:
    result = {
        "left_rows": 0, "right_rows": 0, "common_rows": 0,
        "left_only": 0, "right_only": 0, "different_common_rows": 0,
        "compared_column_indices": columns,
    }
    left_iter, right_iter = rows(left, contig_order), rows(right, contig_order)
    a, b = next(left_iter, None), next(right_iter, None)
    while a is not None or b is not None:
        if b is None or (a is not None and a[0] < b[0]):
            result["left_rows"] += 1
            result["left_only"] += 1
            a = next(left_iter, None)
        elif a is None or b[0] < a[0]:
            result["right_rows"] += 1
            result["right_only"] += 1
            b = next(right_iter, None)
        else:
            result["left_rows"] += 1
            result["right_rows"] += 1
            result["common_rows"] += 1
            result["different_common_rows"] += any(a[1][column] != b[1][column] for column in columns)
            a, b = next(left_iter, None), next(right_iter, None)
    return result


def main() -> None:
    with (ROOT / "inputs/reference.fa.fai").open() as stream:
        contig_order = {line.split("\t", 1)[0]: index for index, line in enumerate(stream)}
    labels = [f"{implementation}-{run}" for implementation in ("rosalind", "pysam", "bcftools") for run in range(1, 4)]
    measurements = {label: measurement(label) for label in labels}
    environment = {
        "rosalind": tool_version("rosalind", ["rosalind", "--version"]),
        "python": tool_version("python", ["python3", "--version"]),
        "pyarrow": tool_version("pyarrow", ["python3", "-c", "import pyarrow; print(pyarrow.__version__)"]),
        "pysam": tool_version("pysam", ["python3", "-c", "import pysam; print(pysam.__version__)"]),
        "bcftools": tool_version("bcftools", ["bcftools", "--version"]),
        "samtools": tool_version("samtools", ["samtools", "--version"]),
    }
    report = {
        "schema": 1,
        "scope": "resource predictability and reproducibility; not a speed or accuracy superiority claim",
        "environment": environment,
        "setup": {
            "common": ["build the digest-pinned environment", "mount one reference and BAM", "prepare coordinate-sorted BAM+BAI"],
            "rosalind": ["reference build", "plan --json before BAM preparation", "features", "verify"],
            "pysam": ["install exact locked Python dependency", "run maintained pysam_features.py"],
            "bcftools": ["install exact Debian bcftools", "run mpileup/query normalizer"],
        },
        "measurements": measurements,
        "repeat_identity": {
            implementation: len({measurements[f"{implementation}-{run}"]["sha256"] for run in range(1, 4)}) == 1
            for implementation in ("rosalind", "pysam", "bcftools")
        },
        "semantic_comparison": {
            "rosalind_vs_pysam_all_fields": compare(ROOT / "outputs/rosalind-1.tsv", ROOT / "outputs/pysam-1.tsv", list(range(2, 19)), contig_order),
            "rosalind_vs_bcftools_ref_and_depth": compare(ROOT / "outputs/rosalind-1.tsv", ROOT / "outputs/bcftools-1.tsv", [2, 3], contig_order),
            "bcftools_unmatched_fields": ["A/C/G/T canonical counts", "strand counts", "mean base quality", "mean mapping quality"],
            "bcftools_policy": "default mpileup filtering and depth sampling differ; recorded differences are not an equivalence certificate",
        },
        "rosalind_contract": json.loads((ROOT / "raw/rosalind-plan.json").read_text()),
        "rosalind_verify": json.loads((ROOT / "raw/rosalind-verify.json").read_text()),
        "arrow": json.loads((ROOT / "raw/arrow-checks.json").read_text()),
        "second_machine_reproduction": {
            "status": "not-run-by-this-harness",
            "command": "benchmarks/platform/run.sh REFERENCE_FASTA COORDINATE_SORTED_BAM OUTPUT_DIR",
            "compare": ["environment-packages.tsv", "environment-python.txt", "report.json", "raw/*.argv.json", "outputs/*"],
        },
    }
    (ROOT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"report": str(ROOT / "report.json"), "repeat_identity": report["repeat_identity"]}, sort_keys=True))

    if not all(report["repeat_identity"].values()):
        raise SystemExit("repeated feature output differs; see report.json")


if __name__ == "__main__":
    main()
