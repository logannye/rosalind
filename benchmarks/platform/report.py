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


def rows(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as stream:
        return {(row[0], row[1]): row for row in csv.reader(stream, delimiter="\t") if row and not row[0].startswith("#")}


def compare(left: Path, right: Path, columns: list[int]) -> dict:
    a, b = rows(left), rows(right)
    common = sorted(a.keys() & b.keys())
    different = sum(any(a[key][column] != b[key][column] for column in columns) for key in common)
    return {
        "left_rows": len(a), "right_rows": len(b), "common_rows": len(common),
        "left_only": len(a.keys() - b.keys()), "right_only": len(b.keys() - a.keys()),
        "different_common_rows": different, "compared_column_indices": columns,
    }


labels = [f"{implementation}-{run}" for implementation in ("rosalind", "pysam", "bcftools") for run in range(1, 4)]
measurements = {label: measurement(label) for label in labels}
environment = {
    "rosalind": subprocess.check_output(["rosalind", "--version"], text=True).strip(),
    "python": subprocess.check_output(["python3", "--version"], text=True).strip(),
    "pyarrow": subprocess.check_output(["python3", "-c", "import pyarrow; print(pyarrow.__version__)"], text=True).strip(),
    "pysam": subprocess.check_output(["python3", "-c", "import pysam; print(pysam.__version__)"], text=True).strip(),
    "bcftools": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
    "samtools": subprocess.check_output(["samtools", "--version"], text=True).splitlines()[0],
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
        "rosalind_vs_pysam_all_fields": compare(ROOT / "outputs/rosalind-1.tsv", ROOT / "outputs/pysam-1.tsv", list(range(2, 19))),
        "rosalind_vs_bcftools_equivalent_fields_only": compare(ROOT / "outputs/rosalind-1.tsv", ROOT / "outputs/bcftools-1.tsv", [2, 3]),
        "bcftools_unmatched_fields": ["A/C/G/T canonical counts", "strand counts", "mean base quality", "mean mapping quality"],
    },
    "rosalind_contract": json.loads((ROOT / "raw/rosalind-plan.json").read_text()),
    "rosalind_verify": json.loads((ROOT / "raw/rosalind-verify.json").read_text()),
    "arrow": json.loads((ROOT / "raw/arrow-checks.json").read_text()),
    "second_machine_reproduction": {
        "command": "benchmarks/platform/run.sh REFERENCE_FASTA COORDINATE_SORTED_BAM OUTPUT_DIR",
        "compare": ["environment-packages.tsv", "environment-python.txt", "report.json", "raw/*.argv.json", "outputs/*"],
    },
}
(ROOT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
print(json.dumps({"report": str(ROOT / "report.json"), "repeat_identity": report["repeat_identity"]}, sort_keys=True))
