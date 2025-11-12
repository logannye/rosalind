#!/usr/bin/env python3
"""Benchmark SAM vs BAM emission using the Rosalind CLI."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Dict, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_DATASET = PROJECT_ROOT / "examples" / "data" / "illumina_toy"
GENERATOR = SCRIPT_DIR / "generate_toy_data.py"


def ensure_dataset(dataset: Path) -> Tuple[Path, Path]:
    reference = dataset / "reference.fa"
    reads = dataset / "reads_R1.fastq"

    if reference.exists() and reads.exists():
        return reference, reads

    print(f"[benchmark] dataset {dataset} missing, generating via {GENERATOR.name} …")
    dataset.mkdir(parents=True, exist_ok=True)
    subprocess.run(["python3", str(GENERATOR), str(dataset)], check=True)
    return reference, reads


def run_command(cmd: list[str], env: Dict[str, str]) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.perf_counter() - start


def human_size(num_bytes: int) -> str:
    for unit in ["B", "KB", "MB", "GB"]:
        if num_bytes < 1024.0 or unit == "GB":
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.1f} GB"


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare SAM vs BAM output performance")
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET, help="Dataset directory produced by generate_toy_data.py")
    parser.add_argument("--iterations", type=int, default=3, help="Number of repetitions per format")
    parser.add_argument("--cargo", default="cargo", help="Path to cargo executable")
    parser.add_argument("--release", action="store_true", help="Use release build (recommended)")
    parser.add_argument("--env", action="append", default=[], help="Additional KEY=VALUE pairs passed to cargo run")
    args = parser.parse_args()

    reference, reads = ensure_dataset(args.dataset)

    cargo_cmd = [args.cargo, "run"]
    if args.release:
        cargo_cmd.append("--release")
    cargo_cmd.extend(["--", "align", "--reference", str(reference), "--reads", str(reads)])

    env = os.environ.copy()
    for assignment in args.env:
        if "=" not in assignment:
            parser.error(f"Invalid env assignment: {assignment}")
        key, value = assignment.split("=", 1)
        env[key] = value

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        results = []
        for format_name in ["sam", "bam"]:
            fmt_cmd = cargo_cmd + ["--format", format_name]
            output_file = tmpdir_path / f"benchmark_output.{format_name}"
            fmt_cmd.extend(["--output", str(output_file)])

            timings = []
            for _ in range(args.iterations):
                try:
                    duration = run_command(fmt_cmd, env)
                except subprocess.CalledProcessError as exc:
                    raise SystemExit(
                        f"Command {' '.join(fmt_cmd)} failed with exit code {exc.returncode}.\n"
                        "If the CLI cannot link against Python, set PYO3_PYTHON to a valid interpreter"
                        " or disable bindings with PYO3_NO_PYTHON=1."
                    ) from exc
                timings.append(duration)

            size = output_file.stat().st_size
            results.append(
                {
                    "format": format_name.upper(),
                    "avg_time": sum(timings) / len(timings),
                    "min_time": min(timings),
                    "max_time": max(timings),
                    "size": size,
                }
            )

    print("\nRosalind format benchmark")
    print("==========================")
    print(f"Reference: {reference}")
    print(f"Reads:     {reads}")
    print(f"Iterations per format: {args.iterations}")
    print()
    print(f"{'Format':<8}{'Avg (s)':>10}{'Min (s)':>10}{'Max (s)':>10}{'Output size':>16}")
    for row in results:
        print(
            f"{row['format']:<8}{row['avg_time']:>10.3f}{row['min_time']:>10.3f}{row['max_time']:>10.3f}{human_size(row['size']):>16}"
        )

    print("\nTip: run with --env PYO3_PYTHON=/path/to/python if the build toolchain needs an explicit interpreter.")


if __name__ == "__main__":
    main()
