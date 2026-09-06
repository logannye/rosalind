#!/usr/bin/env python3
"""Retain a three-repeat evidence budget curve and an independent pysam baseline."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(65536), b""):
            value.update(block)
    return value.hexdigest()


def run_measured(root, label, argv, output):
    raw = root / "raw"
    time_path = raw / (label + ".time.txt")
    arguments = ["/usr/bin/time", "-l" if platform.system() == "Darwin" else "-v",
                 "-o", str(time_path), *map(str, argv)]
    (raw / (label + ".argv.json")).write_text(json.dumps(list(map(str, argv)), indent=2) + "\n")
    started = time.perf_counter()
    with output.open("xb") as stdout, (raw / (label + ".stderr.txt")).open("xb") as stderr:
        code = subprocess.run(arguments, stdout=stdout, stderr=stderr, check=False).returncode
    wall = time.perf_counter() - started
    measured = time_path.read_text()
    if platform.system() == "Darwin":
        match = re.search(r"(\d+)\s+maximum resident set size", measured)
        peak = int(match[1]) if match else None
    else:
        match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", measured)
        peak = int(match[1]) * 1024 if match else None
    return {"argv": list(map(str, argv)), "exit_code": code, "wall_seconds": wall,
            "peak_rss_bytes": peak, "stdout_bytes": output.stat().st_size,
            "stdout_sha256": sha256(output), "raw_time": str(time_path.relative_to(root))}


def compare(left, right):
    different, total, examples = 0, 0, []
    with left.open() as a, right.open() as b:
        for number, (x, y) in enumerate(itertools.zip_longest(a, b), start=1):
            total += 1
            if x != y:
                different += 1
                if len(examples) < 5:
                    examples.append({"line": number, "left": x, "right": y})
    return {"lines_including_header": total, "different_lines": different,
            "equal": different == 0, "first_differences": examples}


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path,
                        help="prepared research-filter directory")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--budgets", default="96,128,256")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cache", action="store_true", help="also measure cold and resumed cache")
    args = parser.parse_args()
    if pysam.__version__ != "0.23.3":
        parser.error("baseline requires pysam==0.23.3")
    budgets = [int(value) for value in args.budgets.split(",")]
    if args.repeats < 1 or any(value <= 0 for value in budgets):
        parser.error("positive budgets and repeats required")
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.output.resolve()
    (root / "raw").mkdir()
    binary, data = args.binary.resolve(), args.data.resolve()
    reference, bam, sites = data / "ex1.fa", data / "sample.bam", data / "candidates.vcf"
    inputs = [reference, Path(str(reference) + ".fai"), bam, Path(str(bam) + ".bai"), sites]
    started = time.perf_counter()
    report = {"schema": 1, "scientific_profile": "shortread-dna-readcount-v1",
              "status": "running", "input_sha256": {str(path): sha256(path) for path in inputs},
              "environment": {"platform": platform.platform(), "machine": platform.machine(),
                              "python": sys.version, "pysam": pysam.__version__,
                              "htslib": pysam.__samtools_version__,
                              "binary_sha256": sha256(binary),
                              "rosalind_version": subprocess.check_output([str(binary), "--version"], text=True).strip()},
              "measurements": {}, "semantic_comparisons": {},
              "limits": ["No OS memory limit is imposed by this harness.",
                         "Rosalind declarations govern only its native process; pysam has no declared-budget contract.",
                         "The pysam baseline extracts equivalent rows but does not seal receipts or verify artifacts.",
                         "Kernel timings are not isolated: analysis_encoding_ms combines those phases.",
                         "This tiny candidate slice is a reproducibility check, not a whole-genome performance result.",
                         "A local source binary is not a published package; its exact SHA and receipt producer identity are retained."]}
    report["harness_input_hashing_seconds"] = time.perf_counter() - started
    preparation = data / "preparation.json"
    if preparation.is_file():
        report["preparation"] = json.loads(preparation.read_text())
    baseline = Path(__file__).with_name("pysam_evidence.py").resolve()
    for repeat in range(1, args.repeats + 1):
        label = f"pysam-{repeat}"
        output = root / (label + ".tsv")
        measured = run_measured(root, label, [sys.executable, baseline, reference, bam, sites], output)
        report["measurements"][label] = measured
    cases = [(f"rosalind-{budget}-{repeat}", budget, []) for budget in budgets
             for repeat in range(1, args.repeats + 1)]
    if args.cache:
        cases += [("rosalind-cache-cold", max(budgets), ["--cache-dir", str(root / "cache")]),
                  ("rosalind-cache-resumed", max(budgets), ["--cache-dir", str(root / "cache"), "--resume"])]
    for label, budget, extra in cases:
        output = root / (label + ".tsv")
        receipt = root / (label + ".manifest.json")
        # TSV travels through stdout for identical measurement plumbing. Persist
        # separately via the native output option to obtain byte-verifiable receipts.
        stdout = root / (label + ".stdout.txt")
        argv = [binary, "analyze", "evidence", "--reference", reference, "--alignments", bam,
                "--sites", sites, "--memory-budget-mb", budget, "--output", output,
                "--manifest", receipt, *extra]
        measured = run_measured(root, label, argv, stdout)
        measured["declared_budget_mib"] = budget
        measured["status"] = {0: "completed", 3: "refused", 4: "resource-failed"}.get(measured["exit_code"], "failed")
        if output.is_file():
            measured.update({"output_bytes": output.stat().st_size, "output_sha256": sha256(output)})
        if receipt.is_file():
            manifest = json.loads(receipt.read_text())
            measured["receipt_measurements"] = manifest.get("measurements", {})
            measured["producer"] = {key: value for key, value in manifest["params"].items()
                                    if key.startswith(("code_", "producer.", "rustc_", "target_", "deps_"))}
        report["measurements"][label] = measured
        if measured["exit_code"] == 0:
            verified = run_measured(root, label + "-verify", [binary, "verify", "--manifest", receipt, "--json"],
                                    root / (label + ".verify.json"))
            measured["verification"] = verified
            report["semantic_comparisons"][label] = compare(root / "pysam-1.tsv", output)
        (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    report["baseline_repeat_identical"] = all(compare(root / "pysam-1.tsv", root / f"pysam-{n}.tsv")["equal"]
                                               for n in range(1, args.repeats + 1))
    ok = report["baseline_repeat_identical"] and all(
        value["exit_code"] == 0 and value.get("verification", {"exit_code": 0})["exit_code"] == 0
        for value in report["measurements"].values()) and all(
        value["equal"] for value in report["semantic_comparisons"].values())
    report["status"] = "passed" if ok else "failed-or-refused"
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{report['status']}: {root / 'report.json'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
