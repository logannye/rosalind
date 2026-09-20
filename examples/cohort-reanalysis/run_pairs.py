#!/usr/bin/env python3
"""Verify explicit paired reports against authored counts and Python Fraction.

Requires a passed run_cohort.py report. Retains commands and failures; synthetic
comparisons are engineering evidence, not independent research adoption.
"""
import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def readable_report(rows):
    """Render validated fixture rows without turning missing evidence into zero."""
    def side(row, prefix):
        if row[prefix + "_status"] == "unmeasured":
            return "unmeasured"
        depth = int(row[prefix + "_callable_depth"])
        counts = f"{row[prefix + '_alt_count']}/{depth}"
        if depth == 0:
            return counts + " (observed; fraction undefined)"
        if row[prefix + "_depth_eligible"] == "false":
            return counts + " (below depth screen)"
        return counts

    lines = ["# Synthetic paired candidate comparisons", "",
             "Direction: **right minus left** observed ALT/read-depth fraction. "
             "Pairs are explicitly supplied; metadata never infers pairing.", "",
             "The depth screen is 10 callable reads. It is a technical screen, "
             "not confidence or a biological response classification.", "",
             "| Pair (left → right) | Position / ALT | Left ALT/depth | Right ALT/depth | Exact difference | Both pass depth screen |",
             "|---|---|---|---|---|---|"]
    for row in rows:
        if row["difference_numerator"] == ".":
            difference = "undefined"
        else:
            difference = Fraction(int(row["difference_numerator"]), int(row["difference_denominator"]))
            if row["difference_negative"] == "true":
                difference = -difference
        eligible = {"true": "yes", "false": "no", ".": "unknown"}[row["both_depth_eligible"]]
        lines.append(f"| {row['pair_id']} ({row['left_member_id']} → {row['right_member_id']}) "
                     f"| {row['pos']} / {row['alt']} | {side(row, 'left')} | {side(row, 'right')} "
                     f"| {difference} | {eligible} |")
    lines += ["", "Unmeasured means no saved observation. Observed zero depth is measured, "
              "but has no defined ALT fraction. Positive low-depth fractions remain mathematically defined.", "",
              "The native TSV and its receipt preserve original integer counts and unreduced difference "
              "components. Fractions above are simplified exactly, without floating-point rounding.", "",
              "[Native paired rows](pairs.tsv) · [Verification and source identity](report.json)", "",
              "These are authored synthetic inputs, not independent research use or clinical validation."]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo-report", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "scope": "authored synthetic pairs; no adoption or clinical claim", "commands": []}
    report_path = args.output / "report.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    save()
    try:
        demo = json.loads(args.demo_report.read_text())
        require(demo["status"] == "passed", "the preceding cohort demonstration must pass")
        binary = args.binary.resolve(strict=True)
        binary_hash = digest(binary)
        require(binary_hash == demo["binary_sha256"], "paired demo must use the same binary as its cohort demonstration")
        report.update(binary_sha256=binary_hash, demo_report_sha256=digest(args.demo_report),
                      script_sha256=digest(Path(__file__)), pair_table_sha256=digest(HERE / "pairs.tsv"),
                      expected_sha256=digest(HERE / "expected.json"))
        cohort = Path(demo["cohort"])
        snapshot = demo["parent_snapshot"]
        candidates = args.demo_report.parent / "first.vcf"
        output = args.output / "pairs.tsv"
        command = [str(binary), "cohort", "compare-pairs", "--cohort", str(cohort), "--snapshot", snapshot,
                   "--pairs", str(HERE / "pairs.tsv"), "--sites", str(candidates), "--missing", "partial",
                   "--memory-budget-mb", "512", "--enforce", "--format", "tsv", "--output", str(output)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=180)
        report["commands"].append({"argv": command, "exit_code": result.returncode,
                                   "stdout": result.stdout, "stderr": result.stderr})
        save()
        require(result.returncode == 0, "native paired report failed")
        with output.open() as stream:
            actual = list(csv.DictReader(stream, delimiter="\t"))
        with (HERE / "pairs.tsv").open() as stream:
            pairs = list(csv.DictReader(stream, delimiter="\t"))
        expected = json.loads((HERE / "expected.json").read_text())
        positions = list(map(int, expected["candidates"]))
        require(len(actual) == len(pairs) * len(positions), "paired row cardinality differs")
        for row, (pair, position) in zip(actual, ((pair, pos) for pair in pairs for pos in positions)):
            require(row["pair_id"] == pair["id"] and int(row["pos"]) == position, "explicit pair/candidate order differs")
            fractions = []
            eligibilities = []
            depths = []
            for side in ("left", "right"):
                member = expected["members"][pair[side]]
                measured = position in member["stored_positions"]
                values = next(values for values in member["full_rows"] if values[0] == position)
                depth = values[3]
                alt = values[4 + "ACGT".index(expected["candidates"][str(position)])]
                require(row[side + "_member_id"] == pair[side], "side identity differs")
                require(row[side + "_status"] == ("observed" if measured else "unmeasured"), "side missingness differs")
                require(row[side + "_callable_depth"] == (str(depth) if measured else "."), "side depth differs")
                require(row[side + "_alt_count"] == (str(alt) if measured else "."), "side ALT count differs")
                eligible = depth >= 10 if measured else None
                require(row[side + "_depth_eligible"] == (str(eligible).lower() if measured else "."), "side eligibility differs")
                require(row[side + "_alt_supported"] == (str(eligible and alt > 0).lower() if measured else "."), "side support differs")
                require(row[side + "_alt_fraction_numerator"] == (str(alt) if measured and depth else "."), "side fraction numerator differs")
                require(row[side + "_alt_fraction_denominator"] == (str(depth) if measured and depth else "."), "side fraction denominator differs")
                fractions.append(Fraction(alt, depth) if measured and depth else None)
                eligibilities.append(eligible)
                depths.append(depth)
            both = None if None in eligibilities else all(eligibilities)
            require(row["both_depth_eligible"] == ("." if both is None else str(both).lower()), "joint eligibility differs")
            if None in fractions:
                require(all(row[key] == "." for key in ("difference_negative", "difference_numerator", "difference_denominator")),
                        "undefined observed fraction became a difference")
            else:
                difference = fractions[1] - fractions[0]
                observed = Fraction(int(row["difference_numerator"]), int(row["difference_denominator"]))
                if row["difference_negative"] == "true":
                    observed = -observed
                require(observed == difference, "exact difference differs from independent Fraction arithmetic")
                require(int(row["difference_denominator"]) == depths[0] * depths[1], "original depth denominator was not retained")
                require(row["difference_negative"] == str(difference < 0).lower(), "difference sign differs")
        receipt_path = Path(str(output) + ".manifest.json")
        receipt = json.loads(receipt_path.read_text())
        require(receipt["measurements"]["execution.original_alignment_records_decoded"] == "0", "paired saved-only query decoded original records")
        verify = [str(binary), "verify", "--manifest", str(receipt_path), "--json"]
        result = subprocess.run(verify, capture_output=True, text=True, timeout=180)
        report["commands"].append({"argv": verify, "exit_code": result.returncode, "stdout": result.stdout, "stderr": result.stderr})
        require(result.returncode == 0 and json.loads(result.stdout)["ok"], "paired receipt verification failed")
        require(digest(binary) == binary_hash, "binary changed during paired demonstration")
        readable = args.output / "REPORT.md"
        readable.write_text(readable_report(actual))
        report.update(status="passed", compared_rows=len(actual), output_sha256=digest(output),
                      receipt_sha256=digest(receipt_path), snapshot_id=snapshot,
                      readable_report_sha256=digest(readable),
                      source_build={key: receipt["params"].get(key) for key in
                                    ("code_git_sha", "code_dirty", "deps_lock_blake3", "target_triple")})
    except BaseException as error:
        report.update(status="failed", failure=f"{type(error).__name__}: {error}")
        raise
    finally:
        save()
    print(f"Verified exact explicit pair comparisons: {report_path}")


if __name__ == "__main__":
    main()
