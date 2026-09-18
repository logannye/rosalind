#!/usr/bin/env python3
"""Check authored expectations against SAM and existing single-sample commands.

No cohort runtime is implemented here. The partial/summary calculations are a
small independent fixture oracle, deliberately bounded to three members/four loci.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            value.update(block)
    return value.hexdigest()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_tsv(path):
    with path.open() as stream:
        return {int(row["pos"]): row for row in csv.DictReader(stream, delimiter="\t")}


def sam_oracle(path, positions):
    """Independent parser for this fixture's unpaired 1M reads, without htslib."""
    counts = {pos: [pos, 0, 0, 0, 0, 0, 0, 0] for pos in positions}
    with path.open() as stream:
        for line in stream:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            require(fields[5] == "1M" and len(fields[9]) == 1, "fixture oracle supports only 1M")
            position = int(fields[3])
            if position not in counts:
                continue
            row = counts[position]
            row[1] += 1
            flag, mapq = int(fields[1]), int(fields[4])
            if flag & (0x100 | 0x800 | 0x200 | 0x400) or mapq == 255 or mapq < 20:
                continue
            row[2] += 1
            base, quality = fields[9], ord(fields[10]) - 33
            if quality < 20 or base not in "ACGT":
                continue
            row[3] += 1
            row[4 + "ACGT".index(base)] += 1
    return counts


def reference_partial(expected):
    rows = []
    for specimen, member in sorted(expected["members"].items()):
        for row in member["full_rows"]:
            pos, _, _, depth, *alleles = row
            if pos not in member["stored_positions"]:
                rows.append([specimen, pos, "unmeasured", None, None, None, None, None, None])
                continue
            alt = expected["candidates"][str(pos)]
            count = alleles["ACGT".index(alt)]
            eligible = depth >= expected["thresholds"]["min_callable_depth"]
            supported = eligible and count >= 1
            rows.append([specimen, pos, "observed", depth, count, eligible, supported,
                         count if depth else None, depth if depth else None])
    return rows


def reference_summary(expected, complete):
    result = []
    for text, alt in expected["candidates"].items():
        pos = int(text)
        observed = [row for member in expected["members"].values()
                    for row in member["full_rows"]
                    if row[0] == pos and (complete or pos in member["stored_positions"])]
        eligible = [row for row in observed if row[3] >= expected["thresholds"]["min_callable_depth"]]
        offset = 4 + "ACGT".index(alt)
        supported = [row for row in eligible if row[offset] >= 1]
        result.append([pos, alt, len(expected["members"]), len(observed), len(eligible), len(supported),
                       sum(row[3] for row in observed), sum(row[offset] for row in observed),
                       sum(row[3] for row in eligible), sum(row[offset] for row in eligible)])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    prepared = args.prepared.resolve(strict=True)
    raw = prepared / "raw"
    output = prepared / "validation"
    output.mkdir(exist_ok=False)
    expected = json.loads((HERE / "expected.json").read_text())
    require(len(expected["members"]) == 3 and len(expected["candidates"]) == 4, "bounded fixture shape")
    require(reference_partial(expected) == expected["partial_rows"], "hand-authored partial expectations")
    require(reference_summary(expected, False) == expected["partial_summary"], "partial summary oracle")
    require(reference_summary(expected, True) == expected["after_fill_summary"], "filled summary oracle")
    commands = []

    def run(arguments, succeeds=True):
        invocation = [str(binary), *map(str, arguments)]
        process = subprocess.run(invocation, text=True, capture_output=True)
        commands.append({"argv": invocation, "exit_code": process.returncode})
        if succeeds:
            require(process.returncode == 0, f"command failed: {invocation}\n{process.stderr}")
        else:
            require(process.returncode != 0, f"expected refusal: {invocation}")
        return process

    def native(specimen, selection, destination, extra=()):
        return run(["analyze", "evidence", "--reference", raw / "reference.fa", "--alignments",
                    raw / (specimen + ".bam"), *selection, "--fields", "depths,alleles",
                    "--memory-budget-mb", "128", "--format", "tsv", "--output", destination, *extra])

    columns = expected["full_row_columns"]
    manifests = {}
    for specimen, member in expected["members"].items():
        parsed = sam_oracle(raw / (specimen + ".sam"), map(int, expected["candidates"]))
        require(list(parsed.values()) == member["full_rows"], f"independent SAM oracle: {specimen}")
        path = output / (specimen + "-stored.tsv")
        native(specimen, ["--regions", HERE / (specimen + ".bed")], path,
               ["--cache-dir", output / (specimen + "-cache")])
        actual = read_tsv(path)
        require(sorted(actual) == member["stored_positions"], f"stored coverage: {specimen}")
        for row in member["full_rows"]:
            if row[0] in actual:
                require([int(actual[row[0]][column]) for column in columns] == row, f"counts: {specimen}/{row[0]}")
        if specimen == "specimen-b":
            for counter, value in expected["specimen_b_position_10_exclusions"].items():
                require(int(actual[10][counter]) == value, f"filter counter: {counter}")
        receipt = json.loads(Path(str(path) + ".manifest.json").read_text())
        source = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
        destination = output / "portable" / specimen
        shutil.copytree(source.parent, destination)
        manifests[specimen] = destination / source.name

    unavailable = prepared / "raw-unavailable"
    raw.rename(unavailable)
    try:
        for specimen, manifest in manifests.items():
            run(["dataset", "verify", "--dataset", manifest])
            offline = output / (specimen + "-offline.tsv")
            run(["dataset", "extract", "--dataset", manifest, "--regions", HERE / (specimen + ".bed"),
                 "--fields", "depths,alleles", "--format", "tsv", "--output", offline])
            require(offline.read_bytes() == (output / (specimen + "-stored.tsv")).read_bytes(), "offline equality")
            if specimen != "specimen-a":
                refused = output / (specimen + "-missing.tsv")
                result = run(["dataset", "extract", "--dataset", manifest, "--sites", HERE / "candidates.vcf",
                              "--format", "tsv", "--output", refused], succeeds=False)
                require("absent" in result.stderr and not refused.exists(), "strict missing-locus refusal")
        missing_fields = output / "missing-fields.tsv"
        result = run(["dataset", "extract", "--dataset", manifests["specimen-a"], "--fields", "strands",
                      "--format", "tsv", "--output", missing_fields], succeeds=False)
        require("fields" in result.stderr and not missing_fields.exists(), "missing-field capability refusal")
    finally:
        unavailable.rename(raw)

    reuse_counts = {}
    for specimen, member in expected["members"].items():
        fresh, reused = output / (specimen + "-fresh.tsv"), output / (specimen + "-reused.tsv")
        native(specimen, ["--sites", HERE / "candidates.vcf"], fresh)
        native(specimen, ["--sites", HERE / "candidates.vcf"], reused,
               ["--reuse-dataset", manifests[specimen], "--workers", "1"])
        require(fresh.read_bytes() == reused.read_bytes(), f"fresh/reuse equality: {specimen}")
        actual = read_tsv(fresh)
        for row in member["full_rows"]:
            require([int(actual[row[0]][column]) for column in columns] == row, "filled observations")
        receipt = json.loads(Path(str(reused) + ".manifest.json").read_text())
        stats = receipt["measurements"]
        reused_count, computed_count = int(stats["execution.reused_loci"]), int(stats["execution.computed_loci"])
        require((reused_count, computed_count) == (len(member["stored_positions"]), 4 - len(member["stored_positions"])),
                "reuse and computation denominators")
        reuse_counts[specimen] = {"reused_loci": reused_count, "computed_loci": computed_count}
        run(["verify", "--manifest", str(reused) + ".manifest.json"])

    report = {
        "fixture": "cohort-contract-synthetic-v1", "passed": True,
        "scope": "Existing single-sample extraction/dataset/reuse checks; no cohort implementation or user validation.",
        "binary_sha256": digest(binary),
        "binary_version": subprocess.check_output([str(binary), "--version"], text=True).strip(),
        "expectations_sha256": digest(HERE / "expected.json"),
        "checks": ["independent-sam-counts", "hand-authored-missingness-and-summary", "exclusive-filter-counters",
                   "relocated-offline-byte-equality", "strict-missing-refusal", "missing-fields-refusal",
                   "fresh-reuse-byte-equality", "reuse-denominators", "reused-artifact-verification"],
        "reuse": reuse_counts, "commands": commands,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Synthetic fixture checks passed; report: {output / 'report.json'}")


if __name__ == "__main__":
    main()
