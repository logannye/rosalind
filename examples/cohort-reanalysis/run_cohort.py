#!/usr/bin/env python3
"""Exercise the source-preview cohort workflow and record its actual costs.

Uses only authored synthetic inputs prepared by prepare.py. All output paths are
create-new; original inputs are moved aside and restored, never deleted. This is
an engineering demonstration, not independent adoption or a performance claim.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import time

from validate import digest, reference_partial, reference_summary, require, sam_oracle

HERE = Path(__file__).resolve().parent


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def size(directory):
    return sum(path.stat().st_size for path in directory.rglob("*") if path.is_file())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--budget-mb", type=int, default=512)
    args = parser.parse_args()
    prepared = args.prepared.resolve(strict=True)
    binary = args.binary.resolve(strict=True)
    def binary_stamp():
        stat = binary.stat()
        return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    initial_binary_stamp = binary_stamp()
    raw = prepared / "raw"
    output = prepared / "cohort-validation"
    output.mkdir(exist_ok=False)
    report = {
        "schema": 1, "status": "running", "scope": "three authored synthetic samples; no performance or adoption claim",
        "binary": str(binary), "binary_sha256": digest(binary),
        "binary_version": subprocess.check_output([str(binary), "--version"], text=True).strip(),
        "budget_mb": args.budget_mb, "commands": [], "costs": {}, "comparisons": {},
        "fixture_sha256": digest(HERE / "expected.json"),
        "preparation": json.loads((prepared / "preparation.json").read_text()),
        "script_sha256": digest(Path(__file__)),
    }
    report_path = output / "report.json"

    def save_report():
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    def run(label, arguments, succeeds=True):
        require(binary_stamp() == initial_binary_stamp, "binary changed during demonstration; use a stable copy")
        invocation = [str(binary), *map(str, arguments)]
        start = time.perf_counter_ns()
        result = subprocess.run(invocation, text=True, capture_output=True, cwd=output)
        entry = {"label": label, "argv": invocation, "expected_success": succeeds, "exit_code": result.returncode,
                 "elapsed_ns": time.perf_counter_ns() - start,
                 "stdout": result.stdout, "stderr": result.stderr}
        if "--output" in invocation:
            artifact = Path(invocation[invocation.index("--output") + 1])
            receipt_path = Path(str(artifact) + ".manifest.json")
            if result.returncode == 0 and receipt_path.exists():
                receipt = json.loads(receipt_path.read_text())
                entry["artifact_bytes"] = artifact.stat().st_size
                entry["receipt_bytes"] = receipt_path.stat().st_size
                entry["measurements"] = receipt.get("measurements", {})
                entry["source_build"] = {key: receipt["params"].get(key)
                    for key in ("code_git_sha", "code_dirty", "deps_lock_blake3", "target_triple")}
        report["commands"].append(entry)
        save_report()
        require(binary_stamp() == initial_binary_stamp, "binary changed during demonstration; use a stable copy")
        require((result.returncode == 0) == succeeds, f"{label}: unexpected exit {result.returncode}\n{result.stderr}")
        return result

    resources = ["--memory-budget-mb", args.budget_mb, "--enforce"]
    cohort = output / "cohort"
    expected = json.loads((HERE / "expected.json").read_text())
    candidates = output / "first.vcf"
    shutil.copyfile(HERE / "candidates.vcf", candidates)
    second = output / "second.vcf"
    second.write_text("##fileformat=VCFv4.3\n##contig=<ID=synthetic,length=64>\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                      "synthetic\t10\t.\tA\tG\t.\tPASS\t.\n"
                      "synthetic\t20\t.\tA\tG\t.\tPASS\t.\n")
    moved_raw = prepared / "raw-cohort-demo-unavailable"
    moved_leaves = output / "leaves-unavailable"
    leaves = output / "leaves"
    leaves.mkdir()
    manifests = {}
    snapshot = None

    def query(action, snapshot_id, selection, destination=None, partial=False, extra=()):
        command = ["cohort", action, "--cohort", cohort, "--snapshot", snapshot_id,
                   "--sites", selection, *resources]
        if partial:
            command += ["--missing", "partial"]
        if destination is not None:
            command += ["--format", "tsv", "--output", destination]
        return command + list(extra)

    def check_partial(path):
        actual = []
        for row in rows(path):
            values = []
            for key in ("callable_depth", "alt_count", "depth_eligible", "alt_supported",
                        "observed_alt_fraction_numerator", "observed_alt_fraction_denominator"):
                value = row[key]
                values.append(None if value == "." else value == "true" if key in
                              ("depth_eligible", "alt_supported") else int(value))
            actual.append([row["sample_id"], int(row["pos"]), row["status"], *values])
        require(actual == expected["partial_rows"], "partial cells differ from hand-authored oracle")

    def check_summary(path, complete):
        counts = ("n_requested", "n_observed", "n_depth_eligible", "n_alt_supported",
                  "callable_total", "alt_total", "eligible_callable_total", "eligible_alt_total")
        actual = [[int(row["pos"]), row["alt"], *(int(row[key]) for key in counts)] for row in rows(path)]
        require(actual == reference_summary(expected, complete), "candidate summary differs from independent arithmetic")

    def compare_fresh(cohort_path, specimen, fresh_path):
        """Compare every shared evidence integer, not just the ALT/depth screen."""
        actual = [row for row in rows(cohort_path) if row["sample_id"] == specimen]
        fresh = rows(fresh_path)
        require(len(actual) == len(fresh), "fresh/cohort candidate cardinality differs")
        for observed, native in zip(actual, fresh):
            require(observed["status"] == "observed", "fresh comparison requires measured cells")
            require(observed["contig"] == native["#contig"] and observed["alt"] == native["requested_alts"],
                    "fresh/cohort candidate identity differs")
            for key, value in native.items():
                if key not in ("#contig", "requested_alts"):
                    require(observed[key] == value, f"fresh/cohort {specimen}:{observed['pos']} {key} differs")
        return len(actual)

    def elapsed(label):
        return next(command["elapsed_ns"] for command in report["commands"] if command["label"] == label)

    try:
        require(args.budget_mb > 0, "budget must be positive")
        require(not moved_raw.exists() and not moved_leaves.exists(), "unavailable-input paths must be new")
        for name, expected_hash in report["preparation"]["prepared_sha256"].items():
            require(digest(raw / name) == expected_hash, f"prepared input identity changed: {name}")
        for name, expected_hash in report["preparation"]["source_sha256"].items():
            require(digest(HERE / name) == expected_hash, f"fixture source identity changed: {name}")
        report["costs"]["original_input_bytes"] = size(raw)
        require(reference_partial(expected) == expected["partial_rows"], "authored missingness oracle")
        for specimen, member in expected["members"].items():
            require(list(sam_oracle(raw / (specimen + ".sam"), map(int, expected["candidates"])).values())
                    == member["full_rows"], "independent single-base SAM oracle")
            stored = output / (specimen + "-stored.tsv")
            run("first-extraction:" + specimen,
                ["analyze", "evidence", "--alignments", raw / (specimen + ".bam"),
                 "--reference", raw / "reference.fa", "--regions", HERE / (specimen + ".bed"),
                 "--fields", "depths,alleles", "--cache-dir", leaves / specimen,
                 "--format", "tsv", "--output", stored, *resources])
            receipt = json.loads(Path(str(stored) + ".manifest.json").read_text())
            manifests[specimen] = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
        members = output / "import.tsv"
        members.write_text("id\tmanifest\tgroup\n" + "".join(
            f"{member}\t{manifest}\tsynthetic\n" for member, manifest in sorted(manifests.items())))
        report["costs"]["original_dataset_bytes"] = size(leaves)
        run("import-plan", ["cohort", "create", "--cohort", cohort, "--members", members, "--plan", *resources])
        created = run("import-copy-and-verify", ["cohort", "create", "--cohort", cohort, "--members", members, *resources])
        snapshot = json.loads(created.stdout)["snapshot_id"]
        report["parent_snapshot"] = snapshot
        report["costs"]["imported_cohort_bytes"] = size(cohort)
        run("full-verification", ["cohort", "verify", "--cohort", cohort, "--snapshot", snapshot, *resources])
        plan = json.loads(run("strict-plan", query("extract", snapshot, candidates, extra=["--plan"])).stdout)
        require(plan["status"] == "blocked", "missing loci must be visible in strict plan")
        refused = output / "strict-refused.tsv"
        run("strict-refusal", query("extract", snapshot, candidates, refused), succeeds=False)
        require(not refused.exists(), "strict refusal published an output")
        partial = output / "partial.tsv"
        run("first-candidate-list", query("extract", snapshot, candidates, partial, True))
        check_partial(partial)
        summary = output / "partial-summary.tsv"
        run("first-summary", query("summarize", snapshot, candidates, summary, True))
        check_summary(summary, False)

        # Physical relocation plus unavailable original alignments AND imported
        # dataset locations tests the portable artifact's complete handoff.
        relocated = output / "relocated-cohort"
        start = time.perf_counter_ns()
        shutil.copytree(cohort, relocated)
        report["costs"]["relocation_copy_ns"] = time.perf_counter_ns() - start
        shutil.rmtree(cohort)  # Only this demonstration's generated copy.
        cohort = relocated
        raw.rename(moved_raw)
        leaves.rename(moved_leaves)
        run("relocated-full-verification", ["cohort", "verify", "--cohort", cohort, "--snapshot", snapshot, *resources])
        repeated = output / "partial-relocated.tsv"
        run("repeat-first-list-offline", query("extract", snapshot, candidates, repeated, True))
        require(repeated.read_bytes() == partial.read_bytes(), "relocation changed primary bytes")
        second_result = output / "second.tsv"
        run("second-question-offline", query("extract", snapshot, second, second_result))
        for row in rows(second_result):
            member = expected["members"][row["sample_id"]]
            original = next(value for value in member["full_rows"] if value[0] == int(row["pos"]))
            require(int(row["callable_depth"]) == original[3] and int(row["alt_count"]) == original[6],
                    "new ALT query differs from stored independent base counts")
        run("replay-offline", ["reproduce", "--manifest", Path(str(second_result) + ".manifest.json"),
                               "--inputs", output, "--binary", binary])
        run("verify-second-offline", ["verify", "--manifest", Path(str(second_result) + ".manifest.json")])
        moved_raw.rename(raw)
        moved_leaves.rename(leaves)

        # Compare matched second-question raw extraction, with startup, source
        # hashing and output finalization included in the recorded wall time.
        compared = 0
        for specimen in expected["members"]:
            fresh = output / (specimen + "-second-fresh.tsv")
            run("second-question-fresh:" + specimen,
                ["analyze", "evidence", "--alignments", raw / (specimen + ".bam"),
                 "--reference", raw / "reference.fa", "--sites", second, "--fields", "depths,alleles",
                 "--format", "tsv", "--output", fresh, *resources])
            compared += compare_fresh(second_result, specimen, fresh)
        report["comparisons"]["second_question"] = {
            "matched_sample_candidate_rows": compared,
            "evidence_fields": "all native depths/alleles columns, including exclusive filter counters",
            "saved_one_process_elapsed_ns": elapsed("second-question-offline"),
            "fresh_three_serial_processes_elapsed_ns": sum(elapsed("second-question-fresh:" + member)
                                                           for member in expected["members"]),
            "caveat": "one tiny ordered run; startup, hashing, verification and finalization included; no repeated trials or speedup claim",
        }
        sources = output / "sources.tsv"
        sources.write_text("id\trole\tpath\n" + "".join(
            f"{member}\t{role}\t{path}\n" for member in ("specimen-b", "specimen-c")
            for role, path in (("alignments", raw / (member + ".bam")),
                               ("alignment-index", raw / (member + ".bam.bai")),
                               ("reference", raw / "reference.fa"),
                               ("reference-fai", raw / "reference.fa.fai"))))
        extended = json.loads(run("extend-missing-only", query("extend", snapshot, candidates,
            extra=["--sources", sources, "--work-dir", output])).stdout)
        child = extended["snapshot_id"]
        require(child != snapshot and extended["parent"] == snapshot, "extension must publish a distinct child")
        require({member["id"]: (member["retained_loci"], member["computed_loci"])
                 for member in extended["members"]}
                == {"specimen-a": (4, 0), "specimen-b": (3, 1), "specimen-c": (3, 1)},
                "extension reused/computed locus counts differ")
        report["extension"] = extended
        report["costs"]["extended_cohort_bytes"] = size(cohort)
        complete = output / "complete.tsv"
        run("extended-query", query("extract", child, candidates, complete))
        for row in rows(complete):
            original = next(value for value in expected["members"][row["sample_id"]]["full_rows"]
                            if value[0] == int(row["pos"]))
            require(row["status"] == "observed" and [int(row[column]) for column in expected["full_row_columns"]]
                    == original, "extended observations disagree with native-source oracle")
        compared = 0
        for specimen in expected["members"]:
            fresh = output / (specimen + "-complete-fresh.tsv")
            run("complete-fresh:" + specimen,
                ["analyze", "evidence", "--alignments", raw / (specimen + ".bam"),
                 "--reference", raw / "reference.fa", "--sites", candidates, "--fields", "depths,alleles",
                 "--format", "tsv", "--output", fresh, *resources])
            compared += compare_fresh(complete, specimen, fresh)
        report["comparisons"]["extended_snapshot"] = {"matched_sample_candidate_rows": compared,
            "evidence_fields": "all native depths/alleles columns, including exclusive filter counters",
            "saved_one_process_elapsed_ns": elapsed("extended-query"),
            "fresh_three_serial_processes_elapsed_ns": sum(elapsed("complete-fresh:" + member)
                                                           for member in expected["members"])}
        complete_summary = output / "complete-summary.tsv"
        run("extended-summary", query("summarize", child, candidates, complete_summary))
        check_summary(complete_summary, True)
        unchanged = output / "parent-still-partial.tsv"
        run("unchanged-parent", query("extract", snapshot, candidates, unchanged, True))
        require(unchanged.read_bytes() == partial.read_bytes(), "extension modified original snapshot results")
        raw.rename(moved_raw)
        leaves.rename(moved_leaves)
        run("verify-extended-offline", ["cohort", "verify", "--cohort", cohort, "--snapshot", child, *resources])
        complete_offline = output / "complete-offline.tsv"
        run("extended-query-offline", query("extract", child, candidates, complete_offline))
        require(complete_offline.read_bytes() == complete.read_bytes(), "child needs original source files")
        run("replay-extended-offline", ["reproduce", "--manifest", Path(str(complete_offline) + ".manifest.json"),
                                       "--inputs", output, "--binary", binary])
        for command in report["commands"]:
            if command["argv"][1] == "cohort" and "measurements" in command:
                require(command["measurements"]["execution.original_alignment_records_decoded"] == "0",
                        "saved-only query reported original alignment decoding")
        report["costs"]["first_extraction_total_ns"] = sum(elapsed("first-extraction:" + member)
                                                         for member in expected["members"])
        report["costs"]["import_copy_and_verify_ns"] = elapsed("import-copy-and-verify")
        report["costs"]["full_verification_ns"] = elapsed("full-verification")
        report["child_snapshot"] = child
        report["cohort"] = str(cohort)
        report["primary_outputs_sha256"] = {path.name: digest(path) for path in
            (partial, summary, second_result, complete, complete_summary)}
        require(digest(binary) == report["binary_sha256"], "binary bytes changed during demonstration")
        report["status"] = "passed"
        lines = ["# Synthetic candidate reanalysis", "", "Read observations; the depth screen is not confidence.", "",
                 "| Sample | Position / ALT | Status | Callable depth | ALT observations |", "|---|---|---|---:|---:|"]
        for row in rows(partial):
            lines.append(f"| {row['sample_id']} | {row['pos']} / {row['alt']} | {row['status']} | {row['callable_depth']} | {row['alt_count']} |")
        lines += ["", "Full commands, failures, hashing/verification/import costs and source identity are in report.json.",
                  "This tiny fixture is dominated by process startup; it does not establish an economic break-even point."]
        (output / "REPORT.md").write_text("\n".join(lines) + "\n")
    except BaseException as error:
        report["status"] = "failed"
        report["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        if moved_raw.exists() and not raw.exists():
            moved_raw.rename(raw)
        if moved_leaves.exists() and not leaves.exists():
            moved_leaves.rename(leaves)
        save_report()
    print(f"Verified cohort workflow and measured costs: {report_path}")


if __name__ == "__main__":
    main()
