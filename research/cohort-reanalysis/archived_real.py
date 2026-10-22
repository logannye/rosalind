#!/usr/bin/env python3
"""Reanalyze one pinned historical HG002 leaf against its retained independent oracle.

This is current saved-evidence validation, not new raw extraction, a fresh oracle
run, a throughput benchmark, or a second biological sample. No downloads/builds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARCHIVE = ROOT / "docs/findings/cram-decoder-validation-2026-09-06/verified-9b8e12f-2026-09-18/matrix-raw.tar.xz"
ARCHIVE_SHA256 = "13a56da1d4281c24dfba1833d41336ec623c470e86a87aabe69feb34d91db464"
LEAF_PREFIX = "matrix/cache/window-bam-sites-workers8-cache-r1/42b2bd2601ffb084eecfff835f6a9bcebf302d81861a5e62bd707c168b8a5863/"
ORACLE = "matrix/artifacts/window-bam-sites-oracle-r1.tsv"
ORACLE_METADATA = "matrix/artifacts/window-bam-sites-oracle-r1.oracle.json"
HISTORICAL_RECEIPT = "matrix/artifacts/window-bam-sites-workers8-cache-r1-cold.tsv.manifest.json"
MAX_COMPRESSED_BYTES = 2 << 20
MAX_SCANNED_BYTES = 400 << 20
MAX_EXTRACTED_BYTES = 16 << 20
MAX_ENTRY_BYTES = 2 << 20
MAX_ENTRIES = 10000
THRESHOLDS = (5, 10, 20)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def extract_selected(archive, output):
    """Stream one allowlisted leaf; never invoke extract/extractall or follow links."""
    if archive.stat().st_size > MAX_COMPRESSED_BYTES or digest(archive) != ARCHIVE_SHA256:
        raise ValueError("archive differs from the pinned repository evidence bundle")
    selected, scanned_bytes, extracted_bytes, entry_count = {}, 0, 0, 0
    special = {ORACLE: "oracle.tsv", ORACLE_METADATA: "oracle-metadata.json",
               HISTORICAL_RECEIPT: "historical-extraction-receipt.json"}
    with tarfile.open(archive, "r|xz") as source:
        for member in source:
            entry_count += 1
            scanned_bytes += member.size
            if entry_count > MAX_ENTRIES or scanned_bytes > MAX_SCANNED_BYTES:
                raise ValueError("archive scan bound exceeded")
            if member.name.startswith(LEAF_PREFIX):
                suffix = member.name[len(LEAF_PREFIX):]
                if member.isdir():
                    continue
                destination = PurePosixPath("leaf") / suffix
            elif member.name in special:
                destination = PurePosixPath(special[member.name])
            else:
                continue
            if (not member.isfile() or destination.is_absolute() or ".." in destination.parts
                    or member.size < 0 or member.size > MAX_ENTRY_BYTES
                    or str(destination) in selected):
                raise ValueError("unsafe, duplicate, or oversized selected archive member")
            extracted_bytes += member.size
            if extracted_bytes > MAX_EXTRACTED_BYTES:
                raise ValueError("selected archive output bound exceeded")
            path = output.joinpath(*destination.parts)
            path.parent.mkdir(parents=True, exist_ok=True)
            stream = source.extractfile(member)
            if stream is None:
                raise ValueError("regular selected archive entry has no body")
            remaining = member.size
            with path.open("xb") as target:
                while remaining:
                    block = stream.read(min(1 << 20, remaining))
                    if not block:
                        raise ValueError("truncated selected archive body")
                    target.write(block)
                    remaining -= len(block)
            selected[str(destination)] = {"archive_member": member.name, "bytes": member.size,
                                           "sha256": digest(path)}
    required = {"leaf/dataset.descriptor.json", "leaf/evidence-dataset.manifest.json", *special.values()}
    if not required <= set(selected):
        raise ValueError("archive lacks required leaf or independent oracle evidence")
    if digest(archive) != ARCHIVE_SHA256:
        raise ValueError("archive changed during extraction")
    return {"archive_sha256": ARCHIVE_SHA256, "archive_bytes": archive.stat().st_size,
            "archive_repository_path": str(ARCHIVE.relative_to(ROOT)), "selected": selected,
            "extracted_bytes": extracted_bytes, "scanned_entry_count": entry_count,
            "declared_archive_bytes_scanned": scanned_bytes}


def prepare(output):
    inventory = extract_selected(ARCHIVE, output)
    descriptor = json.loads((output / "leaf/dataset.descriptor.json").read_text())
    scope = descriptor["sample_scope"]
    if descriptor["version"] != 1 or scope["mode"] != "named" or scope["selected_sample"] != "HG002":
        raise ValueError("pinned leaf is not the historical named HG002 schema-v1 dataset")
    contigs = descriptor["contigs"]
    candidates = descriptor["selection"]["sites"]
    if len(candidates) != 1000 or len(descriptor["partitions"]) != 62:
        raise ValueError("unexpected pinned historical selection/partition shape")
    with (output / "candidates.vcf").open("x") as stream:
        stream.write("##fileformat=VCFv4.2\n")
        for contig in contigs:
            stream.write(f"##contig=<ID={contig['name']},length={contig['length']}>\n")
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for site in candidates:
            alleles = [chr(base) for base in site["alternates"]]
            ref = chr(site["reference"])
            if ref not in "ACGT" or any(alt not in "ACGT" or alt == ref for alt in alleles):
                raise ValueError("historical selection contains unsupported alleles")
            stream.write(f"{contigs[site['contig']]['name']}\t{site['position'] + 1}\t.\t{ref}\t{','.join(alleles)}\t.\t.\t.\n")
    write_json(output / "members.json", {"version": 1, "members": [
        {"specimen_id": "HG002", "dataset": "leaf/evidence-dataset.manifest.json"}]})
    expected = []
    with (output / "oracle.tsv").open() as stream:
        oracle = {(row["#contig"], int(row["pos"])): row for row in csv.DictReader(stream, delimiter="\t")}
    if len(oracle) != 1000:
        raise ValueError("retained independent oracle cardinality changed")
    for site in candidates:
        contig, pos, ref = contigs[site["contig"]]["name"], site["position"] + 1, chr(site["reference"])
        row = oracle[(contig, pos)]
        if row["ref"] != ref:
            raise ValueError("retained oracle/reference annotation disagrees with descriptor")
        depth = int(row["callable_depth"])
        if sum(int(row[base]) for base in "acgt") != depth:
            raise ValueError("retained oracle callable/allele arithmetic is inconsistent")
        for alt in sorted(map(chr, site["alternates"])):
            allele = int(row[alt.lower()])
            for threshold in THRESHOLDS:
                eligible = depth >= threshold
                ed, ea = (depth, allele) if eligible else (0, 0)
                expected.append({"contig": contig, "pos": pos, "ref": ref, "alt": alt,
                    "min_callable_depth": threshold, "n_requested": 1, "n_observed": 1,
                    "n_depth_eligible": int(eligible), "n_alt_supported": int(eligible and allele > 0),
                    "total_callable_depth": depth, "total_alt_count": allele,
                    "eligible_callable_depth": ed, "eligible_alt_count": ea,
                    "alt_fraction_numerator": allele if depth else None, "alt_fraction_denominator": depth or None,
                    "eligible_alt_fraction_numerator": ea if ed else None,
                    "eligible_alt_fraction_denominator": ed or None,
                    "reference_validation": "checked_on_observed_rows"})
    historical = json.loads((output / "historical-extraction-receipt.json").read_text())
    inventory["historical_extraction"] = {key: historical["params"].get(key) for key in
                                        ("code_git_sha", "code_dirty", "producer.version")}
    inventory["historical_oracle"] = json.loads((output / "oracle-metadata.json").read_text())
    inventory["scientific_scope"] = {"named_sample": "HG002", "candidate_loci": len(candidates),
        "physical_partitions": len(descriptor["partitions"]), "physical_fields_mask": descriptor["fields"],
        "projection": ["depths", "alleles"], "thresholds": list(THRESHOLDS), "expected_summary_rows": len(expected)}
    write_json(output / "extracted-inputs.json", inventory)
    write_json(output / "expected-summary.json", expected)
    return inventory, expected


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    inventory, expected = prepare(output)
    if args.prepare_only:
        print(f"Prepared {inventory['extracted_bytes']} historical bytes; native checks not executed")
        return
    if args.binary is None:
        raise ValueError("--binary required for native saved-data validation")
    binary = args.binary.resolve(strict=True)
    report = {"version": 1, "status": "running", "input_archive_sha256": ARCHIVE_SHA256,
        "historical_extraction": inventory["historical_extraction"],
        "current_binary_sha256": digest(binary),
        "current_binary_version": subprocess.check_output([str(binary), "--version"], text=True).strip(),
        "current_source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "scope": "Current saved-only reanalysis of one historical HG002 leaf against retained independent CIGAR oracle outputs; oracle and raw extraction were not rerun.",
        "not_established": ["fresh-versus-saved amortization", "representative throughput", "multiple biological specimens"],
        "commands": [], "configurations": []}
    try:
        command = [str(binary), "dataset", "verify", "--dataset", str(output / "leaf/evidence-dataset.manifest.json"),
                   "--memory-budget-mb", "256"]
        with (output / "native-verify.json").open("x") as stdout, (output / "native-verify.stderr.log").open("x") as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr)
        report["commands"].append({"argv": command, "exit_code": process.returncode})
        if process.returncode:
            raise RuntimeError("current native reader refused historical leaf verification")
        before = None
        for window in (128, 512):
            destination = output / f"window-{window}"
            command = [sys.executable, str(HERE / "analyze.py"), "--binary", str(binary),
                "--members", str(output / "members.json"), "--sites", str(output / "candidates.vcf"),
                "--output", str(destination), "--missingness", "strict", "--window-size", str(window),
                "--native-memory-budget-mb", "256", "--parent-memory-budget-mb", "128"]
            for threshold in THRESHOLDS:
                command += ["--min-callable-depth", str(threshold)]
            started = time.perf_counter()
            with (output / f"window-{window}.stdout.log").open("x") as stdout, (output / f"window-{window}.stderr.log").open("x") as stderr:
                process = subprocess.run(command, stdout=stdout, stderr=stderr)
            report["commands"].append({"argv": command, "exit_code": process.returncode})
            if process.returncode:
                raise RuntimeError(f"archived reanalysis refused at window {window}")
            actual = [json.loads(line) for line in (destination / "summary.jsonl").read_text().splitlines()]
            if actual != expected:
                raise AssertionError("current summary differs from retained independent oracle")
            summary_hash = digest(destination / "summary.jsonl")
            if before is not None and summary_hash != before:
                raise AssertionError("canonical science changes with window size")
            before = summary_hash
            report["configurations"].append({"window_size": window, "rows": len(actual),
                "summary_sha256": summary_hash, "exact_historical_oracle_agreement": True,
                "outer_wall_seconds": time.perf_counter() - started,
                "timing_scope": "single correctness run including oracle comparison; not a performance benchmark"})
        if digest(ARCHIVE) != ARCHIVE_SHA256:
            raise ValueError("source archive changed during validation")
        report["status"] = "passed"
    except BaseException as error:
        report["status"] = "failed"
        report["error"] = repr(error)
        raise
    finally:
        write_json(output / "results.json", report)
    print(output / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true", help="bounded archive/oracle checks only; no native execution")
    run(parser.parse_args())
