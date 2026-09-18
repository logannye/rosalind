#!/usr/bin/env python3
"""Prepare content-locked HG002 chr20 workloads; no simulated or replicated reads."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import urllib.request


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def download(row, cache):
    identifier, name, expected, url = row
    destination = cache / name
    started = time.monotonic()
    if not destination.exists():
        partial = cache / (name + f".download-{os.getpid()}")
        try:
            with urllib.request.urlopen(url, timeout=90) as source, partial.open("xb") as output:
                shutil.copyfileobj(source, output, 1 << 20)
            if digest(partial) != expected:
                raise ValueError(f"download SHA256 mismatch: {identifier}")
            # Preparation has one writer per locked source; never replace a
            # concurrently supplied cache object without checking its bytes.
            if destination.exists():
                if digest(destination) != expected:
                    raise ValueError(f"concurrent cache SHA256 mismatch: {identifier}")
                partial.unlink()
            else:
                partial.rename(destination)
        finally:
            partial.unlink(missing_ok=True)
    if digest(destination) != expected:
        raise ValueError(f"cached SHA256 mismatch: {identifier}")
    return {"id": identifier, "path": str(destination.resolve()), "url": url,
            "sha256": expected, "bytes": destination.stat().st_size,
            "verification_download_seconds": time.monotonic() - started}


def select_sites(fasta, contig, start, end, count):
    """One A/C/G/T site per equal interval; deterministic, without inspecting reads."""
    if not 0 <= start < end <= fasta.get_reference_length(contig) or end - start < count:
        raise ValueError("invalid site interval/count")
    result = []
    for index in range(count):
        left = start + (end - start) * index // count
        right = start + (end - start) * (index + 1) // count
        selected = None
        for window in range(left, right, 4096):
            bases = fasta.fetch(contig, window, min(window + 4096, right)).upper()
            for offset, base in enumerate(bases):
                if base in "ACGT":
                    selected = (window + offset, base)
                    break
            if selected:
                break
        if selected:
            result.append(selected)
    if not result:
        raise ValueError("selection contains no unambiguous reference bases")
    return result


def write_sites(path, contig, length, sites):
    with path.open("x") as output:
        output.write(f"##fileformat=VCFv4.2\n##contig=<ID={contig},length={length}>\n")
        output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for position, reference in sites:
            alternate = next(base for base in "ACGT" if base != reference)
            output.write(f"{contig}\t{position + 1}\t.\t{reference}\t{alternate}\t.\t.\t.\n")


def prepare(args):
    import pysam
    if pysam.__version__ != "0.23.3":
        raise ValueError("preparation requires pysam==0.23.3")
    if args.output.exists():
        raise ValueError("output directory must be new; verified downloads may be reused separately")
    rows = [line.split("\t") for line in args.source_lock.read_text().splitlines()
            if line and not line.startswith("#")]
    required = {"reads_bam", "reads_bai", "reference"}
    rows = [row for row in rows if row[0] in required]
    if len(rows) != 3 or {row[0] for row in rows} != required:
        raise ValueError("source lock requires exactly one BAM, BAI and reference")
    args.download_cache.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        sources = list(pool.map(lambda row: download(row, args.download_cache), rows))
    by_id = {source["id"]: source for source in sources}
    args.output.mkdir(parents=True)
    root = args.output.resolve()
    started = time.monotonic()
    reference = root / "reference.fa"
    with gzip.open(by_id["reference"]["path"], "rb") as source, reference.open("xb") as output:
        shutil.copyfileobj(source, output, 1 << 20)
    pysam.faidx(str(reference))
    bam = Path(by_id["reads_bam"]["path"])
    bai = Path(by_id["reads_bai"]["path"])
    window_bam, window_cram = root / "window.bam", root / "window.cram"
    stats = {"window_records": 0, "window_query_bases": 0, "window_max_read_length": 0}
    with pysam.AlignmentFile(str(bam), "rb", index_filename=str(bai)) as source, \
            pysam.FastaFile(str(reference)) as fasta:
        dictionary = dict(zip(source.references, source.lengths))
        if dictionary != dict(zip(fasta.references, fasta.lengths)):
            raise ValueError("locked source BAM/reference dictionaries disagree")
        if args.contig not in dictionary or not 0 <= args.start < args.end <= dictionary[args.contig]:
            raise ValueError("requested window outside reference dictionary")
        header = source.header.to_dict()
        with pysam.AlignmentFile(str(window_bam), "wb", header=header) as out_bam, \
                pysam.AlignmentFile(str(window_cram), "wc", header=header,
                                    reference_filename=str(reference),
                                    format_options=[b"version=3.0"]) as out_cram:
            for record in source.fetch(args.contig, args.start, args.end):
                out_bam.write(record)
                out_cram.write(record)
                stats["window_records"] += 1
                stats["window_query_bases"] += record.query_length
                stats["window_max_read_length"] = max(stats["window_max_read_length"], record.query_length)
        full_sites = select_sites(fasta, args.contig, 0, dictionary[args.contig], args.sites)
        window_sites = select_sites(fasta, args.contig, args.start, args.end, args.sites)
        write_sites(root / "chromosome-sparse.vcf", args.contig, dictionary[args.contig], full_sites)
        write_sites(root / "window-sparse.vcf", args.contig, dictionary[args.contig], window_sites)
        targets = []
        # Twenty separated500-base targets span the entire1Mb window. The first
        # and last stay150bases inside it, avoiding slice-edge read truncation.
        margin, width, count = 1000, 500, 20
        if args.end - args.start < 2 * margin + count * width:
            raise ValueError("window too short for deterministic panel")
        span = args.end - args.start - 2 * margin - width
        for index in range(count):
            left = args.start + margin + span * index // (count - 1)
            targets.append((left, left + width))
        (root / "panel.bed").write_text("".join(f"{args.contig}\t{left}\t{right}\tT{i+1}\n"
                                               for i, (left, right) in enumerate(targets)))
        samples = sorted({entry.get("SM", "") for entry in header.get("RG", [])})
        if samples != ["HG002"]:
            raise ValueError(f"unexpected locked sample identities: {samples}")
    pysam.index(str(window_bam))
    pysam.index(str(window_cram))
    identities = {}
    def operand(path):
        path = Path(path)
        key = str(path.resolve())
        if key not in identities:
            identities[key] = {"path": os.path.relpath(path, root), "sha256": digest(path)}
        return identities[key].copy()
    cases = [
        {"id": "m96-t512", "budget_mib": 96, "tile_bases": 512, "workers": 1, "cache": "none"},
        {"id": "m128-t4096", "budget_mib": 128, "tile_bases": 4096, "workers": 1, "cache": "none"},
        {"id": "m256-t16384", "budget_mib": 256, "tile_bases": 16384, "workers": 1, "cache": "none"},
        {"id": "workers2", "budget_mib": 512, "tile_bases": 4096, "workers": 2, "cache": "none"},
        {"id": "workers8-cache", "budget_mib": 1024, "tile_bases": 16384, "workers": 8, "cache": "cold-resume"},
    ]
    workloads = []
    def workload(identifier, group, alignment, index, selection, kind, rows, parallel):
        workloads.append({"id": identifier, "equivalence_group": group,
                          "alignments": operand(alignment), "alignment_index": operand(index),
                          "reference": operand(reference), "reference_fai": operand(str(reference) + ".fai"),
                          "selection": {"kind": kind, **operand(root / selection)}, "sample": "HG002",
                          "expected_rows": rows,
                          "cases": [case["id"] for case in (cases if parallel else cases[:3])],
                          "gates": {"min_admitted_budgets": 3, "min_effective_tiles": 2,
                                    "workers": [1, 2, 8] if parallel else [1]}})
    workload("chromosome-bam-sparse", "chromosome-sparse", bam, bai,
             "chromosome-sparse.vcf", "sites", len(full_sites), False)
    for kind, selection, rows in [("sites", "window-sparse.vcf", len(window_sites)),
                                  ("regions", "panel.bed", sum(right-left for left,right in targets))]:
        for encoding, alignment, index in [("bam", window_bam, str(window_bam) + ".bai"),
                                           ("cram", window_cram, str(window_cram) + ".crai")]:
            workload(f"window-{encoding}-{kind}", f"window-{kind}", alignment, index,
                     selection, kind, rows, True)
    provenance = {"source_lock_sha256": digest(args.source_lock), "sources": sources,
                  "preparation_script_sha256": digest(Path(__file__)), "argv": sys.argv,
                  "pysam": pysam.__version__, "samtools": pysam.__samtools_version__,
                  "window": {"contig": args.contig, "start": args.start, "end": args.end},
                  "dictionary_contigs": len(dictionary), "samples": samples, **stats,
                  "preparation_seconds": time.monotonic() - started,
                  "limits": ["Public HG002 short-read data; no reads replicated or simulated.",
                             "Panel selection models an analysis request, not a target-enrichment assay.",
                             "VCF ALT bases are deterministic evidence probes, not validated variants.",
                             "Full source dictionary and all overlapping window reads are retained.",
                             "CRAM3.0 is locally encoded from the same unfiltered source records.",
                             "Three serial budgets and extra parallel/cache cases are not a full Cartesian matrix."]}
    manifest = {"schema": 1, "label": "hg002-chr20-v1", "provenance": provenance,
                "repeats": 3, "seed": 17, "oracle": {"tile_bases": 1024, "max_record_bytes": 1048576},
                "cases": cases, "workloads": workloads}
    (root / "preparation.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (root / "workloads.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(root / "workloads.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-lock", type=Path,
                        default=Path(__file__).resolve().parents[1] / "giab" / "resources.tsv")
    parser.add_argument("--download-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contig", default="chr20")
    parser.add_argument("--start", type=int, default=10000000)
    parser.add_argument("--end", type=int, default=11000000)
    parser.add_argument("--sites", type=int, default=1000)
    args = parser.parse_args()
    if args.sites < 1:
        parser.error("sites must be positive")
    prepare(args)


if __name__ == "__main__":
    main()
