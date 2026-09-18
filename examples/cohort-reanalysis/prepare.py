#!/usr/bin/env python3
"""Prepare authored synthetic single-base reads with pinned pysam; no downloads.

This produces inputs for existing single-sample Rosalind commands, not a cohort API.
Metadata and read order are deterministic; output must be a new directory.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if pysam.__version__ != "0.23.3":
        parser.error("fixture preparation requires pysam==0.23.3")
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()
    with (HERE / "members.tsv").open() as source:
        members = list(csv.DictReader(source, delimiter="\t"))
    with (HERE / "reads.tsv").open() as source:
        records = list(csv.DictReader(source, delimiter="\t"))
    for name in ("reference.fa", "candidates.vcf", *(m["selection"] for m in members)):
        shutil.copyfile(HERE / name, raw / name)
    pysam.faidx(str(raw / "reference.fa"))
    for member in members:
        specimen = member["specimen"]
        group = "rg-" + specimen.removeprefix("specimen-")
        sam = raw / (specimen + ".sam")
        with sam.open("x") as stream:
            stream.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:synthetic\tLN:64\n")
            stream.write(f"@RG\tID:{group}\tSM:{member['sample']}\n")
            selected = [row for row in records if row["specimen"] == specimen]
            selected.sort(key=lambda row: int(row["pos"]))
            ordinal = 0
            for row in selected:
                for _ in range(int(row["copies"])):
                    ordinal += 1
                    name = f"{specimen}-{ordinal:04d}"
                    quality = chr(33 + int(row["base_quality"]))
                    stream.write(f"{name}\t{row['flag']}\tsynthetic\t{row['pos']}\t"
                                 f"{row['mapq']}\t1M\t*\t0\t0\t{row['base']}\t"
                                 f"{quality}\tRG:Z:{group}\n")
        bam = raw / (specimen + ".bam")
        with pysam.AlignmentFile(str(sam), "r") as source:
            with pysam.AlignmentFile(str(bam), "wb", template=source) as output:
                for record in source:
                    output.write(record)
        pysam.index(str(bam))
    metadata = {
        "fixture": "cohort-contract-synthetic-v1",
        "pysam": pysam.__version__,
        "samtools": pysam.__samtools_version__,
        "source_sha256": {p.name: sha256(p) for p in sorted(HERE.iterdir())
                          if p.is_file() and p.suffix in (".tsv", ".bed", ".vcf", ".fa", ".json")},
        "prepared_sha256": {p.name: sha256(p) for p in sorted(raw.iterdir()) if p.is_file()},
    }
    (args.output / "preparation.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"Prepared three synthetic named samples in {args.output}")


if __name__ == "__main__":
    main()
