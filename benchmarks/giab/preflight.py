#!/usr/bin/env python3
"""Verify every manifested prepared GIAB input, including the evaluator's FAI."""
import argparse
import hashlib
import json
from pathlib import Path

REQUIRED_ARTIFACTS = {
    "GRCh38.chr20.fa", "GRCh38.chr20.fa.fai",
    "HG002.chr20.bam", "HG002.chr20.bam.bai",
    "HG002.v5.0q.chr20.vcf", "HG002.v5.0q.chr20.bed", "stratifications.tsv",
    "GRCh38_AllTandemRepeatsandHomopolymers_slop5.chr20.bed.gz",
    "GRCh38_lowmappabilityall.chr20.bed.gz",
    "GRCh38_segdups.chr20.bed.gz", "GRCh38_alldifficultregions.chr20.bed.gz",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_reference_index(reference, index):
    """Derive the single-contig faidx layout with bounded lines, without writes."""
    with reference.open("rb") as stream:
        header = stream.readline(1 << 20)
        if len(header) == 1 << 20 or not header.startswith(b">") or header[1:].split()[:1] != [b"chr20"]:
            raise ValueError("prepared reference must contain exactly chr20")
        offset = stream.tell()
        total, line_bases, line_bytes, previous = 0, None, None, None
        while True:
            line = stream.readline(1 << 20)
            if not line:
                break
            if len(line) == 1 << 20 or line.startswith(b">"):
                raise ValueError("prepared reference has an oversized line or extra contig")
            bases = len(line.rstrip(b"\r\n"))
            if bases == 0:
                raise ValueError("prepared reference has an empty sequence line")
            if previous is not None and previous != (line_bases, line_bytes):
                raise ValueError("prepared reference has inconsistent FASTA wrapping")
            if line_bases is None:
                line_bases, line_bytes = bases, len(line)
            if bases > line_bases:
                raise ValueError("prepared reference has inconsistent FASTA wrapping")
            previous = bases, len(line)
            total += bases
    if total == 0:
        raise ValueError("prepared reference is empty")
    expected = "chr20\t{}\t{}\t{}\t{}".format(total, offset, line_bases, line_bytes)
    with index.open() as stream:
        first = stream.readline(4096)
        if first.rstrip("\r\n") != expected or stream.read(1):
            raise ValueError("prepared reference FAI does not match the FASTA layout")


def verify_prepared(data):
    root = Path(data)
    manifest = json.loads((root / "data-manifest.json").read_text())
    if not isinstance(manifest, dict) or manifest.get("schema") != 1:
        raise ValueError("unsupported GIAB data-manifest schema")
    artifacts = manifest.get("prepared_artifacts", {})
    if not isinstance(artifacts, dict):
        raise ValueError("prepared_artifacts must be an object")
    missing = REQUIRED_ARTIFACTS - artifacts.keys()
    if missing:
        raise ValueError("prepared manifest is missing required artifacts: " + ", ".join(sorted(missing)))
    for name, identity in sorted(artifacts.items()):
        if not isinstance(identity, dict):
            raise ValueError("prepared artifact identity must be an object: " + name)
        if Path(name).name != name or name in {".", ".."}:
            raise ValueError("invalid prepared artifact name: " + name)
        path = root / "prepared" / name
        if not path.is_file():
            raise ValueError("missing prepared artifact: " + name)
        if path.stat().st_size != identity.get("bytes") or sha256(path) != identity.get("sha256"):
            raise ValueError("prepared artifact differs from data manifest: " + name)
    verify_reference_index(root / "prepared/GRCh38.chr20.fa", root / "prepared/GRCh38.chr20.fa.fai")
    return len(artifacts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    args = parser.parse_args()
    try:
        count = verify_prepared(args.data)
    except (OSError, ValueError, TypeError) as error:
        parser.exit(2, "GIAB input preflight failed: {}\n".format(error))
    print("Verified {} prepared GIAB artifacts and reference FAI".format(count))


if __name__ == "__main__":
    main()
