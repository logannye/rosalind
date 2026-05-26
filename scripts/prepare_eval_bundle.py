#!/usr/bin/env python3
"""
Prepare a deterministic evaluation bundle for Rosalind clinical benchmarking.

This script:
- Extracts a single contig from a reference FASTA into `reference.fa`
- Deterministically downsamples paired-end FASTQs (tumor + normal) into the bundle
- Writes SHA256SUMS for reproducibility

Notes:
- This does NOT require network access and does not attempt to fetch public datasets.
- FASTQs are expected to be uncompressed. (Decompress externally if needed.)
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Tuple


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_u64(seed: int, s: str) -> int:
    h = hashlib.blake2b(digest_size=8, person=b"rosalind", key=str(seed).encode("ascii"))
    h.update(s.encode("utf-8", errors="strict"))
    return int.from_bytes(h.digest(), byteorder="little", signed=False)


def normalize_read_name(raw: str) -> str:
    # Strip FASTQ '@' and any Illumina whitespace suffix; drop /1 or /2.
    name = raw.strip()
    if name.startswith("@"):
        name = name[1:]
    name = name.split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


@dataclass(frozen=True)
class FastqRecord:
    name: str
    seq: str
    plus: str
    qual: str


def fastq_iter(path: Path) -> Iterator[FastqRecord]:
    with path.open("r", encoding="ascii") as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not (seq and plus and qual):
                raise ValueError(f"truncated FASTQ at {path}")
            yield FastqRecord(name=name.rstrip("\n"), seq=seq, plus=plus, qual=qual)


def write_fastq_record(out, rec: FastqRecord) -> None:
    out.write(rec.name)
    out.write("\n")
    out.write(rec.seq)
    out.write(rec.plus)
    out.write(rec.qual)


def select_pairs(
    r1: Path,
    r2: Path,
    seed: int,
    max_pairs: Optional[int],
    sample_prob: Optional[float],
) -> Iterator[Tuple[FastqRecord, FastqRecord]]:
    it1 = fastq_iter(r1)
    it2 = fastq_iter(r2)
    kept = 0
    for rec1, rec2 in zip(it1, it2):
        name1 = normalize_read_name(rec1.name)
        name2 = normalize_read_name(rec2.name)
        if name1 != name2:
            raise ValueError(f"mate name mismatch: {name1} vs {name2}")

        keep = True
        if sample_prob is not None:
            # Deterministic sampling by stable hash.
            keep = (stable_u64(seed, name1) / float(2**64)) < sample_prob

        if keep:
            yield rec1, rec2
            kept += 1
            if max_pairs is not None and kept >= max_pairs:
                return


def extract_contig(reference_fasta: Path, contig: str, out_fasta: Path) -> int:
    """
    Extract a single contig sequence from a FASTA.
    Returns length of extracted sequence.
    """
    found = False
    seq_parts = []
    header = None

    with reference_fasta.open("r", encoding="ascii") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if found:
                    break
                if name == contig:
                    found = True
                    header = line
                    continue
            if found and not line.startswith(">"):
                seq_parts.append(line.strip())

    if not found or header is None:
        raise ValueError(f"contig '{contig}' not found in {reference_fasta}")

    seq = "".join(seq_parts).upper()
    if not seq:
        raise ValueError(f"contig '{contig}' has empty sequence in {reference_fasta}")

    with out_fasta.open("w", encoding="ascii") as out:
        out.write(header.split()[0] + "\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i : i + 80] + "\n")

    return len(seq)


def main() -> None:
    ap = argparse.ArgumentParser(description="Prepare deterministic clinical eval bundle")
    ap.add_argument("--reference", type=Path, required=True, help="Input reference FASTA")
    ap.add_argument("--contig", type=str, required=True, help="Contig name to extract (e.g. chr20)")
    ap.add_argument("--tumor-r1", type=Path, required=True, help="Tumor FASTQ R1 (uncompressed)")
    ap.add_argument("--tumor-r2", type=Path, required=True, help="Tumor FASTQ R2 (uncompressed)")
    ap.add_argument("--normal-r1", type=Path, required=True, help="Normal FASTQ R1 (uncompressed)")
    ap.add_argument("--normal-r2", type=Path, required=True, help="Normal FASTQ R2 (uncompressed)")
    ap.add_argument("--out", type=Path, required=True, help="Output directory for bundle")
    ap.add_argument("--seed", type=int, default=1337, help="Sampling seed")
    ap.add_argument("--max-pairs", type=int, default=20000, help="Max read pairs per sample")
    ap.add_argument(
        "--sample-prob",
        type=float,
        default=None,
        help="Optional sampling probability (0..1) per read-pair; deterministic hash-based",
    )
    args = ap.parse_args()

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_out = out_dir / "reference.fa"
    tumor_r1_out = out_dir / "tumor_R1.fastq"
    tumor_r2_out = out_dir / "tumor_R2.fastq"
    normal_r1_out = out_dir / "normal_R1.fastq"
    normal_r2_out = out_dir / "normal_R2.fastq"

    ref_len = extract_contig(args.reference, args.contig, ref_out)

    def write_pairs(src_r1: Path, src_r2: Path, dst_r1: Path, dst_r2: Path) -> int:
        kept = 0
        with dst_r1.open("w", encoding="ascii") as o1, dst_r2.open("w", encoding="ascii") as o2:
            for r1_rec, r2_rec in select_pairs(
                src_r1, src_r2, args.seed, args.max_pairs, args.sample_prob
            ):
                write_fastq_record(o1, r1_rec)
                write_fastq_record(o2, r2_rec)
                kept += 1
        return kept

    tumor_pairs = write_pairs(args.tumor_r1, args.tumor_r2, tumor_r1_out, tumor_r2_out)
    normal_pairs = write_pairs(args.normal_r1, args.normal_r2, normal_r1_out, normal_r2_out)

    manifest = {
        "schema_version": 1,
        "contig": args.contig,
        "reference_len": ref_len,
        "seed": args.seed,
        "max_pairs_per_sample": args.max_pairs,
        "sample_prob": args.sample_prob,
        "tumor_pairs": tumor_pairs,
        "normal_pairs": normal_pairs,
        "files": {
            "reference.fa": "reference.fa",
            "tumor_R1.fastq": "tumor_R1.fastq",
            "tumor_R2.fastq": "tumor_R2.fastq",
            "normal_R1.fastq": "normal_R1.fastq",
            "normal_R2.fastq": "normal_R2.fastq",
        },
    }

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    checksums = {
        "reference.fa": sha256_file(ref_out),
        "tumor_R1.fastq": sha256_file(tumor_r1_out),
        "tumor_R2.fastq": sha256_file(tumor_r2_out),
        "normal_R1.fastq": sha256_file(normal_r1_out),
        "normal_R2.fastq": sha256_file(normal_r2_out),
        "manifest.json": sha256_file(out_dir / "manifest.json"),
    }
    with (out_dir / "SHA256SUMS").open("w", encoding="ascii") as handle:
        for name in sorted(checksums.keys()):
            handle.write(f"{checksums[name]}  {name}\n")

    print(f"Bundle written to: {out_dir}")
    print(f"  reference: {args.contig} ({ref_len} bp)")
    print(f"  tumor pairs: {tumor_pairs}")
    print(f"  normal pairs: {normal_pairs}")
    print(f"  checksums: {out_dir / 'SHA256SUMS'}")


if __name__ == "__main__":
    main()


