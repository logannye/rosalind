#!/usr/bin/env python3
"""Generate a deterministic Illumina-style toy dataset for Rosalind demos."""

from __future__ import annotations

import argparse
import hashlib
import random
from pathlib import Path
from typing import Tuple

NUCLEOTIDES = "ACGT"
DEFAULT_REF_LENGTH = 1_000_000
READ_LENGTH = 150
COVERAGE = 10  # approximate paired depth
ERROR_RATE = 0.01
SEED = 1337


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices(NUCLEOTIDES, k=length))


def introduce_errors(seq: str, rng: random.Random, error_rate: float) -> str:
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < error_rate:
            alternatives = [n for n in NUCLEOTIDES if n != base]
            bases[i] = rng.choice(alternatives)
    return "".join(bases)


def revcomp(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]


def format_fastq(name: str, sequence: str, quality: str) -> str:
    return f"@{name}\n{sequence}\n+\n{quality}\n"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def generate_dataset(output_dir: Path, ref_length: int, coverage: int, seed: int) -> Tuple[Path, Path, Path]:
    rng = random.Random(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference = random_dna(ref_length, rng)
    reference_path = output_dir / "reference.fa"
    with reference_path.open("w", encoding="ascii") as handle:
        handle.write(">chrToy\n")
        for i in range(0, len(reference), 80):
            handle.write(reference[i : i + 80] + "\n")

    total_reads = max(1, int((ref_length * coverage) / READ_LENGTH))
    pair_count = max(1, total_reads // 2)

    r1_path = output_dir / "reads_R1.fastq"
    r2_path = output_dir / "reads_R2.fastq"

    max_start = len(reference) - READ_LENGTH - 1
    quality = "I" * READ_LENGTH

    with r1_path.open("w", encoding="ascii") as r1_handle, r2_path.open("w", encoding="ascii") as r2_handle:
        for idx in range(pair_count):
            start = rng.randint(0, max_start)
            fragment = reference[start : start + READ_LENGTH]
            mate = revcomp(reference[start : start + READ_LENGTH])

            read1_seq = introduce_errors(fragment, rng, ERROR_RATE)
            read2_seq = introduce_errors(mate, rng, ERROR_RATE)

            r1_handle.write(format_fastq(f"toy_read_{idx}/1", read1_seq, quality))
            r2_handle.write(format_fastq(f"toy_read_{idx}/2", read2_seq, quality))

    return reference_path, r1_path, r2_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic Illumina-style dataset")
    parser.add_argument("output", type=Path, help="Directory to place generated files")
    parser.add_argument("--length", type=int, default=DEFAULT_REF_LENGTH, help="Reference length (bp)")
    parser.add_argument("--coverage", type=int, default=COVERAGE, help="Approximate paired depth")
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed")
    args = parser.parse_args()

    reference_path, r1_path, r2_path = generate_dataset(args.output, args.length, args.coverage, args.seed)

    checksums = {
        reference_path.name: sha256_file(reference_path),
        r1_path.name: sha256_file(r1_path),
        r2_path.name: sha256_file(r2_path),
    }

    checksum_path = args.output / "SHA256SUMS"
    with checksum_path.open("w", encoding="ascii") as handle:
        for name, digest in checksums.items():
            handle.write(f"{digest}  {name}\n")

    print("Generated dataset:")
    for name, digest in checksums.items():
        print(f"  {name}: {digest}")


if __name__ == "__main__":
    main()
