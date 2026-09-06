#!/usr/bin/env python3
"""Maintained transparent pysam counterpart to Rosalind's feature extractor."""

import argparse
import csv
import pysam

HEADER = [
    "#contig", "pos", "ref", "depth", "raw_depth", "a", "c", "g", "t",
    "a_fwd", "a_rev", "c_fwd", "c_rev", "g_fwd", "g_rev", "t_fwd", "t_rev",
    "mean_bq", "mean_mapq",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-depth", type=int, default=1000)
    args = parser.parse_args()

    with pysam.FastaFile(args.reference) as reference, pysam.AlignmentFile(args.bam, "rb") as bam, open(
        args.output, "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for contig in reference.references:
            for column in bam.pileup(
                contig,
                0,
                reference.get_reference_length(contig),
                truncate=True,
                stepper="all",
                flag_filter=4 | 256 | 1024 | 2048,
                min_base_quality=0,
                # Observe one overflow read and fail rather than silently sample.
                max_depth=args.max_depth + 1,
                ignore_overlaps=False,
                ignore_orphans=False,
                compute_baq=False,
            ):
                counts = {(base, reverse): 0 for base in "ACGT" for reverse in (False, True)}
                raw_depth, depth, sum_bq, sum_mapq = 0, 0, 0, 0
                for item in column.pileups:
                    read = item.alignment
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                        continue
                    raw_depth += 1
                    if raw_depth > args.max_depth:
                        raise SystemExit("pysam active depth exceeds the comparison capacity; no downsampling performed")
                    if item.is_del or item.is_refskip or item.query_position is None:
                        continue
                    offset = item.query_position
                    base = read.query_sequence[offset].upper()
                    if base not in "ACGT":
                        continue
                    qualities = read.query_qualities
                    counts[(base, read.is_reverse)] += 1
                    depth += 1
                    sum_bq += qualities[offset] if qualities else 0
                    sum_mapq += read.mapping_quality
                if not depth:
                    continue
                alleles = [counts[(base, False)] + counts[(base, True)] for base in "ACGT"]
                strands = [counts[(base, reverse)] for base in "ACGT" for reverse in (False, True)]
                ref = reference.fetch(contig, column.reference_pos, column.reference_pos + 1).upper() or "N"
                writer.writerow([
                    contig, column.reference_pos + 1, ref, depth, raw_depth, *alleles, *strands,
                    f"{((sum_bq * 100 + depth // 2) // depth) / 100:.2f}",
                    f"{((sum_mapq * 100 + depth // 2) // depth) / 100:.2f}",
                ])


if __name__ == "__main__":
    main()
