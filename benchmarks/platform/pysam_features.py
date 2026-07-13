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
                stepper="nofilter",
                min_base_quality=0,
                max_depth=args.max_depth,
                ignore_overlaps=False,
                ignore_orphans=False,
            ):
                observations = []
                for item in column.pileups:
                    read = item.alignment
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                        continue
                    if item.is_del or item.is_refskip or item.query_position is None:
                        continue
                    offset = item.query_position
                    base = read.query_sequence[offset].upper()
                    if base not in "ACGT":
                        continue
                    qualities = read.query_qualities
                    observations.append((base, read.is_reverse, qualities[offset] if qualities else 0, read.mapping_quality))
                if not observations:
                    continue
                counts = {(base, reverse): 0 for base in "ACGT" for reverse in (False, True)}
                for base, reverse, _, _ in observations:
                    counts[(base, reverse)] += 1
                depth = len(observations)
                alleles = [counts[(base, False)] + counts[(base, True)] for base in "ACGT"]
                strands = [counts[(base, reverse)] for base in "ACGT" for reverse in (False, True)]
                ref = reference.fetch(contig, column.reference_pos, column.reference_pos + 1).upper() or "N"
                writer.writerow([
                    contig, column.reference_pos + 1, ref, depth, depth, *alleles, *strands,
                    f"{sum(item[2] for item in observations) / depth:.2f}",
                    f"{sum(item[3] for item in observations) / depth:.2f}",
                ])


if __name__ == "__main__":
    main()
