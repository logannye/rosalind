#!/usr/bin/env python3
"""Normalize the explicitly comparable subset of bcftools mpileup fields."""

import argparse
import csv
import subprocess

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
    mpileup = subprocess.Popen(
        ["bcftools", "mpileup", "--fasta-ref", args.reference, "--max-depth", str(args.max_depth),
         "--annotate", "FORMAT/AD,FORMAT/DP", "--output-type", "u", args.bam],
        stdout=subprocess.PIPE,
    )
    query = subprocess.run(
        ["bcftools", "query", "--format", "%CHROM\\t%POS\\t%REF[\\t%DP\\t%AD]\\n"],
        stdin=mpileup.stdout,
        capture_output=True,
        text=True,
        check=False,
    )
    if mpileup.stdout:
        mpileup.stdout.close()
    mpileup_code = mpileup.wait()
    if query.returncode or mpileup_code:
        raise SystemExit(query.stderr or f"bcftools mpileup exited {mpileup_code}")
    with open(args.output, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for line in query.stdout.splitlines():
            contig, pos, ref, depth, ad = line.split("\t")
            allele_depths = [int(value) for value in ad.split(",") if value not in {".", ""}]
            # AD is REF,ALT-order rather than A/C/G/T. Preserve only DP here and
            # leave non-equivalent fields explicit instead of fabricating equality.
            writer.writerow([contig, pos, ref, depth, depth, *("NA" for _ in range(14))])


if __name__ == "__main__":
    main()
