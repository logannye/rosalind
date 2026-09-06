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
    query = subprocess.Popen(
        ["bcftools", "query", "--format", "%CHROM\\t%POS\\t%REF[\\t%DP\\t%AD]\\n"],
        stdin=mpileup.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
    if mpileup.stdout:
        mpileup.stdout.close()
    try:
        assert query.stdout is not None
        with open(args.output, "w", newline="", encoding="utf-8") as stream:
            write_rows(query.stdout, stream)
        query.stdout.close()
        query_code, mpileup_code = query.wait(), mpileup.wait()
        if query_code or mpileup_code:
            raise SystemExit(
                f"bcftools query exited {query_code}; mpileup exited {mpileup_code}"
            )
    finally:
        # A broken sink or malformed row must not leave either child alive.
        for process in (query, mpileup):
            if process.poll() is None:
                process.terminate()
            process.wait()


def write_rows(lines, stream) -> None:
    """Normalize one row at a time; never retain the chromosome's query output."""
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)
    for line in lines:
        contig, pos, ref, depth, _ad = line.rstrip("\n").split("\t")
        # AD is REF,ALT-order rather than A/C/G/T. Preserve only DP here and
        # leave non-equivalent fields explicit instead of fabricating equality.
        writer.writerow([contig, pos, ref, depth, depth, *("NA" for _ in range(14))])


if __name__ == "__main__":
    main()
