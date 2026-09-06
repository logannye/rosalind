#!/usr/bin/env python3
"""Join a small supplied candidate VCF to evidence; emit an illustrative screen.

This intentionally materializes the small candidate set, not the alignment data
or genome evidence. It is a research triage example, not a calibrated classifier.
"""
import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidates")
    parser.add_argument("evidence")
    args = parser.parse_args()
    candidates = {}
    with open(args.candidates) as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            key = fields[0], fields[1]
            if key in candidates:
                raise ValueError("tutorial expects one candidate VCF record per locus")
            candidates[key] = fields
    output = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    output.writerow(["contig", "pos", "ref", "alt", "callable_depth", "alt_reads",
                     "alt_forward", "alt_reverse", "research_screen"])
    with open(args.evidence) as source:
        for row in csv.DictReader(source, delimiter="\t"):
            key = row["#contig"], row["pos"]
            candidate = candidates.pop(key, None)
            if candidate is None:
                raise ValueError(f"unexpected or duplicate evidence locus: {key}")
            if candidate[3] != row["ref"]:
                raise ValueError(f"candidate REF differs from evidence: {key}")
            depth = int(row["callable_depth"])
            for alt in candidate[4].split(","):
                count = int(row[alt.lower()])
                forward, reverse = int(row[alt.lower() + "_fwd"]), int(row[alt.lower() + "_rev"])
                keep = depth >= 10 and count >= 3 and forward >= 1 and reverse >= 1
                output.writerow([*key, row["ref"], alt, depth, count, forward, reverse,
                                 "review" if keep else "insufficient_example_support"])
    if candidates:
        raise ValueError(f"evidence is missing {len(candidates)} requested loci")


if __name__ == "__main__":
    main()
