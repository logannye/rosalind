#!/usr/bin/env python3
"""Independent, streaming per-locus oracle for the v1 read-count profile.

Only the small supplied candidate VCF is retained; aligned reads are fetched per
locus. Direct CIGAR traversal avoids htslib pileup's depth sampling/BAQ defaults.
This baseline deliberately trades repeated reads for bounded Python state.
"""
import argparse
import csv
import sys

FILTERS = ["secondary", "supplementary", "qc_fail", "duplicate", "unavailable_mapq",
           "low_mapq", "unavailable_base_quality", "low_base_quality", "ambiguous_base"]
SCALARS = ["prefilter_depth", "aligned_depth", "callable_depth", "a", "c", "g", "t",
           "a_fwd", "a_rev", "c_fwd", "c_rev", "g_fwd", "g_rev", "t_fwd", "t_rev",
           "base_quality_sum", "mapping_quality_sum", "read_position_sum", "read_length_sum"]
HEADER = ["#contig", "pos", "ref", "requested_alts", *SCALARS,
          *["filtered_" + name for name in FILTERS], "base_quality_histogram", "mapping_quality_histogram"]


def query_position(read, position):
    reference, query = read.reference_start, 0
    for operation, length in read.cigartuples or ():
        if operation in (0, 7, 8):
            if reference <= position < reference + length:
                return query + position - reference
            reference += length
            query += length
        elif operation in (2, 3):
            if reference <= position < reference + length:
                return None
            reference += length
        elif operation in (1, 4):
            query += length
    return None


def first_failure(read, query, mapq, base_quality):
    for flag, name in [(0x100, "secondary"), (0x800, "supplementary"),
                       (0x200, "qc_fail"), (0x400, "duplicate")]:
        if read.flag & flag:
            return name
    if read.mapping_quality == 255:
        return "unavailable_mapq"
    if read.mapping_quality < mapq:
        return "low_mapq"
    return None


def row_for(bam, fasta, contig, position, ref, alts, mapq=20, base_quality=20):
    if fasta.fetch(contig, position, position + 1).upper() != ref:
        raise ValueError("VCF REF mismatch")
    values = dict.fromkeys(SCALARS + ["filtered_" + name for name in FILTERS], 0)
    bq_hist, mq_hist = [0] * 94, [0] * 255
    for read in bam.fetch(contig, position, position + 1):
        if read.is_unmapped:
            continue
        query = query_position(read, position)
        if query is None:
            continue
        values["prefilter_depth"] += 1
        failure = first_failure(read, query, mapq, base_quality)
        if failure is None:
            values["aligned_depth"] += 1
            qualities = read.query_qualities
            quality = 255 if qualities is None else qualities[query]
            base = read.query_sequence[query].upper()
            if base == "=":
                base = ref
            if quality == 255:
                failure = "unavailable_base_quality"
            elif quality < base_quality:
                failure = "low_base_quality"
            elif base not in "ACGT":
                failure = "ambiguous_base"
        if failure:
            values["filtered_" + failure] += 1
            continue
        if quality > 93:
            raise ValueError("quality exceeds profile envelope")
        values["callable_depth"] += 1
        values[base.lower()] += 1
        values[base.lower() + ("_rev" if read.is_reverse else "_fwd")] += 1
        values["base_quality_sum"] += quality
        values["mapping_quality_sum"] += read.mapping_quality
        values["read_position_sum"] += read.query_length - 1 - query if read.is_reverse else query
        values["read_length_sum"] += read.query_length
        bq_hist[quality] += 1
        mq_hist[read.mapping_quality] += 1
    histogram = lambda bins: ",".join(f"{quality}:{count}" for quality, count in enumerate(bins) if count) or "."
    return [contig, position + 1, ref, "".join(sorted(alts)),
            *[values[name] for name in HEADER[4:-2]], histogram(bq_hist), histogram(mq_hist)]


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("alignments")
    parser.add_argument("sites")
    args = parser.parse_args()
    with pysam.FastaFile(args.reference) as fasta, pysam.AlignmentFile(args.alignments, "rb") as bam:
        ranks = {name: rank for rank, name in enumerate(fasta.references)}
        sites = {}
        with open(args.sites) as source:
            for line in source:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                key = fields[0], int(fields[1]) - 1
                ref, alts = fields[3].upper(), fields[4].upper().split(",")
                if len(ref) != 1 or ref not in "ACGT" or any(len(a) != 1 or a not in "ACGT" for a in alts):
                    raise ValueError("oracle expects SNVs")
                previous = sites.setdefault(key, (ref, set()))
                if previous[0] != ref:
                    raise ValueError("conflicting REF")
                previous[1].update(alts)
        output = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        output.writerow(HEADER)
        for (contig, position), (ref, alts) in sorted(sites.items(), key=lambda item: (ranks[item[0][0]], item[0][1])):
            output.writerow(row_for(bam, fasta, contig, position, ref, alts))


if __name__ == "__main__":
    main()
