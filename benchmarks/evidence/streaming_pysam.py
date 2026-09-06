#!/usr/bin/env python3
"""Independent, bounded direct-CIGAR oracle for Rosalind's legacy 63-field mask.

One selection window and fixed-width integer arrays are retained. Each indexed
window is fetched once; reads, fragments and whole-output rows are never retained.
Native htslib decoding is observed, not protected by a hard allocation limit.
"""
import argparse
from array import array
from bisect import bisect_left
import csv
import json
from pathlib import Path
import sys
import time

FILTERS = ["secondary", "supplementary", "qc_fail", "duplicate", "unavailable_mapq",
           "low_mapq", "unavailable_base_quality", "low_base_quality", "ambiguous_base"]
SCALARS = ["prefilter_depth", "aligned_depth", "callable_depth", "a", "c", "g", "t",
           "a_fwd", "a_rev", "c_fwd", "c_rev", "g_fwd", "g_rev", "t_fwd", "t_rev",
           "base_quality_sum", "mapping_quality_sum", "read_position_sum", "read_length_sum",
           *["filtered_" + name for name in FILTERS]]
HEADER = ["#contig", "pos", "ref", "requested_alts", *SCALARS,
          "base_quality_histogram", "mapping_quality_histogram"]


def selection_loci(kind, path, ranks, lengths):
    """Stream canonical coordinates; reject unsorted input instead of sorting it in RAM."""
    if kind == "regions":
        pending = None
        previous = None
        with open(path) as source:
            for number, line in enumerate(source, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split()
                if len(fields) < 3 or fields[0] not in ranks:
                    raise ValueError(f"invalid BED record {number}")
                name, start, end = fields[0], int(fields[1]), int(fields[2])
                key = ranks[name], start
                if start < 0 or start >= end or end > lengths[name] or (previous and key < previous):
                    raise ValueError("oracle BED must be coordinate sorted and nonempty")
                previous = key
                if pending and pending[0] == name and start <= pending[2]:
                    pending = (name, pending[1], max(end, pending[2]))
                else:
                    if pending:
                        for position in range(pending[1], pending[2]):
                            yield pending[0], position, None, ()
                    pending = (name, start, end)
            if pending:
                for position in range(pending[1], pending[2]):
                    yield pending[0], position, None, ()
        return
    import pysam
    pending = None
    previous = None
    with pysam.VariantFile(str(path)) as source:
        for record in source:
            name, position = record.contig, record.start
            if name not in ranks or position < 0 or position >= lengths[name]:
                raise ValueError("variant coordinate is outside the reference")
            ref, alts = record.ref.upper(), tuple(a.upper() for a in (record.alts or ()))
            if ref not in "ACGT" or len(ref) != 1 or not alts or any(
                    len(a) != 1 or a not in "ACGT" or a == ref for a in alts):
                raise ValueError("oracle selection requires A/C/G/T SNVs")
            key = ranks[name], position
            if previous and key < previous:
                raise ValueError("oracle VCF/BCF must be coordinate sorted")
            if pending and key == previous:
                if ref != pending[2]:
                    raise ValueError("conflicting REF at duplicate selected locus")
                pending = (name, position, ref, tuple(sorted(set(pending[3]) | set(alts))))
            else:
                if pending:
                    yield pending
                pending = (name, position, ref, tuple(sorted(set(alts))))
            previous = key
        if pending:
            yield pending


def selection_windows(loci, width):
    if not 1 <= width <= 16384:
        raise ValueError("oracle tile width must be 1..16384")
    current, key = [], None
    for locus in loci:
        next_key = locus[0], locus[1] // width
        if current and next_key != key:
            yield current
            current = []
        key = next_key
        current.append(locus)
    if current:
        yield current


def sample_filter(header, sample=None, pool=False):
    groups = {}
    for group in header.get("RG", []):
        key = group.get("ID")
        if not key or key in groups:
            raise ValueError("invalid/duplicate alignment read group")
        groups[key] = group.get("SM") or None
    samples = set(value for value in groups.values() if value)
    if pool:
        return lambda read: True
    if sample is None:
        if not samples:
            return lambda read: True
        if len(samples) != 1 or None in groups.values():
            raise ValueError("ambiguous sample scope requires explicit sample or pool")
        sample = next(iter(samples))
    if sample not in samples:
        raise ValueError("requested sample is not declared in alignment header")

    def includes(read):
        try:
            group = read.get_tag("RG")
        except KeyError as error:
            raise ValueError("named sample requires an assignable RG tag") from error
        if group not in groups or groups[group] is None:
            raise ValueError("named sample requires a declared RG/SM")
        return groups[group] == sample
    return includes


def evidence_rows(bam, fasta, loci, tile_bases=1024, max_read_len=250,
                  max_record_bytes=1048576, sample=None, pool=False, stats=None):
    stats = stats if stats is not None else {}
    stats.update(record_visits=0, selected_loci=0, windows=0, sample_filtered_record_visits=0,
                 max_retained_loci=0, max_numeric_array_bytes=0)
    includes = sample_filter(bam.header.to_dict(), sample, pool)
    for window in selection_windows(loci, tile_bases):
        stats["windows"] += 1
        count = len(window)
        stats["max_retained_loci"] = max(stats["max_retained_loci"], count)
        name = window[0][0]
        positions = [locus[1] for locus in window]
        reference = fasta.fetch(name, positions[0], positions[-1] + 1).upper()
        refs = [reference[position - positions[0]] for position in positions]
        for locus, ref in zip(window, refs):
            if locus[2] is not None and locus[2] != ref:
                raise ValueError("selected SNV REF differs from FASTA")
        values = [array("Q", [0]) * count for _ in SCALARS]
        bq_hist, mq_hist = array("Q", [0]) * (count * 94), array("Q", [0]) * (count * 255)
        stats["max_numeric_array_bytes"] = max(stats["max_numeric_array_bytes"], count * (28 + 94 + 255) * 8)
        previous_start = -1
        for read in bam.fetch(name, positions[0], positions[-1] + 1):
            stats["record_visits"] += 1
            if read.query_length > max_read_len:
                raise ValueError("decoded read exceeds oracle max_read_len")
            cigar = read.cigartuples or ()
            # This post-decode lower-bound check never claims to bound native
            # decoder/header/aux allocations; the report labels that limitation.
            if read.query_length * 3 + len(cigar) * 8 + len(read.query_name or "") > max_record_bytes:
                raise ValueError("decoded record exceeds oracle record envelope")
            if read.is_unmapped:
                continue
            if read.reference_start < previous_start:
                raise ValueError("indexed alignments are not coordinate sorted")
            previous_start = read.reference_start
            if not includes(read):
                stats["sample_filtered_record_visits"] += 1
                continue
            failure = next((index for index, flag in enumerate((0x100, 0x800, 0x200, 0x400))
                            if read.flag & flag), None)
            if failure is None:
                failure = 4 if read.mapping_quality == 255 else (5 if read.mapping_quality < 20 else None)
            sequence = (read.query_sequence or "").upper()
            qualities = read.query_qualities
            reference_position, query_position = read.reference_start, 0
            for operation, length in cigar:
                if operation in (0, 7, 8):
                    begin = bisect_left(positions, reference_position)
                    end = bisect_left(positions, reference_position + length, begin)
                    for index in range(begin, end):
                        offset = query_position + positions[index] - reference_position
                        if offset >= len(sequence):
                            raise ValueError("CIGAR exceeds decoded sequence")
                        values[0][index] += 1
                        reason = failure
                        if reason is None:
                            values[1][index] += 1
                            quality = 255 if qualities is None else qualities[offset]
                            base = sequence[offset]
                            if base == "=":
                                base = refs[index]
                            reason = 6 if quality == 255 else (7 if quality < 20 else (8 if base not in "ACGT" else None))
                        if reason is not None:
                            values[19 + reason][index] += 1
                            continue
                        if quality > 93:
                            raise ValueError("base quality exceeds supported profile")
                        allele = "ACGT".index(base)
                        values[2][index] += 1
                        values[3 + allele][index] += 1
                        values[7 + allele * 2 + int(read.is_reverse)][index] += 1
                        values[15][index] += quality
                        values[16][index] += read.mapping_quality
                        values[17][index] += read.query_length - 1 - offset if read.is_reverse else offset
                        values[18][index] += read.query_length
                        bq_hist[index * 94 + quality] += 1
                        mq_hist[index * 255 + read.mapping_quality] += 1
                    reference_position += length
                    query_position += length
                elif operation in (2, 3):
                    reference_position += length
                elif operation in (1, 4):
                    query_position += length
                elif operation not in (5, 6):
                    raise ValueError("unsupported CIGAR operation")
        for index, locus in enumerate(window):
            histogram = lambda bins, width: ",".join(
                f"{quality}:{bins[index * width + quality]}" for quality in range(width)
                if bins[index * width + quality]) or "."
            stats["selected_loci"] += 1
            yield [name, locus[1] + 1, refs[index], ",".join(locus[3]) or ".",
                   *[column[index] for column in values], histogram(bq_hist, 94), histogram(mq_hist, 255)]
        # Drop completed numeric buffers before the next window allocates its
        # arrays. The last short window must not transiently retain both sets.
        del values, bq_hist, mq_hist


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("reference", "reference-fai", "alignments", "alignment-index", "selection", "stats"):
        parser.add_argument("--" + flag, required=True, type=Path)
    parser.add_argument("--selection-kind", choices=("sites", "regions"), required=True)
    parser.add_argument("--tile-bases", type=int, default=1024)
    parser.add_argument("--max-read-len", type=int, default=250)
    parser.add_argument("--max-record-bytes", type=int, default=1048576)
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument("--sample")
    scope.add_argument("--pool-samples", action="store_true")
    args = parser.parse_args()
    if pysam.__version__ != "0.23.3":
        parser.error("oracle is pinned to pysam==0.23.3")
    started, stats = time.perf_counter(), {}
    with pysam.FastaFile(str(args.reference), filepath_index=str(args.reference_fai)) as fasta, pysam.AlignmentFile(
            str(args.alignments), "r", reference_filename=str(args.reference),
            index_filename=str(args.alignment_index)) as bam:
        ranks = {name: rank for rank, name in enumerate(fasta.references)}
        lengths = dict(zip(fasta.references, fasta.lengths))
        loci = selection_loci(args.selection_kind, args.selection, ranks, lengths)
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(evidence_rows(bam, fasta, loci, args.tile_bases, args.max_read_len,
                                      args.max_record_bytes, args.sample, args.pool_samples, stats))
    sys.stdout.flush()
    stats.update(elapsed_seconds=time.perf_counter() - started, pysam=pysam.__version__,
                 htslib=pysam.__samtools_version__, tile_bases=args.tile_bases,
                 assurance="bounded numeric arrays; post-decode envelope; observed native decoder")
    args.stats.write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    main()
