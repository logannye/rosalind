"""Conservative, genotype-aware labels for supplied SNV candidates.

The file-backed constructor imports pysam lazily. ``from_records`` provides an
independent pure-Python entry point for fixtures. VCF records outside selected
intervals and their exclusion neighborhoods are streamed past, never retained.
"""

from __future__ import annotations

import bisect
import gzip
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping

BASES = frozenset("ACGT")
EXCLUSION_PADDING = 10


def _interval(value):
    if isinstance(value, Mapping):
        contig = value.get("contig", value.get("chrom"))
        start, end = value["start"], value["end"]
        extra = {key: value[key] for key in ("partition", "block_id") if key in value}
    else:
        contig, start, end = value[:3]
        extra = {}
    if not isinstance(contig, str) or not contig:
        raise ValueError("interval contig must be a nonempty string")
    if any(isinstance(x, bool) or not isinstance(x, int) for x in (start, end)) or not 0 <= start < end:
        raise ValueError("intervals require nonnegative half-open integer coordinates")
    return {"contig": contig, "start": start, "end": end, **extra}


class _Intervals:
    def __init__(self, intervals):
        grouped = defaultdict(list)
        for raw in intervals:
            item = _interval(raw)
            grouped[item["contig"]].append((item["start"], item["end"]))
        self.ranges = {}
        self.starts = {}
        for contig, ranges in grouped.items():
            merged = []
            for start, end in sorted(ranges):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.ranges[contig] = merged
            self.starts[contig] = [start for start, _ in merged]

    def overlaps(self, contig, start, end):
        ranges = self.ranges.get(contig, ())
        i = bisect.bisect_left(self.starts.get(contig, ()), end) - 1
        return i >= 0 and ranges[i][1] > start

    def contains(self, contig, pos0):
        return self.overlaps(contig, pos0, pos0 + 1)

    def as_dicts(self):
        return [
            {"contig": contig, "start": start, "end": end}
            for contig in sorted(self.ranges) for start, end in self.ranges[contig]
        ]


def _read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError("benchmark BED row has fewer than three columns")
            yield _interval((fields[0], int(fields[1]), int(fields[2])))


class TruthIndex:
    """Labels restricted to confident, selected, unambiguous SNV observations.

    ``positives`` contains assessed-domain ``(contig, pos0, ref, alt)`` keys.
    ``label`` returns ``(1|0|None, reason)``. Unknown rows are never negatives.
    Selected intervals may carry ``partition`` and ``block_id`` for denominator
    accounting. All coordinates passed to this class are zero-based, half-open.
    """

    def __init__(self, vcf_path, bed_path, intervals, reference_path=None):
        import pysam

        selected = [_interval(item) for item in intervals]
        selected_index = _Intervals(selected)
        self._reference_handle = pysam.FastaFile(str(reference_path)) if reference_path else None
        reference_fetch = self._reference_handle.fetch if self._reference_handle else None

        def records():
            with pysam.VariantFile(str(vcf_path)) as vcf:
                if len(vcf.header.samples) != 1:
                    raise ValueError("benchmark VCF must contain exactly one named sample")
                sample = next(iter(vcf.header.samples))
                for record in vcf:
                    end = max(record.start + len(record.ref or "N"), record.stop)
                    if not selected_index.overlaps(
                        record.contig, max(0, record.start - EXCLUSION_PADDING), end + EXCLUSION_PADDING
                    ):
                        continue
                    yield {
                        "contig": record.contig, "pos0": record.start,
                        "ref": record.ref, "alts": record.alts or (),
                        "gt": record.samples[sample].get("GT"), "end": end,
                        "filters": tuple(record.filter),
                    }

        try:
            self._initialize(records(), _read_bed(bed_path), selected, reference_fetch)
        except BaseException:
            self.close()
            raise

    @classmethod
    def from_records(cls, records: Iterable[Mapping], confident_intervals,
                     intervals=None, reference_fetch=None):
        """Build fixtures without pysam or files.

        Record keys: ``contig,pos0,ref,alts,gt``; optional ``end`` is exclusive
        (VCF INFO/END's numeric value), and optional ``filters`` lists FILTER IDs.
        ``reference_fetch(contig,start,end)`` follows pysam's zero-based API.
        """
        instance = cls.__new__(cls)
        instance._reference_handle = None
        confident = [_interval(item) for item in confident_intervals]
        selected = confident if intervals is None else [_interval(item) for item in intervals]
        instance._initialize(records, confident, selected, reference_fetch)
        return instance

    def _initialize(self, records, confident, selected, reference_fetch):
        self.intervals = list(selected)
        self._selected = _Intervals(self.intervals)
        self._confident = _Intervals(confident)
        self._reference_fetch = reference_fetch
        self._truth_refs = {}
        self.exclusion_counts = Counter()
        excluded = []
        proposed = set()
        self.records_considered = 0
        for record in records:
            contig, pos0, ref = record["contig"], record["pos0"], record.get("ref")
            if isinstance(pos0, bool) or not isinstance(pos0, int) or pos0 < 0:
                raise ValueError("benchmark position must be a nonnegative integer")
            span_end = max(pos0 + len(ref or "N"), record.get("end") or pos0 + 1)
            padded_start, padded_end = max(0, pos0 - EXCLUSION_PADDING), span_end + EXCLUSION_PADDING
            if not self._selected.overlaps(contig, padded_start, padded_end):
                continue
            self.records_considered += 1
            alts, gt = tuple(record.get("alts") or ()), record.get("gt")
            reason = None
            if (not isinstance(ref, str) or len(ref) != 1 or ref not in BASES
                    or not alts or any(not isinstance(a, str) or len(a) != 1 or a not in BASES or a == ref for a in alts)):
                reason = "complex_or_unsupported_record"
            elif not gt or any(allele is None for allele in gt):
                reason = "missing_genotype"
            elif any(isinstance(a, bool) or not isinstance(a, int) or a < 0 or a > len(alts) for a in gt):
                reason = "invalid_genotype"
            elif any(value not in ("PASS", ".") for value in record.get("filters", ())):
                reason = "filtered_truth_record"
            elif not self._reference_matches(contig, pos0, ref):
                reason = "unresolved_truth_reference"
            elif (contig, pos0) in self._truth_refs and self._truth_refs[(contig, pos0)] != ref:
                reason = "contradictory_truth_reference"
            if reason:
                excluded.append({"contig": contig, "start": padded_start, "end": padded_end})
                self.exclusion_counts[reason] += 1
                continue
            self._truth_refs[(contig, pos0)] = ref
            for index in set(gt):
                if index:
                    proposed.add((contig, pos0, ref, alts[index - 1]))
        self._excluded = _Intervals(excluded)
        self.excluded_intervals = self._excluded.as_dicts()
        self.positives = {
            key for key in proposed
            if self._selected.contains(key[0], key[1])
            and self._confident.contains(key[0], key[1])
            and not self._excluded.contains(key[0], key[1])
        }

    def _reference_matches(self, contig, pos0, ref):
        if self._reference_fetch is None:
            return True
        try:
            observed = self._reference_fetch(contig, pos0, pos0 + len(ref))
        except (KeyError, ValueError, IndexError):
            return False
        return isinstance(observed, str) and observed.upper() == ref

    def label(self, contig, pos0, ref, alt):
        if (not isinstance(pos0, int) or isinstance(pos0, bool) or pos0 < 0
                or not isinstance(ref, str) or ref not in BASES
                or not isinstance(alt, str) or alt not in BASES or ref == alt):
            return None, "unsupported_candidate"
        if not self._selected.contains(contig, pos0):
            return None, "outside_selected_intervals"
        if not self._confident.contains(contig, pos0):
            return None, "outside_benchmark_regions"
        if self._excluded.contains(contig, pos0):
            return None, "excluded_truth_event"
        if not self._reference_matches(contig, pos0, ref):
            return None, "candidate_reference_mismatch"
        if (contig, pos0) in self._truth_refs and self._truth_refs[(contig, pos0)] != ref:
            return None, "candidate_reference_mismatch"
        if (contig, pos0, ref, alt) in self.positives:
            return 1, "benchmark_supported_alt"
        return 0, "benchmark_region_unmatched_alt"

    def truth_keys(self, partition=None):
        if partition is None:
            return set(self.positives)
        selected = _Intervals(i for i in self.intervals if i.get("partition") == partition)
        return {key for key in self.positives if selected.contains(key[0], key[1])}

    def positive_counts_by_interval(self):
        return [
            {**interval, "n_positive": sum(
                key[0] == interval["contig"] and interval["start"] <= key[1] < interval["end"]
                for key in self.positives
            )}
            for interval in self.intervals
        ]

    @staticmethod
    def _intersections(index, contig, start, end):
        """Yield clipped, disjoint intervals without visiting individual bases."""
        ranges = index.ranges.get(contig, ())
        first = max(0, bisect.bisect_right(index.starts.get(contig, ()), start) - 1)
        for left, right in ranges[first:]:
            if left >= end:
                break
            clipped_start, clipped_end = max(start, left), min(end, right)
            if clipped_start < clipped_end:
                yield clipped_start, clipped_end

    def domain_report(self):
        """Return the geometric label domain and assessed truth denominator.

        Base counts use half-open interval intersections and subtraction. The
        assessed base domain is selected BED intersect benchmark BED minus the
        merged conservative event mask. ``assessed_truth_snvs`` counts distinct
        called ALT keys in that domain, not distinct loci or candidate rows.
        """
        blocks = []
        for interval in self.intervals:
            contig, start, end = interval["contig"], interval["start"], interval["end"]
            confident_bases = excluded_confident_bases = 0
            for left, right in self._intersections(self._confident, contig, start, end):
                confident_bases += right - left
                excluded_confident_bases += sum(
                    excluded_right - excluded_left
                    for excluded_left, excluded_right in self._intersections(
                        self._excluded, contig, left, right
                    )
                )
            blocks.append({
                **interval,
                "selected_bases": end - start,
                "benchmark_confident_bases": confident_bases,
                "excluded_confident_bases": excluded_confident_bases,
                "assessed_bases": confident_bases - excluded_confident_bases,
                "assessed_truth_snvs": sum(
                    key[0] == contig and start <= key[1] < end for key in self.positives
                ),
            })
        return {
            "exclusion_padding_bp": EXCLUSION_PADDING,
            "records_considered": self.records_considered,
            "exclusion_counts": dict(sorted(self.exclusion_counts.items())),
            "excluded_intervals": list(self.excluded_intervals),
            "blocks": blocks,
        }

    def close(self):
        handle = getattr(self, "_reference_handle", None)
        if handle is not None:
            handle.close()
            self._reference_handle = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
