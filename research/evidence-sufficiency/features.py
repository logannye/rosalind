"""Deterministic F0/F1 features from one complete native evidence row.

Inputs use ``src/evidence/encoding.rs`` column names. Integer observations remain
the caller's responsibility to retain: this module deliberately creates model
floats and never mutates the source row. No coordinate or specimen identifier is
included in a feature vector. Omitted evidence is an error, not measured zero.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

BASES = "ACGT"
U64_MAX = (1 << 64) - 1
FILTER_NAMES = (
    "secondary", "supplementary", "qc_fail", "duplicate", "unavailable_mapq",
    "low_mapq", "unavailable_base_quality", "low_base_quality", "ambiguous_base",
)
SUM_NAMES = ("base_quality", "mapping_quality", "read_position", "read_length")


def _base(value: object, name: str) -> str:
    if isinstance(value, int) and not isinstance(value, bool):
        value = chr(value) if 0 <= value <= 127 else None
    if not isinstance(value, str) or value not in BASES or len(value) != 1:
        raise ValueError(f"{name} must be one explicit A/C/G/T base")
    return value


def _u64(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= U64_MAX:
        raise ValueError(f"{name} must be an exact uint64 integer")
    return value


def _required(row: Mapping, name: str) -> int:
    if name not in row:
        raise ValueError(f"required native evidence column is absent: {name}")
    return _u64(row[name], name)


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _nullable(features: dict, name: str, value: float | None) -> None:
    features[name] = value
    features[f"{name}_undefined"] = float(value is None)


def _difference(left: float | None, right: float | None) -> float | None:
    return None if left is None or right is None else left - right


def feature_vectors(raw_row: dict, alt: str) -> dict[str, dict[str, float | None]]:
    """Return stable-schema ``{'f0': ..., 'f1': ...}`` dictionaries.

    A complete F1 native projection is required even when only F0 is subsequently
    selected. Null derived ratios carry explicit ``*_undefined`` indicators;
    missing source fields fail instead of being imputed here. Native ``ref`` may
    be its Arrow one-character string or the equivalent unsigned ASCII value.
    """
    if not isinstance(raw_row, Mapping) or raw_row.get("observed") is False:
        raise ValueError("a measured native evidence row is required")
    ref = _base(raw_row.get("ref"), "ref")
    alt = _base(alt, "alt")
    if ref == alt:
        raise ValueError("ALT must differ from REF")
    depth = _required(raw_row, "callable_depth")
    prefilter = _required(raw_row, "prefilter_depth")
    aligned = _required(raw_row, "aligned_depth")
    counts = [_required(raw_row, base.lower()) for base in BASES]
    if sum(counts) != depth or not prefilter >= aligned >= depth:
        raise ValueError("inconsistent native depth/allele counts")
    filters = [_required(raw_row, "filtered_" + name) for name in FILTER_NAMES]
    if prefilter != aligned + sum(filters[:6]) or aligned != depth + sum(filters[6:]):
        raise ValueError("exclusive filter counters do not reconcile with depths")
    strands = [
        (_required(raw_row, base.lower() + "_fwd"), _required(raw_row, base.lower() + "_rev"))
        for base in BASES
    ]
    if any(sum(pair) != count for pair, count in zip(strands, counts)):
        raise ValueError("strand counts do not reconcile with allele counts")
    global_sums = {name: _required(raw_row, name + "_sum") for name in SUM_NAMES}
    allele_sums = {}
    for name in SUM_NAMES:
        column = "allele_" + name + "_sum"
        values = raw_row.get(column)
        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)) or len(values) != 4:
            raise ValueError(f"{column} must contain four A/C/G/T uint64 values")
        checked = [_u64(value, column) for value in values]
        if sum(checked) != global_sums[name]:
            raise ValueError(f"{column} does not reconcile with its global sum")
        if any(count == 0 and value != 0 for count, value in zip(counts, checked)):
            raise ValueError(f"{column} has a nonzero sum for an unobserved allele")
        allele_sums[name] = checked

    ref_index, alt_index = BASES.index(ref), BASES.index(alt)
    ref_count, alt_count = counts[ref_index], counts[alt_index]
    other_count = depth - ref_count - alt_count
    f0: dict[str, float | None] = {}
    for name, value in (
        ("callable_depth", depth), ("ref_count", ref_count),
        ("alt_count", alt_count), ("other_count", other_count),
    ):
        f0[name] = float(value)
        f0[name + "_log1p"] = math.log1p(value)
    for name, count in (("ref", ref_count), ("alt", alt_count), ("other", other_count)):
        _nullable(f0, name + "_fraction", _ratio(count, depth))
    for base in BASES:
        f0["ref_" + base] = float(ref == base)
        f0["alt_" + base] = float(alt == base)

    f1 = dict(f0)
    for base, pair, count in zip(BASES, strands, counts):
        _nullable(f1, base.lower() + "_forward_fraction", _ratio(pair[0], count))
    ref_strand = _ratio(strands[ref_index][0], ref_count)
    alt_strand = _ratio(strands[alt_index][0], alt_count)
    _nullable(f1, "ref_forward_fraction", ref_strand)
    _nullable(f1, "alt_forward_fraction", alt_strand)
    _nullable(f1, "alt_minus_ref_forward_fraction", _difference(alt_strand, ref_strand))
    for name in SUM_NAMES:
        _nullable(f1, name + "_mean", _ratio(global_sums[name], depth))
        ref_mean = _ratio(allele_sums[name][ref_index], ref_count)
        alt_mean = _ratio(allele_sums[name][alt_index], alt_count)
        _nullable(f1, "ref_" + name + "_mean", ref_mean)
        _nullable(f1, "alt_" + name + "_mean", alt_mean)
        _nullable(f1, "alt_minus_ref_" + name + "_mean", _difference(alt_mean, ref_mean))
    for name, count in zip(FILTER_NAMES, filters):
        _nullable(f1, "filtered_" + name + "_fraction", _ratio(count, prefilter))
    _nullable(f1, "aligned_fraction", _ratio(aligned, prefilter))
    _nullable(f1, "callable_fraction", _ratio(depth, prefilter))
    return {"f0": f0, "f1": f1}
