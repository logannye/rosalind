"""Exact masked reduction for one admitted window, without a member matrix.

Native Arrow columns converted with as_py() supply Python integers, never floats.
One current-member observed-locus set prevents duplicate counting. All counters
and sums use checked uint64 arithmetic. Threshold/ALT rows repeat each locus's
denominator; this module never adds those repeated denominators into a total.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping

from contract import (Candidate, ContractError, MAX_CANDIDATES, MAX_MEMBERS,
                      UINT64_MAX, bounded_text, exact_uint64, normalize_thresholds)


def checked_add(left: int, right: int) -> int:
    result = exact_uint64(left, "sum operand") + exact_uint64(right, "sum operand")
    if result > UINT64_MAX:
        raise ContractError("uint64 reduction overflow")
    return result


def fraction(numerator: int, denominator: int) -> tuple[int | None, int | None]:
    return (numerator, denominator) if denominator else (None, None)


class WindowReducer:
    """Serial member reduction retaining only admitted candidate-window counters.

    Each selected member must have one begin/end pair, including a wholly absent
    member. observe() accepts native output column names: contig, pos (one-based),
    ref, callable_depth, a, c, g, t. Native stream completion and expected covered
    row counts remain the caller's responsibility; strict finish additionally
    rejects any candidate with fewer than n_requested observations.
    """

    def __init__(self, candidates: Iterable[Candidate], thresholds, n_requested: int):
        if type(n_requested) is not int or not 1 <= n_requested <= MAX_MEMBERS:
            raise ContractError(f"n_requested must be between 1 and {MAX_MEMBERS}")
        self.thresholds = normalize_thresholds(thresholds)
        self.n_requested = n_requested
        self.candidates = []
        for candidate in candidates:
            if not isinstance(candidate, Candidate):
                raise ContractError("reducer requires native-normalized Candidate objects")
            self.candidates.append(candidate)
            if len(self.candidates) > MAX_CANDIDATES:
                raise ContractError("candidate window exceeds supported envelope")
        self.candidates.sort()
        self._by_locus = defaultdict(list)
        self._rows = {}
        self._locus_order = {}
        for candidate in self.candidates:
            key = (candidate.contig, candidate.pos)
            entries = self._by_locus[key]
            if any(prior.ref != candidate.ref for prior in entries):
                raise ContractError("conflicting REF annotations at a candidate locus")
            if candidate in entries:
                raise ContractError("duplicate candidate row; normalize candidates first")
            entries.append(candidate)
            self._locus_order[key] = (candidate.contig_id, candidate.position0)
            for threshold in self.thresholds:
                self._rows[(candidate, threshold)] = {
                    **candidate.as_dict(), "min_callable_depth": threshold,
                    "n_requested": n_requested, "n_observed": 0,
                    "n_depth_eligible": 0, "n_alt_supported": 0,
                    "total_callable_depth": 0, "total_alt_count": 0,
                    "eligible_callable_depth": 0, "eligible_alt_count": 0,
                }
        self._members = set()
        self._current = None
        self._observed = set()
        self._last = None
        self._finished = False
        self._failed = False

    def begin_member(self, specimen_id: str) -> None:
        if self._finished or self._failed or self._current is not None:
            raise ContractError("cannot begin member in current reducer state")
        bounded_text(specimen_id, "specimen_id")
        if specimen_id in self._members or len(self._members) >= self.n_requested:
            raise ContractError("duplicate or excess member in reduction")
        self._members.add(specimen_id)
        self._current = specimen_id
        self._observed.clear()
        self._last = None

    def observe(self, row: Mapping) -> None:
        if self._current is None or self._finished or self._failed:
            raise ContractError("begin_member must precede observations")
        # Any rejected native row poisons this reducer. A caller cannot catch an
        # overflow/integrity error and later publish the surviving partial sums.
        self._failed = True
        if not isinstance(row, Mapping):
            raise ContractError("native observation must be a row mapping")
        pos = exact_uint64(row.get("pos"), "native POS")
        key = (row.get("contig"), pos)
        if key not in self._by_locus:
            raise ContractError("native row falls outside the admitted candidate window")
        if key in self._observed:
            raise ContractError("duplicate observed locus within a member")
        order = self._locus_order[key]
        if self._last is not None and order <= self._last:
            raise ContractError("native rows must follow reference dictionary/position order")
        candidates = self._by_locus[key]
        if row.get("ref") != candidates[0].ref:
            raise ContractError("candidate REF does not match observed saved evidence")
        depth = exact_uint64(row.get("callable_depth"), "callable_depth")
        alleles = {base: exact_uint64(row.get(base.lower()), f"{base} count") for base in "ACGT"}
        allele_total = 0
        for count in alleles.values():
            allele_total = checked_add(allele_total, count)
        if allele_total != depth:
            raise ContractError("callable depth must equal the sum of A/C/G/T counts")
        # Compute all updates before committing this observation, so an overflow
        # does not leave a partially updated scientific row.
        updates = []
        for candidate in candidates:
            alt_count = alleles[candidate.alt]
            for threshold in self.thresholds:
                target = self._rows[(candidate, threshold)]
                eligible = depth >= threshold
                deltas = {"n_observed": 1, "total_callable_depth": depth,
                          "total_alt_count": alt_count, "n_depth_eligible": int(eligible),
                          "n_alt_supported": int(eligible and alt_count >= 1),
                          "eligible_callable_depth": depth if eligible else 0,
                          "eligible_alt_count": alt_count if eligible else 0}
                updates.append((target, {key: checked_add(target[key], delta)
                                         for key, delta in deltas.items()}))
        for target, update in updates:
            target.update(update)
        self._observed.add(key)
        self._last = order
        self._failed = False

    def end_member(self) -> int:
        if self._current is None or self._finished or self._failed:
            raise ContractError("no active member to complete")
        observed = len(self._observed)
        self._current = None
        self._observed.clear()
        self._last = None
        return observed

    def finish(self, *, missingness: str = "partial") -> list[dict]:
        if self._failed:
            raise ContractError("failed observation prevents summary publication")
        if missingness not in ("strict", "partial"):
            raise ContractError("missingness must be strict or partial")
        if self._current is not None or len(self._members) != self.n_requested:
            raise ContractError("every requested member must complete before summary publication")
        if missingness == "strict" and any(
                row["n_observed"] != self.n_requested for row in self._rows.values()):
            raise ContractError("strict missingness refuses absent member/candidate evidence")
        result = []
        for candidate in self.candidates:
            for threshold in self.thresholds:
                row = dict(self._rows[(candidate, threshold)])
                row["alt_fraction_numerator"], row["alt_fraction_denominator"] = fraction(
                    row["total_alt_count"], row["total_callable_depth"])
                row["eligible_alt_fraction_numerator"], row["eligible_alt_fraction_denominator"] = fraction(
                    row["eligible_alt_count"], row["eligible_callable_depth"])
                row["reference_validation"] = ("checked_on_observed_rows" if row["n_observed"]
                                               else "unavailable")
                result.append(row)
        self._finished = True
        return result
