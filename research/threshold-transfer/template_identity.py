"""Bounded, source-local template-identity audit for explicit no-RG intake.

Every physical alignment record must pass ``observe`` before downsampling. This
audit never removes records, repairs QNAMEs, or invents biological read groups.
It tests consistency of exact source-local names within the admitted region;
it does not prove that globally distinct molecules always have distinct names.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import struct
from typing import Callable

from preflight import validate_cluster_qname


PAIRED = 0x1
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
REVERSE = 0x10
MATE_REVERSE = 0x20
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0x100
SUPPLEMENTARY = 0x800


@dataclass(slots=True)
class _PrimaryEnd:
    signature: bytes
    reference_id: int
    reference_start: int
    mate_reference_id: int
    mate_reference_start: int
    flag: int


@dataclass(slots=True)
class _Template:
    first: _PrimaryEnd | None = None
    second: _PrimaryEnd | None = None


def _text(digest, value: str | None) -> None:
    """Length-prefix nullable strings so payload boundaries cannot alias."""
    if value is None:
        digest.update(struct.pack("<q", -1))
    else:
        data = value.encode("ascii")
        digest.update(struct.pack("<q", len(data)))
        digest.update(data)


def _primary_signature(record) -> bytes:
    """Hash full primary core fields, excluding optional auxiliary tags.

    QUAL, sequence, CIGAR, full FLAG and all coordinates are included. Tags are
    excluded deliberately: ordering and nonidentity annotations can vary. A
    repeated signature is not described as byte-identical source BAM records.
    """
    digest = hashlib.sha256()
    digest.update(b"rosalind-primary-end-signature-v1\0")
    end = record.reference_end
    digest.update(struct.pack("<8q", int(record.flag), int(record.reference_id),
        int(record.reference_start), -1 if end is None else int(end),
        int(record.mapping_quality), int(record.next_reference_id),
        int(record.next_reference_start), int(record.template_length)))
    cigar = record.cigartuples
    digest.update(struct.pack("<q", -1 if cigar is None else len(cigar)))
    if cigar is not None:
        for operation, length in cigar:
            digest.update(struct.pack("<2q", int(operation), int(length)))
    _text(digest, record.query_sequence)
    qualities = record.query_qualities
    if qualities is None:
        digest.update(struct.pack("<q", -1))
    else:
        values = bytes(qualities)
        digest.update(struct.pack("<q", len(values)))
        digest.update(values)
    return digest.digest()


class TemplateIdentityAudit:
    """Fail-closed audit keyed by exact QNAME in one source/synthetic group.

    ``max_templates`` bounds retained dictionary entries. Compact digests replace
    sequence/quality payloads, and a caller-side process-tree RSS limit remains
    necessary because a template-count bound is not an exact byte bound.
    ``validate_qname`` must return the entire unchanged name when supplied.
    """

    def __init__(self, *, source_identity: str, synthetic_group: str,
                 max_templates: int = 8_000_000,
                 validate_qname: Callable[[str], str] | None = None):
        if not isinstance(source_identity, str) or not source_identity:
            raise ValueError("nonempty source identity required")
        if not isinstance(synthetic_group, str) or not synthetic_group:
            raise ValueError("nonempty technical synthetic group required")
        if type(max_templates) is not int or max_templates < 1:
            raise ValueError("max_templates must be a positive integer")
        self.source_identity = source_identity
        self.synthetic_group = synthetic_group
        self.max_templates = max_templates
        self._validator = validate_qname or validate_cluster_qname
        self._templates: dict[str, _Template] = {}
        self._finished = False
        self._failure = None
        self._records = 0
        self._primary_records = 0
        self._secondary_records = 0
        self._supplementary_records = 0
        self._nonprimary_records = 0
        self._repeated_primary = 0
        self._unique_primary_ends = 0
        self._pairs_observed = 0
        self._mapped_pairs_checked = 0
        self._pairs_with_unmapped_end = 0
        self._stream_digest = hashlib.sha256(b"rosalind-template-audit-record-order-v1\0")
        self._receipt = None

    def _fail(self, reason: str, name: str | None = None):
        suffix = ""
        if name is not None:
            suffix = "; qname_sha256=" + hashlib.sha256(name.encode()).hexdigest()
        self._failure = reason + suffix
        raise ValueError(self._failure)

    def _pair(self, first: _PrimaryEnd, second: _PrimaryEnd, name: str) -> None:
        """Check primary-mate consistency without imposing proper-pair/TLEN rules."""
        first_unmapped, second_unmapped = bool(first.flag & UNMAPPED), bool(second.flag & UNMAPPED)
        if bool(first.flag & MATE_UNMAPPED) != second_unmapped or bool(second.flag & MATE_UNMAPPED) != first_unmapped:
            self._fail("observed primary mate unmapped flags disagree", name)
        if not first_unmapped and bool(second.flag & MATE_REVERSE) != bool(first.flag & REVERSE):
            self._fail("observed first-mate strand disagrees with second mate flag", name)
        if not second_unmapped and bool(first.flag & MATE_REVERSE) != bool(second.flag & REVERSE):
            self._fail("observed second-mate strand disagrees with first mate flag", name)
        if not first_unmapped and not second_unmapped:
            if ((first.mate_reference_id, first.mate_reference_start) != (second.reference_id, second.reference_start)
                    or (second.mate_reference_id, second.mate_reference_start) != (first.reference_id, first.reference_start)):
                self._fail("observed mapped primary mates have nonreciprocal reference/start coordinates", name)
            self._mapped_pairs_checked += 1
        else:
            self._pairs_with_unmapped_end += 1
        self._pairs_observed += 1

    def observe(self, record) -> None:
        """Validate one record, retaining only compact template consistency state."""
        if self._finished:
            raise RuntimeError("template audit already finished")
        if self._failure is not None:
            raise RuntimeError("template audit is invalid after an earlier failure")
        name = record.query_name
        try:
            validated = self._validator(name)
        except (TypeError, ValueError) as exc:
            self._fail("invalid full source QNAME: " + str(exc))
        if validated != name or not isinstance(name, str):
            self._fail("QNAME validator must preserve the complete exact string")
        if record.has_tag("RG"):
            self._fail("no-RG template audit rejects an original RG tag", name)
        flag = int(record.flag)
        if flag < 0 or flag > 0xffff:
            self._fail("record FLAG is outside BAM range", name)
        primary = not (flag & (SECONDARY | SUPPLEMENTARY))
        if primary and (not flag & PAIRED or bool(flag & READ1) == bool(flag & READ2)):
            self._fail("primary record requires paired status and exactly one read-end flag", name)
        if primary and not flag & UNMAPPED:
            if (record.reference_id < 0 or record.reference_start < 0 or record.reference_end is None
                    or record.reference_end <= record.reference_start):
                self._fail("mapped primary end lacks reference/start/end coordinates", name)
        state = self._templates.get(name)
        if state is None:
            if len(self._templates) >= self.max_templates:
                self._fail("template identity audit max_templates exceeded; no entries pruned")
            state = _Template()
            self._templates[name] = state
        self._records += 1
        self._secondary_records += bool(flag & SECONDARY)
        self._supplementary_records += bool(flag & SUPPLEMENTARY)
        # Hash every physical record's name/FLAG in original iteration order.
        # The receipt does not claim to hash every source record byte.
        _text(self._stream_digest, name)
        self._stream_digest.update(struct.pack("<H", flag))
        if not primary:
            self._nonprimary_records += 1
            return
        self._primary_records += 1
        try:
            signature = _primary_signature(record)
        except (TypeError, ValueError, UnicodeError, struct.error) as exc:
            self._fail("primary-end signature fields are invalid: " + str(exc), name)
        existing = state.first if flag & READ1 else state.second
        if existing is not None:
            if existing.signature != signature:
                self._fail("conflicting primary-end identity for exact source-local QNAME", name)
            self._repeated_primary += 1
            return
        current = _PrimaryEnd(signature, int(record.reference_id), int(record.reference_start),
                              int(record.next_reference_id), int(record.next_reference_start), flag)
        self._unique_primary_ends += 1
        if flag & READ1:
            state.first = current
        else:
            state.second = current
        if state.first is not None and state.second is not None:
            self._pair(state.first, state.second, name)

    def finish(self) -> dict:
        """Return a compact audit receipt; an invalid/incomplete audit cannot pass."""
        if self._failure is not None:
            raise RuntimeError("template audit failed: " + self._failure)
        if self._receipt is not None:
            return dict(self._receipt)
        if not self._records:
            self._fail("cannot admit an empty template identity audit")
        self._finished = True
        single_primary = no_primary = 0
        for state in self._templates.values():
            ends = int(state.first is not None) + int(state.second is not None)
            single_primary += ends == 1
            no_primary += ends == 0
        self._receipt = {
            "schema": "rosalind-source-local-template-audit-v1", "passed": True,
            "source_identity": self.source_identity, "technical_synthetic_group": self.synthetic_group,
            "template_key": "exact full unmodified seven-field Illumina QNAME within pinned source and technical group",
            "physical_records_observed": self._records,
            "distinct_exact_qnames": len(self._templates), "max_templates": self.max_templates,
            "primary_records": self._primary_records, "secondary_records": self._secondary_records,
            "supplementary_records": self._supplementary_records,
            "nonprimary_records": self._nonprimary_records,
            "secondary_and_supplementary_counters_may_overlap": True,
            "unique_primary_ends": self._unique_primary_ends,
            "repeated_primary_end_signatures": self._repeated_primary,
            "primary_signature_fields": "Full FLAG, reference/start/end, MAPQ, mate reference/start, TLEN, CIGAR, sequence and qualities; optional tags excluded",
            "templates_with_both_primary_mates": self._pairs_observed,
            "mapped_primary_pairs_reciprocity_checked": self._mapped_pairs_checked,
            "primary_pairs_with_unmapped_end": self._pairs_with_unmapped_end,
            "templates_with_one_primary_mate_not_observed_in_scope": single_primary,
            "templates_with_no_primary_alignment_observed_in_scope": no_primary,
            "physical_records_deduplicated": 0, "qnames_rewritten": 0,
            "original_rg_tags_allowed": False,
            "record_order_qname_flag_sha256": self._stream_digest.hexdigest(),
            "receipt_scope": "Every admitted alignment offered by the caller from its declared bounded region, before downsampling; reconcile this count with caller traversal/exclusion counts; missing mates may be outside scope but that is not inferred",
            "memory_scope": "Bounded template-entry dictionary with compact digests; external process-tree RSS admission remains required",
            "limitations": ["Consistent identical primary signatures cannot exclude physically distinct molecules sharing the same full name",
                "Names are scoped to this source and technical group, not asserted globally unique",
                "Secondary and supplementary payloads are not subjected to primary-end collision checks",
                "Technical grouping supplies no biological lane, library or platform-unit metadata"],
        }
        return dict(self._receipt)
