"""Internal research contracts over metadata validated by native ``inspect``.

No VCF/BED parsing or source/reference reopening occurs here. The native coverage
adapter owns selection parsing and normalization. Source reuse identity and
cross-specimen scientific comparison are deliberately separate relations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import unicodedata
from typing import Any, Iterable, Mapping, Sequence

UINT64_MAX = (1 << 64) - 1
MAX_MEMBERS = 128
MAX_CANDIDATES = 100_000
MAX_THRESHOLDS = 8
REQUIRED_FIELDS = 3
REQUIRED_FIELD_NAMES = ("depths", "alleles")


class ContractError(ValueError):
    """A research input is outside the supported scientific/metadata envelope."""


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False)


def identity(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def exact_uint64(value: Any, name: str) -> int:
    # Do not coerce strings, floats, Arrow scalars, or booleans into counts. Arrow
    # consumers must use as_py(), which preserves native uint64 as Python int.
    if type(value) is not int or not 0 <= value <= UINT64_MAX:
        raise ContractError(f"{name} must be an exact uint64 integer")
    return value


def bounded_text(value: Any, name: str, max_bytes: int = 256) -> str:
    if (not isinstance(value, str) or not value.strip()
            or len(value.encode("utf-8")) > max_bytes
            or any(unicodedata.category(c).startswith("C") for c in value)):
        raise ContractError(f"{name} must be nonempty bounded text without control characters")
    return value


def _positive_limit(value: int, ceiling: int, name: str) -> int:
    if type(value) is not int or not 1 <= value <= ceiling:
        raise ContractError(f"{name} must be between 1 and {ceiling}")
    return value


def _object_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


@dataclass(frozen=True)
class Member:
    specimen_id: str
    dataset: Path
    metadata: Mapping[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"specimen_id": self.specimen_id, "dataset": str(self.dataset),
                **self.metadata}


def load_members(path: str | Path, *, max_manifest_bytes: int = 1 << 20,
                 max_members: int = MAX_MEMBERS) -> list[Member]:
    """Read an explicitly enumerated, byte-bounded one-leaf-per-specimen manifest.

    Paths are resolved against the manifest's directory, never discovered. The
    native reader later verifies each leaf. IDs sort lexically; group, subject,
    and timepoint are bounded assertions without inferred biological meaning.
    """
    _positive_limit(max_members, MAX_MEMBERS, "max_members")
    _positive_limit(max_manifest_bytes, 16 << 20, "max_manifest_bytes")
    path = Path(path).resolve(strict=True)
    with path.open("rb") as stream:
        raw = stream.read(max_manifest_bytes + 1)
    if len(raw) > max_manifest_bytes:
        raise ContractError("member manifest exceeds max_manifest_bytes")
    try:
        document = json.loads(raw, object_pairs_hook=_object_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"invalid member manifest JSON: {error}") from error
    if (not isinstance(document, dict) or set(document) != {"version", "members"}
            or type(document["version"]) is not int or document["version"] != 1):
        raise ContractError("member manifest requires version=1 and members only")
    entries = document["members"]
    if not isinstance(entries, list) or not 1 <= len(entries) <= max_members:
        raise ContractError(f"member manifest requires 1..{max_members} members")
    members, seen = [], set()
    for entry in entries:
        if (not isinstance(entry, dict) or not {"specimen_id", "dataset"} <= set(entry)
                or set(entry) - {"specimen_id", "dataset", "group", "subject", "timepoint"}):
            raise ContractError("each member requires specimen_id and exactly one dataset path")
        specimen = bounded_text(entry["specimen_id"], "specimen_id")
        if specimen in seen:
            raise ContractError(f"duplicate specimen_id: {specimen}")
        seen.add(specimen)
        dataset_text = bounded_text(entry["dataset"], "dataset path", 32768)
        dataset = Path(dataset_text)
        if not dataset.is_absolute():
            dataset = path.parent / dataset
        dataset = dataset.resolve(strict=True)
        if not dataset.is_file():
            raise ContractError("dataset must identify an existing portable leaf manifest file")
        metadata = {key: bounded_text(entry[key], key) for key in
                    ("group", "subject", "timepoint") if key in entry}
        members.append(Member(specimen, dataset, metadata))
    return sorted(members, key=lambda member: member.specimen_id)


def _named_scope(scope: Any) -> dict:
    if not isinstance(scope, dict) or scope.get("version") != 1 or scope.get("mode") != "named":
        raise ContractError("each specimen requires a single named analysis-sample scope")
    selected = bounded_text(scope.get("selected_sample"), "selected sample", 32768)
    samples, groups = scope.get("declared_samples"), scope.get("read_groups")
    if (scope.get("allows_unassigned") is not False or not isinstance(samples, list)
            or not isinstance(groups, list) or selected not in samples):
        raise ContractError("inconsistent named sample scope")
    if samples != sorted(set(samples)):
        raise ContractError("sample scope must preserve native canonical ordering")
    declared, group_ids = set(), []
    for group in groups:
        if not isinstance(group, dict):
            raise ContractError("invalid resolved read group")
        group_ids.append(bounded_text(group.get("id"), "read group ID", 32768))
        if group.get("sample") is not None:
            declared.add(bounded_text(group["sample"], "read group sample", 32768))
    if group_ids != sorted(set(group_ids)) or sorted(declared) != samples:
        raise ContractError("inconsistent resolved read-group/sample metadata")
    return scope


def _contigs(contigs: Any) -> list[dict]:
    if not isinstance(contigs, list) or not contigs:
        raise ContractError("an ordered nonempty reference dictionary is required")
    seen = set()
    for index, contig in enumerate(contigs):
        if not isinstance(contig, dict):
            raise ContractError("invalid reference dictionary entry")
        name = bounded_text(contig.get("name"), "contig name", 4096)
        if (type(contig.get("id")) is not int or contig["id"] != index or name in seen
                or type(contig.get("length")) is not int
                or not 1 <= contig["length"] <= (1 << 32) - 1):
            raise ContractError("invalid ordered reference dictionary")
        seen.add(name)
    return contigs


def _sources(descriptor: Mapping) -> dict[str, dict]:
    values = descriptor.get("sources")
    if not isinstance(values, list) or not 1 <= len(values) <= 256:
        raise ContractError("invalid bounded source identities")
    sources = {}
    for source in values:
        if not isinstance(source, dict):
            raise ContractError("invalid source identity")
        role = bounded_text(source.get("role"), "source role", 128)
        digest = source.get("blake3")
        if (role in sources or not isinstance(digest, str) or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)):
            raise ContractError("invalid or repeated source identity")
        exact_uint64(source.get("bytes"), "source bytes")
        sources[role] = source
    if "alignments" not in sources:
        raise ContractError("alignment content identity is required")
    return sources


def source_scope_identity(descriptor: Mapping) -> str:
    """Duplicate-specimen guard, independent of filters, fields, and selection.

    The alignment content and complete resolved scope identify this relation.
    Different indexes or extraction profiles cannot disguise the same source
    and scope as an independent specimen. Original paths are provenance only.
    """
    source = _sources(descriptor)["alignments"]
    return identity({"alignments": {key: source[key] for key in ("blake3", "bytes")},
                     "sample_scope": _named_scope(descriptor.get("sample_scope"))})


def comparison_profile(descriptor: Mapping) -> dict:
    """Extract comparison semantics; call only after native metadata validation.

    Equal byte identity entails equal representation without guessing a format
    from the recorded pathname. Descriptor v1 does not store a format enum, and
    the engine detects reference formats from bytes. No FASTA/rref equivalence
    is introduced and no original reference path is opened here.
    """
    if (descriptor.get("version") != 1
            or descriptor.get("semantics") != "exact-read-summary-v1"
            or descriptor.get("counting_unit") != "read"
            or descriptor.get("sampling") != "none"
            or descriptor.get("has_reference") is not True):
        raise ContractError("unsupported evidence semantics or unavailable reference")
    fields = descriptor.get("fields")
    if type(fields) is not int or fields & ~127 or fields & REQUIRED_FIELDS != REQUIRED_FIELDS:
        raise ContractError("every member must physically supply depths and alleles")
    if (descriptor.get("schema_version") != (1 if fields == 63 else 2)
            or descriptor.get("fields_version") != (2 if fields & 64 else 1)):
        raise ContractError("unsupported or inconsistent native schema/field-mask version")
    profile = descriptor.get("profile")
    names = {"min_mapq", "min_base_quality", "exclude_secondary",
             "exclude_supplementary", "exclude_qc_fail", "exclude_duplicates"}
    if not isinstance(profile, dict) or set(profile) != names:
        raise ContractError("explicit complete filter profile is required")
    for key, maximum in (("min_mapq", 254), ("min_base_quality", 93)):
        if type(profile[key]) is not int or not 0 <= profile[key] <= maximum:
            raise ContractError(f"invalid explicit {key}")
    if any(type(profile[key]) is not bool for key in names if key.startswith("exclude_")):
        raise ContractError("flag exclusions must be explicit booleans")
    _named_scope(descriptor.get("sample_scope"))
    sources = _sources(descriptor)
    role = "reference" if "reference" in sources else "cram-reference"
    if role not in sources:
        raise ContractError("effective analysis-reference content identity is required")
    reference = {key: sources[role][key] for key in ("role", "blake3", "bytes")}
    reference["representation_policy"] = "exact_source_bytes_no_equivalence"
    return {"version": 1, "semantics": descriptor["semantics"], "counting_unit": "read",
            "sampling": "none", "profile": dict(profile),
            "contigs": _contigs(descriptor.get("contigs")), "effective_reference": reference,
            "required_fields": list(REQUIRED_FIELD_NAMES), "required_field_mask": REQUIRED_FIELDS}


def validate_comparison(members_and_descriptors: Iterable[tuple[Member, Mapping]]) -> dict:
    """Check comparison compatibility and disallow aliasing one source/scope.

    The returned comparison_sha256 is an internal research identity. It does
    not replace or alter native compatibility_blake3, which remains a separate
    same-source reuse key. Caller retains each full descriptor's sources/scope.
    """
    comparison, source_scopes, specimens = None, set(), set()
    for member, descriptor in members_and_descriptors:
        if len(specimens) >= MAX_MEMBERS:
            raise ContractError(f"at most {MAX_MEMBERS} members are supported")
        specimen = bounded_text(member.specimen_id, "specimen_id")
        if specimen in specimens:
            raise ContractError(f"duplicate specimen_id: {specimen}")
        specimens.add(specimen)
        current = comparison_profile(descriptor)
        if comparison is None:
            comparison = current
        elif current != comparison:
            raise ContractError(f"incompatible cross-specimen comparison contract: {specimen}")
        source_scope = source_scope_identity(descriptor)
        if source_scope in source_scopes:
            raise ContractError("same alignment source and sample scope entered under multiple specimen aliases")
        source_scopes.add(source_scope)
    if comparison is None:
        raise ContractError("at least one specimen is required")
    return {**comparison, "comparison_sha256": identity(comparison)}


@dataclass(frozen=True, order=True)
class Candidate:
    contig_id: int
    position0: int
    alt: str
    contig: str = field(compare=False)
    ref: str = field(compare=False)

    @property
    def pos(self) -> int:
        return self.position0 + 1

    def as_dict(self) -> dict:
        return {"contig": self.contig, "pos": self.pos, "ref": self.ref, "alt": self.alt}


def normalized_candidates(selection: Mapping, contigs: list[dict], *,
                          max_candidates: int = MAX_CANDIDATES) -> list[Candidate]:
    """Expand native-normalized SNV sites into canonical distinct ALT rows.

    Input is DescriptorSelection from dataset coverage, not a raw selection
    file. Identical duplicates collapse, multiple ALTs are unioned and sorted,
    and contradictory REF assignments are refused. Empty SNV lists are valid.
    The native position remains zero-based internally; output POS adds one.
    """
    _positive_limit(max_candidates, MAX_CANDIDATES, "max_candidates")
    _contigs(contigs)
    sites = selection.get("sites")
    if not isinstance(sites, list):
        raise ContractError("native normalized SNV sites are required")
    if len(sites) > max_candidates:
        raise ContractError("native candidate metadata exceeds max_candidates")
    result, references = {}, {}
    for site in sites:
        if not isinstance(site, dict):
            raise ContractError("invalid native candidate site")
        contig_id, position = site.get("contig"), site.get("position")
        if (type(contig_id) is not int or not 0 <= contig_id < len(contigs)
                or type(position) is not int or not 0 <= position < contigs[contig_id]["length"]):
            raise ContractError("native candidate outside reference dictionary")
        reference, alternates = site.get("reference"), site.get("alternates")
        if (type(reference) is not int or reference not in b"ACGT"
                or not isinstance(alternates, list) or not 1 <= len(alternates) <= 3):
            raise ContractError("only A/C/G/T SNV candidates are supported")
        locus = (contig_id, position)
        if locus in references and references[locus] != reference:
            raise ContractError("conflicting REF annotations at a candidate locus")
        references[locus] = reference
        for alternate in alternates:
            if type(alternate) is not int or alternate not in b"ACGT" or alternate == reference:
                raise ContractError("ALT must be a distinct A/C/G/T base")
            key = (*locus, alternate)
            result[key] = Candidate(contig_id, position, chr(alternate),
                                    contigs[contig_id]["name"], chr(reference))
            if len(result) > max_candidates:
                raise ContractError("candidate ALT rows exceed max_candidates")
    return sorted(result.values())


def normalize_thresholds(values: Sequence[int]) -> tuple[int, ...]:
    if not isinstance(values, (list, tuple)) or not 1 <= len(values) <= MAX_THRESHOLDS:
        raise ContractError(f"supply 1..{MAX_THRESHOLDS} positive depth thresholds")
    thresholds = [exact_uint64(value, "depth threshold") for value in values]
    if 0 in thresholds:
        raise ContractError("depth thresholds must be positive")
    return tuple(sorted(set(thresholds)))
