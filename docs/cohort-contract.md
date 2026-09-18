# Proposed contract: local cohort candidate reanalysis

**Internal next-minor foundation, not a shipped feature.** This isolated
development branch contains crate-private snapshot import/verification and
comparison contracts. It has no `cohort` command or `open_cohort` Python API.
The implementation stays out of the 0.5 release candidate and its frozen public
contract. Engineering can progress while publication and adoption work continue.
Independent participants and design partners inform usability and product decisions;
their completion is not a software shipping prerequisite. No partner validation,
representative cohort performance or clinical interpretation is claimed.

The first outcome is repeat candidate-SNV questions across explicitly identified
specimens, using saved evidence. Paired/longitudinal analysis is a later feature;
optional subject/timepoint metadata does not imply pairing or longitudinal inference.
The [synthetic fixture](../examples/cohort-reanalysis/README.md) exercises current
single-sample primitives and independently specified cohort expectations.

## Existing implementation seams

| Existing interface | Proposed use |
|---|---|
| `VerifiedEvidenceDataset::open`, `descriptor`, `fields` | Bounded verified metadata and physical capabilities |
| `coverage_split` | Determine represented versus absent loci before decoding |
| `plan`, `visit_batches` | Project and stream covered evidence with native budget checks |
| `source_hashes`, `verified_partition_hashes`, `verify_unchanged` | Bind verified lineage and guard immutable inputs |
| `DatasetQuery`, `EvidenceSelection`, `EvidenceFields` | Reuse selection normalization and capability definitions |
| `VerifiedInputSession::compatibility_key` | Verify unchanged original inputs for same-specimen extension |
| `run_dataset_with_snapshot`, `publish_evidence_dataset` | Publish intact missing-only leaf datasets |
| `plan_reuse`, `run_reusing_dataset` | Existing single-source reuse model and fresh/reused equality oracle |
| `run_evidence_artifact` lifecycle | Admission, cancellation, staged output, receipts and explicit replay |
| Python `EvidenceDataset` / `open_dataset` | Lazy native-process adapter conventions |

These interfaces live in `src/dataset/persisted.rs`, `descriptor.rs`, `reuse.rs`,
`src/evidence/artifact.rs` and `python/rosalind/dataset.py`. The lower-level readers
can feed a cohort source adapter. The current artifact runner handles one source;
its process-wide scope must not be recursively started for every cohort member.
Internal lifecycle reuse must preserve existing single-source APIs and goldens.

## Members and scientific identity

A member is one asserted specimen ID bound to one named analysis sample. IDs must
be unique, nonempty bounded strings with a stable lexical order; reject control
characters and invalid metadata before output creation. Optional group, subject
and timepoint fields are bounded user assertions, carried in snapshot identity.
They are not proof of identity, independence or biological ordering.

Named `@RG SM` scope is required for v1 candidate summaries. Preserve the entire
resolved scope on each leaf; unknown or explicitly pooled scope cannot silently
become an individual. Different members normally have different sample names/read
groups. Do not require equal scope JSON across different specimens. Reject an
identical source/scope entered twice under aliases to prevent accidental double
counting. Distinct BAMs with the same SM are not automatically the same specimen;
the author must supply specimen IDs.

Within a member, every leaf must have the same existing source-reuse compatibility
key. Physical field masks may differ across nonoverlapping leaves; a future query
must find all its requested groups in every consumed leaf. This keeps extensions
tied to unchanged original inputs without equating physical storage with science. Each selected coordinate
belongs to at most one leaf for that member. Reject overlap even when values agree;
never add duplicate counts or apply implicit last-writer precedence.

### Separate comparison compatibility from reuse identity

The existing `compatibility_blake3` includes alignment/index bytes and exact sample
scope. Keep it unchanged for reuse. A new cross-member comparison contract requires:

- the same evidence semantics, read counting unit and exact/no-sampling policy;
- equal explicit mapping/base quality thresholds and flag exclusions;
- equal ordered contig names and lengths, with actual reference sequence available;
- matching effective analysis-reference content identity and representation;
- supported schema/field-mask versions, decoded through their existing validators;
- every selected leaf supplies all requested groups and reducer-required groups.

The effective analysis reference is the descriptor's `reference` source, or its
`cram-reference` source when that supplies the analysis bases. V1 conservatively
requires matching byte identity and source length. Normalize the effective source
role: a BAM using `reference` and a CRAM using `cram-reference` can compare when
the actual reference bytes match. When an effective FAI is present, compare its
byte hash and length too, normalizing `reference-fai`/`cram-reference-fai`. Equal
FASTA bytes and dictionary alone do not establish identical offset mappings.
Paths and filename extensions are provenance,
not comparison keys. Equal dictionary
or assembly label alone is insufficient; a FASTA and an equivalent `.rref` are
not silently equated. A future sequence-based reference identity needs its own
specification and migration evidence.

Alignment/index content and named-sample labels can differ across members.
Physical field masks can differ across members if all cover the requested
projection. Scientific mismatches, missing required groups and unsupported context
are preflight errors, not per-row missingness. V1 supports existing SNV selections;
indels, symbolic alleles, reference-free candidate comparison and filter changes
remain unsupported. A new ALT at a stored locus can use existing A/C/G/T counts.

## Snapshot storage and import

Use one local directory with shared immutable objects:

```text
cohort-root/
  objects/<leaf-manifest-byte-hash>/  # complete original portable dataset layout
  snapshots/<snapshot-content-hash>/
    snapshot.json                   # canonical proposed snapshot metadata
    manifest.json                   # receipt binding snapshot bytes and inputs
```

The internal snapshot schema is version 1, independently versioned from all
existing dataset formats. The leaf object key is lowercase BLAKE3 of the complete
`evidence-dataset.manifest.json` bytes, not its receipt self-hash and not the legacy
`dataset.manifest.json`. The snapshot ID is lowercase BLAKE3 of canonical compact
UTF-8 `snapshot.json`; preserve struct field order and sorted member/leaf arrays.
Readers reject noncanonical encodings, unknown keys and unsupported versions. The snapshot descriptor
pins member metadata, ordered leaf identities, their descriptor hashes, comparison
contract version and optional parent-snapshot identity. Do not include its own
hash inside the bytes being hashed. Paths stored in a snapshot are relative to
`cohort-root`; reject absolute paths, traversal and symlinks inside the store.

`create` receives explicit local dataset locations, copies each declared dataset
inventory into staging with ordinary file copies, verifies all copied metadata,
receipts, ownership and payload bytes, then publishes each object atomically.
Import identical objects once. Do not use hardlinks/symlinks: mutable external
source files must not change the stored object. Existing objects must verify before
reuse. Original dataset wire schemas, partition receipts and producer identities
remain intact; no old partition is relabeled under a new request namespace.

Publish the snapshot directory last with a same-filesystem atomic create-new
directory rename (`RENAME_NOREPLACE` on Linux, `RENAME_EXCL` on macOS). No
replacement fallback is allowed. The directory contains the descriptor and a
self-checking receipt binding its identity, leaf manifests/descriptors and optional
parent descriptor. A failure cannot expose a completed new snapshot. Existing snapshots remain unchanged. Interrupted imports may leave
verified unreferenced objects; automatic garbage collection is out of v1 scope.
Copying the entire root preserves offline access to every referenced snapshot.
Report one-time copying time and additional storage as part of reuse cost.


The copied portable inventory is exactly the portable manifest, descriptor, and
all declared partition `manifest.json`/`evidence.arrow` pairs. Preserve these bytes
without modification. The source directory's legacy cache manifest and unrelated
files are not part of this content-bound portable inventory and are not copied.
Original alignment/reference paths in receipts remain provenance claims; import
and verification never follow them. Existing objects reject undeclared files or
directories, including empty directories. Source inventories reject symlink files
and directory components; ordinary copying creates independent inodes.

Opening a snapshot verifies its metadata, parent chain, leaf metadata and member
ownership. Full verification additionally streams every unique current-snapshot
leaf through the existing bounded Arrow reader, checking every partition hash,
receipt, row count and coordinate. Ancestor descriptors/receipts are checked, but
unreferenced ancestor payloads are not claimed to be consumed. Imported objects
always receive full verification. The store may contain scientifically incompatible
members so they can be inspected and explicitly selected; query planning owns
cross-member compatibility and named-sample requirements.

Snapshot shape (hash values below are explanatory placeholders):

```json
{"version":1,"comparison_version":1,"parent":null,"members":[{"metadata":{"id":"sample-A","group":null,"subject":null,"timepoint":null},"leaves":[{"object_id":"64 lowercase hex characters","descriptor_blake3":"64 lowercase hex characters"}]}]}
```

Members sort by exact UTF-8 ID bytes; leaves sort by object ID. IDs are nonempty,
at most 256 UTF-8 bytes; optional metadata is nonempty when present and at most
1,024 bytes per value. Control characters are rejected and no Unicode normalization
is implicit. Empty cohorts are representable; an included member requires at least
one leaf. Duplicate member IDs, repeated leaf references, overlapping selected
intervals and source/scope aliases are rejected. Disjoint intervals inside the
same canonical partition are valid. Snapshot identity includes asserted metadata
and parent identity, and contains no runtime timestamp or absolute local path.

Default reader/import envelopes are 8 MiB snapshot JSON, 65,536 members, 262,144
leaf references, 262,144 ownership intervals per member, and 1,024 ancestor
snapshots. Snapshot receipts have a 32 MiB envelope. Existing dataset envelopes
remain 32 MiB manifest/descriptor and 64 KiB partition receipt. These are adjustable
operational refusal limits, not claimed scientific scalability. Parsing/encoding,
copy buffers, inventory clones, ownership state and retained mutation guards are
admitted before their associated allocations; leaf decoders run serially. Current
admission uses conservative measured process high-water RSS plus transient
reservations, and can overestimate after earlier allocations are released. It is
not a tight aggregate cohort query planner or an OS allocation cap.

The internal functions obey existing cancellation/governor checkpoints and retain
ordinary-mutation guards through publication. They do not start a nested governor.
The complete outer lifecycle and cross-budget execution study remain C07 work.
File synchronization plus atomic directory visibility does not claim recovery
from power loss. Readers require an immutable store during an operation; integrity
checks do not provide authenticity or defense against deliberate concurrent
metadata restoration. Garbage collection, replacing an object, and mutable
`latest` aliases remain unsupported.

## Query, missingness and candidate summary

Every request identifies an immutable snapshot, member selection, supported SNV
selection and physical fields. Normalize the loci once against the shared
reference dictionary. Plans report each member's covered/missing loci, capability
and compatibility failures, output cardinality and resource reservations. Planning
verifies bounded metadata; execution verifies every consumed partition. Do not
claim untouched partition payloads were checked by an ordinary query.

Default policy is **strict**: any requested member/locus absent from saved evidence
refuses the materialized query before output creation. `partial` is explicit and
emits every requested member/candidate row with an observation status:

| Situation | Counts | Fraction | Eligibility |
|---|---|---|---|
| Not stored (`unmeasured`) | null | null | null |
| Stored, callable depth zero | exact zeros | null | false for positive depth threshold |
| Stored, positive depth below threshold | exact counts | ALT/depth defined | false |
| Stored, depth at/above threshold | exact counts | ALT/depth defined | true |

Use `status=observed` or `status=unmeasured`; derive zero/low/eligible from depth and
recorded threshold rather than multiplying status categories. No field failure or
profile incompatibility is converted into a null row. A zero ALT count is not a
reference genotype and does not prove biological absence.

The first reducer is per locus/ALT, with configurable positive
`min_callable_depth` (default 10). ALT support is fixed at at least one eligible
ALT observation in v1; configurable support thresholds are deferred:

- `n_requested`: selected cohort members;
- `n_observed`: members with a stored row, including zero depth;
- `n_depth_eligible`: observed members at the depth threshold;
- `n_alt_supported`: depth-eligible members with at least one ALT observation;
- total callable and ALT observations across all observed members;
- callable and ALT totals restricted to depth-eligible members.

Report denominators separately; do not call these genotype frequencies or calibrated
confidence. Retain integer fraction numerator/denominator. At a zero denominator,
fraction is null. The initial interoperable output can expose those two integers
without committing to floating-point formatting. A later display adapter may
render a ratio with an explicit policy. Multiallelic rows repeat the per-locus
denominator while using the corresponding ALT count; they do not sum denominators.

Canonical long-form order is member ID, reference dictionary, position, ALT base.
Summary order is dictionary, position, ALT base. Use checked integer addition.
For missingness, call existing `visit_batches` only on the covered selection;
its strict single-dataset behavior must not be weakened.

## Bounded execution and artifact lifecycle

V1 is serial. Long-form extraction streams one member at a time. Summary reduction
uses resource-admitted windows within the existing 16,384-base ownership partitions,
visits members sequentially, and retains only bounded window counters and one
member decoder/projection. Ownership is fixed; computation width is not. A smaller
admitted window may require additional saved-partition decoding. Callback
batches remain at most 1,024 rows. Never materialize a dense members-by-genome
matrix. Include cohort/member metadata, source physical decoders, normalized
selections, output buffers and receipt finalization in admission. Memory is not
independent of metadata cardinality; oversized inventories may refuse explicitly.

Use one outer managed cancellation/governor scope. Native or consumer allocations
can precede checkpoints; a declared budget is not a universal OS allocation cap.
A completed primary output and receipt publish together. Refusal, ordinary failure,
corruption or cancellation cannot publish success; any resource partial is clearly
identified and never accepted as a completed snapshot or query.

## Explicit extension

Queries never discover or reopen original BAM/CRAM/reference files. A separate
`extend` request supplies an explicit member-to-local-source map and new loci.
For each affected member, use a verified input session to check exact original
source/reuse identity, then subtract existing coverage and extract only missing
loci. Preserve the member's profile and scope and explicitly specify the new leaf's
physical fields. Every later query still requires adequate fields at every stored
locus; extension never silently fills missing groups at existing loci. Publish missing-only
leaf objects and a new parent-linked snapshot last. Old snapshots and leaves do
not change. A fully covered no-op creates no new evidence leaf.

New-field backfill, replacement input BAMs, changed filters and overlapping leaves
are out of scope; make a new extraction/snapshot instead. CRAM extension still pays
its complete source-validation cost, and rehashing remains real work. Current
`run_reusing_dataset` can validate fresh-versus-reused equality but does not itself
publish expanded caches. Use the existing native partition publisher for each delta.

## Proposed CLI and Python surface

**The following names are proposals and must not be run against current releases.**

```text
rosalind cohort create --members members.tsv --output cohort-root
rosalind cohort inspect --cohort cohort-root --snapshot SNAPSHOT_ID
rosalind cohort verify --cohort cohort-root --snapshot SNAPSHOT_ID
rosalind cohort extract --cohort cohort-root --snapshot SNAPSHOT_ID \
  --sites candidates.vcf --fields depths,alleles --missing strict --plan
rosalind cohort summarize --cohort cohort-root --snapshot SNAPSHOT_ID \
  --sites candidates.vcf --missing partial --min-callable-depth 10
rosalind cohort extend --cohort cohort-root --snapshot SNAPSHOT_ID \
  --sites expanded.vcf --sources sources.tsv
```

Creation/extension return the immutable snapshot identifier. Require that identifier
for reads; no implicit mutable `latest` pointer. CLI output uses Arrow/TSV and
create-new behavior, with explicit replacement only for derived artifacts where
supported. Snapshot/object replacement is never implicit. `--plan` is a mode of
extract/summarize/extend, not a separate scientific operation.

Proposed Python follows the current lazy dataset adapter:
`open_cohort(root, snapshot=...)`, `.batches(...)`, `.materialize(...)`,
`.summarize(...)`. It delegates to the matching native executable, preserves error
codes and requires stream exhaustion for completed results. Retained Python arrays
are outside native budgeting. Keep a new generic cohort analyzer API internal
until first-party reducers establish the necessary contract.

## Verification, replay and acceptance

Derived receipts bind snapshot bytes, normalized member/query selection, consumed
leaf metadata/partitions, reducer identity/parameters and actual output bytes.
Verification identifies exactly which saved inputs were checked. Add tokenized
Arrow/TSV replay using the intact relocated cohort root and external selection
files. Explicit external binaries remain required. Parquet directory byte replay
is not added by this proposal. A receipt neither authenticates authorship nor
proves biological truth.

Before feature acceptance, require independent per-sample/candidate oracles,
strict/partial behavior, differing compatible masks, incompatible reference/profile
rejection, duplicate members/overlap rejection, uint64 overflow handling, and
admitted budget/tile invariance. Exercise corruption/mutation, bounded metadata,
source-free relocation, cancellation and publication failure. Extension must match
fresh extraction and leave old snapshot hashes unchanged. Existing single-source
schemas, goldens and conformance must remain valid.

The accompanying fixture checks **existing extraction/reuse primitives only**.
It cannot establish a future snapshot implementation, cohort runtime bound, real
assay performance, user demand or partner acceptance. Storage/comparison tests now exercise the internal C02/C03 foundation on this
branch. CLI/Python querying, managed reducers, extension and replay remain C04–C11
work; fixture arithmetic is not partner acceptance or a performance measurement.
