# RFC 0001: receipt schema 6 (proposal only)

Status: open for community feedback. No v0.4 implementation is implied.

## Motivation

Schema 5 deliberately keeps a small, stable envelope and expresses new facts as
additive scalar parameters. That remains compatible, but it cannot natively express
typed artifacts, multiple digest algorithms, format-aware semantic identity, or
explicit provenance edges without conventions.

## Proposed direction

- Replace parallel path/hash records with typed artifact objects carrying a stable
  role, media type, size, and one or more digests.
- Require BLAKE3 and SHA-256 so native content addressing and standards ecosystems
  can share the same artifact identity.
- Add optional semantic digests for formats whose byte encoding can vary while
  preserving records, initially BAM and VCF. Semantic equality must never be
  silently substituted for byte equality.
- Represent parent/child DAG edges explicitly by claim and artifact role rather than
  inferring them from flags or path names.
- Define canonical upgrade and downgrade behavior, with unknown fields preserved.

## Compatibility questions

The RFC must settle canonical hashing, path portability, multi-output replay,
digest-agility downgrade attacks, semantic-digest versioning, and whether
reproduction certificates remain ordinary receipts. Canonical schema 1–5 fixtures
will remain readable. v0.4's in-toto bridge is an export, not a backdoor schema bump.

Community feedback should include real receipts from at least three external
analyzers before the proposal advances.
