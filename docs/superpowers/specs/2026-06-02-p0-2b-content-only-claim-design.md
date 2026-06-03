# P0.2b — Path-normalized (content-only) claim hash

**Status:** design
**Date:** 2026-06-02
**Roadmap:** `docs/ROADMAP.md` §6 Phase 0, item P0.2b
**Predecessor:** P0.2 (claim/measurement split) — surfaced this gap in adversarial review

## The remaining defect

P0.2 removed the machine-dependent measured *cost* from the claim hash. But the claim still
serializes, for every input/output, `{path, blake3}` — and `path` is recorded verbatim from the
CLI (`FileHash.path`, fed unmodified at the four manifest sites in `src/main.rs`). Absolute home
dirs, per-invocation tmp dirs, relative-vs-absolute invocation — all flow straight into the hashed
claim. So **two machines with byte-identical data at different paths still produce different
`manifest_blake3`.** The claim hash is therefore still not a cross-machine content-address — the
exact property chaining / `reproduce` / cohort (Phase 3) depend on.

The content hashes (`blake3`) are already cross-machine stable; only the paths are not. They are
*incidental* to the computation: the same `subcommand` + `params` consuming the same input
*contents* and producing the same output *contents* is the same logical run, wherever the files
live or whatever they are named.

## Design: a content-only claim form (version-gated)

Make the **claim** canonical form (the bytes `content_hash()` hashes) render inputs/outputs as the
**sorted list of their `blake3` digests**, dropping `path`:

```
claim:  {"has_measurements"…,"inputs":["<blake3>",…],"outputs":["<blake3>",…],"params":{…},"subcommand":"…","tool_version":"…"}
```

The **on-disk** receipt (`to_canonical_json`, what `write_manifest` persists and `verify` re-hashes
files from) is **unchanged** — it keeps the full `{path, blake3}` per file, because humans need to
see which files and `verify` re-hashes each file *at its recorded path*. Only the claim form used
for hashing changes. Roles are preserved (inputs and outputs stay separate lists); duplicates are
kept (the sorted blake3 list is a multiset, preserving cardinality).

This makes `manifest_blake3` a genuine content-address: identical input/output contents + identical
params + subcommand + version ⇒ identical claim hash, regardless of path.

### What paths still protect (and why dropping them from the claim is safe)

A path edit in the receipt is still caught — by `verify`'s input/output re-hash loop, which reads
each file *at its recorded path*: a path pointed at a missing or different-content file fails the
re-hash. A path pointed at a same-content file passes both re-hash and claim — which is correct,
because same content = same computation. The claim simply stops committing to the machine-specific
*location*; the content commitment (the `blake3`) is unchanged. Paths become advisory metadata in
the receipt, integrity-checked by re-hash rather than by the portable claim hash.

### Version gating (backward compatibility)

Pre-P0.2b receipts (`schema_version` 1 or 2) recorded `manifest_blake3` over a claim that **included
paths**. Changing the claim form unconditionally would make those receipts fail self-verify. So the
claim form is selected by `schema_version`:

- `schema_version >= 3` → **content-only** claim (path-independent).
- `schema_version < 3` or absent → **path-inclusive** claim (the pre-P0.2b form) — reproduces the
  old hash exactly, so old receipts still `self_hash_ok() == Some(true)`.

Bump `MANIFEST_SCHEMA_VERSION` `2 → 3`. `finalize()` stamps `schema_version` *before* computing the
claim hash (unchanged ordering), so a freshly sealed receipt is schema 3 and commits to the
content-only claim. The version gate lives in one place (`claim_file_render`, consulted by
`content_hash` via `to_canonical_claim_json`); no other code branches on it.

## Implementation (`src/provenance/mod.rs`)

- `MANIFEST_SCHEMA_VERSION` `2 → 3`.
- A `FileRender` enum (`WithPath` | `ContentOnly`) threaded through the canonical serializer.
- `push_canonical(include_measurements, file_render)`: dispatch file arrays to `push_file_hashes`
  (WithPath — existing `[{"blake3","path"}]` sorted by path) or a new `push_blake3_list`
  (ContentOnly — `["<blake3>",…]` sorted by blake3).
- `to_canonical_json()` → `push_canonical(true, WithPath)` (on-disk, unchanged bytes).
- `to_canonical_claim_json()` → `push_canonical(false, self.claim_file_render())`.
- `claim_file_render()` → `ContentOnly` iff `schema_version >= 3`, else `WithPath`.
- `content_hash()` unchanged (still hashes `to_canonical_claim_json()` with `manifest_blake3`
  removed). The parser, `measurements`, `measurement_blake3`, `has_measurements`, and `get_recorded`
  are all unaffected.

`verify` is functionally unchanged: it re-hashes files from the on-disk full form (still has paths)
and checks `self_hash_ok()` (now content-only for schema-3 receipts) + `measurement_hash_ok()`.

## Testing

New `provenance` unit tests:
1. **Path-independence (the keystone)** — two finalized manifests with the *same* input/output
   `blake3` digests but *different* `path`s produce *equal* `content_hash()` and both
   `self_hash_ok() == Some(true)`. (This is the test the review noted would fail before P0.2b.)
2. **Content sensitivity** — changing a `blake3` (different content) *does* change the claim hash.
3. **Claim drops paths, on-disk keeps them** — `to_canonical_claim_json()` does not contain a
   recorded path string; `to_canonical_json()` does.
4. **Pre-P0.2b back-compat** — a hand-built `schema_version = "2"` receipt with paths reproduces its
   path-inclusive `manifest_blake3` (still `self_hash_ok() == Some(true)`), and a schema-3 receipt
   over the same files hashes *differently* (the form genuinely changed).

Integration: bump the `schema_version` assertion `"2" → "3"`
(`tests/plan_enforce.rs::receipt_is_self_hashing_and_schema_versioned`). All existing receipts are
freshly generated, so the rest of the suite exercises the schema-3 path end-to-end (verify still
passes on untampered runs; tamper/strip/consistency tests still fire).

## Out of scope

Signing/attestation (later); build-identity fields (P0.3). The content-only claim is the unit that
P0.3's `code_git_sha`/lockfile and Phase 3's `reproduce` will hash/chain/sign.

## Risks

- **A machine-dependent value other than path still in the claim** → re-breaks cross-machine
  stability. The P0.2 review enumerated the claim fields; only paths remained, and this closes them.
  The path-independence unit test is the guard.
- **Old-receipt regression** → covered by the schema-2 back-compat test (version gate).
