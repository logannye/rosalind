# P0.3 — Build-identity in the receipt

**Status:** design
**Date:** 2026-06-02
**Roadmap:** `docs/ROADMAP.md` §6 Phase 0, item P0.3 (closes Phase 0)
**Predecessor:** P0.2 / P0.2b (the claim is now a cross-machine content-address)

## The defect

The receipt records `tool_version = env!("CARGO_PKG_VERSION") = "0.1.0"` — a static crate version
that is identical across every commit. So "exactly this code" is a lie: two builds from different
commits (different behavior) both record `tool_version: "0.1.0"`. A reproducibility receipt whose
content-address (the P0.2b claim hash) does not commit to *which code, toolchain, and dependencies*
produced the run cannot anchor `reproduce` (Phase 3): "re-run with exactly this code → expect this
output" has no real key.

## Design

Bake five **build-identity** facts into the binary at compile time and stamp them into the receipt's
**claim** (so they are part of the cross-machine content-address and are tamper-protected by the
claim self-hash):

| Field | Source (build time) | Meaning |
|---|---|---|
| `code_git_sha` | `git rev-parse HEAD` | the commit the binary was built from |
| `code_dirty` | `git status --porcelain` non-empty | uncommitted changes at build time (not reproducible from a SHA alone) |
| `rustc_version` | `$RUSTC --version` | the compiler (affects codegen/behavior) |
| `target_triple` | cargo `TARGET` env | the build target |
| `deps_lock_blake3` | BLAKE3 of `Cargo.lock` | the exact resolved dependency versions |

These are deterministic *per binary* — the same binary reports the same values on any machine — so
they fit the P0.2b content-address: two runs of the same binary on the same data hash identically,
while two *different* builds (different commit/toolchain/target/deps) hash differently, which is
exactly the reproduction distinction `reproduce` needs.

### Capture: `build.rs`

A new `build.rs` runs each build, shells out to `git`/`rustc`, reads `Cargo.lock`, and emits the five
values as `cargo:rustc-env=ROSALIND_*` variables; the binary reads them with `env!`. Every fallible
step degrades to the literal `"unknown"` (non-git build, `git` absent, etc.) and the script never
panics, so `env!` always resolves and a build never fails because provenance was unavailable.

`Cargo.lock` is hashed with a `[build-dependencies] blake3` (matching the runtime version) for digest
consistency with the rest of the receipt. Re-run triggers: `Cargo.lock`, `.git/HEAD`, `.git/index`,
and the ref `HEAD` resolves to — so the SHA re-stamps on commit/stage/branch-switch. **Honest
limitation:** `code_dirty` is best-effort — accurate for CI/release builds from a committed tree; a
pure *unstaged* edit between incremental builds may not re-trigger the script, so the flag is
accurate at commit/stage granularity, not at the instant of an unstaged edit.

### Stamp: `finalize()`

`finalize()` stamps the five fields into `params` (the claim) **before** computing the claim hash, so
they are committed to and tamper-evident (editing `code_git_sha` breaks `manifest_blake3` — no
separate presence marker is needed, unlike the measurement block which lives *outside* the claim).
They are stamped in `finalize`, not `new`, so an un-sealed manifest (and the exact-JSON unit test)
stays free of the volatile SHA. Bump `MANIFEST_SCHEMA_VERSION` `3 → 4`: a v4 receipt is guaranteed to
carry build-identity. The `claim_file_render` gate stays `>= 3 → content-only` (v4 ≥ 3), unchanged.

### Check: `verify --expect-code <sha>`

A new `RunManifest::check_expected_code(expected) -> Vec<String>` (the testable core) returns
problems:

- no `code_git_sha` recorded (pre-P0.3 receipt) → cannot check;
- `code_git_sha == "unknown"` (non-git build) → cannot check;
- recorded SHA does not start with `expected` (prefix match, so short SHAs work) → `code mismatch`;
- match **but** `code_dirty == "true"` → matched, but built from a dirty tree → not reproducible
  from a commit SHA alone.

`verify` gains `--expect-code <sha>`; when supplied it extends `problems` with
`check_expected_code` and, on a clean match, prints a positive note. `verify` without the flag is
unchanged. The pure function keeps the match/dirty/absent logic unit-testable without a dirty working
tree (which a CLI e2e cannot reliably produce).

## Testing

Unit (`provenance`):
- `finalize` stamps all five build-identity fields (non-empty) and they are covered by the claim
  hash (editing `code_git_sha` after finalize → `self_hash_ok() == Some(false)`).
- `check_expected_code`: exact match (clean) → no problems; prefix match → no problems; mismatch →
  one problem; recorded dirty + match → one (dirty) problem; absent → one problem; `"unknown"` → one
  problem.

Integration (`tests/`):
- a real receipt's JSON contains `code_git_sha`, `rustc_version`, `target_triple`,
  `deps_lock_blake3`, and `schema_version` `"4"`.
- `verify --expect-code <a-wrong-sha>` exits 5 with `code mismatch` (robust regardless of the build's
  dirty state — the positive/dirty paths are covered by the `check_expected_code` unit tests).

## Out of scope

Signing/attestation of the receipt (Ed25519 `verify-attest`) — Phase 3. Build-identity is the key
that signing and `reproduce` will commit to.

## Risks

- **`env!` fails to resolve** if `build.rs` doesn't emit a var → mitigated by always emitting all five
  (fallback `"unknown"`), so the crate always compiles.
- **Stale `code_dirty`** on incremental local builds → documented as best-effort; the
  reproduction-critical CI/release builds (clean, freshly built) are accurate.
- **`git`/`Cargo.lock` absent** (vendored/published build) → graceful `"unknown"`, no build failure.
