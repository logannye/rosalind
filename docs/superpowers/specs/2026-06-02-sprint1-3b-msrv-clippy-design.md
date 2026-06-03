# Sprint 1.3b — MSRV + Clippy Gate (design)

**Status:** Approved design — 2026-06-02. Increment 1.3b of the engineering roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md), Sprint 1). Picks up the MSRV + clippy items deferred from 1.3.

---

## 1. Problem

Two forkability defects deferred from 1.3, resolved together because they are one decision:

1. **No enforced MSRV.** `Cargo.toml` has no `rust-version`; the prose "1.72" claim is false (the
   committed `Cargo.lock` is format v4 and the dep tree — `url`→`rust-htslib` ICU chain, `proptest` —
   requires rustc **1.83**). A forker on an older toolchain gets an inscrutable failure.
2. **82 pre-existing `clippy --all-targets` lints** and no clippy CI gate, so the backlog grows
   unchecked and a contributor's PR can add more.

## 2. Goals / non-goals

**Goals**
- Declare the **real, verified MSRV (1.83)** and enforce it with a CI job.
- Clear **all 82** `clippy --all-targets` lints — behavior-preserving.
- Add a `clippy … -D warnings` CI gate so the backlog cannot regrow.

**Non-goals**
- Any behavior change. Every lint fix is a semantic no-op; the full test suite staying byte-green is
  the guard (a clippy fix that changes a test result is reverted, not forced).
- Pinning deps down for a lower MSRV (decided: accept 1.83 — zero-maintenance vs a fragile
  `cargo update --precise` treadmill for an audience that runs recent Rust).
- Bumping the crate version or touching the public API.

## 3. Design

### 3.1 MSRV = 1.83

- `Cargo.toml`: `rust-version = "1.83"` in `[package]`. **Verified**: `cargo +1.83.0 check
  --all-targets` builds the lib, bin, and tests against the committed v4 lock.
- `.github/workflows/ci.yml`: a new `msrv` job — `actions-rs/toolchain@v1` pinned to `1.83.0`, then
  `cargo check --all-targets` (matches the existing jobs' toolchain action). Catches any future use of
  a post-1.83 API.

### 3.2 Clearing the 82 clippy lints

Resolved by category (counts from `cargo clippy --all-targets`):

| Lint | × | Resolution |
|---|---|---|
| `manual_div_ceil` | 10 | Replace `(x + n - 1) / n` with `x.div_ceil(n)` — now available at MSRV 1.83 (was the blocker at 1.72). |
| `needless_borrow` | 17 | Drop the redundant `&` (`cargo clippy --fix`). |
| `needless_range_loop` (manual `for`) | 6 | Iterate directly / `enumerate` (`--fix` where safe, else manual). |
| `manual_repeat_n` / `repeat().take()` / manual `str::repeat` | 8 | `std::iter::repeat_n` / `"x".repeat(n)`. |
| `useless_vec` | 2 | `[…]` array instead of `vec![…]`. |
| `bool_assert_comparison` | 2 | `assert!(x)` instead of `assert_eq!(x, true)`. |
| `unnecessary_lazy_evaluations` (`unwrap_or_else(\|\| None-ish)`) | 2 | `unwrap_or(…)`. |
| `manual_is_multiple_of` | 2 | `x.is_multiple_of(n)` (available at 1.83). |
| `ptr_arg` (`&PathBuf`→`&Path`) | 1 | Take `&Path`. |
| `unused_enumerate_index` | 1 | Drop `.enumerate()`. |
| `manual_contains` (`iter().any` vs `contains`) | 1 | `contains`. |
| `too_many_arguments` (8/7, 11/7) | 2 | `#[allow(clippy::too_many_arguments)]` + rationale on the two **streaming-driver** fns (a bounded source + region + params + sinks legitimately needs them; already the pattern in `gvcf.rs`). |
| `type_complexity` | 1 | A `type` alias for the complex `Result<…>`, or a scoped `#[allow]` if the alias hurts readability. |

**Method:** run `cargo clippy --fix --all-targets --allow-dirty` for the mechanically-fixable lints,
then hand-resolve the remainder (the `div_ceil`/`is_multiple_of` rewrites, the two `#[allow]`s, the type
alias). After each pass: `cargo test` (behavior unchanged) + `cargo clippy --all-targets` (0 warnings) +
`cargo fmt`. (Run inline — `cargo clippy --fix` mutating a shared tree is fine for the maintainer, but
must never be delegated to a subagent per the shared-tree hazard.)

### 3.3 The `-D warnings` CI gate

Add a `clippy` job to `ci.yml`: `actions-rs/toolchain@v1` (stable, `components: clippy`), then
`cargo clippy --all-targets --all-features -- -D warnings`. Once the 82 are cleared this stays green and
fails any PR that introduces a new lint.

## 4. File-by-file change list

| File | Change |
|---|---|
| `Cargo.toml` | `rust-version = "1.83"`. |
| `.github/workflows/ci.yml` | New `msrv` job (1.83 `cargo check`); new `clippy` job (`-D warnings`). |
| `src/**`, `tests/**`, `examples/**` | The behavior-preserving lint fixes (many files, small mechanical edits; the test suite is the guard). |

## 5. Testing / verification

- `cargo clippy --all-targets -- -D warnings` → **clean** (the exit criterion).
- `cargo test` → all sections green (behavior preserved by every fix).
- `cargo build --release` (+ examples) → 0 warnings; `cargo fmt --check` clean.
- `cargo +1.83.0 check --all-targets` → builds (MSRV honored).
- CI green on GitHub incl. the new `msrv` + `clippy` jobs.

## 6. References

- The 1.3 spec §3.1/§6 — where MSRV + clippy were deferred, with the dep-floor finding.
- `cargo clippy --all-targets` (2026-06-02) — the 82-lint inventory in §3.2.
- `cargo +1.83.0 check --all-targets` — the verified MSRV.
