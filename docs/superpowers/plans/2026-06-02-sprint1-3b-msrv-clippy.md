# Sprint 1.3b — MSRV + Clippy Gate: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Run INLINE — never delegate `cargo clippy --fix` / mutating git to a subagent (shared-tree hazard).

**Goal:** Declare + enforce the real MSRV (1.83) and clear all 82 `clippy --all-targets` lints behind a `-D warnings` CI gate — behavior-preserving.

**Architecture:** Mechanical lint cleanup (tool-assisted) + two CI jobs. The full test suite staying green is the proof that every fix is a semantic no-op.

**Tech Stack:** `cargo clippy`/`clippy --fix`, Cargo metadata, GitHub Actions YAML, Rust 1.83 (`u64::div_ceil` / `is_multiple_of` now available).

**Spec:** `docs/superpowers/specs/2026-06-02-sprint1-3b-msrv-clippy-design.md`

---

### Task 1: Clear the 82 clippy lints (behavior-preserving)

**Files:** `src/**`, `tests/**`, `examples/**` (small mechanical edits across many files)

- [ ] **Step 1: Snapshot the baseline**

Run: `cargo clippy --all-targets 2>&1 | grep -cE "^warning:"`
Expected: `82` (the baseline to drive to 0).

- [ ] **Step 2: Auto-fix the mechanically-fixable lints**

Run: `cargo clippy --fix --all-targets --allow-dirty --allow-staged 2>&1 | tail -5`
This resolves the machine-applicable lints (needless_borrow, useless_vec, bool_assert_comparison,
unnecessary_lazy_evaluations, ptr_arg, unused_enumerate_index, manual_contains, many of the loop/repeat
ones, and — at MSRV 1.83 — `manual_div_ceil` and `manual_is_multiple_of`).

- [ ] **Step 3: Verify behavior is unchanged after auto-fix**

Run: `cargo test 2>&1 | grep -iE "FAILED" || echo "no failures"`
Expected: `no failures`. **If any test fails, the auto-fix changed behavior — `git diff` the offending
file, revert that specific hunk, and resolve that lint by hand instead.**

- [ ] **Step 4: List the remaining lints**

Run: `cargo clippy --all-targets 2>&1 | grep -E "^warning:" | sort | uniq -c | sort -rn`
These are the ones `--fix` could not apply automatically — expected to be the two `too_many_arguments`
and the one `type_complexity`, plus any loop/repeat rewrites clippy declined.

- [ ] **Step 5: Resolve the `too_many_arguments` lints with a scoped allow**

For each function clippy flags (streaming drivers — a bounded source + region + params + sinks
legitimately need the arity), add directly above the `fn`:

```rust
#[allow(clippy::too_many_arguments)] // a bounded streaming driver: source + region + params + sinks
```

(`gvcf.rs::stream_gvcf_region` already uses this pattern — match it.)

- [ ] **Step 6: Resolve `type_complexity`**

For the flagged complex type (a `Result<(WorkingSet, SkipCounts), CoreError>`-style return), prefer a
local `type` alias near the function if it reads well, e.g.:

```rust
/// (max working set, skip counts) — the bounded-run telemetry a driver returns.
type DriveOutcome = (crate::core::WorkingSet, crate::pileup::SkipCounts);
```

…and use it in the signature. If the alias hurts readability at the single call site, instead add a
scoped `#[allow(clippy::type_complexity)]` above the function with a one-line rationale. Pick one.

- [ ] **Step 7: Resolve any remaining hand-fixable lints**

For each still-listed lint, apply the idiomatic fix clippy suggests (it prints the exact `help:`
rewrite). Common remainders: `needless_range_loop` → iterate the slice / `enumerate`; `manual_repeat_n`
→ `std::iter::repeat_n(x, n)`; manual `str::repeat` → `"x".repeat(n)`. Keep each edit behavior-identical.

- [ ] **Step 8: Drive to zero + verify green**

Run: `cargo clippy --all-targets 2>&1 | grep -cE "^warning:"` → `0`.
Run: `cargo test 2>&1 | grep -iE "FAILED" || echo "no failures"` → `no failures`.
Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"` → `no warnings`.
Run: `cargo fmt && cargo fmt --check && echo "fmt clean"` → `fmt clean`.

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "style(clippy): clear all 82 clippy lints (behavior-preserving; div_ceil/is_multiple_of via MSRV 1.83)"
```

---

### Task 2: Declare the MSRV + add the CI gates

**Files:** `Cargo.toml`, `.github/workflows/ci.yml`

- [ ] **Step 1: Declare the MSRV**

In `Cargo.toml` `[package]`, add after `edition = "2021"`:

```toml
rust-version = "1.83"
```

- [ ] **Step 2: Verify the crate builds on the declared MSRV**

Run: `cargo +1.83.0 check --all-targets 2>&1 | tail -3`
Expected: `Finished` (no errors). If a clippy fix in Task 1 used a >1.83 API, this catches it — fix to
a 1.83-compatible form (do NOT bump the MSRV to paper over it).

- [ ] **Step 3: Add the `msrv` + `clippy` CI jobs**

In `.github/workflows/ci.yml`, append under `jobs:` (match the existing jobs' `actions-rs/toolchain@v1`
structure):

```yaml
  msrv:
    name: MSRV (1.83)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions-rs/toolchain@v1
        with:
          toolchain: 1.83.0
          profile: minimal
          override: true
      - name: cargo check on the declared MSRV
        run: cargo check --all-targets

  clippy:
    name: Clippy (-D warnings)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions-rs/toolchain@v1
        with:
          toolchain: stable
          profile: minimal
          components: clippy
          override: true
      - name: clippy as a gate
        run: cargo clippy --all-targets -- -D warnings
```

- [ ] **Step 4: Verify the gate command passes locally + YAML parses**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: `Finished` (the gate passes — 0 warnings became 0 errors).
Run: `python3 -c "import yaml,sys; yaml.safe_load(open('.github/workflows/ci.yml')); print('yaml ok')"`
Expected: `yaml ok`.

- [ ] **Step 5: Commit**

```bash
git add Cargo.toml Cargo.lock .github/workflows/ci.yml
git commit -m "ci: declare MSRV 1.83 + enforce it and a clippy -D warnings gate"
```

---

## Final verification

- [ ] `cargo clippy --all-targets -- -D warnings` → clean.
- [ ] `cargo test` → all sections green; `cargo build --release` (+ `--examples`) → 0 warnings; `cargo fmt --check` clean.
- [ ] `cargo +1.83.0 check --all-targets` → builds.
- [ ] CI green on GitHub incl. the new `msrv` + `clippy` jobs.
- [ ] Public API unchanged (the lint fixes are internal; confirm no exported signature changed — `git diff main -- src/lib.rs src/call/mod.rs` is empty or only internal).
