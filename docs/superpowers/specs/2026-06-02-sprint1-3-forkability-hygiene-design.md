# Sprint 1.3 — Forkability Hygiene (design)

**Status:** Approved design — 2026-06-02. Increment 1.3 of the engineering roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md), Sprint 1).

---

## 1. Problem

Rosalind's whole value proposition is forkability (~249 stars, 10+ forks), but a new builder hits
verified day-one friction: broken commands in the docs, an unenforced MSRV claim, and no contribution
scaffolding. Cheapest credibility-per-hour in the roadmap.

Verified defects (2026-06-02):
- `README.md:387` tells you to run `cargo test --test space_bounds` — that test was **deleted** in the
  dead-code prune. The very first thing a curious forker copy-pastes errors out.
- `python/README.md:19` prints `ft.array.shape`, but the `FeatureTable` field is **`data`** (`ft.data`)
  — the documented Python example throws `AttributeError`.
- `.github/workflows/ci.yml:55,96` key the pip cache on `hashFiles('pyproject.toml')`, but that file
  was **removed** — `hashFiles` of a missing path is empty, so the key never varies (a silently broken,
  always-cold cache).
- `Cargo.toml` has **no `rust-version`** — the MSRV is asserted in prose but not declared or enforced.
- No `CONTRIBUTING.md`, `ARCHITECTURE.md`, or issue template — a contributor has no paved road.

## 2. Goals / non-goals

**Goals**
- Fix the three verified doc-drift defects.
- Declare + **enforce** the MSRV (`rust-version` + a CI job that builds on it).
- Add `CONTRIBUTING.md`, `ARCHITECTURE.md`, and a deterministic-repro issue template.
- Confirm the crate **packages** cleanly for crates.io (`cargo publish --dry-run`).

**Non-goals (explicitly deferred — see §6)**
- The `clippy -D warnings` CI gate (82 pre-existing lints + an MSRV/`div_ceil` tension → its own
  "1.3b clippy cleanup" increment).
- Actually running `cargo publish` (irreversible; a gated step for the maintainer).
- PyPI / bioconda / Docker channels (the htslib-linking wall — a separate, later effort).
- Running the README quickstart verbatim in CI (`release.yml` already smoke-tests `install.sh` on tag).

## 3. Design

### 3.1 MSRV — declare + enforce

- `Cargo.toml`: add `rust-version = "1.72"` to `[package]` (the documented MSRV).
- `.github/workflows/ci.yml`: add an `msrv` job that installs Rust **1.72** and runs `cargo check
  --all-targets`. This catches any accidental use of a newer-than-1.72 API (the codebase already
  avoids `u64::div_ceil`, a 1.73 API, via manual arithmetic — consistent with 1.72).
- Pre-push local check: `rustup install 1.72.0 && cargo +1.72.0 check --all-targets`. If it does not
  compile on 1.72, **bump `rust-version` to the true minimum** rather than weaken the gate.

### 3.2 Doc-drift fixes

| File:line | Now | Fix |
|---|---|---|
| `README.md:387` | `cargo test --test space_bounds … streaming evaluator` | `cargo test --test plan_enforce      # the memory contract: plan / --enforce / verify / governor` |
| `python/README.md:19` | `print(ft.array.shape)` | `print(ft.data.shape)` |
| `ci.yml:55,96` | `hashFiles('pyproject.toml')` | `hashFiles('.github/workflows/ci.yml')` (busts the cache when the install step changes; an existing path) |

The historical reference in `docs/superpowers/plans/2026-06-01-phase-c1-…md` is an archived plan — left
as-is (it records what was true then).

### 3.3 Forkability scaffolding

- **`CONTRIBUTING.md`** — build/test/fmt commands; the project's `brainstorm → spec → plan → TDD`
  increment flow (with pointers to `docs/superpowers/`); the **two load-bearing invariants** a change
  must preserve (byte-identical determinism per `docs/determinism.md`; the memory contract per
  `CONTRACT.md` — `predicted ≥ realized`, the governor, the self-hashing receipt); the "implement one
  `ColumnAnalyzer`" extension path; branch/PR conventions (a branch per increment, PR for review).
- **`ARCHITECTURE.md`** — a one-screen map of `src/` modules to responsibilities (`core`, `io`,
  `genomics`, `pileup`, `call`, `provenance`, `util`) and the roadmap phases (A–E) to where they live,
  so a contributor knows where a change belongs. Links `docs/ROADMAP.md` + `docs/OPEN_PROBLEMS.md`.
- **`.github/ISSUE_TEMPLATE/bug_report.md`** — a deterministic-repro template: exact command + inputs,
  expected vs actual, the **receipt hash** (`manifest_blake3`) + `rosalind --version`, OS/arch. Steers
  reports toward the reproducible artifacts the engine already produces.

### 3.4 Publish readiness

Run `cargo publish --dry-run` and fix anything that blocks packaging (e.g. excluded files, metadata).
**Do not publish.** The crate metadata is already complete (`description`, `license`, `keywords`,
`categories`, `repository` in `Cargo.toml`); this step confirms the tarball builds and is a no-op if it
already passes. The real `cargo publish` stays a maintainer-gated, irreversible action.

## 4. File-by-file change list

| File | Change |
|---|---|
| `Cargo.toml` | `rust-version = "1.72"`. |
| `.github/workflows/ci.yml` | New `msrv` job (Rust 1.72 `cargo check`); repoint the 2 pip-cache keys off `pyproject.toml`. |
| `README.md` | Replace the `space_bounds` test line. |
| `python/README.md` | `ft.array` → `ft.data`. |
| `CONTRIBUTING.md` (new) | Contribution guide. |
| `ARCHITECTURE.md` (new) | Module + phase map. |
| `.github/ISSUE_TEMPLATE/bug_report.md` (new) | Deterministic-repro bug template. |

## 5. Testing / verification

- `cargo build --release` + `cargo test` stay green; `cargo fmt --check` clean (no code changes, but
  `Cargo.toml` edits must not break the build).
- `cargo +1.72.0 check --all-targets` passes locally (or `rust-version` is bumped to the true min).
- `cargo publish --dry-run` succeeds.
- The fixed doc commands actually run: `cargo test --test plan_enforce` passes; the python snippet
  field (`ft.data`) matches `python/rosalind.py`.
- CI green on GitHub (incl. the new `msrv` job).

## 6. Deferred follow-up (noted, not in this increment)

**1.3b — clippy cleanup + `-D warnings` gate.** Clear all 82 pre-existing `clippy --all-targets` lints
(needless-borrow ×17, `manual_div_ceil` ×10, manual-`for` loops, `repeat().take()`, etc.) and add a
CI clippy gate. The `manual_div_ceil` lints force an MSRV decision: either bump MSRV to **1.73** to use
`u64::div_ceil` (and replace the manual arithmetic, including in `tests/plan_enforce.rs`), or
`#[allow(clippy::manual_div_ceil)]` to stay at 1.72. Resolve there, not here.

## 7. References

- `README.md:387`, `python/README.md:19`, `.github/workflows/ci.yml:55,96` — the verified defects.
- `Cargo.toml` — metadata (publish-ready) + the missing `rust-version`.
- `docs/ROADMAP.md` Sprint 1 (QW-1) — this increment.
