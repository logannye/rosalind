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
- Add `CONTRIBUTING.md`, `ARCHITECTURE.md`, and a deterministic-repro issue template.
- Confirm the crate **packages** cleanly for crates.io (`cargo publish --dry-run`).

**Non-goals (explicitly deferred — see §6)**
- The MSRV declaration + CI gate — moved to **1.3b** after the execution finding in §3.1 (the documented
  1.72 is not buildable; the dep ecosystem forces ≥1.83). It is the same decision as `div_ceil`/clippy.
- The `clippy -D warnings` CI gate (82 pre-existing lints + the MSRV/`div_ceil` tension → **1.3b**).
- Actually running `cargo publish` (irreversible; a gated step for the maintainer).
- PyPI / bioconda / Docker channels (the htslib-linking wall — a separate, later effort).
- Running the README quickstart verbatim in CI (`release.yml` already smoke-tests `install.sh` on tag).

## 3. Design

### 3.1 MSRV — DEFERRED to 1.3b (finding during execution, 2026-06-02)

**Original plan:** declare `rust-version = "1.72"` + a CI job enforcing it.

**Finding:** the documented "1.72" MSRV is stale and not buildable. The committed `Cargo.lock` is
format **v4** (Rust ≥1.78 to parse), and its pinned transitive deps demand an even higher floor that
*climbs* as you test (`proptest 1.9` + its ICU chain → 1.82; `icu_properties_data 2.1.1` → 1.83; …).
Supporting a low MSRV would require pinning many transitive deps via `cargo update --precise` — a
fragile, high-maintenance choice that is a *deliberate dependency policy*, not hygiene. It is also the
**same decision** as the deferred `div_ceil`/clippy item (`div_ceil` needs 1.73). So MSRV moves to the
**1.3b** increment (§6), where the dep-pinning-vs-high-MSRV tradeoff is decided on purpose. No
`rust-version` or MSRV CI job ships in 1.3.

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

Run `cargo publish --dry-run` to confirm the crate packages cleanly. **Do not publish.**

**Result (2026-06-02):** ✅ packages + verifies cleanly — `Packaged 173 files, 1.8 MiB`, the packaged
tarball compiles. BUT a **blocking finding**: the crate name **`rosalind` is already taken on
crates.io** — by an unrelated crate (`github.com/antklim/rosalind`, "solutions of problems published on
Rosalind.info", currently v0.10.0). So `cargo publish` as `rosalind` is **not possible**; publishing to
crates.io requires **renaming the published crate** (e.g. `rosalind-genomics`, `rosalind-engine`) in
`Cargo.toml` `[package].name` — a product/naming decision for the maintainer, deferred out of this
increment. The GitHub-release + `install.sh` + Action distribution path (already shipped) is unaffected.

## 4. File-by-file change list

| File | Change |
|---|---|
| `.github/workflows/ci.yml` | Repoint the 2 pip-cache keys off `pyproject.toml`. (MSRV job deferred to 1.3b.) |
| `README.md` | Replace the `space_bounds` test line. |
| `python/README.md` | `ft.array` → `ft.data`. |
| `CONTRIBUTING.md` (new) | Contribution guide. |
| `ARCHITECTURE.md` (new) | Module + phase map. |
| `.github/ISSUE_TEMPLATE/bug_report.md` (new) | Deterministic-repro bug template. |

## 5. Testing / verification

- `cargo build --release` + `cargo test` stay green; `cargo fmt --check` clean.
- `cargo publish --dry-run` succeeds.
- The fixed doc commands actually run: `cargo test --test plan_enforce` passes; the python snippet
  field (`ft.data`) matches `python/rosalind.py`.
- CI green on GitHub (incl. the new `msrv` job).

## 6. Deferred follow-up (noted, not in this increment)

**1.3b — MSRV policy + clippy cleanup + `-D warnings` gate.** One increment, because MSRV, `div_ceil`,
and the clippy lints are the same decision. Scope: (1) pick a real MSRV policy — either accept a recent
MSRV that the current deps build on (verify by install + `cargo check`), or pin the bleeding-edge
transitive deps down via `cargo update --precise` to support a lower floor — then declare `rust-version`
and add an MSRV CI job that actually passes; (2) clear all 82 `clippy --all-targets` lints (needless-borrow
×17, `manual_div_ceil` ×10 — resolved by the MSRV choice, manual-`for` loops, `repeat().take()`, etc.) and
add the `-D warnings` CI gate. The `Cargo.lock` is format v4 (≥1.78) and dep MSRVs climb to ≥1.83 — that
constraint is the starting point.

## 7. References

- `README.md:387`, `python/README.md:19`, `.github/workflows/ci.yml:55,96` — the verified defects.
- `Cargo.toml` — metadata (publish-ready) + the missing `rust-version`.
- `docs/ROADMAP.md` Sprint 1 (QW-1) — this increment.
