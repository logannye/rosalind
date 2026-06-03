# Sprint 1.3 — Forkability Hygiene: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the verified day-one forkability defects: broken doc commands, an unenforced MSRV, and missing contribution scaffolding; confirm crates.io packaging.

**Architecture:** Config + docs only (no production code). Verification is "the fixed command actually runs" + "CI green incl. a new MSRV job" + "`cargo publish --dry-run` packages."

**Tech Stack:** Cargo metadata, GitHub Actions YAML, Markdown.

**Spec:** `docs/superpowers/specs/2026-06-02-sprint1-3-forkability-hygiene-design.md`

---

### Task 1: Doc-drift fixes (3 verified defects)

**Files:** `README.md`, `python/README.md`, `.github/workflows/ci.yml`

- [ ] **Step 1: Fix the README `space_bounds` line**

In `README.md` (~line 387), replace:

```
cargo test --test space_bounds      # working-set scaling checks for the streaming evaluator
```

with:

```
cargo test --test plan_enforce      # the memory contract: plan / --enforce / verify / governor
```

- [ ] **Step 2: Fix the python README field name**

In `python/README.md` (~line 19), replace `print(ft.array.shape)` with `print(ft.data.shape)`.

- [ ] **Step 3: Repoint the CI pip-cache keys**

In `.github/workflows/ci.yml`, replace BOTH occurrences (lines ~55 and ~96) of:

```
          key: ${{ runner.os }}-pip-${{ hashFiles('pyproject.toml') }}
```

with:

```
          key: ${{ runner.os }}-pip-${{ hashFiles('.github/workflows/ci.yml') }}
```

- [ ] **Step 4: Verify the fixed commands are real**

Run: `cargo test --test plan_enforce 2>&1 | grep "test result"` → PASS (the README command resolves).
Run: `grep -n "ft.data" python/rosalind.py` → the `data` field exists on `FeatureTable`.
Run: `grep -rn "space_bounds\|pyproject.toml\|ft.array" README.md python/README.md .github/workflows/ci.yml` → no matches (all fixed).

- [ ] **Step 5: Commit**

```bash
git add README.md python/README.md .github/workflows/ci.yml
git commit -m "fix(docs): repair day-one doc-drift (space_bounds test, ft.array->ft.data, pip-cache key)"
```

---

### Task 2: MSRV — DEFERRED to 1.3b (execution finding 2026-06-02)

**Skip this task.** The documented 1.72 MSRV is not buildable: `Cargo.lock` is format v4 (≥1.78) and the
pinned transitive deps climb to ≥1.83 (`proptest` ICU chain, `icu_properties_data`). Supporting a low
MSRV needs deliberate `cargo update --precise` dep-pinning — the same decision as the deferred
`div_ceil`/clippy work — so MSRV moves to **1.3b** (see the spec §3.1/§6). No `rust-version` or MSRV CI
job ships in 1.3. The pip-cache-key fix from Task 1 already covers the only `ci.yml` change in 1.3.

**Files:** none (deferred).

- [ ] **Step 1: Add `rust-version` to `Cargo.toml`**

In `Cargo.toml` `[package]`, add after the `edition = "2021"` line:

```toml
rust-version = "1.72"
```

- [ ] **Step 2: Verify the crate actually builds on 1.72 locally**

Run: `rustup install 1.72.0 && cargo +1.72.0 check --all-targets 2>&1 | tail -5`
Expected: PASS. **If it fails to compile on 1.72**, bump `rust-version` in `Cargo.toml` to the lowest
version that compiles (try `1.74.0`, then `1.75.0`) and re-run — do NOT weaken the check. Record the
chosen version.

- [ ] **Step 3: Add the MSRV CI job**

In `.github/workflows/ci.yml`, add a new job (mirror the existing job structure — `runs-on:
ubuntu-latest`, an `actions/checkout@v4` step, then a pinned toolchain). Append this job under `jobs:`
(use the SAME `rust-version` value chosen in Step 2):

```yaml
  msrv:
    name: MSRV (1.72)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions-rs/toolchain@v1
        with:
          toolchain: 1.72.0
          profile: minimal
          override: true
      - name: cargo check on the declared MSRV
        run: cargo check --all-targets
```

(If the repo's other jobs use `dtolnay/rust-toolchain` instead of `actions-rs/toolchain`, match that
action for consistency — check the existing `rust`/`python` jobs and copy their toolchain step.)

- [ ] **Step 4: Verify YAML + the full suite still build**

Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"` → no warnings.
Run: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))" && echo "yaml ok"` (or
`ruby -ryaml -e "YAML.load_file('.github/workflows/ci.yml')"` if no python yaml) → parses.

- [ ] **Step 5: Commit**

```bash
git add Cargo.toml Cargo.lock .github/workflows/ci.yml
git commit -m "ci(msrv): declare rust-version and enforce it with a 1.72 build job"
```

---

### Task 3: Contribution scaffolding

**Files:** `CONTRIBUTING.md` (new), `ARCHITECTURE.md` (new), `.github/ISSUE_TEMPLATE/bug_report.md` (new)

- [ ] **Step 1: Write `CONTRIBUTING.md`**

Create `CONTRIBUTING.md`:

```markdown
# Contributing to Rosalind

Thanks for building on Rosalind. This guide gets you productive fast and explains the two
guarantees every change must preserve.

## Build, test, format

```bash
cargo build --release          # the rosalind binary + library
cargo test                     # full unit + integration suite
cargo fmt --check              # formatting (run `cargo fmt` to fix)
```

The MSRV is **Rust 1.72** (`rust-version` in `Cargo.toml`); CI enforces it.

## Two invariants you must not break

1. **Determinism.** Primary artifacts (VCF, gVCF, feature TSV) are byte-identical across repeated runs
   on identical inputs. See [`docs/determinism.md`](docs/determinism.md). If you touch the calling or
   egress path, add/keep a repeat-run byte-equality test.
2. **The memory contract.** `rosalind plan` predicts a peak that must upper-bound the realized peak;
   `--enforce` honors it (refuse exit 3 / governor exit 4); the receipt is self-hashing. See
   [`CONTRACT.md`](CONTRACT.md). A streaming stage's working set must stay bounded by coverage, not
   input size.

## How we work: brainstorm → spec → plan → TDD

Non-trivial changes go through a short design loop. Specs live in `docs/superpowers/specs/`, plans in
`docs/superpowers/plans/`. Write a failing test first; keep commits small; one branch per increment;
open a PR for review (CI must be green before merge).

## Extending the engine

The fastest way to add your own per-locus analytic is the **ColumnKit SDK**: implement one
`ColumnAnalyzer` method and run it through `run_bounded_whole_genome` to inherit bounded memory,
determinism, and a verifiable receipt — no contract code to re-derive. See
[`examples/columnkit_coverage.rs`](examples/columnkit_coverage.rs) and [`CONTRACT.md`](CONTRACT.md).

## Reporting bugs

Use the bug-report template. Include the exact command, the inputs, and the receipt hash
(`manifest_blake3`) + `rosalind --version` — Rosalind's determinism makes most bugs exactly
reproducible.
```

- [ ] **Step 2: Write `ARCHITECTURE.md`**

Create `ARCHITECTURE.md`:

```markdown
# Rosalind Architecture

A map of the codebase, so you know where a change belongs. Companion to
[`docs/ROADMAP.md`](docs/ROADMAP.md) (where we're going) and
[`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md) (the research thesis).

## The kernel

Rosalind's spine is a **bounded, deterministic `PileupColumn` stream**. Variant calling, gVCF, the
feature substrate, and any ColumnKit analyzer all consume that one stream, so they inherit the same
memory contract and byte-reproducibility.

## Modules (`src/`)

| Module | Responsibility |
|---|---|
| `core/` | Shared types: `Locus`/`ContigSet`/`AlignedRead`, the memory `budget` + `WorkingSet`, the runtime `governor`, the typed `CoreError`. |
| `io/` | Streaming FASTA/FASTQ (gzip/bgzf auto-detect), BAM/VCF read+write (via rust-htslib), decompression. |
| `genomics/` | The FM-index: `suffix_array` (SA-IS), `fm_index`/`fm_backing`, `rank_select`, the persisted `index/` (format, mmap `view`), `compressed_dna`, the `bwt_aligner`, deterministic `sort`, and `eval/` (truth-set comparison). |
| `pileup/` | The streaming pileup `engine` + `column` + read `source` — the bounded kernel. |
| `call/` | Consumers of the kernel: `germline`/`somatic`/`gvcf` calling, the `features` egress, the `columnkit` SDK, the `plan` estimator, fleet `pack`, the `pipeline`/`whole_genome` drivers. |
| `provenance/` | The canonical-JSON, self-hashing BLAKE3 run receipt (`RunManifest`). |
| `util/` | Process RSS (`getrusage`), mmap helpers. |

## Roadmap phases → code

- **A** (streaming pileup + calling): `pileup/`, `call/{germline,somatic,gvcf,pipeline}`.
- **B** (genome-scale index): `genomics/{suffix_array,fm_index,index}`, `io/bam`.
- **C** (the memory contract): `core/{budget,governor}`, `call/plan`, `provenance/`, `verify` in `main.rs`.
- **D** (sublinear-space index *build* — research): `genomics/suffix_array` + a future external-memory constructor. See `docs/OPEN_PROBLEMS.md`.
- **E** (reach + ML substrate): `call/{features,columnkit}`, `python/`, `genomics/eval`.

## CLI

`src/main.rs` is the CLI surface (`index`, `align`, `sort`, `variants`, `somatic`, `features`, `plan`,
`pack`, `verify`, `locate`, `eval-*`), thin over the library.
```

- [ ] **Step 3: Write the bug-report issue template**

Create `.github/ISSUE_TEMPLATE/bug_report.md`:

```markdown
---
name: Bug report
about: A reproducible bug (Rosalind's determinism makes most bugs exactly reproducible)
title: ""
labels: bug
---

**What happened vs. what you expected**

**Exact command**
```sh
rosalind ...
```

**Inputs** (sizes + how to obtain them; a tiny fixture that reproduces it is ideal)

**Receipt** — paste the run's `manifest_blake3` (and `contract_verdict`) from the `.manifest.json`, or
attach the manifest.

**Environment**
- `rosalind --version`:
- OS / arch:
- Installed via (release tarball / `cargo build` / Action):
```

- [ ] **Step 4: Verify the docs render + links resolve**

Run: `ls CONTRIBUTING.md ARCHITECTURE.md .github/ISSUE_TEMPLATE/bug_report.md` → all exist.
Run: `grep -oE "\[.*\]\(([^)]+)\)" CONTRIBUTING.md ARCHITECTURE.md | grep -oE "\(([^)]+)\)" | tr -d '()' | grep -vE "^https?:" | while read p; do [ -e "$p" ] || echo "BROKEN LINK: $p"; done; echo "link check done"`
Expected: no `BROKEN LINK` lines (relative links resolve).

- [ ] **Step 5: Commit**

```bash
git add CONTRIBUTING.md ARCHITECTURE.md .github/ISSUE_TEMPLATE/bug_report.md
git commit -m "docs: add CONTRIBUTING, ARCHITECTURE, and a deterministic-repro issue template"
```

---

### Task 4: Publish readiness (dry-run only)

**Files:** none (verification), or small `Cargo.toml` metadata fix if blocked.

- [ ] **Step 1: Dry-run the crates.io package**

Run: `cargo publish --dry-run --allow-dirty 2>&1 | tail -20`
Expected: `Packaged …` / `Verifying …` with no error. (`--allow-dirty` only because the working tree
may have uncommitted plan edits; the real publish is NOT run.)

- [ ] **Step 2: If it fails, fix the blocker (only if needed)**

Common blockers + fixes (apply ONLY the one that fires; otherwise skip):
- "X files … large" or unwanted files in the package → add an `exclude = ["results/", "target/", ".venv/", "tests/golden/", "tests/snapshots/"]` (or the specific dirs named) to `[package]` in `Cargo.toml`.
- A missing required field → it is already present (`description`/`license`/`repository`); only act on the exact message.

Re-run Step 1 until it packages. **Do NOT run `cargo publish` without `--dry-run`.**

- [ ] **Step 3: Commit (only if Cargo.toml changed)**

```bash
git add Cargo.toml
git commit -m "chore(publish): make the crate package cleanly for crates.io (dry-run verified)"
```

---

## Final verification

- [ ] `cargo test` green; `cargo build --release` no warnings; `cargo fmt --check` clean.
- [ ] `cargo +<msrv>.0 check --all-targets` passes (MSRV honored).
- [ ] `cargo publish --dry-run --allow-dirty` packages cleanly.
- [ ] `grep -rn "space_bounds\|pyproject.toml\|ft.array" README.md python/README.md .github/workflows/ci.yml` → no matches.
- [ ] CI green on GitHub (incl. the new `msrv` job).
