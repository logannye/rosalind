# Phase C2 — `rosalind plan` + `--enforce` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn C1's sound working-set accountant into a contract: `rosalind plan` predicts whether a whole-genome call fits a declared budget *before* committing, and `variants --index --enforce` refuses up front (exit 3) or fails loud post-run (exit 4) instead of silently overrunning.

**Architecture:** A pure estimator (`src/call/plan.rs`) computes the data-dependent working set from the index header + declared `--max-depth`/`--max-read-len`, sharing its cost constants with C1's realized accountant so they cannot drift. The *predicted peak RSS* = a **process baseline measured at call time** (`peak_rss_bytes()` after the index/BAM are open) + that working set — apples-to-apples with the realized `peak_rss` the post-run check uses, with no guessed baseline constant. `rosalind plan` renders the breakdown; `variants --enforce` refuses pre-run (exit 3) when predicted > budget and fails post-run (exit 4) when realized > budget.

**Tech Stack:** Rust 1.72 (MSRV), `clap` derive, `cargo test`/`fmt`/`build`. No new dependencies. Builds on C1 (branch `rosalind/phase-c-contract`).

**Spec:** [`docs/superpowers/specs/2026-06-01-phase-c-contract-design.md`](../specs/2026-06-01-phase-c-contract-design.md) §6.

**Resolved design point (extends spec §6.1):** the budget is RSS (consistent with the existing record-only `budget.admits(peak_rss)`). The pre-run *predicted peak* = `measured baseline RSS (after index+BAM open) + estimated working set`. This avoids a guessed process-baseline constant and keeps predicted/realized comparable; the post-run `peak_rss` check is the hard backstop. `--enforce` requires `--max-depth > 0` (an uncapped active set has no a-priori bound). **The exit-4 (post-run) subprocess test is deferred to C3's CI contract suite** — it is only reachable when the model under-predicts, which is brittle to engineer against coarse process RSS; C2 unit-tests the decision and subprocess-tests exit 3 + the generous-budget pass.

---

## File Structure

- **Modify** `src/core/budget.rs` — add 4 shared `pub const` pileup cost constants.
- **Modify** `src/core/mod.rs` — re-export the 4 constants.
- **Modify** `src/pileup/engine.rs` — `current_working_set` uses the shared constants (behavior-identical).
- **Create** `src/call/plan.rs` — `estimate_variants_working_set`, `predicted_peak_rss_bytes`, `render_variants_plan` (pure, unit-tested).
- **Modify** `src/call/mod.rs` — `pub mod plan;` + re-export.
- **Modify** `src/main.rs` — `Plan` subcommand + `run_plan`; `--enforce`/`--max-depth`/`--max-read-len` on `Variants`; `run_variants_index` gains the cap wiring + pre/post enforcement.
- **Create** `tests/plan_enforce.rs` — subprocess tests for `plan` and `--enforce`.

---

## Task 1: Shared pileup cost constants (no behavior change)

**Files:**
- Modify: `src/core/budget.rs`, `src/core/mod.rs`, `src/pileup/engine.rs`

- [ ] **Step 1: Add the constants.** In `src/core/budget.rs`, after the module doc comment (before `pub struct MemoryBudget`), insert:

```rust
/// Per-base cost of a pileup read's reference→read-offset projection map (one
/// `HashMap<u32, usize>` entry ≈ 16 bytes). Shared by the realized accountant
/// (`PileupEngine::current_working_set`) and the `rosalind plan` estimator so the
/// two cannot drift.
pub const PILEUP_MAP_BYTES_PER_BASE: u64 = 16;
/// Per-base cost of a read's `seq` + `qual` byte buffers (1 byte each).
pub const PILEUP_SEQQUAL_BYTES_PER_BASE: u64 = 2;
/// Fixed per-active-read overhead (handles + struct).
pub const PILEUP_PER_READ_OVERHEAD: u64 = 64;
/// Fixed per-engine overhead.
pub const PILEUP_ENGINE_OVERHEAD: u64 = 256;
```

- [ ] **Step 2: Re-export them.** In `src/core/mod.rs`, replace the budget re-export line:

```rust
pub use budget::{MemoryBudget, WorkingSet};
```
with:
```rust
pub use budget::{
    MemoryBudget, WorkingSet, PILEUP_ENGINE_OVERHEAD, PILEUP_MAP_BYTES_PER_BASE,
    PILEUP_PER_READ_OVERHEAD, PILEUP_SEQQUAL_BYTES_PER_BASE,
};
```

- [ ] **Step 3: Use them in the accountant.** In `src/pileup/engine.rs`, update the import line to add the two constants `current_working_set` needs:

```rust
use crate::core::{
    allele_index, AlignedRead, CoreError, Locus, Position, WorkingSet,
    PILEUP_ENGINE_OVERHEAD, PILEUP_MAP_BYTES_PER_BASE, PILEUP_PER_READ_OVERHEAD,
};
```

Then rewrite the body of `current_working_set` to reference the constants (the numeric values are identical — this is a behavior-preserving refactor):

```rust
    pub fn current_working_set(&self) -> WorkingSet {
        // The decoded reference for this contig is resident in the engine.
        let reference_bytes = self.reference.len() as u64;
        // Each active read holds its projection map plus its seq and qual byte
        // buffers; count all three (the map alone is a large undercount,
        // especially for long reads). Constants are shared with the plan estimator.
        let active_bytes: u64 = self
            .active
            .iter()
            .map(|r| {
                (r.ref_to_read.len() as u64) * PILEUP_MAP_BYTES_PER_BASE
                    + r.seq.len() as u64
                    + r.qual.len() as u64
                    + PILEUP_PER_READ_OVERHEAD
            })
            .sum();
        WorkingSet {
            bytes: reference_bytes + active_bytes + PILEUP_ENGINE_OVERHEAD,
        }
    }
```

- [ ] **Step 4: Run the engine tests (behavior unchanged — the 402 assertion still holds)**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine 2>&1 | tail -6`
Expected: PASS (19 tests; `working_set_counts_reference_and_read_byte_buffers` still computes 402 = `16*4 + 4 + 4 + 64` per read `+ 10` ref `+ 256`).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/core/budget.rs src/core/mod.rs src/pileup/engine.rs && git commit -m "refactor(core): extract shared pileup cost consts; accountant uses them (C2 prep)"
```

---

## Task 2: The pure plan estimator (`src/call/plan.rs`)

**Files:**
- Create: `src/call/plan.rs`
- Modify: `src/call/mod.rs`

- [ ] **Step 1: Create the module with its tests.** Write `src/call/plan.rs`:

```rust
//! Pure, testable planning helpers for `rosalind plan` and `variants --enforce`.
//!
//! The estimator predicts the **data-dependent working set** of a bounded
//! whole-genome germline call from the index header plus the declared depth cap
//! and an assumed max read length. It shares its cost constants with the realized
//! accountant (`PileupEngine::current_working_set`) so a passing plan and the
//! realized receipt cannot silently diverge. The *predicted peak RSS* adds a
//! process baseline the caller measures at runtime (`peak_rss_bytes()`), so the
//! prediction is comparable to the realized `peak_rss` the post-run check uses.

use crate::core::{
    MemoryBudget, WorkingSet, PILEUP_ENGINE_OVERHEAD, PILEUP_MAP_BYTES_PER_BASE,
    PILEUP_PER_READ_OVERHEAD, PILEUP_SEQQUAL_BYTES_PER_BASE,
};

/// Estimate the peak streaming working set of a whole-genome germline call: the
/// largest contig's reference (decoded, 1×) + the depth-capped active read set +
/// the fixed engine overhead. A true upper bound when actual reads do not exceed
/// `max_read_len` and depth is capped at `max_depth` (both enforced at runtime —
/// `max_read_len` is the one assumption, with the post-run check as backstop).
pub fn estimate_variants_working_set(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
) -> WorkingSet {
    let per_read = (max_read_len as u64)
        .saturating_mul(PILEUP_MAP_BYTES_PER_BASE + PILEUP_SEQQUAL_BYTES_PER_BASE)
        .saturating_add(PILEUP_PER_READ_OVERHEAD);
    let active = (max_depth as u64).saturating_mul(per_read);
    WorkingSet {
        bytes: largest_contig_len
            .saturating_add(active)
            .saturating_add(PILEUP_ENGINE_OVERHEAD),
    }
}

/// Predicted peak process RSS = a measured process baseline + the estimated
/// working set. Comparable to the realized `peak_rss` the post-run check uses.
pub fn predicted_peak_rss_bytes(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
    baseline_rss_bytes: u64,
) -> u64 {
    baseline_rss_bytes
        .saturating_add(estimate_variants_working_set(largest_contig_len, max_depth, max_read_len).bytes)
}

/// Render the `rosalind plan --index` breakdown: a measured baseline + the
/// working-set components → predicted peak, tagged `[FITS]`/`[REFUSE]` against the
/// budget (or `[no budget]` when none is declared). Deterministic given inputs.
pub fn render_variants_plan(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
    baseline_rss_bytes: u64,
    budget_mb: Option<u64>,
) -> String {
    const MIB: u64 = 1 << 20;
    let per_read = (max_read_len as u64)
        .saturating_mul(PILEUP_MAP_BYTES_PER_BASE + PILEUP_SEQQUAL_BYTES_PER_BASE)
        .saturating_add(PILEUP_PER_READ_OVERHEAD);
    let active = (max_depth as u64).saturating_mul(per_read);
    let predicted =
        predicted_peak_rss_bytes(largest_contig_len, max_depth, max_read_len, baseline_rss_bytes);
    let verdict = match budget_mb {
        Some(mb) => {
            if MemoryBudget::from_mb(mb).admits(predicted) {
                format!("/ budget {mb} MiB  [FITS]")
            } else {
                format!("/ budget {mb} MiB  [REFUSE]")
            }
        }
        None => "[no budget]".to_string(),
    };
    format!(
        "plan: predicted peak RSS (upper bound)\n  \
         process baseline (measured):        {} MiB\n  \
         reference decode (largest contig):  {} MiB\n  \
         active set @ max-depth {}:           {} MiB\n  \
         engine overhead:                    {} MiB\n  \
         -------------------------------------------------\n  \
         predicted peak: ~{} MiB {}\n",
        baseline_rss_bytes / MIB,
        largest_contig_len / MIB,
        max_depth,
        active / MIB,
        PILEUP_ENGINE_OVERHEAD / MIB,
        predicted / MIB,
        verdict,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_grows_with_inputs_and_does_not_overflow() {
        let small = estimate_variants_working_set(1_000, 100, 150).bytes;
        let deeper = estimate_variants_working_set(1_000, 2_000, 150).bytes;
        let bigger_ref = estimate_variants_working_set(1_000_000, 100, 150).bytes;
        assert!(deeper > small, "more depth → larger working set");
        assert!(bigger_ref > small, "larger contig → larger working set");
        // Active term for D=100, L=150: 100 * (150*18 + 64) = 100 * 2764 = 276_400.
        assert_eq!(small, 1_000 + 276_400 + 256);
        let _ = estimate_variants_working_set(u64::MAX, u32::MAX, u32::MAX); // no panic
    }

    #[test]
    fn predicted_peak_is_baseline_plus_working_set() {
        let ws = estimate_variants_working_set(248_000_000, 1_000, 250).bytes;
        let predicted = predicted_peak_rss_bytes(248_000_000, 1_000, 250, 50_000_000);
        assert_eq!(predicted, 50_000_000 + ws);
    }

    #[test]
    fn render_reports_fits_and_refuse() {
        // Tiny working set; generous budget → FITS.
        let fits = render_variants_plan(1_000, 100, 150, 1_000_000, Some(4096));
        assert!(fits.contains("[FITS]"), "generous budget should fit: {fits}");
        // 248 MiB contig + baseline 50 MiB ≫ 64 MiB budget → REFUSE.
        let refuse = render_variants_plan(248 * (1 << 20), 1000, 250, 50 * (1 << 20), Some(64));
        assert!(refuse.contains("[REFUSE]"), "tight budget should refuse: {refuse}");
        // No budget → advisory.
        assert!(render_variants_plan(1_000, 100, 150, 0, None).contains("[no budget]"));
    }
}
```

- [ ] **Step 2: Wire the module.** In `src/call/mod.rs`, add the module declaration (after `pub mod pipeline;`) and a re-export. Add:

```rust
pub mod plan;
```
and add a re-export line:
```rust
pub use plan::{estimate_variants_working_set, predicted_peak_rss_bytes, render_variants_plan};
```

- [ ] **Step 3: Run the plan tests**

Run: `cd ~/rosalind && cargo test -p rosalind --lib call::plan 2>&1 | tail -8`
Expected: PASS (3 tests).

- [ ] **Step 4: Commit**

```bash
cd ~/rosalind && git add src/call/plan.rs src/call/mod.rs && git commit -m "feat(call): pure variants plan estimator (shares cost consts with the accountant) (C2)"
```

---

## Task 3: `rosalind plan` subcommand

**Files:**
- Modify: `src/main.rs` (`Commands` enum, dispatch, new `run_plan`)
- Create: `tests/plan_enforce.rs`

- [ ] **Step 1: Add the `Plan` variant.** In `src/main.rs`, in `enum Commands`, after the `Locate { … }` variant (before the closing `}` of the enum), add:

```rust
    /// Predict whether a job fits a declared memory budget, before committing.
    Plan {
        /// Persisted index (`rosalind index`): predict the bounded whole-genome
        /// `variants` peak. Mutually exclusive with `--reference`.
        #[arg(
            long,
            conflicts_with = "reference",
            required_unless_present = "reference"
        )]
        index: Option<PathBuf>,
        /// Reference FASTA: predict the index BUILD peak (advisory — build is
        /// O(reference); Phase D enforces). Mutually exclusive with `--index`.
        #[arg(long, required_unless_present = "index")]
        reference: Option<PathBuf>,
        /// Max active depth assumed for the `variants` working-set bound.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed for the `variants` working-set bound.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Declared memory budget (MiB) to check feasibility against.
        #[arg(long)]
        budget_mb: Option<u64>,
    },
```

- [ ] **Step 2: Add the dispatch arm.** In `main()`, after the `Commands::Locate { … } => run_locate(…)?,` arm, add:

```rust
        Commands::Plan {
            index,
            reference,
            max_depth,
            max_read_len,
            budget_mb,
        } => run_plan(index, reference, max_depth, max_read_len, budget_mb)?,
```

- [ ] **Step 3: Implement `run_plan`.** Add this function in `src/main.rs` immediately after `run_index` (after its closing `}`):

```rust
/// Predict whether a job fits a declared budget, before committing. `--index`
/// predicts the bounded whole-genome `variants` peak (largest contig + active set
/// @ the declared cap, atop the measured process baseline). `--reference`
/// predicts the index build peak (advisory; build is O(reference)).
fn run_plan(
    index: Option<PathBuf>,
    reference: Option<PathBuf>,
    max_depth: u32,
    max_read_len: u32,
    budget_mb: Option<u64>,
) -> Result<()> {
    use rosalind::call::plan::render_variants_plan;
    use rosalind::genomics::IndexReader;

    if let Some(index_path) = index {
        let loaded = IndexReader::open(&index_path)
            .with_context(|| format!("failed to open index {}", index_path.display()))?;
        let largest = loaded
            .contigs()
            .iter()
            .map(|c| c.length as u64)
            .max()
            .unwrap_or(0);
        // Measure the process baseline now (binary + libs + index mmap header);
        // the per-contig reference decode + active set are modeled on top.
        let baseline = peak_rss_bytes();
        print!(
            "{}",
            render_variants_plan(largest, max_depth, max_read_len, baseline, budget_mb)
        );
    } else {
        let reference = reference.expect("clap guarantees one of --index/--reference");
        let fasta_reader = open_input(&reference)
            .with_context(|| format!("failed to open reference {}", reference.display()))?;
        let total_bp: u64 = FastaReader::new(fasta_reader)
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("failed to parse FASTA {}", reference.display()))?
            .iter()
            .map(|r| r.sequence.len() as u64)
            .sum();
        let estimate = estimate_build_working_set(total_bp);
        match budget_mb {
            Some(mb) => println!("{}", render_plan_line(estimate, MemoryBudget::from_mb(mb))),
            None => println!(
                "plan: est. build peak ~{} MiB (advisory; build is O(reference)) [no budget]",
                estimate.bytes / (1 << 20)
            ),
        }
    }
    Ok(())
}
```

- [ ] **Step 4: Build to verify it compiles**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -8`
Expected: success, 0 warnings. (`estimate_build_working_set`, `render_plan_line`, `MemoryBudget`, `FastaReader`, `open_input`, `peak_rss_bytes` are already imported/in scope in `main.rs` from the `run_index` path.)

- [ ] **Step 5: Write the subprocess test.** Create `tests/plan_enforce.rs`:

```rust
//! CLI contract surface: `rosalind plan` predicts feasibility, and
//! `variants --index --enforce` refuses up front / passes within budget.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

// Build a tiny 2-contig index in a fresh temp dir; return (dir, index_path).
fn build_index() -> (PathBuf, PathBuf) {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("rosalind-plan-{nanos}"));
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, b">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTT\n").unwrap();
    let idx = dir.join("ref.idx");
    let out = Command::new(bin())
        .args(["index", "--reference"])
        .arg(&fa)
        .arg("--output")
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "index build failed: {out:?}");
    (dir, idx)
}

#[test]
fn plan_index_reports_a_breakdown_and_fits_a_generous_budget() {
    let (dir, idx) = build_index();
    let out = Command::new(bin())
        .args(["plan", "--index"])
        .arg(&idx)
        .args(["--budget-mb", "4096"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("predicted peak"), "missing breakdown: {stdout}");
    assert!(stdout.contains("[FITS]"), "generous budget should FIT: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn plan_reference_reports_build_estimate() {
    let (dir, _idx) = build_index();
    let fa = dir.join("ref.fa");
    let out = Command::new(bin())
        .args(["plan", "--reference"])
        .arg(&fa)
        .args(["--budget-mb", "4096"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("plan:"), "missing build plan line: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 6: Run the plan subprocess test**

Run: `cd ~/rosalind && cargo test --test plan_enforce plan_ 2>&1 | tail -15`
Expected: PASS (2 `plan_*` tests).

- [ ] **Step 7: Commit**

```bash
cd ~/rosalind && git add src/main.rs tests/plan_enforce.rs && git commit -m "feat(cli): rosalind plan — predict variants/build peak vs a declared budget (C2)"
```

---

## Task 4: `--max-depth` / `--max-read-len` / `--enforce` flags + cap wiring

**Files:**
- Modify: `src/main.rs` (`Variants` variant, dispatch, `run_variants_index` signature + cap wiring)

- [ ] **Step 1: Add the flags to the `Variants` variant.** In `enum Commands`, in `Variants { … }`, after the `memory_budget_mb` field, add:

```rust
        /// Cap the active read set per position (deterministic downsampling); the
        /// bound `plan`/`--enforce` rely on. `0` = uncapped.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed by the pre-run `--enforce` estimate.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Honor the budget: refuse up front if predicted peak exceeds it (exit 3),
        /// or fail after the run if the realized peak does (exit 4). Requires
        /// `--memory-budget-mb` and `--max-depth > 0`.
        #[arg(long, default_value_t = false)]
        enforce: bool,
```

- [ ] **Step 2: Update the dispatch.** In `main()`, the `Commands::Variants { … }` destructure — add the three new fields to the pattern, and pass them to `run_variants_index`. Replace the destructure field list to include them and replace the `run_variants_index(…)` call:

```rust
        Commands::Variants {
            index,
            reference,
            alignments,
            chrom,
            region_start,
            mapq_threshold,
            output,
            block_size: _,
            quality_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
        } => {
            if let Some(index) = index {
                if chrom.is_some() || region_start != 0 {
                    bail!("--chrom/--region-start are not valid with --index (the whole index is called)");
                }
                run_variants_index(
                    index,
                    alignments,
                    mapq_threshold,
                    output,
                    quality_threshold,
                    memory_budget_mb,
                    max_depth,
                    max_read_len,
                    enforce,
                )?
            } else {
                let reference = reference.expect("clap guarantees one of --index/--reference");
                run_variants(
                    reference,
                    alignments,
                    chrom,
                    region_start,
                    mapq_threshold,
                    output,
                    1024,
                    quality_threshold,
                )?
            }
        }
```

- [ ] **Step 3: Extend `run_variants_index`'s signature + cap wiring.** Change the signature (add the three params) and the `pileup_params` construction. Replace the function signature:

```rust
fn run_variants_index(
    index_path: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    quality_threshold: f32,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
) -> Result<()> {
```

And replace the `pileup_params` construction (currently `PileupParams { min_mapq: mapq_threshold, ..PileupParams::default() }`) with:

```rust
    let pileup_params = PileupParams {
        min_mapq: mapq_threshold,
        // `--max-depth 0` opts out of the cap (then the working set is unbounded
        // and `--enforce` is rejected below).
        max_depth: if max_depth == 0 { None } else { Some(max_depth) },
        ..PileupParams::default()
    };
```

- [ ] **Step 4: Build to verify it compiles** (enforcement logic comes in Task 5; `max_read_len`/`enforce` are unused for now — silence with a leading `let _ = (max_read_len, enforce);` placeholder is NOT needed because Task 5 lands immediately; if building between tasks, expect an `unused variable` warning only).

Run: `cd ~/rosalind && cargo build 2>&1 | tail -10`
Expected: compiles. (Warnings for unused `max_read_len`/`enforce` are acceptable *only* until Task 5; do not commit Task 4 alone — proceed to Task 5 before committing, or accept the transient warning. To keep commits clean, **commit Task 4 + Task 5 together** after Task 5's tests pass.)

- [ ] **Step 5: Run the existing whole-genome gate (cap default 1000 does not change small-test output)**

Run: `cd ~/rosalind && cargo test --test variants_index 2>&1 | tail -10`
Expected: PASS (5 tests — the tiny test inputs are far below depth 1000).

(No commit here — see Task 5.)

---

## Task 5: `--enforce` — refuse pre-run (exit 3), fail post-run (exit 4)

**Files:**
- Modify: `src/main.rs` (`run_variants_index`: pre-run check after opening the source; post-run check replacing the record-only block)
- Modify: `tests/plan_enforce.rs` (add enforce tests)

- [ ] **Step 1: Add the pre-run refuse check.** In `run_variants_index`, immediately after the `let source = StreamingBamSource::new(…)?;` line (and before the `let max_ws = match &output { … }` streaming block), insert:

```rust
    // `--enforce` contract: predict the peak RSS up front (measured baseline +
    // the depth-capped working set) and refuse cleanly if it won't fit — before
    // doing any work. Never a silent OOM.
    if enforce {
        if memory_budget_mb.is_none() {
            bail!("--enforce requires --memory-budget-mb");
        }
        if max_depth == 0 {
            bail!("--enforce requires --max-depth > 0 (an uncapped active set has no a-priori bound)");
        }
        let mb = memory_budget_mb.unwrap();
        let largest = contigs.iter().map(|c| c.length as u64).max().unwrap_or(0);
        let baseline = peak_rss_bytes();
        let predicted =
            rosalind::call::plan::predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline);
        if !MemoryBudget::from_mb(mb).admits(predicted) {
            eprintln!(
                "contract: REFUSE — declared {} MiB, predicted peak ~{} MiB \
                 (largest contig {} MiB + active @ max-depth {} / max-read-len {} \
                 atop a {} MiB baseline). Raise --memory-budget-mb, lower --max-depth, \
                 or drop --enforce.",
                mb,
                predicted / (1 << 20),
                largest / (1 << 20),
                max_depth,
                max_read_len,
                baseline / (1 << 20),
            );
            std::process::exit(3);
        }
    }
```

- [ ] **Step 2: Replace the post-run record-only block with the enforce-aware version.** Replace the existing budget block at the end of `run_variants_index`:

```rust
    if let Some(mb) = memory_budget_mb {
        let budget = MemoryBudget::from_mb(mb);
        if budget.admits(peak_rss) {
            eprintln!("memory: within budget ({mb} MiB)");
        } else {
            eprintln!(
                "memory: EXCEEDED budget {} MiB (realized peak {} MiB) — record-only, run completed",
                mb,
                peak_rss / (1 << 20)
            );
        }
    }
```
with:
```rust
    if let Some(mb) = memory_budget_mb {
        let budget = MemoryBudget::from_mb(mb);
        let within = budget.admits(peak_rss);
        if enforce {
            if within {
                eprintln!("contract: OK — realized peak {} MiB within declared {mb} MiB", peak_rss / (1 << 20));
            } else {
                eprintln!(
                    "contract: VIOLATED — realized peak {} MiB exceeded declared {mb} MiB (output + receipt written)",
                    peak_rss / (1 << 20)
                );
                std::process::exit(4);
            }
        } else if within {
            eprintln!("memory: within budget ({mb} MiB)");
        } else {
            eprintln!(
                "memory: EXCEEDED budget {mb} MiB (realized peak {} MiB) — record-only, run completed",
                peak_rss / (1 << 20)
            );
        }
    }
```

- [ ] **Step 3: Build (now `max_read_len`/`enforce` are used — 0 warnings)**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -8`
Expected: success, 0 warnings.

- [ ] **Step 4: Add the enforce subprocess tests.** Append to `tests/plan_enforce.rs`. These need a coordinate-sorted BAM; reuse the helper pattern from `tests/variants_index.rs` (it builds a tiny sorted BAM via the `common` test module). Add at the top: `mod common;` is NOT used here — instead drive a refuse that needs no BAM by pointing at a missing BAM is wrong (it would error differently). The deterministic, BAM-free way to hit exit 3 is a budget smaller than the measured baseline alone, which the pre-run check catches before reading the BAM. Add:

```rust
#[test]
fn enforce_refuses_up_front_when_budget_below_predicted() {
    // A 1 MiB budget is below the process baseline alone, so the pre-run check
    // refuses with exit code 3 before touching the (here absent) alignments —
    // the refusal is computed from the index + baseline, not the BAM.
    let (dir, idx) = build_index();
    let bam = dir.join("nonexistent.bam"); // never opened: refuse happens first? No —
    // source opens before the check; use a real (empty-but-valid) BAM instead.
    // Build a minimal sorted BAM by calling `variants --index` is circular; instead
    // assert the refusal via the predicted-peak path using a real sorted BAM from
    // the variants_index fixtures is heavy. Keep this test BAM-light: see note.
    let _ = (bam, idx);
    std::fs::remove_dir_all(&dir).ok();
}
```

**Correction (do this instead of the placeholder above):** the pre-run check runs *after* `StreamingBamSource::new`, so a valid sorted BAM is required even to reach the refusal. Rather than synthesize a BAM here, exercise the refuse path through the **already-present** `tests/variants_index.rs` fixture style. Replace the body above with a test that shells out using the sorted BAM that `tests/variants_index.rs` builds. Since that fixture lives in `tests/common`, add `mod common;` and use it:

```rust
mod common;

#[test]
fn enforce_refuses_up_front_when_budget_below_predicted() {
    // Build a tiny index + a matching coordinate-sorted BAM via the shared
    // fixture, then run with an absurdly small budget under --enforce: the
    // pre-run prediction (baseline + working set) exceeds 1 MiB, so the run
    // refuses with exit code 3 and writes no VCF.
    let fx = common::tiny_sorted_fixture(); // (dir, index_path, sorted_bam_path)
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&fx.index)
        .arg("--alignments")
        .arg(&fx.bam)
        .args(["--memory-budget-mb", "1", "--enforce"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(3), "expected refuse exit 3: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("REFUSE"), "missing refuse message: {stderr}");
    fx.cleanup();
}

#[test]
fn enforce_passes_within_a_generous_budget() {
    let fx = common::tiny_sorted_fixture();
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&fx.index)
        .arg("--alignments")
        .arg(&fx.bam)
        .args(["--memory-budget-mb", "4096", "--enforce"])
        .output()
        .unwrap();
    assert!(out.status.success(), "generous budget should pass: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("contract: OK"), "missing OK line: {stderr}");
    fx.cleanup();
}
```

**Before writing the tests, inspect `tests/common` and `tests/variants_index.rs`** to use the *actual* fixture helper name/shape (the names `tiny_sorted_fixture`/`fx.index`/`fx.bam`/`fx.cleanup` are placeholders for whatever the existing fixture exposes). If `tests/common` has no reusable sorted-BAM fixture, lift the BAM-building code from `tests/variants_index.rs::variants_index_matches_reference_on_the_same_sorted_bam` into this test directly. Match the existing fixture API exactly.

- [ ] **Step 5: Run the enforce tests**

Run: `cd ~/rosalind && cargo test --test plan_enforce 2>&1 | tail -15`
Expected: PASS (the `plan_*` tests + `enforce_refuses_up_front_when_budget_below_predicted` exit 3 + `enforce_passes_within_a_generous_budget`).

- [ ] **Step 6: Commit Tasks 4 + 5 together**

```bash
cd ~/rosalind && git add src/main.rs tests/plan_enforce.rs && git commit -m "feat(cli): variants --index --enforce — refuse pre-run (3) / fail post-run (4); --max-depth default 1000 (C2)"
```

---

## Task 6: Full verification

**Files:** none (verification only)

- [ ] **Step 1: Format**

Run: `cd ~/rosalind && cargo fmt --all && cargo fmt --all -- --check 2>&1 | tail -3`
Expected: clean after applying.

- [ ] **Step 2: Zero-warning builds**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -5 && cargo build --release 2>&1 | tail -5`
Expected: both 0 warnings.

- [ ] **Step 3: Full suite**

Run: `cd ~/rosalind && cargo test 2>&1 | grep -E "test result: FAILED|panicked|[1-9][0-9]* failed" | head; cargo test 2>&1 | grep -cE "test result: ok\."`
Expected: no failures; the count of `ok.` sections ≥ the C1 count + 1 (new `plan_enforce` binary).

- [ ] **Step 4: Commit any fmt fixups** (only if Step 1 changed files)

```bash
cd ~/rosalind && git add -A && git commit -m "style: rustfmt fixups (C2)"
```

---

## Self-Review notes

- **Spec §6 coverage:** §6.1 estimator (shared consts) → Tasks 1+2; §6.2 `plan` subcommand (both modes) → Task 3; §6.3 `--enforce` exit 3/4 → Task 5; §6.4 `--max-depth` default 1000 + `--max-read-len` 250 → Task 4. §6.5 tests → Tasks 2/3/5 (exit-4 subprocess test deferred to C3 per the note up top; exit-4 *path* is implemented and exercised by the generous/within branch).
- **Type consistency:** the four `PILEUP_*` consts are defined once in `core/budget.rs`, re-exported in `core/mod.rs`, and consumed identically in `engine.rs` (Task 1) and `call/plan.rs` (Task 2). `predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline)` has the same argument order in `plan.rs` (Task 2), `render_variants_plan` (Task 2), and the `run_variants_index` pre-run check (Task 5). `run_variants_index` arity is updated in lockstep at its definition (Task 4 Step 3) and its only call site (Task 4 Step 2).
- **Fixture caveat (Task 5 Step 4):** the test fixture helper names are placeholders — inspect `tests/common` + `tests/variants_index.rs` and match the real API before writing, or inline the BAM build. This is the one task requiring a look at existing test code first.
