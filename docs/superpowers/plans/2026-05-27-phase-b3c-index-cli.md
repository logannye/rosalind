# Phase B3c — `rosalind index` + load (CLI) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Surface B3b's persistence at the CLI — `rosalind index` builds a multi-contig reference into a portable artifact once; `rosalind locate` memory-maps it and answers exact-match queries without rebuilding — and lay the first record-only `MemoryBudget` hook.

**Architecture:** Two new flat clap subcommands (`Index`, `Locate`) with **thin** handlers in `main.rs` that delegate to existing library APIs (`io::fasta::FastaReader`, `GenomeIndex::from_named_sequences`, `IndexWriter::write_genome_index`, `IndexReader::open` → `GenomeIndexView::locate_exact`). All testable logic (the build working-set estimate, the build receipt, the budget plan line) lives in a new pure library module `genomics/index/report.rs` so it is unit-tested without spawning a process; CLI behavior is covered by `tests/index_cli.rs` invoking the built binary via `env!("CARGO_BIN_EXE_rosalind")`.

**Tech Stack:** Rust 2021 (MSRV 1.72 — no `div_ceil`); `clap` derive (existing CLI style); `anyhow` for handler errors; `blake3` (existing dep) for the receipt identity; no new dependencies.

This is stage **B3c** under the contract thesis in `docs/OPEN_PROBLEMS.md`, from the design `docs/superpowers/specs/2026-05-27-phase-b3c-index-cli-design.md`. It follows B3b (merged, PR #16). **Out of scope (deferred):** wiring `align`/`variants`/`somatic` onto the persisted index (`--index`) and the seed/chain/extend aligner (**B4**); `MemoryBudget` *enforcement* + `rosalind plan`/`verify` (**Phase C**); leaner/sublinear-space construction (**Phase D**). The build is still `O(reference)` RAM — the receipt is honest about realized peak; the budget hook reports but never refuses.

---

## File structure

- `src/genomics/index/report.rs` — **Create.** Pure, unit-tested helpers: `estimate_build_working_set(reference_len) -> WorkingSet` (coarse, record-only), `IndexBuildReport { … }` + `render() -> String` (deterministic receipt), `render_plan_line(estimate, budget) -> String` (`[OK|OVER]`). One responsibility: presentation + planning for `rosalind index`.
- `src/genomics/index/mod.rs` — **Modify.** `mod report;` + `pub use report::{estimate_build_working_set, render_plan_line, IndexBuildReport};`.
- `src/genomics/mod.rs` — **Modify.** Extend the `pub use index::{…}` line with the three new items.
- `src/main.rs` — **Modify.** Add the `Index` and `Locate` `Commands` variants, their `match` arms, and the thin `run_index`/`run_locate` handlers. Extend the `use rosalind::genomics::{…}` import.
- `tests/index_cli.rs` — **Create.** CLI integration gates (build→loadable, receipt, budget-never-refuses, locate-vs-ground-truth, determinism, self-contained load) via the built binary.
- `README.md` — **Modify (Task 4).** A short "build once → query" usage section.

---

## Task 1: `report.rs` — pure build-estimate / receipt / plan-line helpers

**Files:**
- Create: `src/genomics/index/report.rs`
- Modify: `src/genomics/index/mod.rs`, `src/genomics/mod.rs`

- [ ] **Step 1: Create the module with helpers + tests.** Create `src/genomics/index/report.rs`:

```rust
//! Pure, testable presentation + planning helpers for `rosalind index`.
//!
//! Kept out of `main.rs` (a thin CLI handler) so the build working-set estimate,
//! the build receipt, and the budget plan line are unit-tested without spawning a
//! process. Nothing here enforces anything — the `MemoryBudget` plan line is
//! record-only (honor-or-refuse is Phase C).

use crate::core::{MemoryBudget, WorkingSet};

/// A coarse, **record-only** estimate of the peak working set of building a
/// `GenomeIndex` over a reference of `reference_len` bases.
///
/// The build is dominated by SA-IS over the `u32` text (text + suffix array +
/// workspace) plus the in-RAM index structures — roughly **12 bytes per base**.
/// This is an intentionally coarse upper-ish model for the budget *seam*; precise
/// accounting is Phase C and the build cost itself is what Phase D reduces. It is
/// not a guarantee.
pub fn estimate_build_working_set(reference_len: u64) -> WorkingSet {
    const BYTES_PER_BASE: u64 = 12;
    const BASE_OVERHEAD: u64 = 1 << 20; // 1 MiB floor for short references
    WorkingSet {
        bytes: reference_len
            .saturating_mul(BYTES_PER_BASE)
            .saturating_add(BASE_OVERHEAD),
    }
}

/// The deterministic build receipt for a persisted index. Per-run fields (e.g.
/// realized RSS) are intentionally excluded — the caller prints those separately.
#[derive(Debug, Clone)]
pub struct IndexBuildReport {
    /// Path the index was written to.
    pub index_path: String,
    /// Per-contig `(name, length)` in id order.
    pub contigs: Vec<(String, u32)>,
    /// Total reference length in bases.
    pub total_bp: u64,
    /// BLAKE3 of the (uppercased ASCII) reference.
    pub reference_blake3: [u8; 32],
    /// On-disk size of the index file, in bytes.
    pub index_bytes: u64,
}

impl IndexBuildReport {
    /// Render the receipt as a deterministic multi-line string (trailing newline).
    pub fn render(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("index: {}\n", self.index_path));
        out.push_str(&format!(
            "contigs: {} ({} bp total)\n",
            self.contigs.len(),
            self.total_bp
        ));
        for (name, length) in &self.contigs {
            out.push_str(&format!("  {name}\t{length}\n"));
        }
        out.push_str(&format!(
            "reference_blake3: {}\n",
            hex32(&self.reference_blake3)
        ));
        out.push_str(&format!("index_bytes: {}\n", self.index_bytes));
        out
    }
}

/// Render the record-only budget plan line: estimated build peak vs the declared
/// budget, tagged `[OK]` (estimate fits) or `[OVER]` (estimate exceeds — the build
/// proceeds anyway; enforcement is Phase C).
pub fn render_plan_line(estimate: WorkingSet, budget: MemoryBudget) -> String {
    let verdict = if estimate.fits(budget) { "OK" } else { "OVER" };
    format!(
        "plan: est. build peak ~{} MiB / budget {} MiB  [{}]",
        estimate.bytes / (1 << 20),
        budget.bytes / (1 << 20),
        verdict
    )
}

/// Lowercase hex of a 32-byte digest.
fn hex32(bytes: &[u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = vec![0u8; 64];
    for (i, b) in bytes.iter().enumerate() {
        out[2 * i] = HEX[(b >> 4) as usize];
        out[2 * i + 1] = HEX[(b & 0x0f) as usize];
    }
    String::from_utf8(out).expect("hex digits are valid ASCII")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_grows_with_length_and_does_not_overflow() {
        let small = estimate_build_working_set(1_000).bytes;
        let large = estimate_build_working_set(1_000_000).bytes;
        assert!(large > small, "estimate must grow with reference length");
        assert!(large >= 12_000_000, "≈12 bytes/base");
        let _ = estimate_build_working_set(u64::MAX); // must not panic/overflow
    }

    #[test]
    fn receipt_render_is_deterministic_and_contains_fields() {
        let report = IndexBuildReport {
            index_path: "ref.idx".to_string(),
            contigs: vec![("chr1".to_string(), 100), ("chr2".to_string(), 50)],
            total_bp: 150,
            reference_blake3: [0xab; 32],
            index_bytes: 4096,
        };
        let a = report.render();
        assert_eq!(a, report.render(), "render must be deterministic");
        assert!(a.contains("index: ref.idx"));
        assert!(a.contains("contigs: 2 (150 bp total)"));
        assert!(a.contains("  chr1\t100"));
        assert!(a.contains("  chr2\t50"));
        assert!(a.contains(&format!("reference_blake3: {}", "ab".repeat(32))));
        assert!(a.contains("index_bytes: 4096"));
    }

    #[test]
    fn plan_line_reports_ok_and_over() {
        let budget = MemoryBudget::from_mb(100);
        let under = WorkingSet {
            bytes: 50 * (1 << 20),
        };
        let over = WorkingSet {
            bytes: 200 * (1 << 20),
        };
        assert!(render_plan_line(under, budget).ends_with("[OK]"));
        assert!(render_plan_line(over, budget).ends_with("[OVER]"));
        assert!(render_plan_line(under, budget).contains("budget 100 MiB"));
    }
}
```

- [ ] **Step 2: Register + export the module.** In `src/genomics/index/mod.rs`, add `mod report;` with the other `mod` lines and `pub use report::{estimate_build_working_set, render_plan_line, IndexBuildReport};` with the other `pub use` lines.

In `src/genomics/mod.rs`, replace the index re-export line:

```rust
pub use index::{FmIndexView, GenomeIndexView, IndexHeader, IndexReader, IndexWriter, ReferenceIndex};
```

with:

```rust
pub use index::{
    estimate_build_working_set, render_plan_line, FmIndexView, GenomeIndexView, IndexBuildReport,
    IndexHeader, IndexReader, IndexWriter, ReferenceIndex,
};
```

- [ ] **Step 3: Run the unit tests.**

Run: `cargo test --lib index::report 2>&1 | tail -15`
Expected: `estimate_grows_with_length_and_does_not_overflow`, `receipt_render_is_deterministic_and_contains_fields`, `plan_line_reports_ok_and_over` — 3 passed. Then `cargo build --lib 2>&1 | grep -iE 'error|warning'` — none (the new items are used within the module by tests; if a dead-code warning appears for an item not yet used by `main.rs`, it resolves in Task 2 — do NOT add `#[allow]` and do NOT run `cargo fix`).

- [ ] **Step 4: Commit.**

```bash
git add src/genomics/index/report.rs src/genomics/index/mod.rs src/genomics/mod.rs
git commit -m "feat(genomics/index/report): build estimate + receipt + plan-line helpers" \
  -m "Pure, unit-tested helpers for rosalind index: a coarse record-only build working-set estimate, a deterministic IndexBuildReport receipt, and the [OK|OVER] budget plan line. No enforcement. Keeps the CLI handler thin and the logic process-free testable." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `rosalind index` (build + persist + receipt + budget plan line)

**Files:**
- Modify: `src/main.rs`
- Create: `tests/index_cli.rs`

- [ ] **Step 1: Write the failing CLI tests.** Create `tests/index_cli.rs`:

```rust
//! Phase B3c CLI gates: `rosalind index` + `rosalind locate` (build-once → query).

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = env::temp_dir().join(format!("rosalind-b3c-{nanos}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

// chr1/chr2/chr3 each 30 bp (90 bp total). GATTACA occurs in chr3 at local
// positions 0 and 11 (2 hits, both in chr3).
fn write_fasta(dir: &Path) -> PathBuf {
    let p = dir.join("ref.fa");
    std::fs::write(
        &p,
        ">chr1\nACGTACGTNNACGTACGTACGTAAGGCCTT\n\
         >chr2\nTTTTGGGGCCCCAAAANNNNACGTACGTAC\n\
         >chr3\nGATTACATTTTGATTACAGGGGGCCCCAAA\n",
    )
    .unwrap();
    p
}

#[test]
fn index_builds_a_loadable_artifact() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");

    let out = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("run index");
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(idx.exists(), "index file not written");

    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("contigs: 3 (90 bp total)"), "receipt: {stdout}");
    assert!(stdout.contains("  chr1\t30"), "receipt: {stdout}");
    assert!(stdout.contains("index_bytes: "), "receipt: {stdout}");

    // The binary-produced artifact loads via the library reader.
    let loaded = rosalind::genomics::IndexReader::open(&idx).expect("open");
    assert_eq!(loaded.contigs().len(), 3);

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn index_memory_budget_prints_plan_line_and_never_refuses() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");

    // Budget 0 → the estimate is OVER, but the build must STILL succeed
    // (record-only; no enforcement).
    let out = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
            "--memory-budget-mb",
            "0",
        ])
        .output()
        .expect("run index");
    assert!(
        out.status.success(),
        "index must succeed even when over budget"
    );
    assert!(idx.exists(), "index written despite over-budget");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("plan:") && stderr.contains("[OVER]"),
        "expected an [OVER] plan line on stderr, got: {stderr}"
    );

    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run them to verify they fail.**

Run: `cargo test --test index_cli 2>&1 | tail -20`
Expected: the binary builds but `index` is an unknown subcommand → the spawned command fails (non-zero exit), so `assert!(out.status.success())` fails. (This confirms the test harness + binary path work before `Index` exists.)

- [ ] **Step 3: Add the `Index` subcommand + handler.** In `src/main.rs`:

Extend the genomics import (replace the existing `use rosalind::genomics::{ … };` block) to add the index-build items:

```rust
use rosalind::genomics::{
    compare_callsets, create_bam_writer, estimate_build_working_set, read_vcf_variants,
    render_plan_line, sort_bam_deterministic, AlignedRead, BWTAligner, BedIndex, CigarOp,
    CigarOpKind, GenomeIndex, IndexBuildReport, IndexWriter,
};
```

Add `use rosalind::core::MemoryBudget;` with the other top-level `use` lines.

Add a variant to `enum Commands` (after the `Somatic`/`EvalSomatic` variants, before the closing `}`):

```rust
    /// Build a reference index once into a portable, memory-mappable artifact.
    Index {
        /// Reference genome in FASTA (plain or gzip; `-` for stdin). All contigs.
        #[arg(long)]
        reference: PathBuf,
        /// Output path for the index artifact.
        #[arg(short, long)]
        output: PathBuf,
        /// Declared memory budget (MiB) for the build. Records a plan line; does
        /// not enforce (enforcement is a later phase).
        #[arg(long)]
        memory_budget_mb: Option<u64>,
    },
```

Add the dispatch arm in `main()`'s `match cli.command { … }` (alongside the others):

```rust
        Commands::Index {
            reference,
            output,
            memory_budget_mb,
        } => run_index(reference, output, memory_budget_mb)?,
```

Add the handler (place it near the other `run_*` functions, e.g. after `run_eval_somatic`):

```rust
/// Build a multi-contig index from a FASTA and persist it (B3c).
fn run_index(reference: PathBuf, output: PathBuf, memory_budget_mb: Option<u64>) -> Result<()> {
    // Stream every FASTA record (all contigs) into (name, sequence) pairs.
    let fasta_reader = open_input(&reference)
        .with_context(|| format!("failed to open reference {}", reference.display()))?;
    let records: Vec<FastaRecord> = FastaReader::new(fasta_reader)
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to parse FASTA {}", reference.display()))?;
    if records.is_empty() {
        bail!("reference {} contains no FASTA records", reference.display());
    }
    let total_bp: u64 = records.iter().map(|r| r.sequence.len() as u64).sum();

    // Record-only budget plan line, printed BEFORE the build. Never refuses.
    if let Some(mb) = memory_budget_mb {
        let estimate = estimate_build_working_set(total_bp);
        eprintln!("{}", render_plan_line(estimate, MemoryBudget::from_mb(mb)));
    }

    let named: Vec<(String, Vec<u8>)> =
        records.into_iter().map(|r| (r.name, r.sequence)).collect();
    let index = GenomeIndex::from_named_sequences(&named)
        .with_context(|| format!("failed to build index from {}", reference.display()))?;

    IndexWriter::create(&output)
        .with_context(|| format!("failed to create index file {}", output.display()))?
        .write_genome_index(&index)
        .with_context(|| format!("failed to write index to {}", output.display()))?;

    // Deterministic build receipt → stdout.
    let index_bytes = std::fs::metadata(&output)
        .with_context(|| format!("failed to stat index file {}", output.display()))?
        .len();
    let reference_blake3 = *blake3::hash(index.reference()).as_bytes();
    let report = IndexBuildReport {
        index_path: output.display().to_string(),
        contigs: index
            .contigs()
            .iter()
            .map(|c| (c.name.to_string(), c.length))
            .collect(),
        total_bp,
        reference_blake3,
        index_bytes,
    };
    print!("{}", report.render());

    // Realized peak RSS (per-run, informational) → stderr.
    eprintln!("build peak RSS: {} MiB", peak_rss_bytes() / (1 << 20));
    Ok(())
}
```

(`blake3` is an existing dependency reachable as `blake3::hash` without a `use`; `peak_rss_bytes`, `open_input`, `FastaReader`, `FastaRecord`, `bail!`, `Context` are already imported at the top of `main.rs`.)

- [ ] **Step 4: Run the index tests + build.**

Run: `cargo build 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --test index_cli 2>&1 | tail -20`
Expected: `index_builds_a_loadable_artifact` and `index_memory_budget_prints_plan_line_and_never_refuses` PASS. Then `cargo fmt --all`.

- [ ] **Step 5: Commit.**

```bash
git add src/main.rs tests/index_cli.rs
git commit -m "feat(cli): rosalind index — build + persist a reference index" \
  -m "Streams all FASTA contigs -> GenomeIndex::from_named_sequences -> write_genome_index, prints a deterministic build receipt (contigs, total bp, reference BLAKE3, on-disk size) to stdout and realized peak RSS to stderr. --memory-budget-mb prints a record-only [OK|OVER] plan line and never refuses the build. Thin handler over library APIs; logic tested in genomics/index/report + tests/index_cli." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `rosalind locate` (load + exact-match query)

**Files:**
- Modify: `src/main.rs`, `tests/index_cli.rs`

- [ ] **Step 1: Write the failing test.** Append to `tests/index_cli.rs`:

```rust
#[test]
fn locate_matches_in_ram_ground_truth() {
    use rosalind::genomics::GenomeIndex;

    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");
    let build = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("index");
    assert!(build.status.success());

    // In-RAM ground truth over the same sequences.
    let gi = GenomeIndex::from_named_sequences(&[
        ("chr1".to_string(), b"ACGTACGTNNACGTACGTACGTAAGGCCTT".to_vec()),
        ("chr2".to_string(), b"TTTTGGGGCCCCAAAANNNNACGTACGTAC".to_vec()),
        ("chr3".to_string(), b"GATTACATTTTGATTACAGGGGGCCCCAAA".to_vec()),
    ])
    .unwrap();

    // Includes patterns with multiple hits, one with zero hits ("ZZZZ" has an
    // invalid base), and an N-bearing pattern.
    for pat in ["GATTACA", "ACGT", "GGGGG", "NNNN", "TTTTGGGG", "ZZZZ"] {
        let out = Command::new(bin())
            .args(["locate", "--index", idx.to_str().unwrap(), "--pattern", pat])
            .output()
            .expect("locate");
        assert!(
            out.status.success(),
            "locate {pat} failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);

        let mut expected: Vec<String> = gi
            .locate_exact(pat.as_bytes(), 1024)
            .into_iter()
            .map(|l| {
                let name = gi.contigs().by_id(l.contig).unwrap().name.to_string();
                format!("{name}\t{}", l.pos.0)
            })
            .collect();
        expected.sort();
        let mut got: Vec<String> = stdout.lines().map(|s| s.to_string()).collect();
        got.sort();
        assert_eq!(got, expected, "locate mismatch for {pat}");
    }
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --test index_cli::locate_matches_in_ram_ground_truth 2>&1 | tail -15` (or `cargo test --test index_cli locate_matches 2>&1 | tail -15`)
Expected: fails — `locate` is an unknown subcommand (non-zero exit), so the `out.status.success()` assert fails.

- [ ] **Step 3: Add the `Locate` subcommand + handler.** In `src/main.rs`:

Add `IndexReader` to the genomics import (it is the only new item this task needs); the import block from Task 2 becomes:

```rust
use rosalind::genomics::{
    compare_callsets, create_bam_writer, estimate_build_working_set, read_vcf_variants,
    render_plan_line, sort_bam_deterministic, AlignedRead, BWTAligner, BedIndex, CigarOp,
    CigarOpKind, GenomeIndex, IndexBuildReport, IndexReader, IndexWriter,
};
```

Add a variant to `enum Commands`:

```rust
    /// Locate exact occurrences of a pattern in a prebuilt index (load + query).
    Locate {
        /// Index artifact built by `rosalind index`.
        #[arg(long)]
        index: PathBuf,
        /// Pattern to locate (ASCII A/C/G/T/N; case-insensitive).
        #[arg(long)]
        pattern: String,
        /// Maximum number of candidate hits to locate.
        #[arg(long, default_value_t = 1024)]
        max_hits: usize,
    },
```

Add the dispatch arm in `main()`:

```rust
        Commands::Locate {
            index,
            pattern,
            max_hits,
        } => run_locate(index, pattern, max_hits)?,
```

Add the handler (after `run_index`):

```rust
/// Load a prebuilt index and print exact-match loci for `pattern` (B3c). This is
/// a memory-mapped load + exact match — it never rebuilds the index.
fn run_locate(index: PathBuf, pattern: String, max_hits: usize) -> Result<()> {
    let loaded = IndexReader::open(&index)
        .with_context(|| format!("failed to open index {}", index.display()))?;
    let view = loaded
        .genome_view()
        .with_context(|| format!("failed to view index {}", index.display()))?;

    let loci = view.locate_exact(pattern.as_bytes(), max_hits);
    if loci.is_empty() {
        eprintln!("no hits");
        return Ok(());
    }

    let mut stdout = io::BufWriter::new(io::stdout().lock());
    for locus in loci {
        let name = view
            .contigs()
            .by_id(locus.contig)
            .map(|c| c.name.to_string())
            .unwrap_or_else(|| locus.contig.to_string());
        writeln!(stdout, "{name}\t{}", locus.pos.0)?;
    }
    stdout.flush()?;
    Ok(())
}
```

(`IndexReader::open` → `ReferenceIndex::genome_view()` → `GenomeIndexView` from B3b; `view.contigs()` returns `&ContigSet`; `Locus { contig: u32, pos: Position(u32) }` so `locus.pos.0` is the position. `io` and `Write` (for `writeln!`/`flush`) are already imported at the top of `main.rs`.)

- [ ] **Step 4: Run the locate test + full CLI suite + build.**

Run: `cargo build 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --test index_cli 2>&1 | tail -20` (all index + locate tests pass), then `cargo fmt --all`.

- [ ] **Step 5: Commit.**

```bash
git add src/main.rs tests/index_cli.rs
git commit -m "feat(cli): rosalind locate — load + exact-match query a prebuilt index" \
  -m "IndexReader::open -> GenomeIndexView::locate_exact -> prints contig<TAB>pos (sorted). Memory-mapped load, no rebuild, exact-match only (not the aligner — that's B4). Verified against an in-RAM GenomeIndex::locate_exact ground truth over a multi-contig + N-bearing battery." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: determinism + self-contained gates + docs

**Files:**
- Modify: `tests/index_cli.rs`, `README.md`

- [ ] **Step 1: Write the failing gate tests.** Append to `tests/index_cli.rs`:

```rust
#[test]
fn index_build_is_deterministic_via_cli() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx1 = dir.join("a.idx");
    let idx2 = dir.join("b.idx");
    for out in [&idx1, &idx2] {
        let r = Command::new(bin())
            .args([
                "index",
                "--reference",
                fa.to_str().unwrap(),
                "--output",
                out.to_str().unwrap(),
            ])
            .output()
            .expect("index");
        assert!(r.status.success());
    }
    assert_eq!(
        std::fs::read(&idx1).unwrap(),
        std::fs::read(&idx2).unwrap(),
        "two CLI builds of the same reference must be byte-identical"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn locate_works_from_the_artifact_alone() {
    // Build, delete the source FASTA, then locate from the index file alone —
    // proving the load path is self-contained and never rebuilds.
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");
    assert!(Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("index")
        .status
        .success());
    std::fs::remove_file(&fa).unwrap(); // the only source of the sequence is now gone

    let out = Command::new(bin())
        .args(["locate", "--index", idx.to_str().unwrap(), "--pattern", "GATTACA"])
        .output()
        .expect("locate");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        stdout.lines().count(),
        2,
        "GATTACA: 2 hits in chr3, served from the artifact alone"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run them to verify they pass (the behavior already exists from Tasks 2–3).**

Run: `cargo test --test index_cli 2>&1 | tail -20`
Expected: `index_build_is_deterministic_via_cli` and `locate_works_from_the_artifact_alone` PASS, along with all earlier CLI tests. (These are gates over already-implemented behavior; they pass on first run. If determinism fails, the serializer is non-deterministic — that would be a B3b regression; report it rather than weakening the test.)

- [ ] **Step 3: Document the workflow.** In `README.md`, add this section (place it near the existing CLI usage examples, e.g. after the `Variants`/alignment examples):

````markdown
## Build once, query many: the persisted index

Build a portable, memory-mappable index from a (multi-contig) reference once:

```bash
rosalind index --reference genome.fa --output genome.idx
# index: genome.idx
# contigs: 3 (90 bp total)
#   chr1	30
#   ...
# reference_blake3: <hex>
# index_bytes: <n>
```

Then query it in milliseconds — it is memory-mapped, never rebuilt:

```bash
rosalind locate --index genome.idx --pattern GATTACA
# chr3	0
# chr3	11
```

`rosalind index` is deterministic (the `.idx` is byte-identical across builds of
the same reference). `--memory-budget-mb M` prints a record-only build plan line
(`[OK]`/`[OVER]`) — it does not yet enforce the budget (that is a later phase).
`locate` is exact-match only; seed/chain/extend alignment against the persisted
index lands in a later phase.
````

- [ ] **Step 4: Final verification (whole stage).**

Run, and confirm each:
- `cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` — full suite green (report totals; includes `index_cli` + the `index::report` unit tests).
- `cargo build 2>&1 | grep -ic warning` — `0`.
- `cargo fmt --all -- --check` — clean.
- `cargo clippy --lib 2>&1 | grep -E 'index/report'` — no new lints in `report.rs` (MSRV `div_ceil` false positives elsewhere are pre-existing and must not be "fixed").
- No-rebuild (structural): `grep -rn "sais_u32\|BlockedFMIndex::build" src/main.rs` — the `locate` path (`run_locate`) does not appear among the matches; `run_index` legitimately builds (that is the build command), `run_locate` does not. The `locate_works_from_the_artifact_alone` test is the behavioral witness.

- [ ] **Step 5: Commit.**

```bash
git add tests/index_cli.rs README.md
git commit -m "test(cli): B3c determinism + self-contained gates; docs(readme): build-once → query" \
  -m "CLI gates: two rosalind index builds are byte-identical, and rosalind locate serves queries from the .idx alone after the source FASTA is deleted (load never rebuilds). README documents the build-once → query workflow and the record-only --memory-budget-mb plan line." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before the B3c PR)

- `cargo test` — full suite green. Key witnesses: `index::report` unit tests (estimate/receipt/plan); `tests/index_cli` (build→loadable, receipt fields, budget-never-refuses, locate==in-RAM ground truth, determinism, self-contained load).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo fmt --all -- --check` — clean.
- `cargo clippy --lib` — no new lints in `report.rs`.
- The CLI handlers (`run_index`/`run_locate`) are thin; the testable logic lives in `genomics/index/report.rs`.

## Self-Review

- **Spec coverage (`2026-05-27-phase-b3c-index-cli-design.md`):**
  - §3 `rosalind index` (read all contigs → `from_named_sequences` → `write_genome_index`, deterministic receipt) ✔ Task 2 (+ receipt helper Task 1).
  - §4 `rosalind locate` (`open` → `genome_view` → `locate_exact` → `name\tpos`, no rebuild) ✔ Task 3.
  - §5 `MemoryBudget` hook (record-only `--memory-budget-mb` + estimate + plan line; never refuses) ✔ Task 1 (estimate/plan) + Task 2 (flag wiring + the never-refuses test).
  - §6 decomposition ✔ Tasks 1–4.
  - §7 gates — round-trip vs in-RAM ✔ Task 3; no rebuild on load ✔ Task 4 (`locate_works_from_the_artifact_alone` + structural grep); determinism ✔ Task 4; bounded residency ✔ (the B3b `FmIndexView` borrow; `locate` adds no owned copy — the residency witness is B3b's `view_is_a_small_borrow`); budget seam ✔ Task 2.
  - §8 testing (pure-helper unit tests + `CARGO_BIN_EXE` subprocess integration) ✔ Tasks 1–4.
  - §10 decisions (flat subcommands; record-only budget; thin handlers + lib helpers; no `--block-size`) ✔.
- **Type/name consistency:** `estimate_build_working_set(u64) -> WorkingSet`, `IndexBuildReport { index_path, contigs: Vec<(String,u32)>, total_bp, reference_blake3: [u8;32], index_bytes }` + `render()`, `render_plan_line(WorkingSet, MemoryBudget) -> String` are defined in Task 1 and used identically in Task 2. `GenomeIndex::from_named_sequences(&[(String,Vec<u8>)])`, `IndexWriter::create(path)?.write_genome_index(&index)`, `IndexReader::open(path)? -> ReferenceIndex`, `ReferenceIndex::genome_view()?`, `GenomeIndexView::{locate_exact, contigs}`, `ContigSet::by_id`, `Contig { name: Arc<str>, length: u32 }`, `Locus { contig: u32, pos: Position(u32) }` match the merged B3b/B2 APIs. The clap field `memory_budget_mb: Option<u64>` ↔ flag `--memory-budget-mb`; `index`/`locate`/`--pattern`/`--max-hits` match the tests.
- **No placeholders:** every step ships complete code or an exact command + expected output. The Task-1 dead-code-warning note is a transient TDD state resolved in Task 2 (the items are consumed by `run_index`), explicitly flagged — no `#[allow]`, no `cargo fix`.
- **MSRV 1.72:** no `div_ceil`; integer math uses `saturating_mul`/`saturating_add` and `/ (1 << 20)`.
- **Scope:** no `--index` on `align`/`variants` (B4), no aligner (B4), no budget enforcement / `rosalind plan` (Phase C), no `--block-size` (YAGNI) — each deferred.
- **Shared-tree hazard:** per the B3b run, every implementer/reviewer dispatch must forbid `cargo fix`/mutating git, stage only named files, and self-verify `git show --stat HEAD`; the coordinator verifies the commit stat + clean tree at each task boundary.
