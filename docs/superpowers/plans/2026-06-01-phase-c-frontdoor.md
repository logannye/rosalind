# Move #4 — The Front Door Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Rosalind's shipped breakthrough — the bounded-memory contract — discoverable: a forker's first five minutes (crate-root `use`, `cargo doc`, README, an example) land on the genomics product, not the √t theory layer or the non-bounded plugin path.

**Architecture:** Pure positioning / docs / re-exports / examples — **no code-behavior changes.** Curate a genomics product surface at the crate root; rewrite `lib.rs`'s rustdoc to lead with the contract (with a runnable doctest) and demote √t to an honest research footer; add `CONTRACT.md`; rewrite the README to lead with the contract one-command story and route builders to the bounded `PileupColumn` substrate (plugin lineage demoted+labeled); add a substrate cookbook example; prove the runnable demo with a smoke test; reframe GitHub issue #3 (confirm-first).

**Tech Stack:** Rust 1.72 (MSRV), `cargo test`/`fmt`/`build`/`doc`, Markdown. No new dependencies. Built on `rosalind/phase-c-frontdoor` (stacked on Phase C PR #21).

**Spec:** [`docs/superpowers/specs/2026-06-01-phase-c-frontdoor-design.md`](../specs/2026-06-01-phase-c-frontdoor-design.md).

**§3.5 refinement (flagged during planning):** the in-house aligner is single-contig, `--index` is BAM-only, and there is no SAM→BAM converter — so a *multi-contig* fixture cannot be produced out-of-the-box with Rosalind alone. The runnable bundled demo is therefore **single-contig, fully in-house** on the existing `examples/data/illumina_toy/` fixture; the README documents the **multi-contig flagship** as a command (bring-your-own-aligner BAM). No `generate_toy_data.py` change.

---

## File Structure

- **Modify** `src/lib.rs` — crate-root genomics re-exports; rustdoc rewrite (contract-first + runnable doctest + √t research footer).
- **Modify** `src/pileup/mod.rs` — fix the false "plugins build on" docline.
- **Create** `CONTRACT.md` (repo root) — the authoritative contract doc.
- **Modify** `README.md` — contract-first lead, Extend rewrite, honest brand, demo block, √t footer.
- **Create** `examples/custom_pileup_analytics.rs` — non-caller `PileupColumn`-iterator cookbook.
- **Create** `tests/frontdoor_demo.rs` — smoke test proving the README's in-house demo commands run.

---

## Task 1: Crate-root genomics re-exports + fix the pileup docline

**Files:**
- Modify: `src/lib.rs` (re-export block, lines ~56–61)
- Modify: `src/pileup/mod.rs` (module docstring)

- [ ] **Step 1: Add the genomics product surface to the crate-root re-exports.** In `src/lib.rs`, the current block is:

```rust
// Re-exports for convenience
pub use algebra::{AlgebraicEngine, FiniteField};
pub use blocking::{BlockSummary, MovementLog};
pub use ledger::StreamingLedger;
pub use machine::{Configuration, Move, State, Symbol, Transition, TuringMachine};
pub use tree::{CompressedTree, TreeNode};
```

Replace it with (genomics product surface first; the √t layer kept but regrouped under a research comment):

```rust
// ── Genomics product surface — what builders compose on ───────────────────────
// The bounded streaming substrate:
pub use pileup::{Obs, PileupColumn, PileupEngine, PileupParams, ReadSource, SliceSource};
pub use io::bam::StreamingBamSource;
// The bounded whole-genome germline drive + calls:
pub use call::{
    call_germline_region_streaming, call_germline_whole_genome, GermlineCall, GermlineParams,
};
// The memory contract (declare → plan → honor → verify):
pub use call::{estimate_variants_working_set, predicted_peak_rss_bytes};
pub use core::{MemoryBudget, WorkingSet};
// Build-once → mmap index + the reproducibility receipt:
pub use genomics::{GenomeIndex, IndexReader, ReferenceView};
pub use provenance::RunManifest;

// ── Research layer (√t space-bounded simulation; Phase D — see OPEN_PROBLEMS) ──
pub use algebra::{AlgebraicEngine, FiniteField};
pub use blocking::{BlockSummary, MovementLog};
pub use ledger::StreamingLedger;
pub use machine::{Configuration, Move, State, Symbol, Transition, TuringMachine};
pub use tree::{CompressedTree, TreeNode};
```

- [ ] **Step 2: Rewrite the crate-level rustdoc with a runnable doctest.** In `src/lib.rs`, replace the entire top doc block (lines 1–26, from `//! # O(√t) Space Simulation via Height Compression` through the end of the `//! ```` usage block) with:

```rust
//! # Rosalind — a deterministic, low-memory genomics engine
//!
//! Call variants across a whole genome on a laptop, with memory you can **predict
//! and verify**, and results that are **byte-for-byte reproducible**. Rosalind
//! treats memory as a *contract*: you declare a RAM budget, [`rosalind plan`] tells
//! you up front whether the job fits, the run honors it (fits-or-refuses cleanly —
//! never a silent OOM-kill), and `rosalind verify` re-checks a BLAKE3 receipt
//! proving the realized peak landed inside your budget.
//!
//! The kernel is a streaming, CIGAR-aware **pileup column stream** bounded by local
//! coverage, not input size — a substrate you can compute arbitrary per-locus
//! analytics on. Variant calling is the first consumer, not the whole product.
//!
//! ```
//! use std::sync::Arc;
//! use rosalind::{PileupEngine, PileupParams, SliceSource};
//! use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
//!
//! // One 4bp read "ACGT" aligned at chr0:0 over the reference "ACGT".
//! let read = AlignedRead {
//!     contig: 0,
//!     pos: Position(0),
//!     mapq: 60,
//!     flags: SamFlags(0),
//!     cigar: vec![CigarOp::new(CigarOpKind::Match, 4)],
//!     seq: Arc::from(b"ACGT".to_vec().into_boxed_slice()),
//!     qual: Arc::from(vec![40u8; 4].into_boxed_slice()),
//! };
//! let reference: Arc<[u8]> = Arc::from(b"ACGT".to_vec().into_boxed_slice());
//!
//! // The bounded pileup substrate: one PileupColumn per covered position.
//! let mut engine =
//!     PileupEngine::new(SliceSource::new(vec![read]), reference, 0, 0..4, PileupParams::default());
//! let first = engine.next().unwrap().unwrap();
//! assert_eq!(first.depth(), 1);
//! ```
//!
//! ## Research direction (Phase D)
//!
//! Rosalind is also a research vehicle for **space-bounded genomics**: a `~√t`
//! (square-root-space) evaluation framework (Williams 2025; Cook–Mertz 2024) as a
//! continuous space/time knob, aimed at **sublinear-space index construction**.
//! That layer is future work — not yet load-bearing — tracked in
//! [`docs/OPEN_PROBLEMS.md`](https://github.com/logannye/rosalind/blob/main/docs/OPEN_PROBLEMS.md).
```

(`[`rosalind plan`]` is intentionally plain text in prose — it renders as code; no intra-doc link is implied. Leave the `#![warn(missing_docs, …)]` attribute block that follows untouched.)

- [ ] **Step 3: Fix the false pileup docline.** In `src/pileup/mod.rs`, the docstring says *"the single substrate that variant callers and plugins build on."* Replace that sentence:

```rust
//! `PileupEngine` consumes coordinate-sorted reads and yields one `PileupColumn`
//! per covered reference position. It is the single bounded-memory substrate the
//! germline/somatic callers build on; build your own bounded per-locus analytics
//! over the same stream (see `examples/custom_pileup_analytics.rs`). The legacy
//! `GenomicPlugin`/`framework` lineage is separate and NOT memory-bounded.
//! Reference as `crate::pileup::…`.
```

- [ ] **Step 4: Build + doctest + verify the substrate surface resolves**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -3 && cargo test --doc 2>&1 | tail -8`
Expected: 0 warnings; the doc-test for `src/lib.rs` runs and passes (1 doctest). If `cargo doc` is desired: `cargo doc --no-deps` lands on the genomics intro.

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/lib.rs src/pileup/mod.rs && git commit -m "docs(lib): crate-root genomics surface + contract-first rustdoc (√t demoted) (Move #4)"
```

---

## Task 2: `CONTRACT.md`

**Files:**
- Create: `CONTRACT.md` (repo root)

- [ ] **Step 1: Write `CONTRACT.md`.** Author the file with these sections (concrete content, not placeholders):

  1. **Title + one-line promise:** "Memory is a contract, not a hope." The honest brand line verbatim: *"Rosalind never silently OOM-kills you — it fits, or it tells you up front, and proves the realized peak with a receipt."* (Explicitly NOT "never refuses" — note graceful degrade/spill is Phase D.)
  2. **The four verbs**, each with a copy-pasteable command against the bundled fixture:
     - **Declare** — `--memory-budget-mb N`.
     - **Predict** — `rosalind plan --index g.idx --max-depth 1000 --max-read-len 250 --budget-mb 2048` → the `[FITS]`/`[REFUSE]` breakdown.
     - **Honor** — `rosalind variants --index g.idx --alignments s.sorted.bam --memory-budget-mb 2048 --enforce -o s.vcf` → refuses up front (exit 3) if predicted > budget, fails loud (exit 4) if realized > budget, else completes; without `--enforce` it is record-only.
     - **Verify** — `rosalind verify --manifest s.vcf.manifest.json` → re-hashes inputs/outputs + re-checks the recorded peak vs budget without re-running (exit 5 on mismatch).
  3. **What's bounded (honest scope):** the germline `variants --index` path (peak ≈ largest contig + capped active set, independent of BAM size); somatic is region-bounded; index *build* is O(reference) (Phase D); the engine is single-threaded (no thread-invariance claim).
  4. **Extend — build on the bounded substrate:** route to the `PileupColumn` iterator (`PileupEngine` over a `ReadSource`), pointing at `examples/custom_pileup_analytics.rs`; one short Rust snippet. Then a **Legacy / non-bounded** note: the `GenomicPlugin` trait, `src/framework/`, and the Python `run_rna_seq_plugin` demo still work but do **not** inherit the memory contract — prefer the substrate for bounded work.
  5. **Reproducibility:** the BLAKE3 canonical-JSON receipt; `verify` re-checks it; identical inputs → byte-identical VCF.

- [ ] **Step 2: Commit**

```bash
cd ~/rosalind && git add CONTRACT.md && git commit -m "docs: CONTRACT.md — the memory contract + bounded-substrate extension guide (Move #4)"
```

---

## Task 3: README rewrite

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Rewrite the headline + first command block.** Keep the existing lead paragraph's spirit, but make the first runnable block the **contract loop** (declare → plan → enforce → verify) rather than just `variants`. Add, near the top "headline" section, the four-verb story:

```bash
# Build a portable, mmap-able index of your reference — once.
rosalind index --reference genome.fa --output genome.idx

# Will my whole-genome call fit in 2 GB? Ask before committing a byte.
rosalind plan --index genome.idx --max-depth 1000 --budget-mb 2048

# Call across all contigs from a coordinate-sorted BAM, honoring the budget.
rosalind variants --index genome.idx --alignments sample.sorted.bam \
  --memory-budget-mb 2048 --enforce -o sample.vcf

# Re-check the receipt later — no re-run — to prove it fit and is reproducible.
rosalind verify --manifest sample.vcf.manifest.json
```

Add one sentence: link `CONTRACT.md` ("the full contract: [CONTRACT.md](CONTRACT.md)").

- [ ] **Step 2: Honest brand pass.** Search the README for "never refuses" / over-claims about memory and replace with the honest brand: *"never silently OOM-kills you — it fits, or it tells you up front."* (The current README §"Why it matters" / "memory as a contract" language is close; tighten any absolute claims. `--memory-budget-mb` is no longer only record-only — note `--enforce` makes it honored.)

- [ ] **Step 3: Rewrite the "Extend" section.** Replace the current "Extend" bullets (which lead with `GenomicPlugin`) with substrate-first:

```markdown
## Extend

Rosalind's kernel is a **bounded, deterministic `PileupColumn` stream** — build
your own per-locus analytics (coverage, QC, methylation, ML features) over it and
inherit bounded memory + determinism for free:

- **Rust** — consume the `PileupEngine` iterator over any `ReadSource`. See
  [`examples/custom_pileup_analytics.rs`](examples/custom_pileup_analytics.rs) for a
  non-caller consumer computing per-locus coverage.
- **CLI** — compose subcommands over pipes; `variants --index` is the first consumer.

> **Legacy / non-bounded:** the `GenomicPlugin` trait (`src/plugin/`), the
> `framework/` evaluator, and the Python `run_rna_seq_plugin` demo still work but do
> **not** inherit the memory contract. Prefer the `PileupColumn` substrate for
> bounded work.
```

- [ ] **Step 4: Add the runnable in-house demo + document the multi-contig flagship.** In the "Use it" area, add a **runnable, single-contig, fully-in-house** demo on the bundled fixture, and clearly mark the multi-contig path as the production flagship:

```markdown
### Try the contract end-to-end (bundled data, in-house tools only)

```bash
D=examples/data/illumina_toy
rosalind index  --reference $D/reference.fa --output /tmp/toy.idx
rosalind sort   --input $D/alignments.bam   --output /tmp/toy.sorted.bam
rosalind plan   --index /tmp/toy.idx --budget-mb 512
rosalind variants --index /tmp/toy.idx --alignments /tmp/toy.sorted.bam \
  --memory-budget-mb 512 --enforce -o /tmp/toy.vcf
rosalind verify --manifest /tmp/toy.vcf.manifest.json
```

(This bundled demo is **single-contig** because Rosalind's own aligner is
single-contig. For **whole-genome** calling, align with bwa-mem2/minimap2, sort,
and bring the coordinate-sorted BAM to `variants --index` — which calls *every*
contig in bounded memory.)
```

- [ ] **Step 5: Add the √t research footer.** Ensure the README's space-bounded/√t discussion is a clearly-labeled "Research direction (Phase D)" section near the end (not the lead), honest that it is future work, linking `docs/OPEN_PROBLEMS.md`. (The current "Why it matters" closing paragraph + roadmap already gesture at this — consolidate into one honest footer.)

- [ ] **Step 6: Commit**

```bash
cd ~/rosalind && git add README.md && git commit -m "docs(readme): lead with the contract; substrate-first Extend; honest brand; in-house demo (Move #4)"
```

---

## Task 4: Substrate cookbook example

**Files:**
- Create: `examples/custom_pileup_analytics.rs`

- [ ] **Step 1: Write the example.** Create `examples/custom_pileup_analytics.rs` — a non-caller consumer computing per-locus coverage + a simple low-MAPQ count directly over the `PileupColumn` iterator:

```rust
//! Cookbook: build your own bounded, deterministic per-locus analytics over the
//! `PileupColumn` substrate — no variant calling. Run with:
//!   cargo run --example custom_pileup_analytics

use std::sync::Arc;

use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
use rosalind::{PileupColumn, PileupEngine, PileupParams, SliceSource};

fn read(pos: u32, seq: &[u8]) -> AlignedRead {
    AlignedRead {
        contig: 0,
        pos: Position(pos),
        mapq: 60,
        flags: SamFlags(0),
        cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
        seq: Arc::from(seq.to_vec().into_boxed_slice()),
        qual: Arc::from(vec![40u8; seq.len()].into_boxed_slice()),
    }
}

fn main() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGTACGT".to_vec().into_boxed_slice());
    let reads = vec![read(0, b"ACGT"), read(2, b"GTAC"), read(4, b"ACGT")];

    // The substrate: one PileupColumn per covered position, bounded by coverage.
    let engine = PileupEngine::new(
        SliceSource::new(reads),
        Arc::clone(&reference),
        0,
        0..reference.len() as u32,
        PileupParams::default(),
    );

    // A custom per-locus metric — here, depth — computed without any calling.
    println!("pos\tref\tdepth");
    let mut total_depth = 0u64;
    for column in engine {
        let col: PileupColumn = column.expect("pileup column");
        total_depth += col.depth() as u64;
        println!(
            "{}\t{}\t{}",
            col.locus.pos.0,
            col.ref_base as char,
            col.depth()
        );
    }
    println!("# total observed depth across covered positions: {total_depth}");
}
```

- [ ] **Step 2: Run it**

Run: `cd ~/rosalind && cargo run --example custom_pileup_analytics 2>&1 | tail -8`
Expected: a `pos\tref\tdepth` table (positions 0..11 with depths peaking where reads overlap) + the total line; exit 0.

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add examples/custom_pileup_analytics.rs && git commit -m "docs(example): custom_pileup_analytics — bounded per-locus analytics over the substrate (Move #4)"
```

---

## Task 5: Smoke test for the README in-house demo

**Files:**
- Create: `tests/frontdoor_demo.rs`

- [ ] **Step 1: Write the smoke test.** Mirrors the README's bundled-data demo through the CLI, proving the headline commands actually run (so the README cannot rot). Uses the bundled `examples/data/illumina_toy/` fixture; writes intermediates to a unique temp dir.

```rust
//! Proves the README's in-house contract demo actually runs end-to-end on the
//! bundled single-contig fixture: index → sort → plan → variants --enforce → verify.

use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn unique_dir() -> PathBuf {
    static C: AtomicU64 = AtomicU64::new(0);
    let n = C.fetch_add(1, Ordering::Relaxed);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = std::env::temp_dir().join(format!("rosalind-frontdoor-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().expect("spawn rosalind")
}

#[test]
fn readme_inhouse_contract_demo_runs_end_to_end() {
    // CARGO_MANIFEST_DIR points at the crate root; the fixture is bundled there.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fa = root.join("examples/data/illumina_toy/reference.fa");
    let bam = root.join("examples/data/illumina_toy/alignments.bam");
    assert!(fa.exists(), "bundled reference missing: {}", fa.display());
    assert!(bam.exists(), "bundled alignments missing: {}", bam.display());

    let dir = unique_dir();
    let idx = dir.join("toy.idx");
    let sorted = dir.join("toy.sorted.bam");
    let vcf = dir.join("toy.vcf");
    let manifest = dir.join("toy.vcf.manifest.json");

    let out = run(&["index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()]);
    assert!(out.status.success(), "index: {}", String::from_utf8_lossy(&out.stderr));

    let out = run(&["sort", "--input", bam.to_str().unwrap(), "--output", sorted.to_str().unwrap()]);
    assert!(out.status.success(), "sort: {}", String::from_utf8_lossy(&out.stderr));

    let out = run(&["plan", "--index", idx.to_str().unwrap(), "--budget-mb", "512"]);
    assert!(out.status.success(), "plan: {}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stdout).contains("predicted peak"));

    let out = run(&[
        "variants", "--index", idx.to_str().unwrap(),
        "--alignments", sorted.to_str().unwrap(),
        "--memory-budget-mb", "512", "--enforce",
        "-o", vcf.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "variants --enforce: {}", String::from_utf8_lossy(&out.stderr));
    assert!(manifest.exists(), "receipt sidecar must be written");

    let out = run(&["verify", "--manifest", manifest.to_str().unwrap()]);
    assert!(out.status.success(), "verify: {}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stdout).contains("verify: OK"));

    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run it.** (If the bundled `alignments.bam` turns out unsorted/incompatible with the `--index` monotonicity guard, the `sort` step normalizes it; if `variants` still rejects it, regenerate the fixture's BAM via the existing `examples/data/illumina_toy` pipeline and re-commit — but the bundled BAM is expected to work post-`sort`.)

Run: `cd ~/rosalind && cargo test --test frontdoor_demo 2>&1 | tail -15`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add tests/frontdoor_demo.rs && git commit -m "test: smoke-test the README in-house contract demo end-to-end (Move #4)"
```

---

## Task 6: Full verification

**Files:** none (verification only)

- [ ] **Step 1: Format + zero-warning builds**

Run: `cd ~/rosalind && cargo fmt --all && cargo fmt --all -- --check && echo FMT_CLEAN && cargo build 2>&1 | grep -c warning; cargo build --release 2>&1 | grep -c warning`
Expected: `FMT_CLEAN`; `0` and `0`.

- [ ] **Step 2: Doc + example + full suite**

Run: `cd ~/rosalind && cargo test --doc 2>&1 | tail -5 && cargo run --example custom_pileup_analytics >/dev/null 2>&1 && echo EXAMPLE_OK && cargo test 2>&1 | grep -E "FAILED|panicked|[1-9][0-9]* failed" || echo no failures; cargo test 2>&1 | grep -cE "test result: ok\."`
Expected: doctest passes; `EXAMPLE_OK`; `no failures`; `ok.` count ≥ prior + 1 (new `frontdoor_demo` binary).

- [ ] **Step 3: Commit any fmt fixups**

```bash
cd ~/rosalind && git add -A && git commit -m "style: rustfmt fixups (Move #4)" || true
```

---

## Task 7: issue #3 reframe (GitHub — confirm-first, outward-facing)

**Files:** none (GitHub issue edit)

- [ ] **Step 1: Draft the new issue body** (contract-first thesis; √t as the future Phase-D knob; retire "separate out the theory layer"). Present the full draft to the user.
- [ ] **Step 2: On explicit go-ahead only**, post it: `gh issue edit 3 -R logannye/rosalind --body-file <draft>`. Do **not** edit the public issue without confirmation.

---

## Self-Review notes

- **Spec coverage:** §3.1 re-exports → Task 1; §3.2 lib.rs rustdoc + pileup docline → Task 1; §3.3 CONTRACT.md → Task 2; §3.4 README → Task 3; §3.5 demo → **refined** (single-contig in-house demo in README Task 3 Step 4 + smoke test Task 5; the multi-contig generator change is dropped, with the constraint documented up top); §3.6 cookbook → Task 4; §3.7 issue #3 → Task 7 (confirm-first); §4 testing → Tasks 1/4/5/6.
- **Type consistency:** the re-export paths (Task 1) match the verified module exports (`pileup::{Obs,PileupColumn,PileupEngine,PileupParams,ReadSource,SliceSource}`, `io::bam::StreamingBamSource`, `call::{…}`, `core::{MemoryBudget,WorkingSet}`, `genomics::{GenomeIndex,IndexReader,ReferenceView}`, `provenance::RunManifest`). The doctest + example use `rosalind::{PileupEngine,PileupParams,SliceSource}` (now crate-root) + `rosalind::core::{AlignedRead,CigarOp,CigarOpKind,Position,SamFlags}` (the `AlignedRead` 7-field literal matches `src/pileup/engine.rs` test usage). `PileupColumn::{depth(), locus, ref_base}` are its real accessors/fields.
- **No behavior change:** re-exports + docs + a new example + a new test only; no edits to engine/call/cli logic.
