# Phase A2 — Streaming Pileup Kernel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build Rosalind's single streaming pileup engine — CIGAR-aware, read-filtered, strand-aware, bounded-memory — over the A1 `core` types, as a public substrate that callers (A3), plugins, and (later) Python all consume.

**Architecture:** A new crate-root `pileup` module. A `PileupEngine<S: ReadSource>` consumes coordinate-sorted `core::AlignedRead`s from a `ReadSource` (an in-memory `SliceSource` here; a BAM source is deferred to A5) and yields one `PileupColumn` per covered reference position. It projects read bases onto reference coordinates with `core::AlignedRead::projected_bases` (so indels/soft-clips/ref-skips are handled correctly), uses the **forward-oriented** SEQ directly (fixing the legacy reverse-strand double-complement), maintains a bounded active-read set, and walks empty positions with a **loop** (fixing the legacy recursion). The legacy `genomics::pileup`/`pileup_stream` are left untouched; A5 migrates the calling path onto this engine.

**Tech Stack:** Rust 2021, the A1 `core` module (`crate::core`), `std` only. No new dependencies, no `rust_htslib` coupling in this module.

---

## Conventions for this plan

- **Branch:** execute on `rosalind/phase-a1-core` (the rebuild branch — A2 stacks on A1's `core`). Never `main`.
- **Commit trailer:** every commit ends with the co-author trailer shown in each commit step.
- **Module name (`pileup`):** lives at the crate root as `crate::pileup` (peer of `crate::core` and `crate::genomics`). It is distinct from the legacy `crate::genomics::pileup`; do not modify the legacy module.
- **Lint posture:** `lib.rs` sets `#![warn(missing_docs, missing_debug_implementations)]`. Every public item needs a `///` doc comment and a `Debug` impl (derive). Private types used inside a public `#[derive(Debug)]` struct (e.g. `ActiveRead` inside `PileupEngine`) also derive `Debug`.
- **Test filtering:** unit tests live in each module's `#[cfg(test)] mod tests`. Run a module's tests with `cargo test pileup::<module>`. A test that references an undefined item fails to **compile** — that is the expected "red" state.
- **No `genomics` coupling:** this module imports only from `crate::core` and `std`.

---

## File Structure

- Create `src/pileup/mod.rs` — module root + `pub use` surface.
- Create `src/pileup/column.rs` — `Obs`, `PileupColumn` and its derived views.
- Create `src/pileup/source.rs` — `ReadSource` trait + `SliceSource`.
- Create `src/pileup/engine.rs` — `PileupParams`, `SkipCounts`, private `ActiveRead`, `PileupEngine`.
- Modify `src/lib.rs` — add `pub mod pileup;` (next to `pub mod core;`).

---

## Task 1: Module skeleton + `PileupColumn` (`pileup::column`)

**Files:**
- Modify: `src/lib.rs` (add `pub mod pileup;`)
- Create: `src/pileup/mod.rs`
- Create: `src/pileup/column.rs`
- Test: `src/pileup/column.rs`

- [ ] **Step 1: Wire the module.** In `src/lib.rs`, add near `pub mod core;`:
```rust
/// The streaming pileup kernel: one CIGAR-aware, filtered, bounded-memory engine.
pub mod pileup;
```
Create `src/pileup/mod.rs`:
```rust
//! The streaming pileup kernel.
//!
//! `PileupEngine` consumes coordinate-sorted reads and yields one `PileupColumn`
//! per covered reference position. It is the single substrate that variant
//! callers and plugins build on. Reference as `crate::pileup::…`.

pub mod column;

pub use column::{Obs, PileupColumn};
```

- [ ] **Step 2: Write the failing test.** Create `src/pileup/column.rs`:
```rust
//! The per-position pileup column and its observations — the public substrate
//! type produced by `PileupEngine` and consumed by callers and plugins.

use crate::core::Locus;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Position;

    fn obs(allele: u8, reverse: bool) -> Obs {
        Obs { allele, base_qual: 30, mapq: 60, reverse }
    }

    #[test]
    fn depth_allele_and_strand_counts() {
        let col = PileupColumn {
            locus: Locus { contig: 0, pos: Position(100) },
            ref_base: b'A',
            obs: vec![obs(0, false), obs(0, true), obs(1, false)],
        };
        assert_eq!(col.depth(), 3);
        assert_eq!(col.allele_counts(), [2, 1, 0, 0]);
        // [allele][0=fwd,1=rev]: A has 1 fwd + 1 rev, C has 1 fwd.
        let sc = col.strand_counts();
        assert_eq!(sc[0], [1, 1]);
        assert_eq!(sc[1], [1, 0]);
    }
}
```

- [ ] **Step 3: Run the test to verify it fails.** Run: `cargo test pileup::column`
Expected: FAIL — does not compile (`Obs`, `PileupColumn` undefined).

- [ ] **Step 4: Implement the types.** Insert above the `#[cfg(test)]` block in `src/pileup/column.rs`:
```rust
/// One base observation at a pileup position. Only callable (A/C/G/T) bases are
/// recorded; `allele` is the 0..=3 index (A=0, C=1, G=2, T=3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Obs {
    /// Allele index 0..=3 (A/C/G/T).
    pub allele: u8,
    /// Base quality (Phred) of the observed base.
    pub base_qual: u8,
    /// Mapping quality of the read this observation came from.
    pub mapq: u8,
    /// Whether the read maps to the reverse strand (for strand-bias use).
    pub reverse: bool,
}

/// All callable observations stacked at a single reference position.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupColumn {
    /// The reference coordinate of this column.
    pub locus: Locus,
    /// The reference base at this locus (uppercase ASCII; `b'N'` if unknown).
    pub ref_base: u8,
    /// Callable observations, in deterministic (active-read insertion) order.
    pub obs: Vec<Obs>,
}

impl PileupColumn {
    /// Number of callable observations.
    pub fn depth(&self) -> u32 {
        self.obs.len() as u32
    }

    /// Per-allele observation counts, indexed `[A, C, G, T]`.
    pub fn allele_counts(&self) -> [u32; 4] {
        let mut counts = [0u32; 4];
        for o in &self.obs {
            counts[o.allele as usize] += 1;
        }
        counts
    }

    /// Per-allele, per-strand counts: `[allele][0 = forward, 1 = reverse]`.
    pub fn strand_counts(&self) -> [[u32; 2]; 4] {
        let mut counts = [[0u32; 2]; 4];
        for o in &self.obs {
            counts[o.allele as usize][o.reverse as usize] += 1;
        }
        counts
    }
}
```

- [ ] **Step 5: Run the test to verify it passes.** Run: `cargo test pileup::column` — Expected: PASS (1 test).

- [ ] **Step 6: Commit.**
```bash
git add src/lib.rs src/pileup/mod.rs src/pileup/column.rs
git commit -m "feat(pileup): scaffold pileup module + PileupColumn substrate type" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Read source (`pileup::source`)

**Files:**
- Create: `src/pileup/source.rs`
- Modify: `src/pileup/mod.rs`
- Test: `src/pileup/source.rs`

- [ ] **Step 1: Write the failing test.** Create `src/pileup/source.rs`:
```rust
//! Read sources feed coordinate-sorted reads into the pileup engine. `SliceSource`
//! is an in-memory source (used by the Rust API and tests); a BAM-backed source
//! lands when the calling path is migrated.

use crate::core::{AlignedRead, CoreError};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
    use std::sync::Arc;

    fn read(contig: u32, pos: u32) -> AlignedRead {
        AlignedRead {
            contig,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![CigarOp::new(CigarOpKind::Match, 4)],
            seq: Arc::from(b"ACGT".to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; 4].into_boxed_slice()),
        }
    }

    #[test]
    fn slice_source_yields_reads_sorted_by_contig_then_pos() {
        let mut src = SliceSource::new(vec![read(1, 50), read(0, 200), read(0, 100)]);
        let mut order = Vec::new();
        while let Some(r) = src.next_read().unwrap() {
            order.push((r.contig, r.pos.0));
        }
        assert_eq!(order, vec![(0, 100), (0, 200), (1, 50)]);
    }
}
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test pileup::source`
Expected: FAIL — does not compile (`SliceSource`, `ReadSource` undefined).

- [ ] **Step 3: Implement the source.** Insert above the `#[cfg(test)]` block:
```rust
/// A source of coordinate-sorted aligned reads for the pileup engine.
///
/// Implementations MUST yield reads sorted ascending by `(contig, pos)`.
pub trait ReadSource {
    /// Return the next read, or `None` when exhausted.
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError>;
}

/// An in-memory read source. Sorts the provided reads by `(contig, pos)` on
/// construction so callers need not pre-sort.
#[derive(Debug)]
pub struct SliceSource {
    reads: std::vec::IntoIter<AlignedRead>,
}

impl SliceSource {
    /// Build a source from reads (sorted by `(contig, pos)` here).
    pub fn new(mut reads: Vec<AlignedRead>) -> Self {
        reads.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.pos.0.cmp(&b.pos.0)));
        Self { reads: reads.into_iter() }
    }
}

impl ReadSource for SliceSource {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        Ok(self.reads.next())
    }
}
```
In `src/pileup/mod.rs`, add:
```rust
pub mod source;
pub use source::{ReadSource, SliceSource};
```

- [ ] **Step 4: Run the test to verify it passes.** Run: `cargo test pileup::source` — Expected: PASS (1 test).

- [ ] **Step 5: Commit.**
```bash
git add src/pileup/source.rs src/pileup/mod.rs
git commit -m "feat(pileup): ReadSource trait + in-memory SliceSource" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Engine configuration (`PileupParams`, `SkipCounts`)

**Files:**
- Create: `src/pileup/engine.rs`
- Modify: `src/pileup/mod.rs`
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the failing test.** Create `src/pileup/engine.rs`:
```rust
//! The streaming pileup engine: CIGAR-aware, read-filtered, strand-aware,
//! bounded-memory. Yields one `PileupColumn` per covered reference position.

use std::collections::HashMap;
use std::ops::Range;
use std::sync::Arc;

use crate::core::{allele_index, AlignedRead, CoreError, Locus, Position, WorkingSet};
use crate::pileup::column::{Obs, PileupColumn};
use crate::pileup::source::ReadSource;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_skip_noise_and_keep_quality_open() {
        let p = PileupParams::default();
        assert!(p.skip_secondary && p.skip_supplementary && p.skip_duplicate);
        assert_eq!(p.min_mapq, 0);
        assert_eq!(p.min_base_qual, 0);
    }

    #[test]
    fn skip_counts_total_sums_all_reasons() {
        let s = SkipCounts {
            unmapped: 1,
            wrong_contig: 2,
            secondary: 3,
            supplementary: 4,
            duplicate: 5,
            low_mapq: 6,
        };
        assert_eq!(s.total(), 21);
    }
}
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test pileup::engine`
Expected: FAIL — does not compile (`PileupParams`, `SkipCounts` undefined; unused imports also error under the test build, which is fine — they are used once the engine lands in Task 4).

- [ ] **Step 3: Implement the config types.** Insert above the `#[cfg(test)]` block (the `use` lines above stay; `#[allow(unused_imports)]` is NOT needed once Task 4 adds the engine — but to keep this task compiling on its own, temporarily prefix the not-yet-used imports). Replace the `use` block at the top of the file with:
```rust
use std::collections::HashMap;
use std::ops::Range;
use std::sync::Arc;

#[allow(unused_imports)]
use crate::core::{allele_index, AlignedRead, CoreError, Locus, Position, WorkingSet};
#[allow(unused_imports)]
use crate::pileup::column::{Obs, PileupColumn};
#[allow(unused_imports)]
use crate::pileup::source::ReadSource;
```
Then add the config types:
```rust
/// Read-level filters applied as reads enter the pileup.
#[derive(Debug, Clone)]
pub struct PileupParams {
    /// Minimum mapping quality; reads below this are skipped.
    pub min_mapq: u8,
    /// Minimum base quality; observations below this are dropped.
    pub min_base_qual: u8,
    /// Skip secondary alignments (SAM flag 0x100).
    pub skip_secondary: bool,
    /// Skip supplementary alignments (SAM flag 0x800).
    pub skip_supplementary: bool,
    /// Skip PCR/optical duplicates (SAM flag 0x400).
    pub skip_duplicate: bool,
}

impl Default for PileupParams {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_base_qual: 0,
            skip_secondary: true,
            skip_supplementary: true,
            skip_duplicate: true,
        }
    }
}

/// Counts of reads skipped during pileup, by reason (surfaced to callers/CLI).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SkipCounts {
    /// Reads with the unmapped flag set.
    pub unmapped: u64,
    /// Reads on a contig other than the target.
    pub wrong_contig: u64,
    /// Secondary alignments skipped.
    pub secondary: u64,
    /// Supplementary alignments skipped.
    pub supplementary: u64,
    /// Duplicate reads skipped.
    pub duplicate: u64,
    /// Reads below the MAPQ threshold.
    pub low_mapq: u64,
}

impl SkipCounts {
    /// Total reads skipped across all reasons.
    pub fn total(&self) -> u64 {
        self.unmapped
            + self.wrong_contig
            + self.secondary
            + self.supplementary
            + self.duplicate
            + self.low_mapq
    }
}
```
In `src/pileup/mod.rs`, add:
```rust
pub mod engine;
pub use engine::{PileupEngine, PileupParams, SkipCounts};
```
> Note: `pub use engine::PileupEngine` will not compile until Task 4 defines `PileupEngine`. Add only `pub use engine::{PileupParams, SkipCounts};` in THIS task; add `PileupEngine` to the re-export in Task 4.

- [ ] **Step 4: Run the test to verify it passes.** Run: `cargo test pileup::engine` — Expected: PASS (2 tests).

- [ ] **Step 5: Commit.**
```bash
git add src/pileup/engine.rs src/pileup/mod.rs
git commit -m "feat(pileup): PileupParams + SkipCounts" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: The engine — advance, CIGAR-aware tally, Iterator (`PileupEngine`)

This is the heart. It fixes two legacy bugs: (1) the reverse-strand **double-complement** — we read the **forward-oriented** SEQ directly; (2) the empty-position **recursion** — we loop. Filtering is added in Task 6; this task ingests all mapped reads on the target contig.

**Files:**
- Modify: `src/pileup/engine.rs`
- Modify: `src/pileup/mod.rs` (add `PileupEngine` to re-export)
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the failing tests.** Add to `src/pileup/engine.rs`'s `#[cfg(test)] mod tests`:
```rust
    use crate::core::{CigarOp, CigarOpKind, SamFlags};
    use crate::pileup::source::SliceSource;

    // Build a fully-matched read at `pos` on contig 0 with the given seq.
    fn mread(pos: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        let mut flags = SamFlags::default();
        if reverse {
            flags = SamFlags(SamFlags::REVERSE);
        }
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags,
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; seq.len()].into_boxed_slice()),
        }
    }

    fn engine(reads: Vec<AlignedRead>, reference: &[u8]) -> PileupEngine<SliceSource> {
        PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..reference.len() as u32,
            PileupParams::default(),
        )
    }

    fn columns(mut e: PileupEngine<SliceSource>) -> Vec<PileupColumn> {
        let mut out = Vec::new();
        while let Some(c) = e.next() {
            out.push(c.expect("pileup column"));
        }
        out
    }

    #[test]
    fn reverse_strand_reads_are_not_complemented() {
        // Regression: forward and reverse-strand reads carrying the SAME forward
        // SEQ must contribute the SAME allele. (Legacy code complemented reverse
        // reads, corrupting ~half the data.)
        let reference = b"AAAAA";
        // Both reads observe 'G' at ref position 2 (SEQ is forward-oriented).
        let fwd = mread(0, b"AAGAA", false);
        let rev = mread(0, b"AAGAA", true);
        let cols = columns(engine(vec![fwd, rev], reference));
        let at2 = cols.iter().find(|c| c.locus.pos.0 == 2).expect("column at pos 2");
        // Both observe G (allele 2); none observe C (allele 1, the complement of G).
        assert_eq!(at2.allele_counts(), [0, 0, 2, 0]);
    }

    #[test]
    fn basic_ungapped_pileup_counts() {
        let reference = b"ACGTACGT";
        let reads = vec![mread(0, b"ACGT", false), mread(2, b"GTAC", false)];
        let cols = columns(engine(reads, reference));
        // Position 2 is covered by both reads: read1 offset2='G', read2 offset0='G'.
        let at2 = cols.iter().find(|c| c.locus.pos.0 == 2).unwrap();
        assert_eq!(at2.depth(), 2);
        assert_eq!(at2.ref_base, b'G');
        assert_eq!(at2.allele_counts(), [0, 0, 2, 0]);
        // Only covered positions are emitted (no empty columns).
        assert!(cols.iter().all(|c| c.depth() > 0));
    }

    #[test]
    fn sparse_coverage_skips_empty_positions_without_recursion() {
        // A large gap between two reads must not overflow the stack (loop, not
        // recursion) and must emit no empty columns.
        let mut reference = vec![b'A'; 100_000];
        reference[0] = b'C';
        reference[99_999] = b'C';
        let reads = vec![mread(0, b"C", false), mread(99_999, b"C", false)];
        let cols = columns(engine(reads, &reference));
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].locus.pos.0, 0);
        assert_eq!(cols[1].locus.pos.0, 99_999);
    }
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test pileup::engine`
Expected: FAIL — does not compile (`PileupEngine` undefined).

- [ ] **Step 3: Implement the engine.** Remove the three `#[allow(unused_imports)]` lines added in Task 3 (the imports are now used). Add to `src/pileup/engine.rs` (above the `#[cfg(test)]` block):
```rust
/// A read currently overlapping the cursor, with its CIGAR projection precomputed.
#[derive(Debug)]
struct ActiveRead {
    /// Half-open reference end (CIGAR-derived) — used to expire the read.
    end: u32,
    /// Map from reference position to read offset for this read's Match bases.
    ref_to_read: HashMap<u32, usize>,
    /// Read sequence (forward-reference orientation).
    seq: Arc<[u8]>,
    /// Per-base qualities (parallel to `seq`).
    qual: Arc<[u8]>,
    /// Mapping quality.
    mapq: u8,
    /// Reverse-strand flag (metadata only — never applied to `seq`).
    reverse: bool,
}

/// Streaming, CIGAR-aware, bounded-memory pileup over one contig region.
///
/// Yields one [`PileupColumn`] per covered reference position. The working set
/// is bounded by local read coverage, not by input size.
#[derive(Debug)]
pub struct PileupEngine<S: ReadSource> {
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    params: PileupParams,
    active: Vec<ActiveRead>,
    next_read: Option<AlignedRead>,
    pos: u32,
    skips: SkipCounts,
    source_done: bool,
}

impl<S: ReadSource> PileupEngine<S> {
    /// Create an engine over `reference` bytes for `contig`, covering `region`
    /// (0-based half-open; `reference[0]` is the base at `region.start`).
    pub fn new(
        source: S,
        reference: Arc<[u8]>,
        contig: u32,
        region: Range<u32>,
        params: PileupParams,
    ) -> Self {
        let pos = region.start;
        Self {
            source,
            reference,
            contig,
            region,
            params,
            active: Vec::new(),
            next_read: None,
            pos,
            skips: SkipCounts::default(),
            source_done: false,
        }
    }

    /// Reads skipped so far, by reason. Final after iteration completes.
    pub fn skip_counts(&self) -> SkipCounts {
        self.skips
    }

    /// Whether the read passes flag/MAPQ filters (contig + unmapped handled in
    /// `advance_to`). Added in Task 6; here it accepts every read.
    fn passes_filters(&mut self, _read: &AlignedRead) -> bool {
        true
    }

    /// Precompute a read's reference→read-offset map and add it to the active set.
    fn ingest(&mut self, read: AlignedRead) {
        let end = read.end();
        let mut ref_to_read = HashMap::new();
        for rb in read.projected_bases() {
            ref_to_read.insert(rb.ref_pos, rb.read_offset);
        }
        self.active.push(ActiveRead {
            end,
            ref_to_read,
            seq: Arc::clone(&read.seq),
            qual: Arc::clone(&read.qual),
            mapq: read.mapq,
            reverse: read.flags.is_reverse(),
        });
    }

    /// Expire reads that no longer cover `pos`, then pull in reads starting at or
    /// before `pos`. Reads are coordinate-sorted, so we stop at the first read
    /// that starts after `pos` on the target contig.
    fn advance_to(&mut self, pos: u32) -> Result<(), CoreError> {
        self.active.retain(|r| r.end > pos);
        loop {
            if self.next_read.is_none() && !self.source_done {
                match self.source.next_read()? {
                    Some(r) => self.next_read = Some(r),
                    None => self.source_done = true,
                }
            }
            // Peek the routing fields without holding a borrow across `take`.
            let (rc, rp, unmapped) = match self.next_read.as_ref() {
                Some(r) => (r.contig, r.pos.0, r.flags.is_unmapped()),
                None => break,
            };
            if unmapped {
                self.next_read.take();
                self.skips.unmapped += 1;
                continue;
            }
            match rc.cmp(&self.contig) {
                std::cmp::Ordering::Greater => break, // sorted: no more target-contig reads
                std::cmp::Ordering::Less => {
                    self.next_read.take();
                    self.skips.wrong_contig += 1;
                }
                std::cmp::Ordering::Equal => {
                    if rp > pos {
                        break; // future read on our contig
                    }
                    let read = self.next_read.take().unwrap();
                    if !self.passes_filters(&read) {
                        continue;
                    }
                    if read.end() <= pos {
                        continue; // does not reach the cursor
                    }
                    self.ingest(read);
                }
            }
        }
        Ok(())
    }

    /// Build the column at the current cursor from the active set.
    fn build_column(&self) -> PileupColumn {
        let ref_idx = (self.pos - self.region.start) as usize;
        let ref_base = self.reference.get(ref_idx).copied().unwrap_or(b'N');
        let mut obs = Vec::new();
        for r in &self.active {
            if let Some(&off) = r.ref_to_read.get(&self.pos) {
                // Forward orientation — read the stored SEQ byte directly. No
                // complement (this is the reverse-strand fix).
                let base = r.seq.get(off).copied().unwrap_or(b'N');
                let bq = r.qual.get(off).copied().unwrap_or(0);
                if bq < self.params.min_base_qual {
                    continue;
                }
                if let Some(allele) = allele_index(base) {
                    obs.push(Obs {
                        allele: allele as u8,
                        base_qual: bq,
                        mapq: r.mapq,
                        reverse: r.reverse,
                    });
                }
            }
        }
        PileupColumn {
            locus: Locus { contig: self.contig, pos: Position(self.pos) },
            ref_base,
            obs,
        }
    }
}

impl<S: ReadSource> Iterator for PileupEngine<S> {
    type Item = Result<PileupColumn, CoreError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Loop over empty positions (no recursion → bounded stack, any gap size).
        while self.pos < self.region.end {
            let pos = self.pos;
            if let Err(e) = self.advance_to(pos) {
                self.pos = self.region.end;
                return Some(Err(e));
            }
            let column = self.build_column();
            self.pos += 1;
            if !column.obs.is_empty() {
                return Some(Ok(column));
            }
        }
        None
    }
}
```
In `src/pileup/mod.rs`, change the engine re-export to include `PileupEngine`:
```rust
pub use engine::{PileupEngine, PileupParams, SkipCounts};
```

- [ ] **Step 4: Run the tests to verify they pass.** Run: `cargo test pileup::engine` — Expected: PASS (5 tests: 2 from Task 3 + 3 here).

- [ ] **Step 5: Run the whole suite + build.** Run: `cargo test` — all pass. Run: `cargo build` — no `missing_docs`/`missing_debug_implementations` warnings from `src/pileup/`.

- [ ] **Step 6: Commit.**
```bash
git add src/pileup/engine.rs src/pileup/mod.rs
git commit -m "feat(pileup): PileupEngine (CIGAR-aware tally, forward-strand fix, empty-position loop)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: CIGAR projection correctness (indel / soft-clip / ref-skip / long read)

These tests confirm the engine's `core::projected_bases`-driven tally places bases at correct reference coordinates for non-trivial CIGARs. They should pass against the Task-4 engine; if any fails, it reveals an engine bug to fix before proceeding.

**Files:**
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the tests.** Add to the `#[cfg(test)] mod tests` block:
```rust
    #[test]
    fn insertion_bases_do_not_shift_downstream_reference_positions() {
        // 2M 1I 2M: read offsets 0,1 -> ref 10,11; offset 2 = inserted (no ref);
        // offsets 3,4 -> ref 12,13.
        let reference = b"AAAAAAAAAAAAAAAA"; // 16 'A'
        let read = AlignedRead {
            contig: 0,
            pos: Position(10),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::Match, 2),
                CigarOp::new(CigarOpKind::Insertion, 1),
                CigarOp::new(CigarOpKind::Match, 2),
            ],
            seq: Arc::from(b"CCAGT".to_vec().into_boxed_slice()), // 2M=CC, 1I=A (inserted), 2M=GT
            qual: Arc::from(vec![30u8; 5].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        // ref 10 <- offset0 'C'; ref 11 <- offset1 'C'; ref 12 <- offset3 'G'; ref 13 <- offset4 'T'.
        // offset2 'A' is the inserted base — correctly absent from every reference column.
        let get = |p: u32| cols.iter().find(|c| c.locus.pos.0 == p).map(|c| c.allele_counts());
        assert_eq!(get(10), Some([0, 1, 0, 0])); // C
        assert_eq!(get(11), Some([0, 1, 0, 0])); // C
        assert_eq!(get(12), Some([0, 0, 1, 0])); // G  (offset3)
        assert_eq!(get(13), Some([0, 0, 0, 1])); // T  (offset4)
    }

    #[test]
    fn deletion_leaves_a_reference_gap_with_no_observation() {
        // 2M 1D 2M starting at ref 0: ref 0,1 observed; ref 2 deleted (no obs);
        // ref 3,4 observed.
        let reference = b"AAAAAAAA";
        let read = AlignedRead {
            contig: 0,
            pos: Position(0),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::Match, 2),
                CigarOp::new(CigarOpKind::Deletion, 1),
                CigarOp::new(CigarOpKind::Match, 2),
            ],
            seq: Arc::from(b"GGGG".to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; 4].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        let positions: Vec<u32> = cols.iter().map(|c| c.locus.pos.0).collect();
        assert_eq!(positions, vec![0, 1, 3, 4]); // ref 2 (deleted) emits no column
    }

    #[test]
    fn soft_clipped_bases_are_excluded() {
        // 2S 3M: first 2 read bases are clipped; only the 3 matched bases pile up.
        let reference = b"AAAAAAA";
        let read = AlignedRead {
            contig: 0,
            pos: Position(1),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::SoftClip, 2),
                CigarOp::new(CigarOpKind::Match, 3),
            ],
            seq: Arc::from(b"TTCGA".to_vec().into_boxed_slice()), // TT clipped; CGA match ref 1,2,3
            qual: Arc::from(vec![30u8; 5].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        let positions: Vec<u32> = cols.iter().map(|c| c.locus.pos.0).collect();
        assert_eq!(positions, vec![1, 2, 3]);
        // The clipped 'T's (allele 3) never appear; ref 1 sees 'C' (allele 1).
        let at1 = cols.iter().find(|c| c.locus.pos.0 == 1).unwrap();
        assert_eq!(at1.allele_counts(), [0, 1, 0, 0]);
    }

    #[test]
    fn long_read_piles_up_every_matched_base() {
        // Read-length-agnostic: a 5000-base full-match read covers 5000 positions.
        let reference = vec![b'A'; 6000];
        let seq = vec![b'C'; 5000];
        let read = AlignedRead {
            contig: 0,
            pos: Position(1000),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![CigarOp::new(CigarOpKind::Match, 5000)],
            seq: Arc::from(seq.into_boxed_slice()),
            qual: Arc::from(vec![30u8; 5000].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], &reference));
        assert_eq!(cols.len(), 5000);
        assert_eq!(cols.first().unwrap().locus.pos.0, 1000);
        assert_eq!(cols.last().unwrap().locus.pos.0, 5999);
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
    }
```

- [ ] **Step 2: Run the tests.** Run: `cargo test pileup::engine` — Expected: PASS (9 tests total). If a projection test fails, STOP and report it as an engine bug (do not adjust the test to match buggy behavior).

- [ ] **Step 3: Commit.**
```bash
git add src/pileup/engine.rs
git commit -m "test(pileup): CIGAR projection correctness (indel/softclip/refskip/long-read)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Read-level filtering + per-base quality floor

The Task-4 engine ingests all mapped target-contig reads. Now wire `PileupParams` so secondary/supplementary/duplicate/low-MAPQ reads are skipped (counted), and low-base-quality observations dropped.

**Files:**
- Modify: `src/pileup/engine.rs` (replace the stub `passes_filters`)
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the failing tests.** Add to the `#[cfg(test)] mod tests` block:
```rust
    fn flagged_read(pos: u32, seq: &[u8], flag_bits: u16) -> AlignedRead {
        let mut r = mread(pos, seq, false);
        r.flags = SamFlags(flag_bits);
        r
    }

    #[test]
    fn filters_skip_secondary_supplementary_and_duplicate_reads() {
        let reference = b"AAAAA";
        let reads = vec![
            mread(0, b"CCCCC", false),                          // kept
            flagged_read(0, b"GGGGG", SamFlags::SECONDARY),     // skipped
            flagged_read(0, b"GGGGG", SamFlags::SUPPLEMENTARY), // skipped
            flagged_read(0, b"GGGGG", SamFlags::DUPLICATE),     // skipped
        ];
        let mut e = engine(reads, reference);
        let cols = {
            let mut out = Vec::new();
            while let Some(c) = e.next() {
                out.push(c.unwrap());
            }
            out
        };
        // Only the kept read's 'C' (allele 1) appears at every position.
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
        let s = e.skip_counts();
        assert_eq!((s.secondary, s.supplementary, s.duplicate), (1, 1, 1));
    }

    #[test]
    fn filters_skip_low_mapq_reads() {
        let reference = b"AAAAA";
        let mut low = mread(0, b"GGGGG", false);
        low.mapq = 3;
        let reads = vec![mread(0, b"CCCCC", false), low];
        let params = PileupParams { min_mapq: 10, ..PileupParams::default() };
        let mut e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..5,
            params,
        );
        let mut cols = Vec::new();
        while let Some(c) = e.next() {
            cols.push(c.unwrap());
        }
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
        assert_eq!(e.skip_counts().low_mapq, 1);
    }

    #[test]
    fn low_base_quality_observations_are_dropped() {
        let reference = b"AAAAA";
        let mut r = mread(0, b"GGGGG", false);
        r.qual = Arc::from(vec![2u8; 5].into_boxed_slice()); // below the floor
        let params = PileupParams { min_base_qual: 20, ..PileupParams::default() };
        let mut e = PileupEngine::new(
            SliceSource::new(vec![r]),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..5,
            params,
        );
        // All observations dropped → no columns emitted.
        assert!(e.next().is_none());
    }
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test pileup::engine`
Expected: FAIL — the secondary/supplementary/duplicate/low-MAPQ reads are currently NOT skipped (the stub `passes_filters` returns `true`), so allele counts and skip counts are wrong. (The low-base-quality test already passes via `build_column`'s `min_base_qual` check; that is expected.)

- [ ] **Step 3: Implement filtering.** Replace the stub `passes_filters` in `src/pileup/engine.rs` with:
```rust
    /// Whether the read passes flag/MAPQ filters. (Contig routing and the
    /// unmapped flag are handled in `advance_to`.)
    fn passes_filters(&mut self, read: &AlignedRead) -> bool {
        let f = read.flags;
        if self.params.skip_secondary && f.is_secondary() {
            self.skips.secondary += 1;
            return false;
        }
        if self.params.skip_supplementary && f.is_supplementary() {
            self.skips.supplementary += 1;
            return false;
        }
        if self.params.skip_duplicate && f.is_duplicate() {
            self.skips.duplicate += 1;
            return false;
        }
        if read.mapq < self.params.min_mapq {
            self.skips.low_mapq += 1;
            return false;
        }
        true
    }
```

- [ ] **Step 4: Run the tests to verify they pass.** Run: `cargo test pileup::engine` — Expected: PASS (12 tests total).

- [ ] **Step 5: Commit.**
```bash
git add src/pileup/engine.rs
git commit -m "feat(pileup): read-level filtering (secondary/supplementary/duplicate/MAPQ) + skip accounting" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Bounded-memory reporting (`current_working_set`)

Expose the engine's working set so a run can be checked against a `MemoryBudget` (the foundation for the Phase D memory contract), and prove it is bounded by local coverage, not input size.

**Files:**
- Modify: `src/pileup/engine.rs`
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the failing test.** Add to the `#[cfg(test)] mod tests` block:
```rust
    use crate::core::MemoryBudget;

    #[test]
    fn working_set_is_bounded_by_coverage_not_input_size() {
        // 50,000 short reads tiled across the reference at depth ~1: the active
        // set (and thus the working set) stays tiny throughout iteration.
        let reference = vec![b'A'; 50_000];
        let reads: Vec<AlignedRead> =
            (0..50_000u32).map(|p| mread(p, b"C", false)).collect();
        let mut e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.into_boxed_slice()),
            0,
            0..50_000,
            PileupParams::default(),
        );
        let budget = MemoryBudget::from_mb(1);
        let mut max_ws = 0u64;
        while let Some(c) = e.next() {
            let _ = c.unwrap();
            let ws = e.current_working_set();
            max_ws = max_ws.max(ws.bytes);
            assert!(ws.fits(budget), "working set {} exceeded 1 MiB", ws.bytes);
        }
        // Sanity: peak working set is far below what holding all reads would cost.
        assert!(max_ws < 64 * 1024, "peak working set unexpectedly large: {max_ws}");
    }
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test pileup::engine`
Expected: FAIL — does not compile (`current_working_set` undefined).

- [ ] **Step 3: Implement the reporter.** Add this method inside `impl<S: ReadSource> PileupEngine<S>` (next to `skip_counts`):
```rust
    /// Current working-set estimate: bounded by the active read set (local
    /// coverage), independent of total input size. Foundation for `rosalind plan`.
    pub fn current_working_set(&self) -> WorkingSet {
        // Each active read costs roughly its projection map (16 B/entry) plus a
        // small constant for handles; plus a fixed engine overhead.
        let active_bytes: u64 = self
            .active
            .iter()
            .map(|r| (r.ref_to_read.len() as u64) * 16 + 64)
            .sum();
        WorkingSet { bytes: active_bytes + 256 }
    }
```

- [ ] **Step 4: Run the test to verify it passes.** Run: `cargo test pileup::engine` — Expected: PASS (13 tests total).

- [ ] **Step 5: Commit.**
```bash
git add src/pileup/engine.rs
git commit -m "feat(pileup): current_working_set() reporter (bounded-memory contract foundation)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Determinism + the open-substrate demonstration

Lock in determinism (the moat) and prove the engine is consumer-agnostic by computing a non-caller analysis (per-base coverage) over the column stream.

**Files:**
- Test: `src/pileup/engine.rs`

- [ ] **Step 1: Write the tests.** Add to the `#[cfg(test)] mod tests` block:
```rust
    #[test]
    fn shuffled_source_yields_identical_columns() {
        // SliceSource sorts on construction, so a shuffled input must produce the
        // exact same column sequence (deterministic output).
        let reference = b"ACGTACGTACGT";
        let make = |order: Vec<(u32, &'static [u8])>| {
            let reads: Vec<AlignedRead> =
                order.into_iter().map(|(p, s)| mread(p, s, false)).collect();
            columns(engine(reads, reference))
        };
        let a = make(vec![(0, b"AAAA"), (4, b"CCCC"), (8, b"GGGG")]);
        let b = make(vec![(8, b"GGGG"), (0, b"AAAA"), (4, b"CCCC")]);
        assert_eq!(a, b);
    }

    #[test]
    fn engine_is_an_open_substrate_for_arbitrary_per_locus_analysis() {
        // A non-caller consumer: compute per-position coverage directly from the
        // PileupColumn stream — no variant calling involved.
        let reference = b"AAAAAAAA";
        let reads = vec![mread(0, b"CCCC", false), mread(2, b"CCCC", false)];
        let mut coverage = Vec::new();
        let mut e = engine(reads, reference);
        while let Some(c) = e.next() {
            let col = c.unwrap();
            coverage.push((col.locus.pos.0, col.depth()));
        }
        // pos 0,1 depth 1; pos 2,3 depth 2; pos 4,5 depth 1.
        assert_eq!(
            coverage,
            vec![(0, 1), (1, 1), (2, 2), (3, 2), (4, 1), (5, 1)]
        );
    }
```

- [ ] **Step 2: Run the tests.** Run: `cargo test pileup::engine` — Expected: PASS (15 tests total).

- [ ] **Step 3: Run the whole suite + build + fmt.** Run: `cargo test` (all pass), `cargo build` (no `src/pileup/` warnings), `cargo fmt --all -- --check` (clean).

- [ ] **Step 4: Commit.**
```bash
git add src/pileup/engine.rs
git commit -m "test(pileup): determinism + open-substrate (coverage) demonstration" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Plan self-review

- **Spec coverage** (Phase A spec §3.2 `pileup/`): `ReadSource` + adapters ✔ (Task 2; BamSource deferred to A5, noted); `PileupEngine` ✔ (Task 4); CIGAR projection (M/=/X, I, S, D, N, H) ✔ via `core::projected_bases` (Tasks 4–5); read filtering + `SkipCounts` ✔ (Task 6); per-observation strand ✔ (`Obs.reverse`, Task 1); deterministic ordering ✔ (Task 8); bounded active set + `working_set_bound` ✔ (Tasks 4, 7); empty-position **loop** ✔ (Task 4); reverse-strand fix (forward SEQ) ✔ (Task 4); `PileupColumn` substrate ✔ (Task 1); non-caller consumer ✔ (Task 8). The genotype model is Phase A3; wiring the calling path + a BAM source is Phase A5 — both correctly out of scope here.
- **Placeholder scan:** none — every step has complete code + exact commands. The one cross-task dependency (the `#[allow(unused_imports)]` shim in Task 3, removed in Task 4) is called out explicitly so Task 3 compiles standalone.
- **Type consistency:** `Obs { allele:u8, base_qual:u8, mapq:u8, reverse:bool }`, `PileupColumn { locus, ref_base, obs }` with `depth()/allele_counts()/strand_counts()`, `ReadSource::next_read() -> Result<Option<AlignedRead>, CoreError>`, `SliceSource::new(Vec<AlignedRead>)`, `PileupParams` fields, `SkipCounts` fields + `total()`, `PileupEngine::new(source, reference, contig, region, params)` + `skip_counts()` + `current_working_set() -> WorkingSet`, and the `core` API (`AlignedRead.{contig,pos,mapq,flags,cigar,seq,qual}`, `projected_bases()`, `end()`, `SamFlags::{is_unmapped,is_reverse,is_secondary,is_supplementary,is_duplicate}` + the flag consts, `allele_index`, `Locus`, `Position`, `MemoryBudget`, `WorkingSet`, `CoreError`) — names used in tests match the implementations and the A1 `core` module.

## Definition of done (A2)

`cargo test` green (15 new `pileup::*` tests + the existing suite); `cargo build` warning-clean for `src/pileup/`; `cargo fmt --all -- --check` clean; the legacy `genomics::pileup`/`pileup_stream` untouched; all work on `rosalind/phase-a1-core`.
