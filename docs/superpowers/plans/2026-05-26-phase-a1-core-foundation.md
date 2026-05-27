# Phase A1 — `core` Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create Rosalind's shared core type layer — the contig-aware, 64-bit-safe coordinate model, the canonical CIGAR-aware aligned-read record, the base/allele model, the memory-budget model, and the typed error taxonomy — that every later layer (pileup, call, io, index, align) builds on.

**Architecture:** A new crate-root module `core` with five focused submodules (`error`, `locus`, `sequence`, `record`, `budget`). Self-contained and fully unit-tested; **no existing consumers are migrated in A1** (the legacy `genomics::types` continues to compile unchanged). Later plans (A2, A5) migrate consumers onto these types and delete the legacy ones.

**Tech Stack:** Rust 2021, `thiserror` (already a dependency), `std`. No new dependencies.

---

## Conventions for this plan

- **Branch, never main.** Before Task 1, create a feature branch or worktree (repo policy: branch + PR, never commit to `main`). Example: `git switch -c rosalind/phase-a1-core`.
- **Commit trailer.** Every commit ends with the co-author trailer, shown in each commit step.
- **Module name note (`core`).** This module is named `core` to match the architecture doc. A crate-root `mod core` shadows the std `core` crate for *bare* `core::…` paths within this crate. Rule for contributors: inside the crate, reference these types as `crate::core::…`, and reach the std core crate (rarely needed) as `::core::…`. The codebase uses `std::…`, so this is a non-issue in practice.
- **Lint posture.** `lib.rs` sets `#![warn(missing_docs, missing_debug_implementations)]`. Every public item below has a doc comment and derives `Debug`.
- **Test filtering.** Unit tests live in each module's `#[cfg(test)] mod tests`. Run a module's tests with `cargo test core::<module>`. Before an implementation exists, the test will fail to **compile** (type not yet defined) — that is the expected "red" state for Rust TDD.

---

## File Structure

- Create `src/core/mod.rs` — module root + re-exports (`pub use …`). One responsibility: assemble the core layer's public surface.
- Create `src/core/error.rs` — `CoreError` typed taxonomy.
- Create `src/core/locus.rs` — `Contig`, `Position`, `Locus`, `ContigSet` (per-contig `u32`, 64-bit-safe global offsets).
- Create `src/core/sequence.rs` — `BaseCode`, `allele_index`, IUPAC→N folding.
- Create `src/core/record.rs` — `CigarOpKind` (incl. `RefSkip`), `CigarOp`, `SamFlags`, `AlignedRead`, `RefBase`, CIGAR projection.
- Create `src/core/budget.rs` — `MemoryBudget`, `WorkingSet`.
- Modify `src/lib.rs` — add `pub mod core;` and a doc line.

---

## Task 1: Module skeleton + typed errors (`core::error`)

**Files:**
- Modify: `src/lib.rs` (add `pub mod core;`)
- Create: `src/core/mod.rs`
- Create: `src/core/error.rs`
- Test: `src/core/error.rs` (`#[cfg(test)] mod tests`)

- [ ] **Step 1: Wire the module into the crate**

In `src/lib.rs`, add to the module list (next to the other `pub mod` lines):

```rust
/// Core types: the lingua franca shared by every layer (io, index, align, pileup, call).
pub mod core;
```

Create `src/core/mod.rs`:

```rust
//! Core types — the lingua franca shared by every Rosalind layer.
//!
//! Note: this module is named `core`; inside the crate always reference it as
//! `crate::core::…`. Reach the std `core` crate (rarely needed) as `::core::…`.

pub mod error;

pub use error::CoreError;
// Later tasks add `pub mod locus/budget/sequence/record;` + their re-exports as
// each file is created; declaring them before the files exist would not compile.
```

- [ ] **Step 2: Write the failing test**

Create `src/core/error.rs`:

```rust
//! Typed error taxonomy for the Rosalind core/library layer.
//!
//! The CLI boundary maps these into `anyhow`; library code returns `CoreError`.

use thiserror::Error;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_exceeded_displays_both_numbers() {
        let e = CoreError::BudgetExceeded { needed: 4096, budget: 1024 };
        let msg = e.to_string();
        assert!(msg.contains("4096"), "message should report needed bytes: {msg}");
        assert!(msg.contains("1024"), "message should report budget bytes: {msg}");
    }
}
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cargo test core::error`
Expected: FAIL — does not compile (`cannot find type CoreError in this scope`).

- [ ] **Step 4: Implement `CoreError`**

Insert above the `#[cfg(test)]` block in `src/core/error.rs`:

```rust
/// Errors produced by the Rosalind core/library layer.
#[derive(Debug, Error)]
pub enum CoreError {
    /// A streaming stage's working set would exceed the declared budget.
    #[error("memory budget exceeded: working set {needed} bytes > budget {budget} bytes")]
    BudgetExceeded {
        /// Bytes the stage would require.
        needed: u64,
        /// Bytes the caller permitted.
        budget: u64,
    },
    /// An input record could not be interpreted.
    #[error("malformed record: {0}")]
    MalformedRecord(String),
    /// A contig id was not present in the active `ContigSet`.
    #[error("invalid contig id {0}")]
    InvalidContig(u32),
    /// An underlying I/O failure.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cargo test core::error`
Expected: PASS (1 test).

- [ ] **Step 6: Commit**

```bash
git add src/lib.rs src/core/mod.rs src/core/error.rs
git commit -m "feat(core): scaffold core module + CoreError taxonomy" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Coordinate model (`core::locus`)

**Files:**
- Create: `src/core/locus.rs`
- Modify: `src/core/mod.rs` (re-exports)
- Test: `src/core/locus.rs`

- [ ] **Step 1: Write the failing test**

Create `src/core/locus.rs`:

```rust
//! Contig-aware genomic coordinates. Per-contig positions are `u32`; the
//! concatenated-genome offset used by the multi-contig FM-index (Phase B) is
//! 64-bit-safe, so no single coordinate space is capped at 4.29 Gbp.

use std::sync::Arc;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_set_assigns_ids_and_global_offsets() {
        let mut set = ContigSet::new();
        let chr1 = set.push("chr1", 1_000);
        let chr2 = set.push("chr2", 500);

        assert_eq!(chr1, 0);
        assert_eq!(chr2, 1);
        assert_eq!(set.by_id(chr1).unwrap().global_offset, 0);
        assert_eq!(set.by_id(chr2).unwrap().global_offset, 1_000);
        assert_eq!(set.by_name("chr2").unwrap().id, 1);
        assert_eq!(set.total_length(), 1_500);
    }

    #[test]
    fn loci_order_by_contig_then_position() {
        let a = Locus { contig: 0, pos: Position(100) };
        let b = Locus { contig: 0, pos: Position(200) };
        let c = Locus { contig: 1, pos: Position(0) };
        assert!(a < b && b < c);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test core::locus`
Expected: FAIL — does not compile (`ContigSet`, `Locus`, `Position` undefined).

- [ ] **Step 3: Implement the coordinate types**

Insert above the `#[cfg(test)]` block in `src/core/locus.rs`:

```rust
/// A reference contig (chromosome / sequence) with a stable id.
#[derive(Debug, Clone)]
pub struct Contig {
    /// Dense 0-based id (index into the `ContigSet`).
    pub id: u32,
    /// Contig name as it appears in the reference / SAM `@SQ`.
    pub name: Arc<str>,
    /// Length in bases.
    pub length: u32,
    /// 0-based offset of this contig within the concatenated genome
    /// (64-bit-safe; used by the multi-contig index in Phase B).
    pub global_offset: u64,
}

/// A 0-based position within a single contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Position(pub u32);

/// A canonical genomic coordinate: a contig id plus a position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Locus {
    /// Contig id (index into the originating `ContigSet`).
    pub contig: u32,
    /// 0-based position within the contig.
    pub pos: Position,
}

/// The ordered set of reference contigs: name<->id lookup and global offsets.
#[derive(Debug, Clone, Default)]
pub struct ContigSet {
    contigs: Vec<Contig>,
}

impl ContigSet {
    /// Create an empty contig set.
    pub fn new() -> Self {
        Self { contigs: Vec::new() }
    }

    /// Append a contig, assigning the next id and the running global offset.
    /// Returns the new contig's id.
    pub fn push(&mut self, name: impl Into<Arc<str>>, length: u32) -> u32 {
        let id = self.contigs.len() as u32;
        let global_offset = self
            .contigs
            .last()
            .map_or(0, |c| c.global_offset + c.length as u64);
        self.contigs.push(Contig { id, name: name.into(), length, global_offset });
        id
    }

    /// Look up a contig by id.
    pub fn by_id(&self, id: u32) -> Option<&Contig> {
        self.contigs.get(id as usize)
    }

    /// Look up a contig by name.
    pub fn by_name(&self, name: &str) -> Option<&Contig> {
        self.contigs.iter().find(|c| c.name.as_ref() == name)
    }

    /// Number of contigs.
    pub fn len(&self) -> usize {
        self.contigs.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }

    /// Iterate contigs in id order.
    pub fn iter(&self) -> impl Iterator<Item = &Contig> {
        self.contigs.iter()
    }

    /// Total length of the concatenated genome (64-bit-safe).
    pub fn total_length(&self) -> u64 {
        self.contigs
            .last()
            .map_or(0, |c| c.global_offset + c.length as u64)
    }
}
```

In `src/core/mod.rs`, add the module declaration and its re-exports:

```rust
pub mod locus;
pub use locus::{Contig, ContigSet, Locus, Position};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test core::locus`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/core/locus.rs src/core/mod.rs
git commit -m "feat(core): contig-aware, 64-bit-safe coordinate model" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Memory budget model (`core::budget`)

**Files:**
- Create: `src/core/budget.rs`
- Modify: `src/core/mod.rs`
- Test: `src/core/budget.rs`

- [ ] **Step 1: Write the failing test**

Create `src/core/budget.rs`:

```rust
//! Memory-as-a-contract primitives. Streaming stages report a `WorkingSet`
//! bound so a run can be checked against a `MemoryBudget` *before* it starts
//! (the foundation for `rosalind plan`).

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_admits_within_and_rejects_beyond() {
        let budget = MemoryBudget::from_mb(2);
        assert!(budget.admits(1_000_000));
        assert!(!budget.admits(3 * 1024 * 1024));
    }

    #[test]
    fn unlimited_admits_everything_and_working_set_checks_fit() {
        assert!(MemoryBudget::unlimited().admits(u64::MAX));
        let ws = WorkingSet { bytes: 512 };
        assert!(ws.fits(MemoryBudget::from_mb(1)));
        assert!(!WorkingSet { bytes: 5 * 1024 * 1024 }.fits(MemoryBudget::from_mb(1)));
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test core::budget`
Expected: FAIL — does not compile (`MemoryBudget`, `WorkingSet` undefined).

- [ ] **Step 3: Implement the budget types**

Insert above the `#[cfg(test)]` block:

```rust
/// A declared cap on a streaming stage's working set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryBudget {
    /// Maximum permitted working-set size, in bytes.
    pub bytes: u64,
}

impl MemoryBudget {
    /// A budget expressed in mebibytes.
    pub fn from_mb(mb: u64) -> Self {
        Self { bytes: mb.saturating_mul(1024 * 1024) }
    }

    /// An effectively unbounded budget.
    pub fn unlimited() -> Self {
        Self { bytes: u64::MAX }
    }

    /// Whether a working set of `working_set_bytes` is permitted.
    pub fn admits(self, working_set_bytes: u64) -> bool {
        working_set_bytes <= self.bytes
    }
}

/// A reported working-set bound for a streaming stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkingSet {
    /// Estimated peak working-set size, in bytes.
    pub bytes: u64,
}

impl WorkingSet {
    /// Whether this working set fits within `budget`.
    pub fn fits(self, budget: MemoryBudget) -> bool {
        budget.admits(self.bytes)
    }
}
```

In `src/core/mod.rs`, add:

```rust
pub mod budget;
pub use budget::{MemoryBudget, WorkingSet};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test core::budget`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/core/budget.rs src/core/mod.rs
git commit -m "feat(core): MemoryBudget/WorkingSet (memory-as-a-contract)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Base / allele model (`core::sequence`)

**Files:**
- Create: `src/core/sequence.rs`
- Modify: `src/core/mod.rs`
- Test: `src/core/sequence.rs`

- [ ] **Step 1: Write the failing test**

Create `src/core/sequence.rs`:

```rust
//! Canonical DNA base model. Exact bases are A/C/G/T (U folds to T); every
//! other byte — including IUPAC ambiguity codes — folds to `N` (lossy but
//! explicit). `N` is not a callable allele.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_bases_map_directly_and_u_folds_to_t() {
        assert_eq!(BaseCode::from_ascii_lossy(b'a'), BaseCode::A);
        assert_eq!(BaseCode::from_ascii_lossy(b'C'), BaseCode::C);
        assert_eq!(BaseCode::from_ascii_lossy(b'u'), BaseCode::T);
        assert_eq!(BaseCode::from_ascii_lossy(b'T'), BaseCode::T);
    }

    #[test]
    fn ambiguous_and_unknown_fold_to_n_and_are_not_exact() {
        assert_eq!(BaseCode::from_ascii_lossy(b'R'), BaseCode::N); // IUPAC purine
        assert_eq!(BaseCode::from_ascii_lossy(b'.'), BaseCode::N);
        assert!(!BaseCode::is_exact(b'R'));
        assert!(BaseCode::is_exact(b'g'));
    }

    #[test]
    fn allele_index_is_some_for_acgt_and_none_for_n() {
        assert_eq!(allele_index(b'A'), Some(0));
        assert_eq!(allele_index(b'T'), Some(3));
        assert_eq!(allele_index(b'N'), None);
        assert_eq!(allele_index(b'R'), None);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test core::sequence`
Expected: FAIL — does not compile (`BaseCode`, `allele_index` undefined).

- [ ] **Step 3: Implement the base model**

Insert above the `#[cfg(test)]` block:

```rust
/// A canonical DNA base. `N` represents any non-ACGT / ambiguous base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCode {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine (also the fold target for uracil `U`).
    T,
    /// Any non-ACGT / ambiguous base.
    N,
}

impl BaseCode {
    /// Allele index 0..=3 for A/C/G/T; `None` for `N` (not a callable allele).
    pub fn allele_index(self) -> Option<usize> {
        match self {
            BaseCode::A => Some(0),
            BaseCode::C => Some(1),
            BaseCode::G => Some(2),
            BaseCode::T => Some(3),
            BaseCode::N => None,
        }
    }

    /// Canonical uppercase ASCII byte for this base.
    pub fn to_ascii(self) -> u8 {
        match self {
            BaseCode::A => b'A',
            BaseCode::C => b'C',
            BaseCode::G => b'G',
            BaseCode::T => b'T',
            BaseCode::N => b'N',
        }
    }

    /// Map an ASCII byte to a `BaseCode`. A/C/G/T (any case) map directly,
    /// `U`/`u` fold to `T`, and every other byte folds to `N`.
    pub fn from_ascii_lossy(b: u8) -> BaseCode {
        match b.to_ascii_uppercase() {
            b'A' => BaseCode::A,
            b'C' => BaseCode::C,
            b'G' => BaseCode::G,
            b'T' | b'U' => BaseCode::T,
            _ => BaseCode::N,
        }
    }

    /// Whether `b` is a recognized exact base (A/C/G/T/U). Bytes that fold to
    /// `N` return `false`; callers use this to count ambiguity warnings.
    pub fn is_exact(b: u8) -> bool {
        matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U')
    }
}

/// Convenience: allele index 0..=3 for an ASCII base, or `None` for N/other.
pub fn allele_index(base: u8) -> Option<usize> {
    BaseCode::from_ascii_lossy(base).allele_index()
}
```

In `src/core/mod.rs`, add:

```rust
pub mod sequence;
pub use sequence::{allele_index, BaseCode};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test core::sequence`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/core/sequence.rs src/core/mod.rs
git commit -m "feat(core): BaseCode/allele model with explicit IUPAC->N folding" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Canonical record + CIGAR projection (`core::record`)

This is the most important task in A1: the CIGAR-aware projection here is the
single source of truth the pileup kernel (A2) uses, and it is what makes the
engine correct on indels, soft-clips, ref-skips, and **long reads**.

**Files:**
- Create: `src/core/record.rs`
- Modify: `src/core/mod.rs`
- Test: `src/core/record.rs`

- [ ] **Step 1: Write the failing tests**

Create `src/core/record.rs`:

```rust
//! The canonical aligned-read record and CIGAR projection.
//!
//! SEQ is stored forward-reference-oriented (as in SAM/BAM); strand is
//! metadata and is never applied to the bytes. Records are read-length-
//! agnostic — short Illumina reads and long Nanopore/PacBio reads are handled
//! identically.

use std::sync::Arc;

use crate::core::locus::Position;

#[cfg(test)]
mod tests {
    use super::*;

    fn op(kind: CigarOpKind, len: u32) -> CigarOp {
        CigarOp::new(kind, len)
    }

    fn read(pos: u32, cigar: Vec<CigarOp>, seq: &[u8]) -> AlignedRead {
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags::default(),
            cigar,
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn ref_span_and_end_are_cigar_derived_not_seq_len() {
        // 3M1D2M consumes 6 reference bases from 5 read bases.
        let r = read(100, vec![op(CigarOpKind::Match, 3), op(CigarOpKind::Deletion, 1), op(CigarOpKind::Match, 2)], b"ACGTT");
        assert_eq!(r.ref_span(), 6);
        assert_eq!(r.end(), 106);
    }

    #[test]
    fn projection_handles_softclip_insertion_deletion_refskip() {
        // 2S 3M 1I 2M 1D 2M  over seq positions; ref starts at 50.
        // read offsets: [0,1]=softclip, [2,3,4]=M, [5]=I, [6,7]=M, (1D no read), [8,9]=M
        let r = read(
            50,
            vec![
                op(CigarOpKind::SoftClip, 2),
                op(CigarOpKind::Match, 3),
                op(CigarOpKind::Insertion, 1),
                op(CigarOpKind::Match, 2),
                op(CigarOpKind::Deletion, 1),
                op(CigarOpKind::Match, 2),
            ],
            b"NNACGTTGAA", // 10 bases: 2 clipped + 3 + 1 ins + 2 + 2
        );
        let proj = r.projected_bases();
        let pairs: Vec<(u32, usize)> = proj.iter().map(|p| (p.ref_pos, p.read_offset)).collect();
        assert_eq!(
            pairs,
            vec![
                (50, 2), (51, 3), (52, 4), // first 3M
                (53, 6), (54, 7),          // 2M after the 1I (read offset skips 5)
                (56, 8), (57, 9),          // 2M after the 1D (ref skips 55)
            ]
        );
        // Soft-clipped and inserted read bases never appear in the projection.
        assert!(!pairs.iter().any(|(_, off)| *off == 0 || *off == 1 || *off == 5));
    }

    #[test]
    fn refskip_consumes_reference_only_like_a_long_intron() {
        // 2M 100N 2M (spliced long read): ref advances over the skip, no bases there.
        let r = read(0, vec![op(CigarOpKind::Match, 2), op(CigarOpKind::RefSkip, 100), op(CigarOpKind::Match, 2)], b"ACGT");
        let proj = r.projected_bases();
        assert_eq!(proj.first().unwrap().ref_pos, 0);
        assert_eq!(proj.last().unwrap().ref_pos, 103);
        assert_eq!(proj.len(), 4);
        assert_eq!(r.ref_span(), 104);
    }

    #[test]
    fn long_read_projects_every_match_base() {
        // Read-length-agnostic: a 5000-base full-match "long read".
        let seq = vec![b'A'; 5000];
        let r = read(1000, vec![op(CigarOpKind::Match, 5000)], &seq);
        let proj = r.projected_bases();
        assert_eq!(proj.len(), 5000);
        assert_eq!(proj[0], RefBase { ref_pos: 1000, read_offset: 0 });
        assert_eq!(proj[4999], RefBase { ref_pos: 5999, read_offset: 4999 });
    }

    #[test]
    fn sam_flags_decode_common_bits() {
        let f = SamFlags(SamFlags::REVERSE | SamFlags::DUPLICATE);
        assert!(f.is_reverse());
        assert!(f.is_duplicate());
        assert!(!f.is_secondary());
        assert!(!f.is_unmapped());
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test core::record`
Expected: FAIL — does not compile (`CigarOp`, `CigarOpKind`, `SamFlags`, `AlignedRead`, `RefBase` undefined).

- [ ] **Step 3: Implement the record types and projection**

Insert above the `#[cfg(test)]` block in `src/core/record.rs`:

```rust
/// CIGAR operation kind (SAM). Consuming semantics follow the SAM spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpKind {
    /// Alignment match or mismatch (`M`/`=`/`X`): consumes ref and read.
    Match,
    /// Insertion to the reference (`I`): consumes read only.
    Insertion,
    /// Deletion from the reference (`D`): consumes ref only.
    Deletion,
    /// Skipped reference region (`N`, e.g. an intron): consumes ref only.
    RefSkip,
    /// Soft clip (`S`): read bases present but unaligned; consumes read only.
    SoftClip,
    /// Hard clip (`H`): trimmed bases absent from the read; consumes neither.
    HardClip,
}

impl CigarOpKind {
    /// Whether this op consumes reference bases.
    pub fn consumes_ref(self) -> bool {
        matches!(self, CigarOpKind::Match | CigarOpKind::Deletion | CigarOpKind::RefSkip)
    }

    /// Whether this op consumes read/query bases.
    pub fn consumes_read(self) -> bool {
        matches!(self, CigarOpKind::Match | CigarOpKind::Insertion | CigarOpKind::SoftClip)
    }
}

/// A single CIGAR operation with its run length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// The operation kind.
    pub kind: CigarOpKind,
    /// Number of bases the operation spans.
    pub len: u32,
}

impl CigarOp {
    /// Construct a CIGAR operation.
    pub fn new(kind: CigarOpKind, len: u32) -> Self {
        Self { kind, len }
    }
}

/// SAM flag bitset (the subset Rosalind uses).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SamFlags(pub u16);

impl SamFlags {
    /// Template has multiple segments (paired).
    pub const PAIRED: u16 = 0x1;
    /// Each segment properly aligned (proper pair).
    pub const PROPER_PAIR: u16 = 0x2;
    /// Segment unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// Segment reverse-strand.
    pub const REVERSE: u16 = 0x10;
    /// Secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// PCR or optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// Supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;

    /// Whether the given flag bit is set.
    pub fn contains(self, bit: u16) -> bool {
        self.0 & bit != 0
    }

    /// Segment is unmapped.
    pub fn is_unmapped(self) -> bool {
        self.contains(Self::UNMAPPED)
    }

    /// Segment maps to the reverse strand (metadata only — SEQ stays forward).
    pub fn is_reverse(self) -> bool {
        self.contains(Self::REVERSE)
    }

    /// Alignment is secondary.
    pub fn is_secondary(self) -> bool {
        self.contains(Self::SECONDARY)
    }

    /// Alignment is supplementary.
    pub fn is_supplementary(self) -> bool {
        self.contains(Self::SUPPLEMENTARY)
    }

    /// Read is a PCR/optical duplicate.
    pub fn is_duplicate(self) -> bool {
        self.contains(Self::DUPLICATE)
    }
}

/// One projected base: a read base aligned to a specific reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefBase {
    /// 0-based reference position.
    pub ref_pos: u32,
    /// Offset of the base within the read's `seq`/`qual`.
    pub read_offset: usize,
}

/// A read aligned to a single contig.
///
/// `seq` is uppercase ASCII in forward-reference orientation; `pos` is the
/// 0-based leftmost reference coordinate. Read-length-agnostic.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    /// Contig id (index into the originating `ContigSet`).
    pub contig: u32,
    /// 0-based leftmost reference coordinate.
    pub pos: Position,
    /// Mapping quality (Phred-scaled).
    pub mapq: u8,
    /// SAM flags.
    pub flags: SamFlags,
    /// CIGAR describing the alignment.
    pub cigar: Vec<CigarOp>,
    /// Read sequence, uppercase ASCII, forward orientation.
    pub seq: Arc<[u8]>,
    /// Per-base Phred qualities (parallel to `seq`).
    pub qual: Arc<[u8]>,
}

impl AlignedRead {
    /// Reference bases spanned by the alignment (sum of ref-consuming ops).
    pub fn ref_span(&self) -> u32 {
        self.cigar.iter().filter(|o| o.kind.consumes_ref()).map(|o| o.len).sum()
    }

    /// Half-open reference end coordinate. CIGAR-derived — never `pos + seq_len`.
    pub fn end(&self) -> u32 {
        self.pos.0 + self.ref_span()
    }

    /// Project each `Match` base to its reference position by walking the CIGAR.
    /// Insertions/soft-clips consume read only; deletions/ref-skips consume ref
    /// only; hard-clips consume neither. This is the single CIGAR projection the
    /// pileup kernel consumes.
    pub fn projected_bases(&self) -> Vec<RefBase> {
        let mut out = Vec::new();
        let mut ref_pos = self.pos.0;
        let mut read_off: usize = 0;
        for op in &self.cigar {
            match op.kind {
                CigarOpKind::Match => {
                    for _ in 0..op.len {
                        out.push(RefBase { ref_pos, read_offset: read_off });
                        ref_pos += 1;
                        read_off += 1;
                    }
                }
                CigarOpKind::Insertion | CigarOpKind::SoftClip => {
                    read_off += op.len as usize;
                }
                CigarOpKind::Deletion | CigarOpKind::RefSkip => {
                    ref_pos += op.len;
                }
                CigarOpKind::HardClip => {}
            }
        }
        out
    }
}
```

In `src/core/mod.rs`, add:

```rust
pub mod record;
pub use record::{AlignedRead, CigarOp, CigarOpKind, RefBase, SamFlags};
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test core::record`
Expected: PASS (5 tests).

- [ ] **Step 5: Run the whole suite + lints**

Run: `cargo test` — Expected: all existing tests still pass plus the new `core::*` tests.
Run: `cargo build` — Expected: no `missing_docs`/`missing_debug_implementations` warnings from `src/core/`.

- [ ] **Step 6: Commit**

```bash
git add src/core/record.rs src/core/mod.rs
git commit -m "feat(core): canonical AlignedRead + CIGAR projection (indel/clip/refskip/long-read)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Plan self-review

- **Spec coverage** (Phase A spec §3.1 `core/`): `locus` ✔ (Task 2), `sequence` with IUPAC→N fold ✔ (Task 4), `record` with `RefSkip`, `SamFlags`, CIGAR-derived `end()`, projection, read-length-agnostic ✔ (Task 5), `budget` ✔ (Task 3), `error` ✔ (Task 1). The architecture doc's `core/sequence` 2-bit packing re-home is intentionally **deferred** — the calling vertical doesn't need 2-bit packing; the FM-index keeps using `compressed_dna` until Phase B re-homes it. Noted as out-of-scope for A1.
- **Placeholder scan:** none — every step contains complete code and exact commands.
- **Type consistency:** `Position`, `ContigSet`, `BaseCode`, `allele_index`, `CigarOpKind` (incl. `RefSkip`), `CigarOp::new`, `SamFlags`, `AlignedRead` fields (`contig`, `pos`, `mapq`, `flags`, `cigar`, `seq`, `qual`), `RefBase { ref_pos, read_offset }`, `MemoryBudget`/`WorkingSet`, `CoreError` — names used in tests match the implementations and the re-exports in `core/mod.rs`.

## Definition of done (A1)

`cargo test` green (new `core::*` tests + all pre-existing tests); `cargo build` warning-clean for `src/core/`; legacy `genomics::types` untouched and still compiling; all work on a feature branch, not `main`.
