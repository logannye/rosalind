# Phase B2 — Multi-contig FM-index + global↔local Locus mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make exact-match FM-index lookup **multi-contig** by indexing the *concatenated* genome and resolving each global hit back to a `(contig, position)` `Locus`, rejecting any hit that would straddle a contig boundary — all in-RAM, with the FM-index kernel unchanged.

**Architecture:** Two additive pieces, no change to the FM-index kernel or the aligner. (1) A coordinate **resolver** on `core::locus::ContigSet` (`resolve`, `resolve_span`) that maps a 64-bit-safe global offset to a `Locus` via binary search and enforces the bwa "bns" within-one-contig rule. (2) A new `genomics::GenomeIndex` that pairs a `BlockedFMIndex` (built over the concatenation) with its `ContigSet` and the forward reference, exposing `locate_exact(pattern) -> Vec<Locus>`. Persistence + zero-copy mmap view are **B3**; wiring the seed/chain/extend aligner + whole-genome pileup onto `GenomeIndex` is **B4**.

**Tech Stack:** Rust 2021; the existing `genomics::BlockedFMIndex` (`build`, `backward_search`, `locate_interval`) and `core::{ContigSet, Locus, Position}` (already 64-bit-`global_offset`-aware); `thiserror`.

This is stage **B2** of the Phase B design spec (`docs/superpowers/specs/2026-05-27-phase-b-genome-scale-design.md`, §7 and §9). It lands green and independently and is the correctness-critical core proven *before* persistence (B3). The `u32` suffix-array ceiling (≤ ~4.29 Gbp concatenated; spec §3.4) is enforced here with a typed error.

---

## File structure

- `src/core/locus.rs` — **Modify.** Add `ContigSet::resolve(global: u64) -> Option<Locus>` and `ContigSet::resolve_span(global: u64, len: u32) -> Option<Locus>` (+ unit tests). Pure coordinate logic; no new dependencies.
- `src/genomics/genome_index.rs` — **Create.** `GenomeIndex` (`fm`, `contigs`, `reference`), `GenomeIndexError`, `MAX_GENOME_LEN`, `validate_concat_len`, `build`, `from_named_sequences`, `locate_exact`, accessors (`contigs`/`reference`/`fm`).
- `src/genomics/mod.rs` — **Modify.** Add `mod genome_index;` (between `mod fm_index;` and `mod index;`) and `pub use genome_index::{GenomeIndex, GenomeIndexError};` (after the `fm_index::{…}` re-export).

---

## Task 1: `ContigSet` coordinate resolver (`core::locus`)

**Files:**
- Modify: `src/core/locus.rs`

- [ ] **Step 1: Write the failing tests.** In `src/core/locus.rs`, inside the existing `#[cfg(test)] mod tests { … }` (which already has `use super::*;`), add these four tests:

```rust
    #[test]
    fn resolve_maps_global_offset_to_locus() {
        let mut set = ContigSet::new();
        set.push("chr1", 10); // global 0..10
        set.push("chr2", 5); // global 10..15
        assert_eq!(set.resolve(0), Some(Locus { contig: 0, pos: Position(0) }));
        assert_eq!(set.resolve(9), Some(Locus { contig: 0, pos: Position(9) }));
        assert_eq!(set.resolve(10), Some(Locus { contig: 1, pos: Position(0) }));
        assert_eq!(set.resolve(14), Some(Locus { contig: 1, pos: Position(4) }));
    }

    #[test]
    fn resolve_out_of_range_is_none() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        assert_eq!(set.resolve(10), None); // exactly past chr1's end (the only contig)
        assert_eq!(set.resolve(100), None);
        assert_eq!(ContigSet::new().resolve(0), None); // empty set
    }

    #[test]
    fn resolve_span_within_a_contig_returns_start() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        set.push("chr2", 5);
        assert_eq!(set.resolve_span(8, 2), Some(Locus { contig: 0, pos: Position(8) })); // [8,10) ⊆ chr1
        assert_eq!(set.resolve_span(10, 5), Some(Locus { contig: 1, pos: Position(0) })); // [10,15) ⊆ chr2
    }

    #[test]
    fn resolve_span_crossing_a_boundary_is_rejected() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        set.push("chr2", 5);
        assert_eq!(set.resolve_span(8, 4), None); // [8,12) crosses chr1→chr2
        assert_eq!(set.resolve_span(12, 4), None); // [12,16) runs past genome end
    }
```

- [ ] **Step 2: Run tests to verify they fail.**

Run: `cargo test --lib core::locus 2>&1 | tail -20`
Expected: compile error — `resolve`/`resolve_span` are not methods of `ContigSet`.

- [ ] **Step 3: Implement the resolver.** In `src/core/locus.rs`, add these two methods inside `impl ContigSet { … }` (e.g. just after `total_length`):

```rust
    /// Resolve a 0-based offset in the concatenated genome to a [`Locus`].
    /// Returns `None` if `global` lies at or beyond the end of the last contig.
    /// O(log n) in the number of contigs.
    pub fn resolve(&self, global: u64) -> Option<Locus> {
        // `partition_point` gives the count of contigs whose offset is <= global;
        // the owning contig is the one just before that boundary.
        let after = self.contigs.partition_point(|c| c.global_offset <= global);
        let contig = self.contigs.get(after.checked_sub(1)?)?;
        let pos = global - contig.global_offset;
        if pos < contig.length as u64 {
            Some(Locus {
                contig: contig.id,
                pos: Position(pos as u32),
            })
        } else {
            None
        }
    }

    /// Resolve a `len`-byte span starting at `global`, returning its start
    /// [`Locus`] only if the entire span `[global, global + len)` lies within a
    /// single contig (the bwa "bns" boundary rule). Returns `None` if the span
    /// crosses a contig boundary or runs past the end of the genome.
    pub fn resolve_span(&self, global: u64, len: u32) -> Option<Locus> {
        let start = self.resolve(global)?;
        let contig = self.by_id(start.contig)?;
        let end = global.checked_add(len as u64)?;
        if end <= contig.global_offset + contig.length as u64 {
            Some(start)
        } else {
            None
        }
    }
```

- [ ] **Step 4: Run tests to verify they pass.**

Run: `cargo test --lib core::locus 2>&1 | tail -20`
Expected: all `core::locus` tests pass (the four new ones plus the pre-existing ones).

- [ ] **Step 5: Commit.**

```bash
git add src/core/locus.rs
git commit -m "feat(core/locus): ContigSet global-offset resolver (resolve + boundary-aware resolve_span)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `genomics::GenomeIndex` — multi-contig index + exact-match Locus locate

**Files:**
- Create: `src/genomics/genome_index.rs`
- Modify: `src/genomics/mod.rs`

- [ ] **Step 1: Create the module with types, signatures (unimplemented bodies), and tests.** Create `src/genomics/genome_index.rs`:

```rust
//! Multi-contig genome index: a [`BlockedFMIndex`] over the concatenated
//! reference paired with the [`ContigSet`] that maps global offsets back to
//! `(contig, position)`. Exact-match queries return [`Locus`]es and reject any
//! hit that would straddle a contig boundary (the bwa "bns" rule). The FM-index
//! kernel is unchanged — multi-contig awareness lives entirely in the coordinate
//! mapping. Persistence + a zero-copy mmap view are Phase B3; wiring the
//! seed/chain/extend aligner onto this is Phase B4.

use std::sync::Arc;

use thiserror::Error;

use crate::core::{ContigSet, Locus};
use crate::genomics::{BlockedFMIndex, FMIndexError};

/// Largest concatenated-genome length supported. The suffix array stores
/// positions as `u32` over a text of `len + 1` symbols (one sentinel), so the
/// concatenation must satisfy `len + 1 <= u32::MAX`, i.e. `len <= u32::MAX - 1`
/// (~4.29 Gbp — covers the human genome and the vast majority of edge targets).
pub const MAX_GENOME_LEN: u64 = u32::MAX as u64 - 1;

/// Errors from building or querying a [`GenomeIndex`].
#[derive(Debug, Error)]
pub enum GenomeIndexError {
    /// The concatenated genome exceeds the `u32` suffix-array ceiling.
    #[error("concatenated genome length {len} exceeds the {max}-byte limit (u32 suffix array)")]
    GenomeTooLarge {
        /// Concatenated length in bytes.
        len: u64,
        /// Maximum supported length.
        max: u64,
    },
    /// The reference was empty (no contigs / no bases).
    #[error("genome index requires a non-empty reference")]
    EmptyGenome,
    /// Failure constructing the underlying FM-index (e.g. an unsupported base).
    #[error("fm-index error: {0}")]
    FmIndex(#[from] FMIndexError),
}

/// Validate that a concatenated-genome length fits the `u32` suffix array.
fn validate_concat_len(len: u64) -> Result<(), GenomeIndexError> {
    unimplemented!()
}

/// An FM-index over a concatenated multi-contig reference plus the contig map.
#[derive(Debug)]
pub struct GenomeIndex {
    fm: BlockedFMIndex,
    contigs: ContigSet,
    reference: Arc<[u8]>,
}

impl GenomeIndex {
    /// Build an index over `reference` (the concatenation of every contig in
    /// `contigs`, in id order). `reference` should be uppercase A/C/G/T/N; other
    /// bytes surface as [`FMIndexError`]. Errors if the concatenation is empty
    /// or exceeds [`MAX_GENOME_LEN`].
    pub fn build(contigs: ContigSet, reference: Arc<[u8]>) -> Result<Self, GenomeIndexError> {
        unimplemented!()
    }

    /// The contig map for resolving `Locus`es.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }

    /// The concatenated forward reference (uppercase A/C/G/T/N).
    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    /// The underlying FM-index (consumed by the aligner in Phase B4).
    pub fn fm(&self) -> &BlockedFMIndex {
        &self.fm
    }

    /// Locate exact occurrences of `pattern`, returning up to `max_hits`
    /// [`Locus`]es with boundary-straddling hits removed, sorted by
    /// `(contig, pos)` for a deterministic, kernel-independent order.
    ///
    /// Note: `max_hits` bounds the *candidate* suffixes located before boundary
    /// filtering, so the returned count may be smaller.
    pub fn locate_exact(&self, pattern: &[u8], max_hits: usize) -> Vec<Locus> {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Position;

    // chr1 = ACGTACGT (global 0..8); chr2 = TTTTGGGG (global 8..16).
    // Concatenation: ACGTACGTTTTTGGGG
    fn two_contig_index() -> GenomeIndex {
        let mut contigs = ContigSet::new();
        contigs.push("chr1", 8);
        contigs.push("chr2", 8);
        let reference: Arc<[u8]> = Arc::from(b"ACGTACGTTTTTGGGG".to_vec().into_boxed_slice());
        GenomeIndex::build(contigs, reference).expect("build should succeed")
    }

    // Naive reference scan resolving each within-a-contig match to a Locus.
    fn naive_loci(reference: &[u8], pattern: &[u8], contigs: &ContigSet) -> Vec<Locus> {
        let mut out = Vec::new();
        if pattern.is_empty() || pattern.len() > reference.len() {
            return out;
        }
        for start in 0..=reference.len() - pattern.len() {
            if &reference[start..start + pattern.len()] == pattern {
                if let Some(locus) = contigs.resolve_span(start as u64, pattern.len() as u32) {
                    out.push(locus);
                }
            }
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn stores_contigs_and_reference() {
        let idx = two_contig_index();
        assert_eq!(idx.contigs().len(), 2);
        assert_eq!(idx.contigs().by_name("chr2").unwrap().global_offset, 8);
        assert_eq!(idx.reference(), b"ACGTACGTTTTTGGGG");
    }

    #[test]
    fn locates_a_pattern_unique_to_chr2() {
        let idx = two_contig_index();
        assert_eq!(
            idx.locate_exact(b"GGGG", 16),
            vec![Locus { contig: 1, pos: Position(4) }]
        );
    }

    #[test]
    fn locates_a_pattern_unique_to_chr1() {
        let idx = two_contig_index();
        assert_eq!(
            idx.locate_exact(b"ACGTAC", 16),
            vec![Locus { contig: 0, pos: Position(0) }]
        );
    }

    #[test]
    fn rejects_a_boundary_straddling_match() {
        let idx = two_contig_index();
        // "GTTTTT" matches only at global 6 (chr1[6..8]="GT" + chr2[0..4]="TTTT"),
        // which crosses the chr1/chr2 boundary, so no Locus is returned.
        let loci = idx.locate_exact(b"GTTTTT", 16);
        assert!(loci.is_empty(), "boundary-straddling hit must be rejected, got {loci:?}");
    }

    #[test]
    fn exact_match_matches_a_naive_scan() {
        let idx = two_contig_index();
        let reference = idx.reference().to_vec();
        for pattern in [
            b"ACGT".as_slice(),
            b"GT".as_slice(),
            b"TTTT".as_slice(),
            b"T".as_slice(),
            b"CG".as_slice(),
            b"GTTTTT".as_slice(),
        ] {
            let expected = naive_loci(&reference, pattern, idx.contigs());
            let got = idx.locate_exact(pattern, 1024);
            assert_eq!(
                got,
                expected,
                "mismatch for pattern {:?}",
                std::str::from_utf8(pattern).unwrap()
            );
        }
    }

    #[test]
    fn locate_exact_is_deterministic() {
        let idx = two_contig_index();
        assert_eq!(idx.locate_exact(b"T", 1024), idx.locate_exact(b"T", 1024));
    }

    #[test]
    fn empty_pattern_returns_nothing() {
        let idx = two_contig_index();
        assert!(idx.locate_exact(b"", 16).is_empty());
    }

    #[test]
    fn validate_concat_len_enforces_the_u32_ceiling() {
        assert!(matches!(validate_concat_len(0), Err(GenomeIndexError::EmptyGenome)));
        assert!(validate_concat_len(1).is_ok());
        assert!(validate_concat_len(MAX_GENOME_LEN).is_ok());
        assert!(matches!(
            validate_concat_len(MAX_GENOME_LEN + 1),
            Err(GenomeIndexError::GenomeTooLarge { .. })
        ));
        assert!(matches!(
            validate_concat_len(u32::MAX as u64),
            Err(GenomeIndexError::GenomeTooLarge { .. })
        ));
    }
}
```

- [ ] **Step 2: Register the module.** In `src/genomics/mod.rs`:
  - Add `mod genome_index;` immediately after the `mod fm_index;` line.
  - Add `pub use genome_index::{GenomeIndex, GenomeIndexError};` immediately after the `pub use fm_index::{ … };` block.

- [ ] **Step 3: Run tests to verify they fail.**

Run: `cargo test --lib genome_index 2>&1 | tail -20`
Expected: tests run and FAIL — `not implemented` panics from the `unimplemented!()` bodies (`build`/`locate_exact`/`validate_concat_len`).

- [ ] **Step 4: Implement the three bodies.** In `src/genomics/genome_index.rs`, replace the `unimplemented!()` bodies:

```rust
fn validate_concat_len(len: u64) -> Result<(), GenomeIndexError> {
    if len == 0 {
        return Err(GenomeIndexError::EmptyGenome);
    }
    if len > MAX_GENOME_LEN {
        return Err(GenomeIndexError::GenomeTooLarge { len, max: MAX_GENOME_LEN });
    }
    Ok(())
}
```

```rust
    pub fn build(contigs: ContigSet, reference: Arc<[u8]>) -> Result<Self, GenomeIndexError> {
        validate_concat_len(reference.len() as u64)?;
        // Mirror BWTAligner's block-size heuristic.
        let block_size = ((reference.len() as f64).sqrt().ceil() as usize).max(64);
        let fm = BlockedFMIndex::build(&reference, block_size)?;
        Ok(Self { fm, contigs, reference })
    }
```

```rust
    pub fn locate_exact(&self, pattern: &[u8], max_hits: usize) -> Vec<Locus> {
        if pattern.is_empty() {
            return Vec::new();
        }
        let pattern: Vec<u8> = pattern.iter().map(|b| b.to_ascii_uppercase()).collect();
        let interval = self.fm.backward_search(&pattern);
        let len = pattern.len() as u32;
        let mut loci: Vec<Locus> = self
            .fm
            .locate_interval(interval, max_hits)
            .into_iter()
            .filter_map(|g| self.contigs.resolve_span(g as u64, len))
            .collect();
        loci.sort_unstable();
        loci
    }
```

- [ ] **Step 5: Run tests to verify they pass.**

Run: `cargo test --lib genome_index 2>&1 | tail -20`
Expected: `test result: ok. 8 passed`. Then `cargo build --lib 2>&1 | grep -i warning` (expect none).

- [ ] **Step 6: Commit.**

```bash
git add src/genomics/genome_index.rs src/genomics/mod.rs
git commit -m "feat(genomics/genome_index): multi-contig FM-index + exact-match Locus locate (bns boundary rule)" \
  -m "GenomeIndex pairs a BlockedFMIndex over the concatenated genome with its ContigSet; locate_exact resolves global hits to (contig,pos) and rejects boundary-straddling matches. Enforces the u32 suffix-array ceiling (MAX_GENOME_LEN). FM-index kernel unchanged; aligner wiring is B4." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `GenomeIndex::from_named_sequences` (build from per-contig sequences)

**Files:**
- Modify: `src/genomics/genome_index.rs`

This is the realistic build path: the `rosalind index` subcommand (B3) and the consumers (B4) will map FASTA records into `(name, sequence)` pairs and call this. Keeping the signature on plain `(String, Vec<u8>)` keeps `genomics` independent of the `io` layer.

- [ ] **Step 1: Write the failing test.** In the `#[cfg(test)] mod tests` of `src/genomics/genome_index.rs`, add:

```rust
    #[test]
    fn from_named_sequences_builds_a_multi_contig_index() {
        let idx = GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGT".to_vec()),
            ("chr2".to_string(), b"TTTTGGGG".to_vec()),
            ("chr3".to_string(), b"CCCCAAAA".to_vec()),
        ])
        .expect("build should succeed");

        assert_eq!(idx.contigs().len(), 3);
        assert_eq!(idx.reference(), b"ACGTACGTTTTTGGGGCCCCAAAA");
        // "CCCC" is unique to chr3 at pos 0 (global 16).
        assert_eq!(
            idx.locate_exact(b"CCCC", 16),
            vec![Locus { contig: 2, pos: Position(0) }]
        );
    }

    #[test]
    fn from_named_sequences_rejects_an_empty_genome() {
        assert!(matches!(
            GenomeIndex::from_named_sequences(&[]),
            Err(GenomeIndexError::EmptyGenome)
        ));
    }
```

- [ ] **Step 2: Run tests to verify they fail.**

Run: `cargo test --lib genome_index 2>&1 | tail -20`
Expected: compile error — `from_named_sequences` is not a method of `GenomeIndex`.

- [ ] **Step 3: Implement `from_named_sequences`.** Add this method inside `impl GenomeIndex { … }` (e.g. right after `build`):

```rust
    /// Build from per-contig `(name, sequence)` pairs: assembles the
    /// [`ContigSet`] and the concatenated reference (in the given order), then
    /// delegates to [`GenomeIndex::build`]. Sequences should be uppercase
    /// A/C/G/T/N.
    pub fn from_named_sequences(named: &[(String, Vec<u8>)]) -> Result<Self, GenomeIndexError> {
        let mut contigs = ContigSet::new();
        let total: usize = named.iter().map(|(_, seq)| seq.len()).sum();
        let mut concat = Vec::with_capacity(total);
        for (name, seq) in named {
            contigs.push(name.clone(), seq.len() as u32);
            concat.extend_from_slice(seq);
        }
        Self::build(contigs, Arc::from(concat.into_boxed_slice()))
    }
```

(The empty-genome case falls through to `build` → `validate_concat_len(0)` → `EmptyGenome`.)

- [ ] **Step 4: Run tests to verify they pass.**

Run: `cargo test --lib genome_index 2>&1 | tail -20`
Expected: `test result: ok. 10 passed`. Then `cargo build --lib 2>&1 | grep -i warning` (none) and `cargo fmt --all -- --check`.

- [ ] **Step 5: Commit.**

```bash
git add src/genomics/genome_index.rs
git commit -m "feat(genomics/genome_index): from_named_sequences builder (per-contig (name,seq) → index)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before the B2 PR)

- `cargo test` — full suite green (report totals; the B1 + Phase-A suites must stay green).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo clippy --lib 2>&1 | grep -E 'genome_index|core/locus'` — no new lints in the B2 code.
- `cargo fmt --all -- --check` — clean.

## Self-Review

- **Spec coverage (B2 row of §9 + §7):** multi-contig FM-index over the concatenated genome ✔ (Task 2 `GenomeIndex::build`); global↔local `Locus` map ✔ (Task 1 `resolve`); cross-contig (boundary-straddle) rejection via the bns rule ✔ (Task 1 `resolve_span` + Task 2 `rejects_a_boundary_straddling_match`); exact-match correctness vs naive on a 2-contig reference ✔ (Task 2 `exact_match_matches_a_naive_scan`); determinism ✔ (Task 2 `locate_exact_is_deterministic` + `sort_unstable`); `u32` ceiling enforced with a typed error ✔ (Task 2 `validate_concat_len` + `MAX_GENOME_LEN`, tested at the 2³²−1 boundary). In-RAM only ✔ (no persistence/mmap here).
- **Scope:** no change to the FM-index kernel, the aligner, or `main.rs`. Persistence (B3) and aligner/pileup wiring (B4) are explicitly out.
- **Type consistency:** `resolve(global: u64) -> Option<Locus>` and `resolve_span(global: u64, len: u32) -> Option<Locus>` are used identically in `GenomeIndex::locate_exact` and the naive test helper. `GenomeIndex::build(ContigSet, Arc<[u8]>)`, `from_named_sequences(&[(String, Vec<u8>)])`, `locate_exact(&[u8], usize) -> Vec<Locus>`, and `MAX_GENOME_LEN`/`validate_concat_len` names match across tasks. `Locus`/`Position`/`ContigSet` are the `core` types; `BlockedFMIndex`/`FMIndexError` the `genomics` ones.
- **No placeholders:** every step has complete code or an exact command + expected output; the `unimplemented!()` bodies are a deliberate TDD red-state, replaced in the same task.
- **Naive cross-check rationale:** the `u32` ceiling can't be exercised with a real 4 GB allocation, so it is verified through the pure `validate_concat_len` function at the exact `MAX_GENOME_LEN` / `+1` / `u32::MAX` boundaries instead.
