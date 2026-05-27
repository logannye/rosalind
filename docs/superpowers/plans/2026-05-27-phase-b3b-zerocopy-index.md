# Phase B3b — Zero-copy persisted FM-index Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the reference FM-index a portable, deterministic on-disk artifact that a consumer **memory-maps and queries in place** — no rebuild, no allocation of the bulk index, byte-identical results to the in-RAM index — delivered at the library level (`IndexWriter`/`IndexReader` → `FmIndexView` → `backward_search`/`sa_at`/`locate_exact`).

**Architecture:** One FM query algorithm, two backings. Part 1 extracts the data surface the FM algorithm reads into a `BwtBacking` trait and rewrites `rank`/`symbol_at`/`lf_index`/`backward_search`/`sa_at`/`locate_interval`/`total` as free functions generic over it; `BlockedFMIndex` implements the trait (owned, in-RAM) and its public methods become thin wrappers, so behavior is provably unchanged (existing suites are the witnesses). Part 2 adds a deterministic, 8-aligned little-endian on-disk format that serialises a built `GenomeIndex` faithfully, and a borrowed backing `FmIndexView<'a>` that reads mmap slices via checked `slice::align_to` and implements the same `BwtBacking` — so the view *cannot* diverge on logic, only on backing, which a byte-identical equivalence gate catches loudly.

**Tech Stack:** Rust 2021 (MSRV 1.72 — **no `div_ceil`**, use `(a + b - 1) / b`); `libc::mmap` via the existing `util::mmap::MmapReadOnly`; `blake3` for the reference identity hash; `u64`/`u32` LE word arrays + `count_ones` popcount (mirrors `genomics::rank_select`/`sampled_sa`); **no new dependencies** (zero-copy via `slice::align_to`, not `bytemuck`).

This is stage **B3b** of Phase B's index-persistence pillar — the keystone under the contract thesis in `docs/OPEN_PROBLEMS.md` and the design in `docs/superpowers/specs/2026-05-27-phase-b3b-zerocopy-index-design.md`. It follows B3a (compact `SampledSuffixArray`, merged). **Out of scope (deferred by design):** the `rosalind index` subcommand + CLI `IndexReader` wiring (**B3c**); wiring `align`/`variants`/pileup onto the persisted index (**B4**); `MemoryBudget` *enforcement* + `rosalind plan` (**Phase C**); leaner/faster rank (drop the redundant per-block 2-bit BWT) (**format v2**).

---

## File structure

- `src/genomics/fm_backing.rs` — **Create (Part 1).** `pub(crate) trait BwtBacking` (the data surface the FM algorithm needs) + the generic FM operations as `pub(crate)` free functions (`rank`, `total`, `symbol_at`, `lf_index`, `backward_search`, `sa_at`, `locate_interval`). The single source of truth for the algorithm; both backings call it.
- `src/genomics/fm_index.rs` — **Modify (Part 1 + 3).** Implement `BwtBacking` for `BlockedFMIndex` (delegating to its owned fields); make the public query methods thin wrappers over `fm_backing::*`; remove the now-unused private `lf_index`. Add `pub(crate)` serialization accessors to `BWTBlock`.
- `src/genomics/rank_select.rs` — **Modify (Part 3).** Make `popcount_range` `pub(crate)`; add `pub(crate)` accessors (`bitvectors`, `superblocks`, `len`) to `RankSelectIndex`.
- `src/genomics/sampled_sa.rs` — **Modify (Part 3).** Make `RANK_STRIDE` `pub(crate)`; add `pub(crate)` accessors (`marks`, `superblocks`, `values`) to `SampledSuffixArray`.
- `src/genomics/mod.rs` — **Modify.** `mod fm_backing;`; `pub(crate) use rank_select::popcount_range;`; `pub(crate) use sampled_sa::RANK_STRIDE;`; export the new view types.
- `src/genomics/index/format.rs` — **Modify (Part 2).** Extend `SectionKind` with `Reference2bit`/`FmMeta`/`Boundaries`/`Blocks`; update `from_u32`.
- `src/genomics/index/io.rs` — **Rewrite the body (Part 2).** Replace the scaffold `write_v1`/`ContigInfo` with `IndexWriter::write_genome_index(&GenomeIndex)` (the deterministic serializer) and an `IndexReader::open` that parses the new sections into a `ReferenceIndex` holding a `ContigSet` + parsed section map + mmap.
- `src/genomics/index/view.rs` — **Create (Part 2).** `FmIndexView<'a>` (borrowed backing, `impl BwtBacking`, inherent query wrappers) and `GenomeIndexView<'a>` (`FmIndexView` + `&ContigSet` → `locate_exact -> Vec<Locus>`). The zero-copy slicing helpers (`as_u64_slice`/`as_u32_slice`) live here.
- `src/genomics/index/mod.rs` — **Modify.** `mod view;` + `pub use view::{FmIndexView, GenomeIndexView};` (+ re-export through `genomics/mod.rs`).
- `tests/index_persistence.rs` — **Create (Part 2).** The equivalence battery (view == in-RAM over many patterns, multi-contig + `N` fixtures), determinism (two builds byte-equal), no-rebuild, bounded-residency smoke, corruption rejection, self-contained round-trip.
- `docs/index-format.md` — **Create (Part 2).** The format documented as a versioned, forkable contract (sections, alignment, determinism, the faithful-serialization size note + the v2 lean path).

**Naming contract (used consistently across every task below):**
- `BwtBacking` methods: `bwt_len()`, `block_size()`, `num_blocks()`, `sentinel_pos()`, `sample_rate()`, `c_table() -> [u32; 6]`, `boundary_base(block_idx, base_index) -> u32`, `boundary_sentinel(block_idx) -> u32`, `block_rank(block_idx, FmSymbol, within) -> u32`, `block_symbol(block_idx, within) -> FmSymbol`, `sampled_at(index) -> Option<u32>`.
- Generic ops (free fns in `fm_backing`): `rank(b, FmSymbol, position)`, `total(b, FmSymbol)`, `symbol_at(b, index)`, `lf_index(b, index)`, `backward_search(b, &[u8]) -> FMInterval`, `sa_at(b, index) -> usize`, `locate_interval(b, FMInterval, max_hits) -> Vec<u32>`.
- Serializer entry point: `IndexWriter::write_genome_index(self, &GenomeIndex) -> Result<(), IndexIoError>`.
- Reader entry points: `IndexReader::open(path) -> ReferenceIndex`; `ReferenceIndex::contigs() -> &ContigSet`; `ReferenceIndex::view() -> Result<FmIndexView<'_>, IndexIoError>`; `ReferenceIndex::genome_view() -> Result<GenomeIndexView<'_>, IndexIoError>`.

---

# Part 1 — B3b.1: the `BwtBacking` seam (owned backing, behavior unchanged)

Pure refactor. The trait surface is designed (per spec §3) to serve both backings; Part 1 lands it behind the existing green suites with only the owned backing, so any logic change shows up immediately. The existing FM-index, aligner, and `genome_index` suites — plus `tests/fm_index_props.rs` — are the witnesses that behavior is unchanged.

## Task 1: `BwtBacking` trait + generic FM operations

**Files:**
- Create: `src/genomics/fm_backing.rs`
- Modify: `src/genomics/mod.rs`

- [ ] **Step 1: Create the module with the trait and the generic operations.** Create `src/genomics/fm_backing.rs`:

```rust
//! One FM query algorithm over two backings.
//!
//! The blocked FM-index query algorithm (`rank`/`symbol_at`/`lf_index`/
//! `backward_search`/`sa_at`/`locate_interval`/`total`) reads a small, fixed
//! *data surface*: block-boundary cumulative counts, per-block rank/symbol, the
//! C-table, and the sampled suffix array. [`BwtBacking`] is exactly that surface.
//!
//! Making the algorithm generic over the backing is what lets a persisted,
//! memory-mapped index ([`crate::genomics::FmIndexView`]) answer queries with
//! the *same* code as the in-RAM [`crate::genomics::BlockedFMIndex`] — the view
//! can only differ in how it fetches bytes, never in the query logic, so the
//! equivalence gate (Phase B3b.2) is a complete check.
//!
//! This is an internal seam (`pub(crate)`), not part of the public API.

use crate::genomics::{BaseCode, FMInterval, FmSymbol};

/// The data surface the blocked FM query algorithm reads.
///
/// Boundary accessors are valid for `block_idx` in `0..=num_blocks()` (the
/// terminal boundary at `num_blocks()` carries the whole-BWT cumulative counts);
/// block accessors are valid for `block_idx` in `0..num_blocks()`.
pub(crate) trait BwtBacking {
    /// Length of the BWT string (reference length + 1 sentinel).
    fn bwt_len(&self) -> usize;
    /// Block size used to partition the BWT.
    fn block_size(&self) -> usize;
    /// Number of blocks (`= ceil(bwt_len / block_size)`).
    fn num_blocks(&self) -> usize;
    /// Position of the sentinel (`$`) in the BWT string.
    fn sentinel_pos(&self) -> usize;
    /// Suffix-array sample rate (`>= 1`); used only for the `sa_at` debug guard.
    fn sample_rate(&self) -> usize;
    /// The C table (`$`,A,C,G,T,N cumulative starts) used in LF-mapping.
    fn c_table(&self) -> [u32; 6];
    /// Cumulative count of base `base_index` in `BWT[.. block_idx * block_size]`.
    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32;
    /// Cumulative sentinel count before block `block_idx`.
    fn boundary_sentinel(&self, block_idx: usize) -> u32;
    /// Rank of `symbol` within block `block_idx` over its first `within` positions.
    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32;
    /// The BWT symbol at `within` inside block `block_idx`. The caller guarantees
    /// `block_idx * block_size + within != sentinel_pos` (the sentinel is handled
    /// in [`symbol_at`]), so this never needs to report `FmSymbol::Sentinel`.
    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol;
    /// Sampled suffix-array value at BWT position `index`, if sampled.
    fn sampled_at(&self, index: usize) -> Option<u32>;
}

/// Rank of `symbol` in `BWT[..position)`.
pub(crate) fn rank<B: BwtBacking>(b: &B, symbol: FmSymbol, position: usize) -> u32 {
    let bounded = position.min(b.bwt_len());
    let block_idx = bounded / b.block_size();

    let mut count = match symbol {
        FmSymbol::Sentinel => b.boundary_sentinel(block_idx),
        FmSymbol::Base(code) => b.boundary_base(block_idx, code.index()),
    };

    // The within-block term is skipped only at the terminal position when
    // `bwt_len` is an exact multiple of `block_size` (then `block_idx ==
    // num_blocks`, which has no block — only the terminal boundary).
    if block_idx < b.num_blocks() {
        let within = bounded - block_idx * b.block_size();
        count += b.block_rank(block_idx, symbol, within);
    }

    count
}

/// Total occurrences of `symbol` across the entire BWT string.
pub(crate) fn total<B: BwtBacking>(b: &B, symbol: FmSymbol) -> u32 {
    match symbol {
        FmSymbol::Sentinel => 1,
        FmSymbol::Base(code) => b.boundary_base(b.num_blocks(), code.index()),
    }
}

/// The symbol stored at `index` in the BWT string.
pub(crate) fn symbol_at<B: BwtBacking>(b: &B, index: usize) -> FmSymbol {
    assert!(index < b.bwt_len(), "BWT index out of range");
    if index == b.sentinel_pos() {
        return FmSymbol::Sentinel;
    }
    let block_idx = index / b.block_size();
    let within = index - block_idx * b.block_size();
    b.block_symbol(block_idx, within)
}

/// LF-mapping: the BWT index whose suffix is one position earlier.
pub(crate) fn lf_index<B: BwtBacking>(b: &B, index: usize) -> usize {
    let symbol = symbol_at(b, index);
    let occ_inclusive = rank(b, symbol, index + 1);
    let c_row = b.c_table()[symbol.order()] as usize;
    c_row + occ_inclusive as usize - 1
}

/// Exact FM-index backward search over `pattern` (ASCII A/C/G/T/N).
pub(crate) fn backward_search<B: BwtBacking>(b: &B, pattern: &[u8]) -> FMInterval {
    let mut interval = FMInterval::full(b.bwt_len());

    for &ch in pattern.iter().rev() {
        let base_code = match BaseCode::from_ascii(ch) {
            Some(code) => code,
            None => return FMInterval { lower: 0, upper: 0 },
        };
        let symbol = FmSymbol::Base(base_code);
        let c_row = b.c_table()[symbol.order()];
        let new_lower = c_row + rank(b, symbol, interval.lower as usize);
        let new_upper = c_row + rank(b, symbol, interval.upper as usize);
        interval = FMInterval {
            lower: new_lower,
            upper: new_upper,
        };
        if interval.is_empty() {
            break;
        }
    }

    interval
}

/// The suffix array value corresponding to BWT index `index`.
pub(crate) fn sa_at<B: BwtBacking>(b: &B, index: usize) -> usize {
    assert!(index < b.bwt_len(), "BWT index out of range");

    let mut current = index;
    let mut lf_steps = 0usize;

    // Following LF decreases the SA value by 1 (mod n), so within `rate` steps we
    // must reach a sampled position; the recovered SA value is the sample plus
    // the number of LF steps taken.
    loop {
        if let Some(sampled) = b.sampled_at(current) {
            return sampled as usize + lf_steps;
        }
        current = lf_index(b, current);
        lf_steps += 1;
        debug_assert!(
            lf_steps <= b.sample_rate() + 1,
            "LF steps exceeded sample rate; sampling invariant violated"
        );
    }
}

/// Locate up to `max_hits` 0-based reference positions for `interval`.
pub(crate) fn locate_interval<B: BwtBacking>(
    b: &B,
    interval: FMInterval,
    max_hits: usize,
) -> Vec<u32> {
    let max_hits = max_hits.max(1);
    let mut out = Vec::new();

    let reference_len = b.bwt_len().saturating_sub(1);
    let lower = interval.lower as usize;
    let upper = interval.upper as usize;
    for bwt_idx in lower..upper {
        if out.len() >= max_hits {
            break;
        }
        let sa = sa_at(b, bwt_idx);
        if sa < reference_len {
            out.push(sa as u32);
        }
    }

    out
}
```

- [ ] **Step 2: Register the module.** In `src/genomics/mod.rs`, add `mod fm_backing;` with the other `mod` declarations (alphabetically, after `mod eval;` and before `mod fm_index;`). No `pub use` — the trait and ops are `pub(crate)` and reached via `crate::genomics::fm_backing::*`.

- [ ] **Step 3: Verify it compiles (no implementors yet).**

Run: `cargo build --lib 2>&1 | tail -20`
Expected: compiles clean. The generic functions are not yet called (dead-code is allowed because they are `pub(crate)` items in a library); no warnings about them. If a warning about unused `pub(crate)` items appears, it will resolve in Task 2 when `BlockedFMIndex` routes through them — leave it for now.

- [ ] **Step 4: Commit.**

```bash
git add src/genomics/fm_backing.rs src/genomics/mod.rs
git commit -m "feat(genomics/fm_backing): BwtBacking trait + generic FM operations" \
  -m "Extracts the data surface the blocked FM query algorithm reads into a pub(crate) BwtBacking trait, and rewrites rank/symbol_at/lf_index/backward_search/sa_at/locate_interval/total as free functions generic over it. The single source of truth for the algorithm; the owned BlockedFMIndex (next) and the persisted FmIndexView (B3b.2) both run it, so the borrowed view can only diverge on backing, not logic." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Implement `BwtBacking` for `BlockedFMIndex`; route the public methods through the generic ops

**Files:**
- Modify: `src/genomics/fm_index.rs`

The current `BlockedFMIndex` carries the algorithm as inherent methods over its owned fields. Make it implement `BwtBacking`, and turn the public query methods into thin wrappers over `fm_backing::*`. The inherent `c_table(&self) -> &[u32; 6]` is **kept as-is** (an external caller, `block_alignment::align_within_block`, does `index.c_table()[order]`); the trait's by-value `c_table()` coexists (inherent wins for the concrete type).

- [ ] **Step 1: Add a focused test that pins the trait surface + generic path (failing).** In the `#[cfg(test)] mod tests` of `src/genomics/fm_index.rs`, add:

```rust
    #[test]
    fn owned_backing_matches_generic_ops() {
        use crate::genomics::fm_backing::{self, BwtBacking};

        let reference = b"ACGTNACGTACGTACGTNNACG";
        let index = BlockedFMIndex::build(reference, 5).expect("build");

        // The trait surface and the generic ops agree with the public methods
        // (which now route through them) and with a naive BWT rank.
        let clean = sanitize_reference(reference).unwrap();
        let (bwt, _, _) = build_bwt_and_sa_samples(&clean, 1).unwrap();
        for symbol in [
            FmSymbol::Sentinel,
            FmSymbol::Base(BaseCode::A),
            FmSymbol::Base(BaseCode::C),
            FmSymbol::Base(BaseCode::G),
            FmSymbol::Base(BaseCode::T),
            FmSymbol::Base(BaseCode::N),
        ] {
            for pos in 0..=index.len() {
                let via_generic = fm_backing::rank(&index, symbol, pos);
                let via_public = index.rank(symbol, pos);
                assert_eq!(via_generic, via_public, "generic vs public rank @ {pos}");
            }
        }

        // c_table by value (trait) equals the inherent c_table by reference.
        assert_eq!(BwtBacking::c_table(&index), *index.c_table());

        // sa_at via the generic op matches the public method over the interval.
        let interval = index.backward_search(b"ACGT");
        for bwt_idx in (interval.lower as usize)..(interval.upper as usize) {
            assert_eq!(fm_backing::sa_at(&index, bwt_idx), index.sa_at(bwt_idx));
        }
    }
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --lib fm_index::tests::owned_backing_matches_generic_ops 2>&1 | tail -20`
Expected: compile error — `BlockedFMIndex` does not implement `BwtBacking` (`fm_backing::rank(&index, ...)` requires `B: BwtBacking`).

- [ ] **Step 3: Implement `BwtBacking` for `BlockedFMIndex`.** In `src/genomics/fm_index.rs`, add the import near the top (with the other `use crate::genomics::...` lines):

```rust
use crate::genomics::fm_backing::{self, BwtBacking};
```

Then add the impl block (place it immediately after the `impl BlockedFMIndex { ... }` block):

```rust
impl BwtBacking for BlockedFMIndex {
    fn bwt_len(&self) -> usize {
        self.bwt_len
    }

    fn block_size(&self) -> usize {
        self.block_size
    }

    fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    fn sentinel_pos(&self) -> usize {
        self.sentinel_pos
    }

    fn sample_rate(&self) -> usize {
        self.sampled.rate()
    }

    fn c_table(&self) -> [u32; 6] {
        self.c_table
    }

    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32 {
        self.boundaries.boundary(block_idx).cumulative_counts[base_index]
    }

    fn boundary_sentinel(&self, block_idx: usize) -> u32 {
        self.boundaries.boundary(block_idx).sentinel_count
    }

    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32 {
        self.blocks[block_idx].rank_symbol(symbol, within)
    }

    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol {
        let block = &self.blocks[block_idx];
        let base = block
            .bwt
            .base_at(within)
            .expect("BWT block should contain sequence data");
        let code = BaseCode::from_ascii(base)
            .expect("BWT symbol must be a valid DNA base except sentinel");
        FmSymbol::Base(code)
    }

    fn sampled_at(&self, index: usize) -> Option<u32> {
        self.sampled.sample_at(index)
    }
}
```

(Inside the impl, `self.c_table` is the **field** `[u32; 6]` — field access, returned by value — not the inherent method, so there is no recursion.)

- [ ] **Step 4: Route the public methods through the generic ops; remove the private `lf_index`.** In `impl BlockedFMIndex`, replace the bodies of `backward_search`, `locate_interval`, `rank`, `total`, `symbol_at`, `sa_at` with thin wrappers, and delete the private `lf_index` method (it is now `fm_backing::lf_index`, used only by the generic `sa_at`).

Replace `backward_search`:

```rust
    pub fn backward_search(&self, pattern: &[u8]) -> FMInterval {
        fm_backing::backward_search(self, pattern)
    }
```

Replace `locate_interval`:

```rust
    pub fn locate_interval(&self, interval: FMInterval, max_hits: usize) -> Vec<u32> {
        fm_backing::locate_interval(self, interval, max_hits)
    }
```

Replace `rank`:

```rust
    /// Retrieve rank of `symbol` in `BWT[..position)`.
    pub fn rank(&self, symbol: FmSymbol, position: usize) -> u32 {
        fm_backing::rank(self, symbol, position)
    }
```

Replace `total`:

```rust
    /// Total occurrences of `symbol` across the entire BWT string.
    pub fn total(&self, symbol: FmSymbol) -> u32 {
        fm_backing::total(self, symbol)
    }
```

Replace `symbol_at`:

```rust
    /// Retrieve the symbol stored at `index` in the BWT string.
    pub fn symbol_at(&self, index: usize) -> FmSymbol {
        fm_backing::symbol_at(self, index)
    }
```

Replace `sa_at`:

```rust
    /// Compute the suffix array value corresponding to the provided BWT index.
    pub fn sa_at(&self, index: usize) -> usize {
        fm_backing::sa_at(self, index)
    }
```

Delete the entire private `fn lf_index(&self, index: usize) -> usize { ... }` method.

**Keep unchanged:** the inherent `pub fn c_table(&self) -> &[u32; 6]`, `len`, `block_size`, `num_blocks`, `sentinel_position`, `sa_sample_rate`, `sampled`, `blocks`, `boundaries`, and the private `BWTBlock::rank_symbol`. The `BWTBlock` fields (`bwt`, `occ`, `sentinel_offset`, `start`, `end`) stay private — the `BwtBacking` impl is in the same module, so it reads them directly.

- [ ] **Step 5: Build and run the focused test + the FM/aligner/genome suites.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (expected: none — the generic ops are now used, so no dead-code warning), then
`cargo test --lib fm_index 2>&1 | tail -20` (the new `owned_backing_matches_generic_ops` plus `fm_index_builds_and_ranks`, `sa_at_recovers_reference_position`, `sampled_sa_is_sparse_not_dense` all pass), then
`cargo test --lib genomics 2>&1 | tail -15` (aligner + `genome_index` `locate_exact` tests — which exercise `backward_search`/`sa_at`/`locate_interval` end-to-end — stay green).

- [ ] **Step 6: Run the integration witnesses + format.**

Run: `cargo test --test fm_index_props 2>&1 | tail -15` (the proptest rank/total invariants over the public API — now the generic path — stay green), then
`cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` (full suite green; report totals), `cargo fmt --all -- --check` (clean — run `cargo fmt --all` first if needed), `cargo build 2>&1 | grep -ic warning` (0).

- [ ] **Step 7: Commit.**

```bash
git add src/genomics/fm_index.rs
git commit -m "refactor(genomics/fm_index): route BlockedFMIndex through BwtBacking" \
  -m "BlockedFMIndex implements BwtBacking over its owned fields; backward_search/locate_interval/rank/total/symbol_at/sa_at become thin wrappers over the generic fm_backing ops, and the private lf_index is removed (now fm_backing::lf_index). No observable behavior change — the FM-index/aligner/genome_index suites and tests/fm_index_props are the witnesses. The inherent c_table() -> &[u32;6] is kept for align_within_block; the trait's by-value c_table coexists." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

# Part 2 — B3b.2: format + serializer + `FmIndexView` + the equivalence gate

Extend the format, write the deterministic serializer (`GenomeIndex` → bytes), implement the borrowed `FmIndexView` (the second backing), add the `view == in-RAM` equivalence + byte-identical-build + no-rebuild + bounded-residency gates, and publish `docs/index-format.md`.

## Task 3: serialization accessors + shared rank primitives

Additive, behavior-preserving accessors so the serializer can read each structure's words, and shared `pub(crate)` rank primitives so the borrowed view computes rank with the *same* code as the owned structures.

**Files:**
- Modify: `src/genomics/rank_select.rs`, `src/genomics/sampled_sa.rs`, `src/genomics/fm_index.rs`, `src/genomics/mod.rs`

- [ ] **Step 1: Add a focused test that the new accessors expose the backing words (failing).** In the `#[cfg(test)] mod tests` of `src/genomics/fm_index.rs`, add:

```rust
    #[test]
    fn serialization_accessors_expose_backing() {
        let reference = b"ACGTNACGTACGTACGT";
        let index = BlockedFMIndex::build(reference, 6).expect("build");

        // Boundaries length is num_blocks + 1 (one per block + terminal).
        assert_eq!(index.boundaries().len(), index.num_blocks() + 1);

        // Every block exposes its 2-bit BWT, occ structure, and span.
        for block in index.blocks() {
            assert_eq!(block.end() - block.start(), block.bwt().len());
            // Each of the 5 occ rank bitvectors has ceil(n/64) words.
            let n = block.bwt().len();
            let expected_words = (n + 63) / 64;
            for bv in block.occ().bitvectors() {
                assert_eq!(bv.len(), expected_words);
            }
            // Each of the 5 occ superblock arrays has ceil(n/stride) + 1 entries.
            let expected_sb = (n + block.occ().stride() - 1) / block.occ().stride() + 1;
            for sb in block.occ().superblocks() {
                assert_eq!(sb.len(), expected_sb);
            }
            assert_eq!(block.occ().len(), n);
            let _ = block.sentinel_offset();
        }

        // The sampled SA exposes marks/superblocks/values.
        let s = index.sampled();
        assert_eq!(s.marks().len(), (s.len() + 63) / 64);
        assert_eq!(s.values().len(), s.num_samples());
        assert!(!s.superblocks().is_empty());
    }
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --lib fm_index::tests::serialization_accessors_expose_backing 2>&1 | tail -20`
Expected: compile errors — `BWTBlock::start/end/bwt/occ/sentinel_offset`, `RankSelectIndex::bitvectors/superblocks/len`, `SampledSuffixArray::marks/superblocks/values`, `CompressedBoundaries::len` do not exist.

- [ ] **Step 3: Add `RankSelectIndex` accessors and make `popcount_range` `pub(crate)`.** In `src/genomics/rank_select.rs`:

Change the free function signature from `fn popcount_range` to:

```rust
#[inline]
pub(crate) fn popcount_range(words: &[u64], start: usize, end: usize) -> u32 {
```

Add accessors to `impl RankSelectIndex` (after `totals`):

```rust
    /// The per-base rank bitvectors (`A,C,G,T,N`); bit `i` set iff `BWT[i]` is
    /// that base. Each has `ceil(len/64)` words. (Serialization.)
    pub(crate) fn bitvectors(&self) -> &[Vec<u64>; ALPHABET_SIZE] {
        &self.bitvectors
    }

    /// The per-base superblock prefix-count arrays (each `ceil(len/stride) + 1`
    /// entries). (Serialization.)
    pub(crate) fn superblocks(&self) -> &[Vec<u32>; ALPHABET_SIZE] {
        &self.superblocks
    }

    /// Number of symbols this index was built over (the block length).
    pub(crate) fn len(&self) -> usize {
        self.len
    }
```

- [ ] **Step 4: Add `SampledSuffixArray` accessors and make `RANK_STRIDE` `pub(crate)`.** In `src/genomics/sampled_sa.rs`:

Change `const RANK_STRIDE: usize = 1024;` to:

```rust
/// BWT positions per superblock in the rank structure over `marks`.
pub(crate) const RANK_STRIDE: usize = 1024;
```

Add accessors to `impl SampledSuffixArray` (after `num_samples`):

```rust
    /// The sampled-position bitvector (`ceil(bwt_len/64)` words). (Serialization.)
    pub(crate) fn marks(&self) -> &[u64] {
        &self.marks
    }

    /// The prefix-popcount superblocks over `marks`
    /// (`ceil(bwt_len/RANK_STRIDE) + 1` entries). (Serialization.)
    pub(crate) fn superblocks(&self) -> &[u32] {
        &self.superblocks
    }

    /// The sampled SA values in ascending BWT-position order
    /// (`num_samples` entries). (Serialization.)
    pub(crate) fn values(&self) -> &[u32] {
        &self.values
    }
```

- [ ] **Step 5: Add `BWTBlock` + `CompressedBoundaries` accessors.** In `src/genomics/fm_index.rs`:

Add accessors to `impl BWTBlock` (after the existing `rank_symbol` method, still inside the `impl BWTBlock` block; make them `pub(crate)`):

```rust
    /// Block start offset (inclusive) in the BWT. (Serialization.)
    pub(crate) fn start(&self) -> usize {
        self.start
    }

    /// Block end offset (exclusive) in the BWT. (Serialization.)
    pub(crate) fn end(&self) -> usize {
        self.end
    }

    /// The block's 2-bit BWT payload. (Serialization.)
    pub(crate) fn bwt(&self) -> &CompressedDNA {
        &self.bwt
    }

    /// The block's rank/select (occ) structure. (Serialization.)
    pub(crate) fn occ(&self) -> &RankSelectIndex {
        &self.occ
    }

    /// The sentinel offset within this block, if the sentinel falls here.
    /// (Serialization.)
    pub(crate) fn sentinel_offset(&self) -> Option<usize> {
        self.sentinel_offset
    }
```

Add a `len` accessor to `impl CompressedBoundaries` (after `iter`):

```rust
    /// Number of boundary entries (`num_blocks + 1`). (Serialization.)
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }
```

Confirm `RankSelectIndex` and `CompressedDNA` are imported in `fm_index.rs` — `CompressedDNA` already is; `RankSelectIndex` is used via the `occ` field type. The accessor return types `&CompressedDNA` and `&RankSelectIndex` need both in scope; they are already imported at the top (`use crate::genomics::{... CompressedDNA, RankSelectIndex ...}`). If `RankSelectIndex` is not in the `use` list, add it.

- [ ] **Step 6: Re-export the shared primitives.** In `src/genomics/mod.rs`, add after the existing `pub use` lines:

```rust
pub(crate) use rank_select::popcount_range;
pub(crate) use sampled_sa::RANK_STRIDE;
```

- [ ] **Step 7: Run the test + the rank-select/sampled-sa suites.**

Run: `cargo test --lib fm_index::tests::serialization_accessors_expose_backing 2>&1 | tail -20` (passes), then
`cargo test --lib rank_select 2>&1 | tail -10` and `cargo test --lib sampled_sa 2>&1 | tail -10` (unchanged — green), then
`cargo build --lib 2>&1 | grep -iE 'error|warning'` (none).

- [ ] **Step 8: Commit.**

```bash
git add src/genomics/rank_select.rs src/genomics/sampled_sa.rs src/genomics/fm_index.rs src/genomics/mod.rs
git commit -m "feat(genomics): pub(crate) serialization accessors + shared rank primitives" \
  -m "Adds additive pub(crate) accessors so the B3b serializer can read each structure's backing words (BWTBlock::{start,end,bwt,occ,sentinel_offset}, RankSelectIndex::{bitvectors,superblocks,len}, SampledSuffixArray::{marks,superblocks,values}, CompressedBoundaries::len), and promotes popcount_range + RANK_STRIDE to pub(crate) so the borrowed FmIndexView computes rank with the same code as the owned structures. No behavior change." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: extend the on-disk `SectionKind`

**Files:**
- Modify: `src/genomics/index/format.rs`

- [ ] **Step 1: Add the new section kinds + a discriminant round-trip test (failing).** In `src/genomics/index/format.rs`, add to the `#[cfg(test)]` module at the end of the file (create the module if absent):

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn section_kind_discriminants_round_trip() {
        for (raw, kind) in [
            (1u32, SectionKind::Contigs),
            (2, SectionKind::Reference),
            (3, SectionKind::SaSamples),
            (4, SectionKind::Reference2bit),
            (5, SectionKind::FmMeta),
            (6, SectionKind::Boundaries),
            (7, SectionKind::Blocks),
        ] {
            assert_eq!(SectionKind::from_u32(raw), Some(kind));
            assert_eq!(kind as u32, raw);
        }
        assert_eq!(SectionKind::from_u32(8), None);
        assert_eq!(SectionKind::from_u32(0), None);
    }
}
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --lib index::format 2>&1 | tail -20` (or `cargo test --lib format::tests::section_kind_discriminants_round_trip`)
Expected: compile error — `SectionKind::Reference2bit/FmMeta/Boundaries/Blocks` do not exist.

- [ ] **Step 3: Extend the enum.** In `src/genomics/index/format.rs`, replace the `SectionKind` enum and its `from_u32` with:

```rust
/// A section kind identifier used in the section table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum SectionKind {
    /// Contig metadata table (`name_len:u32, name, length:u64, global_offset:u64`).
    Contigs = 1,
    /// Reserved: the v1 scaffold's raw uppercase-ASCII reference (unused by the
    /// B3b writer, which stores `Reference2bit`).
    Reference = 2,
    /// Compact sampled suffix array (`rate, bwt_len, lengths, marks, superblocks, values`).
    SaSamples = 3,
    /// The 2-bit forward reference (`CompressedDNA`) — self-contained ref-base lookups (B4).
    Reference2bit = 4,
    /// FM-index scalar metadata (`block_size, bwt_len, sentinel_pos, sa_sample_rate, num_blocks, c_table`).
    FmMeta = 5,
    /// Block-boundary cumulative counts (`num_blocks + 1` entries of `[u32;5] + u32`).
    Boundaries = 6,
    /// Per-block BWT + occ payloads (directory of offsets, then block records).
    Blocks = 7,
}

impl SectionKind {
    pub fn from_u32(v: u32) -> Option<Self> {
        match v {
            1 => Some(SectionKind::Contigs),
            2 => Some(SectionKind::Reference),
            3 => Some(SectionKind::SaSamples),
            4 => Some(SectionKind::Reference2bit),
            5 => Some(SectionKind::FmMeta),
            6 => Some(SectionKind::Boundaries),
            7 => Some(SectionKind::Blocks),
            _ => None,
        }
    }
}
```

- [ ] **Step 4: Run the test to verify it passes.**

Run: `cargo test --lib format::tests::section_kind_discriminants_round_trip 2>&1 | tail -10`
Expected: PASS. Then `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none).

- [ ] **Step 5: Commit.**

```bash
git add src/genomics/index/format.rs
git commit -m "feat(genomics/index/format): extend SectionKind for the FM-index sections" \
  -m "Adds Reference2bit/FmMeta/Boundaries/Blocks section kinds (the B3b on-disk format extends the existing versioned v1 header; Reference stays reserved). The fixed header is unchanged." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: the deterministic serializer (`GenomeIndex` → bytes)

Replace the scaffold `write_v1`/`ContigInfo` with `IndexWriter::write_genome_index`, writing every multi-byte array as 8-aligned little-endian words. This task is verifiable on its own: the file is byte-identical across two builds of the same genome, and the section table is well-formed (every expected kind present, extents in-bounds and 8-aligned).

**Files:**
- Modify: `src/genomics/index/io.rs`

- [ ] **Step 1: Replace the IO module body down to (but not including) the reader.** In `src/genomics/index/io.rs`, replace everything from the top of the file through the end of the `impl IndexWriter { ... }` block (i.e. the imports, `IndexIoError`, `ContigInfo`, `ReferenceIndex`, `IndexWriter`, and `write_v1`) with the following. Keep the reader (`IndexReader`, `validate_header`, `read_contigs`, `read_header`, `write_header`, `write_section_table`, `read_section_table`, `read_u16/u32/u64`) for now — Task 6 rewrites the reader; this task only needs the writer to compile and run, so leave the existing reader code in place even though `read_contigs` still parses the old layout (its test is replaced in Task 6).

New top-of-file through the writer:

```rust
//! Index reader/writer for Rosalind's on-disk reference index format.
//!
//! Single-file, little-endian, section-based (see `docs/index-format.md`). Every
//! multi-byte array is written as little-endian words at an 8-byte-aligned file
//! offset so the reader can reinterpret mmap bytes as `&[u64]`/`&[u32]` with a
//! checked `slice::align_to` (zero-copy, no rebuild). Serialization is
//! deterministic: fixed section order, fixed-width LE encoding, zero padding, no
//! timestamps ⇒ byte-identical across repeated builds of the same genome.

use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use blake3::Hasher;
use thiserror::Error;

use crate::core::ContigSet;
use crate::genomics::index::format::{
    IndexHeader, IndexVersion, SectionEntry, SectionKind, ROSALIND_INDEX_MAGIC,
};
use crate::genomics::{CompressedDNA, GenomeIndex};
use crate::util::mmap::MmapReadOnly;

/// Errors for index IO.
#[derive(Debug, Error)]
pub enum IndexIoError {
    #[error("io error: {0}")]
    Io(#[from] io::Error),

    #[error("invalid index: {0}")]
    Invalid(String),
}

/// A loaded reference index: a memory-mapped byte buffer plus the parsed header,
/// section table, and contig set. The FM-index query surface is obtained as a
/// borrowed [`crate::genomics::FmIndexView`] via [`ReferenceIndex::view`].
#[derive(Debug)]
pub struct ReferenceIndex {
    /// Path to the index file.
    pub path: PathBuf,
    /// Parsed file header.
    pub header: IndexHeader,
    /// Parsed contig set (names, lengths, global offsets).
    pub(crate) contigs: ContigSet,
    /// Parsed section table.
    pub sections: Vec<SectionEntry>,
    pub(crate) mmap: MmapReadOnly,
}

impl ReferenceIndex {
    /// Return the raw bytes of the memory-mapped index file.
    pub fn bytes(&self) -> &[u8] {
        self.mmap.as_bytes()
    }

    /// The parsed contig set.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }
}

/// Writes a new index file.
#[derive(Debug)]
pub struct IndexWriter {
    file: File,
}

impl IndexWriter {
    /// Create a new index file for writing (overwriting if it exists).
    pub fn create(path: impl AsRef<Path>) -> Result<Self, IndexIoError> {
        let file = File::create(path)?;
        Ok(Self { file })
    }

    /// Serialise a built [`GenomeIndex`] into the deterministic on-disk format.
    ///
    /// Sections are written in a fixed order, each padded to an 8-byte boundary:
    /// `Contigs`, `Reference2bit`, `FmMeta`, `Boundaries`, `Blocks`, `SaSamples`,
    /// then the section table; the header is rewritten last with the section
    /// table's location. `reference_blake3` is the BLAKE3 of the uppercased ASCII
    /// reference (the stable identity used by the B4 stale-index guard).
    pub fn write_genome_index(self, index: &GenomeIndex) -> Result<(), IndexIoError> {
        let fm = index.fm();
        let contigs = index.contigs();
        if contigs.is_empty() {
            return Err(IndexIoError::Invalid(
                "index must contain at least one contig".to_string(),
            ));
        }

        let reference_blake3 = {
            let mut hasher = Hasher::new();
            hasher.update(index.reference());
            *hasher.finalize().as_bytes()
        };
        let sa_sample_rate = u32::try_from(fm.sa_sample_rate())
            .map_err(|_| IndexIoError::Invalid("sa_sample_rate exceeds u32".to_string()))?;
        let mut header =
            IndexHeader::new_v1(contigs.len() as u32, sa_sample_rate, reference_blake3);

        let mut w = BufWriter::new(self.file);

        // Reserve space for the header (rewritten at the end).
        w.write_all(&vec![0u8; IndexHeader::FIXED_SIZE])?;

        let mut sections: Vec<SectionEntry> = Vec::new();

        // --- Contigs (byte-copy read on load; alignment not required, but we pad
        // every section start to 8 uniformly). ---
        let off = pad_to_8(&mut w)?;
        for contig in contigs.iter() {
            let name = contig.name.as_bytes();
            let name_len = u32::try_from(name.len())
                .map_err(|_| IndexIoError::Invalid("contig name too long".to_string()))?;
            w.write_all(&name_len.to_le_bytes())?;
            w.write_all(name)?;
            w.write_all(&(contig.length as u64).to_le_bytes())?;
            w.write_all(&contig.global_offset.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::Contigs, off, &mut w)?;

        // --- Reference2bit: len, data_words, amb_words, then data[u64], amb[u64]. ---
        let off = pad_to_8(&mut w)?;
        let reference_2bit = CompressedDNA::compress(index.reference())
            .map_err(|e| IndexIoError::Invalid(format!("reference compress: {e}")))?;
        w.write_all(&(reference_2bit.len() as u64).to_le_bytes())?;
        w.write_all(&(reference_2bit.words().len() as u64).to_le_bytes())?;
        w.write_all(&(reference_2bit.ambiguity().bits().len() as u64).to_le_bytes())?;
        write_u64_slice(&mut w, reference_2bit.words())?;
        write_u64_slice(&mut w, reference_2bit.ambiguity().bits())?;
        push_section(&mut sections, SectionKind::Reference2bit, off, &mut w)?;

        // --- FmMeta: 5 u64 then c_table[6 u32]. ---
        let off = pad_to_8(&mut w)?;
        w.write_all(&(fm.block_size() as u64).to_le_bytes())?;
        w.write_all(&(fm.len() as u64).to_le_bytes())?;
        w.write_all(&(fm.sentinel_position() as u64).to_le_bytes())?;
        w.write_all(&(fm.sa_sample_rate() as u64).to_le_bytes())?;
        w.write_all(&(fm.num_blocks() as u64).to_le_bytes())?;
        for v in *fm.c_table() {
            w.write_all(&v.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::FmMeta, off, &mut w)?;

        // --- Boundaries: num_blocks+1 entries of [cumulative_counts; 5] + sentinel. ---
        let off = pad_to_8(&mut w)?;
        debug_assert_eq!(fm.boundaries().len(), fm.num_blocks() + 1);
        for boundary in fm.boundaries().iter() {
            for c in boundary.cumulative_counts {
                w.write_all(&c.to_le_bytes())?;
            }
            w.write_all(&boundary.sentinel_count.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::Boundaries, off, &mut w)?;

        // --- Blocks: directory of num_blocks u64 record-offsets, then records. ---
        let off = pad_to_8(&mut w)?;
        write_blocks_section(&mut w, fm)?;
        push_section(&mut sections, SectionKind::Blocks, off, &mut w)?;

        // --- SaSamples: rate, bwt_len, marks_words, sb_len, values_len, then arrays. ---
        let off = pad_to_8(&mut w)?;
        let sampled = fm.sampled();
        w.write_all(&(sampled.rate() as u64).to_le_bytes())?;
        w.write_all(&(sampled.len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.marks().len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.superblocks().len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.values().len() as u64).to_le_bytes())?;
        write_u64_slice(&mut w, sampled.marks())?;
        write_u32_slice(&mut w, sampled.superblocks())?;
        write_u32_slice(&mut w, sampled.values())?;
        push_section(&mut sections, SectionKind::SaSamples, off, &mut w)?;

        // --- Section table, then rewrite the header. ---
        let section_table_offset = pad_to_8(&mut w)?;
        write_section_table(&mut w, &sections)?;
        let section_table_end = w.stream_position()?;

        header.section_table_offset = section_table_offset;
        header.section_table_bytes = section_table_end - section_table_offset;

        w.seek(SeekFrom::Start(0))?;
        write_header(&mut w, &header)?;
        w.flush()?;
        Ok(())
    }
}

/// Pad `w` with zero bytes up to the next 8-byte boundary; return the (aligned)
/// current offset.
fn pad_to_8(w: &mut (impl Write + Seek)) -> Result<u64, IndexIoError> {
    let pos = w.stream_position()?;
    let rem = pos % 8;
    if rem != 0 {
        let pad = (8 - rem) as usize;
        w.write_all(&[0u8; 8][..pad])?;
    }
    Ok(w.stream_position()?)
}

/// Record a section spanning `[offset, current_position)`.
fn push_section(
    sections: &mut Vec<SectionEntry>,
    kind: SectionKind,
    offset: u64,
    w: &mut (impl Write + Seek),
) -> Result<(), IndexIoError> {
    let end = w.stream_position()?;
    sections.push(SectionEntry {
        kind,
        offset,
        bytes: end - offset,
    });
    Ok(())
}

fn write_u64_slice(w: &mut impl Write, words: &[u64]) -> Result<(), IndexIoError> {
    for &word in words {
        w.write_all(&word.to_le_bytes())?;
    }
    Ok(())
}

fn write_u32_slice(w: &mut impl Write, words: &[u32]) -> Result<(), IndexIoError> {
    for &word in words {
        w.write_all(&word.to_le_bytes())?;
    }
    Ok(())
}

/// Write the `Blocks` section: a `num_blocks × u64` directory of record offsets
/// (relative to the section start, each 8-aligned), then one self-describing
/// record per block. Two passes — sizes first (to lay out the directory), then
/// the records — so the writer streams without materialising the whole section.
fn write_blocks_section(w: &mut impl Write, fm: &crate::genomics::BlockedFMIndex) -> Result<(), IndexIoError> {
    let num_blocks = fm.num_blocks();
    let dir_bytes = (num_blocks as u64) * 8;

    // Pass 1: record offsets (relative to section start).
    let mut offsets = Vec::with_capacity(num_blocks);
    let mut cursor = dir_bytes;
    for block in fm.blocks() {
        offsets.push(cursor);
        cursor += block_record_len(block);
    }

    // Directory.
    for &o in &offsets {
        w.write_all(&o.to_le_bytes())?;
    }

    // Pass 2: records.
    for block in fm.blocks() {
        write_block_record(w, block)?;
    }
    Ok(())
}

/// Byte length of one block record (header + payload), padded to 8.
fn block_record_len(block: &crate::genomics::BWTBlock) -> u64 {
    let bwt_data_words = block.bwt().words().len() as u64;
    let bwt_amb_words = block.bwt().ambiguity().bits().len() as u64;
    let occ_bitvec_words = block.occ().bitvectors()[0].len() as u64;
    let occ_superblock_len = block.occ().superblocks()[0].len() as u64;
    let body = 64 // 8 header u64/i64 fields
        + 8 * bwt_data_words
        + 8 * bwt_amb_words
        + 8 * 5 * occ_bitvec_words
        + 4 * 5 * occ_superblock_len
        + 4 * 5; // totals[5]
    (body + 7) / 8 * 8
}

/// Write one block record: 8 header fields (`start,end,sentinel_offset,stride,
/// bwt_data_words,bwt_amb_words,occ_bitvec_words,occ_superblock_len`), then
/// `bwt_data[u64]`, `bwt_amb[u64]`, 5×`occ_bitvec[u64]`, 5×`occ_superblock[u32]`,
/// `totals[u32;5]`, padded to 8.
fn write_block_record(w: &mut impl Write, block: &crate::genomics::BWTBlock) -> Result<(), IndexIoError> {
    let bwt = block.bwt();
    let occ = block.occ();
    let bwt_data_words = bwt.words().len() as u64;
    let bwt_amb_words = bwt.ambiguity().bits().len() as u64;
    let occ_bitvec_words = occ.bitvectors()[0].len() as u64;
    let occ_superblock_len = occ.superblocks()[0].len() as u64;

    let sentinel_offset: i64 = match block.sentinel_offset() {
        Some(o) => o as i64,
        None => -1,
    };

    w.write_all(&(block.start() as u64).to_le_bytes())?;
    w.write_all(&(block.end() as u64).to_le_bytes())?;
    w.write_all(&sentinel_offset.to_le_bytes())?;
    w.write_all(&(occ.stride() as u64).to_le_bytes())?;
    w.write_all(&bwt_data_words.to_le_bytes())?;
    w.write_all(&bwt_amb_words.to_le_bytes())?;
    w.write_all(&occ_bitvec_words.to_le_bytes())?;
    w.write_all(&occ_superblock_len.to_le_bytes())?;

    write_u64_slice(w, bwt.words())?;
    write_u64_slice(w, bwt.ambiguity().bits())?;
    for bv in occ.bitvectors() {
        write_u64_slice(w, bv)?;
    }
    for sb in occ.superblocks() {
        write_u32_slice(w, sb)?;
    }
    write_u32_slice(w, &occ.totals())?;

    // Pad the record to an 8-byte boundary (the u32 arrays may leave a 4-byte tail).
    let body = 64
        + 8 * bwt_data_words
        + 8 * bwt_amb_words
        + 8 * 5 * occ_bitvec_words
        + 4 * 5 * occ_superblock_len
        + 4 * 5;
    let pad = ((body + 7) / 8 * 8 - body) as usize;
    if pad != 0 {
        w.write_all(&[0u8; 8][..pad])?;
    }
    Ok(())
}
```

Note on imports: `io.rs` needs only `{CompressedDNA, GenomeIndex}` from `crate::genomics` (plus the `core`/`format`/`mmap` imports shown). The writer reaches `BlockedFMIndex`/`BWTBlock` by **full path** in the helper signatures (`fm: &crate::genomics::BlockedFMIndex`, `block: &crate::genomics::BWTBlock`), so they need no `use`. The FM-query imports (`BwtBacking`, `FmSymbol`, `RANK_STRIDE`, `popcount_range`, `BaseCode`, `FMInterval`) are **not** used in `io.rs` at all — they live in `view.rs` (Task 7). Keep this import set exactly as written to stay warning-free.

- [ ] **Step 2: Make the writer compile (the old reader still references `ContigInfo`).** The old reader's `read_contigs` returns `Vec<ContigInfo>` and `IndexReader::open` builds a `ReferenceIndex` with a `contigs: Vec<ContigInfo>` field — but the new `ReferenceIndex` has `contigs: ContigSet`. To keep the crate compiling between Task 5 and Task 6, temporarily stub the reader: replace the body of `IndexReader::open` with `unimplemented!("reader rewritten in Task 6")` and delete `read_contigs` + the `ContigInfo` struct + the old `#[cfg(test)] mod tests` (the `roundtrip_minimal_v1_index` test, which exercised the removed `write_v1`). Keep `validate_header`, `read_header`, `write_header`, `write_section_table`, `read_section_table`, `read_u16/u32/u64` (Task 6 uses them).

Concretely:
- Delete `pub struct ContigInfo { ... }`.
- Replace `impl IndexReader { pub fn open(...) -> ... { ... } }` body with:

```rust
impl IndexReader {
    /// Open and memory-map an existing index file. (Rewritten in B3b.2 Task 6.)
    pub fn open(_path: impl AsRef<Path>) -> Result<ReferenceIndex, IndexIoError> {
        unimplemented!("IndexReader::open is implemented in Task 6")
    }
}
```

- Delete `fn read_contigs(...)`.
- Delete the entire `#[cfg(test)] mod tests { ... }` block at the bottom (replaced by Task 5's determinism test below and Task 6's reader tests).

- [ ] **Step 3: Add the serializer's determinism + structural unit tests.** Append a fresh `#[cfg(test)] mod tests` to `src/genomics/index/io.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::index::format::SectionKind;
    use crate::genomics::GenomeIndex;
    use std::env;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        env::temp_dir().join(format!("rosalind-b3b-{suffix}-{nanos}.idx"))
    }

    fn sample_index() -> GenomeIndex {
        GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGTNNACGTACG".to_vec()),
            ("chr2".to_string(), b"TTTTGGGGCCCCAAAANNNN".to_vec()),
        ])
        .expect("build")
    }

    #[test]
    fn serialized_index_is_byte_identical_across_builds() {
        let idx = sample_index();
        let p1 = temp_path("det1");
        let p2 = temp_path("det2");
        IndexWriter::create(&p1).unwrap().write_genome_index(&idx).unwrap();
        IndexWriter::create(&p2).unwrap().write_genome_index(&idx).unwrap();
        let a = std::fs::read(&p1).unwrap();
        let b = std::fs::read(&p2).unwrap();
        assert_eq!(a, b, "two serializations of the same genome must be byte-identical");
        let _ = std::fs::remove_file(p1);
        let _ = std::fs::remove_file(p2);
    }

    #[test]
    fn serialized_sections_are_present_aligned_and_in_bounds() {
        let idx = sample_index();
        let path = temp_path("sections");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
        let bytes = std::fs::read(&path).unwrap();

        let header = read_header(&bytes).unwrap();
        validate_header(&header).unwrap();
        let sections = read_section_table(&bytes, &header).unwrap();

        for kind in [
            SectionKind::Contigs,
            SectionKind::Reference2bit,
            SectionKind::FmMeta,
            SectionKind::Boundaries,
            SectionKind::Blocks,
            SectionKind::SaSamples,
        ] {
            let s = sections
                .iter()
                .find(|s| s.kind == kind)
                .unwrap_or_else(|| panic!("missing section {kind:?}"));
            assert_eq!(s.offset % 8, 0, "{kind:?} not 8-aligned");
            let end = s.offset + s.bytes;
            assert!(end as usize <= bytes.len(), "{kind:?} out of bounds");
        }
        let _ = std::fs::remove_file(path);
    }
}
```

- [ ] **Step 4: Build and run the serializer tests.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none — `IndexReader::open` is an `unimplemented!()` stub, which compiles), then
`cargo test --lib index::io 2>&1 | tail -20`
Expected: `serialized_index_is_byte_identical_across_builds` and `serialized_sections_are_present_aligned_and_in_bounds` PASS. (Other crates' tests that call `IndexReader::open` — there are none outside this module — would panic; confirm `cargo test --lib 2>&1 | grep -E "test result:"` is otherwise green except any test that *intentionally* hits the stub, of which there are none.)

- [ ] **Step 5: Commit.**

```bash
git add src/genomics/index/io.rs
git commit -m "feat(genomics/index/io): deterministic GenomeIndex serializer" \
  -m "IndexWriter::write_genome_index serialises a built GenomeIndex into the 8-aligned little-endian section format (Contigs+global_offset, Reference2bit, FmMeta, Boundaries, self-describing Blocks directory+records, SaSamples). Byte-identical across repeated builds. The scaffold write_v1/ContigInfo and the old reader are removed; IndexReader::open is stubbed pending Task 6." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: the reader — `IndexReader::open` parses the new sections into a `ReferenceIndex`

Parse the header, section table, and the new `Contigs` layout (with `global_offset`) into a `ContigSet`; validate every section's extent and 8-alignment. The FM sections are *located* here (offsets validated) and *interpreted* by the view in Task 7.

**Files:**
- Modify: `src/genomics/index/io.rs`

- [ ] **Step 1: Add a reader round-trip test (failing).** In the `#[cfg(test)] mod tests` of `src/genomics/index/io.rs`, add:

```rust
    #[test]
    fn open_parses_header_contigs_and_sections() {
        let idx = sample_index();
        let path = temp_path("open");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();

        let loaded = IndexReader::open(&path).expect("open");
        assert_eq!(loaded.header.magic, ROSALIND_INDEX_MAGIC);

        // Contigs round-trip with names, lengths, and global offsets.
        let contigs = loaded.contigs();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs.by_name("chr1").unwrap().length, 17);
        assert_eq!(contigs.by_name("chr2").unwrap().global_offset, 17);

        // reference_blake3 matches the uppercased ASCII reference.
        let mut hasher = blake3::Hasher::new();
        hasher.update(idx.reference());
        assert_eq!(*hasher.finalize().as_bytes(), loaded.header.reference_blake3);

        // The FM sections are all present.
        for kind in [
            SectionKind::FmMeta,
            SectionKind::Boundaries,
            SectionKind::Blocks,
            SectionKind::SaSamples,
        ] {
            assert!(loaded.sections.iter().any(|s| s.kind == kind), "missing {kind:?}");
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn open_rejects_a_truncated_index() {
        let idx = sample_index();
        let path = temp_path("trunc");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.truncate(bytes.len() / 2);
        let truncated = temp_path("trunc-half");
        std::fs::write(&truncated, &bytes).unwrap();
        assert!(IndexReader::open(&truncated).is_err());
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(truncated);
    }

    #[test]
    fn open_rejects_bad_magic() {
        let idx = sample_index();
        let path = temp_path("magic");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        bytes[0] = b'X';
        let corrupt = temp_path("magic-bad");
        std::fs::write(&corrupt, &bytes).unwrap();
        assert!(IndexReader::open(&corrupt).is_err());
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(corrupt);
    }
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --lib index::io::tests::open_parses_header_contigs_and_sections 2>&1 | tail -20`
Expected: panic — `IndexReader::open` is the `unimplemented!()` stub.

- [ ] **Step 3: Implement `IndexReader::open` + `read_contigs`.** In `src/genomics/index/io.rs`, replace the stubbed `impl IndexReader { ... }` with the real reader, and add `read_contigs` (parsing the new layout into a `ContigSet`), plus a section-extent validator. `validate_header`, `read_header`, `read_section_table`, `read_u16/u32/u64` are reused as-is.

```rust
/// Reads an existing index file.
#[derive(Debug)]
pub struct IndexReader;

impl IndexReader {
    /// Open and memory-map an existing index file, parsing the header, section
    /// table, and contig set, and validating every section's extent + 8-alignment.
    /// The FM query surface is obtained via [`ReferenceIndex::view`].
    pub fn open(path: impl AsRef<Path>) -> Result<ReferenceIndex, IndexIoError> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let mmap = MmapReadOnly::map(&file)?;
        let bytes = mmap.as_bytes();

        let header = read_header(bytes)?;
        validate_header(&header)?;

        let version = IndexVersion::from_u16(header.version).ok_or_else(|| {
            IndexIoError::Invalid(format!("unsupported index version {}", header.version))
        })?;
        if version != IndexVersion::V1 {
            return Err(IndexIoError::Invalid(format!(
                "unsupported index version {version:?}"
            )));
        }

        let sections = read_section_table(bytes, &header)?;
        validate_sections(bytes, &sections)?;
        let contigs = read_contigs(bytes, &header, &sections)?;

        Ok(ReferenceIndex {
            path,
            header,
            contigs,
            sections,
            mmap,
        })
    }
}

/// Validate that every section lies within the file and starts 8-aligned.
fn validate_sections(bytes: &[u8], sections: &[SectionEntry]) -> Result<(), IndexIoError> {
    for s in sections {
        if s.offset % 8 != 0 {
            return Err(IndexIoError::Invalid(format!(
                "section {:?} offset {} is not 8-aligned",
                s.kind, s.offset
            )));
        }
        let end = s
            .offset
            .checked_add(s.bytes)
            .ok_or_else(|| IndexIoError::Invalid("section extent overflow".to_string()))?;
        if end as usize > bytes.len() {
            return Err(IndexIoError::Invalid(format!(
                "section {:?} out of bounds",
                s.kind
            )));
        }
    }
    Ok(())
}

/// Parse the `Contigs` section into a `ContigSet`, recomputing global offsets via
/// `push` and asserting they match the stored offsets (an integrity check).
fn read_contigs(
    bytes: &[u8],
    header: &IndexHeader,
    sections: &[SectionEntry],
) -> Result<ContigSet, IndexIoError> {
    let section = sections
        .iter()
        .find(|s| s.kind == SectionKind::Contigs)
        .ok_or_else(|| IndexIoError::Invalid("missing contigs section".to_string()))?;
    let mut offset = section.offset as usize;
    let end = section.offset as usize + section.bytes as usize;

    let mut contigs = ContigSet::new();
    for _ in 0..header.contig_count {
        let name_len = read_u32(bytes, &mut offset)? as usize;
        if offset + name_len > end {
            return Err(IndexIoError::Invalid("contig name out of bounds".to_string()));
        }
        let name = std::str::from_utf8(&bytes[offset..offset + name_len])
            .map_err(|_| IndexIoError::Invalid("contig name not valid utf-8".to_string()))?
            .to_string();
        offset += name_len;
        let length = read_u64(bytes, &mut offset)?;
        let stored_global = read_u64(bytes, &mut offset)?;

        let length_u32 = u32::try_from(length)
            .map_err(|_| IndexIoError::Invalid("contig length exceeds u32".to_string()))?;
        let id = contigs.push(name, length_u32);
        let recomputed = contigs.by_id(id).expect("just pushed").global_offset;
        if recomputed != stored_global {
            return Err(IndexIoError::Invalid(format!(
                "contig {id} global_offset {stored_global} != recomputed {recomputed}"
            )));
        }
    }
    Ok(contigs)
}
```

- [ ] **Step 4: Run the reader tests + full lib suite.**

Run: `cargo test --lib index::io 2>&1 | tail -25`
Expected: `open_parses_header_contigs_and_sections`, `open_rejects_a_truncated_index`, `open_rejects_bad_magic`, plus Task 5's two tests — all PASS. Then `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none).

- [ ] **Step 5: Commit.**

```bash
git add src/genomics/index/io.rs
git commit -m "feat(genomics/index/io): IndexReader::open parses the FM-index sections" \
  -m "open() mmaps the file, validates the header + every section extent/8-alignment, parses the new Contigs layout into a ContigSet (recomputing and integrity-checking global offsets), and locates the FM sections. Rejects truncated/bad-magic files. The borrowed FmIndexView (Task 7) interprets the FM sections." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: `FmIndexView<'a>` (the borrowed backing) + `GenomeIndexView` + the equivalence gate

The second backing: a zero-copy view over mmap slices that implements `BwtBacking` and therefore answers `backward_search`/`sa_at`/`locate_interval` with the *same* generic ops as the in-RAM index. The equivalence gate proves byte-identical results.

**Files:**
- Create: `src/genomics/index/view.rs`
- Modify: `src/genomics/index/mod.rs`, `src/genomics/mod.rs`, `src/genomics/index/io.rs`

- [ ] **Step 1: Create the view module.** Create `src/genomics/index/view.rs`:

```rust
//! Borrowed, zero-copy query view over a memory-mapped index.
//!
//! [`FmIndexView`] reads the FM-index sections directly from the mmap as
//! `&[u64]`/`&[u32]` (via a checked `slice::align_to`) and implements
//! [`crate::genomics::fm_backing::BwtBacking`]. Because the FM query algorithm is
//! generic over that trait, the view runs the *same* code as the in-RAM
//! [`crate::genomics::BlockedFMIndex`] — it can differ only in how it fetches
//! bytes, never in logic. The zero-copy reinterpretation is value-correct only on
//! little-endian hosts (the file is little-endian; [`FmIndexView::new`] rejects
//! big-endian hosts).

use crate::core::{ContigSet, Locus};
use crate::genomics::fm_backing::{self, BwtBacking};
use crate::genomics::index::format::{SectionEntry, SectionKind};
use crate::genomics::index::io::IndexIoError;
use crate::genomics::{popcount_range, BaseCode, FMInterval, FmSymbol, RANK_STRIDE};

/// A borrowed FM-index over memory-mapped bytes.
#[derive(Debug)]
pub struct FmIndexView<'a> {
    bwt_len: usize,
    block_size: usize,
    num_blocks: usize,
    sentinel_pos: usize,
    sample_rate: usize,
    c_table: [u32; 6],
    /// `6 * (num_blocks + 1)` u32: `[c0..c4, sentinel]` per boundary entry.
    boundaries: &'a [u32],
    /// `num_blocks` u64 record offsets, relative to the start of `blocks`.
    block_dir: &'a [u64],
    /// The whole `Blocks` section payload (records are sliced on demand).
    blocks: &'a [u8],
    sampled: SampledView<'a>,
}

#[derive(Debug)]
struct SampledView<'a> {
    rate: usize,
    bwt_len: usize,
    marks: &'a [u64],
    superblocks: &'a [u32],
    values: &'a [u32],
}

/// A parsed block record (computed on demand from the directory; holds only
/// borrowed slices into the mmap).
struct BlockView<'a> {
    start: usize,
    end: usize,
    sentinel_offset: Option<usize>,
    stride: usize,
    bwt_data: &'a [u64],
    bwt_amb: &'a [u64],
    occ_bitvecs: [&'a [u64]; 5],
    occ_superblocks: [&'a [u32]; 5],
}

impl<'a> FmIndexView<'a> {
    /// Build a view from the mmap `bytes` and the parsed `sections`. Validates
    /// section presence, lengths, and the little-endian host requirement.
    pub(crate) fn new(bytes: &'a [u8], sections: &[SectionEntry]) -> Result<Self, IndexIoError> {
        if cfg!(target_endian = "big") {
            return Err(IndexIoError::Invalid(
                "zero-copy index view requires a little-endian host".to_string(),
            ));
        }

        let fm_meta = section_bytes(bytes, sections, SectionKind::FmMeta)?;
        let mut o = 0usize;
        let block_size = read_u64(fm_meta, &mut o)? as usize;
        let bwt_len = read_u64(fm_meta, &mut o)? as usize;
        let sentinel_pos = read_u64(fm_meta, &mut o)? as usize;
        let sample_rate = read_u64(fm_meta, &mut o)? as usize;
        let num_blocks = read_u64(fm_meta, &mut o)? as usize;
        let mut c_table = [0u32; 6];
        for slot in &mut c_table {
            *slot = read_u32(fm_meta, &mut o)?;
        }

        let boundaries_bytes = section_bytes(bytes, sections, SectionKind::Boundaries)?;
        let boundaries = as_u32_slice(boundaries_bytes)?;
        if boundaries.len() != 6 * (num_blocks + 1) {
            return Err(IndexIoError::Invalid("boundaries length mismatch".to_string()));
        }

        let blocks = section_bytes(bytes, sections, SectionKind::Blocks)?;
        let dir_bytes = num_blocks
            .checked_mul(8)
            .ok_or_else(|| IndexIoError::Invalid("block directory overflow".to_string()))?;
        if dir_bytes > blocks.len() {
            return Err(IndexIoError::Invalid("block directory out of bounds".to_string()));
        }
        let block_dir = as_u64_slice(&blocks[..dir_bytes])?;

        // SaSamples.
        let sa = section_bytes(bytes, sections, SectionKind::SaSamples)?;
        let mut so = 0usize;
        let rate = read_u64(sa, &mut so)? as usize;
        let s_bwt_len = read_u64(sa, &mut so)? as usize;
        let marks_words = read_u64(sa, &mut so)? as usize;
        let sb_len = read_u64(sa, &mut so)? as usize;
        let values_len = read_u64(sa, &mut so)? as usize;
        let marks = as_u64_slice(slice_exact(sa, &mut so, marks_words * 8)?)?;
        let superblocks = as_u32_slice(slice_exact(sa, &mut so, sb_len * 4)?)?;
        let values = as_u32_slice(slice_exact(sa, &mut so, values_len * 4)?)?;

        Ok(Self {
            bwt_len,
            block_size,
            num_blocks,
            sentinel_pos,
            sample_rate,
            c_table,
            boundaries,
            block_dir,
            blocks,
            sampled: SampledView {
                rate,
                bwt_len: s_bwt_len,
                marks,
                superblocks,
                values,
            },
        })
    }

    /// Exact FM-index backward search (delegates to the shared generic op).
    pub fn backward_search(&self, pattern: &[u8]) -> FMInterval {
        fm_backing::backward_search(self, pattern)
    }

    /// Suffix array value at BWT index `index` (delegates to the shared generic op).
    pub fn sa_at(&self, index: usize) -> usize {
        fm_backing::sa_at(self, index)
    }

    /// Locate up to `max_hits` 0-based reference positions for `interval`.
    pub fn locate_interval(&self, interval: FMInterval, max_hits: usize) -> Vec<u32> {
        fm_backing::locate_interval(self, interval, max_hits)
    }

    /// Parse block record `block_idx` from the directory (borrowed slices only).
    fn block(&self, block_idx: usize) -> BlockView<'a> {
        let rec = self.block_dir[block_idx] as usize;
        let s = &self.blocks[rec..];
        let mut o = 0usize;
        let start = read_u64_infallible(s, &mut o) as usize;
        let end = read_u64_infallible(s, &mut o) as usize;
        let sentinel_raw = read_i64_infallible(s, &mut o);
        let stride = read_u64_infallible(s, &mut o) as usize;
        let bwt_data_words = read_u64_infallible(s, &mut o) as usize;
        let bwt_amb_words = read_u64_infallible(s, &mut o) as usize;
        let occ_bitvec_words = read_u64_infallible(s, &mut o) as usize;
        let occ_superblock_len = read_u64_infallible(s, &mut o) as usize;

        let bwt_data = as_u64_slice(&s[o..o + bwt_data_words * 8]).expect("aligned");
        o += bwt_data_words * 8;
        let bwt_amb = as_u64_slice(&s[o..o + bwt_amb_words * 8]).expect("aligned");
        o += bwt_amb_words * 8;

        let mut occ_bitvecs: [&[u64]; 5] = [&[]; 5];
        for slot in &mut occ_bitvecs {
            *slot = as_u64_slice(&s[o..o + occ_bitvec_words * 8]).expect("aligned");
            o += occ_bitvec_words * 8;
        }
        let mut occ_superblocks: [&[u32]; 5] = [&[]; 5];
        for slot in &mut occ_superblocks {
            *slot = as_u32_slice(&s[o..o + occ_superblock_len * 4]).expect("aligned");
            o += occ_superblock_len * 4;
        }

        BlockView {
            start,
            end,
            sentinel_offset: if sentinel_raw < 0 {
                None
            } else {
                Some(sentinel_raw as usize)
            },
            stride,
            bwt_data,
            bwt_amb,
            occ_bitvecs,
            occ_superblocks,
        }
    }
}

impl BwtBacking for FmIndexView<'_> {
    fn bwt_len(&self) -> usize {
        self.bwt_len
    }
    fn block_size(&self) -> usize {
        self.block_size
    }
    fn num_blocks(&self) -> usize {
        self.num_blocks
    }
    fn sentinel_pos(&self) -> usize {
        self.sentinel_pos
    }
    fn sample_rate(&self) -> usize {
        self.sample_rate
    }
    fn c_table(&self) -> [u32; 6] {
        self.c_table
    }
    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32 {
        self.boundaries[block_idx * 6 + base_index]
    }
    fn boundary_sentinel(&self, block_idx: usize) -> u32 {
        self.boundaries[block_idx * 6 + 5]
    }
    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32 {
        let block = self.block(block_idx);
        let n = block.end - block.start;
        let bounded = within.min(n);
        match symbol {
            FmSymbol::Sentinel => match block.sentinel_offset {
                Some(off) if off < bounded => 1,
                _ => 0,
            },
            FmSymbol::Base(code) => {
                // Mirror RankSelectIndex::rank over the borrowed occ slices.
                let bi = code.index();
                let sb = bounded / block.stride;
                let within_start = sb * block.stride;
                let mut count = block.occ_superblocks[bi][sb]
                    + popcount_range(block.occ_bitvecs[bi], within_start, bounded);
                // Mirror BWTBlock::rank_symbol's N correction: the sentinel is
                // stored as `N` in the block BWT but is not a real `N`.
                if code == BaseCode::N {
                    if let Some(off) = block.sentinel_offset {
                        if off < bounded {
                            count = count.saturating_sub(1);
                        }
                    }
                }
                count
            }
        }
    }
    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol {
        let block = self.block(block_idx);
        // Mirror CompressedDNA::base_at: ambiguity (1 bit/base) marks `N`, else
        // a 2-bit code (32 bases per u64 word).
        let amb_word = block.bwt_amb[within / 64];
        if amb_word & (1u64 << (within % 64)) != 0 {
            return FmSymbol::Base(BaseCode::N);
        }
        let data_word = block.bwt_data[within / 32];
        let code = ((data_word >> ((within % 32) * 2)) & 0b11) as u8;
        let base = match code {
            0 => BaseCode::A,
            1 => BaseCode::C,
            2 => BaseCode::G,
            _ => BaseCode::T,
        };
        FmSymbol::Base(base)
    }
    fn sampled_at(&self, index: usize) -> Option<u32> {
        debug_assert!(index < self.sampled.bwt_len);
        let word = self.sampled.marks[index / 64];
        if word & (1u64 << (index % 64)) == 0 {
            return None;
        }
        let sb = index / RANK_STRIDE;
        let rank = self.sampled.superblocks[sb]
            + popcount_range(self.sampled.marks, sb * RANK_STRIDE, index);
        Some(self.sampled.values[rank as usize])
    }
}

/// An FM-index view paired with the contig set: exact-match queries to `Locus`es.
#[derive(Debug)]
pub struct GenomeIndexView<'a> {
    fm: FmIndexView<'a>,
    contigs: &'a ContigSet,
}

impl<'a> GenomeIndexView<'a> {
    pub(crate) fn new(fm: FmIndexView<'a>, contigs: &'a ContigSet) -> Self {
        Self { fm, contigs }
    }

    /// The borrowed FM-index.
    pub fn fm(&self) -> &FmIndexView<'a> {
        &self.fm
    }

    /// The contig set.
    pub fn contigs(&self) -> &ContigSet {
        self.contigs
    }

    /// Locate exact occurrences of `pattern`, returning up to `max_hits` `Locus`es
    /// with boundary-straddling hits removed, sorted by `(contig, pos)` — the same
    /// surface (and result) as `GenomeIndex::locate_exact`.
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
}

/// The bytes of section `kind`.
fn section_bytes<'a>(
    bytes: &'a [u8],
    sections: &[SectionEntry],
    kind: SectionKind,
) -> Result<&'a [u8], IndexIoError> {
    let s = sections
        .iter()
        .find(|s| s.kind == kind)
        .ok_or_else(|| IndexIoError::Invalid(format!("missing section {kind:?}")))?;
    let start = s.offset as usize;
    let end = start + s.bytes as usize;
    bytes
        .get(start..end)
        .ok_or_else(|| IndexIoError::Invalid(format!("section {kind:?} out of bounds")))
}

/// Reinterpret a length-8k byte slice (8-aligned start) as `&[u64]`.
fn as_u64_slice(bytes: &[u8]) -> Result<&[u64], IndexIoError> {
    // SAFETY: any byte pattern is a valid `u64`; we require an empty prefix
    // (start is 8-aligned) and empty suffix (length is a multiple of 8). Values
    // are correct on little-endian hosts (checked in `FmIndexView::new`).
    let (prefix, mid, suffix) = unsafe { bytes.align_to::<u64>() };
    if !prefix.is_empty() || !suffix.is_empty() {
        return Err(IndexIoError::Invalid("misaligned u64 array".to_string()));
    }
    Ok(mid)
}

/// Reinterpret a length-4k byte slice (4-aligned start) as `&[u32]`.
fn as_u32_slice(bytes: &[u8]) -> Result<&[u32], IndexIoError> {
    // SAFETY: as `as_u64_slice`, for `u32` (4-aligned, length a multiple of 4).
    let (prefix, mid, suffix) = unsafe { bytes.align_to::<u32>() };
    if !prefix.is_empty() || !suffix.is_empty() {
        return Err(IndexIoError::Invalid("misaligned u32 array".to_string()));
    }
    Ok(mid)
}

/// Take the next `len` bytes from `bytes` starting at `*offset`, advancing it.
fn slice_exact<'a>(
    bytes: &'a [u8],
    offset: &mut usize,
    len: usize,
) -> Result<&'a [u8], IndexIoError> {
    let end = offset
        .checked_add(len)
        .ok_or_else(|| IndexIoError::Invalid("slice overflow".to_string()))?;
    let out = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("slice out of bounds".to_string()))?;
    *offset = end;
    Ok(out)
}

fn read_u32(bytes: &[u8], offset: &mut usize) -> Result<u32, IndexIoError> {
    let end = *offset + 4;
    let buf = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("unexpected EOF".to_string()))?;
    *offset = end;
    Ok(u32::from_le_bytes(buf.try_into().unwrap()))
}

fn read_u64(bytes: &[u8], offset: &mut usize) -> Result<u64, IndexIoError> {
    let end = *offset + 8;
    let buf = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("unexpected EOF".to_string()))?;
    *offset = end;
    Ok(u64::from_le_bytes(buf.try_into().unwrap()))
}

fn read_u64_infallible(bytes: &[u8], offset: &mut usize) -> u64 {
    let v = u64::from_le_bytes(bytes[*offset..*offset + 8].try_into().unwrap());
    *offset += 8;
    v
}

fn read_i64_infallible(bytes: &[u8], offset: &mut usize) -> i64 {
    let v = i64::from_le_bytes(bytes[*offset..*offset + 8].try_into().unwrap());
    *offset += 8;
    v
}
```

- [ ] **Step 2: Wire the view into `ReferenceIndex` and export the types.**

In `src/genomics/index/io.rs`, add view constructors to `impl ReferenceIndex` (after `contigs()`):

```rust
    /// Borrow a zero-copy FM-index view over the memory-mapped sections.
    pub fn view(&self) -> Result<crate::genomics::index::view::FmIndexView<'_>, IndexIoError> {
        crate::genomics::index::view::FmIndexView::new(self.mmap.as_bytes(), &self.sections)
    }

    /// Borrow a zero-copy genome view (FM-index + contig set) for `locate_exact`.
    pub fn genome_view(
        &self,
    ) -> Result<crate::genomics::index::view::GenomeIndexView<'_>, IndexIoError> {
        let fm = self.view()?;
        Ok(crate::genomics::index::view::GenomeIndexView::new(
            fm,
            &self.contigs,
        ))
    }
```

`IndexIoError` must be reachable from `view.rs`; it is `pub` in `io.rs`. Add to `io.rs`'s module the visibility the view needs: `IndexIoError` is already `pub`. The view imports `crate::genomics::index::io::IndexIoError` — ensure `io` is reachable: in `src/genomics/index/mod.rs` the modules are private (`mod io; mod format;`). A sibling module (`view`) can refer to `crate::genomics::index::io::IndexIoError` because `io` is visible to its siblings. Confirm by compiling.

In `src/genomics/index/mod.rs`, register and export the view:

```rust
mod format;
mod io;
mod view;

pub use format::IndexHeader;
pub use io::{IndexReader, IndexWriter, ReferenceIndex};
pub use view::{FmIndexView, GenomeIndexView};
```

In `src/genomics/mod.rs`, extend the index re-export:

```rust
pub use index::{FmIndexView, GenomeIndexView, IndexHeader, IndexReader, IndexWriter, ReferenceIndex};
```

Also make the format module's `SectionEntry`/`SectionKind` reachable from `view.rs`: they are `pub` in `format.rs`, and `format` is visible to its sibling `view` via `crate::genomics::index::format::{SectionEntry, SectionKind}`. Confirm by compiling.

- [ ] **Step 3: Add the equivalence gate (failing until the view compiles + is correct).** Append to the `#[cfg(test)] mod tests` of `src/genomics/index/io.rs`:

```rust
    #[test]
    fn view_is_byte_identical_to_in_ram_index() {
        let idx = sample_index();
        let path = temp_path("equiv");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
        let loaded = IndexReader::open(&path).unwrap();
        let gv = loaded.genome_view().unwrap();
        let fv = gv.fm();
        let fm = idx.fm();

        // backward_search + sa_at over a pattern battery, incl. boundary-straddle
        // and an N-bearing pattern.
        let patterns: &[&[u8]] = &[
            b"A", b"C", b"G", b"T", b"N", b"AC", b"ACG", b"ACGT", b"GT", b"TTTT",
            b"GGGG", b"CCCC", b"AAAA", b"ACGTACG", b"NNAC", b"GTTTTT", b"NNNN",
            b"TACG", b"CGTACG", b"GTACGT",
        ];
        for &p in patterns {
            let a = fm.backward_search(p);
            let b = fv.backward_search(p);
            assert_eq!(a.lower, b.lower, "lower mismatch for {p:?}");
            assert_eq!(a.upper, b.upper, "upper mismatch for {p:?}");
            // sa_at over the whole interval must agree.
            for i in (a.lower as usize)..(a.upper as usize) {
                assert_eq!(fm.sa_at(i), fv.sa_at(i), "sa_at mismatch @ {i} for {p:?}");
            }
            // locate_exact returns identical Locus vectors.
            assert_eq!(
                idx.locate_exact(p, 1024),
                gv.locate_exact(p, 1024),
                "locate_exact mismatch for {p:?}"
            );
        }

        // Exhaustive rank/symbol_at equivalence over every BWT position.
        use crate::genomics::fm_backing::{self, BwtBacking};
        for symbol in [
            FmSymbol::Sentinel,
            FmSymbol::Base(BaseCode::A),
            FmSymbol::Base(BaseCode::C),
            FmSymbol::Base(BaseCode::G),
            FmSymbol::Base(BaseCode::T),
            FmSymbol::Base(BaseCode::N),
        ] {
            for pos in 0..=fm.len() {
                assert_eq!(
                    fm.rank(symbol, pos),
                    fm_backing::rank(fv, symbol, pos),
                    "rank mismatch {symbol:?} @ {pos}"
                );
            }
        }
        for i in 0..fm.len() {
            assert_eq!(fm.symbol_at(i), fm_backing::symbol_at(fv, i), "symbol_at @ {i}");
        }
        let _ = std::fs::remove_file(path);
    }
```

This test needs `BaseCode` in scope — add `use crate::genomics::BaseCode;` and `use crate::genomics::FmSymbol;` to the test module's `use super::*;` imports if not already pulled in (they are re-exported from `crate::genomics`; add explicit `use` lines in the test module).

- [ ] **Step 4: Build and run the equivalence gate + full lib suite.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --lib index 2>&1 | tail -30`
Expected: `view_is_byte_identical_to_in_ram_index` PASS, plus all Task 5/6 io tests + the format test. Then `cargo test --lib 2>&1 | grep -E "test result:"` (green).

- [ ] **Step 5: Commit.**

```bash
git add src/genomics/index/view.rs src/genomics/index/io.rs src/genomics/index/mod.rs src/genomics/mod.rs
git commit -m "feat(genomics/index/view): zero-copy FmIndexView + GenomeIndexView + equivalence gate" \
  -m "FmIndexView reads the FM sections as &[u64]/&[u32] over the mmap (checked align_to, little-endian host) and implements BwtBacking, so backward_search/sa_at/locate_interval run the same generic ops as the in-RAM index. GenomeIndexView pairs it with the ContigSet for locate_exact. The equivalence gate asserts byte-identical backward_search/sa_at/locate_exact/rank/symbol_at over a multi-contig + N-bearing battery." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: integration gates (`tests/index_persistence.rs`) + `docs/index-format.md`

Promote the gates to a public integration test (querying only through the crate's public API, proving the artifact is self-contained), add the determinism / no-rebuild / bounded-residency gates, and publish the format contract.

**Files:**
- Create: `tests/index_persistence.rs`, `docs/index-format.md`

- [ ] **Step 1: Write the integration test.** Create `tests/index_persistence.rs`:

```rust
//! Phase B3b gates: the persisted, memory-mapped FM-index is byte-identical to
//! the in-RAM index, deterministic, self-contained, loads without rebuilding,
//! and keeps the index out of resident memory.

use std::env;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::core::Locus;
use rosalind::genomics::{GenomeIndex, IndexReader, IndexWriter};

fn temp_path(suffix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    env::temp_dir().join(format!("rosalind-b3b-it-{suffix}-{nanos}.idx"))
}

fn multi_contig_index() -> GenomeIndex {
    GenomeIndex::from_named_sequences(&[
        ("chr1".to_string(), b"ACGTACGTNNACGTACGTACGTAAGGCCTT".to_vec()),
        ("chr2".to_string(), b"TTTTGGGGCCCCAAAANNNNACGTACGTAC".to_vec()),
        ("chr3".to_string(), b"GATTACAGATTACANNNNGATTACAGGGGG".to_vec()),
    ])
    .expect("build")
}

#[test]
fn view_equals_in_ram_over_a_pattern_battery() {
    let idx = multi_contig_index();
    let path = temp_path("equiv");
    IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
    let loaded = IndexReader::open(&path).unwrap();
    let gv = loaded.genome_view().unwrap();

    let patterns: &[&[u8]] = &[
        b"A", b"C", b"G", b"T", b"N", b"GATTACA", b"ACGT", b"NNNN", b"GGGGG",
        b"TTTTGGGG", b"AAGGCCTT", b"CGTACGTAC", b"GATTACAGGGGG", b"acgt", b"NnAc",
    ];
    for &p in patterns {
        let expected: Vec<Locus> = idx.locate_exact(p, 4096);
        let got: Vec<Locus> = gv.locate_exact(p, 4096);
        assert_eq!(got, expected, "locate_exact mismatch for {p:?}");
    }
    let _ = std::fs::remove_file(path);
}

#[test]
fn build_is_deterministic() {
    // Two independent builds of the same reference (SA-IS build + serialize) must
    // produce byte-identical files — the end-to-end determinism gate.
    let p1 = temp_path("det1");
    let p2 = temp_path("det2");
    IndexWriter::create(&p1)
        .unwrap()
        .write_genome_index(&multi_contig_index())
        .unwrap();
    IndexWriter::create(&p2)
        .unwrap()
        .write_genome_index(&multi_contig_index())
        .unwrap();
    assert_eq!(std::fs::read(&p1).unwrap(), std::fs::read(&p2).unwrap());
    let _ = std::fs::remove_file(p1);
    let _ = std::fs::remove_file(p2);
}

#[test]
fn index_is_self_contained() {
    // Build, drop the in-RAM index, then query using only the file.
    let path = temp_path("selfcontained");
    {
        let idx = multi_contig_index();
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
    } // `idx` dropped here — the file is the only source.

    let loaded = IndexReader::open(&path).unwrap();
    let gv = loaded.genome_view().unwrap();
    // "GATTACA" occurs in chr3 at pos 0 and pos 18.
    let loci = gv.locate_exact(b"GATTACA", 4096);
    assert_eq!(loci.len(), 2, "expected two GATTACA hits in chr3");
    assert!(loci.iter().all(|l| l.contig == 2));
    let _ = std::fs::remove_file(path);
}

#[test]
fn view_is_a_small_borrow_not_an_owned_copy() {
    // Bounded residency (B3b witness; the enforced RSS gate is Phase C): the view
    // is a handful of slices + scalars, not an owned copy of the index. Its size
    // is independent of genome size.
    use rosalind::genomics::FmIndexView;
    assert!(
        std::mem::size_of::<FmIndexView<'_>>() <= 256,
        "FmIndexView must be a small borrow, got {} bytes",
        std::mem::size_of::<FmIndexView<'_>>()
    );
}
```

Note on `no-rebuild`: it is **enforced structurally** — `IndexReader::open` and `ReferenceIndex::view` never construct a `BlockedFMIndex` and never call `sais_u32` (the only SA-IS entry point); the read path is `mmap` + slice. `index_is_self_contained` is the behavioral witness (querying works after the in-RAM index is dropped, with no rebuild possible). If desired, a reviewer can confirm `grep -n "BlockedFMIndex::build\|sais_u32" src/genomics/index/` returns nothing.

- [ ] **Step 2: Run the integration test.**

Run: `cargo test --test index_persistence 2>&1 | tail -25`
Expected: all four tests PASS.

- [ ] **Step 3: Write the format documentation.** Create `docs/index-format.md`:

```markdown
# Rosalind on-disk index format (v1)

A single-file, little-endian, section-based artifact holding a built
`GenomeIndex`: the contig table, the 2-bit forward reference, and the blocked
FM-index (boundaries + per-block BWT/occ + compact sampled suffix array). It is
designed to be **memory-mapped and queried in place** — the reader never rebuilds
the index and never allocates its bulk arrays (`FmIndexView` borrows mmap slices).

> **Stability:** v1 is pre-1.0 and may change until the 1.0 contract (see
> `docs/OPEN_PROBLEMS.md`). The `IndexVersion` is bumped on any layout change.

## File layout

```
[ fixed header ]
[ Contigs ]        (8-aligned)
[ Reference2bit ]  (8-aligned)
[ FmMeta ]         (8-aligned)
[ Boundaries ]     (8-aligned)
[ Blocks ]         (8-aligned)
[ SaSamples ]      (8-aligned)
[ section table ]  (8-aligned)
```

Every section starts at an 8-byte boundary (zero-padded). All integers are
little-endian. Multi-byte arrays (`[u64]`, `[u32]`) are stored as LE words at
8-aligned offsets so the reader reinterprets them via a checked `slice::align_to`
— zero-copy, dependency-free. The reader requires a little-endian host.

## Header (fixed size)

`magic "ROSALIND"`, `version:u16 (=1)`, `endian:u8 (=1, little)`, `reserved0:u8`,
`flags:u32`, `contig_count:u32`, `sa_sample_rate:u32`, `header_bytes:u64`,
`section_table_offset:u64`, `section_table_bytes:u64`, `reference_blake3:[u8;32]`.
`reference_blake3` is the BLAKE3 of the uppercased ASCII reference — the identity
a consumer checks to confirm an index matches its reference (the B4 stale-index
guard).

## Sections

| Section | Payload |
|---|---|
| `Contigs` | per contig: `name_len:u32, name:[u8], length:u64, global_offset:u64` (read by byte-copy; recomputed offsets are integrity-checked against the stored value) |
| `Reference2bit` | `len:u64, data_words:u64, amb_words:u64`, then `data:[u64]` (2-bit, 32 bases/word) + ambiguity `bits:[u64]` (1 bit/base marks `N`). Self-contained ref-base lookups (B4). |
| `FmMeta` | `block_size:u64, bwt_len:u64, sentinel_pos:u64, sa_sample_rate:u64, num_blocks:u64, c_table:[u32;6]` |
| `Boundaries` | `num_blocks + 1` entries of `cumulative_counts:[u32;5], sentinel_count:u32` |
| `Blocks` | `[directory: u64; num_blocks]` (record offsets relative to the section start), then per block: `start:u64, end:u64, sentinel_offset:i64 (-1 = none), stride:u64, bwt_data_words:u64, bwt_amb_words:u64, occ_bitvec_words:u64, occ_superblock_len:u64`, then `bwt_data:[u64]`, `bwt_amb:[u64]`, `occ_bitvecs:[[u64]; 5]`, `occ_superblocks:[[u32]; 5]`, `totals:[u32; 5]` (record padded to 8). The directory makes block `k` O(1) to locate. |
| `SaSamples` | `rate:u64, bwt_len:u64, marks_words:u64, superblocks_len:u64, values_len:u64`, then `marks:[u64]`, `superblocks:[u32]`, `values:[u32]` (the compact `SampledSuffixArray`; superblock stride is the fixed `RANK_STRIDE`). |

## Determinism

SA-IS is deterministic; sections are written in fixed order with fixed-width LE
encoding and zero padding; no timestamps appear. ⇒ the file is **byte-identical
across repeated builds of the same genome** (gated by `build_is_deterministic`).

## Faithful serialization (size) and the v2 lean path

v1 serialises the blocked structure **faithfully**: both the per-block 2-bit BWT
(for `block_symbol`) and the five occ rank bitvectors (for `block_rank`), ≈7
bits/symbol. This is the honest cost of "one algorithm, two backings" — and RSS
stays bounded regardless, because the index is mmap'd / OS-paged (file size costs
disk, not resident memory). A future **format v2** may derive `block_symbol` from
the rank bitvectors (or rank over the 2-bit BWT) and drop the redundant copy — a
query-speed/size optimization behind a clean `IndexVersion` bump, not a
correctness change.
```

- [ ] **Step 4: Final verification (whole stage).**

Run, and confirm each:
- `cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` — full suite green (report totals).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo clippy --lib 2>&1 | grep -E 'fm_backing|fm_index|index/(io|view|format)'` — no new lints in the touched files.
- `cargo fmt --all -- --check` — clean.
- `grep -rn "BlockedFMIndex::build\|sais_u32" src/genomics/index/` — **empty** (the read path never rebuilds).
- `grep -n "write_v1\|struct ContigInfo" src/genomics/index/io.rs` — **empty** (the scaffold is gone).

- [ ] **Step 5: Commit.**

```bash
git add tests/index_persistence.rs docs/index-format.md
git commit -m "test(index): B3b persistence gates + docs/index-format.md" \
  -m "Integration gates over the public API: view == in-RAM locate_exact over a multi-contig + N battery, byte-identical builds, self-contained query (in-RAM index dropped), and bounded residency (the view is a small borrow, not an owned copy). Documents the versioned on-disk format as a forkable contract, incl. the faithful-serialization size note and the v2 lean-rank path." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before the B3b PR)

- `cargo test` — full suite green. Key witnesses: **B3b.1** — `tests/fm_index_props`, the FM-index/aligner/`genome_index` unit suites (behavior unchanged through the `BwtBacking` refactor); **B3b.2** — `view_is_byte_identical_to_in_ram_index` (unit) and `tests/index_persistence` (equivalence battery, determinism, self-contained, bounded residency).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo clippy --lib` — no new lints in touched files; pay attention to the `unsafe` `align_to` blocks (documented invariants) and `#[allow]`-free imports.
- `cargo fmt --all -- --check` — clean.
- Structural no-rebuild check: `grep -rn "BlockedFMIndex::build\|sais_u32" src/genomics/index/` is empty.

## Self-Review

- **Spec coverage (`2026-05-27-phase-b3b-zerocopy-index-design.md`):**
  - §3 `BwtBacking` seam (one algorithm, two backings) ✔ Tasks 1–2 (trait + generic ops + owned impl) and Task 7 (borrowed impl). The trait adds `num_blocks()`/`sample_rate()` to the spec's listed surface — a deliberate elaboration so the `rank` terminal-block guard and the `sa_at` debug guard are exact (both backings have these natively); noted here as the one refinement of the spec's illustrative list.
  - §4 on-disk format (8-aligned LE, dependency-free `align_to`, faithful v1, BLAKE3) ✔ Tasks 4–5 + `docs/index-format.md` (Task 8). The fixed header is retained; `SectionKind` extended; `Reference2bit` self-contained; v2 lean path documented.
  - §5 `FmIndexView` (mmap → `BwtBacking` → query) + `GenomeIndexView::locate_exact` ✔ Task 7.
  - §6 decomposition (B3b.1 then B3b.2) ✔ Part 1 = Tasks 1–2; Part 2 = Tasks 3–8.
  - §7 gates — equivalence ✔ (Task 7 unit + Task 8 integration battery, multi-contig + boundary-straddle + N), determinism ✔ (Task 5 + Task 8), no rebuild on load ✔ (structural: read path never calls `BlockedFMIndex::build`/`sais_u32`; behavioral: `index_is_self_contained`), bounded residency ✔ (Task 8 `view_is_a_small_borrow_not_an_owned_copy` + the mmap design; enforced RSS gate is Phase C), behavior preserved ✔ (Task 2 witnesses).
  - §8 testing (equivalence battery, determinism, corrupt/truncated rejection, self-contained round-trip) ✔ Tasks 5–8.
- **Type/name consistency:** `BwtBacking` method set and the `fm_backing::*` op signatures are identical in Task 1 (definition), Task 2 (owned impl + wrappers), and Task 7 (borrowed impl + wrappers). Section kinds (`Reference2bit`/`FmMeta`/`Boundaries`/`Blocks`) defined in Task 4 are written in Task 5 and read in Tasks 6–7 with the exact field order documented in Task 8. The serializer's per-block field order (`start,end,sentinel_offset,stride,bwt_data_words,bwt_amb_words,occ_bitvec_words,occ_superblock_len` then payload) matches the view's `block()` parse exactly. `write_genome_index`/`open`/`view`/`genome_view`/`locate_exact` names are used identically across tasks.
- **No placeholders:** every step ships complete code or an exact command + expected output. The one transient `unimplemented!()` (Task 5 reader stub) is a deliberate compile-bridge replaced in Task 6's same-file rewrite; flagged explicitly.
- **MSRV 1.72:** all ceil-divisions use `(a + b - 1) / b` (never `div_ceil`); `align_to`, `BufWriter`, `cfg!(target_endian)` are all ≤1.72.
- **Determinism & safety:** fixed section order, zero padding, no timestamps ⇒ byte-identical builds (gated). The two `unsafe align_to` blocks carry documented invariants (empty prefix/suffix enforced at runtime; LE host enforced in `FmIndexView::new`; `u64`/`u32` have no invalid bit patterns).
- **Scope:** no `rosalind index` CLI (B3c), no `align`/`variants` wiring (B4), no `MemoryBudget` enforcement (Phase C), no lean-rank (format v2) — each named as deferred.
```
