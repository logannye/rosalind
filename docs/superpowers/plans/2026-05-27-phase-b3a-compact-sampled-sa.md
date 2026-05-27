# Phase B3a — Compact sampled suffix array Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the FM-index's dense sampled-suffix-array (`sa_samples: Vec<u32>`, one `u32` per BWT position ≈ 12 GB for a human genome) with a compact sparse representation — a sampled-position bitvector + a packed array of sampled values — making sample storage `O(bwt_len bits + bwt_len/rate values)` instead of `O(bwt_len)` words, with `sa_at` behavior byte-for-byte unchanged.

**Architecture:** A new `SampledSuffixArray` type owns a bitvector `marks` (bit *i* set iff BWT position *i* carries a sampled SA value), superblock prefix-popcounts over `marks` for `O(1)` rank, and a packed `values` array in ascending BWT-position order. `sample_at(i)` returns `Some(values[rank(i)])` when bit *i* is set, else `None`. `BlockedFMIndex` stores one `SampledSuffixArray` instead of the dense vec; `sa_at` walks LF until `sample_at` hits. This is the prerequisite for a bounded-memory persisted index (Phase B3b serialises exactly this structure).

**Tech Stack:** Rust 2021; `u64`-word bitvector + `count_ones` popcount (mirrors the existing `genomics::rank_select` pattern); no new dependencies.

This is stage **B3a** of Phase B's index-persistence work (`docs/superpowers/specs/2026-05-27-phase-b-genome-scale-design.md` §3.6, and the contract thesis in `docs/OPEN_PROBLEMS.md`). It lands green and is independently valuable (it slashes the FM-index's in-RAM footprint). It does **not** add persistence or the zero-copy view (B3b) or `rosalind index` (B3c).

---

## File structure

- `src/genomics/sampled_sa.rs` — **Create.** `SampledSuffixArray` (`marks`/`superblocks`/`values`/`rate`/`bwt_len`), `from_sorted_samples`, `sample_at`, getters, and private `popcount_prefix`/`build_superblocks` helpers. Single responsibility: the compact sampled-SA structure + its rank.
- `src/genomics/mod.rs` — **Modify.** `mod sampled_sa;` + `pub use sampled_sa::SampledSuffixArray;`.
- `src/genomics/fm_index.rs` — **Modify.** `build_bwt_and_sa_samples` returns a `SampledSuffixArray`; `BlockedFMIndex` stores `sampled: SampledSuffixArray` (replacing `sa_sample_rate: usize` + `sa_samples: Vec<u32>`); `sa_at` and `sa_sample_rate()` use it; add a `sampled()` accessor.

---

## Task 1: `SampledSuffixArray` (compact sampled-SA structure)

**Files:**
- Create: `src/genomics/sampled_sa.rs`
- Modify: `src/genomics/mod.rs`

- [ ] **Step 1: Create the module with the type, signatures (unimplemented), and tests.** Create `src/genomics/sampled_sa.rs`:

```rust
//! Compact sampled suffix array for the FM-index.
//!
//! Replaces a dense `Vec<u32>` (one slot per BWT position, ~12 GB for a human
//! genome) with a sparse representation: a sampled-position **bitvector**
//! (`marks`), superblock prefix-popcounts over it for O(1) rank, and a **packed**
//! array of the sampled SA values in ascending BWT-position order. Storage is
//! `O(bwt_len bits + (bwt_len / rate) values)`. This is the structure the on-disk
//! index serialises (Phase B3b).

/// BWT positions per superblock in the rank structure over `marks`.
const RANK_STRIDE: usize = 1024;

/// A compact sampled suffix array.
#[derive(Debug, Clone)]
pub struct SampledSuffixArray {
    rate: usize,
    bwt_len: usize,
    /// Bit `i` set iff BWT position `i` carries a sampled SA value.
    marks: Vec<u64>,
    /// Prefix popcount of `marks` at each `RANK_STRIDE` boundary (`len/stride + 1`).
    superblocks: Vec<u32>,
    /// Sampled SA values, in ascending BWT-position order.
    values: Vec<u32>,
}

impl SampledSuffixArray {
    /// Build from `(bwt_position, sa_value)` pairs of the sampled positions,
    /// which MUST be yielded in **strictly ascending** `bwt_position` order.
    pub fn from_sorted_samples(
        bwt_len: usize,
        rate: usize,
        samples: impl Iterator<Item = (usize, u32)>,
    ) -> Self {
        unimplemented!()
    }

    /// The sampling rate (`>= 1`).
    pub fn rate(&self) -> usize {
        self.rate
    }

    /// Number of BWT positions covered.
    pub fn len(&self) -> usize {
        self.bwt_len
    }

    /// Whether the structure covers zero positions.
    pub fn is_empty(&self) -> bool {
        self.bwt_len == 0
    }

    /// Number of sampled values stored (≈ `bwt_len / rate`).
    pub fn num_samples(&self) -> usize {
        self.values.len()
    }

    /// The sampled SA value at BWT position `bwt_idx`, or `None` if that position
    /// is not sampled. `O(1)` via superblock + word popcount.
    pub fn sample_at(&self, bwt_idx: usize) -> Option<u32> {
        unimplemented!()
    }
}

/// Popcount of set bits in `[start, end)` of the `u64`-word bitvector `words`.
fn popcount_prefix(words: &[u64], start: usize, end: usize) -> u32 {
    if end <= start {
        return 0;
    }
    let start_word = start / 64;
    let end_word = (end - 1) / 64;
    let start_bit = start % 64;
    let end_bit = end % 64;

    if start_word == end_word {
        let mut w = words[start_word];
        w &= !((1u64 << start_bit) - 1);
        if end_bit != 0 {
            w &= (1u64 << end_bit) - 1;
        }
        return w.count_ones();
    }

    let mut count = (words[start_word] & !((1u64 << start_bit) - 1)).count_ones();
    for w in &words[start_word + 1..end_word] {
        count += w.count_ones();
    }
    let mut last = words[end_word];
    if end_bit != 0 {
        last &= (1u64 << end_bit) - 1;
    }
    count + last.count_ones()
}

/// Prefix popcounts of `marks` at each `stride` boundary (entry `k` = set bits in
/// `[0, k*stride)`); length `ceil(bwt_len/stride) + 1`.
fn build_superblocks(marks: &[u64], bwt_len: usize, stride: usize) -> Vec<u32> {
    let num = (bwt_len + stride - 1) / stride;
    let mut sb = Vec::with_capacity(num + 1);
    sb.push(0u32);
    let mut acc = 0u32;
    for i in 0..num {
        let start = i * stride;
        let end = ((i + 1) * stride).min(bwt_len);
        acc += popcount_prefix(marks, start, end);
        sb.push(acc);
    }
    sb
}

#[cfg(test)]
mod tests {
    use super::*;

    // bwt_len 10, rate 3; sample positions 0,3,6,9 with SA values 9,6,3,0.
    fn fixture() -> SampledSuffixArray {
        SampledSuffixArray::from_sorted_samples(
            10,
            3,
            [(0usize, 9u32), (3, 6), (6, 3), (9, 0)].into_iter(),
        )
    }

    #[test]
    fn sample_at_returns_value_at_sampled_positions() {
        let s = fixture();
        assert_eq!(s.sample_at(0), Some(9));
        assert_eq!(s.sample_at(3), Some(6));
        assert_eq!(s.sample_at(6), Some(3));
        assert_eq!(s.sample_at(9), Some(0));
    }

    #[test]
    fn sample_at_returns_none_at_unsampled_positions() {
        let s = fixture();
        for i in [1usize, 2, 4, 5, 7, 8] {
            assert_eq!(s.sample_at(i), None, "position {i} should be unsampled");
        }
    }

    #[test]
    fn metadata_is_correct() {
        let s = fixture();
        assert_eq!(s.rate(), 3);
        assert_eq!(s.len(), 10);
        assert_eq!(s.num_samples(), 4);
        assert!(!s.is_empty());
    }

    #[test]
    fn handles_a_superblock_boundary() {
        // > RANK_STRIDE positions, sampled sparsely across the boundary.
        let n = RANK_STRIDE * 2 + 5;
        let samples: Vec<(usize, u32)> = (0..n).step_by(64).map(|p| (p, p as u32)).collect();
        let s = SampledSuffixArray::from_sorted_samples(n, 64, samples.iter().copied());
        for &(p, v) in &samples {
            assert_eq!(s.sample_at(p), Some(v));
        }
        assert_eq!(s.sample_at(1), None);
        assert_eq!(s.sample_at(RANK_STRIDE + 1), None);
    }

    #[test]
    fn empty_has_no_samples() {
        let s = SampledSuffixArray::from_sorted_samples(0, 1, std::iter::empty());
        assert!(s.is_empty());
        assert_eq!(s.num_samples(), 0);
    }
}
```

- [ ] **Step 2: Register the module.** In `src/genomics/mod.rs`, add `mod sampled_sa;` (with the other `mod` declarations, alphabetically after `mod rank_select;`) and `pub use sampled_sa::SampledSuffixArray;` (with the other `pub use` lines).

- [ ] **Step 3: Run tests to verify they fail.**

Run: `cargo test --lib sampled_sa 2>&1 | tail -20`
Expected: tests run and FAIL — `not implemented` from the two `unimplemented!()` bodies.

- [ ] **Step 4: Implement `from_sorted_samples` and `sample_at`.** Replace the two `unimplemented!()` bodies:

```rust
    pub fn from_sorted_samples(
        bwt_len: usize,
        rate: usize,
        samples: impl Iterator<Item = (usize, u32)>,
    ) -> Self {
        let words = (bwt_len + 63) / 64;
        let mut marks = vec![0u64; words];
        let mut values = Vec::new();
        for (pos, val) in samples {
            debug_assert!(pos < bwt_len, "sample position out of range");
            marks[pos / 64] |= 1u64 << (pos % 64);
            values.push(val);
        }
        let superblocks = build_superblocks(&marks, bwt_len, RANK_STRIDE);
        Self {
            rate: rate.max(1),
            bwt_len,
            marks,
            superblocks,
            values,
        }
    }
```

```rust
    pub fn sample_at(&self, bwt_idx: usize) -> Option<u32> {
        debug_assert!(bwt_idx < self.bwt_len, "bwt_idx out of range");
        let word = self.marks[bwt_idx / 64];
        if word & (1u64 << (bwt_idx % 64)) == 0 {
            return None;
        }
        // rank = number of set bits strictly before bwt_idx = index into `values`.
        let sb = bwt_idx / RANK_STRIDE;
        let rank = self.superblocks[sb] + popcount_prefix(&self.marks, sb * RANK_STRIDE, bwt_idx);
        Some(self.values[rank as usize])
    }
```

- [ ] **Step 5: Run tests to verify they pass.**

Run: `cargo test --lib sampled_sa 2>&1 | tail -20`
Expected: `test result: ok. 5 passed`. Then `cargo build --lib 2>&1 | grep -i warning` (none) and `cargo fmt --all` (apply formatting).

- [ ] **Step 6: Commit.**

```bash
git add src/genomics/sampled_sa.rs src/genomics/mod.rs
git commit -m "feat(genomics/sampled_sa): compact sampled suffix array (bitvector + packed values)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Wire `SampledSuffixArray` into `BlockedFMIndex` (replace the dense vec)

**Files:**
- Modify: `src/genomics/fm_index.rs`

The current `BlockedFMIndex` holds `sa_sample_rate: usize` and `sa_samples: Vec<u32>` (dense, `u32::MAX` = unsampled), and `build_bwt_and_sa_samples` returns that dense vec. Replace both with a single `SampledSuffixArray`, preserving `sa_at` behavior exactly.

- [ ] **Step 1: Add a sparsity-regression test (failing).** In the `#[cfg(test)] mod tests` of `src/genomics/fm_index.rs`, add:

```rust
    #[test]
    fn sampled_sa_is_sparse_not_dense() {
        // A reference long enough that a dense (one-slot-per-position) sample
        // array would be ~bwt_len; the compact structure stores ~bwt_len/rate.
        let reference = vec![b'A'; 4096];
        let index = BlockedFMIndex::build(&reference, 64).expect("build");
        let sampled = index.sampled();
        assert_eq!(sampled.len(), reference.len() + 1); // covers every BWT position
        // Far fewer stored values than positions (sampled at `rate`).
        assert!(
            sampled.num_samples() <= sampled.len() / sampled.rate() + 1,
            "expected ~len/rate samples, got {} for len {} rate {}",
            sampled.num_samples(),
            sampled.len(),
            sampled.rate(),
        );
        assert!(sampled.num_samples() * 4 < sampled.len());
    }
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --lib fm_index::tests::sampled_sa_is_sparse_not_dense 2>&1 | tail -20`
Expected: compile error — `BlockedFMIndex::sampled()` does not exist.

- [ ] **Step 3: Change the struct fields.** In `src/genomics/fm_index.rs`, add the import near the top (with the other `use crate::genomics::...` lines):

```rust
use crate::genomics::SampledSuffixArray;
```

Replace the two fields in `struct BlockedFMIndex`:

```rust
    sa_sample_rate: usize,
    sa_samples: Vec<u32>, // u32::MAX = not sampled
```

with:

```rust
    sampled: SampledSuffixArray,
```

- [ ] **Step 4: Make `build_bwt_and_sa_samples` return a `SampledSuffixArray`.** Change its signature and the sampling tail. The signature becomes:

```rust
fn build_bwt_and_sa_samples(
    reference: &[u8],
    sa_sample_rate: usize,
) -> Result<(Vec<u8>, usize, SampledSuffixArray), FMIndexError> {
```

Inside, **remove** the dense `let mut sa_samples = vec![u32::MAX; text.len()];` line and the `if sa_sample_rate > 0 && (sa_idx % sa_sample_rate == 0) { sa_samples[bwt_idx] = sa_idx_u32; }` block in the BWT loop. After the loop (which still produces `bwt` and `sentinel_pos`), build the compact structure from a lazy pass over `sa` and return it:

```rust
    let rate = sa_sample_rate.max(1);
    let sampled = SampledSuffixArray::from_sorted_samples(
        text.len(),
        rate,
        sa.iter()
            .enumerate()
            .filter_map(|(bwt_idx, &sa_idx)| ((sa_idx as usize) % rate == 0).then_some((bwt_idx, sa_idx))),
    );

    Ok((bwt, sentinel_pos, sampled))
```

(`enumerate()` yields BWT positions in ascending order, satisfying `from_sorted_samples`' ordering contract.)

- [ ] **Step 5: Update `build`, `sa_at`, the accessor.** In `BlockedFMIndex::build`, the call site is `let (bwt, sentinel_pos, sa_samples) = build_bwt_and_sa_samples(&clean, sa_sample_rate)?;` — rename the binding to `sampled`:

```rust
        let (bwt, sentinel_pos, sampled) = build_bwt_and_sa_samples(&clean, sa_sample_rate)?;
```

and in the returned `Ok(Self { … })`, replace the `sa_sample_rate,` and `sa_samples,` fields with `sampled,`.

Replace `sa_sample_rate()` to read from `sampled`:

```rust
    /// Suffix array sample rate used for locating.
    pub fn sa_sample_rate(&self) -> usize {
        self.sampled.rate()
    }
```

Add an accessor (e.g. just after `sa_sample_rate`):

```rust
    /// The compact sampled suffix array backing `sa_at`.
    pub fn sampled(&self) -> &SampledSuffixArray {
        &self.sampled
    }
```

Replace the lookup in `sa_at` — change `let sampled = self.sa_samples[current]; if sampled != u32::MAX { return sampled as usize + lf_steps; }` to use `sample_at`:

```rust
    pub fn sa_at(&self, index: usize) -> usize {
        assert!(index < self.bwt_len, "BWT index out of range");

        let mut current = index;
        let mut lf_steps = 0usize;

        loop {
            if let Some(sampled) = self.sampled.sample_at(current) {
                return sampled as usize + lf_steps;
            }
            current = self.lf_index(current);
            lf_steps += 1;
            debug_assert!(
                lf_steps <= self.sampled.rate() + 1,
                "LF steps exceeded sample rate; sampling invariant violated"
            );
        }
    }
```

- [ ] **Step 6: Build and run the focused + full FM-index/aligner suites.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --lib fm_index 2>&1 | tail -20` (the existing `fm_index_builds_and_ranks`, `sa_at_recovers_reference_position`, plus the new `sampled_sa_is_sparse_not_dense` all pass), then
`cargo test --lib genomics 2>&1 | tail -15` (the aligner + `genome_index` `locate_exact` tests — which exercise `sa_at` end-to-end — stay green, confirming `sa_at` behavior is unchanged).

- [ ] **Step 7: Full suite + format.**

Run: `cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` (all green; report totals), `cargo fmt --all -- --check` (clean — run `cargo fmt --all` first if needed), `cargo build 2>&1 | grep -ic warning` (0).

- [ ] **Step 8: Commit.**

```bash
git add src/genomics/fm_index.rs
git commit -m "refactor(genomics/fm_index): back sa_at with the compact SampledSuffixArray" \
  -m "Replaces the dense sa_samples: Vec<u32> (one slot per BWT position, ~12 GB for a human genome) with a SampledSuffixArray (sampled-position bitvector + packed values, ~bwt_len/rate). sa_at behavior is unchanged; the aligner and genome_index locate tests confirm it end-to-end." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before the B3a PR)

- `cargo test` — full suite green (report totals; the FM-index, aligner, and `genome_index` suites are the key witnesses that `sa_at` is unchanged).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo clippy --lib 2>&1 | grep -E 'sampled_sa|fm_index'` — no new lints in the touched files.
- `cargo fmt --all -- --check` — clean.
- `grep -n "sa_samples\|u32::MAX" src/genomics/fm_index.rs` — the dense `sa_samples` vec and its `u32::MAX` sentinel are gone (no dense per-position array remains).

## Self-Review

- **Spec coverage (§3.6 of the Phase B design + the contract thesis):** dense `Vec<u32>` → compact sparse (bitvector + packed values) ✔ (Task 1 `SampledSuffixArray`, Task 2 wiring); sub-linear sample storage ✔ (`sampled_sa_is_sparse_not_dense` asserts `num_samples*4 < len`); `sa_at` behavior preserved ✔ (existing fm_index/aligner/genome_index tests are the witnesses). Build-side moved to the compact form ✔ (`build_bwt_and_sa_samples` no longer allocates the dense array).
- **Scope:** no persistence, no zero-copy view, no `rosalind index` (those are B3b/B3c). No public-API change beyond the additive `sampled()` accessor (`sa_sample_rate()` semantics preserved).
- **Type consistency:** `SampledSuffixArray::from_sorted_samples(bwt_len: usize, rate: usize, impl Iterator<Item=(usize,u32)>)`, `sample_at(usize) -> Option<u32>`, `rate()/len()/num_samples()/is_empty()` are used identically in Task 1's tests, Task 2's wiring, and the sparsity test. `build_bwt_and_sa_samples` returns `(Vec<u8>, usize, SampledSuffixArray)` and `build` binds `sampled` accordingly.
- **No placeholders:** every step ships complete code or an exact command + expected output; the `unimplemented!()` bodies are a deliberate TDD red-state replaced in the same task.
- **Ordering contract:** `from_sorted_samples` requires ascending BWT positions; the only producer (`build_bwt_and_sa_samples` via `sa.iter().enumerate()`) satisfies it; documented + `debug_assert`-guarded.
