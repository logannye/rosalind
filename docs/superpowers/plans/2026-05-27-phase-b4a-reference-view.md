# Phase B4a — self-contained reference access (`ReferenceView`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Read the reference sequence directly from a persisted index — a borrowed, zero-copy `ReferenceView` over the `Reference2bit` section that decodes bases on demand (no full-reference allocation) — so B4b (aligner DP window) and B4c (variants `ref_base`) can be served from the `.idx` alone.

**Architecture:** A `ReferenceView<'a>` in `genomics/index/view.rs` holding borrowed `&'a [u64]` slices over the mmap'd `Reference2bit` data + ambiguity bits, obtained via the existing checked `as_u64_slice`/`section_bytes`/`slice_exact`/`read_u64` helpers (same zero-copy discipline as `FmIndexView`). `base_at(global)` decodes one ASCII base on demand, mirroring `CompressedDNA::base_at`; `decode_window` fills a bounded caller buffer. Constructed via `ReferenceIndex::reference_view()`.

**Tech Stack:** Rust 2021 (MSRV 1.72 — no `div_ceil`); checked `slice::align_to` (little-endian host); no new dependencies.

This is sub-stage **B4a** of Phase B4 (wire consumers onto the persisted multi-contig index), from the design `docs/superpowers/specs/2026-05-27-phase-b4a-reference-view-design.md`. It follows B3a/B3b/B3c (merged, PRs #15/#16/#17). **Out of scope (deferred):** wiring `ReferenceView` into the aligner (**B4b**) and into `variants`/pileup (**B4c**). B4a delivers and gates only the reader.

---

## File structure

- `src/genomics/index/view.rs` — **Modify.** Add `ReferenceView<'a>` (struct + `pub(crate) fn new` + `len`/`is_empty`/`base_at`/`decode_window`), reusing the file's existing `section_bytes`/`as_u64_slice`/`slice_exact`/`read_u64` helpers.
- `src/genomics/index/io.rs` — **Modify.** Add `ReferenceIndex::reference_view()` (mirrors the existing `view()`/`genome_view()` accessors); add the equivalence-gate test to the existing `#[cfg(test)] mod tests` (reusing `temp_path`/`sample_index`).
- `src/genomics/index/mod.rs` — **Modify.** Re-export `ReferenceView` from the `view` line.
- `src/genomics/mod.rs` — **Modify.** Add `ReferenceView` to the `pub use index::{…}` line.

---

## Task 1: `ReferenceView` + `reference_view()` accessor + equivalence gate

**Files:**
- Modify: `src/genomics/index/view.rs`, `src/genomics/index/io.rs`, `src/genomics/index/mod.rs`, `src/genomics/mod.rs`

- [ ] **Step 1: Write the failing equivalence-gate test.** In `src/genomics/index/io.rs`, append to the existing `#[cfg(test)] mod tests` (it has `temp_path` + `sample_index`, and `sample_index` is multi-contig + `N`-bearing):

```rust
    #[test]
    fn reference_view_decodes_bytes_identical_to_in_ram() {
        let idx = sample_index();
        let path = temp_path("refview");
        IndexWriter::create(&path).unwrap().write_genome_index(&idx).unwrap();
        let loaded = IndexReader::open(&path).unwrap();
        let rv = loaded.reference_view().unwrap();
        let reference = idx.reference();

        // base_at over the whole reference equals the original (incl. N positions).
        assert_eq!(rv.len(), reference.len());
        assert!(!rv.is_empty());
        for i in 0..rv.len() {
            assert_eq!(rv.base_at(i), reference[i], "base_at mismatch @ {i}");
        }

        // decode_window over several ranges == the corresponding slices, including
        // a boundary-spanning range, an N-bearing range, the full reference, and an
        // over-range request (end clamped to len).
        let mut buf = Vec::new();
        let n = reference.len();
        for (s, e) in [(0usize, 10usize), (6, 14), (12, 22), (0, n), (n - 3, n + 5)] {
            rv.decode_window(s, e, &mut buf);
            assert_eq!(buf.as_slice(), &reference[s..e.min(n)], "decode_window {s}..{e}");
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn reference_view_is_a_small_borrow() {
        // Bounded: the view is slices + a scalar, not an owned copy of the reference.
        assert!(
            std::mem::size_of::<crate::genomics::ReferenceView<'_>>() <= 64,
            "ReferenceView must be a small borrow, got {}",
            std::mem::size_of::<crate::genomics::ReferenceView<'_>>()
        );
    }
```

- [ ] **Step 2: Run them to verify they fail.**

Run: `cargo test --lib index::io::tests::reference_view 2>&1 | tail -20`
Expected: compile error — `ReferenceIndex::reference_view` and `crate::genomics::ReferenceView` do not exist.

- [ ] **Step 3: Implement `ReferenceView` in `src/genomics/index/view.rs`.** Add this (place it after the `GenomeIndexView` impl block, before the private free-fn helpers `section_bytes`/`as_u64_slice`/…):

```rust
/// A borrowed, zero-copy view over the persisted 2-bit forward reference
/// (`Reference2bit`). Decodes bases on demand from the mmap — no full-reference
/// allocation — so consumers (the aligner DP window in B4b, variants' `ref_base`
/// in B4c) read the reference from the `.idx` alone. Uses **global** (concatenated)
/// coordinates; `(contig, pos)` mapping stays with `ContigSet`. The little-endian
/// host requirement is the same as `FmIndexView` (checked in `new`).
#[derive(Debug)]
pub struct ReferenceView<'a> {
    len: usize,
    /// 2-bit packed bases, 32 per `u64` word (A/C/G/T = 0/1/2/3).
    data: &'a [u64],
    /// Ambiguity bits, one per base (a set bit marks `N`).
    amb: &'a [u64],
}

impl<'a> ReferenceView<'a> {
    /// Parse + validate the `Reference2bit` section into a borrowed view.
    pub(crate) fn new(bytes: &'a [u8], sections: &[SectionEntry]) -> Result<Self, IndexIoError> {
        if cfg!(target_endian = "big") {
            return Err(IndexIoError::Invalid(
                "zero-copy reference view requires a little-endian host".to_string(),
            ));
        }
        let section = section_bytes(bytes, sections, SectionKind::Reference2bit)?;
        let mut o = 0usize;
        let len = read_u64(section, &mut o)? as usize;
        let data_words = read_u64(section, &mut o)? as usize;
        let amb_words = read_u64(section, &mut o)? as usize;
        let data = as_u64_slice(slice_exact(section, &mut o, data_words.saturating_mul(8))?)?;
        let amb = as_u64_slice(slice_exact(section, &mut o, amb_words.saturating_mul(8))?)?;
        // The slices must cover `len` bases so `base_at` cannot index out of bounds.
        if data.len().saturating_mul(32) < len || amb.len().saturating_mul(64) < len {
            return Err(IndexIoError::Invalid(
                "Reference2bit section too small for the declared length".to_string(),
            ));
        }
        Ok(Self { len, data, amb })
    }

    /// Number of reference bases.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the reference is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The ASCII base at global position `global`. Mirrors `CompressedDNA::base_at`:
    /// a set ambiguity bit decodes to `N`, otherwise the 2-bit code maps to A/C/G/T.
    pub fn base_at(&self, global: usize) -> u8 {
        debug_assert!(global < self.len, "reference index out of range");
        if self.amb[global / 64] & (1u64 << (global % 64)) != 0 {
            return b'N';
        }
        let code = ((self.data[global / 32] >> ((global % 32) * 2)) & 0b11) as u8;
        match code {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }

    /// Decode `[start, end.min(len))` into `out` (cleared first). `end` is clamped
    /// to `len`, so an over-range request never panics; `start <= end` is the
    /// caller's contract. Bounded by the window size — never the whole genome.
    pub fn decode_window(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.clear();
        let end = end.min(self.len);
        for i in start..end {
            out.push(self.base_at(i));
        }
    }
}
```

(`section_bytes`, `as_u64_slice`, `slice_exact`, `read_u64` already exist as private free functions in `view.rs`. `SectionEntry`/`SectionKind` and `IndexIoError` are already imported there. The 2-bit decode is a deliberate 4-line mirror of `CompressedDNA::base_at` in `genomics/compressed_dna.rs`; the equivalence gate is the guarantee against drift.)

- [ ] **Step 4: Add the `reference_view()` accessor in `src/genomics/index/io.rs`.** In `impl ReferenceIndex`, after `genome_view()`, add:

```rust
    /// Borrow a zero-copy view over the persisted 2-bit reference (`Reference2bit`).
    pub fn reference_view(
        &self,
    ) -> Result<crate::genomics::index::view::ReferenceView<'_>, IndexIoError> {
        crate::genomics::index::view::ReferenceView::new(self.mmap.as_bytes(), &self.sections)
    }
```

- [ ] **Step 5: Re-export `ReferenceView`.** In `src/genomics/index/mod.rs`, change:

```rust
pub use view::{FmIndexView, GenomeIndexView};
```

to:

```rust
pub use view::{FmIndexView, GenomeIndexView, ReferenceView};
```

In `src/genomics/mod.rs`, add `ReferenceView` to the `pub use index::{…}` list (keep it alphabetical / consistent with the existing line), e.g.:

```rust
pub use index::{
    estimate_build_working_set, render_plan_line, FmIndexView, GenomeIndexView, IndexBuildReport,
    IndexHeader, IndexReader, IndexWriter, ReferenceIndex, ReferenceView,
};
```

(Match the exact existing item set in that `use` — only ADD `ReferenceView`; do not drop anything.)

- [ ] **Step 6: Run the gate + build + fmt.**

Run: `cargo test --lib index::io::tests::reference_view 2>&1 | tail -20` (both `reference_view_decodes_bytes_identical_to_in_ram` and `reference_view_is_a_small_borrow` PASS), then
`cargo test --lib 2>&1 | grep -E "test result:"` (no failures), then
`cargo build --lib 2>&1 | grep -iE 'error|warning'` (none — `ReferenceView`/`base_at`/`decode_window`/`reference_view` are all used by the tests, so no dead-code warning), then
`cargo fmt --all` then `cargo fmt --all -- --check` (clean).

- [ ] **Step 7: Commit (stage the 4 files).**

```bash
git add src/genomics/index/view.rs src/genomics/index/io.rs src/genomics/index/mod.rs src/genomics/mod.rs
git commit -m "feat(genomics/index/view): ReferenceView — zero-copy reference access from the index" \
  -m "A borrowed ReferenceView over the persisted Reference2bit section (base_at + decode_window, on-demand 2-bit decode mirroring CompressedDNA::base_at, no full-reference Vec), via ReferenceIndex::reference_view(). Gated by base_at/decode_window == the original reference over an N-bearing multi-contig index. The shared reference-access foundation B4b (aligner DP window) and B4c (variants ref_base) consume." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

Then verify: `git show --stat HEAD` lists ONLY the 4 files; `git status --short` is empty.

---

## Final verification (before the B4a PR)

- `cargo test` — full suite green (the new `reference_view_*` tests + the unchanged B3b/B3c suites).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo fmt --all -- --check` — clean.
- No-rebuild (structural): `reference_view()` reads only the mmap; `grep -n "reference_view" src/genomics/index/io.rs` shows it calling `ReferenceView::new`, never `BlockedFMIndex::build`/`sais_u32`.

## Self-Review

- **Spec coverage (`2026-05-27-phase-b4a-reference-view-design.md`):**
  - §4 `ReferenceView` (`new`/`len`/`base_at`/`decode_window`) ✔ Task 1 Step 3.
  - §4 `ReferenceIndex::reference_view()` ✔ Step 4.
  - §5 on-demand zero-copy (borrowed slices, no `Vec`; `decode_window` bounded) ✔ — the struct holds `&'a [u64]` + a scalar; gated by `reference_view_is_a_small_borrow`.
  - §6 gates — equivalence (`base_at`/`decode_window` == original over N-bearing multi-contig) ✔; bounded (`size_of` ≤ 64) ✔; integrity (`new` rejects a too-small section) ✔ (the `data.len()*32 < len` / `amb.len()*64 < len` check); no rebuild ✔ (final verification grep + reads-only-mmap).
  - §7 testing (build→serialize→open→`reference_view`, ranges incl. boundary-spanning + N + over-range clamp) ✔.
- **Type/name consistency:** `ReferenceView<'a> { len: usize, data: &'a [u64], amb: &'a [u64] }`; `new(bytes, sections) -> Result<Self, IndexIoError>`; `len()`, `is_empty()`, `base_at(usize) -> u8`, `decode_window(usize, usize, &mut Vec<u8>)`; `ReferenceIndex::reference_view() -> Result<ReferenceView<'_>, IndexIoError>`; re-exported as `crate::genomics::ReferenceView`. The decode matches `CompressedDNA::base_at` (ambiguity word `i/64` bit `i%64`; data word `i/32` shift `(i%32)*2`).
- **No placeholders:** every step ships complete code or an exact command + expected output.
- **MSRV 1.72:** no `div_ceil`; `saturating_mul` for the length/extent math; `1u64 << (i % 64)` / `(i % 32) * 2` are plain shifts.
- **Reuse/DRY:** reuses `view.rs`'s `section_bytes`/`as_u64_slice`/`slice_exact`/`read_u64`; mirrors the `FmIndexView::new` construction + LE-host pattern; the 2-bit decode mirrors `CompressedDNA::base_at` (drift caught by the gate).
- **Scope:** reader only — no aligner (B4b) or variants/pileup (B4c) wiring; `ReferenceView` is global-coordinate (no contig API).
- **Shared-tree hazard:** the implementer/reviewer dispatches must forbid `cargo fix`/mutating git, stage only the named files, and self-verify `git show --stat HEAD`; the coordinator verifies the commit stat + clean tree at the task boundary.
```
