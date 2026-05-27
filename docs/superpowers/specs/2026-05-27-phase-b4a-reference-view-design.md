# Phase B4a — self-contained reference access (`ReferenceView`) (design)

**Status:** Spec for review — 2026-05-27. The first sub-stage of Phase B4 (wire the consumers onto the persisted multi-contig index), under the contract thesis in [`docs/OPEN_PROBLEMS.md`](../../OPEN_PROBLEMS.md). Follows B3a/B3b/B3c (persisted zero-copy index + `rosalind index`/`locate`, all merged — PRs #15/#16/#17).

## 1. The capability we are shooting for

> **Read the reference sequence directly from a persisted index — no separate `--reference` FASTA, no full-reference allocation — decoding only the small windows a consumer asks for.**

B3b stored the 2-bit forward reference into the index (the `Reference2bit` section) precisely so consumers wouldn't need a separate FASTA — but gave it **no read path**. B4a adds that read path: a borrowed, zero-copy **`ReferenceView<'a>`** that decodes bases on demand from the memory-mapped section. This is the shared foundation the rest of B4 needs — **B4b** (the aligner's banded-DP refinement window) and **B4c** (variants' `ref_base` / pileup) both read reference bases, and after B4a they read them from the `.idx` alone.

## 2. Scope

**In:**
- A borrowed **`ReferenceView<'a>`** over the `Reference2bit` section: `len()`, `base_at(global) -> u8` (on-demand ASCII decode), `decode_window(start, end, &mut Vec<u8>)` (bounded buffer fill).
- **`ReferenceIndex::reference_view() -> Result<ReferenceView<'_>, IndexIoError>`** — parses + validates the section (extent, 8-alignment, little-endian host) and yields the borrowed view.
- The **equivalence gate**: decoded bases (`base_at` over the whole reference, `decode_window` over ranges) are byte-identical to the original `index.reference()`, including an `N`-bearing multi-contig reference.

**Out (deferred, by design):**
- Wiring `ReferenceView` into the aligner's DP window → **B4b**.
- Wiring `ReferenceView` into `variants`/pileup `ref_base` → **B4c**.
- Any contig-coordinate API (`ReferenceView` is **global-coordinate**, matching the concatenated 2-bit reference; `(contig, pos)` resolution stays with `ContigSet`, used by the consumers in B4b/B4c).

## 3. The `Reference2bit` section (already written by B3b)

The B3b serializer (`genomics/index/io.rs`) writes, into the 8-aligned `Reference2bit` section:

```
len:u64, data_words:u64, amb_words:u64, data:[u64; data_words], amb:[u64; amb_words]
```

`data` is the 2-bit packing of `CompressedDNA` (32 bases per `u64`, A/C/G/T = 0/1/2/3); `amb` is the ambiguity bitmap (1 bit per base; a set bit marks `N`). This is exactly `CompressedDNA::compress(index.reference())`'s `words()` + `ambiguity().bits()`.

## 4. `ReferenceView<'a>`

```rust
pub struct ReferenceView<'a> {
    len: usize,
    data: &'a [u64], // 2-bit packed, 32 bases/word
    amb: &'a [u64],  // ambiguity bits, 1/base (set = N)
}
```

- **Construction** (`pub(crate) fn new(bytes, sections)`, called by `ReferenceIndex::reference_view`): reject big-endian hosts; locate the `Reference2bit` section (reuse `section_bytes`); read `len/data_words/amb_words` (reuse `read_u64`); slice `data`/`amb` via the checked `as_u64_slice` (empty `align_to` prefix/suffix — the same zero-copy discipline as `FmIndexView`); **validate the slices cover `len` bases** (`data.len() * 32 >= len` and `amb.len() * 64 >= len`) so `base_at` cannot index out of bounds. Errors (not panics) on any mismatch.
- **`base_at(global) -> u8`** decodes one base, mirroring `CompressedDNA::base_at`: if `amb[global / 64] >> (global % 64) & 1 == 1` → `b'N'`; else `((data[global / 32] >> ((global % 32) * 2)) & 0b11)` → `A/C/G/T`. `debug_assert!(global < len)`.
- **`decode_window(start, end, out: &mut Vec<u8>)`** clears `out` and pushes `base_at(i)` for `i in start..end.min(len)` (`end` is clamped to `len`, so an over-range request never panics; `start <= end` is the caller's contract). Bounded by the window size (the caller decodes only what it needs — e.g. the aligner's band + read length — never the whole genome).
- The 2-bit decode is a 4-line replication of `CompressedDNA::base_at`, annotated as such; the equivalence gate (decoded == original) is the guarantee against drift. (Extracting a shared `pub(crate)` decode helper across `CompressedDNA`/`FmIndexView`/`ReferenceView` is deferred — out of B4a's scope; the gate makes replication safe.)

## 5. Why on-demand zero-copy (not reconstruction)

Reconstructing the full reference into a `Vec<u8>` at open would be simpler but **O(reference) RAM** — defeating the persisted index's bounded-query guarantee (the whole point of B3b). `ReferenceView` instead holds only two borrowed slices + a length and decodes on demand, so reference access stays bounded (a consumer pays only for the windows it reads). This matches the contract thesis and the `FmIndexView` precedent.

## 6. Success criteria (gates)

- **Equivalence:** for a built `GenomeIndex` (multi-contig, `N`-bearing), `reference_view().base_at(i)` equals `index.reference()[i]` for every `i`, and `decode_window(s, e, …)` equals `&index.reference()[s..e]` over several ranges (including ranges spanning contig boundaries in the concatenated coordinate and an all-`N` stretch).
- **No rebuild / self-contained:** `reference_view()` reads only the mmap (never `sais_u32`/`BlockedFMIndex::build`); decoding works after the source FASTA is gone (the bases come from the `.idx`).
- **Bounded:** `ReferenceView` is a small borrow (slices + scalar), not an owned copy; `size_of::<ReferenceView>()` is independent of genome size.
- **Integrity:** a `Reference2bit` section too small for `len` bases is rejected at `reference_view()` (Err, not a later panic).

## 7. Testing

- A unit test in `genomics/index/view.rs`: build → serialize → open → `reference_view()`; assert `base_at` over the whole reference == the original, `decode_window` over ranges == the slices, including an `N`-bearing multi-contig fixture (reuse the B3b/B3c fixture style). Assert `size_of::<ReferenceView>()` is small.
- Round-trip robustness is already covered by B3b's determinism/round-trip; B4a adds the *decode* direction.

## 8. Risks & mitigations

- **Decode divergence from `CompressedDNA`.** *Mitigation:* the decode is a verbatim mirror of `CompressedDNA::base_at`; the equivalence gate (over an `N`-bearing reference) catches any drift loudly.
- **Out-of-bounds on a corrupt/short section.** *Mitigation:* `new` validates `data`/`amb` cover `len` bases and returns `Err`; `base_at` is only reachable after that validation. `decode_window` clamps `end` to `len`, so an over-range request never panics.
- **`amb_words` over-allocation quirk.** B3b stores `amb_words = ceil(len/32)` (the `CompressedDNA` allocation), larger than the `ceil(len/64)` strictly needed; `new`'s `amb.len() * 64 >= len` check accommodates this, and `base_at` indexes `amb[i/64]` which is in-bounds for `i < len`.

## 9. Decisions (this stage)

- **On-demand zero-copy decode** (§5) — bounded reference access is the whole point; full reconstruction is rejected.
- **Global-coordinate view** (§2) — `ReferenceView` matches the concatenated 2-bit reference; `(contig, pos)` mapping stays in `ContigSet` (the consumers compose them in B4b/B4c).
- **Replicate the 2-bit decode, gate against drift** (§4) — a shared decode helper across `CompressedDNA`/`FmIndexView`/`ReferenceView` is deferred (not B4a's job); the equivalence gate makes the 4-line replication safe.
- **Lives in `genomics/index/view.rs`** alongside `FmIndexView`/`GenomeIndexView`, reusing the existing `section_bytes`/`as_u64_slice`/`read_u64` helpers; re-exported via `genomics::index` → `genomics`.
