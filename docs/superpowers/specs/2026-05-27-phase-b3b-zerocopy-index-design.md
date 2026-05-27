# Phase B3b — Zero-copy persisted FM-index (design)

**Status:** Spec for review — 2026-05-27. The keystone of Phase B's persistence pillar, under the
contract thesis in [`docs/OPEN_PROBLEMS.md`](../../OPEN_PROBLEMS.md). Follows B3a (compact
`SampledSuffixArray`, merged).

## 1. The capability we are shooting for

> **Build the reference index once into a portable artifact; thereafter every analysis against it
> starts in milliseconds (it memory-maps the index — it never rebuilds), runs in bounded memory, and
> is byte-identically reproducible.**

B3b is the **keystone**: the persisted, **zero-copy, memory-mapped FM-index** that makes
*instant-start* and *bounded-memory query* real. It ladders directly into **B3c** (`rosalind index`
+ load) and **B4** (wire `align`/`variants` onto it), which turn it into the demonstrable
*build-once → instant, bounded, reproducible analysis* workflow a builder can use. B3b itself
delivers this at the **library** level (`IndexWriter`/`IndexReader` → `FmIndexView` → query).

## 2. Scope

**In:**
- A **`BwtBacking` trait** — the seam that lets one FM query algorithm run over either owned (in-RAM)
  or borrowed (mmap) data.
- A **deterministic on-disk format** (extends the existing versioned `genomics/index/format.rs`) that
  serialises a built `GenomeIndex` (FM-index + contig table + 2-bit reference + compact sampled-SA).
- A **serializer** (`GenomeIndex` → bytes) and a **`FmIndexView<'a>`** (borrowed backing over the
  mmap) that answers `backward_search`/`sa_at`/`locate_exact` with **no allocation and no rebuild**.
- The **equivalence gate**: the view returns byte-identical results to the in-RAM index it was built
  from, over a large pattern battery.
- **`docs/index-format.md`** — the format documented as a versioned, forkable contract.

**Out (deferred, by design):**
- `rosalind index` subcommand + `IndexReader` CLI wiring → **B3c**.
- Wiring `align`/`variants`/`somatic`/pileup onto the persisted index → **B4**.
- `MemoryBudget` *enforcement* + `rosalind plan` → **Phase C** (B3b lays the mmap substrate the budget
  story builds on; it does not add the budget API).
- **Leaner/faster rank** (drop the redundant per-block 2-bit BWT, or rank-over-2-bit) → a future
  **format v2** (query-speed optimization; see §4.3). Not a prerequisite for the capability.

## 3. The `BwtBacking` seam (one algorithm, two backings)

Today `BlockedFMIndex` carries the FM query algorithm as inherent methods over its owned `Vec`-backed
fields (`blocks`, `boundaries`, `c_table`, `sampled`). B3b extracts the *data surface* that algorithm
needs into a trait, and makes the algorithm generic over it:

```rust
pub trait BwtBacking {
    fn bwt_len(&self) -> usize;
    fn block_size(&self) -> usize;
    fn sentinel_pos(&self) -> usize;
    fn c_table(&self) -> [u32; 6];
    /// Cumulative count of `base` in BWT[.. block_idx * block_size] (boundary before the block).
    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32;
    fn boundary_sentinel(&self, block_idx: usize) -> u32;
    /// Rank of `symbol` within block `block_idx` over its first `within` positions.
    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32;
    /// The BWT symbol at `within` inside block `block_idx` (caller handles the sentinel position).
    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol;
    /// Sampled suffix-array value at BWT position `index`, if sampled.
    fn sampled_at(&self, index: usize) -> Option<u32>;
}
```

The FM operations — `rank`, `symbol_at`, `lf_index`, `backward_search`, `sa_at`, `locate_interval`,
`total` — become a single implementation generic over `B: BwtBacking` (free functions, or a blanket
`impl<B: BwtBacking>` on a thin wrapper). **Owned backing** = `BlockedFMIndex` implements `BwtBacking`
by delegating to its existing fields; its public methods become thin wrappers over the generic ops,
so `GenomeIndex`, the aligner, and all existing tests are unaffected. **Borrowed backing** =
`FmIndexView<'a>` implements `BwtBacking` by reading mmap slices.

This is the spec's "one algorithm, two backings": the *algorithm* lives once; the *data access* is the
only thing that differs, which is exactly what makes the borrowed view trustworthy (it can't diverge
on logic — only on backing).

## 4. The on-disk format

### 4.1 Mechanism (zero-copy)
Every multi-byte array is stored as **little-endian bytes at an 8-byte-aligned file offset**. The view
obtains a `&[u64]`/`&[u32]` over the mmap via a checked `slice::align_to` — the format guarantees
alignment (so the unaligned prefix is empty) and the header validates **little-endian-only** (so the
reinterpretation is value-correct on supported hosts). A small `unsafe` with a documented invariant,
in the spirit of the existing `util::mmap` (raw `libc::mmap`). **No new dependency** (no `bytemuck`).
Result: rank/`sa_at` run over borrowed slices — no allocation, no rebuild, OS-paged residency.

### 4.2 Sections (extending `genomics/index/format.rs`)
The existing fixed header (magic `ROSALIND`, version, endian marker, `contig_count`, section-table
location, reference BLAKE3) is retained; the `SectionKind` set is extended. All section payloads are
flat, 8-aligned, little-endian:

| Section | Payload |
|---|---|
| `Contigs` | per contig: `name_len:u32, name:[u8], length:u64, global_offset:u64` |
| `Reference2bit` | the 2-bit forward reference (`CompressedDNA`): `len:u64`, packed `data:[u64]`, ambiguity `bits:[u64]` — self-contained, for `ref_base` lookups (B4) without a separate `.fa` |
| `FmMeta` | `block_size:u64, bwt_len:u64, sentinel_pos:u64, sa_sample_rate:u64, num_blocks:u64, c_table:[u32;6]` |
| `Boundaries` | `num_blocks + 1` entries of `{ cumulative_counts:[u32;5], sentinel_count:u32 }` |
| `Blocks` | a `num_blocks × u64` offset directory, then per block: `{ start:u64, end:u64, sentinel_offset:i64 (-1 = none), stride:u64 }` + the 2-bit BWT (`data:[u64]`, ambiguity `bits:[u64]`) + the occ's five rank bitvectors `[[u64];5]` + five superblock arrays `[[u32];5]` + `totals:[u32;5]` |
| `SaSamples` | the compact `SampledSuffixArray` (B3a): `rate:u64, bwt_len:u64, marks:[u64], superblocks:[u32], values:[u32]` |

The per-block directory makes the `Blocks` section self-describing (the view seeks block `k` by its
offset), robust to the last block's odd size.

### 4.3 Faithful serialization (the size characteristic, documented)
We serialise the blocked structure **faithfully**: both the per-block 2-bit BWT (used by
`block_symbol`) and the five rank bitvectors (used by `block_rank`) — ≈7 bits/symbol. This is the
honest cost of "one algorithm, two backings" in v1. **RSS stays bounded regardless** (mmap'd /
OS-paged — file size costs disk, not resident memory), so the capability holds. `docs/index-format.md`
documents the size and names the **v2 lean path** (derive `block_symbol` from the rank bitvectors and
drop the 2-bit BWT, or rank-over-2-bit) as a clean `IndexVersion` bump — a query-speed/size
optimization, not a correctness change.

### 4.4 Determinism + integrity
SA-IS is deterministic; sections are written in a fixed order with fixed-width LE encoding; no
timestamps appear in the file ⇒ **byte-identical across repeated builds**. The header's
`reference_blake3` lets a consumer verify an index came from exactly its reference (the stale-index
guard, used in B4).

## 5. `FmIndexView<'a>`
`IndexReader::open(path)` mmaps the file, validates the header/sections, parses the `ContigSet`, and
yields a `ReferenceIndex` owning the `MmapReadOnly`. `ReferenceIndex::view(&self) -> FmIndexView<'_>`
borrows mmap slices for the FM data and implements `BwtBacking`. Because the FM ops are generic over
`BwtBacking`, the view gets `backward_search`/`sa_at`/`locate_exact` for free — querying the mmap in
place. A `GenomeIndexView` (or the existing `GenomeIndex` made generic) pairs the `FmIndexView` with
the parsed `ContigSet` to expose `locate_exact -> Vec<Locus>` (boundary-aware, sorted), the same
surface B2 shipped.

## 6. Decomposition (two green sub-stages)

- **B3b.1 — `BwtBacking` trait + route `BlockedFMIndex` through it (owned backing).** Pure refactor:
  define the trait, move the FM ops to a generic implementation, implement `BwtBacking` for the owned
  index, keep its public methods as thin wrappers. Existing FM-index, aligner, and `genome_index`
  tests prove behavior is unchanged. The trait surface is designed (per §3) to serve both backings.
- **B3b.2 — format + serializer + `FmIndexView` + equivalence gate + `docs/index-format.md`.** Extend
  `format.rs` sections; write the deterministic serializer (`GenomeIndex` → bytes); implement
  `IndexReader::open` → `ReferenceIndex` → `FmIndexView` (the borrowed backing); add the
  `view == in-RAM` equivalence + byte-identical-build + no-rebuild gates; publish the format doc.

## 7. Success criteria (gates)

- **Equivalence:** build a `GenomeIndex` in RAM, serialise it, mmap it as an `FmIndexView`, and assert
  `backward_search`/`sa_at`/`locate_exact` are **byte-identical** to the in-RAM index over a large,
  fixed pattern battery (incl. multi-contig + boundary-straddle cases).
- **Determinism:** the index file is byte-identical across two builds of the same reference.
- **No rebuild on load:** the read path never calls SA-IS (`sais_u32`) — enforced structurally (the
  loader constructs a view, never `BlockedFMIndex::build*`) and asserted behaviorally.
- **Bounded residency:** opening + querying the view does not make the index resident (it is mmap'd;
  a basic residency/working-set check — the enforced RSS gate is Phase C).
- **Behavior preserved (B3b.1):** the owned-backing refactor changes no observable behavior (existing
  suites green).

## 8. Testing
- **B3b.1:** the existing FM-index/aligner/`genome_index` suites are the witnesses (owned-via-trait ==
  owned-direct); add a focused test that the generic ops over the owned backing match expected ranks
  / `sa_at` for a known reference.
- **B3b.2:** the equivalence battery (view vs in-RAM, many patterns over a multi-contig reference incl.
  an `N`-bearing fixture); determinism (two builds byte-equal); a corrupt/truncated index is rejected
  (section-extent + magic/version/endian validation); round-trip over a self-contained index (query
  with only the index, no separate reference).

## 9. Risks & mitigations
- **The `BwtBacking` refactor is the subtle part.** *Mitigation:* B3b.1 lands it behind the existing
  green suites (owned backing only, behavior unchanged) before any format work; the generic ops are
  one source of truth, so the borrowed view can only diverge on backing, which the equivalence gate
  catches loudly.
- **`align_to` UB if a section is misaligned.** *Mitigation:* the serializer pads every multi-byte
  section to 8 bytes; `IndexReader::open` validates alignment (and that section extents lie in the
  file); the view asserts the `align_to` prefix is empty. Endianness is header-validated LE-only.
- **Format churn before forkers depend on it.** *Mitigation:* `IndexVersion` is bumped on any layout
  change; the format is documented as v1, stability-deferred-to-1.0 (per OPEN_PROBLEMS), so the v2 lean
  path (§4.3) is a clean migration, not a break.

## 10. Decisions (this stage)
- **Faithful serialization in v1** (§4.3) — the capability is mmap-delivered (bounded RSS) regardless;
  lean-rank is a versioned v2 query-speed win.
- **Dependency-free zero-copy** via aligned LE sections + checked `align_to` (§4.1) — dep-light *is* the
  portability/reliability the capability promises.
- **One algorithm, two backings** via `BwtBacking` (§3) — the view's trustworthiness comes from sharing
  the algorithm, not reimplementing it.
