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
| `Blocks` | `[directory: u64; num_blocks]` (record offsets relative to the section start), then per block: `start:u64, end:u64, sentinel_offset:i64 (-1 = none), stride:u64, bwt_data_words:u64, bwt_amb_words:u64, occ_bitvec_words:u64, occ_superblock_len:u64`, then `bwt_data:[u64]`, `bwt_amb:[u64]`, `occ_bitvecs:[[u64]; 5]`, `occ_superblocks:[[u32]; 5]`, `totals:[u32; 5]` (record padded to 8). The directory makes block `k` O(1) to locate; the reader validates every record's extent at `open`. |
| `SaSamples` | `rate:u64, bwt_len:u64, marks_words:u64, superblocks_len:u64, values_len:u64`, then `marks:[u64]`, `superblocks:[u32]`, `values:[u32]` (the compact `SampledSuffixArray`; superblock stride is the fixed `RANK_STRIDE`; `rate` is checked to equal the header `sa_sample_rate`). |

## Determinism

SA-IS is deterministic; sections are written in fixed order with fixed-width LE
encoding and zero padding; no timestamps appear. ⇒ the file is **byte-identical
across repeated builds of the same genome** (gated by `build_is_deterministic`).

## Integrity

`IndexReader::open` validates the magic/version/endian, every section's
8-alignment and in-file extent, the contig global-offset consistency, every block
record's extent (via checked arithmetic — a corrupt record is rejected at open,
never silently mis-read), and `SaSamples.rate == sa_sample_rate`. A corrupt or
truncated file is rejected rather than mis-queried.

## Faithful serialization (size) and the v2 lean path

v1 serialises the blocked structure **faithfully**: both the per-block 2-bit BWT
(for `block_symbol`) and the five occ rank bitvectors (for `block_rank`), ≈7
bits/symbol. This is the honest cost of "one algorithm, two backings" — and RSS
stays bounded regardless, because the index is mmap'd / OS-paged (file size costs
disk, not resident memory). A future **format v2** may derive `block_symbol` from
the rank bitvectors (or rank over the 2-bit BWT) and drop the redundant copy — a
query-speed/size optimization behind a clean `IndexVersion` bump, not a
correctness change.
