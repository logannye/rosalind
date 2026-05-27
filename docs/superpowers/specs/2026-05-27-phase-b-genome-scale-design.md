# Phase B — Make it real on a genome (multi-contig + persisted index + streaming ingestion)

**Status:** Spec for review — 2026-05-27. Second vertical of the target architecture
(`2026-05-26-rosalind-target-architecture.md`, §9.2). Builds directly on the Phase A kernel
(`core/`, `pileup/`, `call/`, `io/vcf`, `provenance/`).

**One-line goal.** Take the Phase A calling vertical from a single-contig, rebuild-the-index-every-run,
uncompressed-input toy to a tool that runs on a **real, multi-contig genome from a build-once,
memory-mapped index**, reads **gzipped** FASTA/FASTQ as a streaming library, and **composes over
pipes** — without surrendering any of the five promises (deterministic, memory-as-contract,
standards-compliant, verifiable, embeddable).

**Why this phase, and why this shape.** A builder forking Rosalind today is hard-blocked: `align`
indexes only the first FASTA record, the FM-index is rebuilt in RAM on every run, and only
uncompressed input is accepted. The Phase A `PileupColumn` substrate is clean but starved — it can
only ever see one contig of one toy. Phase B raises the ceiling on *all* downstream builder value by
making the engine work on a real genome. We deliberately do the roadmap-literal infrastructure (not a
builder-experience detour) because the substrate is worthless on data it cannot ingest.

---

## 1. Success criteria (definition of done)

Phase B is done when, on a multi-contig reference with gzipped reads, all of the following hold:

- **Multi-contig:** `align` indexes and maps reads across *every* contig; emitted `core::AlignedRead`s
  carry the correct `(contig, pos)`; the "first FASTA record only" warning path is gone. Reads whose
  alignment would straddle a contig boundary are rejected, not mis-mapped.
- **Build-once, mmap-load:** `rosalind index ref.fa[.gz]` produces a single index artifact; `align` /
  `variants` / `somatic` load it via mmap and **never rebuild the FM-index** (no SA-IS on the read
  path). Querying the loaded index returns **byte-identical** results to the in-RAM index it was built
  from.
- **Streaming ingestion:** FASTA and FASTQ readers live in `io/`, stream multi-record input, and read
  plain or gzip/bgzf transparently (auto-detected); a `.gz` reference/reads produce identical results
  to their decompressed equivalents.
- **Self-contained index:** the index stores the 2-bit reference, so `variants`/`somatic` need only
  `--index` (not also `--reference`). A stale index (reference BLAKE3 mismatch) is rejected with a
  clear error.
- **Deterministic:** the index file is byte-identical across repeated builds (no timestamps, fixed
  section order, fixed-width LE encoding); VCF + receipt remain byte-identical across runs and across
  the file-based vs piped paths.
- **Memory honest:** loading the index does not make it resident (near-zero declared RSS; OS-paged).
  The `MemoryBudget` contract continues to govern the streaming working set; the index's
  page-cache residency is documented as OS-managed and separate.
- **Pipe-native:** `align --reads - | rosalind sort - | rosalind variants --alignments -` produces the
  same VCF as the file-based path; the receipt records the index artifact's hash.
- **Embeddable:** the readers, the index reader/writer, and `FmIndexView` are library-first — no
  `rust_htslib` or `anyhow` types in public signatures; the on-disk format is documented as a
  public, forkable artifact (`docs/index-format.md`).
- **Standards:** multi-contig VCF carries a `##contig` line per contig and a multi-contig BAM carries an
  `@SQ` line per contig, both from the `ContigSet`; `bcftools`/`samtools` validate them in CI when
  available.

---

## 2. Scope

**In:**
1. `io/fasta` + `io/fastq` — streaming, multi-record, gzip/bgzf-transparent readers; ingestion leaves
   `main.rs`.
2. Multi-contig FM-index build over the concatenated genome + global↔local `Locus` mapping.
3. Index persistence — a zero-copy on-disk layout, a serializer, an mmap-backed `FmIndexView`, and a
   `rosalind index` subcommand. Self-contained (stores the 2-bit reference).
4. Consumer wiring — `align`/`variants`/`somatic` load the persisted multi-contig index; the pileup
   walks the whole genome contig-by-contig; `@SQ`/`##contig` from the `ContigSet`.
5. Pipe-native `-`/stdin/stdout composition; receipt extended to cover the index artifact.
6. Minimal per-pillar CI for the surface B introduces + a docs-truth pass (README scope +
   `docs/index-format.md`).

**Out (deferred, by design — do not let these creep in):**
- **Alignment *algorithm* rewrite** (real MAPQ calibration, soft-clip scoring, better chaining). Phase B
  keeps the existing seed→chain→extend math and only makes it multi-contig + index-backed. Algorithm
  quality is a later `align/` phase.
- **rayon / real RSS gate / `rosalind plan` / `rosalind verify`** → Phase D. (B exposes the honest
  memory story and a basic load-residency check; the *enforced* RSS gate is D.)
- **`.csi` / bgzf random-access region `fetch`** → Phase D. B reads bgzf by sequential full
  decompression only.
- **Germline indels** → Phase C. **Somatic indels** remain Phase C (already deferred in Phase A).
- **Python binding** → Phase P. (The readers/view/index API are shaped to be bindable, but no PyO3
  surface is added here.)
- **Theory-layer feature-gating, clippy `-D warnings` across the legacy tree, `lib.rs` crate-doc
  honesty pass** → Phase E. (Noted as a known builder-legibility wart; intentionally untouched here.)

---

## 3. Locked decisions (2026-05-27)

These were settled in the brainstorming dialogue and are not re-opened during implementation:

1. **Phase B mission = genome-scale infrastructure, roadmap-literal.** The four roadmap pillars, not a
   builder-experience reprioritization.
2. **Index memory model = mmap read-only, OS-paged.** The `MemoryBudget` contract governs the
   *streaming working set* (pileup/sort, coverage-bounded). The index is mmap'd: near-zero declared
   RSS, demand-paged, evictable, shareable across processes. Index *build* is O(reference) and
   documented separately. No `--load-index`/pin option now (YAGNI; revisit in D with the RSS gate).
3. **Index load = zero-copy, via a build-side/read-side split.** The existing `BlockedFMIndex` stays the
   *build-side* type (builds in RAM, then serializes). A new *read-side* `FmIndexView<'a>` borrows
   `&'a [u8]` slices from the mmap and runs the FM algorithm directly over them. **On-disk layout =
   query layout.** One source of truth for the FM algorithm via a `BwtBacking` trait implemented by
   both; an equivalence property test asserts the view returns identical results to the in-RAM index.
   No throwaway owned-deserialization path.
4. **`sa_samples` stays `u32`** → a documented ≤4.29 Gbp concatenated-genome ceiling (covers human and
   the overwhelming majority of edge targets; keeps SA-sample memory lean). `build_multi` errors
   clearly if the concatenated length reaches 2³². The `u64` coordinate API in `ContigSet` is
   preserved.
5. **Self-contained index** — the index stores the 2-bit forward reference, so `variants`/`somatic`
   require only `--index`. A stale index is rejected via the header's reference BLAKE3.
6. **SA samples are stored compactly (sparse), not as a dense full-length array.** The current
   `BlockedFMIndex` holds `sa_samples: Vec<u32>` of length `bwt_len` with `u32::MAX` holes
   (`fm_index.rs:469`) — so the sample *rate* only changes density, never storage size (~12 GB for a
   human genome; a persisted index would be ~13 GB). The on-disk `SaSamples` section and the read-side
   `sa_at` use a compact representation — a sampled-position **bitvector** (reusing `rank_select`) plus a
   **packed array of the sampled values** — making SA-sample storage sub-linear (`O(bwt_len / rate)`).
   The build-side is moved to the same compact form, eliminating the dense array. This is a prerequisite
   for honestly persisting a bounded-memory index, not an optimization.

---

## 4. Module map (changes)

```
io/
  fasta        NEW  streaming, multi-record, gz/bgzf-transparent FASTA reader → core types
  fastq        NEW  streaming, gz/bgzf-transparent FASTQ reader → core types
  decompress   NEW  magic-sniffing Read wrapper (plain / gzip / bgzf-sequential) + open_input/open_output (- = std)
  bam, vcf     KEEP (vcf: emit @SQ/##contig for all contigs from ContigSet)
genomics/
  fm_index     EXTEND  build_multi over the concatenated genome; factor FM ops over a `BwtBacking` trait
  index/
    format     EXTEND  add section kinds (BwtBlocks, RankCheckpoints, Boundaries, CTable, Reference2bit, real SaSamples); contig global_offset
    io         REWRITE  serialize the queryable FM-index; FmIndexView<'a> over the mmap; stale-index guard
  bwt_aligner  REWRITE  generic over `B: BwtBacking` (owned index for build/tests; FmIndexView for production); contig-aware locate → Locus; reject boundary-straddling hits
  suffix_array, rank_select, compressed_dna  KEEP kernels (rank_select/compressed_dna gain borrowed-slice read accessors used by FmIndexView)
core/
  locus        KEEP  ContigSet already 64-bit-offset-ready; add global_offset→Locus resolve (binary search) + boundary-cross check helper
provenance/
  mod          EXTEND  RunManifest records the index artifact path + BLAKE3
main.rs        REWRITE (slim)  `rosalind index` subcommand; align/variants/somatic take --index; - conventions; reading removed (now in io/)
docs/
  index-format.md  NEW  the on-disk format as a public, forkable artifact
```

**Public-API discipline (unchanged hard rule):** no `rust_htslib` types in public signatures; readers
and the index API return/borrow `core` types; `anyhow` only at the CLI.

**New dependency:** `flate2` (default `miniz_oxide` backend — pure Rust, portable/edge-friendly) for
gzip/bgzf-sequential decompression. bgzf *random access* is not pulled in here (Phase D).

---

## 5. The on-disk index format (public, forkable artifact)

A single little-endian file. The header (fixed size, already implemented) is readable without parsing
the body and carries magic, version, endian marker, `contig_count`, `sa_sample_rate`, the section-table
location, and the **reference BLAKE3**. The body is a set of length-delimited sections located by a
section table; **section payloads are flat, naturally-aligned, fixed-width LE arrays so the read-side
view can index into them directly from the mmap (zero-copy).**

Sections (extending today's `Contigs` / `Reference` / `SaSamples`):

| Section | Payload (v1) |
|---|---|
| `Contigs` | per contig: `name_len:u32, name:[u8], length:u64, global_offset:u64` (offset stored, not re-derived, so the view is self-describing) |
| `Reference2bit` | the concatenated forward reference, 2-bit packed (A/C/G/T) + an N/ambiguity bitmask (lossless), for `ref_base` lookups without a separate `.fa` |
| `BwtBlocks` | the 2-bit BWT, block-partitioned (`block_size` from the header), exactly the bytes `CompressedDNA` already produces |
| `RankCheckpoints` | per-block `Occ`/rank-select payload (`RankSelectIndex` data as flat arrays) |
| `Boundaries` | per-block `{start:u64, cumulative_counts:[u32;5], sentinel_count:u32}` (the `CompressedBoundaries` table) + `c_table:[u32;6]` + `sentinel_pos:u64` |
| `SaSamples` | **compact/sparse** (decision §3.6): a sampled-position bitvector over BWT indices + its rank checkpoints + a packed array of the sampled `u32` SA values — *not* the dense `bwt_len × u32` array, and not the current scaffold's zero samples |

**Determinism:** SA-IS is deterministic; sections are written in a fixed order with fixed-width LE
encoding; no timestamps appear anywhere in the file. ⇒ byte-identical across builds of the same
reference. **Integrity:** the header already stores the reference BLAKE3; `IndexReader::open` validates
magic/version/endian/header-size/section-table and (cheaply) that section extents lie within the file.

`docs/index-format.md` documents all of the above as a stable, forkable contract so a third party can
write a compatible builder or reader.

---

## 6. Data flow after Phase B

```
# build once
rosalind index ref.fa[.gz] [-o ref.fa.rosalind.idx]
  io::fasta (stream, multi-record) -> ContigSet + concatenated text (err if len >= 2^32)
  -> BlockedFMIndex::build_multi (SA-IS, single terminal sentinel, real SA samples)
  -> genomics::index serializer (deterministic, flat LE sections, BLAKE3) -> ref.fa.rosalind.idx

# align (streaming reads -> SAM/BAM)
rosalind align --index ref.fa.rosalind.idx --reads R.fq[.gz] [--output - | out.bam]
  IndexReader::open -> mmap -> FmIndexView (queryable; NO rebuild)
  io::fastq (stream, gz) -> seed -> chain -> extend -> resolve Locus -> reject cross-contig
  -> core::AlignedRead{contig,pos,..} -> SAM/BAM (@SQ per contig) -> stdout/-

# call (whole-genome pileup -> VCF)
rosalind variants --index ref.fa.rosalind.idx [--alignments - | sorted.bam]
  BamSource -> for each contig in ContigSet: PileupEngine(region = 0..contig.len)
            -> call::germline -> io::vcf (##contig per contig)
  ref_base served from the index's Reference2bit (no separate --reference)
  + RunManifest records index path + BLAKE3
```

---

## 7. Multi-contig FM-index design

- **Concatenation, single sentinel.** All contig sequences are concatenated (in `ContigSet` id order)
  into one text; a single terminal `$` sentinel is appended (exactly as today); SA-IS runs once over
  the whole text. **No per-contig separators** — keeping the kernel unchanged.
- **Coordinate mapping.** A located global SA position `p` resolves to `Locus{contig, pos}` by binary
  search over the contigs' `global_offset` (`pos = p - global_offset`). `ContigSet` already stores
  `global_offset: u64`.
- **Boundary rejection (the bwa "bns" model).** A hit of read length `L` at global `p` is **rejected**
  if `[p, p+L)` crosses a contig boundary (i.e., spans into the next contig). This is the only thing
  preventing phantom cross-contig matches; the SA itself may contain boundary-straddling suffixes,
  which are simply never reported. Covered by a dedicated test (a read placed to straddle a boundary
  yields no hit).
- **u32 ceiling.** Concatenated length must be `< 2^32`; `build_multi` returns a typed error otherwise.
  Documented in README + `index-format.md`.
- **Reference bases.** The forward 2-bit reference (with N/ambiguity mask) is stored so the pileup can
  read `ref_base` at any locus without a separate FASTA.

---

## 8. Zero-copy loading design

- **`trait BwtBacking`** — the minimal surface the FM algorithm needs: BWT symbol/word access, per-block
  rank checkpoint access, boundary/c_table access, and compact SA-sample access (the sampled-position
  bitvector + its rank + the packed sampled values, per §3.6). The FM operations (`backward_search`,
  `rank`, `lf_index`, `sa_at`, `locate_interval`) are written **once** against this trait; `sa_at` walks
  LF until it lands on a sampled BWT position (rank over the sampled bitvector → index into the packed
  values).
- **Build-side backing** = the existing owned `BlockedFMIndex` data (Vecs in RAM). Used by `rosalind
  index` and by unit/property tests.
- **Read-side backing** = `FmIndexView<'a>`, holding `&'a [u8]` subslices of the mmap (BWT words, rank
  checkpoints, boundaries, c_table, sa_samples) + the parsed `ContigSet`. `ReferenceIndex` owns the
  `MmapReadOnly`; `ReferenceIndex::view(&self) -> FmIndexView<'_>` borrows from it. **No allocation of
  the big arrays; the OS pages them.**
- **`BWTAligner` becomes generic** over `B: BwtBacking`: production wires it to `FmIndexView`;
  build/tests wire it to the owned index. Contig-aware locate returns `Locus`.
- **Equivalence gate (the headline correctness test):** build an index in RAM, serialize it, mmap it as
  an `FmIndexView`, and assert that `backward_search` / `rank` / `sa_at` / `locate_interval` return
  **identical** results for a large, fixed battery of patterns. This simultaneously proves zero-copy
  correctness *and* serialization round-trip fidelity.
- **No-rebuild gate:** the read path must not call SA-IS. Enforced structurally (the loader constructs a
  view, never a `BlockedFMIndex::build*`) and asserted in an e2e timing/behavioral check.

---

## 9. Decomposition — six stages, each lands green, dependency-ordered

Each stage is independently compilable, testable, and green; later stages depend only on earlier ones.
Each becomes its own implementation plan (writing-plans) and PR, mirroring the Phase A A1–A7 cadence.

### B1 — Streaming readers (`io/fasta`, `io/fastq`, `io/decompress`)
- Streaming, multi-record FASTA/FASTQ over any `Read`; magic-sniffing transparent gz/bgzf-sequential
  decompression; `-` = stdin. Library-first (typed errors, `core` types, no htslib/anyhow). Ingestion
  removed from `main.rs`.
- **DoD:** a 3-record FASTA yields a 3-contig `ContigSet`; `.gz` and plain read identically (byte-for-
  byte equal records); reading from `-` works; round-trip + malformed-input tests; CI reads a gz toy.
- **Why first:** hard prerequisite for index build and for gz reads.

### B2 — Multi-contig FM-index build + `Locus` mapping (in-RAM)
- `BlockedFMIndex::build_multi(ContigSet, &concat)`; global→`Locus` resolve; cross-boundary rejection;
  contig-aware locate. Still in-RAM (no persistence yet) to **isolate multi-contig correctness from
  serialization.**
- **DoD:** 2-contig reference — a read maps to the correct `(contig,pos)`; a boundary-straddling read is
  rejected; exact-match correctness vs naive on both contigs; deterministic; `build_multi` errors at
  `len >= 2^32`.

### B3 — Zero-copy persistence: format + serializer + `FmIndexView` + `rosalind index`
- Extend the section set (§5); serialize the queryable FM-index deterministically; implement
  `BwtBacking` + `FmIndexView<'a>`; add the `rosalind index` subcommand.
- **DoD (headline gates):** view == in-RAM index (equivalence battery); index file byte-identical across
  builds; load performs no SA-IS; basic load-residency check (index not resident); **SA samples stored
  compactly** (the index file scales as `O(reference + bwt_len/rate)`, not `O(bwt_len)` in the sample
  array). Round-trip test over a multi-contig reference.

### B4 — Wire consumers + whole-genome pileup walk
- `align`/`variants`/`somatic` take `--index`; load `FmIndexView`; drop the single-contig blocker and
  the per-run `BWTAligner::new(reference)` rebuild. Pileup iterates per contig across the genome.
  `@SQ`/`##contig` for all contigs. Stale-index BLAKE3 guard. `variants`/`somatic` drop the
  `--reference` requirement (served from the index).
- **DoD:** multi-contig FASTA → index → align → sort → variants → VCF with calls on ≥2 contigs, correct
  `##contig`, deterministic; stale index rejected with a clear message.

### B5 — Pipe-native composition + receipt covers the index
- `-` conventions for align reads/output and variants alignments; `open_input`/`open_output` helpers.
  Sort stays spill-to-disk internally but reads/writes streams at its ends; `variants` documents that
  stdin must be coordinate-sorted. `RunManifest` records the index path + BLAKE3.
- **DoD:** `align - | sort - | variants -` equals the file-based VCF; receipt includes the index hash;
  determinism holds across the pipe.

### B6 — Minimal per-pillar CI + docs-truth (scoped to what B ships)
- CI: index round-trip/determinism job; multi-contig gz e2e (index→align→sort→variants on a ≥2-contig
  toy); `bcftools view`/`samtools quickcheck` when available. README "Current scope" updated
  (multi-contig, persisted index, gzip); `docs/index-format.md` published.
- **DoD:** new CI jobs green; README no longer claims single-contig/uncompressed-only; format
  documented. (Clippy `-D` + theory gating remain Phase E.)

---

## 10. Testing strategy (TDD — failing test first per stage)

- **Reader fidelity (B1):** multi-record FASTA/FASTQ parsing; gz vs plain byte-equality; `-`/stdin;
  truncated/malformed records produce typed errors, not panics.
- **Multi-contig correctness (B2):** map-to-correct-contig; boundary-straddle rejection; exact-match vs
  naive over a 2-contig reference; `len >= 2^32` error path (small synthetic via a stubbed length
  check).
- **Zero-copy equivalence + determinism (B3):** the equivalence battery (view vs in-RAM, many
  patterns); index file byte-identical across two builds; section-extent validation; corrupt/truncated
  index rejected.
- **End-to-end multi-contig (B4):** 2-contig gz reference + reads spanning both → align → sort →
  variants → expected multi-`##contig` VCF; stale-index guard fires on reference mismatch.
- **Pipe equivalence (B5):** piped path VCF == file-based path VCF; receipt determinism incl. the index
  hash.
- **Property/regression carried from A:** existing determinism/space-bounds/fm-index-props/golden-vcf
  suites stay green; `golden_vcf` extended for multi-contig.

---

## 11. Risks & mitigations

- **Zero-copy refactor is the largest piece.** *Mitigation:* the build-side/read-side split avoids a
  risky in-place rewrite of `BlockedFMIndex`; the `BwtBacking` trait keeps one FM-algorithm source of
  truth; the equivalence test makes any divergence loud. B2 proves multi-contig correctness *before* B3
  touches serialization.
- **Alignment-quality scope creep.** *Mitigation:* §2 hard boundary — keep the existing seed→chain→
  extend math; only make it multi-contig + index-backed. MAPQ/soft-clip quality is explicitly a later
  phase.
- **bgzf vs gzip confusion.** *Mitigation:* B reads both *sequentially* via `flate2::MultiGzDecoder`
  (bgzf is concatenated gzip members); random-access bgzf (virtual offsets, `.gzi`/`.csi`) is Phase D
  and explicitly out.
- **u32 SA ceiling surprising a user with a >4.29 Gbp target.** *Mitigation:* a typed build error +
  documentation in README and `index-format.md`; the `u64` coordinate API is preserved so a future u64
  SA is a format-version bump, not an API break.
- **mmap memory story overclaiming.** *Mitigation:* document that the index is OS page-cache (not the
  budgeted working set); the load-residency check guards against accidental full-read; the enforced RSS
  gate lands in D.
- **`--reference` removal is a CLI change.** *Mitigation:* note in CHANGELOG; `variants`/`somatic` read
  the reference from the self-contained index; keep a clear error if neither an index nor a reference is
  available.
- **Index *build* is memory-heavy on large genomes** (SA-IS workspace + the SA itself are `O(reference)`;
  the dense sample array is removed per §3.6, but SA-IS over a human genome is still a big-machine, not a
  laptop, operation). *Mitigation:* this is the accepted "build is O(reference), documented separately"
  caveat — the promise is that *using* (mmap-load + query) the index is bounded. Document the build
  footprint in `index-format.md`; a streaming/low-memory SA construction is explicitly out of Phase B
  scope.

---

## 12. CI additions (minimal, per-pillar — full hardening is Phase E)

- **Index determinism job:** build the toy index twice; assert byte-identical files; assert
  `IndexReader::open` round-trips.
- **Multi-contig gz e2e:** generate a ≥2-contig gzipped toy; `index` → `align` (from `-`) → `sort` →
  `variants`; assert calls on ≥2 contigs and a `##contig` per contig.
- **Standards validation (best-effort):** `samtools quickcheck` the BAM and `bcftools view` the VCF when
  the tools are present on the runner.

Clippy `-D warnings` across the legacy tree, theory-layer feature-gating, and the `lib.rs` crate-doc
honesty pass remain **Phase E**.
