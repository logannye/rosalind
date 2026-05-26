# Phase A — The pileup kernel + calling vertical + reproducibility receipt

**Status:** Spec for review — 2026-05-26. First vertical slice of the target architecture
(`2026-05-26-rosalind-target-architecture.md`).

**One-line goal.** Stand up Rosalind's defining primitive — a single CIGAR-aware, filtered,
strand-aware, **budget-honoring, read-length-agnostic streaming pileup kernel** over
target-shaped `core` types — make variant calling its *first consumer* with **calibrated,
abstention-aware** confidence, emit **spec-valid** VCF, and wrap every run in a **verifiable
reproducibility receipt**. On real coordinate-sorted reads (reverse-strand, indel-bearing,
soft-clipped, short or long), Rosalind produces correct, deterministic, trustworthy output that
builders can also consume as a library.

This slice is deliberately broader than a bug-fix sprint: it lays the `core` types and the
public kernel API that every later phase builds on, and it demonstrates four of the five product
promises (deterministic, standards, verifiable reproducibility, embeddable) end-to-end on day
one. Guided by the core principle: *most uniquely valuable to edge builders.*

---

## 1. Success criteria (definition of done)

- **Correct:** reverse-strand reads contribute correct alleles; indel/soft-clip/refskip reads
  contribute matched bases at correct reference coordinates; **long reads** (indel-rich) handled
  identically to short reads. No silently-dropped reads — skips are counted by reason.
- **Kernel/substrate:** the pileup engine is a public, documented, consumer-agnostic API
  yielding `PileupColumn`s; no htslib types leak into public signatures; the same stream that
  feeds the caller can feed a plugin. A worked example consumes columns without the caller.
- **Calibrated & honest:** germline calls carry true Phred `QUAL`/`GQ` from a genotype-likelihood
  model and an explicit FILTER vocabulary; below-evidence sites are flagged/abstained, never
  emitted as confident calls.
- **Standards:** germline + somatic VCF pass `bcftools view` (CI when available), with
  `##fileformat`, `##contig`, typed `##INFO`/`##FORMAT`/`##FILTER`, and `FORMAT` + sample
  column(s) (germline 1 sample; somatic `TUMOR`/`NORMAL`).
- **Verifiable:** every `variants`/`somatic` run emits a JSON receipt (tool version, normalized
  input paths, BLAKE3 content hashes of inputs + params + outputs); a second run reproduces it
  byte-for-byte.
- **Budget-aware:** the engine exposes `working_set_bound()` and accepts a `MemoryBudget`; the
  active-set bound is documented (foundation for `rosalind plan`).
- **Deterministic:** identical inputs → byte-identical VCF + receipt; shuffled reads → identical
  calls. Public `call_variants(Vec<AlignedRead>)` API preserved.
- `CompressedEvaluator`/`PileupProcessor` off the calling path (retained for plugins).

---

## 2. Scope

**In:** `core/` (locus, sequence, record, **budget**, error); the public `pileup/` engine;
`call/germline` (calibrated diploid biallelic SNV GL model) + re-homed `call/somatic`;
`io/vcf` (spec-valid writer); a **minimal `provenance/`** receipt for `variants`/`somatic`; the
correctness fixes (§7); the migration (§8); a worked "plugin consumes the column stream"
example; long-read correctness + tests.

**Out (deferred, by design):** germline **indel calling** (C) — indel reads are *handled* by
the pileup, not *called*; multi-allelic decomposition (biallelic only); MAPQ *calibration*
(D/align) — A applies the filter with existing MAPQ; full `rosalind plan`/`verify` subcommands
and region/`.csi` fetch + whole-genome multi-contig iteration (B/D) — A uses the single-region
model on contig-aware types; the Python binding (Phase P) — but the API is shaped for it now;
mate-overlap dedup, strand-bias *filter*, BAQ (C) — A collects per-strand counts so the filter
is a drop-in; rayon (D).

---

## 3. New & changed modules

### 3.1 `core/` (NEW — the lingua franca)
- **`locus`**: `Contig {id, name, length}`, `Position(u32)`, `Locus {contig, pos}`; `ContigSet`
  (name↔id, 64-bit-safe global offsets for B). Phase A is single-contig but uses these types.
- **`sequence`**: 2-bit DNA (from `compressed_dna`); non-ACGT folds to `N` with a *counted
  warning* (not a hard error); ambiguity mask preserved (lossless serialize). `BaseCode` is the
  single source of truth (today duplicated in 4 files).
- **`record`**: canonical `AlignedRead` — `contig, pos, mapq, flags (SAM bitflags), cigar:
  Vec<CigarOp>, seq, qual`. `CigarOpKind` gains `RefSkip`. `ref_span()/end()` CIGAR-derived.
  SEQ forward-oriented; strand is metadata, never applied to bytes. **Read-length-agnostic.**
- **`budget`**: `MemoryBudget` (declared cap) + a `WorkingSet` model streaming stages report.
- **`error`**: `thiserror` taxonomy.

### 3.2 `pileup/` (REWRITE — the public kernel)
- **`PileupParams { min_mapq, min_base_qual, skip_secondary, skip_supplementary,
  skip_duplicate }`**; optional `MemoryBudget`.
- **`ReadSource`** trait yielding coordinate-sorted `AlignedRead`s; adapters: `BamSource` (lazy
  `bam::Record → AlignedRead`, bounded — conversion happens *here*, htslib stays internal) and
  `SliceSource` (pre-sorted `&[AlignedRead]`).
- **`PileupEngine`** — `Iterator<Item = Result<PileupColumn>>`. Responsibilities: CIGAR
  projection (Match/Equal/Diff → ref+read/emit base; Ins/SoftClip → read-only; Del/RefSkip →
  ref-only/no obs; HardClip → neither); read-level filtering with a `SkipCounts` accumulator
  (unmapped/wrong_contig/secondary/supplementary/duplicate/low_mapq/malformed); per-observation
  strand; deterministic active-set + intra-column observation order; **bounded active set**
  (expired by CIGAR-derived `end`); **empty stretches advance by loop, not recursion** (fixes
  `pileup_stream.rs:192`); `working_set_bound()` reporting; budget honored or fail-fast.
- **`PileupColumn { locus, ref_base, obs: Vec<Obs> }`**, `Obs { allele: u8, bq: u8, reverse:
  bool, mapq: u8 }`; derived `allele_counts()`, `strand_counts()`, `depth()`. **This is the
  public substrate type** — stable, documented, designed for callers, plugins, and (Phase P)
  Python/ML consumers. A `tests/`-level example shows a non-caller consumer (e.g. a coverage/
  per-base-feature reducer) over the same stream.

### 3.3 `call/` (NEW germline; re-home somatic) — calibrated & abstention-aware
- **`call::germline(column, params) -> GermlineSite`** — the GL model (§5). Returns a *site
  decision* including possible abstention: a confident variant (`PASS`), a low-confidence
  variant (`LowQual`/`LowDepth`), or no-call. Honest by default: thin evidence is flagged, not
  emitted as a confident call.
- **`call::somatic(tumor, normal, params) -> Option<SomaticCall>`** — existing exact f64
  binomial LLR moved here verbatim, fed corrected columns; tumor/normal co-walk driver over two
  `PileupEngine`s.
- `GermlineCall`/`SomaticCall` carry GT/GQ/PL/AD/DP/QUAL/FILTER (+ tumor/normal AD/DP/AF).

### 3.4 `io/vcf` (REWRITE — one spec-valid emitter)
- `VcfHeader` builder: `##fileformat=VCFv4.2`, `##contig` from `ContigSet`, typed
  `##INFO=<DP,AF>`, `##FORMAT=<GT,GQ,DP,AD,PL>`, `##FILTER=<PASS,LowQual,LowDepth>`, then
  `#CHROM…FORMAT\t<sample(s)>`. `VcfRecord` model + writer; germline 1 sample, somatic
  `TUMOR`/`NORMAL`. Deterministic ordering preserved.

### 3.5 `provenance/` (NEW — minimal receipt)
- `RunManifest { tool_version, subcommand, inputs:[{path_normalized, blake3}], params,
  outputs:[{path, blake3}], created_utc? }` serialized as canonical JSON (sorted keys, no
  timestamps in the hashed core). `variants`/`somatic` write `<output>.manifest.json`. BLAKE3
  is already a dependency. Full `provenance` module + `rosalind verify` is Phase D; A emits the
  receipt and a determinism test hashes it.

---

## 4. Data flow
- **Germline, BAM:** `BamSource → PileupEngine → call::germline → VcfWriter (+ manifest)`.
- **Germline, in-memory/SAM:** load reads → **sort by pos** → `SliceSource` → same engine →
  same writer. `StreamingVariantCaller` re-implemented over the engine; public
  `call_variants(Vec<AlignedRead>)` preserved.
- **Somatic:** two engines (tumor, normal) co-walked by position → `call::somatic` → writer with
  two sample columns (+ manifest).
- **Substrate demo:** `BamSource/SliceSource → PileupEngine → custom reducer` (no caller),
  proving the kernel is consumer-agnostic.

---

## 5. The genotype model (minimal, standard, calibrated, abstention-aware)

Reference `R`, candidate alt `A` = most-supported non-ref allele (biallelic). Per observation
base `bᵢ`, base-quality `qᵢ`, error `εᵢ = 10^(−qᵢ/10)`:
```
P(bᵢ | X)      = 1 − εᵢ  if bᵢ == X else εᵢ/3
P(bᵢ | hom R)  = P(bᵢ | R)
P(bᵢ | hom A)  = P(bᵢ | A)
P(bᵢ | het RA) = ½·P(bᵢ | R) + ½·P(bᵢ | A)
```
Accumulate log-likelihoods (f64, fixed order) for `{0/0, 0/1, 1/1}`; apply a configurable prior
(default heterozygosity θ = 1e-3: `P(0/1)=θ`, `P(1/1)=θ/2`, `P(0/0)=1−1.5θ`).
- **PL**: Phred `−10·log10 L(g)`, normalized min=0, cap 255.
- **GT** = argmax posterior (argmin PL). **GQ** = second-smallest PL, cap 99. **QUAL** =
  `−10·log10 P(0/0 | data)` (true Phred site-is-variant). **AD/DP** as counts.

**Abstention/honesty rules (the unique-value posture):**
- Emit a record only when `GT ≠ 0/0`. `FILTER = LowDepth` if `DP < min_depth`; `LowQual` if
  `QUAL < min_qual`; else `PASS`.
- Defaults tuned to *not overcall*: when evidence is insufficient the site is filtered, not
  confidently called. Confidence is calibrated (PL/QUAL derived from the model, not a heuristic),
  so a downstream/field user can trust the numbers and threshold for their risk tolerance.
- Determinism: f64, fixed summation order; identical across runs.

---

## 6. VCF output contract (germline example)
```
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">
##FILTER=<ID=LowQual,Description="QUAL below threshold">
##FILTER=<ID=LowDepth,Description="Depth below threshold">
#CHROM POS ID REF ALT QUAL FILTER INFO        FORMAT         SAMPLE
chr1   101 .  A   G   48.0 PASS    DP=30;AF=0.5 GT:GQ:DP:AD:PL 0/1:48:30:15,15:48,0,49
```
Somatic adds `TUMOR`/`NORMAL` columns and retains tumor/normal AF/DP.

---

## 7. Correctness fixes folded in
1. **Reverse strand**: use forward SEQ directly; delete complement+offset-flip
   (`pileup_stream.rs:33-42`, `somatic/model.rs:325-327`).
2. **CIGAR-aware projection** (§3.2): stop dropping non-simple-match reads; `end()` = ref span.
3. **Read filters** in-engine.
4. **Empty-position recursion → loop** (`pileup_stream.rs:192`).
5. **`checked_sub`** on offsets (`variant_caller.rs:125`).
6. **Skip-and-count** replaces the SAM-path `bail!` and silent `continue`s; summary logged.

---

## 8. Migration / deprecation
- `StreamingVariantCaller` re-implemented over `PileupEngine`; `call_variants(Vec<AlignedRead>)`
  preserved.
- `statistics::bayesian_variant_caller` **deleted** → `call::germline`.
- `main.rs::parse_cigar` generalized → full CIGAR parser.
- `somatic/model.rs` LLR + co-walk → `call/`; reverse-complement bug removed.
- `PileupProcessor`/`PileupSummary`/`PileupWorkload` + `CompressedEvaluator` retained,
  re-exported, **documented plugin-facing**, off the calling path.
- `vcf.rs` + `somatic/vcf.rs` string writers replaced by `io/vcf`.

---

## 9. Error handling
Typed library errors; `anyhow` only at the CLI. **Lenient by default** (malformed/filtered reads
counted + skipped, never abort); IO/parse failures propagate. (`--strict` later.)

---

## 10. Testing (TDD — failing test first)
- **Reverse-strand regression** (± strand over a known SNV → identical correct counts) — written
  first, proves the bug before any change.
- **CIGAR projection** (I/D/S/N land correctly; clipped/inserted excluded; previously-dropped
  reads now contribute).
- **Long-read correctness** (a long, indel-rich read piles up identically to its short-read
  equivalents) — guards read-length-agnosticism.
- **Filtering** (secondary/supplementary/duplicate/low-MAPQ excluded; `SkipCounts` exact).
- **Genotype/calibration** (hand-computed PL/GT/GQ/QUAL for hom-ref (20,0), het (~15,15),
  hom-alt (0,20); **abstention** at low depth (2,1) → LowDepth not a confident call;
  QUAL monotonic with evidence).
- **VCF validity** (parse + header tags + FORMAT/sample well-formed; `bcftools view` round-trip
  in CI when available).
- **Substrate** (the non-caller column-stream example produces expected per-base output).
- **Determinism + receipt** (twice → byte-identical VCF *and* manifest; shuffled reads →
  identical calls).
- **End-to-end smoke** (synthetic ref + ± strand + indel + a long read → align → sort →
  variants → expected VCF + manifest). Update `golden_vcf`.

---

## 11. Risks & mitigations
- *Reverse-strand fix rests on "BAM SEQ is forward-oriented."* True per SAM spec; **the failing
  regression test is written first** to prove it on this code before changing anything.
- *Public API stability* — preserve `call_variants`; cover with a test; public-API surface test
  asserts no htslib types leak.
- *Somatic two-sample VCF is downstream-visible* — note in CHANGELOG (more conformant).
- *Scope creep into C* — hard boundary: biallelic SNV genotyping only; indels/multi-allelic/
  SB-filter explicitly out. The pileup *handles* indel/long reads; the caller does not *emit*
  indels.
- *Manifest determinism* — canonical JSON (sorted keys), no timestamps in the hashed core;
  covered by the receipt determinism test.
