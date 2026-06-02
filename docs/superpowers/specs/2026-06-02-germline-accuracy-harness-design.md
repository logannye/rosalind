# Simulated-Truth Germline Accuracy Harness (design)

**Status:** DESIGN SPEC — 2026-06-02. **Branch:** `rosalind/germline-accuracy` (off `main` `95c4ba8`).
Addresses the reflection audit's #1 blind spot: the memory contract is validated, but the variant
**calls** have never been measured against ground truth (only simulated E. coli with no injected
variants). User chose "simulated diploid truth" as the first accuracy increment.

## 1. Goal

Measure the germline caller's **detection precision / recall / F1** against a known truth set,
deterministically and in-house (no downloads), and gate it in CI. Reuse the existing GIAB-shaped
comparator (`compare_callsets` + `BedIndex` + `VcfVariant`) so a real GIAB BAM slice can plug into the
same interface later. Specifically validate the depth-cap fix in an accuracy frame (a deep het site
whose alt reads start at the variant must be recovered).

## 2. What already exists (reuse, don't rebuild)

- `compare_callsets(reference: &[u8], calls: &[VcfVariant], truth: &[VcfVariant], bed) ->
  ComparisonReport` with `true_positive`/`false_positive`/`false_negative`, `precision()`,
  `recall()`, and a per-`VariantType` breakdown. Matching is **detection by normalized (pos, ref,
  alt)** — genotype is not compared (the standard primary metric; genotype concordance is a later
  refinement).
- `VcfVariant { chrom: String, pos0: u32, reference: Vec<u8>, alternate: Vec<u8>, qual, filter }`,
  `read_vcf_variants`, `BedIndex::from_str` — `src/genomics/eval/`.
- `run_eval_somatic` (`src/main.rs:400`) is already VCF-agnostic comparison logic.
- BAM construction from hand-built `Record`s (the pattern in `tests/truthset_validation.rs`).
- `call_germline_whole_genome` (library) returns `(Locus, ref_base, GermlineCall)` rows;
  `GermlineCall { genotype, alt_base: u8, qual, filter, … }`.

## 3. Deliverables

### 3a. `eval-germline` CLI (the GIAB-ready interface)

Generalize `run_eval_somatic` → a shared `run_eval(reference, calls, truth, regions)`; keep
`Commands::EvalSomatic` calling it; add `Commands::EvalGermline { reference, calls, truth, regions }`
calling the same `run_eval`. Add **F1** to the printed report:
`f1 = 2·P·R / (P + R)` (0 when `P + R == 0`). A real GIAB germline run uses exactly this command
(`rosalind eval-germline --reference grch38.fa --calls out.vcf --truth giab.vcf --regions
highconf.bed`).

### 3b. `tests/germline_accuracy.rs` — the gated accuracy harness

A deterministic, download-free integration test. Bypasses the toy aligner (builds the BAM directly
from known coordinates) so it measures the **caller**, not alignment.

**Simulation (fixed-seed LCG — no `rand` dependency):**
- Reference: `N = 20_000` bp of deterministic pseudo-random DNA (LCG over `ACGT`). Single contig
  `chr1`.
- Truth variants: `~40` SNVs at positions spaced ~`N/45` apart, avoiding the first/last 200 bp. Each
  is `Het` or `HomAlt` (assigned by the LCG). `alt_base = next base after ref_base in "ACGT"`
  (guarantees alt ≠ ref).
- Two haplotype arrays built once: `hap1 = reference + HomAlt variants`; `hap2 = reference + HomAlt +
  Het variants`. (So Het sites carry alt on hap2 only; HomAlt on both.)
- Reads: length `L = 100`, target coverage `~40×` via tiled starts with small LCG jitter. Each read is
  sampled from `hap1`/`hap2` 50/50; `SEQ = hap[start..start+L]` with each base flipped to a random
  other base with probability `E = 0.005` (sequencing error). Record: tid 0, `pos = start`, CIGAR
  `L M`, qual `40`. All forward strand (strand is metadata; the reverse path is covered elsewhere).
- **Deep-het probe (the depth-cap fix in an accuracy frame):** pick one Het truth site; additionally
  emit `~30` ref-bearing reads that START upstream and span it, and `~30` alt-bearing reads that START
  AT the site — so with a depth cap below the local depth, the biased (old) cap would drop the
  at-variant alt reads and miss the het, while the unbiased reservoir recovers it. The test runs with
  `PileupParams.max_depth = Some(50)` so normal ~40× sites do not cap but the deep-het site (~60×+)
  does.

**Call + compare:**
- Write the reads to a BAM, `sort_bam_deterministic` → sorted BAM. Build a Rosalind index for the
  reference (`IndexWriter`/`IndexReader`), get `ReferenceView` + `ContigSet`.
- `call_germline_whole_genome(StreamingBamSource, …, PileupParams{ max_depth: Some(50), .. },
  GermlineParams::default(), sink)` collecting rows.
- Convert each emitted row → `VcfVariant { chrom: "chr1", pos0: locus.pos.0, reference:
  vec![ref_base], alternate: vec![call.alt_base], qual: Some(call.qual as f32), filter: "." }`.
  (Compare ALL emitted variant records = detection recall; the caller already abstains on hom-ref.)
- Truth → `VcfVariant` from the injected list.
- `report = compare_callsets(&reference_bytes, &calls, &truth, None)`.

**Assertions (thresholds CALIBRATED to the measured run, set a margin below — report the actuals):**
- `report.recall() >= R_MIN` and `report.precision() >= P_MIN` (initial targets `R_MIN = 0.90`,
  `P_MIN = 0.90`; adjust to the measured numbers, never paper over a low result — a low number is a
  real finding about the caller, to be surfaced).
- The deep-het site is in the called set (a TP) — the cap fix, proven on accuracy.
- A clear panic message printing the full report on failure (so CI shows the numbers).

## 4. Honesty + scope

- **Report the real numbers.** If precision is dragged down by error-driven FPs, that is a genuine
  finding about the diploid model's error-robustness; surface it (and `min_qual` is a legitimate knob
  to discuss), do not silently tune until green.
- Detection metrics only (pos+ref+alt); **genotype concordance** (het-vs-hom) is a documented future
  refinement, not in scope.
- Real **GIAB HG002** is out of scope here (needs a downloaded pre-aligned BAM slice; this env has no
  bwa/samtools and GIAB BAMs are large) — but the `eval-germline` CLI + the `VcfVariant`/BED interface
  make it a drop-in second increment.
- A short **findings doc** (`docs/findings/2026-06-02-germline-accuracy.md`) records the measured
  precision/recall/F1 and the simulation parameters.

## 5. Risks

- **FPs from sequencing errors (medium).** At 0.5%/40× the diploid θ=1e-3 prior should abstain on thin
  alt, but clustered errors could produce a few FPs. Mitigation: measure; if precision is low, report
  it honestly and consider a `min_qual`/`min_depth` in the harness (documented), not a silent fudge.
- **Edge effects (low).** Variants near contig ends or reads spilling past the end — keep variants ≥
  200 bp from both ends; reads clamp at the contig end (handled).
- **Threshold brittleness (low).** Fixed-seed determinism makes the measured numbers stable;
  thresholds are set a margin below the measured value, so normal runs don't flake.

## 6. Self-review

- **Coverage:** the comparator already exists; new = `eval-germline` CLI (3a) + the gated simulated
  harness (3b) + a findings doc. ✓
- **Type consistency:** `VcfVariant` fields match `compare_callsets`; `GermlineCall.alt_base: u8` and
  the row's `ref_base: u8` map to `Vec<u8>` alleles. ✓
- **No placeholders:** the one calibration step (R_MIN/P_MIN) is explicitly "measure then set a margin
  below, report actuals." ✓
- **Determinism:** fixed-seed LCG for reference, variants, haplotype/read/error sampling →
  reproducible numbers, non-flaky gate. ✓
- **Scope:** detection-only, simulated-only; GIAB + genotype concordance explicitly deferred. ✓
