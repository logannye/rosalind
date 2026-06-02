# Finding — germline caller detection accuracy (simulated diploid truth)

**Date:** 2026-06-02. **Harness:** `tests/germline_accuracy.rs` (deterministic, fixed-seed, gated in
CI). **Comparator:** the existing `compare_callsets` (detection by normalized pos+ref+alt) + the new
`rosalind eval-germline` CLI. Addresses the reflection audit's #1 blind spot: the memory contract was
validated, but the **calls themselves** had never been measured against ground truth.

## What was measured

A diploid donor is built by injecting ~40 known het/hom SNVs into a 20 kbp reference; reads are
sampled 50/50 from the two haplotypes (with sequencing error) and written **directly as a sorted
BAM** — bypassing the toy aligner, so this isolates the **caller's** detection accuracy (genotype
model + the unbiased depth cap), not alignment. PASS = records passing the default filters
(`min_qual 30`, `min_depth 8`); `all` = every emitted variant record.

| Regime | metric | truth | calls | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|---|---|
| 40× / 0.5% err (clean) | **PASS** | 40 | 40 | 40 | 0 | 0 | **1.000** | **1.000** | **1.000** |
| 40× / 0.5% err (clean) | all | 40 | 40 | 40 | 0 | 0 | 1.000 | 1.000 | 1.000 |
| 12× / 1.5% err (stress) | **PASS** | 40 | 53 | 36 | 17 | 4 | **0.679** | **0.900** | **0.774** |
| 12× / 1.5% err (stress) | all | 40 | 138 | 40 | 98 | 0 | 0.290 | 1.000 | 0.449 |

## Honest reading

- **On clean data (40× / 0.5%) the caller is exact**: every het and hom-alt SNV is recovered with zero
  false positives. The diploid GL model + the θ=1e-3 prior correctly abstain on error-driven thin alt.
- **At low coverage + high error (12× / 1.5%) detection degrades, and the quality filter matters.**
  Unfiltered, the caller finds *every* true variant (recall 1.000) but emits 98 error-driven false
  positives (precision 0.290). The default PASS filter (`min_qual 30` / `min_depth 8`) cuts those FPs
  to 17 (precision → 0.679) at the cost of 4 true variants flagged low-quality (recall → 0.900). So
  the caller has a **working precision/recall knob**, and the honest characterization is: *robust
  recall, precision that is filter- and coverage-dependent.*
- **The depth-cap fix holds in an accuracy frame.** A deep het site (alt reads starting *at* the
  variant, under a cap below local depth) is recovered as a PASS call — the unbiased reservoir, which
  the old leftmost-arrival cap would have silently dropped.

## Scope + honest limitations

- **Simulated, not GIAB.** Reads are simulated and the BAM is built from known coordinates (no real
  alignment). This measures the caller on clean, indel-free, single-contig data — it does **not** cover
  alignment error, mapping ambiguity, repeats, or real error profiles. A real **GIAB HG002** run is the
  natural next step and plugs into the *same* comparator via `rosalind eval-germline --reference
  --calls --truth --regions <highconf.bed>`; it needs a downloaded pre-aligned BAM slice (out of scope
  here — this environment has no bwa/samtools and GIAB BAMs are large).
- **Detection only.** Matching is pos+ref+alt; **genotype concordance** (het vs hom) is not scored — a
  documented future refinement.
- **SNV only.** The caller is biallelic-SNV; indels are out of scope (they need local realignment,
  which is not a per-column streaming operation).

## Reproduce

```
cargo test --test germline_accuracy -- --nocapture
```

The numbers above are deterministic (fixed-seed LCG). The CI gate asserts a margin below the measured
PASS values (clean ≥ 0.97/0.97; stress recall ≥ 0.85, precision ≥ 0.60) so a real regression fails
loudly without flaking.
