# Genotype-aware germline accuracy (the GIAB-grade comparator)

**Date:** 2026-06-02. Sprint 2.1. Companion to
[`2026-06-02-germline-accuracy.md`](2026-06-02-germline-accuracy.md) (the detection-only predecessor).

## What changed

`eval-germline`'s comparator is now **GIAB-grade**: beyond detecting *whether* a variant is present
(position + ref + alt), it scores **genotype concordance** — whether a het was called het and a hom
called hom. It also **decomposes multi-allelic records** (`A,C` with `GT=1/2` → two biallelic records)
instead of dropping them, and parses `GT` from the FORMAT/sample columns. A het called as hom-alt is now
counted as a *detection* true positive **and** a *genotype error* — the distinction the old detection-only
comparator could not see, and the metric GIAB tools (hap.py, vcfeval) report.

## Measured (synthetic diploid truth, the caller in isolation)

Known het/hom SNVs injected into a reference, reads sampled from both haplotypes, the BAM built directly
(so this isolates the **caller**, not alignment); deterministic (fixed-seed LCG). Numbers from
`cargo test --test germline_accuracy -- --nocapture`:

| Regime | Detection P / R / F1 | **Genotype concordance** | concordant / discordant |
|---|---|---|---|
| **40× / 0.5% error (clean), PASS** | 1.0000 / 1.0000 / 1.0000 | **1.0000** | 40 / 0 |
| 12× / 1.5% error (stress), PASS | 0.6792 / 0.9000 / 0.7742 | **0.9722** | 35 / 1 |
| 12× / 1.5% error (stress), all | 0.2899 / 1.0000 / 0.4494 | **0.9750** | 39 / 1 |

On clean data every recovered SNV is **genotyped correctly** (40/40). Under stress (12×, 1.5% error) a
single site's zygosity is mis-assigned (35/36 concordant among PASS TPs) — the honest graded behavior
of the diploid likelihood as evidence thins.

## Reproducibility

The harness is deterministic: identical seeds → identical truth, reads, calls, and therefore identical
metrics. In a CLI run, the call set carries a canonical-JSON BLAKE3 receipt (`--manifest`), so an
accuracy card and the run that produced it are **content-addressed together** — the artifact pairing no
other caller publishes.

## Honest scope

Synthetic (not yet GIAB), SNV-focused, zygosity (not phase) concordance. The comparator is now the
GIAB-grade interface; the real **HG002** number is the drop-in next step (this environment lacks the
`samtools`/aligner tooling to run it here):

```sh
rosalind eval-germline \
  --reference GRCh38.fa \
  --calls sample.vcf \
  --truth HG002_GRCh38_1_22_v4.2.1_benchmark.vcf \
  --regions HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

It prints `precision`, `recall`, `f1`, and the new `genotype_concordant` / `genotype_discordant` /
`genotype_unknown` / `genotype_concordance` lines. Folding **MAPQ** into the likelihood (relevant for
real externally-mapped BAMs, where MAPQ varies) is the deferred Sprint-2.2 half.
