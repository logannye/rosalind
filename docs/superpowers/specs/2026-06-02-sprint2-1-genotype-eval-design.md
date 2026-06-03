# Sprint 2.1 — Genotype-Aware Evaluation (design)

**Status:** Approved design — 2026-06-02. Increment 2.1 of the engineering roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md), Sprint 2). The GIAB credibility unlock — the *comparator*
upgrade that makes `eval-germline` GIAB-grade; the real HG002 run plugs into it.

---

## 1. Problem

The accuracy story is the biggest remaining objection ("I trust the contract, but can I trust the
*caller*?"). `eval-germline` exists, but the comparator is below the GIAB bar in two concrete ways:

1. **Detection-only.** `compare_callsets` (`eval/compare.rs`) matches variants by locus+allele
   (`BTreeSet<NormalizedVariant>`) and reports precision/recall/F1 — but **never checks genotype**. A
   het called as hom-alt counts as a full true positive. GIAB-grade tools (hap.py, vcfeval) score
   *genotype concordance*; a het-as-hom is a genotype error.
2. **Multi-allelic records are dropped.** `read_vcf_variants` (`eval/vcf.rs:69`) does
   `if alt_field.contains(',') { continue; }` — so every multi-allelic truth/call site is silently
   skipped, which on a real GIAB VCF discards a meaningful fraction of variants.

This increment closes both, plus the GT parsing they require, and emits a genotype-aware accuracy card
paired with the memory receipt. It is **fully testable here** against synthetic diploid truth (the
`germline_accuracy` harness already injects known het/hom SNVs). Rosalind's own germline VCF already
emits `GT` (`0/1`/`1/1` — `io/vcf.rs:23`), so its genotype concordance is directly scorable.

## 2. Goals / non-goals

**Goals**
- Parse `GT` from the VCF FORMAT/sample columns into the eval types.
- **Decompose** multi-allelic records into biallelic ones (per-ALT) with the correct decomposed
  genotype — stop dropping them.
- Score **genotype concordance** alongside detection P/R/F1 in `compare_callsets` /
  `ComparisonReport` / the `eval-germline` CLI output.
- Emit a **genotype-aware accuracy card** (`docs/findings/…`) from the synthetic harness, paired with
  the BLAKE3 receipt — the artifact no other caller publishes.

**Non-goals (deferred)**
- Folding **MAPQ** into the germline likelihood — 2.2 (behavior-changing; pairs with the real GIAB run,
  since our own aligner emits a constant MAPQ).
- Running the **real GIAB HG002** benchmark here (needs a downloaded BAM + samtools; this env lacks the
  tooling). It plugs into the upgraded comparator via the existing
  `eval-germline --reference … --calls … --truth … --regions highconf.bed` — documented, not run.
- Phasing (`0|1` vs `0/1` is treated as the same zygosity; we score zygosity concordance, not phase).
- Any change to the *caller* — this increment touches only the evaluator + tests + the card.

## 3. Design

### 3.1 Zygosity (`eval`, new)

```rust
/// Diploid zygosity of a (biallelic) variant call, parsed from `GT`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity { Het, Hom }
```

A small enum local to `genomics::eval` (no coupling to `call::types::Genotype`). Phase-insensitive:
`0/1`, `1/0`, `0|1` → `Het`; `1/1`, `1|1` → `Hom`. Missing/`.`/`./.` GT → `None` (scored as
genotype-*unknown*, never an error — real truth VCFs may lack per-sample GT).

### 3.2 `VcfVariant` gains an optional genotype

Add `pub genotype: Option<Zygosity>` to `VcfVariant`. Parsed from the FORMAT (`cols[8]`) + first sample
(`cols[9]`) columns when present: find the `GT` sub-field by its index in the colon-delimited FORMAT
key list, read it from the sample. Absent FORMAT/sample → `None` (back-compat: an 8-column VCF still
parses).

### 3.3 Multi-allelic decomposition (`read_vcf_variants`)

Replace the `contains(',') { continue }` skip with decomposition: for `ALT = a1,a2,…`, emit one
`VcfVariant` per ALT allele `ai` (1-based index `i`). The decomposed genotype counts how many times `i`
appears in the GT allele indices:

| GT (indices) | ALT | record for `a1` | record for `a2` |
|---|---|---|---|
| `0/1` | `a1` | `Het` | — |
| `1/1` | `a1` | `Hom` | — |
| `1/2` | `a1,a2` | `Het` | `Het` |
| `2/2` | `a1,a2` | (index 1 count 0 → **dropped**) | `Hom` |

So a record for allele `i` is emitted iff `i` appears in GT; `Het` if it appears once, `Hom` if twice.
When there is no GT, all comma-ALTs are emitted as separate biallelic records with `genotype = None`
(decomposed for detection, unscored for genotype). This matches the `bcftools norm -m-` convention
closely enough for SNV/indel scoring (the v1 scope).

### 3.4 Genotype-aware comparison (`compare_callsets` / `ComparisonReport`)

Switch the two `BTreeSet<NormalizedVariant>` to **`BTreeMap<NormalizedVariant, Option<Zygosity>>`** —
the key is the normalized locus+allele *identity* (unchanged: detection still matches by key), the
value is the decomposed zygosity. `NormalizedVariant` is **not** modified (genotype is *not* part of
identity, so a het-vs-hom at the same locus still *matches* for detection and is then compared).

`ComparisonReport` gains, scored over the detection true-positives:

```rust
pub genotype_concordant: usize,   // TP where both GTs known and equal
pub genotype_discordant: usize,   // TP where both GTs known and differ (a genotype error)
pub genotype_unknown: usize,      // TP where either GT is absent (unscored)
```

plus `pub fn genotype_concordance(&self) -> f64` = `concordant / (concordant + discordant)` (0.0 when
the denominator is 0). Detection precision/recall/F1 are unchanged. When the same key appears in both
maps (a TP), compare the call zygosity to the truth zygosity and increment the right counter.

### 3.5 CLI output + the accuracy card

- `run_eval` / `eval-germline` (`main.rs`) prints the existing detection line plus a genotype line, e.g.
  `genotype: concordance 0.98 (concordant 49, discordant 1, unknown 0 of 50 TP)`.
- A committed findings card: `docs/findings/2026-06-02-genotype-aware-eval.md`, generated by the
  synthetic harness — detection P/R/F1 + genotype concordance on injected het/hom truth — explicitly
  paired with the run's BLAKE3 receipt, and stating the real-HG002 command as the next step.

## 4. File-by-file change list

| File | Change |
|---|---|
| `src/genomics/eval/vcf.rs` | `Zygosity` enum; `VcfVariant.genotype`; GT parsing from FORMAT/sample; multi-allelic decomposition (replace the skip). |
| `src/genomics/eval/compare.rs` | `BTreeMap<NormalizedVariant, Option<Zygosity>>`; genotype-concordance counters + `genotype_concordance()` on `ComparisonReport`. |
| `src/genomics/eval/mod.rs` | Re-export `Zygosity`. |
| `src/main.rs` | `run_eval`: print the genotype-concordance line. |
| `tests/germline_accuracy.rs` (or a new `tests/genotype_eval.rs`) | Synthetic het+hom truth → assert genotype concordance; het-called-as-hom = detection-TP + genotype-discordant; a multi-allelic `1/2` truth decomposes to two biallelic records; a golden left-align/decomposition case. |
| `docs/findings/2026-06-02-genotype-aware-eval.md` (new) | The accuracy card + the real-GIAB command. |
| `README.md` (Accuracy section) | Note genotype-concordance is now scored; the GIAB command is the drop-in. |

## 5. Testing

1. **GT parsing** — `0/1`/`0|1`→`Het`, `1/1`→`Hom`, `./.`/absent→`None`; an 8-column (no-sample) VCF
   still parses (back-compat).
2. **Multi-allelic decomposition** — `ALT=A,C` with `GT=1/2` → two biallelic records (A:Het, C:Het);
   `GT=2/2` → only the C record (Hom); detection count reflects both alleles.
3. **Genotype concordance** — synthetic truth with known het + hom SNVs; Rosalind's calls scored:
   assert detection P/R and genotype concordance; then a deliberately mis-genotyped call (het truth,
   hom call) is a **detection TP but genotype-discordant** (the new capability the old comparator
   missed).
4. **Determinism / regression** — the existing `germline_accuracy` detection numbers are unchanged
   (genotype scoring is additive); `eval-germline` exit codes unchanged.
5. **The card** — generated deterministically; its numbers match the test assertions.

## 6. References

- `src/genomics/eval/{vcf,compare,normalize}.rs` — the comparator being upgraded.
- `src/io/vcf.rs:23` — Rosalind emits `GT` (`0/1`/`1/1`), so its genotype is scorable.
- `tests/germline_accuracy.rs` — the synthetic diploid-truth harness this extends.
- `docs/ROADMAP.md` Sprint 2 (FB-1) — this increment; MAPQ folding is FB-1's deferred half (2.2).
