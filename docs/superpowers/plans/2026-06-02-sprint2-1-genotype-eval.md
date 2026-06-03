# Sprint 2.1 — Genotype-Aware Evaluation: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Run INLINE.

**Goal:** Upgrade `eval-germline`'s comparator to GIAB-grade — parse `GT`, decompose multi-allelic records, and score genotype concordance alongside detection P/R/F1 — validated on synthetic diploid truth.

**Architecture:** Add a `Zygosity` enum + optional genotype to `VcfVariant`; decompose multi-allelic ALTs in the parser; switch `compare_callsets` to a `BTreeMap<NormalizedVariant, Option<Zygosity>>` so genotype is scored over detection true-positives without changing the detection key. The synthetic `germline_accuracy` harness (which already knows het/hom truth) asserts the concordance end-to-end.

**Tech Stack:** Rust, the existing `src/genomics/eval/` module + `tests/germline_accuracy.rs`.

**Spec:** `docs/superpowers/specs/2026-06-02-sprint2-1-genotype-eval-design.md`

---

### Task 1: `Zygosity` + GT parsing + multi-allelic decomposition (`eval/vcf.rs`)

**Files:** `src/genomics/eval/vcf.rs`, `src/genomics/eval/mod.rs`

- [ ] **Step 1: Write the failing unit tests** (append to `src/genomics/eval/vcf.rs`'s `#[cfg(test)] mod tests`, or add one)

```rust
#[cfg(test)]
mod gt_tests {
    use super::*;

    #[test]
    fn parses_gt_zygosity_from_format_sample() {
        let vcf = "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:30\n\
                   chr1\t20\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\n\
                   chr1\t30\t.\tT\tA\t.\tPASS\t.\tGT\t0|1\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v[0].genotype, Some(Zygosity::Het));
        assert_eq!(v[1].genotype, Some(Zygosity::Hom));
        assert_eq!(v[2].genotype, Some(Zygosity::Het)); // phase-insensitive
    }

    #[test]
    fn missing_sample_is_genotype_none_and_still_parses() {
        // An 8-column VCF (no FORMAT/sample) still parses, genotype unknown.
        let vcf = "chr1\t10\t.\tA\tC\t.\tPASS\t.\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].genotype, None);
    }

    #[test]
    fn multiallelic_decomposes_per_alt_with_decomposed_genotype() {
        // ALT A,C with GT 1/2 -> two biallelic records, each Het.
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\tGT\t1/2\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].alternate, b"A");
        assert_eq!(v[0].genotype, Some(Zygosity::Het));
        assert_eq!(v[1].alternate, b"C");
        assert_eq!(v[1].genotype, Some(Zygosity::Het));
    }

    #[test]
    fn multiallelic_hom_second_allele_drops_the_absent_first() {
        // ALT A,C with GT 2/2 -> only the C record, Hom.
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\tGT\t2/2\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].alternate, b"C");
        assert_eq!(v[0].genotype, Some(Zygosity::Hom));
    }

    #[test]
    fn multiallelic_without_gt_emits_all_alleles_unscored() {
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].genotype, None);
        assert_eq!(v[1].genotype, None);
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib gt_tests 2>&1 | grep -E "error|test result"`
Expected: FAIL — `Zygosity`, `VcfVariant.genotype`, and decomposition do not exist.

- [ ] **Step 3: Add `Zygosity` + `genotype` field + GT parse + decomposition**

In `src/genomics/eval/vcf.rs`, add the enum near the top (after the imports):

```rust
/// Diploid zygosity of a biallelic call, parsed from a VCF `GT`. Phase-insensitive
/// (`0|1` and `0/1` are both `Het`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    /// One alt copy (`0/1`).
    Het,
    /// Two alt copies (`1/1`).
    Hom,
}
```

Add the field to `VcfVariant` (after `info`):

```rust
    /// Diploid zygosity from the sample `GT`, if a FORMAT/sample column is present.
    pub genotype: Option<Zygosity>,
```

Replace the body of `read_vcf_variants`'s per-line ALT handling. Replace this:

```rust
        let alt_field = cols[4];
        // Single-allele only for v1 harness.
        if alt_field.contains(',') {
            continue;
        }
        let alternate = alt_field.as_bytes().to_vec();
        let qual = if cols[5] == "." {
            None
        } else {
            Some(cols[5].parse().map_err(|_| VcfParseError::InvalidLine {
                line: line_no,
                msg: "QUAL is not a float".to_string(),
            })?)
        };
        let filter = cols[6].to_string();
        let info = parse_info(cols[7]);

        out.push(VcfVariant {
            chrom,
            pos0,
            reference,
            alternate,
            qual,
            filter,
            info,
        });
```

with:

```rust
        let qual = if cols[5] == "." {
            None
        } else {
            Some(cols[5].parse().map_err(|_| VcfParseError::InvalidLine {
                line: line_no,
                msg: "QUAL is not a float".to_string(),
            })?)
        };
        let filter = cols[6].to_string();
        let info = parse_info(cols[7]);

        // GT allele indices from the FORMAT (cols[8]) + first sample (cols[9]), if present.
        let gt_indices = parse_gt_indices(cols.get(8).copied(), cols.get(9).copied());

        // Decompose multi-allelic ALT into one biallelic record per ALT allele.
        for (alt_idx, alt) in cols[4].split(',').enumerate() {
            let allele_index = (alt_idx + 1) as u8; // GT indices: ref=0, first alt=1, ...
            let genotype = match &gt_indices {
                // No GT column: emit the allele for detection, unscored for genotype.
                None => None,
                // GT present: zygosity = how many times this allele index appears in it.
                Some(idx) => match idx.iter().filter(|&&a| a == allele_index).count() {
                    0 => continue, // allele absent from the genotype -> drop this record
                    1 => Some(Zygosity::Het),
                    _ => Some(Zygosity::Hom),
                },
            };
            out.push(VcfVariant {
                chrom: chrom.clone(),
                pos0,
                reference: reference.clone(),
                alternate: alt.as_bytes().to_vec(),
                qual,
                filter: filter.clone(),
                info: info.clone(),
                genotype,
            });
        }
```

Add the GT parser helper (near `parse_info`):

```rust
/// Parse the GT allele indices (e.g. `0/1` -> [0,1], `1|2` -> [1,2]) from the
/// FORMAT + first sample columns. Returns `None` if no GT field is present.
fn parse_gt_indices(format: Option<&str>, sample: Option<&str>) -> Option<Vec<u8>> {
    let (format, sample) = (format?, sample?);
    let gt_pos = format.split(':').position(|k| k == "GT")?;
    let gt = sample.split(':').nth(gt_pos)?;
    let indices: Vec<u8> = gt
        .split(['/', '|'])
        .filter_map(|a| a.parse::<u8>().ok()) // '.' (no-call) is filtered out
        .collect();
    if indices.is_empty() {
        None
    } else {
        Some(indices)
    }
}
```

> Note: `reference`/`chrom`/`filter`/`info` are now cloned per ALT (cheap; multi-allelic sites are
> rare). The single-allele case (`enumerate` yields one item) is unchanged in behavior except for the
> added `genotype`.

- [ ] **Step 4: Re-export `Zygosity`**

In `src/genomics/eval/mod.rs`, add `Zygosity` to the `pub use vcf::{…}` line.

- [ ] **Step 5: Run to verify the GT tests pass**

Run: `cargo test --lib gt_tests 2>&1 | grep -E "test result"` → PASS.
(Other crates/tests that construct `VcfVariant` now fail to compile — fixed in Task 3. That is expected
at this checkpoint; do not commit yet.)

- [ ] **Step 6: Hold the commit until Task 3** (adding the field breaks the 2 harness constructors; commit Tasks 1–3 together once green).

---

### Task 2: Genotype-concordance scoring (`eval/compare.rs`)

**Files:** `src/genomics/eval/compare.rs`

- [ ] **Step 1: Write the failing unit test** (append to `compare.rs` `#[cfg(test)] mod tests`, or add)

```rust
#[cfg(test)]
mod gt_concordance_tests {
    use super::*;
    use crate::genomics::eval::Zygosity;
    use std::collections::BTreeMap;

    fn vv(pos1: u32, refb: &str, alt: &str, gt: Option<Zygosity>) -> VcfVariant {
        VcfVariant {
            chrom: "chr1".to_string(),
            pos0: pos1 - 1,
            reference: refb.as_bytes().to_vec(),
            alternate: alt.as_bytes().to_vec(),
            qual: None,
            filter: ".".to_string(),
            info: BTreeMap::new(),
            genotype: gt,
        }
    }

    #[test]
    fn het_called_as_hom_is_a_detection_tp_but_a_genotype_error() {
        let refs = BTreeMap::from([("chr1".to_string(), b"ACGTACGTAC".to_vec())]);
        let truth = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Hom))]; // right site, wrong genotype
        let r = compare_callsets(&refs, &calls, &truth, None).unwrap();
        assert_eq!(r.true_positive, 1, "detection still matches");
        assert_eq!(r.false_positive, 0);
        assert_eq!(r.false_negative, 0);
        assert_eq!(r.genotype_concordant, 0);
        assert_eq!(r.genotype_discordant, 1);
        assert_eq!(r.genotype_concordance(), 0.0);
    }

    #[test]
    fn matching_genotype_is_concordant() {
        let refs = BTreeMap::from([("chr1".to_string(), b"ACGTACGTAC".to_vec())]);
        let truth = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let r = compare_callsets(&refs, &calls, &truth, None).unwrap();
        assert_eq!((r.genotype_concordant, r.genotype_discordant), (1, 0));
        assert_eq!(r.genotype_concordance(), 1.0);
    }

    #[test]
    fn unknown_genotype_is_counted_separately_not_as_error() {
        let refs = BTreeMap::from([("chr1".to_string(), b"ACGTACGTAC".to_vec())]);
        let truth = vec![vv(5, "A", "C", None)]; // truth without GT
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let r = compare_callsets(&refs, &calls, &truth, None).unwrap();
        assert_eq!(r.true_positive, 1);
        assert_eq!((r.genotype_concordant, r.genotype_discordant, r.genotype_unknown), (0, 0, 1));
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib gt_concordance 2>&1 | grep -E "error|test result"`
Expected: FAIL — the report has no genotype fields and the comparator uses sets, not maps.

- [ ] **Step 3: Add the genotype fields + map-based scoring**

In `ComparisonReport`, add (after `by_type`):

```rust
    /// True positives whose genotype matched truth (both genotypes known).
    pub genotype_concordant: usize,
    /// True positives whose genotype disagreed with truth (both known) — a genotype error.
    pub genotype_discordant: usize,
    /// True positives where either genotype was absent (unscored for concordance).
    pub genotype_unknown: usize,
```

Add the method to the `impl ComparisonReport`:

```rust
    /// Genotype concordance = concordant / (concordant + discordant), over the scorable
    /// true positives. `0.0` when none are scorable.
    pub fn genotype_concordance(&self) -> f64 {
        let scorable = self.genotype_concordant + self.genotype_discordant;
        if scorable == 0 {
            0.0
        } else {
            self.genotype_concordant as f64 / scorable as f64
        }
    }
```

Replace the two `BTreeSet<NormalizedVariant>` with maps and the scoring loop. Change:

```rust
    let mut calls_set: BTreeSet<NormalizedVariant> = BTreeSet::new();
    let mut truth_set: BTreeSet<NormalizedVariant> = BTreeSet::new();

    for v in calls {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        calls_set.insert(normalize_variant(reference_for(references, v)?, v)?);
    }
    for v in truth {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        truth_set.insert(normalize_variant(reference_for(references, v)?, v)?);
    }
```

to:

```rust
    // Key = normalized locus+allele (the detection identity); value = decomposed zygosity.
    let mut calls_map: BTreeMap<NormalizedVariant, Option<Zygosity>> = BTreeMap::new();
    let mut truth_map: BTreeMap<NormalizedVariant, Option<Zygosity>> = BTreeMap::new();

    for v in calls {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        calls_map.insert(normalize_variant(reference_for(references, v)?, v)?, v.genotype);
    }
    for v in truth {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        truth_map.insert(normalize_variant(reference_for(references, v)?, v)?, v.genotype);
    }
```

Then update the scoring loops to use the maps and accumulate genotype counters. Replace the
`for v in calls_set.iter()` / `for v in truth_set.iter()` blocks and the `Ok(ComparisonReport { … })`
with:

```rust
    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut fn_ = 0usize;
    let mut gt_concordant = 0usize;
    let mut gt_discordant = 0usize;
    let mut gt_unknown = 0usize;

    let mut by_type: BTreeMap<VariantType, (usize, usize, usize)> = BTreeMap::new();

    for (v, call_gt) in calls_map.iter() {
        let vt = variant_type(v);
        if let Some(truth_gt) = truth_map.get(v) {
            tp += 1;
            by_type.entry(vt).or_insert((0, 0, 0)).0 += 1;
            match (call_gt, truth_gt) {
                (Some(c), Some(t)) if c == t => gt_concordant += 1,
                (Some(_), Some(_)) => gt_discordant += 1,
                _ => gt_unknown += 1,
            }
        } else {
            fp += 1;
            by_type.entry(vt).or_insert((0, 0, 0)).1 += 1;
        }
    }

    for v in truth_map.keys() {
        if !calls_map.contains_key(v) {
            fn_ += 1;
            by_type.entry(variant_type(v)).or_insert((0, 0, 0)).2 += 1;
        }
    }

    Ok(ComparisonReport {
        total_truth: truth_map.len(),
        total_calls: calls_map.len(),
        true_positive: tp,
        false_positive: fp,
        false_negative: fn_,
        by_type,
        genotype_concordant: gt_concordant,
        genotype_discordant: gt_discordant,
        genotype_unknown: gt_unknown,
    })
```

Update the imports at the top of `compare.rs`: `use std::collections::BTreeMap;` (drop `BTreeSet`),
and add `Zygosity` to the `use super::{…}` line.

- [ ] **Step 4: Run to verify the concordance tests pass**

Run: `cargo test --lib gt_concordance 2>&1 | grep -E "test result"` → PASS (the lib compiles; `tests/`
binaries still fail until Task 3).

- [ ] **Step 5: Hold the commit until Task 3.**

---

### Task 3: Wire-up — CLI, harness, card, README; commit Tasks 1–3

**Files:** `src/main.rs`, `tests/germline_accuracy.rs`, `docs/findings/2026-06-02-genotype-aware-eval.md` (new), `README.md`

- [ ] **Step 1: Fix the 2 harness `VcfVariant` constructions + assert concordance** (`tests/germline_accuracy.rs`)

The calls construction (~`:223`) — add `genotype` from the call's genotype. Add to the struct literal
(after `info: BTreeMap::new(),`):

```rust
            genotype: Some(match call.genotype {
                Genotype::Het => Zygosity::Het,
                _ => Zygosity::Hom, // HomAlt; HomRef is never emitted as a variant call
            }),
```

(Add `Genotype` and `Zygosity` to the test's imports — `use rosalind::...::{Genotype, Zygosity}`; the
existing `GermlineCall` import path shows where.)

The truth construction (~`:243`) — add `genotype` from the truth het flag (after `info: BTreeMap::new(),`):

```rust
            genotype: Some(if v.het { Zygosity::Het } else { Zygosity::Hom }),
```

After the existing detection assertions in the clean-baseline test, add a genotype-concordance assertion
(use the `report_pass` or `report_all` already computed; on the clean 40× baseline every recovered SNV
should be genotype-concordant):

```rust
    // Genotype concordance: on the clean baseline, recovered SNVs match zygosity.
    assert!(
        o.pass.genotype_concordant > 0,
        "expected some genotype-concordant TPs"
    );
    assert_eq!(
        o.pass.genotype_discordant, 0,
        "clean baseline should have no genotype errors: {:?}",
        o.pass
    );
```

> Adapt the field access (`o.pass` / `report_pass`) to the harness's actual `AccuracyOutcome` shape —
> the point is to assert `genotype_discordant == 0` and `genotype_concordant > 0` on the clean run.

- [ ] **Step 2: Add the genotype line to `run_eval`** (`src/main.rs`)

After the `println!("f1={f1:.6}");` line, add:

```rust
    println!("genotype_concordant={}", report.genotype_concordant);
    println!("genotype_discordant={}", report.genotype_discordant);
    println!("genotype_unknown={}", report.genotype_unknown);
    println!("genotype_concordance={:.6}", report.genotype_concordance());
```

- [ ] **Step 3: Run the full suite green**

Run: `cargo test 2>&1 | grep -iE "FAILED" || echo "no failures"` → `no failures`.
Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -2` → clean (the new gate must stay green).
Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"; cargo fmt && cargo fmt --check && echo "fmt clean"`.

- [ ] **Step 4: Capture the accuracy numbers + write the card**

Run the accuracy test with output to read the real numbers:
`cargo test --test germline_accuracy -- --nocapture 2>&1 | grep -iE "accuracy|genotype|precision|recall"`

Write `docs/findings/2026-06-02-genotype-aware-eval.md` recording: the harness setup (synthetic diploid
truth, het+hom), the measured detection P/R/F1 AND genotype concordance from the run above, the explicit
pairing with the run's reproducibility receipt (the determinism guarantee), and the one-command
real-GIAB next step:

```
rosalind eval-germline --reference GRCh38.fa --calls sample.vcf \
  --truth HG002_GRCh38_1_22_v4.2.1_benchmark.vcf --regions HG002_..._highconf.bed
```

State the honest scope (synthetic, SNV-focused; the real HG002 number is the next step) — mirror the
framing in the existing `docs/findings/2026-06-02-germline-accuracy.md`.

- [ ] **Step 5: Update the README Accuracy section** (`README.md`)

Add one sentence: the comparator now scores **genotype concordance** (not just detection), and the
GIAB-grade `eval-germline` is the drop-in for a real HG002 run. Point at the new findings card.

- [ ] **Step 6: Commit Tasks 1–3 together**

```bash
git add src/genomics/eval/ src/main.rs tests/germline_accuracy.rs docs/findings/2026-06-02-genotype-aware-eval.md README.md
git commit -m "feat(eval): genotype-aware comparison — GT parsing, multi-allelic decomposition, concordance scoring"
```

---

## Final verification

- [ ] `cargo test` green (incl. the new GT + concordance + harness assertions).
- [ ] `cargo clippy --all-targets -- -D warnings` clean; `cargo build --release` 0 warnings; `cargo fmt --check` clean.
- [ ] `eval-germline` on a small synthetic pair prints the `genotype_concordance=…` line; a het-vs-hom mismatch shows up as `genotype_discordant`.
- [ ] The findings card's numbers match the test output and name the real-GIAB command.
- [ ] CI green on GitHub.
