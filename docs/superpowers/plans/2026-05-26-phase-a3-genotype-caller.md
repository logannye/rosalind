# Phase A3 — Calibrated, abstention-aware genotype caller Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `src/call/` — a real diploid biallelic-SNV genotype-likelihood model that turns the `PileupColumn` stream into calibrated germline calls (GT/GQ/PL/AD/DP/QUAL/FILTER) with honest abstention, plus a re-homed tumor/normal somatic LLR — replacing the legacy heuristic that discarded its posterior.

**Architecture:** A new crate-root `call` module (peer of `core` and `pileup`), built on `crate::core` + `crate::pileup` + `std` only. Pure, column-in/decision-out functions: `call::germline(&PileupColumn, &GermlineParams) -> Option<GermlineCall>` and `call::somatic(&PileupColumn, &PileupColumn, &SomaticParams) -> Option<SomaticCall>`. No VCF writing, no CLI wiring, no two-engine co-walk (those are A4/A5). No `genomics` coupling. The somatic likelihood math (`ln_factorial`/`ln_choose`/`ln_binom_pmf`/`somatic_llr`) is re-homed **verbatim** from `genomics/somatic/model.rs`; it is fed corrected columns, so the reverse-complement bug never enters.

**Tech Stack:** Rust, `std` f64 math (no new deps), TDD with `cargo test`.

**Design reference:** `docs/superpowers/specs/2026-05-26-phase-a-unified-pileup-genotype-design.md` §5 (genotype model), §3.3 (call module), §10 (tests).

---

## The genotype model (from spec §5, with the numerics pinned)

Reference allele `R = allele_index(column.ref_base)`; alt allele `A` = most-supported non-reference allele (biallelic; ties broken by lowest index). Per observation with base-quality `q`, error `ε = min(0.75, 10^(−q/10))` (capped at 0.75 so `1−ε > 0` even at `q=0` — no base is worse than a uniform random draw):

```
P(b | X)      = 1 − ε   if b == X   else   ε / 3
log L(0/0)   += ln P(b | R)
log L(0/1)   += ln( ½·P(b | R) + ½·P(b | A) )
log L(1/1)   += ln P(b | A)
```

Prior (default heterozygosity θ = 1e-3): `P(0/0)=1−1.5θ`, `P(0/1)=θ`, `P(1/1)=θ/2`. Log-posterior `= log L + ln prior`.

- **GT** = argmax log-posterior (ties → lower index, i.e. conservative). If GT is `0/0`, **return `None`** (abstain — emit no record).
- **PL** = `round(−10·log10 L(g))` re-normalized so the most-likely-by-likelihood genotype is 0, each capped at 255. Order `[0/0, 0/1, 1/1]` (VCF genotype order).
- **GQ** = second-smallest PL, capped at 99.
- **QUAL** = `−10·log10 P(0/0 | data)` (posterior, via stable log-sum-exp), floored at 0.
- **AD** = `[ref_count, alt_count]`; **DP** = callable depth (`column.depth()`).
- **FILTER** = `LowDepth` if `DP < min_depth`; else `LowQual` if `QUAL < min_qual`; else `PASS`.

**Note on honest abstention (principled deviation from spec §10's literal example):** Spec §10 lists "(2,1) → LowDepth". Under the calibrated θ=1e-3 model, a single alt read among 2–4 reads cannot overcome the hom-ref prior, so the model **abstains (`None`)** there — which spec §3.3 explicitly sanctions as a valid "no-call" outcome and which is the stronger, more honest form of "not a confident call." The `LowDepth` *flag* is demonstrated on a shallow-but-unambiguous variant (e.g. `(0,3)`), where a variant is genuinely called but flagged as too shallow to trust. Both behaviors are tested.

---

## File Structure

- `src/call/mod.rs` — module banner + re-exports (the public surface).
- `src/call/types.rs` — `Genotype`, `Filter`, `GermlineParams`, `GermlineCall`, `SomaticParams`, `SomaticCall`, `pub(crate) const ACGT`.
- `src/call/germline.rs` — `site_likelihoods` (private) + `call_germline` (public) + germline tests.
- `src/call/somatic.rs` — re-homed `ln_factorial`/`ln_choose`/`ln_binom_pmf`/`somatic_llr` (private) + `call_somatic` (public) + somatic tests.
- `src/lib.rs` — register `pub mod call;`.

Each file has one responsibility: `types` is data only, `germline` and `somatic` are the two decision functions. Tests live beside the code they cover.

---

### Task 1: `call` module scaffold + types

**Files:**
- Create: `src/call/types.rs`
- Create: `src/call/mod.rs`
- Modify: `src/lib.rs` (add `pub mod call;` next to `pub mod pileup;`)
- Test: in `src/call/types.rs`

- [ ] **Step 1: Write the failing test.** Create `src/call/types.rs` with the test at the bottom (types not yet defined → fails to compile):

```rust
//! Result and parameter types for the calling layer — pure data, no logic.

/// Diploid genotype over the biallelic {reference, alt} pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    /// Homozygous reference (0/0).
    HomRef,
    /// Heterozygous (0/1).
    Het,
    /// Homozygous alternate (1/1).
    HomAlt,
}

/// VCF FILTER status for an emitted call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    /// Passed all thresholds.
    Pass,
    /// QUAL below the configured minimum.
    LowQual,
    /// Depth below the configured minimum.
    LowDepth,
}

/// ACGT index (0..=3) → uppercase ASCII base. Shared by the calling functions.
pub(crate) const ACGT: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Tunable thresholds for germline genotyping and honest filtering.
#[derive(Debug, Clone)]
pub struct GermlineParams {
    /// Prior heterozygosity θ (default 1e-3).
    pub heterozygosity: f64,
    /// Calls with DP below this are flagged `LowDepth` (default 8).
    pub min_depth: u32,
    /// Calls with QUAL below this are flagged `LowQual` (default 30.0).
    pub min_qual: f64,
}

impl Default for GermlineParams {
    fn default() -> Self {
        Self {
            heterozygosity: 1e-3,
            min_depth: 8,
            min_qual: 30.0,
        }
    }
}

/// A germline variant call. Emitted only when GT ≠ HomRef. Locus-free — the
/// consumer pairs it with the originating column's `Locus`.
#[derive(Debug, Clone, PartialEq)]
pub struct GermlineCall {
    /// Called genotype.
    pub genotype: Genotype,
    /// Alternate base (uppercase ASCII).
    pub alt_base: u8,
    /// Phred site-is-variant quality, −10·log10 P(0/0 | data).
    pub qual: f64,
    /// Genotype quality (second-smallest PL, capped 99).
    pub gq: u8,
    /// Phred genotype likelihoods, ordered [0/0, 0/1, 1/1], min 0, cap 255.
    pub pl: [u32; 3],
    /// Allelic depths [ref, alt].
    pub ad: [u32; 2],
    /// Callable read depth.
    pub dp: u32,
    /// FILTER status.
    pub filter: Filter,
}

/// Tunable thresholds for tumor/normal somatic SNV calling.
#[derive(Debug, Clone)]
pub struct SomaticParams {
    /// Minimum tumor depth to consider a site (default 8).
    pub min_tumor_depth: u32,
    /// Minimum normal depth to consider a site (default 8).
    pub min_normal_depth: u32,
    /// Minimum tumor alt allele fraction (default 0.05).
    pub min_tumor_af: f32,
    /// Maximum tolerated normal alt allele fraction (default 0.01).
    pub max_normal_af: f32,
    /// Minimum LLR-derived quality to emit (default 10.0).
    pub min_quality: f64,
    /// Sequencing error rate for the noise model (default 1e-3).
    pub seq_error_rate: f64,
}

impl Default for SomaticParams {
    fn default() -> Self {
        Self {
            min_tumor_depth: 8,
            min_normal_depth: 8,
            min_tumor_af: 0.05,
            max_normal_af: 0.01,
            min_quality: 10.0,
            seq_error_rate: 1e-3,
        }
    }
}

/// A tumor/normal somatic SNV call. Locus-free.
#[derive(Debug, Clone, PartialEq)]
pub struct SomaticCall {
    /// Reference base (uppercase ASCII).
    pub ref_base: u8,
    /// Alternate base (uppercase ASCII).
    pub alt_base: u8,
    /// Tumor alt count / depth.
    pub tumor_alt: u32,
    /// Tumor depth.
    pub tumor_depth: u32,
    /// Normal alt count.
    pub normal_alt: u32,
    /// Normal depth.
    pub normal_depth: u32,
    /// Tumor alt allele fraction.
    pub tumor_af: f32,
    /// Normal alt allele fraction.
    pub normal_af: f32,
    /// LLR-derived Phred-scaled quality (clamped 0..200).
    pub quality: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_spec() {
        let g = GermlineParams::default();
        assert_eq!(g.heterozygosity, 1e-3);
        assert_eq!(g.min_depth, 8);
        assert_eq!(g.min_qual, 30.0);

        let s = SomaticParams::default();
        assert_eq!(s.min_tumor_depth, 8);
        assert_eq!(s.min_tumor_af, 0.05);
        assert_eq!(s.max_normal_af, 0.01);
        assert_eq!(s.seq_error_rate, 1e-3);

        assert_eq!(ACGT[0], b'A');
        assert_eq!(ACGT[3], b'T');
        assert_ne!(Genotype::Het, Genotype::HomAlt);
        assert_ne!(Filter::Pass, Filter::LowDepth);
    }
}
```

- [ ] **Step 2: Create the module file.** Create `src/call/mod.rs`:

```rust
//! The calling layer: turn the `PileupColumn` stream into calibrated,
//! abstention-aware variant calls. Built on `crate::core` + `crate::pileup`
//! only; no VCF writing or CLI wiring (those are later phases).

pub mod germline;
pub mod somatic;
pub mod types;

pub use germline::call_germline;
pub use somatic::call_somatic;
pub use types::{
    Filter, GermlineCall, GermlineParams, Genotype, SomaticCall, SomaticParams,
};
```

- [ ] **Step 3: Create placeholder submodules so `mod.rs` compiles.** Create `src/call/germline.rs`:

```rust
//! The germline diploid genotype-likelihood model (filled in Task 2–3).
```

Create `src/call/somatic.rs`:

```rust
//! The tumor/normal somatic LLR caller (filled in Task 4).
```

- [ ] **Step 4: Register the module.** In `src/lib.rs`, add `pub mod call;` immediately after the existing `pub mod pileup;` line (keep modules grouped; alphabetical among the new crate-root modules `call`, `core`, `pileup` is fine — place it before `core` or after `pileup`, matching the file's existing ordering style).

- [ ] **Step 5: Run the test.** Run: `cargo test call::types`
Expected: PASS (`defaults_match_spec`). Also run `cargo build 2>&1 | grep -i warning` — expect no `src/call/` warnings (the empty `germline`/`somatic` banners are fine; `call_germline`/`call_somatic` re-exports resolve once Tasks 2–4 land — if the `pub use` lines error because the functions don't exist yet, temporarily comment the two `pub use germline::...`/`somatic::...` lines and restore them in Task 3/4; note this in your commit).

- [ ] **Step 6: Commit.**

```bash
git add src/call/mod.rs src/call/types.rs src/call/germline.rs src/call/somatic.rs src/lib.rs
git commit -m "feat(call): scaffold call module + germline/somatic result & param types" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Germline genotype log-likelihoods

**Files:**
- Modify: `src/call/germline.rs`
- Test: in `src/call/germline.rs`

- [ ] **Step 1: Write the failing test.** Replace the contents of `src/call/germline.rs` with the imports, a test helper, and a test that exercises `site_likelihoods` (not yet defined → fails to compile):

```rust
//! The germline diploid genotype-likelihood model: per-observation base-quality
//! error integrated into [0/0, 0/1, 1/1] log-likelihoods, then a prior, then a
//! calibrated, abstention-aware call.

use crate::call::types::{Filter, GermlineCall, GermlineParams, Genotype, ACGT};
use crate::core::allele_index;
use crate::pileup::PileupColumn;

const LN10: f64 = std::f64::consts::LN_10;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Locus, Position};
    use crate::pileup::Obs;

    /// Build a column at chr0:100 with `ref_base` and observations given as
    /// `(allele_index, base_qual)` pairs (all forward strand, mapq 60).
    fn col(ref_base: u8, obs_spec: &[(u8, u8)]) -> PileupColumn {
        let obs: Vec<Obs> = obs_spec
            .iter()
            .map(|&(allele, base_qual)| Obs {
                allele,
                base_qual,
                mapq: 60,
                reverse: false,
            })
            .collect();
        PileupColumn {
            locus: Locus {
                contig: 0,
                pos: Position(100),
            },
            ref_base,
            raw_depth: obs.len() as u32,
            obs,
        }
    }

    #[test]
    fn likelihoods_rank_genotypes_correctly() {
        // 15 ref (A) + 15 alt (C), all bq30 → het is the most-likely genotype.
        let het: Vec<(u8, u8)> = (0..15)
            .map(|_| (0u8, 30u8))
            .chain((0..15).map(|_| (1u8, 30u8)))
            .collect();
        let sl = site_likelihoods(&col(b'A', &het)).unwrap();
        assert_eq!(sl.alt_idx, 1); // C
        assert_eq!(sl.ad, [15, 15]);
        assert_eq!(sl.dp, 30);
        // het (index 1) has the largest log-likelihood
        assert!(sl.log_l[1] > sl.log_l[0]);
        assert!(sl.log_l[1] > sl.log_l[2]);

        // 20 alt only → hom-alt dominates.
        let homalt: Vec<(u8, u8)> = (0..20).map(|_| (1u8, 30u8)).collect();
        let sl = site_likelihoods(&col(b'A', &homalt)).unwrap();
        assert!(sl.log_l[2] > sl.log_l[1]);
        assert!(sl.log_l[2] > sl.log_l[0]);
    }

    #[test]
    fn no_alt_or_non_callable_ref_yields_none() {
        // All reference → nothing to call.
        let allref: Vec<(u8, u8)> = (0..20).map(|_| (0u8, 30u8)).collect();
        assert!(site_likelihoods(&col(b'A', &allref)).is_none());
        // Non-callable reference base.
        assert!(site_likelihoods(&col(b'N', &[(1, 30), (1, 30)])).is_none());
    }
}
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test call::germline`
Expected: FAIL — does not compile (`site_likelihoods` undefined).

- [ ] **Step 3: Implement `site_likelihoods`.** Add above the `#[cfg(test)]` block in `src/call/germline.rs`:

```rust
/// Per-genotype log-likelihoods plus the biallelic context for a column.
struct SiteLikelihoods {
    /// Natural-log likelihoods, ordered [0/0, 0/1, 1/1].
    log_l: [f64; 3],
    /// Index (0..=3) of the alt allele.
    alt_idx: usize,
    /// Allelic depths [ref, alt].
    ad: [u32; 2],
    /// Callable depth.
    dp: u32,
}

/// Accumulate diploid genotype log-likelihoods from per-observation base
/// qualities. Returns `None` when the reference base is non-callable or no
/// non-reference allele is observed (nothing to call).
fn site_likelihoods(column: &PileupColumn) -> Option<SiteLikelihoods> {
    let ref_idx = allele_index(column.ref_base)?;
    let counts = column.allele_counts();

    // Alt = most-supported non-reference allele; ties → lowest index.
    let mut alt_idx: Option<usize> = None;
    let mut best = 0u32;
    for i in 0..4 {
        if i == ref_idx {
            continue;
        }
        if counts[i] > best {
            best = counts[i];
            alt_idx = Some(i);
        }
    }
    let alt_idx = alt_idx?;
    if best == 0 {
        return None;
    }

    let mut log_l = [0.0f64; 3];
    for o in &column.obs {
        // ε from Phred base quality, capped at 0.75 so 1−ε stays positive (and
        // the log finite) even at q=0 — no base is worse than a random draw.
        let eps = (10f64.powf(-(o.base_qual as f64) / 10.0)).min(0.75);
        let a = o.allele as usize;
        let p_ref = if a == ref_idx { 1.0 - eps } else { eps / 3.0 };
        let p_alt = if a == alt_idx { 1.0 - eps } else { eps / 3.0 };
        log_l[0] += p_ref.ln();
        log_l[1] += (0.5 * p_ref + 0.5 * p_alt).ln();
        log_l[2] += p_alt.ln();
    }

    Some(SiteLikelihoods {
        log_l,
        alt_idx,
        ad: [counts[ref_idx], counts[alt_idx]],
        dp: column.depth(),
    })
}
```

- [ ] **Step 4: Run the test to verify it passes.** Run: `cargo test call::germline`
Expected: PASS (2 tests). Also `cargo build 2>&1 | grep -i warning` — `call_germline`/`LN10`/`ACGT`/`Filter`/`GermlineCall`/`GermlineParams`/`Genotype` may be reported unused until Task 3; that is expected mid-task and resolved in Task 3. If `mod.rs`'s `pub use germline::call_germline` was commented out in Task 1, leave it commented until Task 3.

- [ ] **Step 5: Commit.**

```bash
git add src/call/germline.rs
git commit -m "feat(call): germline genotype log-likelihoods (per-base-quality ε model)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Germline call — posterior, PL/GT/GQ/QUAL, FILTER, abstention

**Files:**
- Modify: `src/call/germline.rs`
- Modify: `src/call/mod.rs` (restore the `pub use germline::call_germline;` if it was commented in Task 1)
- Test: in `src/call/germline.rs`

- [ ] **Step 1: Write the failing tests.** Add these tests inside the existing `#[cfg(test)] mod tests` block in `src/call/germline.rs` (reuse the `col` helper):

```rust
    #[test]
    fn hom_ref_site_abstains() {
        // No alt observed → no variant record.
        let allref: Vec<(u8, u8)> = (0..20).map(|_| (0u8, 30u8)).collect();
        assert!(call_germline(&col(b'A', &allref), &GermlineParams::default()).is_none());
    }

    #[test]
    fn thin_alt_evidence_abstains() {
        // 3 ref + 1 alt at bq30: a lone alt read cannot overcome the θ=1e-3
        // prior, so the calibrated model abstains (no-call) rather than overcall.
        let thin = [(0u8, 30u8), (0, 30), (0, 30), (1, 30)];
        assert!(call_germline(&col(b'A', &thin), &GermlineParams::default()).is_none());
    }

    #[test]
    fn clear_het_passes() {
        let het: Vec<(u8, u8)> = (0..15)
            .map(|_| (0u8, 30u8))
            .chain((0..15).map(|_| (1u8, 30u8)))
            .collect();
        let call = call_germline(&col(b'A', &het), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::Het);
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.ad, [15, 15]);
        assert_eq!(call.dp, 30);
        assert_eq!(call.pl[1], 0); // het is the most-likely-by-likelihood genotype
        assert!(call.pl[0] > 0 && call.pl[2] > 0);
        assert_eq!(call.filter, Filter::Pass);
        assert!(call.gq > 0);
    }

    #[test]
    fn clear_hom_alt_passes() {
        let homalt: Vec<(u8, u8)> = (0..20).map(|_| (1u8, 30u8)).collect();
        let call = call_germline(&col(b'A', &homalt), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::HomAlt);
        assert_eq!(call.pl[2], 0);
        assert_eq!(call.filter, Filter::Pass);
    }

    #[test]
    fn shallow_variant_is_flagged_low_depth_not_dropped() {
        // 3 clean alt reads, 0 ref: unambiguously hom-alt, but DP=3 < 8 → emitted
        // with the LowDepth flag (a real variant, too shallow to fully trust).
        let shallow: Vec<(u8, u8)> = (0..3).map(|_| (1u8, 30u8)).collect();
        let call = call_germline(&col(b'A', &shallow), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::HomAlt);
        assert_eq!(call.dp, 3);
        assert_eq!(call.filter, Filter::LowDepth);
    }

    #[test]
    fn qual_increases_with_evidence() {
        let mk = |n: usize| -> f64 {
            let obs: Vec<(u8, u8)> = (0..n)
                .map(|_| (0u8, 30u8))
                .chain((0..n).map(|_| (1u8, 30u8)))
                .collect();
            call_germline(&col(b'A', &obs), &GermlineParams::default())
                .unwrap()
                .qual
        };
        // Same 0.5 alt fraction, deeper site → higher site-is-variant QUAL.
        assert!(mk(20) > mk(8));
    }
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test call::germline`
Expected: FAIL — does not compile (`call_germline` undefined).

- [ ] **Step 3: Implement `call_germline` and `argmax3`.** Add above the `#[cfg(test)]` block in `src/call/germline.rs`:

```rust
/// Call a germline genotype from a pileup column. Returns `None` (abstains) when
/// the most-probable genotype is homozygous reference — honest by default: thin
/// or absent evidence is a no-call, not a confident reference assertion turned
/// variant.
pub fn call_germline(column: &PileupColumn, params: &GermlineParams) -> Option<GermlineCall> {
    let sl = site_likelihoods(column)?;

    let theta = params.heterozygosity;
    let log_prior = [
        (1.0 - 1.5 * theta).ln(),
        theta.ln(),
        (theta / 2.0).ln(),
    ];
    let log_post = [
        sl.log_l[0] + log_prior[0],
        sl.log_l[1] + log_prior[1],
        sl.log_l[2] + log_prior[2],
    ];

    // GT = argmax posterior; hom-ref → abstain.
    let gt_idx = argmax3(&log_post);
    if gt_idx == 0 {
        return None;
    }

    // PL: −10·log10 L(g), re-normalized so the max-likelihood genotype is 0,
    // each capped at 255. Order [0/0, 0/1, 1/1].
    let max_log_l = sl.log_l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0u32; 3];
    for g in 0..3 {
        let phred = -10.0 * (sl.log_l[g] - max_log_l) / LN10;
        pl[g] = phred.round().min(255.0) as u32;
    }

    // GQ = second-smallest PL, capped 99.
    let mut sorted = pl;
    sorted.sort_unstable();
    let gq = sorted[1].min(99) as u8;

    // QUAL = −10·log10 P(0/0 | data), via a numerically stable log-sum-exp.
    let max_lp = log_post.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let log_z = max_lp
        + log_post
            .iter()
            .map(|lp| (lp - max_lp).exp())
            .sum::<f64>()
            .ln();
    let log_p_homref = log_post[0] - log_z; // ≤ 0
    let qual = (-10.0 * log_p_homref / LN10).max(0.0);

    let genotype = if gt_idx == 1 {
        Genotype::Het
    } else {
        Genotype::HomAlt
    };
    let filter = if sl.dp < params.min_depth {
        Filter::LowDepth
    } else if qual < params.min_qual {
        Filter::LowQual
    } else {
        Filter::Pass
    };

    Some(GermlineCall {
        genotype,
        alt_base: ACGT[sl.alt_idx],
        qual,
        gq,
        pl,
        ad: sl.ad,
        dp: sl.dp,
        filter,
    })
}

/// Index of the maximum of three values; ties resolve to the lowest index.
fn argmax3(v: &[f64; 3]) -> usize {
    let mut best = 0;
    for i in 1..3 {
        if v[i] > v[best] {
            best = i;
        }
    }
    best
}
```

- [ ] **Step 4: Restore the re-export.** If Task 1 Step 5 commented out `pub use germline::call_germline;` in `src/call/mod.rs`, uncomment it now.

- [ ] **Step 5: Run the tests to verify they pass.** Run: `cargo test call::germline`
Expected: PASS (all germline tests — the 2 from Task 2 plus the 6 added here). Then `cargo test` (full suite green), `cargo build 2>&1 | grep -i warning` (no `src/call/` warnings), `cargo clippy --lib 2>&1 | grep -A2 'src/call'` (report any lints — fix obvious ones like needless casts).

- [ ] **Step 6: Commit.**

```bash
git add src/call/germline.rs src/call/mod.rs
git commit -m "feat(call): calibrated abstention-aware germline call (GT/GQ/PL/QUAL/FILTER)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Somatic LLR re-home

**Files:**
- Modify: `src/call/somatic.rs`
- Modify: `src/call/mod.rs` (restore `pub use somatic::call_somatic;` if commented in Task 1)
- Test: in `src/call/somatic.rs`

- [ ] **Step 1: Write the failing tests.** Replace the contents of `src/call/somatic.rs` with imports, the test helper, and tests (the functions don't exist yet → fails to compile):

```rust
//! The tumor/normal somatic SNV caller. The likelihood math
//! (`ln_factorial`/`ln_choose`/`ln_binom_pmf`/`somatic_llr`) is re-homed
//! verbatim from the legacy `genomics::somatic::model`; here it is fed corrected
//! `PileupColumn`s, so the legacy reverse-complement bug never enters.

use crate::call::types::{SomaticCall, SomaticParams, ACGT};
use crate::core::allele_index;
use crate::pileup::PileupColumn;

const LN10: f64 = std::f64::consts::LN_10;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Locus, Position};
    use crate::pileup::Obs;

    fn col(ref_base: u8, obs_spec: &[(u8, u8)]) -> PileupColumn {
        let obs: Vec<Obs> = obs_spec
            .iter()
            .map(|&(allele, base_qual)| Obs {
                allele,
                base_qual,
                mapq: 60,
                reverse: false,
            })
            .collect();
        PileupColumn {
            locus: Locus {
                contig: 0,
                pos: Position(100),
            },
            ref_base,
            raw_depth: obs.len() as u32,
            obs,
        }
    }

    /// `n` observations of a single allele index, all bq30.
    fn reads(allele: u8, n: usize) -> Vec<(u8, u8)> {
        (0..n).map(|_| (allele, 30u8)).collect()
    }

    #[test]
    fn filters_normal_contamination() {
        // Tumor 15A+5C (alt C, AF 0.25); normal 18A+2C (AF 0.10 > max 0.01) → None.
        let mut t = reads(0, 15);
        t.extend(reads(1, 5));
        let mut n = reads(0, 18);
        n.extend(reads(1, 2));
        let call = call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default());
        assert!(call.is_none());
    }

    #[test]
    fn calls_clean_somatic_snv() {
        // Tumor 10A+10C (AF 0.5); normal 20A (AF 0) → somatic C call.
        let mut t = reads(0, 10);
        t.extend(reads(1, 10));
        let n = reads(0, 20);
        let call =
            call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default()).unwrap();
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.tumor_alt, 10);
        assert_eq!(call.normal_alt, 0);
        assert!((call.tumor_af - 0.5).abs() < 1e-6);
        assert!(call.quality >= SomaticParams::default().min_quality);
    }

    #[test]
    fn respects_min_depth() {
        // Below min_tumor_depth → None even with a clean signal.
        let t = reads(1, 4); // 4 alt reads, depth 4 < 8
        let n = reads(0, 20);
        assert!(call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default()).is_none());
    }

    #[test]
    fn llr_math_matches_legacy_reference() {
        // Pin the re-homed math against hand-computed expectations.
        assert_eq!(ln_choose(5, 0), 0.0);
        assert!((ln_choose(5, 2) - 10f64.ln()).abs() < 1e-9); // C(5,2)=10
        assert_eq!(ln_choose(2, 5), f64::NEG_INFINITY);
        // Strong tumor signal, clean normal → positive LLR (somatic favored).
        assert!(somatic_llr(10, 20, 0, 20, 1e-3) > 0.0);
        // No tumor signal → non-positive LLR (noise favored).
        assert!(somatic_llr(0, 20, 0, 20, 1e-3) <= 0.0);
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail.** Run: `cargo test call::somatic`
Expected: FAIL — does not compile (`call_somatic`, `ln_choose`, `somatic_llr` undefined).

- [ ] **Step 3: Implement the re-homed math + `call_somatic`.** Add above the `#[cfg(test)]` block in `src/call/somatic.rs`:

```rust
/// Call a somatic SNV from a tumor column and a position-matched normal column.
/// Filter order matches the legacy caller: depth → alt-exists → AF gates →
/// LLR/quality. Returns `None` when any gate fails. Locus-free.
pub fn call_somatic(
    tumor: &PileupColumn,
    normal: &PileupColumn,
    params: &SomaticParams,
) -> Option<SomaticCall> {
    let t_depth = tumor.depth();
    let n_depth = normal.depth();
    if t_depth < params.min_tumor_depth || n_depth < params.min_normal_depth {
        return None;
    }

    let ref_idx = allele_index(tumor.ref_base)?;
    let t_counts = tumor.allele_counts();
    let n_counts = normal.allele_counts();

    // Alt = tumor's most-supported non-reference allele; ties → lowest index.
    let mut alt_idx: Option<usize> = None;
    let mut best = 0u32;
    for i in 0..4 {
        if i == ref_idx {
            continue;
        }
        if t_counts[i] > best {
            best = t_counts[i];
            alt_idx = Some(i);
        }
    }
    let alt_idx = alt_idx?;
    let t_alt = t_counts[alt_idx];
    if t_alt == 0 {
        return None;
    }
    let n_alt = n_counts[alt_idx];

    let t_af = t_alt as f32 / t_depth as f32;
    let n_af = n_alt as f32 / n_depth as f32;
    if t_af < params.min_tumor_af || n_af > params.max_normal_af {
        return None;
    }

    let llr = somatic_llr(t_alt, t_depth, n_alt, n_depth, params.seq_error_rate);
    let quality = (llr / LN10 * 10.0).clamp(0.0, 200.0);
    if quality < params.min_quality {
        return None;
    }

    Some(SomaticCall {
        ref_base: tumor.ref_base.to_ascii_uppercase(),
        alt_base: ACGT[alt_idx],
        tumor_alt: t_alt,
        tumor_depth: t_depth,
        normal_alt: n_alt,
        normal_depth: n_depth,
        tumor_af: t_af,
        normal_af: n_af,
        quality,
    })
}

// --- Likelihood math, re-homed verbatim from genomics::somatic::model ---

fn ln_factorial(n: u32) -> f64 {
    // For sequencing depths, n is small; compute exactly and deterministically.
    let mut acc = 0.0f64;
    for k in 2..=n {
        acc += (k as f64).ln();
    }
    acc
}

fn ln_choose(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

fn ln_binom_pmf(k: u32, n: u32, p: f64) -> f64 {
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    ln_choose(n, k) + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln()
}

fn somatic_llr(t_alt: u32, t_depth: u32, n_alt: u32, n_depth: u32, err: f64) -> f64 {
    // Somatic vs noise-only: under somatic the tumor alt fraction is the MLE and
    // the normal alt fraction is ~err; under noise both are ~err.
    let t_p_hat = (t_alt as f64 / t_depth.max(1) as f64).clamp(err, 1.0 - 1e-6);
    let somatic = ln_binom_pmf(t_alt, t_depth, t_p_hat) + ln_binom_pmf(n_alt, n_depth, err);
    let noise = ln_binom_pmf(t_alt, t_depth, err) + ln_binom_pmf(n_alt, n_depth, err);
    somatic - noise
}
```

- [ ] **Step 4: Restore the re-export.** If Task 1 Step 5 commented out `pub use somatic::call_somatic;` in `src/call/mod.rs`, uncomment it now.

- [ ] **Step 5: Run the tests to verify they pass.** Run: `cargo test call::somatic`
Expected: PASS (4 somatic tests). Then `cargo test` (full suite green), `cargo build 2>&1 | grep -i warning` (no `src/call/` warnings), `cargo fmt --all -- --check` (clean; run `cargo fmt --all` if needed), `cargo clippy --lib 2>&1 | grep -A2 'src/call'` (report/fix obvious lints).

- [ ] **Step 6: Commit.**

```bash
git add src/call/somatic.rs src/call/mod.rs
git commit -m "feat(call): re-home tumor/normal somatic LLR onto corrected pileup columns" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage (§5, §3.3, §10):**
- §5 genotype model — Task 2 (likelihoods) + Task 3 (prior, PL, GT, GQ, QUAL, FILTER, abstention). ✔
- §3.3 `call::germline -> site decision (PASS / LowQual / LowDepth / no-call)` — Task 3 returns `Option<GermlineCall>` with `Filter`; `None` = no-call. ✔
- §3.3 `call::somatic -> Option<SomaticCall>` with tumor/normal AD/DP/AF — Task 4. ✔ (Two-engine co-walk + VCF are explicitly A4/A5, per spec scope.)
- §10 germline calibration tests (hom-ref/het/hom-alt; abstention; QUAL monotonic) — Task 3. The (2,1)→LowDepth example is implemented as the principled "(thin)→no-call" + "(shallow clear)→LowDepth" pair, documented above and sanctioned by §3.3's no-call option. ✔
- §10 somatic — Task 4 re-homes the normal-contamination filter test + adds a positive call + depth gate + math pin. ✔

**Placeholder scan:** none — every step has complete code and exact commands. The one cross-task dependency (the `pub use` re-exports in `mod.rs` referencing functions added in Tasks 3–4) is handled explicitly with a comment-out/restore note.

**Type consistency:** `GermlineCall { genotype, alt_base:u8, qual:f64, gq:u8, pl:[u32;3], ad:[u32;2], dp:u32, filter }`; `SomaticCall { ref_base, alt_base, tumor_alt, tumor_depth, normal_alt, normal_depth, tumor_af:f32, normal_af:f32, quality:f64 }`; `call_germline(&PileupColumn, &GermlineParams) -> Option<GermlineCall>`; `call_somatic(&PileupColumn, &PileupColumn, &SomaticParams) -> Option<SomaticCall>`. Substrate API used: `PileupColumn { ref_base, obs, depth(), allele_counts() }`, `Obs { allele, base_qual }`, `core::allele_index(u8) -> Option<usize>`, `core::{Locus, Position}` (tests). All match the shipped A1/A2 code (verified against `src/core/sequence.rs`, `src/core/locus.rs`, `src/pileup/column.rs`).

**Scope:** germline GL model + somatic decision only; no VCF, no CLI, no co-walk, no `genomics` deletion (all A4/A5). Focused and self-contained.
