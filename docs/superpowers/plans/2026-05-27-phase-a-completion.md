# Phase A Completion Plan — somatic on the new engine + legacy retirement

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Finish Phase A: (A6) route the *somatic* CLI through the new engine (SNV co-walk → `call::somatic` → `io::vcf` + BLAKE3 receipt), then (A7) delete the entire legacy calling path the spec §8 named for removal and migrate its tests onto the new path. Somatic **indels are deferred to Phase C** (the legacy indel path is reverse-strand-buggy; Phase C rebuilds indels — germline + somatic — correctly on the new engine).

**Architecture:** A6 adds a tumor/normal co-walk to `call::pipeline` (two `PileupEngine`s, position-joined, `call::somatic` per shared position) and rewrites `run_somatic`. A7 deletes the legacy callers + string VCF writers + the legacy somatic model, updates `genomics` exports + `main.rs` imports, and migrates the four legacy-coupled tests to the new path (dropping the deferred-indel test). `PileupProcessor`/`PileupSummary`/`PileupWorkload` + `CompressedEvaluator` are **kept** (the spec's plugin substrate; used by `plugin/*` + `python_bindings`). `pileup_stream::BamPileupStream` becomes off-calling-path legacy — **kept** for now (flag for Phase-E cleanup).

**Tech Stack:** Rust; the A1–A5 modules; `rust_htslib` (already used in `io/bam`).

**Branch:** `rosalind/phase-a-completion` (off `main` @ the merged Phase-A head). One PR at the end.

---

## Stage A6 — Somatic SNV calling on the new engine

### Task 1: `call::pipeline::call_somatic_region` (tumor/normal co-walk)

**Files:** Modify `src/call/pipeline.rs` (add the function + a test). 

- [ ] **Step 1: Write the failing test.** Add to the `#[cfg(test)] mod tests` in `src/call/pipeline.rs` (reuse the existing `read` helper; add `use crate::call::SomaticParams;` and `use crate::pileup::SliceSource;` if not already imported):

```rust
    #[test]
    fn calls_a_somatic_snv_from_tumor_normal_cowalk() {
        // Reference AAAA. At position 1: tumor has C (alt) in 6/12 reads; normal is
        // all A. Expect one somatic SNV at pos 1, alt C.
        let reference: Arc<[u8]> = Arc::from(b"AAAA".to_vec().into_boxed_slice());
        let tumor: Vec<AlignedRead> = (0..6)
            .map(|_| read(0, b"ACAA", false))
            .chain((0..6).map(|_| read(0, b"AAAA", false)))
            .collect();
        let normal: Vec<AlignedRead> = (0..12).map(|_| read(0, b"AAAA", false)).collect();
        let params = SomaticParams {
            min_tumor_depth: 4,
            min_normal_depth: 4,
            min_tumor_af: 0.1,
            max_normal_af: 0.05,
            min_quality: 0.0,
            seq_error_rate: 1e-3,
        };
        let calls = call_somatic_region(
            SliceSource::new(tumor),
            SliceSource::new(normal),
            reference,
            0,
            0..4,
            PileupParams::default(),
            &params,
        )
        .unwrap();
        assert_eq!(calls.len(), 1);
        let (locus, call) = &calls[0];
        assert_eq!(locus.pos, Position(1));
        assert_eq!(call.ref_base, b'A');
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.normal_alt, 0);
    }
```

- [ ] **Step 2: Run it to confirm it fails.** `cargo test call::pipeline` → FAIL (`call_somatic_region` undefined).

- [ ] **Step 3: Implement the co-walk.** Add to `src/call/pipeline.rs` (extend the `use` line with `call_somatic, SomaticCall, SomaticParams` from `crate::call`, and `PileupColumn` from `crate::pileup`):

```rust
/// Call somatic SNVs by co-walking a tumor and a normal pileup over `region` of
/// `contig`. A call is attempted only at positions covered in BOTH samples
/// (a somatic call needs a normal baseline). Returns `(locus, call)` per emitted
/// site. Collects both column streams first (bounded by region; a bounded-memory
/// streaming co-walk is Phase-D work alongside `.csi` fetch).
pub fn call_somatic_region<T: ReadSource, N: ReadSource>(
    tumor: T,
    normal: N,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    somatic_params: &SomaticParams,
) -> Result<Vec<(Locus, SomaticCall)>, CoreError> {
    let tumor_cols: Vec<PileupColumn> =
        PileupEngine::new(tumor, Arc::clone(&reference), contig, region.clone(), pileup_params.clone())
            .collect::<Result<_, _>>()?;
    let normal_cols: Vec<PileupColumn> =
        PileupEngine::new(normal, reference, contig, region, pileup_params).collect::<Result<_, _>>()?;

    // Merge-join by position (both streams are ascending in `pos`).
    let mut out = Vec::new();
    let mut n = 0usize;
    for t in &tumor_cols {
        while n < normal_cols.len() && normal_cols[n].locus.pos < t.locus.pos {
            n += 1;
        }
        if n < normal_cols.len() && normal_cols[n].locus.pos == t.locus.pos {
            if let Some(call) = call_somatic(t, &normal_cols[n], somatic_params) {
                out.push((t.locus, call));
            }
        }
    }
    Ok(out)
}
```

- [ ] **Step 4: Run it to confirm it passes.** `cargo test call::pipeline` → PASS (the prior pipeline tests + this one). Then `cargo build 2>&1 | grep -i warning` (none new), `cargo fmt --all -- --check`.

- [ ] **Step 5: Commit.**
```bash
git add src/call/pipeline.rs
git commit -m "feat(call/pipeline): somatic SNV tumor/normal co-walk over the pileup engine" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 2: Rewire `run_somatic` onto the new vertical

**Files:** Modify `src/main.rs`.

- [ ] **Step 1: Replace the somatic calling + output section of `run_somatic`.** Keep everything up to and including the BAM sort (`sort_bam_deterministic` of tumor/normal). Replace the block from `// Call somatic SNVs and indels.` through the end of the manifest writing with:

```rust
    // Call somatic SNVs on the new engine (indels deferred to Phase C).
    use rosalind::call::{call_somatic_region, SomaticParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::write_somatic_vcf;
    use rosalind::pileup::PileupParams;
    use rosalind::provenance::{blake3_file, write_manifest, FileHash, RunManifest};

    let start_call = Instant::now();
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(fasta.name.clone(), reference.len() as u32);
    let region = 0..(reference.len() as u32);

    let tumor_src = BamSource::new(&tumor_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read tumor BAM {}: {e}", tumor_sorted.display()))?;
    let normal_src = BamSource::new(&normal_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read normal BAM {}: {e}", normal_sorted.display()))?;
    let calls = call_somatic_region(
        tumor_src,
        normal_src,
        Arc::clone(&reference),
        contig_id,
        region,
        PileupParams::default(),
        &SomaticParams::default(),
    )
    .map_err(|e| anyhow!("somatic calling failed: {e}"))?;
    let dur_call = start_call.elapsed();

    // Write spec-valid somatic VCF (TUMOR/NORMAL).
    {
        let file = File::create(&output_vcf)
            .with_context(|| format!("failed to create somatic VCF {}", output_vcf.display()))?;
        let mut writer = io::BufWriter::new(file);
        write_somatic_vcf(&mut writer, &contigs, &calls)?;
        writer.flush()?;
    }

    // Reproducibility receipt (BLAKE3, canonical JSON).
    let tumor_inputs: Vec<PathBuf> = match (&tumor_fastq, &tumor_r1, &tumor_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let normal_inputs: Vec<PathBuf> = match (&normal_fastq, &normal_r1, &normal_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let mut manifest = RunManifest::new("somatic");
    manifest.inputs.push(FileHash {
        path: reference_path.display().to_string(),
        blake3: blake3_file(&reference_path)?,
    });
    for p in tumor_inputs.iter().chain(normal_inputs.iter()) {
        if p.exists() {
            manifest.inputs.push(FileHash {
                path: p.display().to_string(),
                blake3: blake3_file(p)?,
            });
        }
    }
    manifest.outputs.push(FileHash {
        path: output_vcf.display().to_string(),
        blake3: blake3_file(&output_vcf)?,
    });
    manifest.params.insert("somatic_snv_only".to_string(), "true".to_string());
    let manifest_path = write_manifest(&output_vcf, &manifest)?;
    eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
```

Adapt the surrounding timing/logging lines (`dur_align_tumor`, `dur_sort`, `dur_call`, `start_total`) to whatever the function already prints — keep the existing summary log working; just drop any references to `indels`/`snvs.len()+indels.len()` (report `calls.len()` instead).

- [ ] **Step 2: Drop the now-unused legacy somatic imports + helpers used only here.** In `run_somatic`, `SomaticCaller`/`SomaticCallerConfig` and `write_combined_somatic_vcf`/`write_somatic_manifest` are no longer called. Leave the `fn write_combined_somatic_vcf` / `fn write_somatic_manifest` definitions for Task 4 to delete (they may still be referenced until then — if the compiler now flags them `dead_code`, add a temporary `#[allow(dead_code)]` with a `// removed in A7` note, OR delete them now if nothing else uses them). Do NOT remove the `SomaticCaller` import from the top-of-file `use` yet (Task 4 handles imports after the somatic module is deleted) unless the compiler errors — if it only *warns*, leave it; if a name is now genuinely unused and warns, you may remove just that name.

- [ ] **Step 3: Build + smoke test.** `cargo build 2>&1 | grep -iE 'error|warning'`. Then replicate a somatic run if toy somatic data exists, or construct a minimal one: confirm `run_somatic` produces a `##fileformat=VCFv4.2` somatic VCF with `TUMOR` + `NORMAL` columns and a `.manifest.json` sidecar. Report commands + output. `cargo test` — full suite still green (the legacy somatic tests still pass; they're migrated in A7).

- [ ] **Step 4: Commit.**
```bash
git add src/main.rs
git commit -m "feat(cli): route 'somatic' through the new engine (SNV co-walk + VCF + receipt)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Stage A7 — Retire the legacy calling path + migrate tests

**Files:** delete `src/genomics/{statistics.rs, variant_caller.rs, vcf.rs}` + `src/genomics/somatic/` (whole dir); modify `src/genomics/mod.rs`, `src/main.rs`; migrate `tests/{variant_pipeline,determinism,golden_vcf,truthset_validation}.rs`; delete `tests/somatic_indels.rs`.

### Task 3: Delete legacy germline + migrate germline tests

- [ ] **Step 1: Migrate the germline tests FIRST (so deletion doesn't leave them dangling).**

`tests/variant_pipeline.rs` — replace the whole file with a test that drives the NEW germline pipeline:
```rust
use std::sync::Arc;

use rosalind::call::{call_germline_region, GermlineParams, Genotype};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
use rosalind::pileup::{PileupParams, SliceSource};

fn read(pos: u32, seq: &[u8]) -> AlignedRead {
    AlignedRead {
        contig: 0,
        pos: Position(pos),
        mapq: 60,
        flags: SamFlags::default(),
        cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
        seq: Arc::from(seq.to_vec().into_boxed_slice()),
        qual: Arc::from(vec![35u8; seq.len()].into_boxed_slice()),
    }
}

#[test]
fn germline_pipeline_detects_a_het_snv() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGT".to_vec().into_boxed_slice());
    // Position 3 (ref T): half the reads carry A.
    let reads = vec![
        read(0, b"ACGTACGT"),
        read(0, b"ACGAACGT"),
        read(0, b"ACGTACGT"),
        read(0, b"ACGAACGT"),
    ];
    let sites = call_germline_region(
        SliceSource::new(reads),
        reference,
        0,
        0..8,
        PileupParams::default(),
        &GermlineParams::default(),
    )
    .unwrap();
    let site = sites.iter().find(|(l, _, _)| l.pos == Position(3)).expect("variant at pos 3");
    assert_eq!(site.1, b'T'); // ref base
    assert_eq!(site.2.alt_base, b'A');
    assert_eq!(site.2.genotype, Genotype::Het);
}
```

`tests/determinism.rs` — replace with a determinism test over the new path + writer:
```rust
use std::collections::HashSet;
use std::sync::Arc;

use blake3::hash;
use rosalind::call::{call_germline_region, GermlineParams};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, ContigSet, Position, SamFlags};
use rosalind::io::vcf::{render_germline_vcf, GermlineRow};
use rosalind::pileup::{PileupParams, SliceSource};

fn read(pos: u32, seq: &[u8]) -> AlignedRead {
    AlignedRead {
        contig: 0,
        pos: Position(pos),
        mapq: 60,
        flags: SamFlags::default(),
        cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
        seq: Arc::from(seq.to_vec().into_boxed_slice()),
        qual: Arc::from(vec![35u8; seq.len()].into_boxed_slice()),
    }
}

#[test]
fn germline_calling_and_vcf_are_deterministic() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGTACGTACGT".to_vec().into_boxed_slice());
    let reads = vec![read(0, b"ACGTACGT"), read(2, b"GTAATCGT"), read(0, b"ACAATCGT")];
    let mut contigs = ContigSet::new();
    contigs.push("chrDet", reference.len() as u32);

    let mut fingerprints = HashSet::new();
    for _ in 0..5 {
        let sites = call_germline_region(
            SliceSource::new(reads.clone()),
            Arc::clone(&reference),
            0,
            0..reference.len() as u32,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();
        let rows: Vec<GermlineRow> = sites
            .into_iter()
            .map(|(locus, ref_base, call)| GermlineRow { locus, ref_base, call })
            .collect();
        let vcf = render_germline_vcf(&contigs, "S", &rows).unwrap();
        fingerprints.insert(hash(vcf.as_bytes()));
    }
    assert_eq!(fingerprints.len(), 1, "germline VCF diverged across runs");
}
```

`tests/golden_vcf.rs` — replace with a golden snapshot of the NEW germline VCF, and **regenerate the snapshot file**:
```rust
#[path = "common/mod.rs"]
mod common;
use common::assert_snapshot;
use rosalind::call::{Filter, GermlineCall, Genotype};
use rosalind::core::{ContigSet, Locus, Position};
use rosalind::io::vcf::{render_germline_vcf, GermlineRow};

#[test]
fn germline_vcf_matches_golden() {
    let mut contigs = ContigSet::new();
    contigs.push("chr1", 100_000);
    let rows = vec![
        GermlineRow {
            locus: Locus { contig: 0, pos: Position(99) },
            ref_base: b'T',
            call: GermlineCall {
                genotype: Genotype::Het,
                alt_base: b'A',
                qual: 42.0,
                gq: 40,
                pl: [42, 0, 60],
                ad: [6, 6],
                dp: 12,
                filter: Filter::Pass,
            },
        },
        GermlineRow {
            locus: Locus { contig: 0, pos: Position(199) },
            ref_base: b'G',
            call: GermlineCall {
                genotype: Genotype::HomAlt,
                alt_base: b'C',
                qual: 88.0,
                gq: 60,
                pl: [120, 60, 0],
                ad: [0, 8],
                dp: 8,
                filter: Filter::Pass,
            },
        },
    ];
    let actual = render_germline_vcf(&contigs, "SAMPLE", &rows).unwrap();
    assert_snapshot("variants/simple.vcf", &actual);
}
```
Then **regenerate `tests/snapshots/variants/simple.vcf`**: read `tests/common/mod.rs` to learn how `assert_snapshot` updates (it likely supports an `UPDATE_SNAPSHOTS`/`UPDATE_EXPECT` env var, or writes on mismatch). If it has an update mode, run it (e.g. `UPDATE_SNAPSHOTS=1 cargo test golden_vcf`); otherwise compute the expected VCF by running the test once, capture `actual`, and write it to the snapshot file. The new snapshot must be the exact `render_germline_vcf` output for the two rows above. Verify `cargo test golden_vcf` passes against the committed snapshot.

- [ ] **Step 2: Delete the legacy germline source.** Delete `src/genomics/statistics.rs`, `src/genomics/variant_caller.rs`, `src/genomics/vcf.rs`. In `src/genomics/mod.rs` remove the `mod statistics;`, `mod variant_caller;`, `mod vcf;` declarations and their `pub use` lines (`bayesian_variant_caller, VariantCall`; `StreamingVariantCaller, Variant, VariantCallerError`; `render_vcf, write_vcf`).

- [ ] **Step 3: Fix `main.rs` fallout.** Remove any now-dangling references. The germline path already uses the new vertical (A5). Remove `StreamingVariantCaller`/`write_vcf`/`Variant`/`VariantCall` from imports if present (A5 may have already removed some). Build and fix only legacy-germline-related errors.

- [ ] **Step 4: Build + test germline migration.** `cargo build 2>&1 | grep -iE 'error|warning'` (none). `cargo test variant_pipeline determinism golden_vcf` (3 migrated tests pass). `cargo test` (full suite — somatic legacy tests still pass; A7 Task 4 handles those).

- [ ] **Step 5: Commit.**
```bash
git add -A
git commit -m "refactor(genomics): delete legacy germline caller + VCF writer (migrated to new vertical)" \
  -m "Removes bayesian_variant_caller, StreamingVariantCaller, and the legacy VCF string writer; their tests now drive call::pipeline + io::vcf. call_variants(Vec<AlignedRead>) is superseded by call::pipeline::call_germline_region." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 4: Delete legacy somatic + migrate/drop somatic tests

- [ ] **Step 1: Migrate the somatic SNV test.** In `tests/truthset_validation.rs`, replace the body of `somatic_snv_truth_mixture_simple` so it calls the NEW path instead of `SomaticCaller`. Keep the BAM construction (12 tumor-alt C + 8 tumor-ref A at pos 50; 20 normal-ref A; `sort_bam_deterministic`), then:
```rust
    use rosalind::call::{call_somatic_region, SomaticParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::pileup::PileupParams;

    let mut contigs = ContigSet::new();
    contigs.push("chr1", 1000);
    let calls = call_somatic_region(
        BamSource::new(&tumor_sorted, &contigs).unwrap(),
        BamSource::new(&normal_sorted, &contigs).unwrap(),
        Arc::clone(&reference),
        0,
        0..200,
        PileupParams::default(),
        &SomaticParams { min_tumor_depth: 10, min_normal_depth: 10, min_tumor_af: 0.2, max_normal_af: 0.01, min_quality: 0.0, seq_error_rate: 1e-3 },
    )
    .unwrap();
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].0.pos, Position(50));   // import core::Position
    assert_eq!(calls[0].1.ref_base, b'A');
    assert_eq!(calls[0].1.alt_base, b'C');
```
Update the imports at the top of the file (drop `SomaticCaller, SomaticCallerConfig`; add what the new path needs). The `alignment_truth_exact_matches_map_to_expected_position` test (uses `BWTAligner`) stays unchanged.

- [ ] **Step 2: Delete the deferred-indel test.** `git rm tests/somatic_indels.rs` (somatic indel calling is deferred to Phase C; document in the commit message).

- [ ] **Step 3: Delete the legacy somatic source.** Delete the whole `src/genomics/somatic/` directory (`model.rs`, `vcf.rs`, `mod.rs`). In `src/genomics/mod.rs` remove `mod somatic;` and its `pub use somatic::{...}` line.

- [ ] **Step 4: Fix `main.rs` fallout.** Remove `SomaticCaller, SomaticCallerConfig, SomaticIndel, SomaticVariant` from the `use rosalind::genomics::{…}` import. Delete the now-unused `fn write_combined_somatic_vcf` and `fn write_somatic_manifest` (and any `write_somatic_manifest`/text-manifest helpers used only by the old somatic path). Build and fix only legacy-somatic fallout. If `EvalSomatic`/`run_eval_somatic` references nothing deleted (it uses `compare_callsets`/`read_vcf_variants`/`BedIndex` — all in `genomics::eval`, retained), leave it untouched.

- [ ] **Step 5: Build + full verification.** `cargo build 2>&1 | grep -iE 'error|warning'` (none). `cargo test` — FULL suite green. `cargo clippy --lib 2>&1 | grep -iE 'src/(call|io|pileup|core|provenance)'` (no new lints in the new modules; pre-existing legacy lints OK). `cargo fmt --all -- --check`. `grep -rn "bayesian_variant_caller\|StreamingVariantCaller\|SomaticCaller\|render_vcf\|write_vcf\b" src/ tests/` — should be EMPTY except the new `io::vcf` `write_somatic_vcf`/`render_somatic_vcf` (which are the NEW ones — distinguish by path `src/io/vcf.rs`).

- [ ] **Step 6: Commit.**
```bash
git add -A
git commit -m "refactor(genomics): delete legacy somatic caller + VCF writer (SNV path on new engine; indels deferred to Phase C)" \
  -m "Removes SomaticCaller/SomaticVariant/SomaticIndel + the legacy somatic VCF writer. Somatic SNV calling now runs through call::somatic; somatic indel calling is deferred to Phase C (the legacy indel path was reverse-strand-buggy). Drops tests/somatic_indels.rs accordingly." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before PR)
- `cargo test` — full suite green (report totals).
- `cargo build` — 0 warnings; `cargo clippy` — no new lints in `core/pileup/call/io/provenance`; `cargo fmt --all -- --check` clean.
- CLI smoke: replicate the germline e2e (still green) AND a somatic run (new VCF + manifest).
- Confirm retained: `PileupProcessor`/`PileupSummary`/`PileupWorkload` + `CompressedEvaluator` (plugin substrate) still compile + exported; `plugin/*` + `python_bindings` unaffected; `pileup_stream::BamPileupStream` + `tests/pileup_stream.rs` retained (off-path legacy, flagged for Phase E).

## Self-Review
- **Spec §8 coverage:** somatic LLR + co-walk → `call/` ✔ (A3 LLR + A6 co-walk); `bayesian_variant_caller` deleted ✔; legacy `StreamingVariantCaller` deleted ✔ (`call_variants` superseded by `call::pipeline` — full-retirement decision); `vcf.rs` + `somatic/vcf.rs` replaced by `io/vcf` ✔; `PileupProcessor`/`CompressedEvaluator` retained as plugin-facing ✔.
- **Deferred (documented):** somatic indels → Phase C; `pileup_stream::BamPileupStream` off-path legacy → Phase E; bounded-streaming co-walk/BamSource → Phase D.
- **No placeholders;** test migrations specify the new-API calls + assertions; deletions are explicit file lists.
- **Safety:** germline tests migrated before germline deletion; somatic test migrated before somatic deletion; full suite is the gate at each step.
