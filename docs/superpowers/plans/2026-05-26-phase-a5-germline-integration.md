# Phase A5 — Germline calling integration onto the new engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route the `rosalind variants` CLI through the new vertical — `BamSource`/`SliceSource` → `PileupEngine` → `call::germline` → `io::vcf` + a BLAKE3 provenance receipt — so the shipped command produces calibrated, spec-valid VCF from the bug-free engine (reverse-strand + CIGAR correct, MAPQ filter actually honored) instead of the legacy heuristic.

**Architecture:** Two new small pieces — `io::bam::BamSource` (the htslib→`core::AlignedRead` boundary, implementing `pileup::ReadSource`) and `call::pipeline::call_germline_region` (the engine→caller orchestration) — then a surgical rewrite of `main.rs::run_variants` to use them. The legacy `StreamingVariantCaller`, `bayesian_variant_caller`, `write_vcf`/`Variant`, and `SomaticCaller` are **left intact** (still compiled and tested) so every existing test stays green; this PR migrates the germline CLI path only. Full legacy deletion and somatic-CLI rewiring are explicit follow-ons (see "Deferred", below).

**Tech Stack:** Rust, `rust_htslib` (already a dep, stays inside `io/`), the A1–A4 modules (`core`, `pileup`, `call`, `io::vcf`, `provenance`).

**Design reference:** `docs/superpowers/specs/2026-05-26-phase-a-unified-pileup-genotype-design.md` §4 (data flow), §8 (migration). **Deferred from §8 (documented, not done here):** deleting `statistics::bayesian_variant_caller` + the legacy `StreamingVariantCaller`/`vcf.rs` (they back existing tests — `variant_pipeline.rs`, `golden_vcf.rs`; deleting requires migrating those tests, best done under review); somatic CLI rewiring onto `call::somatic` + `io::vcf` somatic (the somatic LLR + writer already exist & are tested as library code from A3/A4 — only the `run_somatic` wiring remains); a true bounded-memory streaming `BamSource` (this one pre-loads; bounded BAM streaming belongs with Phase-D `.csi` fetch + RSS gating).

**Consumes (shipped & green):** `core::{AlignedRead{contig:u32,pos:Position,mapq:u8,flags:SamFlags,cigar:Vec<CigarOp>,seq:Arc<[u8]>,qual:Arc<[u8]>}, CigarOp, CigarOpKind{Match,Insertion,Deletion,RefSkip,SoftClip,HardClip,Pad}, SamFlags(pub u16), ContigSet, Locus, Position, CoreError{MalformedRecord(String),InvalidContig,Io}}`; `pileup::{PileupEngine::new(source, reference:Arc<[u8]>, contig:u32, region:Range<u32>, params:PileupParams), PileupParams{min_mapq,min_base_qual,skip_secondary,skip_supplementary,skip_duplicate}, ReadSource, SliceSource}`; `call::{call_germline(&PileupColumn,&GermlineParams)->Option<GermlineCall>, GermlineCall, GermlineParams}`; `io::vcf::{GermlineRow{locus,ref_base,call}, write_germline_vcf}`; `provenance::{RunManifest, FileHash, blake3_file, write_manifest}`.

---

## File Structure
- `src/io/bam.rs` — `read_bam_as_core_reads` + `BamSource` (ReadSource). The only htslib→core boundary.
- `src/io/mod.rs` — add `pub mod bam;`.
- `src/call/pipeline.rs` — `call_germline_region<S: ReadSource>`.
- `src/call/mod.rs` — add `pub mod pipeline;` + re-export `call_germline_region`.
- `src/main.rs` — rewrite `run_variants`; add a private `legacy_read_to_core` helper.

---

### Task 1: `io::bam::BamSource` — the htslib → `core::AlignedRead` boundary

**Files:**
- Create: `src/io/bam.rs`
- Modify: `src/io/mod.rs` (add `pub mod bam;`)
- Test: in `src/io/bam.rs`

- [ ] **Step 1: Write the failing test.** Create `src/io/bam.rs`:

```rust
//! The BAM → `core::AlignedRead` boundary. `rust_htslib` lives here and never
//! leaks past this module: `BamSource` yields canonical `core::AlignedRead`s and
//! implements `pileup::ReadSource`, so the kernel stays htslib-free.
//!
//! This adapter pre-loads + sorts the records (it is not yet bounded-memory; a
//! streaming `.csi`-fetch BAM source is a later phase). SEQ is taken
//! forward-oriented per the SAM spec — no reverse-complement is applied.

use std::path::Path;
use std::sync::Arc;

use rust_htslib::bam::record::Cigar as BamCigar;
use rust_htslib::bam::{self, Read as BamRead};

use crate::core::{AlignedRead, CigarOp, CigarOpKind, ContigSet, CoreError, Position, SamFlags};
use crate::pileup::ReadSource;

/// Read all mapped records of a BAM into canonical `core::AlignedRead`s, mapping
/// each record's reference name to a contig id via `contigs`. Records that are
/// unmapped, have no tid, or whose reference is absent from `contigs` are
/// skipped. SEQ is uppercased and kept forward-oriented; strand is recorded in
/// `flags` only.
pub fn read_bam_as_core_reads(
    path: &Path,
    contigs: &ContigSet,
) -> Result<Vec<AlignedRead>, CoreError> {
    let mut reader = bam::Reader::from_path(path)
        .map_err(|e| CoreError::MalformedRecord(format!("open BAM {}: {e}", path.display())))?;
    let header = reader.header().to_owned();

    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::MalformedRecord(e.to_string()))?;
        if rec.is_unmapped() {
            continue;
        }
        let tid = rec.tid();
        if tid < 0 {
            continue;
        }
        let name = std::str::from_utf8(header.tid2name(tid as u32))
            .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
        let contig = match contigs.by_name(name) {
            Some(c) => c.id,
            None => continue,
        };
        let pos0 = rec.pos();
        if pos0 < 0 {
            continue;
        }

        let mut cigar = Vec::new();
        for c in rec.cigar().iter() {
            let (kind, len) = match *c {
                BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => {
                    (CigarOpKind::Match, l)
                }
                BamCigar::Ins(l) => (CigarOpKind::Insertion, l),
                BamCigar::Del(l) => (CigarOpKind::Deletion, l),
                BamCigar::RefSkip(l) => (CigarOpKind::RefSkip, l),
                BamCigar::SoftClip(l) => (CigarOpKind::SoftClip, l),
                BamCigar::HardClip(l) => (CigarOpKind::HardClip, l),
                BamCigar::Pad(l) => (CigarOpKind::Pad, l),
            };
            cigar.push(CigarOp::new(kind, len));
        }

        let seq: Vec<u8> = rec.seq().as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
        let qual: Vec<u8> = rec.qual().to_vec();

        out.push(AlignedRead {
            contig,
            pos: Position(pos0 as u32),
            mapq: rec.mapq(),
            flags: SamFlags(rec.flags()),
            cigar,
            seq: Arc::from(seq.into_boxed_slice()),
            qual: Arc::from(qual.into_boxed_slice()),
        });
    }
    Ok(out)
}

/// A `ReadSource` over a BAM file. Pre-loads and coordinate-sorts the records on
/// construction, then yields them in `(contig, pos)` order.
#[derive(Debug)]
pub struct BamSource {
    reads: std::vec::IntoIter<AlignedRead>,
}

impl BamSource {
    /// Open `path`, convert its mapped records to `core::AlignedRead`s (mapping
    /// reference names via `contigs`), and sort by `(contig, pos)`.
    pub fn new(path: &Path, contigs: &ContigSet) -> Result<Self, CoreError> {
        let mut reads = read_bam_as_core_reads(path, contigs)?;
        reads.sort_by_key(|r| (r.contig, r.pos));
        Ok(Self {
            reads: reads.into_iter(),
        })
    }
}

impl ReadSource for BamSource {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        Ok(self.reads.next())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp(name: &str) -> std::path::PathBuf {
        let ts = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        std::env::temp_dir().join(format!("rosalind-bamsrc-{name}-{ts}.bam"))
    }

    fn one_contig() -> ContigSet {
        let mut c = ContigSet::new();
        c.push("chr1", 1000);
        c
    }

    #[test]
    fn reads_bam_into_core_reads_with_flags_and_cigar() {
        let path = tmp("basic");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );
        {
            let mut w = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            // Forward read at pos 10, 3M.
            let mut fwd = Record::new();
            fwd.set(
                b"fwd",
                Some(&CigarString::from(vec![Cigar::Match(3)])),
                b"ACG",
                b"III",
            );
            fwd.set_tid(0);
            fwd.set_pos(10);
            fwd.set_flags(0);
            fwd.set_mapq(60);
            w.write(&fwd).unwrap();
            // Reverse read at pos 20, 2M.
            let mut rev = Record::new();
            rev.set(
                b"rev",
                Some(&CigarString::from(vec![Cigar::Match(2)])),
                b"TT",
                b"II",
            );
            rev.set_tid(0);
            rev.set_pos(20);
            rev.set_flags(0x10); // REVERSE
            rev.set_mapq(40);
            w.write(&rev).unwrap();
        }

        let reads = read_bam_as_core_reads(&path, &one_contig()).unwrap();
        assert_eq!(reads.len(), 2);
        let fwd = reads.iter().find(|r| r.pos.0 == 10).unwrap();
        assert_eq!(fwd.contig, 0);
        assert_eq!(fwd.mapq, 60);
        assert!(!fwd.flags.is_reverse());
        assert_eq!(fwd.cigar, vec![CigarOp::new(CigarOpKind::Match, 3)]);
        assert_eq!(&fwd.seq[..], b"ACG");
        let rev = reads.iter().find(|r| r.pos.0 == 20).unwrap();
        assert!(rev.flags.is_reverse());
        // SEQ stays forward-oriented (no reverse-complement applied).
        assert_eq!(&rev.seq[..], b"TT");

        // BamSource yields the same reads, coordinate-sorted.
        let mut src = BamSource::new(&path, &one_contig()).unwrap();
        let first = src.next_read().unwrap().unwrap();
        assert_eq!(first.pos.0, 10);
        let second = src.next_read().unwrap().unwrap();
        assert_eq!(second.pos.0, 20);
        assert!(src.next_read().unwrap().is_none());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn unknown_contig_records_are_skipped() {
        let path = tmp("unknown");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chrX")
                .push_tag(b"LN", &1000),
        );
        {
            let mut w = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let mut r = Record::new();
            r.set(b"x", Some(&CigarString::from(vec![Cigar::Match(2)])), b"AC", b"II");
            r.set_tid(0);
            r.set_pos(5);
            r.set_flags(0);
            r.set_mapq(60);
            w.write(&r).unwrap();
        }
        // ContigSet has chr1 only → the chrX record is skipped.
        let reads = read_bam_as_core_reads(&path, &one_contig()).unwrap();
        assert!(reads.is_empty());
        std::fs::remove_file(&path).ok();
    }
}
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test io::bam`
Expected: FAIL — `src/io/bam.rs` not declared (module unknown) / functions undefined.

- [ ] **Step 3: Register the module.** In `src/io/mod.rs`, add after `pub mod vcf;`:

```rust
pub mod bam;
```

- [ ] **Step 4: Run the tests to verify they pass.** Run: `cargo test io::bam`
Expected: PASS (2 tests). If the `match *c` over `BamCigar` errors with "non-exhaustive" (the enum is `#[non_exhaustive]` in this `rust_htslib` version), add a final `_ => continue,` arm. Then `cargo build 2>&1 | grep -i warning` (no `src/io/` warnings), `cargo fmt --all -- --check`.

- [ ] **Step 5: Commit.**

```bash
git add src/io/bam.rs src/io/mod.rs
git commit -m "feat(io/bam): BamSource — htslib to core::AlignedRead boundary (ReadSource)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `call::pipeline::call_germline_region` — engine → caller orchestration

**Files:**
- Create: `src/call/pipeline.rs`
- Modify: `src/call/mod.rs` (add `pub mod pipeline;` + re-export)
- Test: in `src/call/pipeline.rs`

- [ ] **Step 1: Write the failing test.** Create `src/call/pipeline.rs`:

```rust
//! Orchestration: drive a `ReadSource` through the `PileupEngine` and call a
//! germline genotype at every covered position, collecting the variant sites
//! (hom-ref columns are abstained on by `call_germline` and never appear).

use std::ops::Range;
use std::sync::Arc;

use crate::call::{call_germline, GermlineCall, GermlineParams};
use crate::core::{CoreError, Locus};
use crate::pileup::{PileupEngine, PileupParams, ReadSource};

/// Call germline variants across `region` of `contig`. Returns one entry per
/// emitted site as `(locus, ref_base, call)`; the caller pairs these into VCF
/// rows. Hom-ref / no-evidence positions are abstained on (absent from output).
pub fn call_germline_region<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<Vec<(Locus, u8, GermlineCall)>, CoreError> {
    let engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut out = Vec::new();
    for column in engine {
        let column = column?;
        if let Some(call) = call_germline(&column, germline_params) {
            out.push((column.locus, column.ref_base, call));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::Genotype;
    use crate::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
    use crate::pileup::SliceSource;

    fn read(pos: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        let flags = if reverse {
            SamFlags(SamFlags::REVERSE)
        } else {
            SamFlags::default()
        };
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags,
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![35u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn calls_a_het_snv_including_a_reverse_strand_read() {
        // Reference AAAA over [0,4). At position 1, reads split A (ref) / C (alt).
        // One alt read is reverse-strand — proving the end-to-end reverse-strand
        // fix (it must contribute a C, not a complemented G).
        let reference: Arc<[u8]> = Arc::from(b"AAAA".to_vec().into_boxed_slice());
        let reads = vec![
            read(0, b"AAAA", false), // all ref
            read(0, b"ACAA", false), // alt C at pos 1 (forward)
            read(0, b"ACAA", true),  // alt C at pos 1 (reverse — SEQ already forward)
            read(0, b"ACAA", false),
            read(0, b"AAAA", false),
            read(0, b"ACAA", true),
        ];
        let source = SliceSource::new(reads);
        let calls = call_germline_region(
            source,
            reference,
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();

        // Exactly one variant site, at position 1, alt C, heterozygous.
        assert_eq!(calls.len(), 1);
        let (locus, ref_base, call) = &calls[0];
        assert_eq!(locus.pos, Position(1));
        assert_eq!(*ref_base, b'A');
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.genotype, Genotype::Het);
        assert!(call.dp >= 6);
    }

    #[test]
    fn pure_reference_yields_no_calls() {
        let reference: Arc<[u8]> = Arc::from(b"ACGT".to_vec().into_boxed_slice());
        let reads = vec![read(0, b"ACGT", false), read(0, b"ACGT", false)];
        let source = SliceSource::new(reads);
        let calls = call_germline_region(
            source,
            reference,
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();
        assert!(calls.is_empty());
    }
}
```

- [ ] **Step 2: Run the test to verify it fails.** Run: `cargo test call::pipeline`
Expected: FAIL — `pipeline` module not declared / `call_germline_region` undefined.

- [ ] **Step 3: Register the module.** In `src/call/mod.rs`, add `pub mod pipeline;` (next to the other `pub mod` lines) and add to the re-export block:

```rust
pub use pipeline::call_germline_region;
```

- [ ] **Step 4: Run the tests to verify they pass.** Run: `cargo test call::pipeline`
Expected: PASS (2 tests). Then `cargo test` (full suite green), `cargo build 2>&1 | grep -i warning` (none new), `cargo fmt --all -- --check`, `cargo clippy --lib 2>&1 | grep -A2 'src/call/pipeline'`.

- [ ] **Step 5: Commit.**

```bash
git add src/call/pipeline.rs src/call/mod.rs
git commit -m "feat(call/pipeline): germline calling over the pileup engine (ReadSource -> calls)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Rewire `main.rs::run_variants` onto the new vertical + emit a manifest

**Files:**
- Modify: `src/main.rs` (rewrite `run_variants`; add `legacy_read_to_core`)
- Test: the CLI e2e (replicate `.github/workflows/ci.yml`'s `cli-e2e` steps locally) + full suite

**Context:** `main.rs` imports legacy types via `use rosalind::genomics::{… AlignedRead, CigarOp, CigarOpKind …}`, so those identifiers are the **legacy** ones. New code must use fully-qualified `rosalind::core::…`, `rosalind::call::…`, `rosalind::io::…`, `rosalind::pileup::…`, `rosalind::provenance::…`. Keep the `run_variants` signature identical so its callers (the `Commands::Variants` match arm) are unaffected; `block_size` is now unused by the streaming engine (prefix it `_block_size` in the signature or `let _ = block_size;`).

- [ ] **Step 1: Add the legacy→core read converter.** Add this private helper near `read_bam_alignment_file` in `src/main.rs`:

```rust
/// Convert a legacy `genomics::AlignedRead` into a canonical `core::AlignedRead`
/// for the new calling vertical. The legacy type has no RefSkip/Pad CIGAR ops.
fn legacy_read_to_core(r: &AlignedRead, contig: u32) -> rosalind::core::AlignedRead {
    use rosalind::core::{CigarOp as CoreOp, CigarOpKind as CoreKind};
    let cigar = r
        .cigar
        .iter()
        .map(|op| {
            let kind = match op.kind {
                CigarOpKind::Match => CoreKind::Match,
                CigarOpKind::Insertion => CoreKind::Insertion,
                CigarOpKind::Deletion => CoreKind::Deletion,
                CigarOpKind::SoftClip => CoreKind::SoftClip,
                CigarOpKind::HardClip => CoreKind::HardClip,
            };
            CoreOp::new(kind, op.len)
        })
        .collect();
    let flags = if r.is_reverse {
        rosalind::core::SamFlags(rosalind::core::SamFlags::REVERSE)
    } else {
        rosalind::core::SamFlags::default()
    };
    rosalind::core::AlignedRead {
        contig,
        pos: rosalind::core::Position(r.pos),
        mapq: r.mapq,
        flags,
        cigar,
        seq: std::sync::Arc::clone(&r.sequence),
        qual: std::sync::Arc::clone(&r.qualities),
    }
}
```

(If the legacy `genomics::CigarOpKind` has variants beyond those five, add the missing arms — map any reference-skipping op to `CoreKind::RefSkip`, padding to `CoreKind::Pad`.)

- [ ] **Step 2: Replace the body of `run_variants`.** Replace the entire `run_variants` function body (keep the signature; rename `block_size` → `_block_size`) with:

```rust
fn run_variants(
    reference_path: PathBuf,
    alignments_path: PathBuf,
    chrom: Option<String>,
    region_start: u32,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    _block_size: usize,
    quality_threshold: f32,
) -> Result<()> {
    use rosalind::call::{call_germline_region, GermlineParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::{write_germline_vcf, GermlineRow};
    use rosalind::pileup::{PileupParams, SliceSource};
    use rosalind::provenance::{blake3_file, write_manifest, FileHash, RunManifest};

    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let chrom_name = chrom.unwrap_or_else(|| fasta.name.clone());
    let chrom_arc: Arc<str> = chrom_name.clone().into();
    let reference: Arc<[u8]> = Arc::from(fasta.sequence.into_boxed_slice());

    // Single-contig run: one contig spanning the reference window.
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(chrom_name.clone(), region_start + reference.len() as u32);

    let region = region_start..(region_start + reference.len() as u32);
    let pileup_params = PileupParams {
        min_mapq: mapq_threshold, // now actually honored by the engine
        ..PileupParams::default()
    };
    let germline_params = GermlineParams {
        min_qual: quality_threshold as f64,
        ..GermlineParams::default()
    };

    let is_bam = alignments_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("bam"))
        .unwrap_or(false);

    let sites = if is_bam {
        let source = BamSource::new(&alignments_path, &contigs)
            .map_err(|e| anyhow!("failed to read BAM {}: {e}", alignments_path.display()))?;
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (BAM): {e}"))?
    } else {
        let legacy = read_alignment_file(&alignments_path, Some(&chrom_arc)).with_context(|| {
            format!("failed to read SAM alignments from {}", alignments_path.display())
        })?;
        let core_reads: Vec<rosalind::core::AlignedRead> = legacy
            .iter()
            .map(|r| legacy_read_to_core(r, contig_id))
            .collect();
        let source = SliceSource::new(core_reads);
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (SAM): {e}"))?
    };

    let rows: Vec<GermlineRow> = sites
        .into_iter()
        .map(|(locus, ref_base, call)| GermlineRow {
            locus,
            ref_base,
            call,
        })
        .collect();

    match output {
        Some(path) => {
            let file = File::create(&path)
                .with_context(|| format!("failed to create VCF file {}", path.display()))?;
            let mut writer = io::BufWriter::new(file);
            write_germline_vcf(&mut writer, &contigs, &chrom_name, &rows)?;
            writer.flush()?;
            drop(writer);

            // Reproducibility receipt next to the VCF.
            let mut manifest = RunManifest::new("variants");
            manifest.inputs.push(FileHash {
                path: reference_path.display().to_string(),
                blake3: blake3_file(&reference_path)?,
            });
            manifest.inputs.push(FileHash {
                path: alignments_path.display().to_string(),
                blake3: blake3_file(&alignments_path)?,
            });
            manifest.outputs.push(FileHash {
                path: path.display().to_string(),
                blake3: blake3_file(&path)?,
            });
            manifest
                .params
                .insert("mapq_threshold".to_string(), mapq_threshold.to_string());
            manifest
                .params
                .insert("min_qual".to_string(), (quality_threshold as f64).to_string());
            manifest
                .params
                .insert("region_start".to_string(), region_start.to_string());
            let manifest_path = write_manifest(&path, &manifest)?;
            eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            write_germline_vcf(&mut handle, &contigs, &chrom_name, &rows)?;
        }
    }

    Ok(())
}
```

- [ ] **Step 3: Build and run the full suite.** Run: `cargo build 2>&1 | grep -iE 'error|warning'` (expect none), then `cargo test` (full suite green — the legacy tests are untouched and still pass). Fix any compile errors (most likely: an unused legacy import now that `StreamingVariantCaller`/`write_vcf` may no longer be used by `run_variants` — if `cargo build` warns that `StreamingVariantCaller` or `write_vcf` is now unused in `main.rs`, remove only those names from the `use rosalind::genomics::{…}` import list; do NOT remove anything still used elsewhere in `main.rs`).

- [ ] **Step 4: Replicate the CLI e2e locally.** Read `.github/workflows/ci.yml`'s `cli-e2e` job and run the same steps against `examples/data/illumina_toy` (generate toy data → index/align → sort → `cargo run -- variants … --output …/variants.vcf`). Confirm:
  - `grep -q '^#CHROM' …/variants.vcf` succeeds and the file now starts with `##fileformat=VCFv4.2`.
  - There is at least one non-`#` record line (or the documented "no variant lines" warning path still holds).
  - A `…/variants.vcf.manifest.json` sidecar was written and is valid JSON (`python3 -c "import json;json.load(open('…/variants.vcf.manifest.json'))"`).
Report the exact commands you ran and their output.

- [ ] **Step 5: Run fmt + clippy.** `cargo fmt --all -- --check` (run `cargo fmt --all` + amend if it drifts), `cargo clippy --bin rosalind 2>&1 | grep -A2 run_variants` (fix obvious lints).

- [ ] **Step 6: Commit.**

```bash
git add src/main.rs
git commit -m "feat(cli): route 'variants' through the new engine + calibrated caller + VCF + receipt" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage (§4, §8 — germline subset):**
- §4 germline data flow (BAM: BamSource→engine→call::germline→VCF+manifest; SAM: load→sort→SliceSource→same) — Tasks 1–3. ✔
- §8 "StreamingVariantCaller re-implemented over PileupEngine" — realized as `call::pipeline::call_germline_region` + CLI wiring; the legacy struct is retained (not deleted) so its tests stay green (deletion deferred + documented). ✔ (partial, by design)
- §8 "parse_cigar generalized" — already done in the earlier CI fix. ✔
- The reverse-strand + CIGAR + MAPQ-honoring fixes now reach the CLI (the engine applies them; `min_mapq` is wired from `mapq_threshold`). ✔

**Deferred (documented, NOT in this plan):** delete `statistics::bayesian_variant_caller` + legacy `StreamingVariantCaller`/`vcf.rs` (needs migrating `variant_pipeline.rs`/`golden_vcf.rs`); somatic CLI rewiring onto `call::somatic` + `io::vcf` somatic (+ `truthset_validation.rs` stays on the legacy `SomaticCaller`); bounded-memory streaming `BamSource` (Phase D, with `.csi` fetch + RSS gate).

**Placeholder scan:** none — complete code for Tasks 1–2; exact replacement function for Task 3.

**Type consistency:** uses `core::AlignedRead`/`SamFlags(u16)`/`CigarOpKind` (incl. RefSkip/Pad), `pileup::{PileupEngine::new(source, Arc<[u8]>, u32, Range<u32>, PileupParams), PileupParams, ReadSource, SliceSource}`, `call::{call_germline_region, GermlineParams}`, `io::vcf::{GermlineRow, write_germline_vcf}`, `provenance::{RunManifest, FileHash, blake3_file, write_manifest}` — all as shipped in A1–A4. Legacy `genomics::{AlignedRead, CigarOpKind}` used only in `legacy_read_to_core`.

**Safety:** legacy calling path retained → `variant_pipeline.rs`, `golden_vcf.rs`, `truthset_validation.rs`, and the `variant_caller`/`statistics` unit tests stay green; the e2e is format-agnostic (greps `^#CHROM` + a non-comment line) so the new spec-valid VCF passes. No existing public API removed.
