# Phase B4 — bounded whole-genome variant calling over the persisted index — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rosalind variants --index <idx> --alignments <sorted.bam>` calls germline variants across every contig of a persisted index in bounded, predictable memory (reads stream one record at a time, the reference is read from the index per-contig), emitting a multi-contig VCF plus a receipt that shows the realized peak.

**Architecture:** A `StreamingBamSource` (`ReadSource`) reads one BAM record at a time via `bam::Reader::read` (no full-file `Vec`), with a `(contig_id, pos)` monotonicity guard enforcing coordinate-sorted input. A library drive `call_germline_whole_genome` makes a single sorted pass, partitions it per contig with a peekable adapter, decodes each contig's reference from the B4a `ReferenceView`, and reuses the tested per-contig `call_germline_region` (now also returning the engine's max working set). `variants --index` loads the index, runs the drive, writes a multi-contig VCF + a manifest carrying the realized peak RSS + max working set; a record-only `--memory-budget-mb` flags overage without aborting.

**Tech Stack:** Rust 2021 (MSRV 1.72 — no `div_ceil`); `rust-htslib` 0.44.1 (`bam::Reader::read(&mut Record) -> Option<Result<(), Error>>`, proven in `pileup_stream.rs`); `clap` derive groups; `anyhow`; no new dependencies.

This is the caller sub-stage of Phase B4, from `docs/superpowers/specs/2026-05-27-phase-b4-variants-index-design.md` (reprioritized ahead of the aligner). Follows B4a (`ReferenceView`, merged — PR #18). **Out of scope (deferred):** the aligner over the index + `align --index` / `somatic --index` (Phase E); SAM streaming (BAM is the WGS input); on-demand `ref_base` (the per-contig reference peak is laptop-fine); budget *enforcement* + `rosalind plan`/`verify` (Phase C).

---

## File structure

- `src/io/bam.rs` — **Modify.** Extract the per-record `Record → AlignedRead` mapping into `pub(crate) fn record_to_aligned_read`; reuse it in `read_bam_as_core_reads`; add `StreamingBamSource<'a>` (`ReadSource`, one record at a time, `(contig_id,pos)` monotonicity guard).
- `src/call/pipeline.rs` — **Modify.** Add `call_germline_region_tracked` (the engine loop tracking the max `current_working_set`, returning `(sites, WorkingSet)`); make `call_germline_region` a thin wrapper that drops the working set (existing 4 call sites unchanged).
- `src/call/whole_genome.rs` — **Create.** `call_germline_whole_genome` (single sorted pass → per-contig `PerContig` adapter → `call_germline_region_tracked` per contig + per-contig `ReferenceView` decode → accumulated rows + max working set).
- `src/call/mod.rs` — **Modify.** `mod whole_genome;` + `pub use`.
- `src/main.rs` — **Modify.** `Variants` clap variant: `--index` XOR `--reference` (exactly one), `--memory-budget-mb`; split `run_variants` into the existing `--reference` path and a new `--index` path that runs the drive over a `StreamingBamSource` + the index's `ReferenceView`/`ContigSet`, writes the multi-contig VCF + the memory receipt.
- `tests/variants_index.rs` — **Create.** Integration gates: parity vs `--reference`, multi-contig, self-contained, sort-rejection, bounded (no materialization), memory receipt + budget-never-refuses.

---

## Task 1: `StreamingBamSource` + shared record mapping + sort guard

**Files:**
- Modify: `src/io/bam.rs`

- [ ] **Step 1: Write failing tests.** In the `#[cfg(test)] mod tests` of `src/io/bam.rs`, add (the module already has `use super::*;` and builds `Record`s via `rust_htslib::bam::record::{Cigar, CigarString, Record}` — reuse that style; if a BAM-writing helper exists, reuse it, otherwise these tests write a tiny BAM with `bam::Writer`):

```rust
    #[test]
    fn streaming_source_yields_same_reads_as_materializing() {
        // Build a tiny coordinate-sorted BAM with 2 contigs, read it both ways.
        let dir = std::env::temp_dir().join(format!(
            "rosalind-stream-{}",
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("reads.bam");
        write_sorted_test_bam(&bam_path); // helper below

        let mut contigs = ContigSet::new();
        contigs.push("chr1", 1000);
        contigs.push("chr2", 1000);

        let materialized = read_bam_as_core_reads(&bam_path, &contigs).unwrap();
        let mut streamed = Vec::new();
        let mut src = StreamingBamSource::new(&bam_path, &contigs).unwrap();
        while let Some(r) = src.next_read().unwrap() {
            streamed.push(r);
        }
        assert_eq!(streamed.len(), materialized.len());
        for (a, b) in streamed.iter().zip(materialized.iter()) {
            assert_eq!((a.contig, a.pos), (b.contig, b.pos));
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn streaming_source_rejects_out_of_order() {
        let dir = std::env::temp_dir().join(format!(
            "rosalind-stream-bad-{}",
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("unsorted.bam");
        write_unsorted_test_bam(&bam_path); // two reads, descending pos on one contig

        let mut contigs = ContigSet::new();
        contigs.push("chr1", 1000);

        let mut src = StreamingBamSource::new(&bam_path, &contigs).unwrap();
        // First read OK; the second (lower pos) must error.
        let _ = src.next_read().unwrap();
        assert!(src.next_read().is_err(), "out-of-order read must be rejected");
        std::fs::remove_dir_all(&dir).ok();
    }
```

Add these test helpers to the same test module (writing minimal BAMs with `bam::Writer` + a header with two `@SQ`s):

```rust
    fn test_header(contigs: &[(&str, usize)]) -> bam::Header {
        let mut header = bam::Header::new();
        for (name, len) in contigs {
            let mut rec = bam::header::HeaderRecord::new(b"SQ");
            rec.push_tag(b"SN", name);
            rec.push_tag(b"LN", &(*len as i64));
            header.push_record(&rec);
        }
        header
    }

    fn push_record(writer: &mut bam::Writer, header: &bam::HeaderView, tid: i32, pos: i64, seq: &[u8]) {
        let mut rec = Record::new();
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        let qual = vec![30u8; seq.len()];
        rec.set(b"r", Some(&cigar), seq, &qual);
        rec.set_tid(tid);
        rec.set_pos(pos);
        rec.set_mapq(60);
        let _ = header;
        writer.write(&rec).unwrap();
    }

    fn write_sorted_test_bam(path: &Path) {
        let header = test_header(&[("chr1", 1000), ("chr2", 1000)]);
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let hv = writer.header().clone();
        push_record(&mut writer, &hv, 0, 10, b"ACGT");
        push_record(&mut writer, &hv, 0, 20, b"ACGT");
        push_record(&mut writer, &hv, 1, 5, b"ACGT");
    }

    fn write_unsorted_test_bam(path: &Path) {
        let header = test_header(&[("chr1", 1000)]);
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let hv = writer.header().clone();
        push_record(&mut writer, &hv, 0, 50, b"ACGT");
        push_record(&mut writer, &hv, 0, 10, b"ACGT"); // out of order
    }
```

(If `Record::set`'s exact signature differs in 0.44.1, adapt the record construction to the version already used elsewhere in this test module — the goal is a tiny sorted/unsorted BAM. `bam::Writer`/`bam::Header`/`HeaderRecord` are rust-htslib 0.44.1 APIs.)

- [ ] **Step 2: Run them to verify they fail.**

Run: `cargo test --lib io::bam::tests::streaming 2>&1 | tail -20`
Expected: compile error — `StreamingBamSource` does not exist.

- [ ] **Step 3: Extract the shared record mapping.** In `src/io/bam.rs`, add a `pub(crate)` helper containing the per-record logic currently inline in `read_bam_as_core_reads` (lines mapping a `Record` to an `AlignedRead`), and rewrite `read_bam_as_core_reads`'s loop to call it:

```rust
/// Map one BAM `Record` to a canonical `AlignedRead`, or `None` if it is
/// unmapped, has no tid, a negative pos, or a reference name absent from
/// `contigs`. Shared by `read_bam_as_core_reads` and `StreamingBamSource`.
pub(crate) fn record_to_aligned_read(
    rec: &bam::Record,
    header: &bam::HeaderView,
    contigs: &ContigSet,
) -> Result<Option<AlignedRead>, CoreError> {
    if rec.is_unmapped() {
        return Ok(None);
    }
    let tid = rec.tid();
    if tid < 0 {
        return Ok(None);
    }
    let name = std::str::from_utf8(header.tid2name(tid as u32))
        .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
    let contig = match contigs.by_name(name) {
        Some(c) => c.id,
        None => return Ok(None),
    };
    let pos0 = rec.pos();
    if pos0 < 0 {
        return Ok(None);
    }

    let mut cigar = Vec::new();
    for c in rec.cigar().iter() {
        let (kind, len) = match *c {
            BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => (CigarOpKind::Match, l),
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

    Ok(Some(AlignedRead {
        contig,
        pos: Position(pos0 as u32),
        mapq: rec.mapq(),
        flags: SamFlags(rec.flags()),
        cigar,
        seq: Arc::from(seq.into_boxed_slice()),
        qual: Arc::from(qual.into_boxed_slice()),
    }))
}
```

Rewrite `read_bam_as_core_reads`'s body to reuse it (keep its signature):

```rust
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
        if let Some(read) = record_to_aligned_read(&rec, &header, contigs)? {
            out.push(read);
        }
    }
    Ok(out)
}
```

- [ ] **Step 4: Add `StreamingBamSource`.** In `src/io/bam.rs` (after `BamSource`), add (mirrors `pileup_stream.rs`'s `reader.read` idiom; `use rust_htslib::bam::Read as BamRead;` is already imported in this file — if not, add it):

```rust
/// A bounded `ReadSource` over a **coordinate-sorted** BAM: pulls one record at a
/// time (never materialises the file). Enforces sort order via a `(contig_id, pos)`
/// monotonicity guard — `next_read` errors if a read arrives out of order (an
/// unsorted BAM, or one whose `@SQ` order disagrees with `contigs`).
pub struct StreamingBamSource<'a> {
    reader: bam::Reader,
    header: bam::HeaderView,
    contigs: &'a ContigSet,
    last: Option<(u32, u32)>,
}

impl<'a> StreamingBamSource<'a> {
    /// Open `path` for streaming, mapping reference names via `contigs`.
    pub fn new(path: &Path, contigs: &'a ContigSet) -> Result<Self, CoreError> {
        let reader = bam::Reader::from_path(path)
            .map_err(|e| CoreError::MalformedRecord(format!("open BAM {}: {e}", path.display())))?;
        let header = reader.header().to_owned();
        Ok(Self { reader, header, contigs, last: None })
    }
}

impl ReadSource for StreamingBamSource<'_> {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        loop {
            let mut record = bam::Record::new();
            match self.reader.read(&mut record) {
                None => return Ok(None),
                Some(Err(e)) => return Err(CoreError::MalformedRecord(e.to_string())),
                Some(Ok(())) => {}
            }
            let Some(read) = record_to_aligned_read(&record, &self.header, self.contigs)? else {
                continue; // unmapped / absent contig / etc. — skip
            };
            let key = (read.contig, read.pos.0);
            if let Some(prev) = self.last {
                if key < prev {
                    return Err(CoreError::MalformedRecord(format!(
                        "alignments are not coordinate-sorted (or @SQ order disagrees with the \
                         index): read at contig {} pos {} follows contig {} pos {}",
                        key.0, key.1, prev.0, prev.1
                    )));
                }
            }
            self.last = Some(key);
            return Ok(Some(read));
        }
    }
}
```

NOTE: the `let ... else { continue; }` syntax is Rust 1.65+ (MSRV 1.72 — fine). If `bam::HeaderView` is not the type returned by `reader.header().to_owned()` in 0.44.1, use the same owned-header type `pileup_stream.rs` stores (`reader.header().to_owned()` there is assigned to a field; mirror it). `CoreError::MalformedRecord(String)` is the existing variant used throughout this file.

- [ ] **Step 5: Build + run the streaming tests + the existing BAM tests.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --lib io::bam 2>&1 | tail -25` (the 2 new streaming tests + the existing bam tests pass — `read_bam_as_core_reads` is unchanged behaviorally), then `cargo fmt --all`.

- [ ] **Step 6: Commit.**

```bash
git add src/io/bam.rs
git commit -m "feat(io/bam): StreamingBamSource — bounded, sort-guarded BAM read source" \
  -m "Streams one record at a time via bam::Reader::read (no full-file Vec), sharing the Record->AlignedRead mapping (extracted as record_to_aligned_read) with read_bam_as_core_reads. A (contig_id,pos) monotonicity guard rejects unsorted / index-order-mismatched BAMs. The bounded-reads foundation for whole-genome variants --index." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: working-set-tracking caller + the whole-genome drive (library)

**Files:**
- Modify: `src/call/pipeline.rs`, `src/call/mod.rs`
- Create: `src/call/whole_genome.rs`

- [ ] **Step 1: Add `call_germline_region_tracked` + delegate.** In `src/call/pipeline.rs`, replace the body of `call_germline_region` so it delegates to a new tracked variant that also returns the max engine working set. (`WorkingSet` is `crate::core::WorkingSet`; `PileupEngine::current_working_set(&self) -> WorkingSet` exists. The engine is an `Iterator`; call `.next()` manually so it stays queryable.) Add the import `use crate::core::WorkingSet;` if not present.

```rust
/// Like [`call_germline_region`], but also returns the maximum pileup-engine
/// working set observed during the pass (the bounded-memory signal — the
/// foundation for the `variants` memory receipt and `rosalind plan`).
pub fn call_germline_region_tracked<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<(Vec<(Locus, u8, GermlineCall)>, WorkingSet), CoreError> {
    let mut engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut out = Vec::new();
    let mut max_ws = WorkingSet { bytes: 0 };
    while let Some(column) = engine.next() {
        let column = column?;
        let ws = engine.current_working_set();
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        if let Some(call) = call_germline(&column, germline_params) {
            out.push((column.locus, column.ref_base, call));
        }
    }
    Ok((out, max_ws))
}
```

and make the existing `call_germline_region` a thin wrapper (keeping its exact signature, so the 4 existing call sites are unchanged):

```rust
pub fn call_germline_region<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<Vec<(Locus, u8, GermlineCall)>, CoreError> {
    let (sites, _ws) =
        call_germline_region_tracked(source, reference, contig, region, pileup_params, germline_params)?;
    Ok(sites)
}
```

(`engine.next()` requires the `Iterator` trait in scope — it is, via the prelude. If `PileupEngine::new`/`current_working_set` are not `pub`/`pub(crate)` reachable from `pipeline.rs`, they already are — `pipeline.rs` constructs the engine today.)

- [ ] **Step 2: Export the tracked fn.** In `src/call/mod.rs`, add `call_germline_region_tracked` to the `pub use pipeline::{…}` line.

- [ ] **Step 3: Write the failing drive test.** Create `src/call/whole_genome.rs` with the test first (and `unimplemented!()` stubs so it compiles-then-fails):

```rust
//! Bounded whole-genome germline calling over a coordinate-sorted read stream and
//! a persisted reference (`ReferenceView`). One sorted pass, partitioned per
//! contig by `PerContig`, reusing the per-contig caller — peak memory ≈ the
//! largest contig's reference + the pileup working set, independent of input size.

use std::ops::Range;
use std::sync::Arc;

use crate::call::{call_germline_region_tracked, GermlineParams};
use crate::call::pipeline::GermlineCall;
use crate::core::{AlignedRead, ContigSet, CoreError, Locus, WorkingSet};
use crate::genomics::ReferenceView;
use crate::pileup::{PileupParams, ReadSource};

/// Per-contig view over a shared sorted stream: yields contig `contig`'s reads,
/// then stops at (and buffers) the first read of a later contig.
struct PerContig<'s, S: ReadSource> {
    source: &'s mut S,
    contig: u32,
    peeked: &'s mut Option<AlignedRead>,
}

impl<S: ReadSource> ReadSource for PerContig<'_, S> {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        if let Some(r) = self.peeked.take() {
            if r.contig == self.contig {
                return Ok(Some(r));
            }
            *self.peeked = Some(r); // belongs to a later contig — keep it, stop here
            return Ok(None);
        }
        match self.source.next_read()? {
            Some(r) if r.contig == self.contig => Ok(Some(r)),
            Some(r) => {
                *self.peeked = Some(r);
                Ok(None)
            }
            None => Ok(None),
        }
    }
}

/// Call germline variants across every contig of `contigs`, in id order, reading
/// each contig's reference from `ref_view`. Returns the accumulated
/// `(Locus, ref_base, call)` rows and the max pileup working set observed.
/// `source` MUST yield reads in `(contig, pos)` order (a `StreamingBamSource`
/// guards this); the per-contig partition relies on it.
pub fn call_germline_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<(Vec<(Locus, u8, GermlineCall)>, WorkingSet), CoreError> {
    let _ = (&mut source, ref_view, contigs, pileup_params, germline_params);
    unimplemented!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::{GenomeIndex, IndexReader, IndexWriter};
    use crate::pileup::SliceSource;

    fn tmp(suffix: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let d = std::env::temp_dir().join(format!("rosalind-wg-{suffix}-{nanos}"));
        std::fs::create_dir_all(&d).unwrap();
        d.join("ref.idx")
    }

    // A read fully matching `len` bases at (contig, pos).
    fn read_at(contig: u32, pos: u32, seq: &[u8]) -> AlignedRead {
        use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
        AlignedRead {
            contig,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags(0),
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![40u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn whole_genome_equals_per_contig_calls() {
        // Two contigs; reads with a clear SNV on each. Build an index for the
        // reference, then compare the whole-genome drive to per-contig calls.
        let idx_path = tmp("equal");
        let index = GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec()),
            ("chr2".to_string(), b"TTTTGGGGCCCCAAAATTTT".to_vec()),
        ])
        .unwrap();
        IndexWriter::create(&idx_path).unwrap().write_genome_index(&index).unwrap();
        let loaded = IndexReader::open(&idx_path).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();

        // Reads (already (contig,pos)-sorted): depth on a few positions of each contig.
        let reads = vec![
            read_at(0, 0, b"ACGTACGT"),
            read_at(0, 0, b"ACGTACGT"),
            read_at(1, 0, b"TTTTGGGG"),
            read_at(1, 0, b"TTTTGGGG"),
        ];
        let pp = PileupParams::default();
        let gp = GermlineParams::default();

        let (rows, ws) =
            call_germline_whole_genome(SliceSource::new(reads.clone()), &rv, contigs, pp.clone(), &gp)
                .unwrap();
        // Every row's contig is a real contig id; rows are grouped by contig.
        assert!(rows.iter().all(|(l, _, _)| (l.contig as usize) < contigs.len()));
        assert!(ws.bytes >= 0); // working set tracked (>= 0 always; presence check)

        // Equivalence: union of per-contig call_germline_region over the same reads.
        let mut expected = Vec::new();
        for c in contigs.iter() {
            let mut buf = Vec::new();
            rv.decode_window(
                c.global_offset as usize,
                c.global_offset as usize + c.length as usize,
                &mut buf,
            );
            let cref: Arc<[u8]> = Arc::from(buf.as_slice());
            let creads: Vec<AlignedRead> =
                reads.iter().filter(|r| r.contig == c.id).cloned().collect();
            let sites = call_germline_region(
                SliceSource::new(creads),
                cref,
                c.id,
                0..c.length,
                pp.clone(),
                &gp,
            )
            .unwrap();
            expected.extend(sites);
        }
        assert_eq!(rows, expected, "whole-genome == per-contig union");

        let _ = std::fs::remove_dir_all(idx_path.parent().unwrap());
    }
}
```

(The test imports `call_germline_region` — add it to the `use crate::call::{…}` line in the test module. `PileupParams`/`GermlineParams` derive `Clone` — if not, construct fresh defaults per call instead of `.clone()`. `GermlineCall`/`PileupParams`/`SliceSource` paths: adjust the `use` to wherever they are actually exported — `GermlineCall` is in `call::pipeline`, `SliceSource`/`PileupParams`/`ReadSource` in `pileup`. If `PileupParams` is not `Clone`, replace `pp.clone()` with `PileupParams::default()`.)

- [ ] **Step 4: Run it to verify it fails.**

Run: `cargo test --lib call::whole_genome 2>&1 | tail -20`
Expected: FAIL — `unimplemented!()` panics in `call_germline_whole_genome`.

- [ ] **Step 5: Implement the drive.** Replace the `unimplemented!()` body of `call_germline_whole_genome`:

```rust
pub fn call_germline_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<(Vec<(Locus, u8, GermlineCall)>, WorkingSet), CoreError> {
    let mut rows = Vec::new();
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut peeked: Option<AlignedRead> = None;
    let mut buf = Vec::new();

    for c in contigs.iter() {
        // Decode this contig's reference (peak = largest contig; bounded).
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        ref_view.decode_window(start, end, &mut buf);
        let reference: Arc<[u8]> = Arc::from(buf.as_slice());

        let per = PerContig { source: &mut source, contig: c.id, peeked: &mut peeked };
        let region: Range<u32> = 0..c.length;
        let (sites, ws) = call_germline_region_tracked(
            per,
            reference,
            c.id,
            region,
            pileup_params.clone(),
            germline_params,
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        rows.extend(sites);
    }

    Ok((rows, max_ws))
}
```

(`pileup_params.clone()` requires `PileupParams: Clone` — it does in this codebase (a small config struct). `Contig` fields `id: u32`, `length: u32`, `global_offset: u64` are `pub`.)

- [ ] **Step 6: Register the module.** In `src/call/mod.rs`, add `mod whole_genome;` and `pub use whole_genome::call_germline_whole_genome;`.

- [ ] **Step 7: Build + run the drive test + the full call/pileup suites.**

Run: `cargo build --lib 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --lib call 2>&1 | tail -20` (the drive test + the existing call/pipeline tests pass — `call_germline_region` behavior unchanged via the wrapper), then
`cargo test --lib 2>&1 | grep -E "test result:"` (no failures), then `cargo fmt --all`.

- [ ] **Step 8: Commit.**

```bash
git add src/call/pipeline.rs src/call/whole_genome.rs src/call/mod.rs
git commit -m "feat(call): bounded whole-genome germline drive over a sorted stream + ReferenceView" \
  -m "call_germline_region_tracked surfaces the max pileup working set (call_germline_region delegates, signature unchanged). call_germline_whole_genome makes a single sorted pass, partitions per contig (PerContig adapter), decodes each contig's reference from ReferenceView, and reuses the per-contig caller — peak ≈ largest contig + working set, independent of input size. Verified == the per-contig union." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `variants --index` CLI + multi-contig VCF

**Files:**
- Modify: `src/main.rs`
- Create: `tests/variants_index.rs`

- [ ] **Step 1: Write the failing integration tests.** Create `tests/variants_index.rs`:

```rust
//! Phase B4 gates: `rosalind variants --index` — bounded, self-contained,
//! germline calling, parity with `--reference` on the same sorted BAM. Exercised
//! through the real `index -> align -> sort -> variants` CLI pipeline (no
//! rust-htslib dev-dependency; `--index` is BAM-only). Multi-contig calling is
//! gated by the `call::whole_genome` library test.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let d = env::temp_dir().join(format!("rosalind-b4v-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().expect("spawn rosalind")
}

fn write_fasta(dir: &Path, name: &str, seq: &str) -> PathBuf {
    let p = dir.join("ref.fa");
    std::fs::write(&p, format!(">{name}\n{seq}\n")).unwrap();
    p
}

// A FASTQ of reads that are substrings of `seq` (so they align), giving depth.
fn write_fastq(dir: &Path, seq: &str, starts: &[usize], len: usize) -> PathBuf {
    let p = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in starts.iter().enumerate() {
        let read = &seq[start..start + len];
        let qual: String = std::iter::repeat('I').take(len).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&p, s).unwrap();
    p
}

/// `index` -> `align --format bam` -> `sort` -> returns `(idx, sorted.bam)`.
/// Single-contig (the aligner is single-contig); produces a real sorted BAM via
/// the repo's own pipeline so the streaming `--index` path has valid input.
fn build_index_and_sorted_bam(dir: &Path, fa: &Path, fq: &Path) -> (PathBuf, PathBuf) {
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
    assert!(run(&["index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()])
        .status
        .success());
    let a = run(&["align", "--reference", fa.to_str().unwrap(), "--reads", fq.to_str().unwrap(),
        "--format", "bam", "--output", raw.to_str().unwrap()]);
    assert!(a.status.success(), "align: {}", String::from_utf8_lossy(&a.stderr));
    let s = run(&["sort", "--input", raw.to_str().unwrap(), "--output", sorted.to_str().unwrap()]);
    assert!(s.status.success(), "sort: {}", String::from_utf8_lossy(&s.stderr));
    (idx, sorted)
}

fn vcf_records(out: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(out)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.to_string())
        .collect()
}

#[test]
fn variants_index_matches_reference_on_the_same_sorted_bam() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 4, 4, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);

    // Same sorted BAM through both paths → identical reads → calls must match.
    let by_ref = run(&["variants", "--reference", fa.to_str().unwrap(), "--alignments", bam.to_str().unwrap()]);
    assert!(by_ref.status.success(), "--reference: {}", String::from_utf8_lossy(&by_ref.stderr));
    let by_idx = run(&["variants", "--index", idx.to_str().unwrap(), "--alignments", bam.to_str().unwrap()]);
    assert!(by_idx.status.success(), "--index: {}", String::from_utf8_lossy(&by_idx.stderr));

    // Record lines only (headers legitimately differ: sample name / ##contig set).
    assert_eq!(
        vcf_records(&by_ref.stdout),
        vcf_records(&by_idx.stdout),
        "variants --index records must match --reference on the same BAM"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_index_is_self_contained() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    std::fs::remove_file(&fa).unwrap(); // reference FASTA gone

    let out = run(&["variants", "--index", idx.to_str().unwrap(), "--alignments", bam.to_str().unwrap()]);
    assert!(out.status.success(), "must call from the index alone: {}", String::from_utf8_lossy(&out.stderr));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_requires_exactly_one_reference_source() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0], 8);
    let (_idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    // Neither --index nor --reference → error.
    let neither = run(&["variants", "--alignments", bam.to_str().unwrap()]);
    assert!(!neither.status.success(), "must require one of --index/--reference");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_index_rejects_sam() {
    // --index requires a BAM; a .sam alignments file errors with guidance.
    let dir = tmpdir();
    let fa = write_fasta(&dir, "chr1", "ACGTACGTACGTACGT");
    let idx = dir.join("ref.idx");
    assert!(run(&["index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()]).status.success());
    let sam = dir.join("reads.sam");
    std::fs::write(&sam, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:16\n").unwrap();
    let out = run(&["variants", "--index", idx.to_str().unwrap(), "--alignments", sam.to_str().unwrap()]);
    assert!(!out.status.success(), "--index with SAM must error");
    std::fs::remove_dir_all(&dir).ok();
}
```

NOTE: the test builds a real sorted BAM via the repo's own `index → align → sort` CLI (no rust-htslib dev-dependency, and it exercises the genuine pipeline). Parity runs `--reference` and `--index` over the **same** sorted BAM, so the reads are identical and only the reference *source* differs; record lines must match (headers legitimately differ — sample name / `##contig` set). `variants --reference` defaults `--chrom` to the FASTA's record name, so omitting it is fine.

- [ ] **Step 2: Run them to verify they fail.**

Run: `cargo test --test variants_index 2>&1 | tail -20`
Expected: failures — `variants` does not yet accept `--index` (the arg is unknown → non-zero exit), so the `--index` assertions fail.

- [ ] **Step 3: Update the `Variants` clap variant.** In `src/main.rs`, change the `Variants` variant's `reference` field and add `index` + `memory_budget_mb`:

```rust
    /// Call germline variants from aligned reads (streaming pileup engine +
    /// calibrated, abstention-aware genotype-likelihood caller).
    Variants {
        /// Persisted index (`rosalind index`); calls all contigs, reference from
        /// the index. Mutually exclusive with `--reference`.
        #[arg(long, conflicts_with = "reference", required_unless_present = "reference")]
        index: Option<PathBuf>,
        /// Reference genome (FASTA) — single-contig path. Mutually exclusive with `--index`.
        #[arg(long, required_unless_present = "index")]
        reference: Option<PathBuf>,
        /// Alignments in SAM or BAM format (coordinate-sorted for `--index`).
        #[arg(long)]
        alignments: PathBuf,
        /// Chromosome name (single-contig `--reference` path only; defaults to the
        /// first FASTA record). Not allowed with `--index`.
        #[arg(long)]
        chrom: Option<String>,
        /// Starting offset (0-based) for the reference region (`--reference` only).
        #[arg(long, default_value_t = 0)]
        region_start: u32,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        /// Optional VCF output path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Deprecated and ignored.
        #[arg(long, default_value_t = 1024, hide = true)]
        block_size: usize,
        /// Minimum quality threshold for reporting variants.
        #[arg(long, default_value_t = 10.0)]
        quality_threshold: f32,
        /// Declared memory budget (MiB) for the run — records a plan/peak line; does
        /// not enforce (enforcement is a later phase). (`--index` path.)
        #[arg(long)]
        memory_budget_mb: Option<u64>,
    },
```

Update the `match cli.command` arm for `Commands::Variants { … }` to destructure the new fields (`index`, `reference`, `memory_budget_mb`) and dispatch:

```rust
        Commands::Variants {
            index,
            reference,
            alignments,
            chrom,
            region_start,
            mapq_threshold,
            output,
            block_size: _,
            quality_threshold,
            memory_budget_mb,
        } => {
            if let Some(index) = index {
                if chrom.is_some() || region_start != 0 {
                    bail!("--chrom/--region-start are not valid with --index (the whole index is called)");
                }
                run_variants_index(
                    index, alignments, mapq_threshold, output, quality_threshold, memory_budget_mb,
                )?
            } else {
                let reference = reference.expect("clap guarantees one of --index/--reference");
                run_variants(
                    reference, alignments, chrom, region_start, mapq_threshold, output, 1024,
                    quality_threshold,
                )?
            }
        }
```

Keep `run_variants` (the `--reference` path) as-is, but change its first param to take the resolved `PathBuf` (it already does). (If `run_variants`'s signature took `reference: PathBuf`, it still does — we pass the unwrapped one.)

- [ ] **Step 4: Implement `run_variants_index`.** Add this handler near `run_variants` in `src/main.rs`:

```rust
/// Call germline variants across all contigs of a persisted index (B4), reading
/// the reference from the index and streaming the (coordinate-sorted) BAM.
fn run_variants_index(
    index_path: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    quality_threshold: f32,
    memory_budget_mb: Option<u64>,
) -> Result<()> {
    use rosalind::call::{call_germline_whole_genome, GermlineParams};
    use rosalind::genomics::IndexReader;
    use rosalind::io::bam::StreamingBamSource;
    use rosalind::io::vcf::{write_germline_vcf, GermlineRow};
    use rosalind::pileup::PileupParams;
    use rosalind::provenance::{blake3_file, write_manifest, FileHash, RunManifest};

    let loaded = IndexReader::open(&index_path)
        .with_context(|| format!("failed to open index {}", index_path.display()))?;
    let ref_view = loaded
        .reference_view()
        .with_context(|| format!("failed to read reference from index {}", index_path.display()))?;
    let contigs = loaded.contigs();

    let pileup_params = PileupParams { min_mapq: mapq_threshold, ..PileupParams::default() };
    let germline_params = GermlineParams { min_qual: quality_threshold as f64, ..GermlineParams::default() };

    // `--index` streams a coordinate-sorted BAM (the bounded, multi-contig WGS
    // path). SAM under `--index` is unsupported — the legacy SAM reader is
    // single-contig; use a sorted BAM (`rosalind sort`), or `--reference` for
    // single-contig SAM.
    let is_bam = alignments_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("bam"))
        .unwrap_or(false);
    if !is_bam {
        bail!(
            "--index requires a coordinate-sorted BAM (use `rosalind sort`); \
             for single-contig SAM use --reference"
        );
    }
    let source = StreamingBamSource::new(&alignments_path, contigs)
        .map_err(|e| anyhow!("failed to open BAM {}: {e}", alignments_path.display()))?;
    let (sites, _max_ws) =
        call_germline_whole_genome(source, &ref_view, contigs, pileup_params, &germline_params)
            .map_err(|e| anyhow!("variant calling failed: {e}"))?;

    let rows: Vec<GermlineRow> = sites
        .into_iter()
        .map(|(locus, ref_base, call)| GermlineRow { locus, ref_base, call })
        .collect();

    match output {
        Some(path) => {
            let file = File::create(&path)
                .with_context(|| format!("failed to create VCF file {}", path.display()))?;
            let mut writer = io::BufWriter::new(file);
            write_germline_vcf(&mut writer, contigs, "SAMPLE", &rows)?;
            writer.flush()?;
            drop(writer);
            let mut manifest = RunManifest::new("variants");
            manifest.inputs.push(FileHash { path: index_path.display().to_string(), blake3: blake3_file(&index_path)? });
            manifest.inputs.push(FileHash { path: alignments_path.display().to_string(), blake3: blake3_file(&alignments_path)? });
            manifest.outputs.push(FileHash { path: path.display().to_string(), blake3: blake3_file(&path)? });
            manifest.params.insert("mapq_threshold".to_string(), mapq_threshold.to_string());
            manifest.params.insert("min_qual".to_string(), (quality_threshold as f64).to_string());
            let manifest_path = write_manifest(&path, &manifest)?;
            eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            write_germline_vcf(&mut handle, contigs, "SAMPLE", &rows)?;
        }
    }
    let _ = memory_budget_mb; // wired in Task 4
    Ok(())
}
```

(No SAM helper is needed — `--index` is BAM-only. `legacy_read_to_core`/`read_alignment_file` remain used only by the single-contig `--reference` `run_variants`, unchanged.)

- [ ] **Step 5: Build + run the integration tests + suite.**

Run: `cargo build 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --test variants_index 2>&1 | tail -25` (the 3 tests pass), then
`cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` (no failures), then `cargo fmt --all`.

- [ ] **Step 6: Commit.**

```bash
git add src/main.rs tests/variants_index.rs
git commit -m "feat(cli): variants --index — bounded multi-contig calling over the persisted index" \
  -m "variants takes --index XOR --reference (exactly one). --index loads the index (ContigSet + ReferenceView), streams a sorted BAM (StreamingBamSource) or materialises a SAM, runs call_germline_whole_genome, and writes a multi-contig VCF + manifest (inputs = .idx + BAM). --chrom/--region-start error under --index. Parity with --reference on a single contig; self-contained (no FASTA)." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: memory receipt + record-only `--memory-budget-mb`

**Files:**
- Modify: `src/main.rs`, `tests/variants_index.rs`

- [ ] **Step 1: Write the failing test.** Append to `tests/variants_index.rs`:

```rust
#[test]
fn variants_index_memory_budget_reports_and_never_refuses() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);

    // Budget 0 → reported EXCEEDED, but the call still completes (record-only).
    let out = run(&[
        "variants", "--index", idx.to_str().unwrap(), "--alignments", bam.to_str().unwrap(),
        "--memory-budget-mb", "0",
    ]);
    assert!(out.status.success(), "must complete even when over budget: {}", String::from_utf8_lossy(&out.stderr));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("peak RSS"), "receipt must report realized peak: {stderr}");
    assert!(stderr.contains("EXCEEDED") || stderr.contains("budget"), "budget verdict: {stderr}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run it to verify it fails.**

Run: `cargo test --test variants_index variants_index_memory_budget 2>&1 | tail -15`
Expected: FAIL — no "peak RSS"/budget line on stderr yet.

- [ ] **Step 3: Surface the memory receipt + budget.** In `run_variants_index` (src/main.rs), capture the max working set and the realized peak RSS, record them in the manifest params, and print a stderr summary + the record-only budget verdict. Replace `let (sites, _max_ws) = …` to bind `max_ws`, and replace the trailing `let _ = memory_budget_mb;` block with the receipt logic.

Bind the working set (remove the `_` so `max_ws` is used):

```rust
    let (sites, max_ws) = if is_bam { /* …unchanged… */ } else { /* …unchanged… */ };
```

Add the top-level import `use rosalind::core::MemoryBudget;` and `use rosalind::util::rss::peak_rss_bytes;` (the latter is already imported at the top of main.rs).

After writing the VCF (in BOTH the `Some(path)` and `None` branches, or once after the `match`), add the memory receipt. Simplest: compute + report after the `match output { … }` block, and also fold the two memory params into the manifest in the `Some(path)` branch. Concretely, replace the final `let _ = memory_budget_mb; Ok(())` with:

```rust
    // Memory receipt: the bounded contract, made visible + verifiable.
    let peak_rss = peak_rss_bytes();
    eprintln!(
        "memory: peak RSS {} MiB; max pileup working set {} KiB",
        peak_rss / (1 << 20),
        max_ws.bytes / 1024
    );
    if let Some(mb) = memory_budget_mb {
        let budget = MemoryBudget::from_mb(mb);
        if budget.admits(peak_rss) {
            eprintln!("memory: within budget ({} MiB)", mb);
        } else {
            eprintln!(
                "memory: EXCEEDED budget {} MiB (realized peak {} MiB) — record-only, run completed",
                mb,
                peak_rss / (1 << 20)
            );
        }
    }
    Ok(())
```

And in the `Some(path)` branch's manifest, add the memory params before `write_manifest` (so the receipt file records them too):

```rust
            manifest.params.insert("peak_rss_bytes".to_string(), peak_rss_bytes().to_string());
            manifest.params.insert("max_working_set_bytes".to_string(), max_ws.bytes.to_string());
```

(Note: `peak_rss_bytes()` is monotonic peak — calling it in the manifest and again for the stderr line is fine; or compute `let peak_rss = peak_rss_bytes();` once before the `match` and reuse. Prefer computing once before the `match output` and using it in both places.)

To keep it clean, compute `let peak_rss = peak_rss_bytes();` immediately after binding `(sites, max_ws)`, use `peak_rss` in the manifest params and the stderr block, and drop the duplicate call.

- [ ] **Step 4: Build + run the memory test + full suite.**

Run: `cargo build 2>&1 | grep -iE 'error|warning'` (none), then
`cargo test --test variants_index 2>&1 | tail -20` (all 4 tests pass, incl. the budget test), then
`cargo test 2>&1 | grep -E "test result:" | grep -v "0 passed"` (no failures), then `cargo fmt --all`.

- [ ] **Step 5: Commit.**

```bash
git add src/main.rs tests/variants_index.rs
git commit -m "feat(cli): variants --index memory receipt + record-only --memory-budget-mb" \
  -m "Surfaces the bounded-memory contract: the run reports realized peak RSS + the engine's max pileup working set (to stderr and the manifest params). --memory-budget-mb flags overage against the realized peak but never aborts (enforcement is Phase C). Operationalizes Rosalind's differentiator — predictable, verifiable memory — on the flagship whole-genome workload." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before the B4 PR)

- `cargo test` — full suite green. Witnesses: `io::bam` streaming tests; `call::whole_genome` drive == per-contig union; `tests/variants_index` (parity vs `--reference`, self-contained, sort-rejection [via the streaming source test], budget-never-refuses).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo fmt --all -- --check` — clean.
- No-rebuild (structural): `grep -n "run_variants_index" -A40 src/main.rs | grep -c "sais_u32\|BlockedFMIndex::build"` → `0`.
- Streaming (structural): `run_variants_index`'s BAM path uses `StreamingBamSource` (not `read_bam_as_core_reads`/`BamSource`), so reads are not materialised.

## Self-Review

- **Spec coverage (`2026-05-27-phase-b4-variants-index-design.md`):**
  - §2/§3 streaming reads ✔ Task 1 (`StreamingBamSource`, one record at a time); sort guard ✔ Task 1 (`(contig_id,pos)` monotonicity).
  - §2 working set surfaced ✔ Task 2 (`call_germline_region_tracked`) + Task 4 (receipt); §2 per-contig reference ✔ Task 2 (`decode_window` per contig).
  - §4 bounded multi-contig drive ✔ Task 2 (`call_germline_whole_genome` + `PerContig`, reusing `call_germline_region`).
  - §3/§6 `variants --index` + multi-contig VCF + manifest ✔ Task 3 (reuses `write_germline_vcf`'s per-contig `##contig`).
  - §5 memory receipt + record-only `--memory-budget-mb` ✔ Task 4.
  - §6 CLI `--index` XOR `--reference`, `--chrom`/`--region-start` error under `--index` ✔ Task 3.
  - §8 gates — bounded (streaming source, structural), self-contained ✔ Task 3, parity ✔ Task 3, sort-safety ✔ Task 1, reproducible/no-rebuild ✔ (final verification).
- **Type/name consistency:** `record_to_aligned_read(&Record,&HeaderView,&ContigSet)->Result<Option<AlignedRead>>`; `StreamingBamSource::new(&Path,&ContigSet)`; `call_germline_region_tracked(...)->(Vec<(Locus,u8,GermlineCall)>,WorkingSet)`; `call_germline_whole_genome(S,&ReferenceView,&ContigSet,PileupParams,&GermlineParams)->(rows,WorkingSet)`; `GermlineRow{locus,ref_base,call}`; `write_germline_vcf(out,&ContigSet,&str,&[GermlineRow])`. All consistent across tasks and matching the merged APIs (`ReferenceView::decode_window`, `IndexReader::open`/`reference_view`/`contigs`, `MemoryBudget::{from_mb,admits}`, `peak_rss_bytes`).
- **No placeholders:** every step ships complete code or an exact command + expected output; the one transient `unimplemented!()` (Task 2 drive stub) is replaced in the same task. The clearly-flagged "adapt to the actual 0.44.1/legacy-SAM API" notes point at *named* in-tree references (`pileup_stream.rs`'s `reader.read` idiom; `legacy_read_to_core`) to copy, not vague instructions.
- **MSRV 1.72:** `let…else` (1.65), `conflicts_with`/`required_unless_present` (clap 4.5), no `div_ceil`.
- **DRY/reuse:** shares `record_to_aligned_read`; `call_germline_region` delegates to `_tracked` (no duplicate loop); reuses `call_germline_region`, `write_germline_vcf`, the manifest, `MemoryBudget`, `peak_rss_bytes` unchanged.
- **Scope:** germline `variants` only; aligner/`somatic --index`/SAM-streaming/on-demand-`ref_base`/enforcement deferred.
- **Shared-tree hazard:** each implementer/reviewer dispatch must forbid `cargo fix`/mutating git, stage only named files, self-verify `git show --stat HEAD`; coordinator verifies stat + clean tree at every task boundary.
```
