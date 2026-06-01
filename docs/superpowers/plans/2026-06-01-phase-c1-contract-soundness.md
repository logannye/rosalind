# Phase C1 — Working-Set Soundness Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work) or superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the bounded `variants --index` path's reported working set a *true conservative upper bound* on realized peak — the hard correctness prerequisite for `rosalind plan` (C2) — without changing default output.

**Architecture:** Four changes, all on the streaming germline path: (1) fix `PileupEngine::current_working_set()` to count the resident reference + per-read `seq`/`qual` bytes it omits today; (2) add a deterministic `--max-depth` cap to the pileup engine (mechanism only — defaults off in C1); (3) split the VCF writer into header + per-row functions so calls can stream; (4) change `call_germline_whole_genome` to drive a row sink (no genome-wide `Vec`) and decode-then-**move** each contig's reference (kills the persistent `buf`+`Arc` double-copy). Then wire `main.rs::run_variants_index` to write the header once and stream rows through the sink.

**Tech Stack:** Rust 1.72 (MSRV — no `div_ceil`/`is_none_or`), `cargo test`/`fmt`/`build`, the existing `core`/`pileup`/`call`/`io`/`provenance` modules. No new dependencies.

**Spec:** [`docs/superpowers/specs/2026-06-01-phase-c-contract-design.md`](../specs/2026-06-01-phase-c-contract-design.md) §5.

**Note on the C1 soundness proof (deliberate refinement of spec §5.5):** the rigorous proof that the accountant is sound is a **library test** on `call_germline_whole_genome` with an explicit `max_depth: Some(D)` — it asserts the returned `WorkingSet` (a) counts the reference term, (b) is bounded by `reference + D·active_cost`, and (c) is **flat as the read count grows**. This is more robust than a process-RSS assertion. The subprocess/CLI **real-RSS CI gate** lands in C3 (where the `--max-depth` CLI default and the contract suite live).

---

## File Structure

- **Modify** `src/pileup/engine.rs` — `PileupParams.max_depth` field; `SkipCounts.over_max_depth` field + `total()`; `current_working_set()` accounting fix; deterministic cap in `advance_to`. Tests in-file.
- **Modify** `src/io/vcf.rs` — add `write_germline_header` + `write_germline_row`; rewrite `write_germline_vcf` as a (sorting) wrapper over them. Tests in-file.
- **Modify** `src/call/pipeline.rs` — add `call_germline_region_streaming` (sink, returns `WorkingSet`); make `call_germline_region_tracked` a wrapper over it. Tests in-file.
- **Modify** `src/call/mod.rs` — export `call_germline_region_streaming`.
- **Modify** `src/call/whole_genome.rs` — `call_germline_whole_genome` takes a row sink, returns `WorkingSet`; decode-then-move reference. Update its in-file test. Add the soundness test.
- **Modify** `src/main.rs` — `run_variants_index`: write header once, stream rows via the sink, keep the file-path manifest (peak_rss + max_ws) as today.

---

## Task 1: Add the `max_depth` param and `over_max_depth` skip counter (fields only, no behavior)

**Files:**
- Modify: `src/pileup/engine.rs` (`PileupParams`, `PileupParams::default`, `SkipCounts`, `SkipCounts::total`)

- [ ] **Step 1: Write the failing test** — append to the `tests` module in `src/pileup/engine.rs`:

```rust
    #[test]
    fn params_default_max_depth_is_none_and_skipcounts_total_includes_over_max_depth() {
        assert_eq!(PileupParams::default().max_depth, None);
        let s = SkipCounts {
            unmapped: 1,
            wrong_contig: 2,
            secondary: 3,
            supplementary: 4,
            duplicate: 5,
            low_mapq: 6,
            over_max_depth: 7,
        };
        assert_eq!(s.total(), 28);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine::tests::params_default_max_depth 2>&1 | tail -20`
Expected: FAIL — compile error (`PileupParams` has no field `max_depth`; `SkipCounts` has no field `over_max_depth`).

- [ ] **Step 3: Add the fields.** In `PileupParams` (after `skip_duplicate`):

```rust
    /// Skip PCR/optical duplicates (SAM flag 0x400).
    pub skip_duplicate: bool,
    /// Cap on the active read set per position (deterministic downsampling).
    /// `None` = uncapped (default). When `Some(d)`, reads arriving at a position
    /// already covered by `d` active reads are dropped (counted `over_max_depth`).
    pub max_depth: Option<u32>,
```

In `impl Default for PileupParams` (after `skip_duplicate: true,`):

```rust
            skip_duplicate: true,
            max_depth: None,
```

In `SkipCounts` (after `low_mapq`):

```rust
    /// Reads below the MAPQ threshold.
    pub low_mapq: u64,
    /// Reads dropped because the position was already at `max_depth`.
    pub over_max_depth: u64,
```

In `SkipCounts::total` (add the term):

```rust
    pub fn total(&self) -> u64 {
        self.unmapped
            + self.wrong_contig
            + self.secondary
            + self.supplementary
            + self.duplicate
            + self.low_mapq
            + self.over_max_depth
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine::tests::params_default_max_depth 2>&1 | tail -20`
Expected: PASS. (Other `engine` tests that construct `SkipCounts` with all fields — `skip_counts_total_sums_all_reasons` — need the new field; see Step 5.)

- [ ] **Step 5: Fix the existing `SkipCounts` literal.** The test `skip_counts_total_sums_all_reasons` builds a `SkipCounts { ... }` without `over_max_depth` — add `over_max_depth: 0,` to that literal and leave its `assert_eq!(s.total(), 21)` unchanged (0 added).

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine 2>&1 | tail -20`
Expected: PASS (all engine tests).

- [ ] **Step 6: Commit**

```bash
cd ~/rosalind && git add src/pileup/engine.rs && git commit -m "feat(pileup): add max_depth param + over_max_depth skip counter (C1, fields only)"
```

---

## Task 2: Fix `current_working_set()` to count the resident reference + per-read seq/qual

**Files:**
- Modify: `src/pileup/engine.rs` (`current_working_set`)

- [ ] **Step 1: Write the failing test** — append to the `tests` module:

```rust
    #[test]
    fn working_set_counts_reference_and_read_byte_buffers() {
        // One 4-base read fully covering a 10-base reference. After advancing to
        // pos 0 the active set holds that read; the working set must include the
        // reference bytes (10) AND the read's seq+qual buffers (4+4), not just the
        // projection map.
        let reference = b"ACGTACGTAC"; // 10 bytes
        let mut e = engine(vec![mread(0, b"ACGT", false)], reference);
        let first = e.next().expect("a column").expect("ok"); // drives advance_to(0)
        assert_eq!(first.locus.pos.0, 0);
        let ws = e.current_working_set().bytes;
        // reference (10) + map(4*16=64) + seq(4) + qual(4) + per-read(64) + fixed(256)
        // = 10 + 64 + 4 + 4 + 64 + 256 = 402.
        assert_eq!(ws, 402);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine::tests::working_set_counts_reference 2>&1 | tail -20`
Expected: FAIL — `assert_eq!` mismatch (current value omits reference + seq + qual; it reports `4*16 + 64 + 256 = 384`... actually the current formula is `len*16 + 64` per read summed `+ 256` = `64 + 64 + 256 = 384`, no reference, no seq/qual).

- [ ] **Step 3: Implement the corrected accountant.** Replace the body of `current_working_set`:

```rust
    pub fn current_working_set(&self) -> WorkingSet {
        // The decoded reference for this contig is resident in the engine.
        let reference_bytes = self.reference.len() as u64;
        // Each active read holds its projection map (16 B/entry) plus its seq and
        // qual byte buffers; count all three (the map alone is a large undercount,
        // especially for long reads).
        let active_bytes: u64 = self
            .active
            .iter()
            .map(|r| {
                (r.ref_to_read.len() as u64) * 16
                    + r.seq.len() as u64
                    + r.qual.len() as u64
                    + 64
            })
            .sum();
        WorkingSet {
            bytes: reference_bytes + active_bytes + 256,
        }
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine::tests::working_set_counts_reference 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Confirm the existing coverage-bound test still holds.** `working_set_is_bounded_by_coverage_not_input_size` uses a 50,000-byte reference at depth ~1; the new accounting reports ~50,338 bytes, still `< 64*1024` and `fits(1 MiB)`.

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine 2>&1 | tail -20`
Expected: PASS (all engine tests). If `working_set_is_bounded_by_coverage_not_input_size`'s `max_ws < 64 * 1024` now fails, raise that literal bound to `< 128 * 1024` and update the comment to note the reference term is now counted.

- [ ] **Step 6: Commit**

```bash
cd ~/rosalind && git add src/pileup/engine.rs && git commit -m "fix(pileup): current_working_set counts resident reference + read seq/qual (C1 soundness core)"
```

---

## Task 3: Deterministic `max_depth` cap in the pileup engine

**Files:**
- Modify: `src/pileup/engine.rs` (`advance_to`, the `Ordering::Equal` arm)

- [ ] **Step 1: Write the failing test** — append to the `tests` module:

```rust
    #[test]
    fn max_depth_caps_active_set_deterministically() {
        // 5 reads all covering pos 0..4; cap at 2. Only the first 2 (arrival order)
        // are kept; the other 3 are counted over_max_depth. Capped output is
        // identical regardless of input order (SliceSource sorts on construction).
        let reference = b"AAAA";
        let params = PileupParams {
            max_depth: Some(2),
            ..PileupParams::default()
        };
        let run = |reads: Vec<AlignedRead>| -> (Vec<u32>, u64) {
            let mut e = PileupEngine::new(
                SliceSource::new(reads),
                Arc::from(reference.to_vec().into_boxed_slice()),
                0,
                0..4,
                params.clone(),
            );
            let mut depths = Vec::new();
            while let Some(c) = e.next() {
                depths.push(c.unwrap().raw_depth);
            }
            (depths, e.skip_counts().over_max_depth)
        };
        let reads_a = vec![
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
        ];
        let mut reads_b = reads_a.clone();
        reads_b.reverse();
        let (depths_a, over_a) = run(reads_a);
        let (depths_b, over_b) = run(reads_b);
        // Capped: every position sees at most 2 reads.
        assert!(depths_a.iter().all(|&d| d <= 2), "raw depth must be capped at 2");
        assert_eq!(over_a, 3, "3 of 5 reads dropped over max_depth");
        // Deterministic regardless of input order.
        assert_eq!(depths_a, depths_b);
        assert_eq!(over_a, over_b);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine::tests::max_depth_caps 2>&1 | tail -20`
Expected: FAIL — `over_a` is 0 and depths reach 5 (no cap applied yet).

- [ ] **Step 3: Implement the cap.** In `advance_to`, the `std::cmp::Ordering::Equal` arm, insert the cap check between the `read.end() <= pos` guard and `self.ingest(read)`:

```rust
                std::cmp::Ordering::Equal => {
                    if rp > pos {
                        break; // future read on our contig
                    }
                    let read = self.next_read.take().unwrap();
                    if !self.passes_filters(&read) {
                        continue;
                    }
                    if read.end() <= pos {
                        continue; // does not reach the cursor
                    }
                    // Deterministic max-depth cap: once `max_depth` reads already
                    // cover the cursor, drop arrivals (counted) so the active set —
                    // and thus the working set — is bounded by the declared depth.
                    if let Some(max) = self.params.max_depth {
                        if self.active.len() as u32 >= max {
                            self.skips.over_max_depth += 1;
                            continue;
                        }
                    }
                    self.ingest(read);
                }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd ~/rosalind && cargo test -p rosalind --lib pileup::engine 2>&1 | tail -20`
Expected: PASS (new test + all existing engine tests — default `max_depth: None` leaves behavior unchanged).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/pileup/engine.rs && git commit -m "feat(pileup): deterministic max_depth cap on the active set (C1)"
```

---

## Task 4: Split the germline VCF writer into header + row functions

**Files:**
- Modify: `src/io/vcf.rs` (add `write_germline_header`, `write_germline_row`; rewrite `write_germline_vcf` as a wrapper)

- [ ] **Step 1: Write the failing test** — append to the `tests` module in `src/io/vcf.rs`:

```rust
    #[test]
    fn header_then_streamed_rows_equals_batch_write() {
        let r1 = row(0, 100, b'A', het_call());
        let r2 = row(0, 50, b'A', het_call());
        let r3 = row(1, 10, b'A', het_call());
        // Batch writer (sorts internally).
        let batch = render_germline_vcf(&contigs(), "S", &[r1.clone(), r2.clone(), r3.clone()]).unwrap();
        // Streaming: header once, then rows in already-sorted (contig,pos) order.
        let mut buf = Vec::new();
        write_germline_header(&mut buf, &contigs(), "S").unwrap();
        for r in [&r2, &r1, &r3] {
            write_germline_row(&mut buf, &contigs(), r).unwrap();
        }
        let streamed = String::from_utf8(buf).unwrap();
        assert_eq!(streamed, batch, "streamed header+rows must equal the batch write byte-for-byte");
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib io::vcf::tests::header_then_streamed_rows 2>&1 | tail -20`
Expected: FAIL — `write_germline_header` / `write_germline_row` not found.

- [ ] **Step 3: Implement the split.** In `src/io/vcf.rs`, add the two functions and rewrite `write_germline_vcf`. Replace the existing `write_germline_vcf` (lines 48–132) with:

```rust
/// Write the germline VCFv4.2 header (everything up to and including the
/// `#CHROM` line). Pair with `write_germline_row` to stream records.
pub fn write_germline_header<W: Write>(out: &mut W, contigs: &ContigSet, sample: &str) -> io::Result<()> {
    write_fileformat_and_contigs(out, contigs)?;
    writeln!(
        out,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">"#
    )?;
    writeln!(
        out,
        r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">"#
    )?;
    writeln!(
        out,
        r#"##FILTER=<ID=LowQual,Description="QUAL below threshold">"#
    )?;
    writeln!(
        out,
        r#"##FILTER=<ID=LowDepth,Description="Depth below threshold">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )
}

/// Write one germline record line. The caller supplies rows in canonical
/// (contig, pos, ref, alt) order — this does not sort (use `write_germline_vcf`
/// for an unordered batch).
pub fn write_germline_row<W: Write>(out: &mut W, contigs: &ContigSet, r: &GermlineRow) -> io::Result<()> {
    let chrom = contigs
        .by_id(r.locus.contig)
        .map(|c| c.name.as_ref())
        .unwrap_or(".");
    let pos = r.locus.pos.0 as u64 + 1;
    let total = (r.call.ad[0] + r.call.ad[1]).max(1);
    let af = r.call.ad[1] as f64 / total as f64;
    writeln!(
        out,
        "{chrom}\t{pos}\t.\t{ref_b}\t{alt}\t{qual:.1}\t{filt}\tDP={dp};AF={af:.3}\tGT:GQ:DP:AD:PL\t{gt}:{gq}:{dp}:{ad0},{ad1}:{pl0},{pl1},{pl2}",
        ref_b = r.ref_base as char,
        alt = r.call.alt_base as char,
        qual = r.call.qual,
        filt = filter_str(r.call.filter),
        dp = r.call.dp,
        gt = genotype_str(r.call.genotype),
        gq = r.call.gq,
        ad0 = r.call.ad[0],
        ad1 = r.call.ad[1],
        pl0 = r.call.pl[0],
        pl1 = r.call.pl[1],
        pl2 = r.call.pl[2],
    )
}

/// Write a spec-valid germline (single-sample) VCFv4.2 to `out`. Records are
/// emitted in canonical (contig, pos, ref, alt) order regardless of input order.
pub fn write_germline_vcf<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    sample: &str,
    rows: &[GermlineRow],
) -> io::Result<()> {
    write_germline_header(out, contigs, sample)?;
    let mut ordered: Vec<&GermlineRow> = rows.iter().collect();
    ordered.sort_by(|a, b| {
        a.locus
            .cmp(&b.locus)
            .then_with(|| a.ref_base.cmp(&b.ref_base))
            .then_with(|| a.call.alt_base.cmp(&b.call.alt_base))
    });
    for r in ordered {
        write_germline_row(out, contigs, r)?;
    }
    out.flush()
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd ~/rosalind && cargo test -p rosalind --lib io::vcf 2>&1 | tail -25`
Expected: PASS (new test + all existing vcf tests — `write_germline_vcf` output is unchanged: header then sorted rows).

- [ ] **Step 5: Run the golden VCF snapshot to confirm byte-identity**

Run: `cd ~/rosalind && cargo test --test golden_vcf 2>&1 | tail -20`
Expected: PASS (no snapshot drift).

- [ ] **Step 6: Commit**

```bash
cd ~/rosalind && git add src/io/vcf.rs && git commit -m "refactor(vcf): split germline writer into header + row; write_germline_vcf wraps them (C1)"
```

---

## Task 5: Streaming germline caller (`call_germline_region_streaming`)

**Files:**
- Modify: `src/call/pipeline.rs` (add `call_germline_region_streaming`; make `call_germline_region_tracked` a wrapper)
- Modify: `src/call/mod.rs` (export)

- [ ] **Step 1: Write the failing test** — append to the `tests` module in `src/call/pipeline.rs`:

```rust
    #[test]
    fn streaming_emits_same_sites_as_tracked_and_returns_working_set() {
        let reference: Arc<[u8]> = Arc::from(b"AAAA".to_vec().into_boxed_slice());
        let reads = vec![
            read(0, b"ACAA", false),
            read(0, b"ACAA", false),
            read(0, b"AAAA", false),
            read(0, b"ACAA", false),
        ];
        // Reference path: the tracked collector.
        let (collected, _ws) = call_germline_region_tracked(
            SliceSource::new(reads.clone()),
            Arc::clone(&reference),
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();
        // Streaming path: push into a Vec via the sink, capture the working set.
        let mut streamed = Vec::new();
        let ws = call_germline_region_streaming(
            SliceSource::new(reads),
            reference,
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
            &mut |row| {
                streamed.push(row);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(streamed, collected, "streaming sites == tracked sites");
        assert!(ws.bytes > 0, "working set tracked and non-zero");
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib call::pipeline::tests::streaming_emits_same 2>&1 | tail -20`
Expected: FAIL — `call_germline_region_streaming` not found.

- [ ] **Step 3: Implement.** In `src/call/pipeline.rs`, add the streaming primitive and rewrite `call_germline_region_tracked` as a wrapper over it. Replace the current `call_germline_region_tracked` (lines 17–39) with:

```rust
/// Stream germline calls over `region` of `contig` to a sink, returning the
/// maximum pileup-engine working set observed (the bounded-memory signal behind
/// the `variants` receipt and `rosalind plan`). The sink receives each emitted
/// site as `(locus, ref_base, call)` in ascending position order; no
/// genome-wide buffer accumulates. Hom-ref / no-evidence positions are abstained
/// on (the sink is not called for them).
pub fn call_germline_region_streaming<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<WorkingSet, CoreError> {
    let mut engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut max_ws = WorkingSet { bytes: 0 };
    while let Some(column) = engine.next() {
        let column = column?;
        let ws = engine.current_working_set();
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        if let Some(call) = call_germline(&column, germline_params) {
            on_row((column.locus, column.ref_base, call))?;
        }
    }
    Ok(max_ws)
}

/// Like [`call_germline_region`], but also returns the maximum pileup-engine
/// working set observed during the pass. Collects sites into a `Vec` via
/// [`call_germline_region_streaming`].
pub fn call_germline_region_tracked<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<(Vec<(Locus, u8, GermlineCall)>, WorkingSet), CoreError> {
    let mut out = Vec::new();
    let ws = call_germline_region_streaming(
        source,
        reference,
        contig,
        region,
        pileup_params,
        germline_params,
        &mut |row| {
            out.push(row);
            Ok(())
        },
    )?;
    Ok((out, ws))
}
```

In `src/call/mod.rs`, update the `pipeline` re-export line:

```rust
pub use pipeline::{
    call_germline_region, call_germline_region_streaming, call_germline_region_tracked,
    call_somatic_region,
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd ~/rosalind && cargo test -p rosalind --lib call::pipeline 2>&1 | tail -20`
Expected: PASS (new test + the two existing pipeline tests, which go through the unchanged `call_germline_region` → `_tracked` → `_streaming`).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/call/pipeline.rs src/call/mod.rs && git commit -m "feat(call): call_germline_region_streaming (sink, returns WorkingSet); tracked wraps it (C1)"
```

---

## Task 6: `call_germline_whole_genome` drives a row sink + decode-then-move

**Files:**
- Modify: `src/call/whole_genome.rs` (signature → sink + `WorkingSet`; decode-then-move; update in-file test; add the soundness test)

- [ ] **Step 1: Write the failing tests** — replace the in-file `whole_genome_equals_per_contig_calls` test body's call site and add the soundness test. First, update the existing test to use the sink (replace the `call_germline_whole_genome(...)` call and the `(rows, ws)` binding, lines 148–155, with):

```rust
        let mut rows: Vec<(crate::core::Locus, u8, GermlineCall)> = Vec::new();
        let ws = call_germline_whole_genome(
            SliceSource::new(reads.clone()),
            &rv,
            contigs,
            pp.clone(),
            &gp,
            &mut |row| {
                rows.push(row);
                Ok(())
            },
        )
        .unwrap();
```

Then append a new soundness test to the `tests` module:

```rust
    #[test]
    fn working_set_is_bounded_by_reference_and_capped_depth_not_read_count() {
        // One small contig; pour in increasing numbers of reads at the SAME few
        // positions with a depth cap. The returned working set must (a) count the
        // reference, (b) stay bounded by reference + capped active, and (c) NOT
        // grow with the number of input reads.
        let idx_path = tmp("bounded");
        let refseq = vec![b'A'; 2000];
        let index = GenomeIndex::from_named_sequences(&[(
            "chr1".to_string(),
            refseq.clone(),
        )])
        .unwrap();
        IndexWriter::create(&idx_path)
            .unwrap()
            .write_genome_index(&index)
            .unwrap();
        let loaded = IndexReader::open(&idx_path).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();

        let params = PileupParams {
            max_depth: Some(8),
            ..PileupParams::default()
        };
        let gp = GermlineParams::default();

        let run = |n: usize| -> u64 {
            // n reads, each 100bp, all starting at pos 0 (depth would be n without
            // the cap; capped at 8).
            let reads: Vec<AlignedRead> = (0..n)
                .map(|_| read_at(0, 0, &vec![b'C'; 100]))
                .collect();
            let mut sink_calls = 0u64;
            let ws = call_germline_whole_genome(
                SliceSource::new(reads),
                &rv,
                contigs,
                params.clone(),
                &gp,
                &mut |_row| {
                    sink_calls += 1;
                    Ok(())
                },
            )
            .unwrap();
            ws.bytes
        };

        let ws_small = run(20);
        let ws_large = run(2000);
        // (a) reference is counted: bound exceeds the 2000-byte reference.
        assert!(ws_small > 2000, "working set must include the reference bytes");
        // (c) flat in read count: 100x more reads, same bounded working set.
        assert_eq!(
            ws_small, ws_large,
            "working set must not grow with the number of input reads"
        );
        // (b) bounded by reference + capped active: ref(2000) + 8 reads *
        // (map 100*16 + seq 100 + qual 100 + 64) + 256, generously bounded.
        let bound = 2000 + 8 * (100 * 16 + 100 + 100 + 64) + 256;
        assert!(
            ws_small <= bound,
            "working set {ws_small} exceeded the analytic bound {bound}"
        );

        let _ = std::fs::remove_dir_all(idx_path.parent().unwrap());
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd ~/rosalind && cargo test -p rosalind --lib call::whole_genome 2>&1 | tail -25`
Expected: FAIL — `call_germline_whole_genome` arity/return mismatch (no sink param; returns a tuple).

- [ ] **Step 3: Implement the sink + decode-then-move.** Replace the body of `call_germline_whole_genome` (the function at lines 47–87) with:

```rust
pub fn call_germline_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<WorkingSet, CoreError> {
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut peeked: Option<AlignedRead> = None;

    for c in contigs.iter() {
        // Decode this contig's reference into a fresh Vec and MOVE it into the
        // Arc — no persistent second copy (the steady-state reference resident is
        // one contig, not two). Peak = the largest contig; bounded.
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        let mut decoded = Vec::new();
        ref_view.decode_window(start, end, &mut decoded);
        let reference: Arc<[u8]> = Arc::from(decoded);

        let per = PerContig {
            source: &mut source,
            contig: c.id,
            peeked: &mut peeked,
        };
        let region: Range<u32> = 0..c.length;
        let ws = call_germline_region_streaming(
            per,
            reference,
            c.id,
            region,
            pileup_params.clone(),
            germline_params,
            on_row,
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
    }

    Ok(max_ws)
}
```

Update the imports at the top of `whole_genome.rs`: change
`use crate::call::{call_germline_region_tracked, GermlineCall, GermlineParams};`
to
`use crate::call::{call_germline_region_streaming, GermlineCall, GermlineParams};`
(the doc-comment on the function should also drop the "accumulated rows" wording — it now returns only the max working set and streams rows to `on_row`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd ~/rosalind && cargo test -p rosalind --lib call::whole_genome 2>&1 | tail -25`
Expected: PASS (the updated equivalence test + the new soundness test).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/call/whole_genome.rs && git commit -m "feat(call): whole-genome drive streams rows to a sink + decode-then-move reference (C1)"
```

---

## Task 7: Wire `main.rs::run_variants_index` to stream via the sink

**Files:**
- Modify: `src/main.rs` (`run_variants_index`, lines ~1001–1102: imports + the call + the two output branches)

- [ ] **Step 1: Update the `use` line.** In `run_variants_index`, change:

```rust
    use rosalind::io::vcf::{write_germline_vcf, GermlineRow};
```
to:
```rust
    use rosalind::core::CoreError;
    use rosalind::io::vcf::{write_germline_header, write_germline_row, GermlineRow};
```

- [ ] **Step 2: Replace the call + output section.** Replace the block from `let (sites, max_ws) =` (line ~1044) through the end of the `match output { … }` block (line ~1102) with the streaming form:

```rust
    let source = StreamingBamSource::new(&alignments_path, contigs)
        .map_err(|e| anyhow!("failed to open BAM {}: {e}", alignments_path.display()))?;

    // Stream calls straight to the VCF writer (header once, then one row per
    // emitted call) so no genome-wide row buffer accumulates. The returned
    // WorkingSet is the high-water (reference + active set), captured per contig.
    let max_ws = match &output {
        Some(path) => {
            let file = File::create(path)
                .with_context(|| format!("failed to create VCF file {}", path.display()))?;
            let mut writer = io::BufWriter::new(file);
            write_germline_header(&mut writer, contigs, "SAMPLE")?;
            let ws = call_germline_whole_genome(
                source,
                &ref_view,
                contigs,
                pileup_params,
                &germline_params,
                &mut |(locus, ref_base, call)| {
                    write_germline_row(&mut writer, contigs, &GermlineRow { locus, ref_base, call })
                        .map_err(CoreError::from)
                },
            )
            .map_err(|e| anyhow!("variant calling failed: {e}"))?;
            writer.flush()?;
            ws
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            write_germline_header(&mut handle, contigs, "SAMPLE")?;
            let ws = call_germline_whole_genome(
                source,
                &ref_view,
                contigs,
                pileup_params,
                &germline_params,
                &mut |(locus, ref_base, call)| {
                    write_germline_row(&mut handle, contigs, &GermlineRow { locus, ref_base, call })
                        .map_err(CoreError::from)
                },
            )
            .map_err(|e| anyhow!("variant calling failed: {e}"))?;
            handle.flush()?;
            ws
        }
    };
    // Realized peak (monotonic high-water mark) captured after the calling pass.
    let peak_rss = peak_rss_bytes();

    // Reproducibility + memory receipt (file output only; stdout receipt is C3).
    if let Some(path) = &output {
        let mut manifest = RunManifest::new("variants");
        manifest.inputs.push(FileHash {
            path: index_path.display().to_string(),
            blake3: blake3_file(&index_path)?,
        });
        manifest.inputs.push(FileHash {
            path: alignments_path.display().to_string(),
            blake3: blake3_file(&alignments_path)?,
        });
        manifest.outputs.push(FileHash {
            path: path.display().to_string(),
            blake3: blake3_file(path)?,
        });
        manifest
            .params
            .insert("mapq_threshold".to_string(), mapq_threshold.to_string());
        manifest.params.insert(
            "min_qual".to_string(),
            (quality_threshold as f64).to_string(),
        );
        manifest
            .params
            .insert("peak_rss_bytes".to_string(), peak_rss.to_string());
        manifest.params.insert(
            "max_working_set_bytes".to_string(),
            max_ws.bytes.to_string(),
        );
        let manifest_path = write_manifest(path, &manifest)?;
        eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
    }
```

(The `let rows: Vec<GermlineRow> = sites …collect();` block that previously sat between the call and the `match` is **removed** — rows are streamed, not collected.)

The existing memory-receipt + record-only budget block (`eprintln!("memory: peak RSS …")` … through the `record-only, run completed` branch, lines ~1103–1120) is **unchanged** and remains after this block.

- [ ] **Step 3: Build to verify it compiles**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -25`
Expected: success, 0 warnings. (If `write_germline_vcf`/`GermlineRow` become unused elsewhere in `main.rs`, the build will warn — confirm `run_variants` (the single-contig `--reference` path, line ~857) still uses `write_germline_vcf`; it does, so the import stays valid there.)

- [ ] **Step 4: Run the whole-genome CLI gate**

Run: `cd ~/rosalind && cargo test --test variants_index 2>&1 | tail -25`
Expected: PASS (the multi-contig `variants --index` output is byte-identical — header then sorted rows, same as before).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/main.rs && git commit -m "feat(cli): variants --index streams rows to the VCF writer via the sink (C1)"
```

---

## Task 8: Full-suite verification + format/warning gates

**Files:** none (verification only)

- [ ] **Step 1: Format check**

Run: `cd ~/rosalind && cargo fmt --all -- --check 2>&1 | tail -20`
Expected: no output (clean). If it reports diffs, run `cargo fmt --all` and re-commit the touched files.

- [ ] **Step 2: Zero-warning build (debug + release)**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -10 && cargo build --release 2>&1 | tail -10`
Expected: both finish with 0 warnings.

- [ ] **Step 3: Full test suite**

Run: `cd ~/rosalind && cargo test 2>&1 | tail -40`
Expected: all binaries pass — in particular `determinism`, `golden_vcf`, `variants_index`, `space_bounds`, `pileup_stream`, and the in-lib `pileup`/`call`/`io` tests.

- [ ] **Step 4: Commit any formatting fixups** (only if Step 1 required changes)

```bash
cd ~/rosalind && git add -A && git commit -m "style: rustfmt fixups (C1)"
```

---

## Self-Review notes (filled during writing)

- **Spec coverage (§5):** §5.1 accountant → Task 2; §5.2 cap → Tasks 1+3; §5.3 incremental writer → Task 4; §5.4 sink + decode-then-move → Tasks 5+6+7; §5.5 proof → Task 6 soundness test (library-level; the process-RSS CI gate is deliberately deferred to C3, noted at the top).
- **Type consistency:** the sink type `&mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>` is identical in `call_germline_region_streaming` (Task 5), `call_germline_whole_genome` (Task 6), and both `main.rs` call sites (Task 7). `write_germline_row(out, contigs, row)` arity matches across Task 4 (definition) and Tasks 6/7 callers. `WorkingSet`/`CoreError` are the existing `core` types.
- **No behavior change by default:** `max_depth` defaults to `None` (Task 1); `current_working_set` only grows the reported number (Task 2); the VCF wrapper still sorts (Task 4); streamed rows arrive pre-sorted so output is byte-identical (Tasks 6/7, gated by `golden_vcf` + `variants_index`).
