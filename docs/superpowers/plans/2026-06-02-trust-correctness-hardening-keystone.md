# Trust + Correctness Hardening Keystone — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:executing-plans (inline) to implement
> this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make the shipped memory contract's correctness + trust guarantees actually hold and be
tested: stop the depth cap from silently dropping real variants, keep output byte-identical, surface
when the cap engages, enforce `--max-read-len`, and regression-guard the exit-4 breach path + a CI
memory gate.

**Architecture:** All changes are on the call/contract path. The depth cap becomes a fixed-seed
min-hash bounded reservoir (unbiased, deterministic, bounded). Obs are emitted in a canonical total
order so the f64 likelihood sum is order-independent. `SkipCounts` is threaded out to the CLI +
receipt. `--max-read-len` becomes an ingest precondition under `--enforce`. A post-run-only RSS test
seam makes the exit-4 branch testable; CI exercises the bounded `--index` path.

**Tech stack:** Rust (MSRV 1.72), `cargo test`/`fmt`, GitHub Actions. Spec:
`docs/superpowers/specs/2026-06-02-trust-correctness-hardening-keystone-design.md`.

**Conventions:** TDD (failing test first), per-task commit, 0 warnings (debug + release),
`cargo fmt --all -- --check` clean, full `cargo test` green at each task boundary. Branch:
`rosalind/trust-correctness-keystone`.

---

## File map

- `src/pileup/engine.rs` — Item A (reservoir + `read_priority` + `ActiveRead.priority`), Item B
  (canonical obs in `build_column`), Item D (`PileupParams.max_read_len` + ingest check).
- `src/call/pipeline.rs`, `src/call/whole_genome.rs` — Item C (return `SkipCounts`).
- `src/core/error.rs` — Item D (`ReadExceedsDeclaredLength` variant).
- `src/main.rs` — Item C (stderr summary + receipt fields), Item D (pass `max_read_len` under
  `--enforce`, help-text fixes), Item E (force-RSS seam at the realized-peak capture).
- `src/call/plan.rs` — Item D (fix the contradictory comment).
- `tests/plan_enforce.rs` — Items C/D/E integration tests.
- `.github/workflows/ci.yml` — Item F (contract gate on `--index`).

---

## Task 1: Canonical obs order in `build_column` (Item B)

**Files:** Modify `src/pileup/engine.rs` (`build_column`, ~262-302) + a unit test.

Rationale: Item A reorders the active set; the germline f64 likelihood sum
(`src/call/germline.rs:48-58`) is over `column.obs` in order and is non-associative. Emitting obs in
a canonical total order makes downstream QUAL/PL byte-identical regardless of active-set order.

- [ ] **Step 1: Write the failing test** (append to `engine.rs` `#[cfg(test)] mod tests`)

```rust
    #[test]
    fn obs_are_emitted_in_canonical_order() {
        // A column with mixed alleles/quals must emit obs sorted by
        // (allele, base_qual, mapq, reverse) — independent of read arrival order,
        // so the downstream f64 likelihood sum is order-stable.
        let reference = b"AAAA";
        // Three reads covering pos 0 with different alleles at offset 0:
        //  C (allele 1), A (allele 0, ref), G (allele 2). Arrival order C,A,G.
        let reads = vec![
            mread(0, b"C", false),
            mread(0, b"A", false),
            mread(0, b"G", false),
        ];
        let cols = columns(engine(reads, reference));
        let at0 = cols.iter().find(|c| c.locus.pos.0 == 0).unwrap();
        let alleles: Vec<u8> = at0.obs.iter().map(|o| o.allele).collect();
        // Canonical: sorted ascending by allele (0,1,2) regardless of C,A,G arrival.
        assert_eq!(alleles, vec![0, 1, 2]);
    }
```

- [ ] **Step 2: Run it — expect FAIL** (current order is arrival order C,A,G → `[1,0,2]`)

Run: `cargo test -p rosalind --lib pileup::engine::tests::obs_are_emitted_in_canonical_order`
Expected: assertion failure `[1, 0, 2]` vs `[0, 1, 2]`.

- [ ] **Step 3: Implement — sort obs canonically in `build_column`**

In `build_column`, after the `for r in &self.active { … obs.push(Obs { … }) … }` loop and before
constructing `PileupColumn`, insert:

```rust
        // Emit observations in a canonical total order so the downstream diploid
        // log-likelihood accumulation (germline.rs) is order-independent and the
        // VCF stays byte-identical regardless of active-set internal order.
        obs.sort_by(|a, b| {
            (a.allele, a.base_qual, a.mapq, a.reverse as u8).cmp(&(
                b.allele,
                b.base_qual,
                b.mapq,
                b.reverse as u8,
            ))
        });
```

- [ ] **Step 4: Run the new test + the golden/determinism suite**

Run: `cargo test -p rosalind --lib pileup::engine::tests::obs_are_emitted_in_canonical_order`
Expected: PASS.
Run: `cargo test --test golden_vcf --test determinism`
Expected: PASS (uniform-quality fixtures ⇒ reordering identical terms doesn't change QUAL/PL). If
`golden_vcf` fails, inspect the diff — it must be at most a ±1 rounding artifact from reordered
summation; regenerate with `ROSALIND_UPDATE_SNAPSHOTS=1 cargo test --test golden_vcf` and note it in
the commit. Any larger/structural change is a bug — stop and investigate.

- [ ] **Step 5: Commit**

```bash
git add src/pileup/engine.rs tests/snapshots 2>/dev/null; git add -A
git commit -m "fix(pileup): emit obs in canonical order (order-stable QUAL/PL)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Unbiased hash-priority depth-cap reservoir (Item A)

**Files:** Modify `src/pileup/engine.rs` (`ActiveRead`, `ingest`, `advance_to`, a new `read_priority`
free fn) + tests.

- [ ] **Step 1: Write the failing regression test** (append to `engine.rs` tests)

```rust
    #[test]
    fn unbiased_cap_keeps_downstream_starting_reads() {
        // The biased (leftmost-arrival) cap dropped reads that START at a deep
        // variant, zeroing the alt allele. The min-hash reservoir keeps a
        // start-position-independent sample, so a variant carried only by reads
        // that begin AT the variant column survives. Reads must be DISTINCT
        // (real reads are) — identical content shares one hash and degenerates.
        let reference = vec![b'A'; 64];
        let v = 32u32; // variant column
        let cap = 30u32;
        let mut reads = Vec::new();
        // 30 ref reads starting upstream (pos 0), distinct lengths -> distinct
        // hashes, all carrying A (ref) at v.
        for i in 0..30u32 {
            let len = 33 + i as usize; // covers v=32 (len>32); 33..62 < 64
            reads.push(mread(0, &vec![b'A'; len], false));
        }
        // 30 alt reads starting AT v, distinct lengths, carrying C at v (offset 0).
        for j in 0..30u32 {
            let len = 1 + j as usize; // start v, covers v
            let mut seq = vec![b'A'; len];
            seq[0] = b'C';
            reads.push(mread(v, &seq, false));
        }
        let params = PileupParams {
            max_depth: Some(cap),
            ..PileupParams::default()
        };
        let mut e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.into_boxed_slice()),
            0,
            0..64,
            params,
        );
        let mut at_v = None;
        while let Some(c) = e.next() {
            let col = c.unwrap();
            if col.locus.pos.0 == v {
                at_v = Some(col);
            }
        }
        let counts = at_v.expect("a column at v").allele_counts();
        // Both alleles present: ref (A=0) AND alt (C=1) — the biased cap gave alt==0.
        assert!(counts[1] > 0, "alt allele must survive the cap: {counts:?}");
        assert!(counts[0] > 0, "ref allele must remain too: {counts:?}");
        // Capped: total observed at v <= cap.
        assert!(counts.iter().sum::<u32>() <= cap, "depth must be capped");
    }
```

- [ ] **Step 2: Run it — expect FAIL** (current biased cap: `counts[1] == 0`)

Run: `cargo test -p rosalind --lib pileup::engine::tests::unbiased_cap_keeps_downstream_starting_reads`
Expected: assertion `alt allele must survive the cap` fails (alt == 0).

- [ ] **Step 3: Implement the reservoir**

(3a) Add a fixed-seed FNV-1a `read_priority` free fn near the top of `engine.rs` (after imports):

```rust
/// Fixed-seed FNV-1a-64 over a read's identity (position, end, flags, CIGAR, seq,
/// qual). Deterministic across runs/processes (NOT `DefaultHasher`, whose seed is
/// per-process). Uncorrelated with the allele at any single column — the basis of
/// the unbiased depth-cap sample.
fn read_priority(read: &AlignedRead) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    #[inline]
    fn fold(mut h: u64, bytes: &[u8]) -> u64 {
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(FNV_PRIME);
        }
        h
    }
    let mut h = FNV_OFFSET;
    h = fold(h, &read.pos.0.to_le_bytes());
    h = fold(h, &read.end().to_le_bytes());
    h = fold(h, &read.flags.0.to_le_bytes());
    for op in &read.cigar {
        h = fold(h, &[op.kind as u8]);
        h = fold(h, &op.len.to_le_bytes());
    }
    h = fold(h, &read.seq);
    h = fold(h, &read.qual);
    h
}
```

(3b) Add `priority: u64` to `ActiveRead` (after `reverse: bool`):

```rust
    /// Fixed-seed content hash — the reservoir admission/eviction key.
    priority: u64,
```

(3c) Change `ingest` to take the precomputed priority:

```rust
    /// Precompute a read's reference→read-offset map and add it to the active set.
    fn ingest(&mut self, read: AlignedRead, priority: u64) {
        let end = read.end();
        let mut ref_to_read = HashMap::new();
        for rb in read.projected_bases() {
            ref_to_read.insert(rb.ref_pos, rb.read_offset);
        }
        self.active.push(ActiveRead {
            end,
            ref_to_read,
            seq: Arc::clone(&read.seq),
            qual: Arc::clone(&read.qual),
            mapq: read.mapq,
            reverse: read.flags.is_reverse(),
            priority,
        });
    }
```

(3d) In `advance_to`, replace the cap block (the `if let Some(max) = self.params.max_depth { if
self.active.len() as u32 >= max { self.skips.over_max_depth += 1; continue; } }` and the following
`self.ingest(read);`) with:

```rust
                    let prio = read_priority(&read);
                    if let Some(max) = self.params.max_depth {
                        if self.active.len() as u32 >= max {
                            // Unbiased min-hash reservoir: keep the `max`
                            // smallest-priority reads covering the cursor. Evict the
                            // greatest-priority resident iff the newcomer ranks below
                            // it; otherwise refuse. Either way the cap removed one
                            // read (counted). Priority is a fixed-seed content hash,
                            // so selection is independent of start position / allele.
                            let mut worst = 0usize;
                            for i in 1..self.active.len() {
                                if self.active[i].priority > self.active[worst].priority {
                                    worst = i;
                                }
                            }
                            self.skips.over_max_depth += 1;
                            if prio < self.active[worst].priority {
                                self.active.swap_remove(worst);
                                self.ingest(read, prio);
                            }
                            continue;
                        }
                    }
                    self.ingest(read, prio);
```

- [ ] **Step 4: Run the new test + the existing cap/bounded tests**

Run: `cargo test -p rosalind --lib pileup::engine::tests`
Expected: PASS, including `unbiased_cap_keeps_downstream_starting_reads`,
`max_depth_caps_active_set_deterministically` (5 identical reads share a hash ⇒ ties ⇒ no eviction ⇒
first-2 kept, over=3 — unchanged), and `working_set_is_bounded_by_coverage_not_input_size`.

- [ ] **Step 5: Full suite + warnings + fmt**

Run: `cargo test` ; `cargo build --release 2>&1 | grep -c warning: || true` (expect 0) ;
`cargo fmt --all -- --check`
Expected: all green, 0 warnings.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "fix(pileup): unbiased min-hash depth-cap reservoir (no silent variant drops)

Replaces leftmost-arrival truncation with a fixed-seed content-hash bounded
reservoir + eviction. Keeps a start-position/allele-independent sample at the
declared cap; a deep het variant carried by reads starting at the variant
survives. Deterministic and bounded at max_depth. Regression test included.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Surface `SkipCounts` on the `--index` path + receipt (Item C)

**Files:** Modify `src/pileup/engine.rs` (`SkipCounts::accumulate`), `src/call/pipeline.rs`,
`src/call/whole_genome.rs`, `src/main.rs` + tests.

- [ ] **Step 1: Add a `SkipCounts` accumulator + test** (engine.rs)

In `impl SkipCounts` (after `total`):

```rust
    /// Field-wise sum (across per-contig engines on the whole-genome path).
    pub fn accumulate(&mut self, other: &SkipCounts) {
        self.unmapped += other.unmapped;
        self.wrong_contig += other.wrong_contig;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.duplicate += other.duplicate;
        self.low_mapq += other.low_mapq;
        self.over_max_depth += other.over_max_depth;
    }
```

Test (engine.rs tests):

```rust
    #[test]
    fn skip_counts_accumulate_sums_fieldwise() {
        let mut a = SkipCounts {
            over_max_depth: 2,
            low_mapq: 1,
            ..SkipCounts::default()
        };
        let b = SkipCounts {
            over_max_depth: 3,
            duplicate: 4,
            ..SkipCounts::default()
        };
        a.accumulate(&b);
        assert_eq!(a.over_max_depth, 5);
        assert_eq!(a.low_mapq, 1);
        assert_eq!(a.duplicate, 4);
    }
```

Run: `cargo test -p rosalind --lib pileup::engine::tests::skip_counts_accumulate_sums_fieldwise` →
PASS (pure addition).

- [ ] **Step 2: Change `call_germline_region_streaming` to return `(WorkingSet, SkipCounts)`**

In `src/call/pipeline.rs`, change the signature + body:

```rust
pub fn call_germline_region_streaming<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
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
    Ok((max_ws, engine.skip_counts()))
}
```

Add `SkipCounts` to the `use crate::pileup::{…}` import line. Update `call_germline_region_tracked`
(it calls `_streaming`): destructure `(ws, _skips)` and keep returning `(Vec, WorkingSet)`:

```rust
    let (ws, _skips) = call_germline_region_streaming( … same args … )?;
    Ok((out, ws))
```

(`call_germline_region` is unchanged — it calls `_tracked`.)

- [ ] **Step 3: Change `call_germline_whole_genome` to return `(WorkingSet, SkipCounts)`**

In `src/call/whole_genome.rs`: add `SkipCounts` to the `use crate::pileup::{…}` import; change the
return type and accumulate per contig:

```rust
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    let mut peeked: Option<AlignedRead> = None;

    for c in contigs.iter() {
        // … decode reference, build `per`, region (unchanged) …
        let (ws, contig_skips) = call_germline_region_streaming(
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
        skips.accumulate(&contig_skips);
    }

    Ok((max_ws, skips))
}
```

Update the in-file test `whole_genome_equals_per_contig_calls`: it binds `let ws =
call_germline_whole_genome(…)` → change to `let (ws, _skips) = call_germline_whole_genome(…)`.
Update `working_set_is_bounded_by_reference_and_capped_depth_not_read_count`: the `run` closure ends
`.unwrap().bytes` → change to `.unwrap().0.bytes`.

- [ ] **Step 4: Update `main.rs` call sites + stderr summary + receipt fields**

In `run_variants_index` (`src/main.rs:1294-1350`), both match arms bind `let ws =
call_germline_whole_genome(…)?;` and `ws` at the arm tail. Change the outer binding to capture both:
make each arm return `(ws, skips)` and bind `let (max_ws, skips) = match &output { … };`.
Concretely, in each arm replace `let ws = call_germline_whole_genome(…)…?; writer.flush()?; ws`
with `let (ws, sk) = call_germline_whole_genome(…)…?; writer.flush()?; (ws, sk)` and change the
outer `let max_ws = match &output {` to `let (max_ws, skips) = match &output {`.

After the realized-peak line (`eprintln!("memory: peak RSS …")`, ~1430-1434), add the skip summary:

```rust
    if skips.over_max_depth > 0 {
        eprintln!(
            "pileup: dropped {} reads at --max-depth {} (deep-site downsampling — \
             calls at those sites use a bounded unbiased sample)",
            skips.over_max_depth, max_depth
        );
    }
    let other_skipped = skips.total() - skips.over_max_depth;
    if other_skipped > 0 {
        eprintln!(
            "pileup: skipped {other_skipped} reads by filter (unmapped/wrong-contig/\
             secondary/supplementary/duplicate/low-mapq)"
        );
    }
```

In the receipt block (after `max_working_set_bytes`, ~1412), add:

```rust
        manifest
            .params
            .insert("over_max_depth".to_string(), skips.over_max_depth.to_string());
        manifest
            .params
            .insert("reads_skipped_total".to_string(), skips.total().to_string());
```

- [ ] **Step 5: Add an integration test** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn receipt_records_skip_counts() {
    // The standard fixture (5 reads, max depth way below default 1000) drops
    // nothing → over_max_depth 0. A tight --max-depth 1 forces drops → > 0.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let manifest = dir.join("run.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--max-depth", "1", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");
    let json = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        json.contains("\"over_max_depth\":") && json.contains("\"reads_skipped_total\":"),
        "manifest missing skip fields: {json}"
    );
    // At pos 0 the fixture stacks >1 read; --max-depth 1 must drop at least one.
    let m = rosalind::provenance::RunManifest::from_canonical_json(&json).unwrap();
    let over: u64 = m.params.get("over_max_depth").unwrap().parse().unwrap();
    assert!(over > 0, "tight cap should drop reads, got {over}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 6: Run, verify, commit**

Run: `cargo test` ; `cargo build --release 2>&1 | grep -c warning: || true` ;
`cargo fmt --all -- --check`
Expected: all green, 0 warnings. (Also confirm `stdout_run_persists_a_self_describing_receipt` and
`estimator_upper_bounds_the_realized_working_set` still pass — they don't read the new fields.)

```bash
git add -A
git commit -m "feat(variants): surface depth-cap + filter skip counts (stderr + receipt)

Thread SkipCounts out of the streaming + whole-genome callers; print a summary
when the cap engages; record over_max_depth + reads_skipped_total in the
manifest so the receipt is honest about whether downsampling was load-bearing.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Enforce `--max-read-len` at ingest under `--enforce` (Item D)

**Files:** `src/core/error.rs`, `src/pileup/engine.rs` (`PileupParams` + `advance_to`),
`src/main.rs`, `src/call/plan.rs` (comment) + tests.

- [ ] **Step 1: Add the error variant** (`src/core/error.rs`, inside `enum CoreError`)

```rust
    /// A read longer than the declared `--max-read-len` voids the predicted
    /// memory envelope (raised only under `--enforce`).
    #[error("read length {len} exceeds declared --max-read-len {declared}; raise --max-read-len or drop --enforce")]
    ReadExceedsDeclaredLength {
        /// The offending read's sequence length.
        len: u32,
        /// The declared cap.
        declared: u32,
    },
```

- [ ] **Step 2: Add `max_read_len` to `PileupParams` + failing test** (`src/pileup/engine.rs`)

Add field (after `max_depth`):

```rust
    /// When `Some(m)` (set only under `--enforce`), a read whose `seq.len()`
    /// exceeds `m` aborts the run (the predicted envelope assumes `<= m`).
    /// `None` = no check (default).
    pub max_read_len: Option<u32>,
```

Add `max_read_len: None` to the `Default` impl.

Failing test (engine.rs tests):

```rust
    #[test]
    fn read_exceeding_declared_max_read_len_errors() {
        let reference = b"AAAAAAAA";
        let params = PileupParams {
            max_read_len: Some(4),
            ..PileupParams::default()
        };
        let mut e = PileupEngine::new(
            SliceSource::new(vec![mread(0, b"CCCCCC", false)]), // len 6 > 4
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..8,
            params,
        );
        let err = loop {
            match e.next() {
                Some(Ok(_)) => continue,
                Some(Err(err)) => break err,
                None => panic!("expected an error, got clean end"),
            }
        };
        assert!(matches!(
            err,
            crate::core::CoreError::ReadExceedsDeclaredLength { len: 6, declared: 4 }
        ));
    }
```

Run: `cargo test -p rosalind --lib pileup::engine::tests::read_exceeding_declared_max_read_len_errors`
→ FAIL (no check yet; `CoreError` variant must compile first — add Step 1 before running).

- [ ] **Step 3: Implement the ingest check** (`advance_to`, just before `let prio = read_priority`)

```rust
                    if let Some(maxlen) = self.params.max_read_len {
                        if read.seq.len() as u32 > maxlen {
                            return Err(CoreError::ReadExceedsDeclaredLength {
                                len: read.seq.len() as u32,
                                declared: maxlen,
                            });
                        }
                    }
```

Run the test → PASS. Confirm `long_read_piles_up_every_matched_base` (max_read_len `None`) still
passes.

- [ ] **Step 4: Wire it under `--enforce` in `main.rs` + fix help/comment text**

(4a) `run_variants_index` builds `pileup_params`. Find where `PileupParams { … max_depth … }` is
constructed (the `max_depth == 0 → None` logic) and set `max_read_len`:

```rust
        max_read_len: if enforce { Some(max_read_len) } else { None },
```

(4b) Fix the stale `--memory-budget-mb` help (`src/main.rs:96-97`):

```rust
        /// Declared memory budget (MiB) for the run — records a plan/peak line.
        /// With `--enforce`, it is honored (exit 3 refuse / exit 4 breach). (`--index` path.)
```

(4c) Update the `--max-read-len` help (`src/main.rs:104-106`):

```rust
        /// Max read length assumed by the pre-run `--enforce` estimate AND
        /// enforced at ingest under `--enforce` (a longer read aborts the run).
```

(4d) Fix the contradictory comment in `src/call/plan.rs:19-20` (it claims `max_read_len` is
"enforced at runtime"): reword to state the estimate assumes `<= --max-read-len`, which `--enforce`
now checks at ingest (engine.rs).

- [ ] **Step 5: Integration test** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn enforce_aborts_on_a_read_longer_than_declared_max_read_len() {
    // Standard fixture has 16 bp reads. Declare --max-read-len 8 under --enforce
    // (with a generous budget so we hit the ingest check, not the exit-3 gate).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args([
            "--memory-budget-mb", "4096",
            "--max-read-len", "8",
            "--max-depth", "1000",
            "--enforce",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "over-long read must abort: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("exceeds declared --max-read-len"),
        "missing clear message: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

Note: confirm the budget passes the exit-3 gate at `--max-read-len 8` (predicted uses 8 → smaller →
fits 4096), so the run proceeds to ingest and hits the abort. If exit-3 intervenes, raise the budget.

- [ ] **Step 6: Run, verify, commit**

Run: `cargo test` ; `cargo build --release 2>&1 | grep -c warning: || true` ;
`cargo fmt --all -- --check` → all green, 0 warnings.

```bash
git add -A
git commit -m "feat(enforce): check --max-read-len at ingest (close the silent-OOM path)

Under --enforce, a read whose seq.len() exceeds the declared --max-read-len
aborts loudly (CoreError::ReadExceedsDeclaredLength) instead of silently voiding
the predicted envelope. Non-enforced runs are unaffected. Fixes the contradictory
plan.rs comment and the stale --memory-budget-mb / --max-read-len help text.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Test the exit-4 realized-breach path (Item E)

**Files:** `src/main.rs` (force-RSS seam at the realized-peak capture), `tests/plan_enforce.rs`.

- [ ] **Step 1: Add the post-run-only force seam** (`src/main.rs`, the realized-peak capture ~1352)

Replace `let peak_rss = peak_rss_bytes();` with:

```rust
    // Realized peak (monotonic high-water mark) captured after the calling pass.
    // Test-only seam: ROSALIND_FORCE_PEAK_RSS_BYTES overrides ONLY the post-run
    // realized peak (never the pre-run baseline at the --enforce gate), so the
    // exit-4 breach branch can be exercised deterministically without allocating.
    let peak_rss = std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or_else(peak_rss_bytes);
```

- [ ] **Step 2: Write the exit-4 test** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn enforce_breach_exits_4_after_writing_output_and_receipt() {
    // Budget 4096 MiB passes the pre-run exit-3 gate (tiny fixture), but a forced
    // realized peak of 8 GiB trips the post-run breach -> exit 4, with the VCF +
    // receipt still written (the documented "output + receipt written" semantics).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .env("ROSALIND_FORCE_PEAK_RSS_BYTES", "8589934592") // 8 GiB
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(4), "expected breach exit 4: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("VIOLATED"), "missing VIOLATED line: {stderr}");
    // Output + receipt were still written before the breach exit.
    assert!(vcf.exists(), "VCF must be written before exit 4");
    let json = std::fs::read_to_string(&manifest).expect("receipt written");
    assert!(
        json.contains("\"contract_verdict\":\"over\""),
        "verdict should be over: {json}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 3: Run + commit**

Run: `cargo test --test plan_enforce` (incl. the new breach test) ; then `cargo test` full ;
`cargo build --release 2>&1 | grep -c warning: || true` ; `cargo fmt --all -- --check`
Expected: all green, 0 warnings.

```bash
git add -A
git commit -m "test(enforce): cover the exit-4 realized-breach backstop

Adds a post-run-only ROSALIND_FORCE_PEAK_RSS_BYTES seam to deterministically
force predicted-fits-yet-realized-overruns, and asserts exit 4 + VIOLATED +
VCF/receipt still written. The contract's most safety-critical branch was
previously untested.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: CI memory-envelope gate on the bounded `--index` path (Item F)

**Files:** `.github/workflows/ci.yml`.

- [ ] **Step 1: Add a contract-gate step to the `cli-e2e` job**

After the existing "Call variants" step (which uses the legacy `--reference` path), append steps that
exercise the bounded `--index` contract on the toy data. (The toy `alignments.bam` is produced by the
existing "Align to BAM" step; build an index from `reference.fa`, sort the BAM, then enforce.)

```yaml
      - name: Build index for the contract path
        run: |
          cargo run --release -- index \
            --reference examples/data/illumina_toy/reference.fa \
            --output examples/data/illumina_toy/reference.idx
      - name: Sort the BAM
        run: |
          cargo run --release -- sort \
            --input examples/data/illumina_toy/alignments.bam \
            --output examples/data/illumina_toy/sorted.bam
      - name: Contract gate — fits, breaches, refuses
        run: |
          set -e
          # (1) Generous budget fits -> exit 0 + "contract: OK".
          cargo run --release -- variants --index examples/data/illumina_toy/reference.idx \
            --alignments examples/data/illumina_toy/sorted.bam \
            --memory-budget-mb 4096 --enforce \
            -o examples/data/illumina_toy/contract.vcf 2> ok.log
          grep -q "contract: OK" ok.log
          # (2) verify the receipt without re-running -> exit 0.
          cargo run --release -- verify \
            --manifest examples/data/illumina_toy/contract.vcf.manifest.json | grep -q "verify: OK"
          # (3) 1 MiB budget refuses up front -> exit 3, no VCF.
          rm -f examples/data/illumina_toy/refused.vcf
          set +e
          cargo run --release -- variants --index examples/data/illumina_toy/reference.idx \
            --alignments examples/data/illumina_toy/sorted.bam \
            --memory-budget-mb 1 --enforce \
            -o examples/data/illumina_toy/refused.vcf 2> refuse.log
          code=$?
          set -e
          test "$code" -eq 3
          grep -q "REFUSE" refuse.log
          test ! -f examples/data/illumina_toy/refused.vcf
```

- [ ] **Step 2: Validate the YAML locally (best-effort) + commit**

Run (best-effort lint; skip if `yamllint` absent):
`yamllint .github/workflows/ci.yml || echo "yamllint not installed — skipping"`
Manually re-read the indentation against the existing steps in the file.

```bash
git add .github/workflows/ci.yml
git commit -m "ci: exercise the bounded variants --index contract (fits / refuse / verify)

Adds a contract gate to cli-e2e: build index -> sort -> variants --index
--enforce at a generous budget (exit 0 + contract: OK), verify the receipt,
and a 1 MiB budget (exit 3 refuse, no VCF). Makes 'the bounded contract is
exercised in CI' literally true, via deterministic exit codes (no RSS-noise
flakiness; exit-4 is covered deterministically by the force-seam unit test).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Finish — verify, push, PR (no merge)

**Files:** none (release hygiene).

- [ ] **Step 1: Full green gate**

Run: `cargo fmt --all -- --check` ; `cargo test` ;
`cargo build 2>&1 | grep -c warning: || true` ; `cargo build --release 2>&1 | grep -c warning: || true`
Expected: fmt clean, full suite green, 0 warnings (debug + release).

- [ ] **Step 2: Skim the diff for honesty + scope**

Run: `git log --oneline main..HEAD` and `git diff --stat main..HEAD`. Confirm only the 6 keystone
items changed; no cross-cutting (sort.rs, bam.rs @SQ, rss.rs unit fix) leaked in.

- [ ] **Step 3: Push the branch + open the PR (DO NOT MERGE)**

```bash
git push -u origin rosalind/trust-correctness-keystone
gh pr create --title "Trust + correctness hardening keystone" --body "$(cat <<'EOF'
## Summary
The 6 call/contract-path items from the 2026-06-02 reflection audit:
- **Unbiased depth cap** — replace leftmost-arrival truncation with a fixed-seed min-hash bounded reservoir (no more silent het drops at deep sites); deterministic + bounded at `max_depth`.
- **Canonical obs order** — `build_column` emits obs sorted by `(allele, base_qual, mapq, reverse)` so the f64 QUAL/PL sum is order-independent (byte-identical VCF).
- **Surface skip counts** — thread `SkipCounts` out to a stderr summary + `over_max_depth`/`reads_skipped_total` in the receipt.
- **`--max-read-len` enforcement** — under `--enforce`, an over-long read aborts loudly instead of voiding the predicted envelope; fixes stale help/comment text.
- **Exit-4 breach test** — a post-run-only `ROSALIND_FORCE_PEAK_RSS_BYTES` seam makes the "never a silent overrun" backstop testable; asserts exit 4 + `VIOLATED` + output/receipt written.
- **CI memory gate** — `cli-e2e` now exercises the bounded `variants --index` contract (fits / refuse / verify).

Deferred to a fast-follow (out of scope here): sort tie-break, `@SQ` length check, BSD `ru_maxrss` unit fix + cgroup-awareness, front-door honesty pass, merging PR #23.

## Test plan
- [ ] `cargo test` green (incl. new unbiased-cap, obs-order, skip-count, max-read-len, and exit-4 tests)
- [ ] `cargo fmt --all -- --check` clean; 0 warnings (debug + release)
- [ ] CI contract gate passes (fits exit 0 / refuse exit 3 / verify OK)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 4: Report the PR URL to the user; await merge authorization.** Do NOT merge to `main`.
