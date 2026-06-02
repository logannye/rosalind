# Act 1 — Contract Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Execute INLINE (not via subagents) — this repo is a single shared tree; subagents running `cargo`/git concurrently corrupt it ([[feedback_subagent_shared_tree_no_cargo_fix]]).

**Goal:** Make the memory contract *true and trusted on real genomes* — close the gaps where the contract holds on the synthetic/single-contig demo path but breaks on a real multi-contig human reference, and where the advertised first-60-seconds isn't actually verifiable.

**Architecture:** Three focused PRs in sequence, each independently shippable. This plan details **PR-A (soundness)** fully; **PR-B (real-genome correctness)** and **PR-C (trust on-ramp + hygiene)** are scoped at the end and will each get their own detailed plan when reached.

**Tech Stack:** Rust 2021 (MSRV 1.72), rust-htslib (BAM/VCF), BLAKE3 receipts, clap CLI, subprocess integration tests.

**Source of truth:** the 60-agent reflection audit (`/private/tmp/.../tasks/audit.json`), findings triaged in [[project_rosalind_strategy]]. Every fix below traces to an adversarially-verified gap.

---

## The Act-1 sequence (header — full detail for PR-A only)

- **PR-A — Soundness: make the predicted peak a true upper bound.** (§A) The keystone. The `Arc::from(decoded)` reallocation in `whole_genome.rs:69` makes two ASCII-decoded copies of the largest contig briefly co-resident, but the estimator (`plan.rs:33`) and accountant (`engine.rs:195`) both count the reference once — so a `[FITS]` plan can breach at runtime, caught only post-hoc by exit 4. And **no test compares predicted *peak RSS* to realized *peak RSS***. This PR records the prediction in the receipt, adds the missing real-peak soundness test, kills the transient, and gives the estimator an honest RSS-level margin.
- **PR-B — Real-genome correctness.** (§B) Multi-contig `eval` (load all FASTA records + per-variant reference lookup), contig-naming guard (error when BAM `@SQ` ∩ index names is empty; warn when 0 reads contributed), IUPAC ambiguity codes → `N` at index build. (Somatic streaming-co-walk + contract surface is a larger separate item — note, don't bundle.)
- **PR-C — Trust on-ramp + hygiene.** (§C/§D) Fix the README/CHANGELOG Action snippet (`logannye/rosalind-budget@v1` → `logannye/rosalind@v0.1.0`); make `install.sh` actually verify the `.sha256`; fix the CI `cli-e2e` checksum tautology (`--length 4000`); add `verify` manifest self-hash + cross-consistency; address the `IndexFreeIterator` public-API `unimplemented!()` panic; track `Cargo.lock`.

---

## PR-A: Soundness — predicted peak is a true upper bound

**Files:**
- Modify: `src/main.rs` (record prediction in receipt; ~`run_variants_index` near line 1634/1722, and the features path ~line 1373)
- Modify: `src/call/plan.rs` (RSS-level margin in `predicted_peak_rss_bytes` + `render_variants_plan`)
- Modify: `src/core/budget.rs` (new `PILEUP_IO_RSS_OVERHEAD` constant + re-export)
- Modify: `src/genomics/index/view.rs` (new `decode_window_arc`)
- Modify: `src/call/whole_genome.rs` + `src/call/features.rs` (use `decode_window_arc`; kill the transient)
- Test: `tests/plan_enforce.rs` (new real-peak soundness test), `src/genomics/index/view.rs` (unit test), `src/call/plan.rs` (update value test)

### Task 1: Record the predicted peak in the receipt

The receipt records `peak_rss_bytes` (realized) but not the prediction. Recording it (a) enables the soundness assertion to read both numbers from one manifest, and (b) is exactly what the Act-2 `plan --fleet` scheduler needs. Compute it unconditionally (not only under `--enforce`).

- [ ] **Step 1: Locate the realized-peak recording in `run_variants_index`.** In `src/main.rs`, the post-run block records `peak_rss_bytes` (~line 1781) and `max_working_set_bytes` (~line 1783). The `--enforce` gate (~line 1634) already computes `let baseline = peak_rss_bytes();` then `let predicted = predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline);`. Currently that lives *inside* the `if enforce {…}` block.

- [ ] **Step 2: Lift baseline+predicted out of the enforce block so it's always computed.** Just before calling begins (before the `if enforce` gate), add:
  ```rust
  // Measured baseline (binary + libs + index/BAM open) → predicted peak RSS.
  // Recorded in the receipt unconditionally: it is the contract's up-front
  // claim, and the post-run check + `verify` assert the realized peak honors it.
  let predicted_peak = rosalind::call::plan::predicted_peak_rss_bytes(
      largest, max_depth, max_read_len, peak_rss_bytes(),
  );
  ```
  Have the existing `--enforce` gate reuse `predicted_peak` instead of recomputing (keep its REFUSE message). Confirm `largest`, `max_depth`, `max_read_len` are in scope at that point (the enforce gate uses them, so they are).

- [ ] **Step 3: Record it in the receipt params** alongside `peak_rss_bytes`:
  ```rust
  params.insert("predicted_peak_rss_bytes".to_string(), predicted_peak.to_string());
  ```

- [ ] **Step 4: Mirror the same recording in the features path** (`run_features`, ~line 1373) so `features --index` receipts also carry the prediction (it shares the engine and contract).

- [ ] **Step 5: Build.** Run: `cargo build` — Expected: compiles, 0 warnings.

### Task 2: The failing real-peak soundness test (TDD red)

- [ ] **Step 1: Add the test to `tests/plan_enforce.rs`.** It builds an index with ONE multi-MiB contig and SHALLOW coverage (so the reference decode, not the active set, is the high-water — this is exactly where the transient dominates), runs the real binary, and asserts the recorded prediction bounds the recorded realized peak:
  ```rust
  #[test]
  fn predicted_peak_rss_upper_bounds_realized_peak() {
      // A multi-MiB contig with SHALLOW coverage: the reference decode is the
      // RSS high-water, so the Arc::from transient (two copies of the contig)
      // shows up in realized peak_rss but NOT in a reference-counted-once model.
      // This is the contract's core inequality (predicted peak >= realized peak),
      // which no existing test exercises (the others compare working-set vs
      // working-set, both modeling the reference exactly once).
      let dir = unique_dir("peak-soundness");
      let idx = dir.join("ref.idx");

      // ~8 MiB single contig of 'A'.
      let bases = vec![b'A'; 8 * 1024 * 1024];
      let index = rosalind::genomics::GenomeIndex::from_named_sequences(&[(
          "chr1".to_string(),
          bases,
      )])
      .unwrap();
      rosalind::genomics::IndexWriter::create(&idx)
          .unwrap()
          .write_genome_index(&index)
          .unwrap();

      // A handful of short reads → shallow coverage; reference dominates RSS.
      let bam = dir.join("reads.bam");
      write_sorted_bam(&bam, &idx, &shallow_reads()); // helper below

      let vcf = dir.join("calls.vcf");
      let manifest = dir.join("calls.vcf.manifest.json");
      let out = Command::new(bin())
          .args(["variants", "--index"])
          .arg(&idx)
          .arg("--alignments")
          .arg(&bam)
          .args(["--max-depth", "1000", "--max-read-len", "250", "-o"])
          .arg(&vcf)
          .output()
          .unwrap();
      assert!(out.status.success(), "run failed: {out:?}");

      let text = std::fs::read_to_string(&manifest).unwrap();
      let m = rosalind::provenance::RunManifest::from_canonical_json(&text).unwrap();
      let predicted: u64 = m.params.get("predicted_peak_rss_bytes").unwrap().parse().unwrap();
      let realized: u64 = m.params.get("peak_rss_bytes").unwrap().parse().unwrap();
      assert!(
          predicted >= realized,
          "predicted peak RSS {predicted} must be >= realized peak RSS {realized} \
           (gap = the unmodeled reference-decode transient)"
      );
      std::fs::remove_dir_all(&dir).ok();
  }
  ```
  Reuse the fixture helpers already in the file (`unique_dir`, `bin`, the sorted-BAM builder used by `build_sorted_bam_fixture`). If a small standalone `write_sorted_bam`/`shallow_reads` helper is needed, factor it from the existing fixture builder — do NOT duplicate htslib glue.

- [ ] **Step 2: Run it — expect RED.** Run: `cargo test --test plan_enforce predicted_peak_rss_upper_bounds_realized_peak -- --nocapture`. Expected: **FAIL** — realized peak ≈ baseline + ~16 MiB (8 MiB Arc + 8 MiB transient Vec), predicted ≈ baseline + ~8 MiB. This pins the soundness gap. (If it unexpectedly passes, the contig is too small for the transient to clear measurement noise — raise to 16–32 MiB.)

### Task 3: Kill the reference-decode transient (green, part 1)

- [ ] **Step 1: Add `decode_window_arc` to `ReferenceView`** in `src/genomics/index/view.rs`, building the `Arc<[u8]>` directly from a `TrustedLen` range iterator so there is no full intermediate ASCII `Vec` (no two-copies-resident spike):
  ```rust
  /// Decode `[start, end.min(len))` directly into an `Arc<[u8]>` without
  /// materializing a separate owned `Vec` first — the collect writes into the
  /// Arc's allocation in place (the range iterator is `TrustedLen`), so peak
  /// resident reference memory is ONE copy, not the transient two that
  /// `Arc::from(vec)` (which reallocates) would hold. Same bytes as
  /// `decode_window`; bounded by the window size.
  pub fn decode_window_arc(&self, start: usize, end: usize) -> std::sync::Arc<[u8]> {
      let end = end.min(self.len);
      (start..end).map(|i| self.base_at(i)).collect()
  }
  ```

- [ ] **Step 2: Unit test byte-equivalence with `decode_window`** (in `view.rs` tests): decode the same window both ways on a small fixture, assert the `Arc` slice equals the `Vec`. This guarantees zero behavior change in the reference bytes the caller sees.

- [ ] **Step 3: Use it in `call_germline_whole_genome`** (`src/call/whole_genome.rs:67-69`). Replace:
  ```rust
  let mut decoded = Vec::new();
  ref_view.decode_window(start, end, &mut decoded);
  let reference: Arc<[u8]> = Arc::from(decoded);
  ```
  with:
  ```rust
  let reference: Arc<[u8]> = ref_view.decode_window_arc(start, end);
  ```
  Update the comment block (lines 62-64) to state the transient is gone.

- [ ] **Step 4: Use it in the features whole-genome driver** too — grep `decode_window` + `Arc::from` in `src/call/features.rs` and apply the same replacement (the features `--enforce` path inherits the identical contract).

- [ ] **Step 5: Run the existing equivalence + bounded tests.** Run: `cargo test --test '*' whole_genome && cargo test -p rosalind --lib call::features`. Expected: PASS (byte-identical calls; `whole_genome_equals_per_contig_calls` still green).

### Task 4: Honest RSS-level margin (green, part 2)

The working-set estimate stays the pure pileup working set (comparable to the engine's accountant). The *RSS* prediction adds a named margin for what realized peak RSS includes but the working set omits: the BGZF decompression input buffer, the VCF `BufWriter`, and allocator slack.

- [ ] **Step 1: Add the constant to `src/core/budget.rs`** next to the other `PILEUP_*` constants, and re-export it from the crate's core prelude the way the others are:
  ```rust
  /// Peak-RSS overhead the streaming *working set* does not model: BGZF input
  /// buffer + VCF BufWriter + allocator slack between RSS and live bytes. Added
  /// to the predicted peak RSS (not the working-set estimate) so the contract's
  /// up-front claim is a conservative upper bound, not a steady-state count.
  pub const PILEUP_IO_RSS_OVERHEAD: u64 = 16 * 1024 * 1024; // 16 MiB; calibrated by the soundness test
  ```
  (Start at 16 MiB; Task 4 Step 3 calibrates it against the measured gap.)

- [ ] **Step 2: Add it in `predicted_peak_rss_bytes` and `render_variants_plan`** (`src/call/plan.rs`) — at the RSS level only, NOT in `estimate_variants_working_set`:
  ```rust
  pub fn predicted_peak_rss_bytes(/* …unchanged args… */) -> u64 {
      baseline_rss_bytes
          .saturating_add(estimate_variants_working_set(largest_contig_len, max_depth, max_read_len).bytes)
          .saturating_add(PILEUP_IO_RSS_OVERHEAD)
  }
  ```
  Add a matching `io overhead (RSS): N MiB` line to the `render_variants_plan` breakdown so `rosalind plan` stays transparent about the margin.

- [ ] **Step 3: Update the value-level unit test** `predicted_peak_is_baseline_plus_working_set` in `plan.rs` to `assert_eq!(predicted, 50_000_000 + ws + PILEUP_IO_RSS_OVERHEAD)`. The working-set exact-value test (`estimate_grows_…`) is unchanged (the working-set estimate did not move). Leave `estimator_upper_bounds_the_realized_working_set` unchanged (still holds).

- [ ] **Step 4: Run the soundness test — expect GREEN.** Run: `cargo test --test plan_enforce predicted_peak_rss_upper_bounds_realized_peak -- --nocapture`. Expected: PASS, with the transient eliminated (Task 3) and the margin (Task 4) covering io/allocator slack. If it still fails, read the printed predicted/realized gap and raise `PILEUP_IO_RSS_OVERHEAD` to cover it with headroom (this is the calibration step — the test is the spec).

### Task 5: Full verification + commit

- [ ] **Step 1: Whole suite.** Run: `cargo test` — Expected: all sections green (the new test + all existing, including `frontdoor_demo`, `germline_accuracy`, golden VCF).
- [ ] **Step 2: Lints + format.** Run: `cargo fmt --all && cargo clippy --all-targets -- -D warnings && cargo build --release`. Expected: clean, 0 warnings debug+release.
- [ ] **Step 3: Re-run the flagship contract path locally** on the toy fixture to confirm `plan`/`--enforce`/`verify` still behave and the receipt now carries `predicted_peak_rss_bytes`. Eyeball that predicted ≥ realized in the receipt.
- [ ] **Step 4: Commit** (use `git commit -F -` with a quoted heredoc — backticks in `-m` trigger shell substitution, [[project_rosalind_strategy]] lesson):
  ```
  fix(contract): predicted peak RSS is a true upper bound on real genomes

  - kill the Arc::from(Vec) reference-decode transient (decode_window_arc)
  - record predicted_peak_rss_bytes in the receipt (unconditional)
  - honest PILEUP_IO_RSS_OVERHEAD margin at the RSS level
  - the missing real peak-RSS soundness test (predicted >= realized)
  ```

---

## Self-review (run before opening the PR)

1. **Spec coverage:** §A's two halves — eliminate the transient AND make the estimator honest — both land, and the real-peak inequality is now a test (the audit's named gap "the contract's core soundness claim is NOT tested" is closed).
2. **Soundness direction:** the margin only ever makes `predicted` larger → the contract can over-refuse but never under-predict. Correct safe direction.
3. **No behavior change in calls:** `decode_window_arc` is byte-equivalent to `decode_window` (unit-tested); `whole_genome_equals_per_contig_calls` and the golden VCF pin output bytes.
4. **Type consistency:** `decode_window_arc` returns `Arc<[u8]>` — the exact type `call_germline_region_streaming` already takes; no signature ripple.

## Follow-on (PR-B, PR-C) — to be detailed when reached

- **PR-B (correctness):** `read_fasta` → load ALL records into a `name→seq` map for `eval`; `compare_callsets` looks up the per-variant reference by `chrom` (error on absent contig, naming it); `StreamingBamSource::new` errors when no `@SQ` name resolves into the index and `run_variants_index` warns loudly when `reads_used == 0`; `sanitize_reference` maps IUPAC codes to `N` (report contig:coord on the strict path).
- **PR-C (trust/hygiene):** README/CHANGELOG `rosalind-budget@v1` → `logannye/rosalind@v0.1.0` (+ a CI lint that the snippet ref resolves); `install.sh` fetches and checks the `.sha256`; CI `cli-e2e` passes `--length 4000` and stops overwriting the committed fixture; `verify` adds a manifest self-hash + `max_working_set_bytes <= peak_rss_bytes` / verdict-vs-budget consistency checks; doc-or-remove the `IndexFreeIterator::next_item()` `unimplemented!()` from the public API; track `Cargo.lock`.
