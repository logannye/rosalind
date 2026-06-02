# Trust + Correctness Hardening — Keystone (design)

**Status:** DESIGN SPEC — 2026-06-02. Approved decisions captured below; ready for the
implementation plan after user review. **Branch:** `rosalind/trust-correctness-keystone` (off `main`
`b6fdfda`). **Source:** the 2026-06-02 reflection audit (16-agent workflow) — see the conversation
record; the audit's `strongest_gap` (depth-cap silent variant drops) and `strongest_next_step`
(exit-4 + CI gate untested) are the two anchors of this increment.

## 1. Goal

Make the shipped memory contract's **correctness** and **trust** guarantees actually hold and be
tested, before we deepen or grow the system. Concretely: stop the variant caller from silently
dropping real variants, and regression-guard the two claims most central to the pitch
("honor-or-breach / never a silent overrun" and "CI-enforced bound"), which the audit found are the
two **least** covered.

This is the **keystone** increment. It is scoped to the call/contract path so it is one coherent,
reviewable PR. Cross-cutting hardening (sort tie-break, `@SQ` length check, BSD `ru_maxrss` unit
fix + cgroup-awareness, manifest provenance, front-door honesty pass, merging PR #23) is sequenced
to a **fast-follow** increment (§7) and is explicitly out of scope here.

## 2. Approved decisions (locked)

1. **Depth-cap algorithm:** unbiased **hash-priority bounded reservoir** (min-hash with eviction).
2. **Default `--max-depth`:** keep **1000** (now unbiased), plus a stderr warning + a receipt field
   when the cap engages.
3. **First-PR scope:** the **6 keystone items** (§4). Guards + cross-cutting → fast-follow.

## 3. Background: the two verified problems

- **Silent variant drops (correctness blocker).** `PileupEngine::advance_to`
  (`src/pileup/engine.rs:249-254`) refuses *new arrivals* once `self.active.len() >= max`, retaining
  reads that started upstream. A read's start position correlates with which allele it carries near
  that start, so the cap is **allele-biased**: an agent reproduced a true 30:30 het SNV becoming
  **zero variant calls** at `--max-depth 30` (the alt reads, which start at the variant, are all
  dropped → the site reads as hom-ref → abstain). The default cap is 1000, so this bites at any site
  deeper than 1000× (panels/amplicon/deep-tumor). The `over_max_depth` skip count is computed
  (`engine.rs:251`) but **never surfaced** on the `--index` path or in the receipt.
- **f64 order-sensitivity (propagation).** `site_likelihoods` (`src/call/germline.rs:48-58`)
  accumulates `log_l[0/1/2] += …` by summing over `column.obs` **in active-set order**. f64 addition
  is non-associative, so any fix that reorders the active set would change QUAL/PL bytes. `build_column`
  (`engine.rs:271-292`) emits `obs` in `self.active` Vec order. **The cap fix and an obs-order
  canonicalization must land together** to preserve byte-identical output.

## 4. The six keystone items

### Item A — unbiased hash-priority depth-cap reservoir

**Files:** `src/pileup/engine.rs` (modify).

Replace leftmost-arrival truncation with a **min-hash bounded reservoir**:

- Add a fixed-seed content hash to each active read. `priority(read) = fnv1a64(pos.0, end, flags.0,
  cigar bytes, &seq, &qual)`. **Fixed seed** (a hand-rolled FNV-1a, *not* `std::collections`'
  `RandomState`/`DefaultHasher`, whose per-process seed would break determinism). The hash is
  independent of arrival order and uncorrelated with the allele at any single column — this is what
  removes the bias.
- `ActiveRead` gains a `priority: u64` field (computed once in `ingest`).
- In `advance_to`, when a read passes filters and reaches the cursor:
  - if `self.active.len() < max` → ingest (admit);
  - else find the resident with the **maximum** priority; if `new.priority < worst.priority` →
    evict that resident (`self.skips.over_max_depth += 1`) and ingest the newcomer; else refuse the
    newcomer (`self.skips.over_max_depth += 1`).
  - `over_max_depth` counts **every read removed by the cap** (refused arrival OR evicted resident);
    documented as such.
- Tie semantics: equal priorities ⇒ `new.priority < worst.priority` is false ⇒ no eviction ⇒
  first-`max` retained. So genuinely identical reads (same content ⇒ same hash) keep the old
  first-N behavior — the existing `max_depth_caps_active_set_deterministically` test still passes.
- Finding the max-priority resident: a linear scan over `self.active` (≤ `max` entries) is simple and
  deterministic; acceptable for the first correct version (optimization to a heap is a later perf
  fast-follow, not needed for correctness).

**Determinism argument:** given a canonical (coordinate-sorted) input stream, "keep the `max`
smallest-priority reads among those covering the cursor that we still hold" is order-independent
(verified by case analysis in the design discussion: arrival orders R1,R2,R3 / R3,R2,R1 / R2,R1,R3
all converge to the same kept set). Memory stays ≤ `max` active reads (the bounded-memory invariant
is preserved exactly).

**Tests (engine.rs):**
- NEW `unbiased_cap_keeps_both_alleles_at_a_deep_het`: reference `A…`; `max_depth = 30`; 30 ref reads
  whose spans start upstream of the variant column and 30 alt reads starting at the variant column
  (so the biased cap would drop all alt). Assert the variant column observes **both** alleles
  (alt count > 0) — i.e. the reproduced bug is fixed.
- NEW `cap_is_order_independent_with_distinct_reads`: distinct-content reads in two input orders
  produce identical kept allele counts and identical `over_max_depth`.
- KEEP `max_depth_caps_active_set_deterministically` (identical reads ⇒ ties ⇒ unchanged) and
  `working_set_is_bounded_by_coverage_not_input_size` green.

### Item B — obs-order canonicalization (fixes the f64 QUAL/PL order-sensitivity)

**Files:** `src/pileup/engine.rs` (`build_column`).

After collecting `obs` for a column, **sort** it by a total key before returning:
`(allele, base_qual, mapq, reverse as u8)`. Two observations equal on that key contribute identical
terms to the log-likelihood sum, so the f64 accumulation in `germline.rs` becomes order-independent
and byte-identical regardless of active-set order. Cost: sort of ≤ `max_depth` items per column —
negligible.

**Tests:**
- NEW (germline.rs): `qual_pl_independent_of_obs_order` — same multiset of observations fed in two
  orders yields byte-identical `qual`/`pl`. (Construct two `PileupColumn`s with shuffled `obs` and
  assert equal `GermlineCall`.)
- The existing determinism test (`tests/determinism.rs`) and golden VCF (`tests/golden_vcf.rs`) must
  stay green (the canonical order must not change the *current* golden output — verify; if the golden
  ordering differs, regenerate the snapshot under `ROSALIND_UPDATE_SNAPSHOTS` and note it in the PR).

### Item C — surface SkipCounts on the `--index` path + in the receipt

**Files:** `src/call/pipeline.rs`, `src/call/whole_genome.rs`, `src/main.rs` (modify).

The engine already tracks `SkipCounts` and exposes `skip_counts()`. Plumb it out:

- `call_germline_region_streaming` returns `Result<(WorkingSet, SkipCounts), CoreError>` (read
  `engine.skip_counts()` after the iteration loop). Update `call_germline_region_tracked` and
  `call_germline_region` accordingly (they discard the new field).
- `call_germline_whole_genome` returns `Result<(WorkingSet, SkipCounts), CoreError>`, **summing**
  `SkipCounts` across contigs (add a `SkipCounts` accumulator; needs field-wise add — add an
  `impl std::ops::Add` or a `merge`/`accumulate` method on `SkipCounts`).
- `run_variants_index` (`src/main.rs:1294-1350`): capture the returned `SkipCounts`; after the memory
  line, print a stderr summary **only when any cap drop occurred**, e.g.
  `pileup: dropped <over_max_depth> reads at --max-depth <d> (deep-site downsampling)`; print the
  other skip reasons (low_mapq, secondary, …) if non-zero too (one concise line).
- Manifest (`src/main.rs:1374-1428`): add params `over_max_depth` and `reads_skipped_total`
  (`SkipCounts::total()`), so the receipt is honest about whether the cap was load-bearing.

**Tests:**
- Extend the whole-genome equivalence test (`whole_genome.rs`) to assert the summed `SkipCounts`
  equals the per-contig sum.
- `tests/plan_enforce.rs` (or `tests/variants_index.rs`): a `--max-depth`-engaging run records
  `over_max_depth > 0` in the manifest; a non-engaging run records `0`.

### Item D — `--max-read-len` ingest enforcement (close the live silent-OOM path under `--enforce`)

**Files:** `src/pileup/engine.rs` (`PileupParams`, `advance_to`), `src/main.rs`, `src/call/plan.rs`
(comment), `src/core` (error variant if needed).

`--max-read-len` is currently prediction-only (`plan.rs`) and never checked at ingest, so a BAM with
reads longer than the declared cap breaks the predicted envelope by ~40–200× in the under-estimation
direction (the live "plan says FITS → OOM" path).

- `PileupParams` gains `max_read_len: Option<u32>` (default `None` = no check). It is set to
  `Some(max_read_len)` **only under `--enforce`** (so non-enforced runs stay backward-compatible and
  accept long reads, just unbudgeted).
- In `advance_to`, before ingesting a read that reaches the cursor, if `Some(m) = max_read_len` and
  the read's `seq.len() as u32 > m`, return `Err(CoreError::…)` with a clear message:
  `read length <N> exceeds declared --max-read-len <M>; raise --max-read-len or drop --enforce`.
  Under `--enforce` this aborts the run loudly (propagates to a non-zero exit) rather than silently
  blowing the envelope. (A reference variant of `CoreError` already exists; reuse the closest generic
  one or add `ReadExceedsDeclaredLength`. Decide in the plan after reading `src/core/error.rs`.)
- Honesty fixes (same item): correct the contradictory `plan.rs:19-20` comment (it claims
  `max_read_len` is "enforced at runtime"); fix the stale `main.rs:96-97` `--memory-budget-mb` help
  ("does not enforce (enforcement is a later phase)" — `--enforce` now exists); update the
  `--max-read-len` flag help (`main.rs:104-106`) to state it is enforced at ingest under `--enforce`.

**Tests:**
- engine.rs: a read longer than `max_read_len` (when `Some`) yields `Err`; when `None`, long reads
  pile up normally (the existing `long_read_piles_up_every_matched_base` stays green).
- `tests/plan_enforce.rs`: `variants --index --enforce` on an alignment set containing an over-long
  read exits non-zero with the clear message.

### Item E — test the exit-4 realized-breach path

**Files:** `src/util/rss.rs` (test seam) or `src/main.rs` (capture site), `tests/plan_enforce.rs`.

The exit-4 branch (`src/main.rs:1449`, the literal "never a silent overrun" backstop) has zero tests.
To force *predicted-fits-yet-realized-overruns* deterministically (without allocating gigabytes):

- Add a **post-run-only** test seam: at the realized-peak capture (`main.rs:1352`), read
  `ROSALIND_FORCE_PEAK_RSS_BYTES`; if set and parseable, use it as the realized peak instead of
  `peak_rss_bytes()`. **It must NOT affect the pre-run baseline** (`main.rs:1267`), so the up-front
  predicted peak still uses the real (small) baseline and the run passes the exit-3 gate. Document the
  env var as a test-only seam in code.
- Test (`tests/plan_enforce.rs`): run `variants --index --enforce --memory-budget-mb B -o <vcf>` on
  the bundled fixture with `B` chosen above the real predicted peak (so exit-3 passes) and
  `ROSALIND_FORCE_PEAK_RSS_BYTES` set above `B` (so the realized check fails). Assert: **exit 4**,
  stderr contains `VIOLATED`, and the **VCF + manifest were still written** (the documented
  "output + receipt written" semantics). Add a companion assertion that the manifest's
  `contract_verdict` is `over`.

### Item F — CI memory-envelope gate on the bounded `--index` path

**Files:** `.github/workflows/ci.yml` (modify).

CI has no RSS gate today and its e2e job exercises only the legacy `--reference` path. Add a job (or
extend `cli-e2e`) that exercises the **bounded `--index` contract** on a fixture, asserting the
contract fires — using deterministic, non-flaky signals (not a tight RSS threshold on a noisy
runner):

- index → sort → `variants --index --enforce --memory-budget-mb <fits> -o out.vcf` → assert exit 0
  and stderr `contract: OK — realized peak … within`.
- `variants --index --enforce --memory-budget-mb 1` → assert **exit 3** (`REFUSE`, pre-run,
  deterministic) and that **no VCF was written**.
- `verify` the receipt → assert exit 0.
- Optionally assert the receipt's `max_working_set_bytes` (deterministic) is below a generous bound.

This makes "the bounded contract is exercised in CI" literally true via the `--index` path and the
exit-3 refuse, without depending on run-to-run RSS variance. (The exit-4 branch is covered
deterministically by Item E's force-seam test, which runs under `cargo test` in CI.)

## 5. Cross-cutting requirements

- **Determinism is the non-negotiable invariant.** After Items A+B, the `variants --index` VCF must be
  byte-identical regardless of input read order and of the (content-based) cap selection. The golden
  VCF + determinism tests are the guard.
- **No new warnings** (debug + release); `cargo fmt --check` clean; full `cargo test` green.
- **Backward compatibility:** non-`--enforce` runs behave as before except for the (now unbiased) cap
  and the canonical obs order. The cap default stays 1000.
- **Honesty:** every place a number or guarantee is surfaced (stderr, receipt, help text) must match
  the code. The Item-D help-text fixes are part of this.

## 6. Risks

- **Golden-snapshot churn (low):** obs-order canonicalization could change the byte order of the
  current golden VCF if the existing arrival order happened to differ from the canonical key order.
  Mitigation: verify; regenerate the snapshot deliberately under `ROSALIND_UPDATE_SNAPSHOTS` if needed
  and call it out in the PR (the *values* must be unchanged — only obs accumulation order, which must
  not change emitted QUAL/PL beyond the documented canonicalization).
- **Eviction accounting (low):** counting both refused arrivals and evicted residents as
  `over_max_depth` slightly overcounts "reads not present in any column," but it honestly reports
  "reads the cap removed." Documented; the point is to surface that the cap engaged, not exact
  bookkeeping.
- **`--max-read-len` abort leaves a partial VCF (low/acceptable):** under `--enforce`, an over-long
  read aborts mid-run after some rows are written. This is the *safe* direction (loud failure, not a
  silent OOM); documented. A clean pre-scan would defeat streaming, so we accept the loud mid-run
  abort.

## 7. Out of scope (fast-follow increment 2)

Deferred, with the audit severities noted — to be specced separately after this lands:
- **sort k-way-merge tie-break** (finding #5, high) — `src/genomics/sort.rs`: add `source_idx`
  (global arrival index) as the final tie-break in `HeapItem::cmp` so sort output is budget-invariant.
  (Note: Item B already defuses the *VCF* propagation of this; the remaining issue is sorted-BAM
  byte-identity across `--memory-mb`.)
- **`@SQ` contig-length cross-check** (finding #7, high) — `src/io/bam.rs`: reject a BAM whose header
  contig lengths disagree with the index.
- **BSD `ru_maxrss` unit fix + cgroup-awareness + manifest os/arch provenance** (findings #6, medium)
  — `src/util/rss.rs`, `src/provenance/mod.rs`.
- **Front-door honesty pass** (low) — kill the `O(√t)` overclaim in `Cargo.toml` description, the CLI
  `--help` about-string, and `scale_test_results.txt`; fix `CONTRACT.md`'s "byte-identical manifest".
- **Merge PR #23** (the 41 B/base `BuildMemoryModel`) so `plan --reference` stops under-predicting
  build RAM ~3.4× — a merge decision reserved for the user.

## 8. Self-review

- **Spec coverage:** all 6 approved items have a file list, a design, and tests. The two anchors
  (depth-cap blocker; exit-4 + CI) are Items A/C and E/F. ✓
- **Type consistency:** `SkipCounts` return-type change is threaded through all three call-path
  functions (`_streaming`, `_tracked`, `_region`, `_whole_genome`) and the `main.rs` call sites;
  `PileupParams` gains one optional field used in `advance_to`. ✓
- **No placeholders:** the one deliberately-deferred detail is the exact `CoreError` variant for Item
  D (decide after reading `src/core/error.rs` in the plan) — flagged, not hidden. ✓
- **Scope:** single coherent PR on the call/contract path; cross-cutting items explicitly deferred
  with severities (§7). ✓
- **Determinism:** Items A+B are designed together precisely because the cap reorders the active set
  and the f64 sum is order-sensitive; the golden/determinism tests guard byte-identity. ✓
