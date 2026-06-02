# Sprint 1.1 — The Unbreakable Bound (design)

**Status:** Approved design — 2026-06-02. Increment 1.1 of the engineering roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md), Sprint 1). Companion to [`CONTRACT.md`](../../../CONTRACT.md).

---

## 1. Problem

Rosalind's headline promise is *"never a silent OOM — it fits, or it tells you up front, and it
proves the realized peak with a receipt."* Today that promise has a hole and an unproven assumption:

1. **The silent-OOM window.** Under `--enforce`, the budget is checked **before** the run
   (`main.rs:1623`, refuse exit 3) and **after** the run (`main.rs:1681`, fail-loud exit 4). Between
   them, the BAM streams (`call_germline_whole_genome` / `run_bounded_whole_genome`) with **no live
   guard**. If the realized peak crosses the budget mid-run — e.g. the prediction was wrong — the
   kernel OOM-killer can fire *before* the post-run check is ever reached. The one failure mode the
   thesis explicitly forbids is currently reachable.

2. **The unproven margin.** The entire "predicted peak is an upper bound" claim rests on
   `PILEUP_IO_RSS_OVERHEAD = 8 MiB` (`core/budget.rs:24`) — a hard-coded constant standing in for
   {htslib/BGZF buffers, the VCF `BufWriter`, allocator slack}. There is no evidence that
   `realized_peak − working_set − baseline ≤ 8 MiB` across platforms, allocators, and htslib
   versions, and the only upper-bound test (`tests/plan_enforce.rs`, `predicted_peak_rss_upper_bounds_realized_peak`)
   uses a 4 MiB contig where the margin is never stressed near saturation.

This increment closes the window with a runtime governor and starts measuring the real margin, so a
future increment can re-tune the constant with evidence rather than a guess.

## 2. Goals / non-goals

**Goals**
- A runtime **memory governor** that, under `--enforce`, fails **loud** (exit 4, with output +
  receipt written) at the moment a breach is detected — never a silent kernel OOM.
- Record the **measured RSS residual** (`peak − working_set − baseline`) and the assumed margin in
  every receipt, so the 8 MiB constant becomes an evidence-backed decision later.
- A **near-saturation** soundness test that stresses the margin at the budget boundary on a contig
  large enough that reference-decode dominates — not the 4 MiB toy.

**Non-goals (explicitly deferred)**
- Re-tuning / replacing the `PILEUP_IO_RSS_OVERHEAD` value (future, evidence-driven — we record first).
- Receipt self-hash / tamper-evidence / `SCHEMA_VERSION` (Sprint 1.2).
- `cargo publish`, clippy/MSRV gate, doc-drift fixes (Sprint 1.3).
- Any change to the prediction formula or the exit-3 pre-run refuse path.

## 3. Design

### 3.1 The governor (`src/core/governor.rs`, new)

```text
MemoryGovernor::start(budget_bytes: u64, poll: Duration, rss_source: fn() -> u64) -> MemoryGovernor
  - spawns ONE background thread:
      loop { if rss_source() > budget_bytes { tripped.store(true, Release); break }
             if stop.load(Acquire) { break }
             sleep(poll) }
  - holds: Arc<AtomicBool> tripped, Arc<AtomicBool> stop, JoinHandle
MemoryGovernor::tripped() -> &AtomicBool        // the cancel token the drivers observe
MemoryGovernor::stop(self)                       // sets stop, joins the thread
impl Drop                                        // stop-and-join if not already stopped
```

- `rss_source` defaults to `rosalind::util::rss::peak_rss_bytes` (getrusage `ru_maxrss`, a
  monotonic high-water mark — exactly "did peak ever cross the budget"). It is injectable purely as
  the **test seam** (§3.5): a test feeds a rising sequence with no real memory pressure.
- `poll` defaults to **100 ms**; `run_variants_index`/`run_features` read an override from
  `ROSALIND_GOVERNOR_POLL_MS` (tests use a small value for responsiveness).
- The governor is **started only under `--enforce`**. Record-only runs construct **no** governor and
  pass `None` as the cancel token (telemetry in §3.3 is still recorded — it is pure measurement).

### 3.2 Cooperative cancellation (the abort path)

The guard never interrupts the main thread mid-allocation (that would corrupt the
*"you keep the data and the proof it overran"* property). It sets `tripped`; the bounded drivers
observe it cooperatively and return an error that flows through the **existing** exit-4 machinery.

- **Signature change:** add `cancel: Option<&AtomicBool>` to:
  - `call_germline_whole_genome` (`call/whole_genome.rs:49`) and the per-region caller it drives,
    `call_germline_region_streaming` (`call/pipeline.rs:20`).
  - `run_bounded_whole_genome` (`call/columnkit.rs:61`).
  - the features region helpers `stream_features_whole_genome` (`call/features.rs:100`) /
    `stream_features_region` (`call/features.rs:74`).
  - `None` preserves today's behavior; **all existing call sites and tests pass `None`** and are
    unchanged.
- **Two check points** (each a single `Relaxed` atomic load — negligible):
  1. **Top of the per-contig loop, before `decode_window_arc`** (`whole_genome.rs:61`,
     `features.rs:111`, the columnkit driver loop). This catches the *dominant* breach mode — decoding
     the next large contig's reference — **before** the allocation happens.
  2. **Once per column** in the region loop (`call_germline_region_streaming` `pipeline.rs:31`;
     `stream_features_region` `features.rs:84`). Catches within-contig drift (I/O buffers, allocator
     growth).
- **On `tripped`** → the driver returns `Err(CoreError::BudgetExceeded { .. })`, reusing the existing
  (currently unraised) variant (`core/error.rs:12`). The driver needs to know only *that* it was
  cancelled, not the budget; `run_variants_index` / `run_features` **catch the error and populate the
  authoritative numbers** for the loud message + receipt — `budget` = the declared budget, `needed` =
  the post-trip `peak_rss_bytes()`. (Field population by the caller keeps the driver decoupled from
  the budget it does not otherwise hold.)

### 3.3 Receipt — the measured residual (telemetry)

Recorded in the `params` map of the `variants` **and** `features` receipts, in **all** modes
(record-only and enforce — it is measurement, not enforcement). All values are decimal strings, as
the receipt's existing convention:

| key | value |
|---|---|
| `baseline_rss_bytes` | the baseline measured at `main.rs:1610` (already computed, just not recorded) |
| `rss_residual_bytes` | `peak_rss.saturating_sub(max_working_set).saturating_sub(baseline)` — the real I/O+slack overhead this run incurred |
| `io_rss_overhead_assumed_bytes` | the `PILEUP_IO_RSS_OVERHEAD` constant (assumed-vs-realized, side by side) |
| `governor` | `"enforced"` when `--enforce`, else `"record-only"` |

> Determinism note: `baseline_rss_bytes` and `rss_residual_bytes` are machine-dependent
> measurements, exactly like the already-present `peak_rss_bytes`. The VCF/feature output stays
> byte-identical; the manifest differs only in these measured fields (the same caveat
> `CONTRACT.md` already documents for `peak_rss_bytes`). `verify`'s existing internal-consistency
> checks are unaffected; the new fields are recorded, not re-derived.

### 3.4 The constant + prediction reframe (no value change)

`PILEUP_IO_RSS_OVERHEAD` **stays 8 MiB**. You cannot measure a run before running it, so the up-front
prediction keeps a fixed conservative margin. What changes is that **its correctness is no longer
load-bearing for safety**: if it ever under-predicts, the governor catches the live breach and fails
loud (exit 4) instead of a silent kernel OOM. We deliberately do **not** re-tune it here — we record
`rss_residual_bytes` first and tune later with data. This reframe is documented in the `budget.rs`
doc-comment and the `CONTRACT.md` "honor" section.

### 3.5 Error handling, exit codes, determinism

- `Err(CoreError::BudgetExceeded)` propagates to `run_variants_index` / `run_features`, which route it
  into the **existing exit-4 path**: write the partial output, write the receipt with
  `contract_verdict = "over"`, print the loud breach message, `std::process::exit(4)`. A live breach
  now yields the *same* artifacts a post-run breach does.
- The **post-run check stays** as a backstop, covering the ≤ `poll` sampling gap between the last
  poll and a breach.
- **Determinism for fitting runs is preserved.** The guard only reads RSS and sets a flag; it never
  touches output bytes, and in a run that fits it never fires. A breach is by nature a
  nondeterministic abort (which column trips depends on timing) — but a breach exits non-zero and is
  a *failure*, not a reproducible artifact, so no determinism guarantee is weakened. Stated
  explicitly in `docs/determinism.md`.

### 3.6 Testing

1. **Governor abort (integration, deterministic).** Inject an `rss_source` that returns a rising
   sequence crossing the budget after a few polls; run a small `variants --index --enforce`; assert:
   exit 4, `contract_verdict=over` in the receipt, and the partial VCF **and** manifest exist. No real
   memory pressure → CI-safe and fast.
2. **Governor unit test.** `MemoryGovernor` with a stub `rss_source` trips `tripped` once the stub
   crosses the budget and not before; `stop()` joins cleanly; below-budget never trips.
3. **Near-saturation soundness.** Extend `predicted_peak_rss_upper_bounds_realized_peak`: build a
   larger synthetic contig (a few MB — large enough that reference-decode dominates), run with
   `budget == predicted_peak` (zero extra headroom), assert the run completes within budget (the
   margin holds at the boundary), and `predicted_peak ≥ realized_peak`.
4. **Residual recorded.** Assert `baseline_rss_bytes`, `rss_residual_bytes`,
   `io_rss_overhead_assumed_bytes`, `governor` appear in the receipt and round-trip through
   `RunManifest::from_canonical_json` / `verify`.
5. **Regression.** Existing exit-3, exit-4 (`ROSALIND_FORCE_PEAK_RSS_BYTES`), and determinism tests
   stay green; all non-enforce call sites compile with `None`.

## 4. File-by-file change list

| File | Change |
|---|---|
| `src/core/governor.rs` (new) | `MemoryGovernor` (poll thread, injectable `rss_source`, `tripped`/`stop`/`Drop`). |
| `src/core/mod.rs` | Export `governor` + `MemoryGovernor`. |
| `src/core/budget.rs` | Doc-comment: the prediction-margin reframe (governor is the safety net, not the constant). No value change. |
| `src/call/whole_genome.rs` | `cancel: Option<&AtomicBool>` param; check at per-contig top (pre-decode); thread into the region caller. |
| `src/call/pipeline.rs` | `cancel` param on `call_germline_region_streaming` (+ `_tracked` wrapper); per-column check in the `while let Some(column)` loop. |
| `src/call/columnkit.rs` | `cancel: Option<&AtomicBool>` param on `run_bounded_whole_genome`; per-contig + per-column checks. |
| `src/call/features.rs` | `cancel` threaded into `stream_features_region` / `stream_features_whole_genome`; per-column check. |
| `src/main.rs` | `run_variants_index` + `run_features`: under `--enforce` start a `MemoryGovernor`, pass `governor.tripped()` as `cancel`, stop on completion; map `Err(BudgetExceeded)` to the exit-4 path; record the four new receipt fields (`baseline_rss_bytes`, `rss_residual_bytes`, `io_rss_overhead_assumed_bytes`, `governor`). |
| `tests/plan_enforce.rs` | Governor-abort integration test; near-saturation soundness; residual-recorded assertions. |
| `CONTRACT.md`, `docs/determinism.md` | Document the live governor + the breach-determinism note. |

## 5. Risks & mitigations

- **Sampling gap.** A breach faster than one `poll` interval could still OOM. *Mitigation:* the
  working set grows gradually with coverage (depth-capped), the dominant breach (reference decode) is
  caught pre-allocation at the contig boundary, and the post-run check remains a backstop. 100 ms is
  conservative for genome-scale runs that take seconds-to-minutes per contig.
- **Thread + determinism.** Covered in §3.5 — guard is read-only w.r.t. output; fitting runs never
  fire.
- **Signature churn.** Mitigated by `Option<&AtomicBool>` defaulting to `None`, so every existing
  caller/test is a one-token change with identical behavior.
- **getrusage cost.** A cheap syscall every 100 ms is negligible.

## 6. References

- `CONTRACT.md` — the four verbs; this hardens "honor".
- `core/budget.rs:24`, `call/plan.rs` — the constant and the predictor.
- `main.rs:1610` (baseline), `:1623` (exit 3), `:1681` (post-run peak), `:1925` (germline driver),
  `:1651/:1666` (features driver) — the integration points.
- `core/error.rs:12` — the reused `BudgetExceeded` variant.
- `tests/plan_enforce.rs:predicted_peak_rss_upper_bounds_realized_peak` — the test extended in §3.6(3).
