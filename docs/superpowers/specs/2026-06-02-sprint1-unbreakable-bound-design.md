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

### 3.1 The governor (`src/core/governor.rs`, new) — process-global

RSS is a **process-global** resource (the budget is on the whole process, not one call), so the
governor's cancellation is process-global state, not a value threaded through every call. This keeps
the entire public streaming API — including the `run_bounded_whole_genome` ColumnKit SDK entry point
and every crate-root re-export — **signature-stable**: a breach is observed by a one-line
`governor::checkpoint()?` at the hot-loop sites, with no parameter changes anywhere.

Process-global cancellation state (private to `core::governor`), read on the hot path with a single
relaxed atomic load:

```text
static ARMED: AtomicBool         // is a governor active? (one per process)
static TRIPPED: AtomicBool       // has the budget been breached?
static PEAK_AT_TRIP: AtomicU64   // realized peak at the breach (for the error / receipt)
static BUDGET: AtomicU64         // declared budget in bytes (for the error)

pub fn checkpoint() -> Result<(), CoreError>
  // if TRIPPED (relaxed load) -> Err(CoreError::BudgetExceeded { needed: PEAK_AT_TRIP, budget: BUDGET })
  // else Ok(())   — a single relaxed load when no governor is armed (library callers, record-only)
```

The RAII handle owns the poll thread and arms/disarms the global state:

```text
MemoryGovernor::start(budget_bytes, poll, rss_source: impl Fn() -> u64 + Send + 'static)
    -> Result<MemoryGovernor, GovernorError>
  - ARMED.swap(true) == true  -> Err(AlreadyActive)   (one governor per process)
  - reset TRIPPED=false, PEAK_AT_TRIP=0, BUDGET=budget_bytes
  - spawn ONE thread: loop { let r = rss_source();
                             if r > budget_bytes { PEAK_AT_TRIP=r; TRIPPED=true; break }
                             if stop { break } sleep(poll) }
impl Drop  // stop + join the thread, then ARMED=false, TRIPPED=false (clean for the next run)
```

- `rss_source` is `impl Fn() -> u64 + Send + 'static`. Production passes `|| peak_rss_bytes()`
  (getrusage `ru_maxrss`, a monotonic high-water mark — exactly "did peak ever cross the budget").
  It is injectable purely as the **unit-test seam** (§3.6): a closure capturing a rising counter
  trips the governor with no real memory pressure.
- `poll` defaults to **100 ms**; `run_variants_index`/`run_features` read an override from
  `ROSALIND_GOVERNOR_POLL_MS` (tests use a small value for responsiveness).
- The governor is **constructed only under `--enforce`** (held in a local; dropped at function end).
  Record-only runs never arm it, so `checkpoint()` is always `Ok` — the telemetry in §3.3 is still
  recorded (it is pure measurement, independent of the governor).

### 3.2 Cooperative cancellation (the abort path)

The guard never interrupts the main thread mid-allocation (that would corrupt the
*"you keep the data and the proof it overran"* property). The poll thread sets `TRIPPED`; the bounded
drivers observe it cooperatively at `checkpoint()` calls and return an error that the CLI routes into
the exit-4 machinery.

- **No signature changes.** Insert `governor::checkpoint()?;` at six hot-loop sites (one line each):
  1. **Top of each per-contig loop, before `decode_window_arc`** — `call_germline_whole_genome`
     (`whole_genome.rs:61`), `stream_features_whole_genome` (`features.rs:111`),
     `stream_gvcf_whole_genome` (`gvcf.rs:259`). Catches the *dominant* breach mode — decoding the
     next large contig's reference — **before** the allocation happens.
  2. **Once per column** in each region loop — `call_germline_region_streaming` (`pipeline.rs:31`),
     `stream_features_region` (`features.rs:84`), `stream_gvcf_region` (`gvcf.rs:226`). Catches
     within-contig drift (I/O buffers, allocator growth). `run_bounded_whole_genome` inherits the
     features checks (it delegates to `stream_features_whole_genome`).
  - Each call is a single relaxed atomic load when no governor is armed — negligible, and a no-op for
    every library caller and existing test (no behavior change).
- **The gVCF path is covered too** — `variants --index --gvcf --enforce` is a real path; omitting it
  would leave a governor gap in the very feature being built.
- **On `TRIPPED`** → `checkpoint()` returns `Err(CoreError::BudgetExceeded { needed, budget })` with
  **authoritative numbers** carried in the global state (`needed` = `PEAK_AT_TRIP`, `budget` =
  `BUDGET`), so the driver stays decoupled (it raises the error without holding the budget). The CLI
  catches it (§3.5). The existing (currently unraised) `BudgetExceeded` variant (`core/error.rs:12`)
  is reused unchanged.

### 3.3 Receipt — the measured residual (telemetry)

Recorded in the `params` map of the `variants` **and** `features` receipts, in **all** modes
(record-only and enforce — it is measurement, not enforcement). All values are decimal strings, as
the receipt's existing convention:

| key | value |
|---|---|
| `baseline_rss_bytes` | the baseline measured at `main.rs:1610` (already computed, just not recorded) |
| `rss_residual_bytes` | `peak_rss.saturating_sub(max_working_set).saturating_sub(baseline)` — the real I/O+slack overhead this run incurred |
| `io_rss_overhead_assumed_bytes` | the `PILEUP_IO_RSS_OVERHEAD` constant (assumed-vs-realized, side by side) |
| `governor` | `"enforced"` (ran under `--enforce`, fit) · `"tripped"` (governor aborted the run mid-stream) · `"record-only"` (no `--enforce`) |

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

- The driver call in `run_variants_index` / `run_features` is wrapped so its `Result` is **inspected,
  not `?`-propagated**, and the writer is **flushed unconditionally** (so partial output survives a
  breach):
  - `Ok((max_ws, skips))` → today's path unchanged (post-run peak, verdict, receipt, backstop exit 4).
  - `Err(CoreError::BudgetExceeded { needed, .. })` (governor abort) → record `peak_rss = needed`,
    `governor = "tripped"`, `contract_verdict = "over"`, `max_working_set_bytes = 0` (the run did not
    complete, so no working-set high-water was returned — `0` is the honest sentinel), the residual
    fields best-effort; **write the receipt**, print the loud `VIOLATED` message, `std::process::exit(4)`.
  - `Err(other)` → `bail!` as today.
- This keeps the promise: a live breach yields the **proof** (a `verdict=over`, `governor=tripped`
  receipt) plus whatever partial output was flushed — never a silent OOM.
- The **post-run check stays** as a backstop on the `Ok` path, covering the ≤ `poll` sampling gap
  between the last poll and a breach that the kernel did not preempt.
- **Determinism for fitting runs is preserved.** The poll thread only reads RSS and sets a flag; it
  never touches output bytes, and in a run that fits it never fires (every `checkpoint()` is `Ok`). A
  breach is by nature a nondeterministic abort (which column trips depends on timing) — but a breach
  exits non-zero and is a *failure*, not a reproducible artifact, so no determinism guarantee is
  weakened. Stated explicitly in `docs/determinism.md`.

### 3.6 Testing

1. **Governor abort (integration, deterministic).** The CLI's `rss_source` closure reads an env seam
   `ROSALIND_FORCE_LIVE_RSS_BYTES` (a fixed value standing in for live RSS; falls back to
   `peak_rss_bytes()`). The test runs `variants --index --enforce` with a budget *above* the predicted
   peak (so the exit-3 gate passes) but with `ROSALIND_FORCE_LIVE_RSS_BYTES` set *above* the budget, so
   the governor trips on the first poll; assert: exit 4, `governor=tripped` + `contract_verdict=over`
   in the receipt, and the manifest exists. No real memory pressure → CI-safe and fast. (This seam is
   distinct from `ROSALIND_FORCE_PEAK_RSS_BYTES`, which overrides only the *post-run* peak.)
2. **Governor unit test.** `MemoryGovernor::start` with a closure `rss_source` capturing a rising
   counter: `governor::checkpoint()` returns `Ok` before the counter crosses the budget and
   `Err(BudgetExceeded { needed, budget })` after; `Drop` joins cleanly and disarms (a second
   `start` then succeeds); a below-budget source never trips. A second `start` while one is active
   returns `GovernorError::AlreadyActive`.
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
| `src/core/governor.rs` (new) | Process-global `ARMED`/`TRIPPED`/`PEAK_AT_TRIP`/`BUDGET` statics; `pub fn checkpoint() -> Result<(), CoreError>`; `MemoryGovernor::start(budget, poll, rss_source)` (poll thread) + `Drop` (stop/join/disarm); `GovernorError::AlreadyActive`; unit tests. **No signature changes anywhere else.** |
| `src/core/mod.rs` | `pub mod governor;` + re-export `checkpoint`, `MemoryGovernor`, `GovernorError`. |
| `src/core/budget.rs` | Doc-comment: the prediction-margin reframe (the governor is the safety net, not the constant). No value change. |
| `src/call/whole_genome.rs` | One line: `governor::checkpoint()?;` at the top of the `for c in contigs.iter()` loop (`:61`), before `decode_window_arc`. |
| `src/call/pipeline.rs` | One line: `governor::checkpoint()?;` at the top of the `while let Some(column)` loop in `call_germline_region_streaming` (`:31`). |
| `src/call/features.rs` | Two lines: `checkpoint()?` at the per-contig loop top (`:111`) and the per-column loop top in `stream_features_region` (`:84`). |
| `src/call/gvcf.rs` | Two lines: `checkpoint()?` at the per-contig loop top (`:259`) and the per-column loop top in `stream_gvcf_region` (`:226`). |
| `src/call/columnkit.rs` | No change — `run_bounded_whole_genome` delegates to `stream_features_whole_genome`, inheriting its checks. Public SDK signature **unchanged**. |
| `src/main.rs` | `run_variants_index` + `run_features`: under `--enforce`, hold a `MemoryGovernor` (rss_source = the `ROSALIND_FORCE_LIVE_RSS_BYTES`-aware closure); flush the writer unconditionally and **inspect** the driver `Result` instead of `?`-propagating it; on `Err(BudgetExceeded)` write the `governor=tripped`/`verdict=over` receipt and `exit(4)`; record the four new receipt fields (`baseline_rss_bytes`, `rss_residual_bytes`, `io_rss_overhead_assumed_bytes`, `governor`) on both paths. |
| `tests/plan_enforce.rs` | Governor-abort integration test (env seam); near-saturation soundness; residual-recorded assertions. |
| `CONTRACT.md`, `docs/determinism.md` | Document the live governor + the breach-determinism note. |

## 5. Risks & mitigations

- **Sampling gap.** A breach faster than one `poll` interval could still OOM. *Mitigation:* the
  working set grows gradually with coverage (depth-capped), the dominant breach (reference decode) is
  caught pre-allocation at the contig boundary, and the post-run check remains a backstop. 100 ms is
  conservative for genome-scale runs that take seconds-to-minutes per contig.
- **Thread + determinism.** Covered in §3.5 — the poll thread is read-only w.r.t. output; fitting runs
  never fire.
- **Global mutable state.** The cancellation state is process-global because RSS *is* process-global.
  It is fully encapsulated in `core::governor` behind `checkpoint()` + the RAII `MemoryGovernor`, armed
  only while a governor lives, with an `AlreadyActive` guard against two concurrent governors (the CLI
  is one-job-per-process). This buys a **signature-stable public API** — the ColumnKit SDK and every
  re-export are untouched — which is worth more than avoiding one well-scoped global.
- **`checkpoint()` hot-path cost.** A single relaxed atomic load per column; a no-op (always `Ok`) when
  no governor is armed, so library callers and existing tests pay effectively nothing.
- **getrusage cost.** A cheap syscall every 100 ms is negligible.

## 6. References

- `CONTRACT.md` — the four verbs; this hardens "honor".
- `core/budget.rs:24`, `call/plan.rs` — the constant and the predictor.
- `main.rs:1610` (baseline), `:1623` (exit 3), `:1681` (post-run peak), `:1925` (germline driver),
  `:1651/:1666` (features driver) — the integration points.
- `core/error.rs:12` — the reused `BudgetExceeded` variant.
- `tests/plan_enforce.rs:predicted_peak_rss_upper_bounds_realized_peak` — the test extended in §3.6(3).
