# Sprint 1.1 — The Unbreakable Bound: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a process-global runtime memory governor so an `--enforce`'d run fails **loud** (exit 4, output + receipt written) the instant realized peak RSS crosses the declared budget — closing the silent-kernel-OOM window between the pre-run refuse (exit 3) and the post-run check — and record the measured RSS residual in every receipt so the 8 MiB prediction margin becomes evidence-backed.

**Architecture:** A new `src/core/governor.rs` holds process-global cancellation state (RSS *is* process-global) behind a `checkpoint()` function and an RAII `MemoryGovernor` (a background poll thread). The bounded streaming drivers gain a one-line `governor::checkpoint()?` at six hot-loop sites — **no public-API signature changes**, so the ColumnKit SDK stays intact. The CLI (`run_variants_index`/`run_features`) holds a `MemoryGovernor` under `--enforce`, inspects the driver result instead of `?`-propagating it, and routes a governor trip into the existing exit-4 path. Four telemetry fields (`baseline_rss_bytes`, `rss_residual_bytes`, `io_rss_overhead_assumed_bytes`, `governor`) are added to the receipt.

**Tech Stack:** Rust (edition 2021, MSRV 1.72), `std::sync::atomic` + `std::thread`, the existing `CoreError::BudgetExceeded` variant, `getrusage`-based `peak_rss_bytes()`.

**Spec:** `docs/superpowers/specs/2026-06-02-sprint1-unbreakable-bound-design.md`

---

### Task 1: The governor module (`core::governor`)

**Files:**
- Create: `src/core/governor.rs`
- Modify: `src/core/mod.rs`
- Test: `tests/governor.rs` (new — a separate test binary so arming the process-global state cannot contaminate the in-process streaming unit tests)

- [ ] **Step 1: Write the failing library test** (`tests/governor.rs`)

```rust
//! Library-API tests for the process-global memory governor. These run in their
//! OWN test binary (separate process) so arming the global breach state cannot
//! contaminate the in-process streaming unit tests; a module-local lock serializes
//! the (process-global) governor tests against each other.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use rosalind::core::governor::{checkpoint, GovernorError, MemoryGovernor};
use rosalind::core::CoreError;

// The governor state is process-global; serialize these tests against each other.
static LOCK: Mutex<()> = Mutex::new(());

#[test]
fn checkpoint_is_ok_when_no_governor_is_armed() {
    let _l = LOCK.lock().unwrap();
    assert!(checkpoint().is_ok());
}

#[test]
fn governor_trips_when_rss_crosses_budget_and_disarms_on_drop() {
    let _l = LOCK.lock().unwrap();
    // rss rises 100 bytes per call; budget 250 -> trips on the 3rd sample (300).
    let counter = Arc::new(AtomicU64::new(0));
    let c = Arc::clone(&counter);
    let gov = MemoryGovernor::start(250, Duration::from_millis(1), move || {
        c.fetch_add(100, Ordering::SeqCst) + 100
    })
    .expect("start");

    // Spin (bounded) until checkpoint observes the breach.
    let mut tripped = false;
    for _ in 0..2000 {
        if checkpoint().is_err() {
            tripped = true;
            break;
        }
        std::thread::sleep(Duration::from_millis(1));
    }
    assert!(tripped, "checkpoint should observe the breach");
    match checkpoint() {
        Err(CoreError::BudgetExceeded { needed, budget }) => {
            assert!(needed > 250, "needed {needed} should exceed the budget");
            assert_eq!(budget, 250);
        }
        other => panic!("expected BudgetExceeded, got {other:?}"),
    }

    drop(gov);
    // Disarmed on drop -> checkpoint is Ok again, and a fresh governor can start.
    assert!(checkpoint().is_ok(), "drop must disarm the global state");
    let again = MemoryGovernor::start(1000, Duration::from_millis(1), || 0).expect("re-start");
    drop(again);
}

#[test]
fn a_below_budget_source_never_trips() {
    let _l = LOCK.lock().unwrap();
    let gov = MemoryGovernor::start(1_000_000, Duration::from_millis(1), || 100).expect("start");
    std::thread::sleep(Duration::from_millis(20));
    assert!(checkpoint().is_ok(), "a source below budget must never trip");
    drop(gov);
}

#[test]
fn a_second_governor_while_one_is_active_is_rejected() {
    let _l = LOCK.lock().unwrap();
    let first = MemoryGovernor::start(1_000_000, Duration::from_millis(50), || 0).expect("first");
    match MemoryGovernor::start(1_000_000, Duration::from_millis(50), || 0) {
        Err(GovernorError::AlreadyActive) => {}
        _ => panic!("a second concurrent governor must be rejected"),
    }
    drop(first);
}
```

- [ ] **Step 2: Run the test to verify it fails to compile**

Run: `cargo test --test governor`
Expected: FAIL — `unresolved import rosalind::core::governor` (the module does not exist yet).

- [ ] **Step 3: Create `src/core/governor.rs`**

```rust
//! Process-global runtime memory governor.
//!
//! Closes the window between the pre-run refuse (exit 3) and the post-run check
//! where an `--enforce`'d run that mispredicts its peak could be silently
//! OOM-killed by the kernel. While a [`MemoryGovernor`] is alive, a background
//! thread polls process RSS; the first sample above the declared budget arms a
//! process-global breach flag that the bounded streaming drivers observe via
//! [`checkpoint`] and turn into a loud `Err(CoreError::BudgetExceeded)`.
//!
//! The state is process-global because RSS *is* a process-global resource — the
//! budget is on the whole process, not one call. This keeps the public streaming
//! API (the ColumnKit SDK) signature-stable: a driver checks the budget with a
//! one-line [`checkpoint`] call, not a threaded parameter.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::Duration;

use crate::core::error::CoreError;

/// Is a governor active in this process? Guards against two concurrent governors.
static ARMED: AtomicBool = AtomicBool::new(false);
/// Has the active governor observed a breach?
static TRIPPED: AtomicBool = AtomicBool::new(false);
/// Realized peak RSS (bytes) at the breach — carried into the error + receipt.
static PEAK_AT_TRIP: AtomicU64 = AtomicU64::new(0);
/// Declared budget (bytes) of the active governor — carried into the error.
static BUDGET: AtomicU64 = AtomicU64::new(0);

/// A cooperative cancellation point for the bounded streaming drivers. Returns
/// `Err(CoreError::BudgetExceeded)` once the active governor has observed a breach,
/// else `Ok(())`. A single relaxed atomic load on the hot path; a no-op (always
/// `Ok`) when no governor is armed (library callers, record-only runs).
pub fn checkpoint() -> Result<(), CoreError> {
    if TRIPPED.load(Ordering::Relaxed) {
        Err(CoreError::BudgetExceeded {
            needed: PEAK_AT_TRIP.load(Ordering::Acquire),
            budget: BUDGET.load(Ordering::Acquire),
        })
    } else {
        Ok(())
    }
}

/// Failure starting a [`MemoryGovernor`].
#[derive(Debug)]
pub enum GovernorError {
    /// A governor is already active in this process (the CLI is one-job-per-process).
    AlreadyActive,
}

impl std::fmt::Display for GovernorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GovernorError::AlreadyActive => write!(f, "a memory governor is already active"),
        }
    }
}

impl std::error::Error for GovernorError {}

/// An RAII runtime memory governor. While alive, a background thread polls
/// `rss_source` every `poll`; the first sample above `budget_bytes` arms the
/// process-global breach state [`checkpoint`] observes. Dropping it stops the
/// thread and disarms the state.
pub struct MemoryGovernor {
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl MemoryGovernor {
    /// Start a governor for `budget_bytes`, polling `rss_source` every `poll`.
    /// Performs one synchronous check before spawning (so an already-breached
    /// source trips deterministically, before the first `checkpoint`). Errors if a
    /// governor is already active in this process.
    pub fn start<F>(budget_bytes: u64, poll: Duration, rss_source: F) -> Result<Self, GovernorError>
    where
        F: Fn() -> u64 + Send + 'static,
    {
        if ARMED.swap(true, Ordering::AcqRel) {
            return Err(GovernorError::AlreadyActive);
        }
        TRIPPED.store(false, Ordering::Release);
        PEAK_AT_TRIP.store(0, Ordering::Release);
        BUDGET.store(budget_bytes, Ordering::Release);

        // Synchronous initial check: a source already above budget trips now, so
        // the very first driver checkpoint fails (no race with the poll thread).
        let initial = rss_source();
        if initial > budget_bytes {
            PEAK_AT_TRIP.store(initial, Ordering::Release);
            TRIPPED.store(true, Ordering::Release);
        }

        let stop = Arc::new(AtomicBool::new(false));
        let stop_thread = Arc::clone(&stop);
        let handle = std::thread::spawn(move || loop {
            let rss = rss_source();
            if rss > budget_bytes {
                PEAK_AT_TRIP.store(rss, Ordering::Release);
                TRIPPED.store(true, Ordering::Release);
                break;
            }
            if stop_thread.load(Ordering::Acquire) {
                break;
            }
            std::thread::sleep(poll);
        });
        Ok(Self {
            stop,
            handle: Some(handle),
        })
    }
}

impl Drop for MemoryGovernor {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Release);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
        TRIPPED.store(false, Ordering::Release);
        PEAK_AT_TRIP.store(0, Ordering::Release);
        ARMED.store(false, Ordering::Release);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Non-arming check only (arming tests live in tests/governor.rs, a separate
    // process, so they cannot trip the global flag during in-process unit tests).
    #[test]
    fn checkpoint_ok_when_unarmed_and_governor_error_displays() {
        assert!(checkpoint().is_ok());
        assert!(GovernorError::AlreadyActive.to_string().contains("already active"));
    }
}
```

- [ ] **Step 4: Wire the module into `src/core/mod.rs`**

Find the module declarations and the re-export block. Add the module declaration alongside the others (e.g. after `pub mod error;`):

```rust
pub mod governor;
```

And add a re-export next to the existing `pub use` lines:

```rust
pub use governor::{checkpoint, GovernorError, MemoryGovernor};
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cargo test --test governor && cargo test --lib governor`
Expected: PASS — all four `tests/governor.rs` tests and the in-crate `checkpoint_ok_when_unarmed_and_governor_error_displays`.

- [ ] **Step 6: Verify no warnings and commit**

Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"`
Expected: `no warnings`

```bash
git add src/core/governor.rs src/core/mod.rs tests/governor.rs
git commit -m "feat(governor): process-global memory governor + checkpoint (core::governor)"
```

---

### Task 2: Insert `checkpoint()` at the six hot-loop sites

No behavior change when no governor is armed (every `checkpoint()` is `Ok`), so the entire existing suite must stay green. The observable effect lands in Task 3.

**Files:**
- Modify: `src/call/whole_genome.rs`, `src/call/pipeline.rs`, `src/call/features.rs`, `src/call/gvcf.rs`

- [ ] **Step 1: `src/call/whole_genome.rs` — per-contig pre-decode check**

Add the import near the top (with the other `use crate::...` lines):

```rust
use crate::core::governor;
```

In `call_germline_whole_genome`, at the very top of the `for c in contigs.iter() {` loop body (before the `let start = ...` decode), insert:

```rust
        // Cooperative budget check (no-op unless an --enforce governor is armed):
        // catch a breach before decoding the next contig's reference.
        governor::checkpoint()?;
```

- [ ] **Step 2: `src/call/pipeline.rs` — per-column germline check**

Add the import near the top:

```rust
use crate::core::governor;
```

In `call_germline_region_streaming`, inside `while let Some(column) = engine.next() {`, immediately after `let column = column?;`, insert:

```rust
        governor::checkpoint()?;
```

- [ ] **Step 3: `src/call/features.rs` — per-contig and per-column checks**

Add the import near the top:

```rust
use crate::core::governor;
```

In `stream_features_whole_genome`, at the top of the `for c in contigs.iter() {` loop body (before `let start = ...`), insert:

```rust
        governor::checkpoint()?;
```

In `stream_features_region`, inside `while let Some(column) = engine.next() {`, immediately after `let column = column?;`, insert:

```rust
        governor::checkpoint()?;
```

- [ ] **Step 4: `src/call/gvcf.rs` — per-contig and per-column checks**

Add the import near the top:

```rust
use crate::core::governor;
```

In `stream_gvcf_whole_genome`, at the top of the `for c in contigs.iter() {` loop body (before `let start = ...`), insert:

```rust
        governor::checkpoint()?;
```

In `stream_gvcf_region`, inside `while let Some(column) = engine.next() {`, immediately after `let column = column?;`, insert:

```rust
        governor::checkpoint()?;
```

- [ ] **Step 5: Run the full suite to verify no regression**

Run: `cargo test 2>&1 | tail -20`
Expected: all tests PASS (the checkpoints are no-ops with no governor armed).

- [ ] **Step 6: Verify no warnings and commit**

Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"`
Expected: `no warnings`

```bash
git add src/call/whole_genome.rs src/call/pipeline.rs src/call/features.rs src/call/gvcf.rs
git commit -m "feat(governor): cooperative checkpoint() at the 6 bounded-driver hot-loop sites"
```

---

### Task 3: CLI governor wiring + breach path in `run_variants_index`

**Files:**
- Modify: `src/main.rs` (`run_variants_index`, ~1798–2106)
- Test: `tests/plan_enforce.rs`

- [ ] **Step 1: Write the failing integration test** (append to `tests/plan_enforce.rs`)

Use this file's existing helpers (verified present): `bin()` (the binary path), `build_sorted_bam_fixture() -> (dir, idx, bam)` (single-contig index + sorted BAM via index/align/sort), and `Command::new(bin())`. With `-o calls.vcf` the receipt sidecar is `calls.vcf.manifest.json`.

```rust
#[test]
fn governor_aborts_loud_when_live_rss_exceeds_budget() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    // 4096 MiB passes the pre-run exit-3 gate, but the live-RSS seam reports
    // 5000 MiB > budget, so the governor trips mid-run -> exit 4.
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .env("ROSALIND_FORCE_LIVE_RSS_BYTES", (5_000u64 * 1024 * 1024).to_string())
        .env("ROSALIND_GOVERNOR_POLL_MS", "1")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(4), "governor breach must exit 4: {out:?}");
    let m = std::fs::read_to_string(&manifest).expect("manifest written on breach");
    assert!(m.contains("\"governor\":\"tripped\""), "manifest: {m}");
    assert!(m.contains("\"contract_verdict\":\"over\""), "manifest: {m}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --test plan_enforce governor_aborts_loud -- --nocapture`
Expected: FAIL — today the run completes (exit 0) and the `ROSALIND_FORCE_LIVE_RSS_BYTES` seam + `governor` field do not exist.

- [ ] **Step 3: Add the governor lifecycle + live-RSS seam in `run_variants_index`**

In `src/main.rs`, add the import inside `run_variants_index`'s `use` block (near `use rosalind::provenance::...`):

```rust
    use rosalind::core::governor::MemoryGovernor;
```

Immediately AFTER the `--enforce` pre-run gate (after the closing `}` of `if enforce { ... }` at ~line 1902, before the `macro_rules! drive!`), insert the governor start + seam:

```rust
    // Live RSS source for the governor: a test seam (ROSALIND_FORCE_LIVE_RSS_BYTES)
    // standing in for live RSS, else the real getrusage high-water mark. Distinct
    // from ROSALIND_FORCE_PEAK_RSS_BYTES, which overrides only the POST-run peak.
    let live_rss = || {
        std::env::var("ROSALIND_FORCE_LIVE_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let poll_ms = std::env::var("ROSALIND_GOVERNOR_POLL_MS")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(100);
    // Under --enforce, the governor fails the run LOUD the moment live RSS crosses
    // the budget (exit 4 with output + receipt) — never a silent kernel OOM. Held
    // for the duration of the calling pass; dropped (thread stopped) at scope end.
    let _governor_guard = if enforce {
        let mb = memory_budget_mb.expect("--enforce requires --memory-budget-mb (checked above)");
        Some(
            MemoryGovernor::start(
                MemoryBudget::from_mb(mb).bytes,
                std::time::Duration::from_millis(poll_ms),
                live_rss,
            )
            .map_err(|e| anyhow!("failed to start memory governor: {e}"))?,
        )
    } else {
        None
    };
```

- [ ] **Step 4: Make `drive!` flush unconditionally and stop `?`-propagating the driver error**

Replace the existing `drive!` macro body so the driver `Result` is returned (not `?`-mapped) and the writer flushes on both paths:

```rust
    macro_rules! drive {
        ($writer:expr) => {{
            let w = $writer;
            let r = if gvcf {
                write_gvcf_header(&mut *w, contigs, "SAMPLE")?;
                stream_gvcf_whole_genome(
                    source,
                    &ref_view,
                    contigs,
                    pileup_params,
                    &germline_params,
                    &mut *w,
                )
            } else {
                write_germline_header(&mut *w, contigs, "SAMPLE")?;
                call_germline_whole_genome(
                    source,
                    &ref_view,
                    contigs,
                    pileup_params,
                    &germline_params,
                    &mut |(locus, ref_base, call)| {
                        write_germline_row(
                            &mut *w,
                            contigs,
                            &GermlineRow {
                                locus,
                                ref_base,
                                call,
                            },
                        )
                        .map_err(rosalind::core::CoreError::from)
                    },
                )
            };
            // Flush even on a governed abort so partial output survives; surface a
            // genuine flush failure only when the calling pass itself succeeded.
            if r.is_ok() {
                w.flush()?;
            } else {
                let _ = w.flush();
            }
            r
        }};
    }
```

- [ ] **Step 5: Inspect the driver result and derive the breach state**

Replace the `let (max_ws, skips) = match &output { ... };` block AND the subsequent `let peak_rss = ...` / `let verdict = ...` block with:

```rust
    let drive_result: Result<(rosalind::core::WorkingSet, rosalind::pileup::SkipCounts), rosalind::core::CoreError> =
        match &output {
            Some(path) => {
                let file = File::create(path)
                    .with_context(|| format!("failed to create VCF file {}", path.display()))?;
                let mut writer = io::BufWriter::new(file);
                drive!(&mut writer)
            }
            None => {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                drive!(&mut handle)
            }
        };

    // A governor trip is the one error we do NOT bail on: we still write the proof
    // receipt (verdict=over, governor=tripped) and exit 4 via the existing post-run
    // check. Any other error is a genuine failure.
    let (max_ws, skips, breached, breach_peak) = match drive_result {
        Ok((ws, sk)) => (ws, sk, false, 0u64),
        Err(rosalind::core::CoreError::BudgetExceeded { needed, .. }) => (
            rosalind::core::WorkingSet { bytes: 0 },
            rosalind::pileup::SkipCounts::default(),
            true,
            needed,
        ),
        Err(e) => return Err(anyhow!("variant calling failed: {e}")),
    };

    // Realized peak: the governor's tripping peak on a breach, else the post-run
    // high-water (ROSALIND_FORCE_PEAK_RSS_BYTES overrides only this post-run value).
    let peak_rss = if breached {
        breach_peak
    } else {
        std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let verdict = match memory_budget_mb.map(|mb| MemoryBudget::from_mb(mb).admits(peak_rss)) {
        None => "unset",
        Some(true) => "within",
        Some(false) => "over",
    };
    let governor_state = if breached {
        "tripped"
    } else if enforce {
        "enforced"
    } else {
        "record-only"
    };
```

(Note: `breached` ⇒ `verdict == "over"` because `breach_peak > budget`; the existing post-run check at the end then prints `VIOLATED` and exits 4 — no new exit branch needed.)

- [ ] **Step 6: Insert the `governor` field into the receipt block**

In the receipt-writing block (`if let Some(dest) = receipt_dest { ... }`), immediately after the existing `max_working_set_bytes` insert (~line 2034), add:

```rust
        manifest
            .params
            .insert("governor".to_string(), governor_state.to_string());
```

- [ ] **Step 7: Run the test to verify it passes, then commit Task 3**

Run: `cargo test --test plan_enforce governor_aborts_loud -- --nocapture && cargo build --release 2>&1 | grep -i warning || echo "no warnings"`
Expected: PASS (exit 4 + `governor=tripped` + `contract_verdict=over` on breach), `no warnings`.

```bash
git add src/main.rs tests/plan_enforce.rs
git commit -m "feat(governor): wire the runtime governor into variants --index (loud exit-4 on breach)"
```

---

### Task 4: Residual telemetry fields in `run_variants_index`

**Files:**
- Modify: `src/main.rs` (`run_variants_index` receipt block, ~2008–2049)
- Test: `tests/plan_enforce.rs`

- [ ] **Step 1: Write the failing residual-recorded test** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn receipt_records_residual_and_governor_fields_on_a_fitting_run() {
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
        .output()
        .unwrap();
    assert!(out.status.success(), "a fitting run should exit 0: {out:?}");

    let text = std::fs::read_to_string(&manifest).expect("manifest");
    for key in [
        "\"baseline_rss_bytes\":",
        "\"rss_residual_bytes\":",
        "\"io_rss_overhead_assumed_bytes\":",
        "\"governor\":\"enforced\"",
    ] {
        assert!(text.contains(key), "manifest missing {key}: {text}");
    }
    // The receipt round-trips through the canonical parser.
    let m = rosalind::provenance::RunManifest::from_canonical_json(&text).expect("parse");
    assert_eq!(m.params.get("governor").map(String::as_str), Some("enforced"));
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --test plan_enforce receipt_records_residual -- --nocapture`
Expected: FAIL — `governor=enforced` is present (from Task 3) but the three residual fields are not in the manifest yet.

- [ ] **Step 3: Derive the residual and add the three telemetry fields**

In `run_variants_index`, add the import for the constant inside the `use` block:

```rust
    use rosalind::core::PILEUP_IO_RSS_OVERHEAD;
```

Just before the receipt-writing block (`if let Some(dest) = receipt_dest {`), derive the residual locals (all in scope: `baseline` from ~line 1874, `peak_rss`/`max_ws` from Task 3):

```rust
    let baseline_rss_bytes = baseline;
    let rss_residual_bytes = peak_rss
        .saturating_sub(max_ws.bytes)
        .saturating_sub(baseline_rss_bytes);
```

In the receipt-writing block, immediately after the `governor` insert added in Task 3, add the three residual fields:

```rust
        manifest.params.insert(
            "baseline_rss_bytes".to_string(),
            baseline_rss_bytes.to_string(),
        );
        manifest.params.insert(
            "rss_residual_bytes".to_string(),
            rss_residual_bytes.to_string(),
        );
        manifest.params.insert(
            "io_rss_overhead_assumed_bytes".to_string(),
            PILEUP_IO_RSS_OVERHEAD.to_string(),
        );
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test --test plan_enforce receipt_records_residual -- --nocapture`
Expected: PASS — the three residual fields + `governor=enforced` are present and the receipt round-trips.

- [ ] **Step 5: Run the full suite + warnings check, then commit Task 4**

Run: `cargo test 2>&1 | tail -15 && cargo build --release 2>&1 | grep -i warning || echo "no warnings"`
Expected: all PASS, `no warnings`.

```bash
git add src/main.rs tests/plan_enforce.rs
git commit -m "feat(governor): record baseline + RSS residual telemetry in the variants receipt"
```

---

### Task 5: Mirror the governor + residual fields in `run_features`

**Files:**
- Modify: `src/main.rs` (`run_features`, ~1557–1796)
- Test: `tests/plan_enforce.rs` (or `tests/features.rs` if that is where feature CLI tests live)

- [ ] **Step 1: Write the failing test** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn features_governor_aborts_loud_and_records_residual() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let tsv = dir.join("feats.tsv");
    let manifest = dir.join("feats.tsv.manifest.json");
    let out = Command::new(bin())
        .args(["features", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&tsv)
        .env("ROSALIND_FORCE_LIVE_RSS_BYTES", (5_000u64 * 1024 * 1024).to_string())
        .env("ROSALIND_GOVERNOR_POLL_MS", "1")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(4), "features breach must exit 4: {out:?}");
    let m = std::fs::read_to_string(&manifest).expect("manifest on breach");
    assert!(m.contains("\"governor\":\"tripped\""), "manifest: {m}");
    assert!(m.contains("\"contract_verdict\":\"over\""), "manifest: {m}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --test plan_enforce features_governor_aborts_loud -- --nocapture`
Expected: FAIL — no seam / `governor` field in `run_features` yet.

- [ ] **Step 3: Add the governor lifecycle in `run_features`**

Add the imports inside `run_features`'s `use` block:

```rust
    use rosalind::core::governor::MemoryGovernor;
    use rosalind::core::PILEUP_IO_RSS_OVERHEAD;
```

Immediately AFTER the `--enforce` pre-run gate (after the closing `}` of `if enforce { ... }` at ~line 1637, before the `let mut analyzer = ...`), insert the SAME governor block as Task 3 Step 3:

```rust
    let live_rss = || {
        std::env::var("ROSALIND_FORCE_LIVE_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let poll_ms = std::env::var("ROSALIND_GOVERNOR_POLL_MS")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(100);
    let _governor_guard = if enforce {
        let mb = memory_budget_mb.expect("--enforce requires --memory-budget-mb (checked above)");
        Some(
            MemoryGovernor::start(
                MemoryBudget::from_mb(mb).bytes,
                std::time::Duration::from_millis(poll_ms),
                live_rss,
            )
            .map_err(|e| anyhow!("failed to start memory governor: {e}"))?,
        )
    } else {
        None
    };
```

- [ ] **Step 4: Inspect the driver result (flush unconditionally) in `run_features`**

Replace the `let (max_ws, skips) = match &output { ... };` block (the two `run_bounded_whole_genome(...)` arms) with an inspecting version that flushes on both paths:

```rust
    let mut analyzer = rosalind::call::FeatureAnalyzer::default();
    // Inline both arms (match arms are exclusive, so moving `source`/`pileup_params`
    // in each is fine — the same pattern the original code used). Flush on BOTH
    // paths so partial output survives a governed abort; inspect the Result rather
    // than `?`-propagating it.
    let drive_result: Result<(rosalind::core::WorkingSet, rosalind::pileup::SkipCounts), rosalind::core::CoreError> =
        match &output {
            Some(path) => {
                let file = File::create(path)
                    .with_context(|| format!("failed to create features file {}", path.display()))?;
                let mut writer = io::BufWriter::new(file);
                let r = run_bounded_whole_genome(
                    &mut analyzer,
                    source,
                    &ref_view,
                    contigs,
                    pileup_params,
                    &mut writer,
                );
                if r.is_ok() {
                    writer.flush()?;
                } else {
                    let _ = writer.flush();
                }
                r
            }
            None => {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                let r = run_bounded_whole_genome(
                    &mut analyzer,
                    source,
                    &ref_view,
                    contigs,
                    pileup_params,
                    &mut handle,
                );
                if r.is_ok() {
                    handle.flush()?;
                } else {
                    let _ = handle.flush();
                }
                r
            }
        };
    let (max_ws, skips, breached, breach_peak) = match drive_result {
        Ok((ws, sk)) => (ws, sk, false, 0u64),
        Err(rosalind::core::CoreError::BudgetExceeded { needed, .. }) => (
            rosalind::core::WorkingSet { bytes: 0 },
            rosalind::pileup::SkipCounts::default(),
            true,
            needed,
        ),
        Err(e) => return Err(anyhow!("feature streaming failed: {e}")),
    };
    let feature_rows = analyzer.rows();
```

Then replace the `let peak_rss = ...` / `let verdict = ...` block with the full derivation (`run_features` is a single task, so it derives the residual locals and inserts all four receipt fields together):

```rust
    let peak_rss = if breached {
        breach_peak
    } else {
        std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let verdict = match memory_budget_mb.map(|mb| MemoryBudget::from_mb(mb).admits(peak_rss)) {
        None => "unset",
        Some(true) => "within",
        Some(false) => "over",
    };
    let governor_state = if breached {
        "tripped"
    } else if enforce {
        "enforced"
    } else {
        "record-only"
    };
    let baseline_rss_bytes = baseline;
    let rss_residual_bytes = peak_rss
        .saturating_sub(max_ws.bytes)
        .saturating_sub(baseline_rss_bytes);
```

- [ ] **Step 5: Add the four receipt fields in `run_features`**

In the `run_features` receipt block, immediately after the existing `max_working_set_bytes` insert (~line 1738), add:

```rust
        manifest.params.insert(
            "baseline_rss_bytes".to_string(),
            baseline_rss_bytes.to_string(),
        );
        manifest.params.insert(
            "rss_residual_bytes".to_string(),
            rss_residual_bytes.to_string(),
        );
        manifest.params.insert(
            "io_rss_overhead_assumed_bytes".to_string(),
            PILEUP_IO_RSS_OVERHEAD.to_string(),
        );
        manifest
            .params
            .insert("governor".to_string(), governor_state.to_string());
```

- [ ] **Step 6: Run the test + full suite, then commit**

Run: `cargo test --test plan_enforce features_governor -- --nocapture && cargo test 2>&1 | tail -15`
Expected: PASS. Then `cargo build --release 2>&1 | grep -i warning || echo "no warnings"` → `no warnings`.

```bash
git add src/main.rs tests/plan_enforce.rs
git commit -m "feat(governor): wire the governor + residual telemetry into features --index"
```

---

### Task 6: Near-saturation soundness test

Validate that the 8 MiB prediction margin actually covers the real residual at the budget boundary on a contig large enough that reference-decode dominates — and that the governor does NOT spuriously fire on a run that genuinely fits.

**Files:**
- Test: `tests/plan_enforce.rs`

- [ ] **Step 1: Write the near-saturation soundness test**

Uses the file's existing `build_big_contig_fixture(ref_len) -> (dir, idx, bam)` (the same fixture the 4 MiB `predicted_peak_rss_upper_bounds_realized_peak` test uses) and the canonical-manifest parser. Single-process (no cross-process baseline jitter), and it validates the new `rss_residual_bytes` field directly against the assumed margin.

```rust
#[test]
fn near_saturation_margin_holds_on_a_larger_contig() {
    // Soundness of the 8 MiB prediction margin at a bigger scale than the 4 MiB
    // sibling test: on an ~8 MB contig (reference-decode dominates), the predicted
    // peak must still upper-bound the realized peak — the fixed I/O+slack margin
    // covers the real residual — AND the recorded residual is within that margin.
    // (Governor no-false-fire on fitting runs is covered by the existing --enforce
    // tests, which now run with the governor armed and still exit 0.)
    let (dir, idx, bam) = build_big_contig_fixture(8 * 1024 * 1024);
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

    let m = rosalind::provenance::RunManifest::from_canonical_json(
        &std::fs::read_to_string(&manifest).unwrap(),
    )
    .unwrap();
    let predicted: u64 = m.params.get("predicted_peak_rss_bytes").unwrap().parse().unwrap();
    let realized: u64 = m.params.get("peak_rss_bytes").unwrap().parse().unwrap();
    let residual: u64 = m.params.get("rss_residual_bytes").unwrap().parse().unwrap();
    let assumed: u64 = m.params.get("io_rss_overhead_assumed_bytes").unwrap().parse().unwrap();
    assert!(
        predicted >= realized,
        "8 MB contig: predicted {predicted} must be >= realized {realized}"
    );
    assert!(
        residual <= assumed,
        "recorded residual {residual} should be within the assumed margin {assumed} \
         (if not, the 8 MiB constant is too small for this workload — a REAL finding, \
          not a test to silence)"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run the test**

Run: `cargo test --test plan_enforce near_saturation_margin -- --nocapture`
Expected: PASS — predicted ≥ realized and the recorded residual is within the assumed margin. (If it FAILS, that is a real finding — the 8 MiB margin is too small for this workload; record it and raise the issue rather than silencing the test.)

- [ ] **Step 3: Commit**

```bash
git add tests/plan_enforce.rs
git commit -m "test(governor): near-saturation soundness — budget == predicted peak still fits"
```

---

### Task 7: Documentation

**Files:**
- Modify: `src/core/budget.rs` (doc-comment), `CONTRACT.md`, `docs/determinism.md`

- [ ] **Step 1: `src/core/budget.rs` — reframe the constant's role**

Extend the doc-comment on `PILEUP_IO_RSS_OVERHEAD` (currently ~lines 16–24) with a sentence on the governor reframe (no value change):

```rust
/// ...existing text...
///
/// Reframe (Sprint 1.1): this margin's correctness is no longer load-bearing for
/// *safety*. Under `--enforce` a runtime `MemoryGovernor` (`core::governor`) fails
/// the run loud (exit 4, output + receipt written) the moment realized peak RSS
/// crosses the budget, so an under-prediction here is caught live rather than
/// silently OOM-killed. The realized residual (`peak − working_set − baseline`) is
/// now recorded in every receipt (`rss_residual_bytes`) to make a future, evidence-
/// based re-tuning of this constant possible.
```

- [ ] **Step 2: `CONTRACT.md` — document the live governor in the "Honor" section**

In the `### 3. Honor` section, after the bullet describing `realized peak > budget → fail loud (exit 4)`, add:

```markdown
- **runtime governor** — under `--enforce`, a background guard polls process RSS and
  fails the run **loud at the moment of breach** (exit 4, partial output + a
  `governor=tripped`, `contract_verdict=over` receipt), so a misprediction is caught
  *during* the run, not by a silent kernel OOM. The receipt also records the realized
  RSS residual (`rss_residual_bytes`) against the assumed margin
  (`io_rss_overhead_assumed_bytes`).
```

- [ ] **Step 3: `docs/determinism.md` — the breach-determinism note**

Add a short paragraph:

```markdown
## The memory governor and determinism

Under `--enforce`, a background governor thread reads process RSS and, on a breach,
cooperatively aborts the run (exit 4). This does **not** weaken determinism: the
guard never touches output bytes, and a run that fits never fires it — identical
inputs still produce a byte-identical VCF/feature table. A *breach* is a
non-deterministic abort (which locus trips depends on timing), but a breach exits
non-zero and is a failure, not a reproducible artifact.
```

- [ ] **Step 4: Verify the build + full suite, then commit**

Run: `cargo test 2>&1 | tail -15 && cargo build --release 2>&1 | grep -i warning || echo "no warnings"`
Expected: all PASS, `no warnings`.

```bash
git add src/core/budget.rs CONTRACT.md docs/determinism.md
git commit -m "docs(governor): document the runtime governor + residual telemetry + determinism"
```

---

## Final verification

- [ ] **Full suite green:** `cargo test 2>&1 | tail -20` — all sections PASS.
- [ ] **Zero warnings:** `cargo build --release 2>&1 | grep -i warning || echo "no warnings"` and `cargo build --release --examples 2>&1 | grep -i warning || echo "no warnings"`.
- [ ] **Format clean:** `cargo fmt --check`.
- [ ] **Smoke the contract end-to-end** on the bundled fixture: `cargo run --release -- variants --index <toy.idx> --alignments <toy.sorted.bam> --memory-budget-mb 512 --enforce -o /tmp/c.vcf` → `contract: OK`, and the manifest contains `governor`, `baseline_rss_bytes`, `rss_residual_bytes`, `io_rss_overhead_assumed_bytes`.
- [ ] **ColumnKit SDK signature unchanged:** confirm `run_bounded_whole_genome` and the re-exports in `src/lib.rs` / `src/call/mod.rs` are byte-for-byte the same as before (no parameter added).
- [ ] CI green on GitHub before reporting done (local-green ≠ CI-green).
