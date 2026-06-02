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
#[derive(Debug)]
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
        assert!(GovernorError::AlreadyActive
            .to_string()
            .contains("already active"));
    }
}
