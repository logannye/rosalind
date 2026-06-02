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
    assert!(
        checkpoint().is_ok(),
        "a source below budget must never trip"
    );
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
