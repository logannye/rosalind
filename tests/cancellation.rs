//! Process-wide cancellation tests run in isolated child invocations.
use rosalind::core::cancellation::{CancellationScope, CancellationToken, SignalCancellationGuard};
use rosalind::core::governor::checkpoint;
use rosalind::core::CoreError;

fn child(name: &str) -> bool {
    if std::env::var_os("ROSALIND_CANCELLATION_CHILD").is_some() {
        return true;
    }
    let status = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", name, "--nocapture"])
        .env("ROSALIND_CANCELLATION_CHILD", "1")
        .status()
        .unwrap();
    assert!(status.success());
    false
}
#[test]
fn token_scope_cancels_workers_and_old_tokens_cannot_cancel_a_new_job() {
    if !child("token_scope_cancels_workers_and_old_tokens_cannot_cancel_a_new_job") {
        return;
    }
    let token = CancellationToken::new();
    let scope = CancellationScope::start(token.clone()).unwrap();
    assert!(checkpoint().is_ok());
    assert!(CancellationScope::start(CancellationToken::new()).is_err());
    let cloned = token.clone();
    std::thread::spawn(move || cloned.cancel()).join().unwrap();
    assert!(matches!(checkpoint(), Err(CoreError::Cancelled)));
    assert!(token.is_cancelled());
    drop(scope);
    assert!(checkpoint().is_ok());
    let next = CancellationToken::new();
    let _scope = CancellationScope::start(next.clone()).unwrap();
    token.cancel();
    assert!(checkpoint().is_ok());
    assert!(!next.is_cancelled());
}
#[test]
fn pre_cancelled_tokens_are_observed_before_any_work_and_drop_disarms() {
    if !child("pre_cancelled_tokens_are_observed_before_any_work_and_drop_disarms") {
        return;
    }
    let token = CancellationToken::new();
    token.cancel();
    let scope = CancellationScope::start(token).unwrap();
    assert!(matches!(checkpoint(), Err(CoreError::Cancelled)));
    drop(scope);
    assert!(checkpoint().is_ok());
}
#[cfg(unix)]
#[test]
fn opt_in_signals_cancel_the_active_scope_and_restore_existing_handlers() {
    if !child("opt_in_signals_cancel_the_active_scope_and_restore_existing_handlers") {
        return;
    }
    for signal in [libc::SIGINT, libc::SIGTERM] {
        // SAFETY: sigaction writes valid prior-handler state into this live object.
        let mut previous: libc::sigaction = unsafe { std::mem::zeroed() };
        assert_eq!(
            unsafe { libc::sigaction(signal, std::ptr::null(), &mut previous) },
            0
        );
        let token = CancellationToken::new();
        let scope = CancellationScope::start(token.clone()).unwrap();
        let guard = SignalCancellationGuard::install(&scope).unwrap();
        assert!(SignalCancellationGuard::install(&scope).is_err());
        // SAFETY: this isolated child installed its own cancellation handler.
        assert_eq!(unsafe { libc::raise(signal) }, 0);
        assert!(matches!(checkpoint(), Err(CoreError::Cancelled)));
        assert!(token.is_cancelled());
        drop(guard);
        let mut restored: libc::sigaction = unsafe { std::mem::zeroed() };
        assert_eq!(
            unsafe { libc::sigaction(signal, std::ptr::null(), &mut restored) },
            0
        );
        assert_eq!(previous.sa_sigaction, restored.sa_sigaction);
        drop(scope);
        assert!(token.is_cancelled());
        assert!(checkpoint().is_ok());
    }
}
