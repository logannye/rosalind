//! Cooperative cancellation shared by native hashing, decoding and batch drivers.
//!
//! Scope installation and cancellation take a short mutex; hot checkpoints use
//! only atomic loads. A scope is process-wide, like the memory governor. Separate
//! simultaneous jobs should use separate native processes. Native C calls remain
//! cooperative: cancellation is observed when control returns to a checkpoint.

use super::CoreError;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

static ACTIVE: Mutex<Option<Arc<AtomicBool>>> = Mutex::new(None);
static GENERATION: AtomicU64 = AtomicU64::new(0);
static NEXT_GENERATION: AtomicU64 = AtomicU64::new(1);
static SIGNAL_GENERATION: AtomicU64 = AtomicU64::new(0);
static SIGNAL_HANDLERS_INSTALLED: AtomicBool = AtomicBool::new(false);
static REQUESTED: AtomicBool = AtomicBool::new(false);

/// Cloneable one-way cancellation request. Cancelling an inactive or old token
/// cannot cancel a different job. Tokens may be cancelled before a scope starts.
#[derive(Debug, Clone, Default)]
pub struct CancellationToken {
    cancelled: Arc<AtomicBool>,
}
impl CancellationToken {
    /// Create an uncancelled token.
    pub fn new() -> Self {
        Self::default()
    }
    /// Request cancellation from an ordinary application thread.
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Release);
        let active = ACTIVE.lock().unwrap_or_else(|error| error.into_inner());
        if active
            .as_ref()
            .is_some_and(|token| Arc::ptr_eq(token, &self.cancelled))
        {
            REQUESTED.store(true, Ordering::Release);
        }
    }
    /// Whether this token, including a signal delivered to its active scope, was cancelled.
    pub fn is_cancelled(&self) -> bool {
        if self.cancelled.load(Ordering::Acquire) {
            return true;
        }
        let active = ACTIVE.lock().unwrap_or_else(|error| error.into_inner());
        active
            .as_ref()
            .is_some_and(|token| Arc::ptr_eq(token, &self.cancelled))
            && cancellation_requested()
    }
    /// Return a typed cancellation error when cancellation has been requested.
    pub fn check(&self) -> Result<(), CoreError> {
        if self.is_cancelled() {
            Err(CoreError::Cancelled)
        } else {
            Ok(())
        }
    }
}

/// Another process-wide cancellation scope is already active.
#[derive(Debug, thiserror::Error)]
#[error("a process-wide cancellation scope is already active")]
pub struct CancellationScopeActive;

/// Activates a token for shared checkpoints and clears that binding on drop.
/// Only one scope may exist per process. It does not install signal handlers.
#[derive(Debug)]
pub struct CancellationScope {
    token: CancellationToken,
    generation: u64,
}
impl CancellationScope {
    /// Bind a token before input validation/hashing. A pre-cancelled token is
    /// observed by the first checkpoint; installation itself remains reversible.
    pub fn start(token: CancellationToken) -> Result<Self, CancellationScopeActive> {
        let mut active = ACTIVE.lock().unwrap_or_else(|error| error.into_inner());
        if active.is_some() {
            return Err(CancellationScopeActive);
        }
        let generation = NEXT_GENERATION
            .fetch_update(Ordering::AcqRel, Ordering::Acquire, |value| {
                value.checked_add(1)
            })
            .map_err(|_| CancellationScopeActive)?;
        REQUESTED.store(token.cancelled.load(Ordering::Acquire), Ordering::Release);
        *active = Some(Arc::clone(&token.cancelled));
        GENERATION.store(generation, Ordering::Release);
        Ok(Self { token, generation })
    }
}
impl Drop for CancellationScope {
    fn drop(&mut self) {
        let mut active = ACTIVE.lock().unwrap_or_else(|error| error.into_inner());
        if active
            .as_ref()
            .is_some_and(|token| Arc::ptr_eq(token, &self.token.cancelled))
        {
            if cancellation_requested() {
                self.token.cancelled.store(true, Ordering::Release);
            }
            GENERATION.store(0, Ordering::Release);
            *active = None;
            REQUESTED.store(false, Ordering::Release);
        }
    }
}

pub(crate) fn checkpoint() -> Result<(), CoreError> {
    if cancellation_requested() {
        Err(CoreError::Cancelled)
    } else {
        Ok(())
    }
}

fn cancellation_requested() -> bool {
    let generation = GENERATION.load(Ordering::Acquire);
    generation != 0
        && (REQUESTED.load(Ordering::Acquire)
            || SIGNAL_GENERATION.load(Ordering::Acquire) == generation)
}

/// Opt-in SIGINT/SIGTERM handling for a native CLI's active cancellation scope.
/// Existing process handlers are restored on drop. Library runners should leave
/// this opt-in disabled unless the embedding application delegates signal ownership.
pub struct SignalCancellationGuard<'scope> {
    scope: std::marker::PhantomData<&'scope CancellationScope>,
    #[cfg(unix)]
    previous: [(libc::c_int, libc::sigaction); 2],
}
impl std::fmt::Debug for SignalCancellationGuard<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SignalCancellationGuard")
            .finish_non_exhaustive()
    }
}
#[cfg(unix)]
extern "C" fn cancel_signal(_: libc::c_int) {
    // AtomicBool is lock-free on the supported Linux/macOS targets. Do not lock,
    // allocate, log or call into an application callback from a signal handler.
    let generation = GENERATION.load(Ordering::Acquire);
    if generation != 0 {
        SIGNAL_GENERATION.fetch_max(generation, Ordering::AcqRel);
    }
}
impl<'scope> SignalCancellationGuard<'scope> {
    /// Install handlers for an active scope; unsupported platforms or a second
    /// live guard refuse. The guard cannot outlive its scope.
    ///
    /// ```compile_fail
    /// use rosalind::core::cancellation::*;
    /// let scope = CancellationScope::start(CancellationToken::new()).unwrap();
    /// let signals = SignalCancellationGuard::install(&scope).unwrap();
    /// drop(scope);
    /// drop(signals);
    /// ```
    pub fn install(scope: &'scope CancellationScope) -> std::io::Result<Self> {
        if GENERATION.load(Ordering::Acquire) != scope.generation {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "signal cancellation requires an active cancellation scope",
            ));
        }
        #[cfg(unix)]
        {
            if SIGNAL_HANDLERS_INSTALLED
                .compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
                .is_err()
            {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::AlreadyExists,
                    "cancellation signal handlers are already installed",
                ));
            }
            // SAFETY: zero is a valid starting representation for sigaction;
            // sigemptyset initializes its mask before it is installed.
            let mut action: libc::sigaction = unsafe { std::mem::zeroed() };
            action.sa_sigaction = cancel_signal as *const () as libc::sighandler_t;
            action.sa_flags = libc::SA_RESTART;
            // SAFETY: the mask points into our live, writable sigaction structure.
            unsafe {
                libc::sigemptyset(&mut action.sa_mask);
            }
            // SAFETY: output storage is initialized before any restoration call.
            let mut previous: [(libc::c_int, libc::sigaction); 2] = unsafe { std::mem::zeroed() };
            for (index, signal) in [libc::SIGINT, libc::SIGTERM].into_iter().enumerate() {
                previous[index].0 = signal;
                // SAFETY: all pointers remain valid for this synchronous libc call.
                if unsafe { libc::sigaction(signal, &action, &mut previous[index].1) } != 0 {
                    let error = std::io::Error::last_os_error();
                    for (signal, old) in &previous[..index] {
                        // SAFETY: restore only successfully initialized prior actions.
                        unsafe {
                            libc::sigaction(*signal, old, std::ptr::null_mut());
                        }
                    }
                    SIGNAL_HANDLERS_INSTALLED.store(false, Ordering::Release);
                    return Err(error);
                }
            }
            Ok(Self {
                previous,
                scope: std::marker::PhantomData,
            })
        }
        #[cfg(not(unix))]
        {
            Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "signal cancellation requires Unix",
            ))
        }
    }
}
impl Drop for SignalCancellationGuard<'_> {
    fn drop(&mut self) {
        #[cfg(unix)]
        {
            for (signal, old) in &self.previous {
                // SAFETY: restore the valid actions captured during installation.
                unsafe {
                    libc::sigaction(*signal, old, std::ptr::null_mut());
                }
            }
            SIGNAL_HANDLERS_INSTALLED.store(false, Ordering::Release);
        }
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    #[test]
    fn delayed_old_signal_cannot_cancel_a_new_job_or_erase_its_cancellation() {
        const NAME: &str = "core::cancellation::tests::delayed_old_signal_cannot_cancel_a_new_job_or_erase_its_cancellation";
        if std::env::var_os("ROSALIND_SIGNAL_GENERATION_CHILD").is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", NAME, "--nocapture"])
                .env("ROSALIND_SIGNAL_GENERATION_CHILD", "1")
                .status()
                .unwrap();
            assert!(status.success());
            return;
        }
        let first = CancellationScope::start(CancellationToken::new()).unwrap();
        let delayed = first.generation;
        drop(first);
        let next = CancellationScope::start(CancellationToken::new()).unwrap();
        // Model a handler paused after reading the previous scope's generation.
        SIGNAL_GENERATION.fetch_max(delayed, Ordering::AcqRel);
        assert!(checkpoint().is_ok());
        cancel_signal(libc::SIGINT);
        SIGNAL_GENERATION.fetch_max(delayed, Ordering::AcqRel);
        assert!(matches!(checkpoint(), Err(CoreError::Cancelled)));
        assert!(next.token.is_cancelled());
    }
}
