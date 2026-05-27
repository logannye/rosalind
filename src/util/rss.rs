//! RSS / peak-memory measurement utilities.
//!
//! This module provides a **real** process RSS signal to complement the logical
//! `SpaceTracker` counters used in other parts of the codebase.

use libc::rusage;

/// Return the process peak RSS in bytes (best-effort, platform-dependent).
///
/// - On Linux, `ru_maxrss` is reported in KiB.
/// - On macOS, `ru_maxrss` is reported in bytes.
pub fn peak_rss_bytes() -> u64 {
    let mut usage: rusage = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage as *mut rusage) };
    if rc != 0 {
        return 0;
    }
    let raw = usage.ru_maxrss as u64;

    #[cfg(target_os = "linux")]
    {
        raw.saturating_mul(1024)
    }
    #[cfg(not(target_os = "linux"))]
    {
        raw
    }
}
