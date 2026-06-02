//! RSS / peak-memory measurement utilities.
//!
//! This module provides a **real** process RSS signal to complement the logical
//! `SpaceTracker` counters used in other parts of the codebase.

use libc::rusage;

/// Return the process peak RSS in bytes (best-effort, platform-dependent).
///
/// `ru_maxrss` units differ by platform: **bytes** on Darwin (macOS/iOS),
/// **KiB** on Linux and the BSDs (FreeBSD/NetBSD/OpenBSD/DragonFly). We treat
/// Darwin as bytes and everything else as KiB (×1024). For an unenumerated
/// target this defaults to KiB — the common case for `getrusage`, and the safe
/// (never-undercount) direction for the memory contract: under-reporting peak
/// RSS would let `--enforce`/`verify` pass a job that actually breached its
/// budget, which is exactly the failure the contract forbids.
pub fn peak_rss_bytes() -> u64 {
    let mut usage: rusage = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage as *mut rusage) };
    if rc != 0 {
        return 0;
    }
    let raw = usage.ru_maxrss as u64;

    #[cfg(any(target_os = "macos", target_os = "ios"))]
    {
        raw // already bytes
    }
    #[cfg(not(any(target_os = "macos", target_os = "ios")))]
    {
        raw.saturating_mul(1024) // KiB -> bytes (Linux, the BSDs, other unix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_rss_reflects_a_real_allocation_in_bytes() {
        // Allocate ~64 MiB and touch one byte per 4 KiB page to force residency,
        // then assert the reported peak is on the order of MiB (>= 32 MiB). This
        // catches a dropped/incorrect KiB->bytes conversion on the run platform:
        // a 1024x undercount would report ~64 KiB, far below the bound. (Lower
        // bound only — peak RSS is a process-global high-water mark, so other
        // tests can only raise it.)
        let n = 64usize * 1024 * 1024;
        let mut buf = vec![0u8; n];
        let mut i = 0;
        while i < n {
            buf[i] = 1;
            i += 4096;
        }
        std::hint::black_box(&buf);
        let peak = peak_rss_bytes();
        assert!(
            peak >= 32 * 1024 * 1024,
            "peak RSS {peak} bytes implausibly small after a 64 MiB allocation \
             (a dropped KiB->bytes conversion would land here)"
        );
    }
}
