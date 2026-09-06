//! RSS / peak-memory measurement utilities.
//!
//! This module provides a **real** process RSS signal to complement the logical
//! `SpaceTracker` counters used in other parts of the codebase.

use libc::rusage;

/// Return the process peak RSS in bytes (best-effort, platform-dependent).
///
/// Linux uses `/proc/self/status`'s `VmHWM`, which measures the current executable's
/// resident high-water mark. `getrusage` can retain memory inherited from a parent
/// across `exec`, incorrectly charging a Python caller's memory to its native
/// worker. When procfs is unavailable or malformed, fall back conservatively to
/// `getrusage` instead of undercounting.
///
/// `ru_maxrss` units differ by platform: **bytes** on Darwin (macOS/iOS),
/// **KiB** on Linux and the BSDs (FreeBSD/NetBSD/OpenBSD/DragonFly). We treat
/// Darwin as bytes and everything else as KiB (×1024). For an unenumerated
/// target this defaults to KiB — the common case for `getrusage`, and the safe
/// (never-undercount) direction for the memory contract: under-reporting peak
/// RSS would let `--enforce`/`verify` pass a job that actually breached its
/// budget, which is exactly the failure the contract forbids.
pub fn peak_rss_bytes() -> u64 {
    #[cfg(target_os = "linux")]
    {
        use std::sync::atomic::{AtomicU64, Ordering};

        // Linux's batched RSS accounting can settle between samples, including
        // after unmapping an allocation. Never forget a peak we already saw.
        // This state starts at zero in a newly exec'd native worker.
        static OBSERVED_PEAK: AtomicU64 = AtomicU64::new(0);
        let bytes = linux_peak_rss_bytes().unwrap_or_else(rusage_peak_rss_bytes);
        OBSERVED_PEAK.fetch_max(bytes, Ordering::Relaxed).max(bytes)
    }
    #[cfg(not(target_os = "linux"))]
    {
        rusage_peak_rss_bytes()
    }
}

fn rusage_peak_rss_bytes() -> u64 {
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

#[cfg(target_os = "linux")]
fn linux_peak_rss_bytes() -> Option<u64> {
    use std::io::Read;

    // VmHWM normally appears near the start of this small virtual file. Limit
    // observation scratch even for unusual process metadata; a missing or
    // truncated line uses the conservative getrusage fallback above.
    let mut status = [0u8; 16 * 1024];
    let mut file = std::fs::File::open("/proc/self/status").ok()?;
    let mut used = 0;
    while used < status.len() {
        let count = file.read(&mut status[used..]).ok()?;
        if count == 0 {
            break;
        }
        used += count;
        if let Some(bytes) = parse_linux_peak_rss(&status[..used]) {
            return Some(bytes);
        }
    }
    None
}

#[cfg(target_os = "linux")]
fn parse_linux_peak_rss(status: &[u8]) -> Option<u64> {
    for line in status.split_inclusive(|byte| *byte == b'\n') {
        let Some(value) = line.strip_prefix(b"VmHWM:") else {
            continue;
        };
        if !line.ends_with(b"\n") {
            return None;
        }
        let mut fields = value
            .split(u8::is_ascii_whitespace)
            .filter(|field| !field.is_empty());
        let number = fields.next()?;
        if number.is_empty()
            || !number.iter().all(u8::is_ascii_digit)
            || fields.next()? != b"kB"
            || fields.next().is_some()
        {
            return None;
        }
        return std::str::from_utf8(number)
            .ok()?
            .parse::<u64>()
            .ok()?
            .checked_mul(1024);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_os = "linux")]
    #[test]
    fn linux_peak_parser_requires_a_complete_valid_kib_value() {
        assert_eq!(
            parse_linux_peak_rss(b"Name:\tworker\nVmHWM:\t123 kB\nVmRSS:\t100 kB\n"),
            Some(123 * 1024)
        );
        for malformed in [
            b"VmRSS:\t123 kB\n".as_slice(),
            b"VmHWM:\t123 kB",
            b"VmHWM:\t123 MB\n",
            b"VmHWM:\t-1 kB\n",
            b"VmHWM:\t123 kB extra\n",
            b"VmHWM:\t18446744073709551615 kB\n",
        ] {
            assert_eq!(parse_linux_peak_rss(malformed), None);
        }
    }

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
