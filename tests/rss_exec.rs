//! Native worker accounting must exclude allocations retained by its caller.

#![cfg(target_os = "linux")]

use std::process::Command;

use rosalind::util::rss::peak_rss_bytes;

const PARENT_BYTES: usize = 128 << 20;
const CHILD_BYTES: usize = 32 << 20;

fn resident_allocation(bytes: usize) -> Vec<u8> {
    let mut allocation = vec![0; bytes];
    for offset in (0..bytes).step_by(4096) {
        allocation[offset] = 1;
    }
    allocation
}

#[test]
fn native_worker_peak_excludes_parent_memory_and_retains_its_own_peak() {
    const CHILD: &str = "ROSALIND_TEST_RSS_EXEC_CHILD";
    if std::env::var_os(CHILD).is_none() {
        let parent_memory = resident_allocation(PARENT_BYTES);
        let output = Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "native_worker_peak_excludes_parent_memory_and_retains_its_own_peak",
                "--nocapture",
            ])
            .env(CHILD, "1")
            .output()
            .unwrap();
        std::hint::black_box(&parent_memory);
        assert!(
            output.status.success(),
            "native worker RSS regression failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    assert_eq!(unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) }, 0);
    let inherited_peak = (usage.ru_maxrss as u64) * 1024;
    assert!(inherited_peak >= PARENT_BYTES as u64);
    let initial = peak_rss_bytes();
    assert!(
        initial + (PARENT_BYTES as u64 / 2) < inherited_peak,
        "worker peak {initial} must exclude inherited parent peak {inherited_peak}"
    );

    let child_memory = resident_allocation(CHILD_BYTES);
    std::hint::black_box(&child_memory);
    let touched = peak_rss_bytes();
    assert!(touched >= initial + (CHILD_BYTES as u64 / 2));
    assert!(touched >= CHILD_BYTES as u64);
    drop(child_memory);
    assert!(
        peak_rss_bytes() >= touched,
        "released native allocations must remain in the worker high-water mark"
    );
}
