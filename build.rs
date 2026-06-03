//! Bake build-identity into the binary so the reproducibility receipt records exactly
//! which code, toolchain, and dependencies produced a run. Every value degrades to
//! "unknown" (non-git build, `git` absent, …) and the script never panics, so the
//! crate always compiles and `env!("ROSALIND_*")` always resolves.

use std::path::Path;
use std::process::Command;

fn main() {
    emit_rerun_triggers();

    let git_sha = git(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
    let git_dirty = match git(&["status", "--porcelain"]) {
        Some(s) => if s.is_empty() { "false" } else { "true" }.to_string(),
        None => "unknown".to_string(),
    };
    let rustc_version = rustc_version().unwrap_or_else(|| "unknown".to_string());
    let target = std::env::var("TARGET").unwrap_or_else(|_| "unknown".to_string());
    let deps_lock_blake3 = lockfile_blake3().unwrap_or_else(|| "unknown".to_string());

    println!("cargo:rustc-env=ROSALIND_GIT_SHA={git_sha}");
    println!("cargo:rustc-env=ROSALIND_GIT_DIRTY={git_dirty}");
    println!("cargo:rustc-env=ROSALIND_RUSTC_VERSION={rustc_version}");
    println!("cargo:rustc-env=ROSALIND_TARGET={target}");
    println!("cargo:rustc-env=ROSALIND_DEPS_LOCK_BLAKE3={deps_lock_blake3}");
}

/// Re-run when the commit, staged set, or lockfile changes so the baked identity stays
/// current. Best-effort: a pure unstaged edit between builds may not re-trigger this
/// (so `code_dirty` is accurate at commit/stage granularity).
fn emit_rerun_triggers() {
    println!("cargo:rerun-if-changed=Cargo.lock");
    if Path::new(".git/HEAD").exists() {
        println!("cargo:rerun-if-changed=.git/HEAD");
        println!("cargo:rerun-if-changed=.git/index");
        if let Ok(head) = std::fs::read_to_string(".git/HEAD") {
            if let Some(r) = head.strip_prefix("ref: ") {
                println!("cargo:rerun-if-changed=.git/{}", r.trim());
            }
        }
    }
}

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn rustc_version() -> Option<String> {
    let rustc = std::env::var("RUSTC").ok()?;
    let out = Command::new(rustc).arg("--version").output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn lockfile_blake3() -> Option<String> {
    let dir = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let bytes = std::fs::read(Path::new(&dir).join("Cargo.lock")).ok()?;
    Some(blake3::hash(&bytes).to_hex().to_string())
}
