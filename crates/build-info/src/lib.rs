//! Build-script support for Rosalind-compatible producer identity.
//!
//! Add this crate as a build dependency, then call [`emit`] from `build.rs`.
//! The runtime binary can install the resulting values into
//! `rosalind::provenance::set_build_identity`.

use std::path::Path;
use std::process::Command;

/// Emit `ROSALIND_*` compile-time environment values for the package currently
/// being built. Every value degrades to `unknown`; this function never panics.
pub fn emit() {
    emit_rerun_triggers();
    let git_sha = git(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
    let git_dirty = match git(&["status", "--porcelain"]) {
        Some(status) if status.is_empty() => "false".to_string(),
        Some(_) => "true".to_string(),
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

fn emit_rerun_triggers() {
    println!("cargo:rerun-if-changed=Cargo.lock");
    if Path::new(".git/HEAD").exists() {
        println!("cargo:rerun-if-changed=.git/HEAD");
        println!("cargo:rerun-if-changed=.git/index");
        if let Ok(head) = std::fs::read_to_string(".git/HEAD") {
            if let Some(reference) = head.strip_prefix("ref: ") {
                println!("cargo:rerun-if-changed=.git/{}", reference.trim());
            }
        }
    }
}

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn rustc_version() -> Option<String> {
    let rustc = std::env::var("RUSTC").ok()?;
    let output = Command::new(rustc).arg("--version").output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn lockfile_blake3() -> Option<String> {
    let directory = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let bytes = std::fs::read(Path::new(&directory).join("Cargo.lock")).ok()?;
    Some(blake3::hash(&bytes).to_hex().to_string())
}
