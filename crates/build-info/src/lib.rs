//! Build-script support for Rosalind-compatible producer identity.
//!
//! Add this crate as a build dependency, then call [`emit`] from `build.rs`.
//! The runtime binary can install the resulting values into
//! `rosalind::provenance::set_build_identity`.

use std::path::{Path, PathBuf};
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
    for path in git_rerun_paths(Path::new(".")) {
        println!("cargo:rerun-if-changed={}", path.display());
    }
}

/// Ask Git to resolve administrative files: linked worktrees have a `.git` file,
/// and a workspace member's working directory need not contain `.git` at all.
/// Watch loose refs even while absent, since a commit can create one from a
/// packed ref without changing HEAD or packed-refs. Cargo may rerun while an
/// optional path is absent; that is preferable to retaining stale provenance.
fn git_rerun_paths(directory: &Path) -> Vec<PathBuf> {
    let mut paths = Vec::new();
    for name in ["HEAD", "index", "packed-refs"] {
        if let Some(path) = git_in(directory, &["rev-parse", "--git-path", name]) {
            paths.push(directory.join(path));
        }
    }
    if let Some(reference) = git_in(directory, &["symbolic-ref", "--quiet", "HEAD"]) {
        if let Some(path) = git_in(directory, &["rev-parse", "--git-path", &reference]) {
            paths.push(directory.join(path));
        }
    }
    if let Some(root) = git_in(directory, &["rev-parse", "--show-toplevel"]) {
        let pointer = Path::new(&root).join(".git");
        if pointer.is_file() {
            paths.push(pointer);
        }
    }
    paths
}

fn git(args: &[&str]) -> Option<String> {
    git_in(Path::new("."), args)
}

fn git_in(directory: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .current_dir(directory)
        .args(args)
        .output()
        .ok()?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Temporary(PathBuf);
    impl Temporary {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "rosalind-build-info-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(fs::canonicalize(path).unwrap())
        }
    }
    impl Drop for Temporary {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn run(directory: &Path, command: &str, args: &[&str]) -> String {
        let output = Command::new(command)
            .current_dir(directory)
            .args(args)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{command} {args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout).unwrap().trim().into()
    }
    fn repository(path: &Path) {
        fs::create_dir(path).unwrap();
        run(path, "git", &["init", "--quiet"]);
        run(path, "git", &["symbolic-ref", "HEAD", "refs/heads/main"]);
        run(path, "git", &["config", "user.name", "Build-info fixture"]);
        run(
            path,
            "git",
            &["config", "user.email", "fixture@example.invalid"],
        );
        run(path, "git", &["config", "commit.gpgsign", "false"]);
        run(path, "git", &["config", "core.hooksPath", "/dev/null"]);
        fs::write(path.join(".gitignore"), "target/\n").unwrap();
        run(path, "git", &["add", ".gitignore"]);
        run(path, "git", &["commit", "--quiet", "-m", "initial"]);
    }
    fn assert_watches(directory: &Path, expected: &Path) {
        let expected = fs::canonicalize(expected).unwrap();
        assert!(git_rerun_paths(directory)
            .into_iter()
            .filter_map(|path| fs::canonicalize(path).ok())
            .any(|path| path == expected));
    }

    #[test]
    fn resolves_normal_nested_linked_and_packed_reference_paths() {
        let temporary = Temporary::new();
        let root = temporary.0.join("repository");
        repository(&root);
        let nested = root.join("crates/member");
        fs::create_dir_all(&nested).unwrap();
        for directory in [&root, &nested] {
            assert_watches(directory, &root.join(".git/HEAD"));
            assert_watches(directory, &root.join(".git/index"));
            assert_watches(directory, &root.join(".git/refs/heads/main"));
        }
        let linked = temporary.0.join("linked");
        run(
            &root,
            "git",
            &[
                "worktree",
                "add",
                "--quiet",
                "-b",
                "linked",
                linked.to_str().unwrap(),
            ],
        );
        let linked_nested = linked.join("crates/member");
        fs::create_dir_all(&linked_nested).unwrap();
        let admin = PathBuf::from(git_in(&linked, &["rev-parse", "--absolute-git-dir"]).unwrap());
        for directory in [&linked, &linked_nested] {
            assert_watches(directory, &linked.join(".git"));
            assert_watches(directory, &admin.join("HEAD"));
            assert_watches(directory, &admin.join("index"));
            assert_watches(directory, &root.join(".git/refs/heads/linked"));
        }
        run(&root, "git", &["pack-refs", "--all", "--prune"]);
        assert!(!root.join(".git/refs/heads/linked").exists());
        assert_watches(&linked_nested, &root.join(".git/packed-refs"));
        assert!(git_rerun_paths(&linked_nested)
            .iter()
            .any(|path| path.ends_with("refs/heads/linked")));
        run(&linked, "git", &["checkout", "--quiet", "--detach"]);
        assert_watches(&linked_nested, &admin.join("HEAD"));
        assert!(!git_rerun_paths(&linked_nested)
            .iter()
            .any(|path| path.ends_with("refs/heads/linked")));
    }

    #[test]
    fn no_repository_has_no_git_triggers_or_identity() {
        let temporary = Temporary::new();
        assert!(git_rerun_paths(&temporary.0).is_empty());
        assert_eq!(git_in(&temporary.0, &["rev-parse", "HEAD"]), None);
    }

    #[test]
    fn cargo_refreshes_embedded_identity_after_linked_worktree_commits() {
        let temporary = Temporary::new();
        let root = temporary.0.join("repository");
        repository(&root);
        fs::create_dir(root.join("src")).unwrap();
        fs::write(root.join("Cargo.toml"), format!(
            "[package]\nname = \"producer-probe\"\nversion = \"0.0.0\"\nedition = \"2021\"\n[build-dependencies]\nrosalind-build-info = {{ path = {:?} }}\n",
            env!("CARGO_MANIFEST_DIR")
        )).unwrap();
        fs::write(
            root.join("build.rs"),
            "fn main() { rosalind_build_info::emit(); }\n",
        )
        .unwrap();
        fs::write(
            root.join("src/main.rs"),
            "fn main() { println!(\"{}\", env!(\"ROSALIND_GIT_SHA\")); }\n",
        )
        .unwrap();
        let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());
        run(&root, &cargo, &["generate-lockfile", "--offline"]);
        run(&root, "git", &["add", "."]);
        run(
            &root,
            "git",
            &["commit", "--quiet", "-m", "producer fixture"],
        );
        let linked = temporary.0.join("linked");
        run(
            &root,
            "git",
            &[
                "worktree",
                "add",
                "--quiet",
                "-b",
                "linked",
                linked.to_str().unwrap(),
            ],
        );
        // A stable packed-refs file ensures the regression cannot accidentally
        // pass merely because Cargo reruns for a missing optional trigger.
        run(&root, "git", &["pack-refs", "--all", "--no-prune"]);
        for directory in [&root, &linked] {
            let target = directory.join("target");
            let build = || {
                run(
                    directory,
                    &cargo,
                    &[
                        "build",
                        "--locked",
                        "--offline",
                        "--quiet",
                        "--target-dir",
                        target.to_str().unwrap(),
                    ],
                )
            };
            let binary = target
                .join("debug")
                .join(format!("producer-probe{}", std::env::consts::EXE_SUFFIX));
            build();
            let before = run(directory, binary.to_str().unwrap(), &[]);
            assert_eq!(
                Some(before.clone()),
                git_in(directory, &["rev-parse", "HEAD"])
            );
            run(
                directory,
                "git",
                &["commit", "--quiet", "--allow-empty", "-m", "advance"],
            );
            build();
            let after = run(directory, binary.to_str().unwrap(), &[]);
            assert_ne!(before, after);
            assert_eq!(Some(after), git_in(directory, &["rev-parse", "HEAD"]));
        }
    }
}
