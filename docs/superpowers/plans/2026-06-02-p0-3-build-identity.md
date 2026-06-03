# P0.3 — Build-identity in the receipt — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans. Steps use `- [ ]`.

**Goal:** Bake `code_git_sha` + `code_dirty` + `rustc_version` + `target_triple` + `deps_lock_blake3`
into the binary at build time, stamp them into the receipt's claim, and add `verify --expect-code`.

**Architecture:** A `build.rs` emits the five facts as `ROSALIND_*` rustc-env vars (graceful
`"unknown"` fallbacks). `finalize()` stamps them into `params` (claim, before the self-hash). A pure
`check_expected_code` powers `verify --expect-code`. Schema 3 → 4.

**Reference:** spec `docs/superpowers/specs/2026-06-02-p0-3-build-identity-design.md`.

---

### Task 1: `build.rs` + the blake3 build-dependency

**Files:** Create `build.rs`; modify `Cargo.toml`.

- [ ] **Step 1: Add the build-dependency** — in `Cargo.toml`, after the `[dev-dependencies]` block:

```toml
[build-dependencies]
# Same hash as the runtime receipt digests — for deps_lock_blake3 in build.rs.
blake3 = "1.5"
```

- [ ] **Step 2: Write `build.rs`**

```rust
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
```

- [ ] **Step 3: Verify it builds + emits** — `cargo build 2>&1 | tail -3` (clean). Sanity:
  `./target/debug/rosalind --help` still runs.

---

### Task 2: Stamp build-identity in `finalize()` (TDD)

**Files:** Modify `src/provenance/mod.rs`.

- [ ] **Step 1: Write failing unit tests** (append to `mod tests`):

```rust
    #[test]
    fn finalize_stamps_build_identity_into_the_claim() {
        let mut m = RunManifest::new("variants");
        m.finalize();
        for k in [
            "code_git_sha",
            "code_dirty",
            "rustc_version",
            "target_triple",
            "deps_lock_blake3",
        ] {
            assert!(
                m.params.get(k).is_some_and(|v| !v.is_empty()),
                "finalize must stamp {k}"
            );
        }
        // Build-identity is in the claim → tampering it breaks the claim self-hash.
        assert_eq!(m.self_hash_ok(), Some(true));
        m.params
            .insert("code_git_sha".to_string(), "tampered".to_string());
        assert_eq!(m.self_hash_ok(), Some(false));
    }

    #[test]
    fn check_expected_code_matches_mismatches_and_flags_dirty() {
        let mk = |sha: &str, dirty: &str| {
            let mut m = RunManifest::new("variants");
            m.params.insert("code_git_sha".to_string(), sha.to_string());
            m.params.insert("code_dirty".to_string(), dirty.to_string());
            m
        };
        // Exact + prefix match on a clean build → no problems.
        assert!(mk("abc123def456", "false")
            .check_expected_code("abc123def456")
            .is_empty());
        assert!(mk("abc123def456", "false")
            .check_expected_code("abc123")
            .is_empty());
        // Mismatch → one problem mentioning "mismatch".
        let mm = mk("abc123", "false").check_expected_code("deadbeef");
        assert_eq!(mm.len(), 1);
        assert!(mm[0].contains("mismatch"));
        // Match but dirty → one problem mentioning "dirty".
        let dirty = mk("abc123", "true").check_expected_code("abc123");
        assert_eq!(dirty.len(), 1);
        assert!(dirty[0].contains("dirty"));
        // Absent / unknown → cannot check (one problem each).
        assert_eq!(
            RunManifest::new("variants")
                .check_expected_code("abc")
                .len(),
            1
        );
        assert_eq!(mk("unknown", "false").check_expected_code("abc").len(), 1);
    }
```

- [ ] **Step 2: Confirm RED** — `cargo test -p rosalind --lib provenance 2>&1 | tail -12`
  (compile error: `check_expected_code` missing; `finalize_stamps_build_identity` fails).

- [ ] **Step 3: Bump the schema version**

```rust
/// … v3: content-only claim (paths dropped). v4: the claim records build-identity
/// (`code_git_sha`/`code_dirty`/`rustc_version`/`target_triple`/`deps_lock_blake3`),
/// so it commits to exactly which code, toolchain, and deps produced the run.
pub const MANIFEST_SCHEMA_VERSION: u32 = 4;
```

- [ ] **Step 4: Stamp build-identity in `finalize`** — in `finalize()`, immediately before the
  `schema_version` insert (step 3 of the existing body), add:

```rust
        // Build-identity (baked at compile time by build.rs) — part of the claim, so it
        // is committed to by the self-hash and forms the reproduction key.
        for (k, v) in [
            ("code_git_sha", env!("ROSALIND_GIT_SHA")),
            ("code_dirty", env!("ROSALIND_GIT_DIRTY")),
            ("rustc_version", env!("ROSALIND_RUSTC_VERSION")),
            ("target_triple", env!("ROSALIND_TARGET")),
            ("deps_lock_blake3", env!("ROSALIND_DEPS_LOCK_BLAKE3")),
        ] {
            self.params.insert(k.to_string(), v.to_string());
        }
```

- [ ] **Step 5: Add `check_expected_code`** (near `claims_measurements`):

```rust
    /// Check the recorded `code_git_sha` against an expected commit (prefix match, so
    /// short SHAs work). Returns the problems found: a mismatch, a clean match from a
    /// DIRTY tree (not reproducible from a SHA alone), or an inability to check
    /// (no/`unknown` SHA). An empty vec means a clean, matching build.
    pub fn check_expected_code(&self, expected: &str) -> Vec<String> {
        match self.params.get("code_git_sha").map(String::as_str) {
            None => vec![
                "cannot check --expect-code: the receipt records no code_git_sha (a pre-P0.3 receipt)"
                    .to_string(),
            ],
            Some("unknown") => vec![
                "cannot check --expect-code: the receipt's code_git_sha is 'unknown' (a non-git build)"
                    .to_string(),
            ],
            Some(sha) if !sha.starts_with(expected) => vec![format!(
                "code mismatch: receipt was built from {sha}, expected {expected}"
            )],
            Some(_) => {
                if self.params.get("code_dirty").map(String::as_str) == Some("true") {
                    vec![format!(
                        "code matches {expected} but the receipt was built from a DIRTY tree \
                         (uncommitted changes) — not reproducible from a commit SHA alone"
                    )]
                } else {
                    Vec::new()
                }
            }
        }
    }
```

- [ ] **Step 6: GREEN** — `cargo test -p rosalind --lib provenance 2>&1 | tail -12` (all pass,
  including every P0.2/P0.2b test — build-identity is identical for both manifests in the
  cross-machine tests, so they still match).

- [ ] **Step 7: `cargo fmt` + commit**

```bash
cargo fmt
git add build.rs Cargo.toml Cargo.lock src/provenance/mod.rs \
  docs/superpowers/specs/2026-06-02-p0-3-build-identity-design.md \
  docs/superpowers/plans/2026-06-02-p0-3-build-identity.md
git commit -m "feat(provenance): build-identity in the receipt — code/toolchain/deps reproduction key"
```

---

### Task 3: `verify --expect-code`

**Files:** Modify `src/main.rs`.

- [ ] **Step 1: Add the CLI arg** — in the `Verify` variant (after `budget_mb`):

```rust
        /// Assert the receipt was built from exactly this commit SHA (prefix ok).
        /// Fails verify on a mismatch, or on a clean match from a dirty build.
        #[arg(long)]
        expect_code: Option<String>,
```

- [ ] **Step 2: Thread it to `run_verify`** — update the match arm and the signature:

```rust
        Commands::Verify {
            manifest,
            budget_mb,
            expect_code,
        } => run_verify(manifest, budget_mb, expect_code)?,
```

and `fn run_verify(manifest_path: PathBuf, budget_mb: Option<u64>, expect_code: Option<String>) -> Result<()>`.

- [ ] **Step 3: Apply the check** — in `run_verify`, after the measurement-hash match block (before
  the `if problems.is_empty()`):

```rust
    // Build-identity: assert the receipt came from exactly the expected commit.
    if let Some(expected) = &expect_code {
        let code_problems = manifest.check_expected_code(expected);
        if code_problems.is_empty() {
            println!("verify: code matches {expected} (clean build)");
        } else {
            problems.extend(code_problems);
        }
    }
```

- [ ] **Step 4: Build + the verify-related integration tests compile** —
  `cargo test -p rosalind --test plan_enforce verify 2>&1 | tail -10`.

- [ ] **Step 5: `cargo fmt` + commit**

```bash
cargo fmt
git add src/main.rs
git commit -m "feat(verify): --expect-code <sha> gate against the receipt's build-identity"
```

---

### Task 4: Integration tests + full verification

**Files:** Modify `tests/plan_enforce.rs`.

- [ ] **Step 1: Bump the schema assertion** — `"\"schema_version\":\"3\""` → `"\"schema_version\":\"4\""`.

- [ ] **Step 2: Add a build-identity + expect-code e2e** (append to `tests/plan_enforce.rs`):

```rust
#[test]
fn receipt_records_build_identity_and_verify_expect_code_catches_a_mismatch() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    let text = std::fs::read_to_string(&manifest).unwrap();
    for needle in [
        "\"code_git_sha\":",
        "\"code_dirty\":",
        "\"rustc_version\":",
        "\"target_triple\":",
        "\"deps_lock_blake3\":",
    ] {
        assert!(text.contains(needle), "receipt missing {needle}: {text}");
    }

    // A wrong commit SHA must fail verify (robust regardless of the build's dirty state).
    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .args(["--expect-code", "0000000000000000000000000000000000000000"])
        .output()
        .unwrap();
    assert_eq!(
        v.status.code(),
        Some(5),
        "a wrong --expect-code must fail verify: {v:?}"
    );
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        stderr.contains("code mismatch"),
        "expected a code-mismatch error: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 3: Full suite + clippy + MSRV**

```bash
cargo test 2>&1 | grep -E "test result: FAILED|FAILED|^error"   # expect none
cargo clippy --all-targets -- -D warnings 2>&1 | tail -3        # clean
cargo +1.83 build 2>&1 | tail -3                                # MSRV (CI covers if absent)
```

- [ ] **Step 4: `cargo fmt` + commit, push, PR, watch CI, merge on green**

```bash
cargo fmt
git add tests/plan_enforce.rs
git commit -m "test: receipt schema 4 + build-identity / --expect-code e2e"
git push -u origin rosalind/p0-3-build-identity
gh pr create --title "feat(provenance): build-identity in the receipt + verify --expect-code (P0.3)" --body "<summary>"
gh pr checks --watch
```
