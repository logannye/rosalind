# `rosalind diff <a> <b>` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `rosalind diff <a> <b>`, a claim-level localizer that buckets two run receipts' differences by causal role (inputs/code-identity/science-params as causes, outputs as effect, measurements as noise) and renders a causal verdict.

**Architecture:** A pure, `std`-only `diff_receipts(a, b) -> ReceiptDiff` in the `rosalind-receipt` crate (wasm-clean; reuses `RunManifest`). A new `BUILD_IDENTITY_KEYS` const lets it segregate code-identity. A thin `rosalind diff` CLI loads two receipts, calls it, renders, and exits `0`/`1`/`2`.

**Tech Stack:** Rust 1.83, the `rosalind-receipt` crate (hand-rolled canonical JSON), clap. Tests: `#[cfg(test)]` units + a stdlib `Command` integration test.

**Spec:** `docs/superpowers/specs/2026-06-09-receipt-diff-design.md`.

**Key fact (drives the tests):** `RunManifest::finalize()` OVERWRITES the 5 build-identity claim keys (to `"unknown"` in tests) and stamps `manifest_blake3`/`schema_version`/`has_measurements`. So a test needing a specific `code_git_sha` must set it AFTER `finalize()`. `content_hash()` recomputes from current params (ignoring the stored `manifest_blake3`), so a post-finalize override is reflected in `claims_identical`.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `crates/receipt/src/lib.rs` | Add `pub const BUILD_IDENTITY_KEYS`; `mod diff; pub use diff::{…}`; a drift-guard test | Const + wiring |
| `crates/receipt/src/diff.rs` | Create | `diff_receipts` + `ReceiptDiff`/`FieldChange`/`OperandChange` + verdict/exit/json + unit tests |
| `src/main.rs` | Add `Diff` clap variant + dispatch + `run_diff` | The CLI |
| `tests/diff.rs` | Create | `rosalind diff` end-to-end |

---

## Task 1: `BUILD_IDENTITY_KEYS` const + drift guard

**Files:** Modify `crates/receipt/src/lib.rs`

- [ ] **Step 1: Write the failing drift-guard test**

In `crates/receipt/src/lib.rs`'s `#[cfg(test)] mod tests`, add:
```rust
    #[test]
    fn build_identity_keys_match_the_pairs_helper() {
        let keys: Vec<&str> = build_identity_pairs().iter().map(|(k, _)| *k).collect();
        assert_eq!(keys, BUILD_IDENTITY_KEYS, "BUILD_IDENTITY_KEYS must track build_identity_pairs");
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p rosalind-receipt build_identity_keys_match 2>&1 | tail -6`
Expected: FAIL — `cannot find value BUILD_IDENTITY_KEYS in this scope`.

- [ ] **Step 3: Add the const** (in `lib.rs`, immediately after the `MEASUREMENT_KEYS` const):
```rust
/// The claim keys that record build identity (code / toolchain / deps). Segregated by
/// `rosalind diff` as a distinct cause bucket. Must track [`build_identity_pairs`].
pub const BUILD_IDENTITY_KEYS: &[&str] = &[
    "code_git_sha",
    "code_dirty",
    "rustc_version",
    "target_triple",
    "deps_lock_blake3",
];
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p rosalind-receipt build_identity_keys_match 2>&1 | tail -6`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/receipt/src/lib.rs
git commit -m "feat(receipt): BUILD_IDENTITY_KEYS const (the code-identity claim keys)

A pub const naming the 5 build-identity claim keys, guarded by a test that it
tracks build_identity_pairs(). rosalind diff uses it to segregate code-identity.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: the pure `diff_receipts` core

**Files:** Create `crates/receipt/src/diff.rs`; Modify `crates/receipt/src/lib.rs` (wire the module)

- [ ] **Step 1: Wire the module** — in `crates/receipt/src/lib.rs`, after `pub use chain::{…};`, add:
```rust
mod diff;
pub use diff::{diff_receipts, FieldChange, OperandChange, ReceiptDiff};
```

- [ ] **Step 2: Write `crates/receipt/src/diff.rs` with types, an unimplemented core, and failing tests**

```rust
//! `rosalind diff` — a pure, claim-level divergence localizer over two run receipts.
//! It names WHICH hashed field differs, bucketed by causal role: inputs / code-identity /
//! science-params (causes), outputs (effect), measurements (machine-dependent noise).
//! No I/O, no htslib — wasm-portable. Compares two receipts to each other (distinct from
//! `reproduce`, which compares a receipt's recorded outputs against freshly-produced ones).

use std::collections::BTreeMap;

use crate::{RunManifest, BUILD_IDENTITY_KEYS};

/// Params that are derived or redundant, excluded from the science-params bucket:
/// `manifest_blake3` is the self-hash (a function of everything else); `command` is the
/// recipe whose operands/opts are already surfaced by the input/output/param buckets.
const SKIP_PARAMS: &[&str] = &["manifest_blake3", "command"];

/// One differing scalar claim/measurement field. `None` = absent on that side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldChange {
    pub key: String,
    pub a: Option<String>,
    pub b: Option<String>,
}

/// One differing input/output operand, labeled by its CLI flag (recovered from `command`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperandChange {
    pub flag: String,
    pub a: Option<String>,
    pub b: Option<String>,
}

/// The bucketed difference between two receipts' claims.
#[derive(Debug, Clone)]
pub struct ReceiptDiff {
    pub subcommand: Option<(String, String)>,
    pub inputs: Vec<OperandChange>,
    pub outputs: Vec<OperandChange>,
    pub code_identity: Vec<FieldChange>,
    pub science_params: Vec<FieldChange>,
    pub measurements: Vec<FieldChange>,
    /// `content_hash(a) == content_hash(b)` — the cross-machine claim addresses match.
    pub claims_identical: bool,
}

/// Recover `(flag, hash)` operand pairs of `marker` (`"@in:"` or `"@out:"`) from a command.
fn operands(command: &str, marker: &str) -> BTreeMap<String, String> {
    let toks: Vec<&str> = command.split(' ').collect();
    let mut out = BTreeMap::new();
    for (i, t) in toks.iter().enumerate() {
        if let Some(h) = t.strip_prefix(marker) {
            let flag = if i > 0 { toks[i - 1].to_string() } else { "?".to_string() };
            out.insert(flag, h.to_string());
        }
    }
    out
}

fn operand_changes(a: &RunManifest, b: &RunManifest, marker: &str) -> Vec<OperandChange> {
    let ma = operands(a.params.get("command").map(String::as_str).unwrap_or(""), marker);
    let mb = operands(b.params.get("command").map(String::as_str).unwrap_or(""), marker);
    let mut flags: Vec<String> = ma.keys().chain(mb.keys()).cloned().collect();
    flags.sort();
    flags.dedup();
    flags
        .into_iter()
        .filter(|f| ma.get(f) != mb.get(f))
        .map(|f| OperandChange {
            a: ma.get(&f).cloned(),
            b: mb.get(&f).cloned(),
            flag: f,
        })
        .collect()
}

fn map_changes(
    a: &BTreeMap<String, String>,
    b: &BTreeMap<String, String>,
    keep: impl Fn(&str) -> bool,
) -> Vec<FieldChange> {
    let mut keys: Vec<String> = a.keys().chain(b.keys()).cloned().collect();
    keys.sort();
    keys.dedup();
    keys.into_iter()
        .filter(|k| keep(k))
        .filter(|k| a.get(k) != b.get(k))
        .map(|k| FieldChange {
            a: a.get(&k).cloned(),
            b: b.get(&k).cloned(),
            key: k,
        })
        .collect()
}

/// Bucket the difference between two receipts' claims by causal role. Pure.
pub fn diff_receipts(a: &RunManifest, b: &RunManifest) -> ReceiptDiff {
    let subcommand = if a.subcommand != b.subcommand {
        Some((a.subcommand.clone(), b.subcommand.clone()))
    } else {
        None
    };
    ReceiptDiff {
        subcommand,
        inputs: operand_changes(a, b, "@in:"),
        outputs: operand_changes(a, b, "@out:"),
        code_identity: map_changes(&a.params, &b.params, |k| BUILD_IDENTITY_KEYS.contains(&k)),
        science_params: map_changes(&a.params, &b.params, |k| {
            !BUILD_IDENTITY_KEYS.contains(&k) && !SKIP_PARAMS.contains(&k)
        }),
        measurements: map_changes(&a.measurements, &b.measurements, |k| k != "measurement_blake3"),
        claims_identical: a.content_hash() == b.content_hash(),
    }
}

impl ReceiptDiff {
    /// `0` if the claims are identical, else `1` (read/parse errors are `2`, in the CLI).
    pub fn exit_code(&self) -> i32 {
        if self.claims_identical {
            0
        } else {
            1
        }
    }

    /// Whether an upstream CAUSE differs (vs only the output effect / measurement noise).
    pub fn has_cause(&self) -> bool {
        self.subcommand.is_some()
            || !self.inputs.is_empty()
            || !self.science_params.is_empty()
            || !self.code_identity.is_empty()
    }

    /// A one-line causal localization.
    pub fn verdict(&self) -> String {
        if self.claims_identical {
            if self.measurements.is_empty() {
                return "IDENTICAL claims".to_string();
            }
            let keys: Vec<&str> = self.measurements.iter().map(|c| c.key.as_str()).collect();
            return format!(
                "IDENTICAL claims — only machine-dependent measurements differ ({})",
                keys.join(", ")
            );
        }
        if self.has_cause() {
            let mut causes: Vec<String> = Vec::new();
            if self.subcommand.is_some() {
                causes.push("subcommand".to_string());
            }
            if !self.code_identity.is_empty() {
                causes.push(format!(
                    "code-identity ({})",
                    self.code_identity.iter().map(|c| c.key.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
            if !self.inputs.is_empty() {
                causes.push(format!(
                    "{} input(s) ({})",
                    self.inputs.len(),
                    self.inputs.iter().map(|c| c.flag.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
            if !self.science_params.is_empty() {
                causes.push(format!(
                    "{} science param(s) ({})",
                    self.science_params.len(),
                    self.science_params.iter().map(|c| c.key.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
            let effect = if self.outputs.is_empty() {
                String::new()
            } else {
                format!("; effect: {} output(s)", self.outputs.len())
            };
            return format!("claims DIFFER — cause: {}{}", causes.join(", "), effect);
        }
        if !self.outputs.is_empty() {
            return "claims DIFFER — outputs differ with identical inputs/params/code → nondeterminism or corruption".to_string();
        }
        "claims DIFFER".to_string()
    }

    /// A compact, dependency-free JSON summary for `--json`.
    pub fn to_json(&self) -> String {
        format!(
            "{{\"claims_identical\":{},\"inputs\":{},\"outputs\":{},\"code_identity\":{},\"science_params\":{},\"measurements\":{},\"exit_code\":{}}}",
            self.claims_identical,
            self.inputs.len(),
            self.outputs.len(),
            self.code_identity.len(),
            self.science_params.len(),
            self.measurements.len(),
            self.exit_code()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build + finalize a manifest from a command, extra params, and measurements.
    /// (Empty inputs[]/outputs[] — the operand diff reads the `command` param, and
    /// `command` is itself a claim param so content_hash reflects it.)
    fn mk(sub: &str, command: &str, params: &[(&str, &str)], measurements: &[(&str, &str)]) -> RunManifest {
        let mut m = RunManifest::new(sub);
        m.params.insert("command".to_string(), command.to_string());
        for (k, v) in params {
            m.params.insert(k.to_string(), v.to_string());
        }
        for (k, v) in measurements {
            m.record_measurement(*k, *v);
        }
        m.finalize();
        m
    }

    #[test]
    fn input_change_is_a_cause_with_the_output_as_effect() {
        let a = mk("variants", "variants --index @in:IDX --alignments @in:B1 -o @out:O1", &[], &[]);
        let b = mk("variants", "variants --index @in:IDX --alignments @in:B2 -o @out:O2", &[], &[]);
        let d = diff_receipts(&a, &b);
        assert!(!d.claims_identical);
        assert_eq!(d.exit_code(), 1);
        assert_eq!(d.inputs.len(), 1);
        assert_eq!(d.inputs[0].flag, "--alignments");
        assert_eq!(d.inputs[0].a.as_deref(), Some("B1"));
        assert_eq!(d.inputs[0].b.as_deref(), Some("B2"));
        assert_eq!(d.outputs.len(), 1, "the -o change is the effect");
        assert!(d.has_cause());
        assert!(d.verdict().contains("--alignments"), "{}", d.verdict());
    }

    #[test]
    fn science_param_change_is_isolated_from_code_identity() {
        let a = mk("features", "features --index @in:I -o @out:O", &[("max_depth", "1000")], &[]);
        let b = mk("features", "features --index @in:I -o @out:O", &[("max_depth", "500")], &[]);
        let d = diff_receipts(&a, &b);
        assert_eq!(d.science_params.len(), 1);
        assert_eq!(d.science_params[0].key, "max_depth");
        assert!(d.code_identity.is_empty());
        assert!(d.inputs.is_empty());
        assert!(d.outputs.is_empty(), "same @out:O");
        assert_eq!(d.exit_code(), 1);
    }

    #[test]
    fn code_identity_drift_is_its_own_bucket() {
        // Set code_git_sha AFTER finalize (finalize overwrites build-identity keys).
        let mut a = mk("variants", "variants --index @in:I -o @out:O", &[], &[]);
        a.params.insert("code_git_sha".to_string(), "aaaa111".to_string());
        let mut b = mk("variants", "variants --index @in:I -o @out:O", &[], &[]);
        b.params.insert("code_git_sha".to_string(), "bbbb222".to_string());
        let d = diff_receipts(&a, &b);
        assert_eq!(d.code_identity.len(), 1);
        assert_eq!(d.code_identity[0].key, "code_git_sha");
        assert!(d.science_params.is_empty(), "code-identity is segregated from science params");
        assert!(d.verdict().contains("code-identity"), "{}", d.verdict());
    }

    #[test]
    fn outputs_differ_with_identical_cause_flags_nondeterminism() {
        let a = mk("variants", "variants --index @in:I -o @out:O1", &[], &[]);
        let b = mk("variants", "variants --index @in:I -o @out:O2", &[], &[]);
        let d = diff_receipts(&a, &b);
        assert!(!d.claims_identical);
        assert!(d.inputs.is_empty());
        assert!(d.science_params.is_empty());
        assert!(d.code_identity.is_empty());
        assert!(!d.has_cause());
        assert_eq!(d.outputs.len(), 1);
        assert!(d.verdict().contains("nondeterminism"), "{}", d.verdict());
    }

    #[test]
    fn identical_claims_with_only_a_measurement_diff_exit_zero() {
        let a = mk("features", "features --index @in:I -o @out:O", &[], &[("peak_rss_bytes", "100")]);
        let b = mk("features", "features --index @in:I -o @out:O", &[], &[("peak_rss_bytes", "200")]);
        let d = diff_receipts(&a, &b);
        assert!(d.claims_identical, "measurements are excluded from the claim hash");
        assert_eq!(d.exit_code(), 0);
        assert_eq!(d.measurements.len(), 1);
        assert_eq!(d.measurements[0].key, "peak_rss_bytes");
        assert!(d.verdict().contains("measurements differ"), "{}", d.verdict());
        assert_eq!(
            d.to_json(),
            "{\"claims_identical\":true,\"inputs\":0,\"outputs\":0,\"code_identity\":0,\"science_params\":0,\"measurements\":1,\"exit_code\":0}"
        );
    }
}
```

- [ ] **Step 3: Run the unit tests**

Run: `cargo test -p rosalind-receipt diff:: 2>&1 | tail -12`
Expected: PASS — all five `diff::tests::*` pass. (If you staged the core as `unimplemented!()` first, restore the implementation above.)

- [ ] **Step 4: Crate builds clean (wasm-portability guard)**

Run: `cargo build -p rosalind-receipt`
Expected: PASS — `diff.rs` uses only `std::collections::BTreeMap` + the crate's own types.

- [ ] **Step 5: Commit**

```bash
git add crates/receipt/src/diff.rs crates/receipt/src/lib.rs
git commit -m "feat(receipt): pure diff_receipts claim-level localizer

std-only, wasm-portable. Buckets two receipts' claim differences by causal role
(inputs/code-identity/science-params as causes; outputs as effect; measurements
as noise), with a causal verdict and 0/1 exit code. Compares receipt-vs-receipt
(distinct from reproduce's recorded-vs-produced compare).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: the `rosalind diff` CLI

**Files:** Modify `src/main.rs` (`Diff` clap variant + dispatch + `run_diff`); Test: `tests/diff.rs` (create)

- [ ] **Step 1: Write the failing integration test** (`tests/diff.rs`):
```rust
//! `rosalind diff <a> <b>`: claim-level divergence localization between two receipts.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let d = env::temp_dir().join(format!("rosalind-diff-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().expect("spawn rosalind")
}

fn build_index_and_sorted_bam(dir: &Path) -> (PathBuf, PathBuf) {
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 8].iter().enumerate() {
        let read = &seq[start..start + 16];
        let qual: String = std::iter::repeat_n('I', 16).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, s).unwrap();
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let bam = dir.join("sorted.bam");
    assert!(run(&["index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()]).status.success());
    assert!(run(&["align", "--reference", fa.to_str().unwrap(), "--reads", fq.to_str().unwrap(), "--format", "bam", "--output", raw.to_str().unwrap()]).status.success());
    assert!(run(&["sort", "--input", raw.to_str().unwrap(), "--output", bam.to_str().unwrap()]).status.success());
    (idx, bam)
}

fn features(dir: &Path, idx: &Path, bam: &Path, out: &str, max_depth: &str) -> PathBuf {
    let tsv = dir.join(out);
    assert!(run(&[
        "features", "--index", idx.to_str().unwrap(), "--alignments", bam.to_str().unwrap(),
        "--max-depth", max_depth, "-o", tsv.to_str().unwrap(),
    ])
    .status
    .success());
    PathBuf::from(format!("{}.manifest.json", tsv.display()))
}

#[test]
fn diff_localizes_a_science_param_change_and_exits_one() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let a = features(&d, &idx, &bam, "a.tsv", "1000");
    let b = features(&d, &idx, &bam, "b.tsv", "500");

    let out = run(&["diff", a.to_str().unwrap(), b.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(out.status.code(), Some(1), "differing claims must exit 1. stdout:\n{stdout}");
    assert!(stdout.contains("max_depth"), "must localize the param change: {stdout}");
    assert!(stdout.contains("DIFFER"), "{stdout}");

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn diff_of_a_receipt_against_itself_is_identical_and_exits_zero() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let a = features(&d, &idx, &bam, "a.tsv", "1000");

    let out = run(&["diff", a.to_str().unwrap(), a.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(out.status.code(), Some(0), "identical claims must exit 0: {stdout}");
    assert!(stdout.contains("IDENTICAL"), "{stdout}");

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --test diff 2>&1 | tail -10`
Expected: FAIL — `clap` rejects the unknown `diff` subcommand (the exit-code assertion fails).

- [ ] **Step 3: Add the `Diff` clap variant** (in the `Commands` enum — e.g. after the `Verify { … }` variant):
```rust
    /// Localize how two run receipts' claims differ, bucketed by causal role: inputs /
    /// code-identity / params (causes), outputs (effect), measurements (noise). Exit
    /// 0 = identical claims, 1 = claims differ, 2 = read/parse error.
    Diff {
        /// First receipt (`*.manifest.json`).
        a: PathBuf,
        /// Second receipt (`*.manifest.json`).
        b: PathBuf,
        /// Emit a compact JSON summary instead of the human report.
        #[arg(long)]
        json: bool,
    },
```

- [ ] **Step 4: Add the dispatch arm** (after the `Commands::Verify { … } => run_verify(…)` arm):
```rust
        Commands::Diff { a, b, json } => run_diff(a, b, json)?,
```

- [ ] **Step 5: Implement `run_diff`** (e.g. just after `run_verify`):
```rust
/// Localize how two receipts' claims differ (claim-level only; it does not re-derive or
/// diff output bytes). Exit 0 = identical, 1 = differ, 2 = read/parse error.
fn run_diff(a: PathBuf, b: PathBuf, json: bool) -> Result<()> {
    use rosalind::provenance::{diff_receipts, RunManifest};

    let parse = |p: &PathBuf| -> RunManifest {
        let text = match std::fs::read_to_string(p) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("diff: cannot read {}: {e}", p.display());
                std::process::exit(2);
            }
        };
        match RunManifest::from_canonical_json(&text) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("diff: cannot parse {}: {e}", p.display());
                std::process::exit(2);
            }
        }
    };
    let ma = parse(&a);
    let mb = parse(&b);
    let report = diff_receipts(&ma, &mb);

    let short = |o: &Option<String>| -> String {
        match o {
            Some(s) if s.len() > 12 => format!("{}…", &s[..12]),
            Some(s) => s.clone(),
            None => "(absent)".to_string(),
        }
    };

    if json {
        println!("{}", report.to_json());
    } else {
        println!("diff: {}  vs  {}", ma.subcommand, mb.subcommand);
        for c in &report.code_identity {
            println!("CAUSE  — code-identity {}  {} → {}", c.key, short(&c.a), short(&c.b));
        }
        for c in &report.inputs {
            println!("CAUSE  — input  {}  {} → {}", c.flag, short(&c.a), short(&c.b));
        }
        for c in &report.science_params {
            println!("CAUSE  — param  {}  {} → {}", c.key, short(&c.a), short(&c.b));
        }
        for c in &report.outputs {
            println!("EFFECT — output {}  {} → {}", c.flag, short(&c.a), short(&c.b));
        }
        for c in &report.measurements {
            println!("noise  — measurement {}  {} → {}", c.key, short(&c.a), short(&c.b));
        }
        println!("VERDICT: {}", report.verdict());
    }
    std::process::exit(report.exit_code());
}
```

- [ ] **Step 6: Run to verify it passes**

Run: `cargo test --test diff 2>&1 | tail -8`
Expected: PASS — `diff` exits 1 localizing `max_depth`; a receipt vs itself exits 0 (`IDENTICAL`).

- [ ] **Step 7: Commit**

```bash
git add src/main.rs tests/diff.rs
git commit -m "feat(cli): rosalind diff <a> <b> — claim-level divergence localizer

Loads two receipts, buckets their claim differences by causal role, prints a
one-line localization (cause vs effect vs noise), exits 0=identical/1=differ/
2=error. --json for a compact summary.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: full gate

**Files:** none (verification only)

- [ ] **Step 1: Full suite + lint + format**

Run:
```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```
Expected: PASS on all three. (If `cargo fmt --check` reports diffs, run `cargo fmt` and amend the relevant commit.)

- [ ] **Step 2: End-to-end smoke**

Run:
```bash
BIN=target/debug/rosalind; D=$(mktemp -d)
printf '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$D/ref.fa"
printf '@r0\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n' > "$D/reads.fq"
$BIN index --reference "$D/ref.fa" --output "$D/i.idx" >/dev/null 2>&1
$BIN align --reference "$D/ref.fa" --reads "$D/reads.fq" --format bam --output "$D/raw.bam" >/dev/null 2>&1
$BIN sort --input "$D/raw.bam" --output "$D/s.bam" >/dev/null 2>&1
$BIN features --index "$D/i.idx" --alignments "$D/s.bam" --max-depth 1000 -o "$D/a.tsv" >/dev/null 2>&1
$BIN features --index "$D/i.idx" --alignments "$D/s.bam" --max-depth 500  -o "$D/b.tsv" >/dev/null 2>&1
echo "--- diff (expect localize max_depth, exit 1) ---"
$BIN diff "$D/a.tsv.manifest.json" "$D/b.tsv.manifest.json"; echo "exit=$?"
echo "--- diff self (expect IDENTICAL, exit 0) ---"
$BIN diff "$D/a.tsv.manifest.json" "$D/a.tsv.manifest.json"; echo "exit=$?"
rm -rf "$D"
```
Expected: the first prints a `param max_depth 1000 → 500` line + `VERDICT: claims DIFFER` and `exit=1`; the second prints `VERDICT: IDENTICAL claims` and `exit=0`.

- [ ] **Step 3: No commit** (verification only).

---

## Self-Review

**1. Spec coverage:**
- §3 `diff_receipts` + `ReceiptDiff`/`FieldChange`/`OperandChange`; bucketing (operand-flag inputs/outputs, code-identity via `BUILD_IDENTITY_KEYS`, science-params excluding `manifest_blake3`/`command`, measurements excluding `measurement_blake3`) → **Task 2**. `BUILD_IDENTITY_KEYS` + drift guard → **Task 1**.
- §4 verdict + exit codes (`0`/`1`, `2` in the CLI) → **Task 2** (`verdict`/`exit_code`) + **Task 3** (`run_diff` exits `2` on read/parse error).
- §5 CLI (`Diff` variant + `run_diff` + `--json`) → **Task 3**.
- §6 tests: unit input/param/code-identity/nondeterminism/identical+measurements (Task 2); integration localize-a-param + diff-self-identical (Task 3); gates (Task 4).
- §7 non-goals (claim-level only, no re-derivation, no `classify_outputs` overlap, pure-std) → respected (no VCF/BAM bytes touched; `diff.rs` is `std`-only).

**2. Placeholder scan:** No TBD/TODO. Every step has literal code or an exact command + expected output. (Task 3 Step 5 flags the exact import line to use.)

**3. Type consistency:** `diff_receipts(&RunManifest, &RunManifest) -> ReceiptDiff`; `ReceiptDiff` fields (`inputs`/`outputs`/`code_identity`/`science_params`/`measurements`/`claims_identical`) + methods (`exit_code`/`has_cause`/`verdict`/`to_json`) are defined in Task 2 and used identically in Task 2's tests and Task 3's `run_diff`. `FieldChange { key, a, b }` / `OperandChange { flag, a, b }` field names match across the core, the renderer, and the tests. `BUILD_IDENTITY_KEYS` (Task 1) is consumed in `diff.rs` (Task 2). The CLI imports are `rosalind::provenance::{diff_receipts, RunManifest}`.
```
