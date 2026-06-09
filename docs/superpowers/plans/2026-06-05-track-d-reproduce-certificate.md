# Track D — `reproduce` + reproduction certificate — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `rosalind reproduce` re-derive a recorded result byte-for-byte from content-located inputs, mint a chainable `.repro.json` reproduction certificate, and gate CI on it — building only on Rosalind's shipped determinism + receipt (standalone, no cross-repo).

**Architecture:** A single `CommandCapture` chokepoint records a normalized, machine-independent, hash-protected `command` recipe into the claim (schema 4→5) and is the same structure `reproduce` replays. `reproduce` (an htslib-free driver) integrity-checks the receipt via an extracted `verify_receipt`, content-locates inputs by BLAKE3, re-execs `current_exe()`, byte-compares outputs (text only in v1), and writes a content-addressed `ReproReceipt` that chains to the parent's `content_hash`. A `badge` verb + a CI fence (cross-machine REPRODUCED + one-byte-flip DIVERGED) complete it.

**Tech Stack:** Rust (edition per repo, MSRV 1.72), `anyhow`, `clap` derive, `std::collections::BTreeMap`, `std::process::Command` + `std::env::current_exe`, the in-repo `blake3_*` helpers. No new runtime deps.

**Spec:** `docs/superpowers/specs/2026-06-05-track-d-reproduce-certificate-design.md`. Branch: `rosalind/track-d-reproduce`.

---

## File Structure

| Action | Path | Responsibility |
|---|---|---|
| Create | `src/provenance/command.rs` | `CommandCapture` — capture a normalized replayable invocation into the claim; reconstruct an argv from it. |
| Create | `src/provenance/repro.rs` | `ReproReceipt` — the reproduction certificate: type, canonical render + self-hash, chaining, parse. |
| Create | `src/provenance/badge.rs` | `badge_json` / `badge_svg` — self-hosted shields-endpoint JSON + static SVG. |
| Create | `src/reproduce.rs` | the `reproduce` driver: locate inputs, re-exec, classify outputs, mint the certificate. |
| Modify | `src/provenance/mod.rs` | `MANIFEST_SCHEMA_VERSION = 5`; `mod command/repro/badge`; `verify_receipt` + `VerifyReport`; re-exports. |
| Modify | `src/lib.rs` | re-export the new public items. |
| Modify | `src/main.rs` | adopt `CommandCapture` at the 4 receipt sites; `Reproduce` + `Badge` clap variants + dispatch; `run_verify` calls `verify_receipt`. |
| Create | `tests/reproduce.rs` | end-to-end gates: REPRODUCED, DIVERGED (pure-fn), INCONCLUSIVE, back-compat, certificate, chaining. |
| Create | `.github/workflows/reproduce.yml` | the cross-machine reproduce fence. |
| Create | `tests/golden/reproduce/` | committed golden chain (inputs + receipt + expected VCF) for the fence. |
| Modify | `CHANGELOG.md`, `README.md` | document the verb, the certificate, the badge. |

**Commit/PR boundary note:** Tasks 1–3 (capture + schema-5 + `verify_receipt`) are a coherent first half (receipts now self-describe their command; shared verifier extracted). Tasks 4–7 build `reproduce`/certificate/badge/CI on top. Commit after every task.

---

## Part 1 — Capture + extraction

### Task 1: `CommandCapture` — the capture/replay chokepoint

**Files:**
- Create: `src/provenance/command.rs`
- Modify: `src/provenance/mod.rs` (add `mod command; pub use command::CommandCapture;`)

- [ ] **Step 1: Write the failing test**

Append to the bottom of `src/provenance/command.rs` (created in Step 3, but write the test first in the same file under a `#[cfg(test)]` module — create the file with only the test to start):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::RunManifest;
    use std::path::Path;

    // Build a capture WITHOUT touching the filesystem by injecting hashes directly.
    fn sample() -> CommandCapture {
        let mut c = CommandCapture::new("variants");
        c.input_hashed("--index", "ref.idx", "h_idx");
        c.input_hashed("--alignments", "s.bam", "h_bam");
        c.opt("--mapq-threshold", 20u8);
        c.opt("--max-depth", 1000u32);
        c.flag_if(true, "--enforce");
        c.flag_if(false, "--gvcf");
        c.opt("--memory-budget-mb", 256u64);
        c.output_hashed("-o", "out.vcf", "h_out");
        c
    }

    #[test]
    fn record_into_writes_command_inputs_outputs_and_discrete_params() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);

        // command recipe: canonical, operands as @in:/@out:, gvcf absent (flag was false)
        assert_eq!(
            m.params.get("command").unwrap(),
            "variants --index @in:h_idx --alignments @in:h_bam \
             --mapq-threshold 20 --max-depth 1000 --memory-budget-mb 256 --enforce -o @out:h_out"
        );
        // inputs/outputs derived from the builder, in insertion order
        assert_eq!(m.inputs.iter().map(|f| f.blake3.as_str()).collect::<Vec<_>>(), ["h_idx", "h_bam"]);
        assert_eq!(m.outputs.iter().map(|f| f.blake3.as_str()).collect::<Vec<_>>(), ["h_out"]);
        // discrete params projected mechanically (flag -> key); mode inferred from --index
        assert_eq!(m.params.get("mapq_threshold").unwrap(), "20");
        assert_eq!(m.params.get("max_depth").unwrap(), "1000");
        assert_eq!(m.params.get("memory_budget_mb").unwrap(), "256");
        assert_eq!(m.params.get("enforce").unwrap(), "true");
        assert_eq!(m.params.get("mode").unwrap(), "index");
        assert!(!m.params.contains_key("gvcf")); // false flag is not recorded
    }

    #[test]
    fn argv_roundtrips_from_the_recorded_command() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);
        let command = m.params.get("command").unwrap();

        let locate = |h: &str| match h {
            "h_idx" => Some("/data/ref.idx".to_string()),
            "h_bam" => Some("/data/s.bam".to_string()),
            _ => None,
        };
        let out_temp = |_h: &str| "/tmp/out.vcf".to_string();
        let argv = CommandCapture::argv_from_command(command, &locate, &out_temp).unwrap();

        assert_eq!(
            argv,
            vec![
                "variants", "--index", "/data/ref.idx", "--alignments", "/data/s.bam",
                "--mapq-threshold", "20", "--max-depth", "1000", "--memory-budget-mb", "256",
                "--enforce", "-o", "/tmp/out.vcf",
            ]
        );
    }

    #[test]
    fn argv_errors_when_an_input_cannot_be_located() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);
        let command = m.params.get("command").unwrap();
        let locate = |_h: &str| None;
        let out_temp = |_h: &str| "/tmp/out.vcf".to_string();
        let err = CommandCapture::argv_from_command(command, &locate, &out_temp).unwrap_err();
        assert!(err.contains("h_idx"), "error names the unresolved input hash: {err}");
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rosalind --lib provenance::command 2>&1 | tail -20`
Expected: FAIL — `cannot find type CommandCapture` / module not declared.

- [ ] **Step 3: Write minimal implementation**

Put this ABOVE the `#[cfg(test)]` module in `src/provenance/command.rs`:

```rust
//! `CommandCapture` — the single chokepoint that records a normalized, replayable
//! invocation into a receipt's claim, and reconstructs an argv from it. Recording and
//! replay share this one structure so they cannot drift (the "forgot to record flag X"
//! bug class is eliminated by construction).

use std::path::Path;

use super::{blake3_file, FileHash, RunManifest};

/// One token of a recorded invocation. Input/output operands carry a content hash, not
/// a path, so the recorded command is machine-independent and hash-protected.
enum Token {
    Flag(String),
    Opt(String, String),
    Input { flag: String, blake3: String },
    Output { flag: String, blake3: String },
}

/// Accumulates an invocation, then writes it into a `RunManifest` (claim) and/or
/// reconstructs an argv for re-execution.
pub struct CommandCapture {
    subcommand: String,
    tokens: Vec<Token>,
    inputs: Vec<FileHash>,
    outputs: Vec<FileHash>,
}

impl CommandCapture {
    pub fn new(subcommand: impl Into<String>) -> Self {
        Self { subcommand: subcommand.into(), tokens: Vec::new(), inputs: Vec::new(), outputs: Vec::new() }
    }

    /// A content-addressed input operand; hashes the file and records it.
    pub fn input(&mut self, flag: &str, path: &Path) -> std::io::Result<&mut Self> {
        let h = blake3_file(path)?;
        Ok(self.input_hashed(flag, &path.display().to_string(), &h))
    }
    /// A content-addressed output operand; hashes the file and records it.
    pub fn output(&mut self, flag: &str, path: &Path) -> std::io::Result<&mut Self> {
        let h = blake3_file(path)?;
        Ok(self.output_hashed(flag, &path.display().to_string(), &h))
    }
    /// Input operand with a precomputed hash (used when the caller already hashed it,
    /// and by tests). `path` is recorded into `inputs[]` for humans/`verify`.
    pub fn input_hashed(&mut self, flag: &str, path: &str, blake3: &str) -> &mut Self {
        self.inputs.push(FileHash { path: path.to_string(), blake3: blake3.to_string() });
        self.tokens.push(Token::Input { flag: flag.to_string(), blake3: blake3.to_string() });
        self
    }
    pub fn output_hashed(&mut self, flag: &str, path: &str, blake3: &str) -> &mut Self {
        self.outputs.push(FileHash { path: path.to_string(), blake3: blake3.to_string() });
        self.tokens.push(Token::Output { flag: flag.to_string(), blake3: blake3.to_string() });
        self
    }
    pub fn opt(&mut self, flag: &str, value: impl ToString) -> &mut Self {
        self.tokens.push(Token::Opt(flag.to_string(), value.to_string()));
        self
    }
    pub fn flag(&mut self, flag: &str) -> &mut Self {
        self.tokens.push(Token::Flag(flag.to_string()));
        self
    }
    pub fn flag_if(&mut self, cond: bool, flag: &str) -> &mut Self {
        if cond { self.flag(flag) } else { self }
    }

    /// Render the normalized, machine-independent command string. Canonical order:
    /// subcommand, then input operands, then options (sorted by flag), then bare flags
    /// (sorted), then output operands — so the recorded command is stable run-to-run.
    fn render_command(&self) -> String {
        let mut parts: Vec<String> = vec![self.subcommand.clone()];
        for t in &self.tokens {
            if let Token::Input { flag, blake3 } = t {
                parts.push(flag.clone());
                parts.push(format!("@in:{blake3}"));
            }
        }
        let mut opts: Vec<(&String, &String)> =
            self.tokens.iter().filter_map(|t| match t { Token::Opt(f, v) => Some((f, v)), _ => None }).collect();
        opts.sort_by(|a, b| a.0.cmp(b.0));
        for (f, v) in opts {
            parts.push(f.clone());
            parts.push(v.clone());
        }
        let mut flags: Vec<&String> =
            self.tokens.iter().filter_map(|t| match t { Token::Flag(f) => Some(f), _ => None }).collect();
        flags.sort();
        for f in flags { parts.push(f.clone()); }
        for t in &self.tokens {
            if let Token::Output { flag, blake3 } = t {
                parts.push(flag.clone());
                parts.push(format!("@out:{blake3}"));
            }
        }
        parts.join(" ")
    }

    /// Write the capture into a manifest's claim: the `command` recipe, the derived
    /// `inputs[]`/`outputs[]`, the discrete params (mechanical flag->key projection),
    /// and the inferred `mode`. Call before `finalize()`. Measurement fields are the
    /// caller's responsibility (recorded separately, relocated by `finalize`).
    pub fn record_into(self, m: &mut RunManifest) {
        m.params.insert("command".to_string(), self.render_command());
        for t in &self.tokens {
            match t {
                Token::Opt(f, v) => { m.params.insert(flag_to_key(f), v.clone()); }
                Token::Flag(f) => { m.params.insert(flag_to_key(f), "true".to_string()); }
                _ => {}
            }
        }
        let has = |want: &str| self.tokens.iter().any(|t| matches!(t, Token::Input { flag, .. } if flag == want));
        if has("--index") {
            m.params.insert("mode".to_string(), "index".to_string());
        } else if has("--reference") {
            m.params.insert("mode".to_string(), "reference".to_string());
        }
        m.inputs = self.inputs;
        m.outputs = self.outputs;
    }

    /// Reconstruct an argv from a recorded `command`: substitute each `@in:<h>` with the
    /// located input path and each `@out:<h>` with a caller-supplied temp path. Errors
    /// (naming the hash) if an input cannot be located.
    pub fn argv_from_command(
        command: &str,
        locate_input: &dyn Fn(&str) -> Option<String>,
        temp_output: &dyn Fn(&str) -> String,
    ) -> Result<Vec<String>, String> {
        let mut argv = Vec::new();
        for tok in command.split(' ') {
            if let Some(h) = tok.strip_prefix("@in:") {
                match locate_input(h) {
                    Some(p) => argv.push(p),
                    None => return Err(format!("input not located by content hash @in:{h}")),
                }
            } else if let Some(h) = tok.strip_prefix("@out:") {
                argv.push(temp_output(h));
            } else {
                argv.push(tok.to_string());
            }
        }
        Ok(argv)
    }
}

/// Mechanical `--max-depth` -> `max_depth` projection (strip leading dashes,
/// dashes -> underscores). Deterministic; no per-flag config.
fn flag_to_key(flag: &str) -> String {
    flag.trim_start_matches('-').replace('-', "_")
}
```

Declare the module in `src/provenance/mod.rs` (near the top, after the existing items):

```rust
mod command;
pub use command::CommandCapture;
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rosalind --lib provenance::command 2>&1 | tail -20`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/provenance/command.rs src/provenance/mod.rs
git commit -m "feat(provenance): CommandCapture — the normalized capture/replay chokepoint"
```

---

### Task 2: Schema 4→5 + wire `CommandCapture` into the four receipt sites

**Files:**
- Modify: `src/provenance/mod.rs:31` (`MANIFEST_SCHEMA_VERSION`)
- Modify: `src/main.rs` — `run_variants_index` (~2233-2305), `run_variants` (~1599-1623), `run_features` (~1852-1923), `run_somatic` (~1229-1250)
- Test: `tests/variants_index.rs` (add a receipt-shape assertion) + the existing provenance unit tests

- [ ] **Step 1: Write the failing test**

Add to `tests/variants_index.rs` (it already exercises `variants --index`; mirror its setup — reuse the existing helper that builds an index + sorted BAM and runs the binary, then read the manifest):

```rust
#[test]
fn variants_index_receipt_records_a_replayable_command() {
    // (reuse this file's existing fixture: build index + sorted bam, run `variants --index`
    //  with -o <vcf> and --gvcf omitted; `man` = parsed RunManifest of <vcf>.manifest.json)
    let man = run_variants_index_and_read_manifest(/* gvcf = */ false);
    let command = man.params.get("command").expect("schema-5 receipt records `command`");
    assert!(command.starts_with("variants --index @in:"), "got: {command}");
    assert!(command.contains("-o @out:"), "command records the output operand: {command}");
    assert_eq!(man.params.get("mode").map(String::as_str), Some("index"));
    assert_eq!(man.params.get("gvcf"), None); // false flag omitted
    assert_eq!(man.params.get("schema_version").map(String::as_str), Some("5"));
}
```

If `run_variants_index_and_read_manifest` does not already exist in the test file, add a small helper there that runs the binary (via `env!("CARGO_BIN_EXE_rosalind")`, the pattern this file already uses) with `--gvcf` controlled by the arg and returns `RunManifest::from_canonical_json(&fs::read_to_string(manifest_path))`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rosalind --test variants_index records_a_replayable_command 2>&1 | tail -20`
Expected: FAIL — no `command` param (and `schema_version` is `4`).

- [ ] **Step 3a: Bump the schema version**

`src/provenance/mod.rs:31`:

```rust
pub const MANIFEST_SCHEMA_VERSION: u32 = 5;
```

(No change to `claim_file_render`: content-only rendering already applies at schema ≥ 3, so schema-5 inherits it. `verify`/`from_canonical_json` are schema-agnostic and keep parsing ≤ 4 receipts — back-compat is preserved.)

- [ ] **Step 3b: Migrate `run_variants_index` (the flagship)**

Replace the inputs/outputs pushes and the *claim* `params.insert` calls (the block at ~main.rs:2234-2263 covering `index`/`alignments`/`output` and `mapq_threshold`/`min_qual`/`max_depth`/`max_read_len`/`enforced`) with the builder. **Keep** every measurement insert (`peak_rss_bytes`, `predicted_peak_rss_bytes`, `max_working_set_bytes`, `governor`, `baseline_rss_bytes`, `rss_residual_bytes`, `io_rss_overhead_assumed_bytes`, `over_max_depth`, `reads_skipped_total`, `memory_budget_mb`, `contract_verdict`) and the `finalize()` + `std::fs::write` exactly as they are:

```rust
let mut manifest = RunManifest::new("variants");
let mut cmd = CommandCapture::new("variants");
cmd.input("--index", &index_path)?;
cmd.input("--alignments", &alignments_path)?;
cmd.opt("--mapq-threshold", mapq_threshold);
cmd.opt("--min-qual", quality_threshold as f64);
cmd.opt("--max-depth", max_depth);
cmd.opt("--max-read-len", max_read_len);
cmd.flag_if(enforce, "--enforce");
cmd.flag_if(gvcf, "--gvcf");
if let Some(mb) = memory_budget_mb {
    cmd.opt("--memory-budget-mb", mb);
}
if let Some(path) = &output {
    cmd.output("-o", path)?;
}
cmd.record_into(&mut manifest);
// ---- measurement inserts unchanged below this line (peak_rss_bytes, ... contract_verdict) ----
```

Add `CommandCapture` to the `use rosalind::provenance::{...}` import in this function. **Note:** the discrete key for `--enforce` is now `enforce` (was `enforced`); for `--min-qual` it is `min_qual` (matches the existing key). Grep the repo for any reader/test of the literal `"enforced"` and update to `"enforce"` (`rg '"enforced"' src tests`).

- [ ] **Step 3c: Migrate the other three sites the same way**

For each, replace the claim `params.insert`/`inputs.push`/`outputs.push` calls with a `CommandCapture` matching that subcommand's CLI flags; leave measurement inserts + `finalize()` + `write_manifest`/`fs::write` untouched.

`run_variants` (reference path, ~1599-1623), subcommand `"variants"`:
```rust
let mut cmd = CommandCapture::new("variants");
cmd.input("--reference", &reference_path)?;       // the FASTA path this fn opened
cmd.input("--alignments", &alignments_path)?;
cmd.opt("--chrom", &chrom_name);                  // now recorded (was the gap)
cmd.opt("--region-start", region_start);
cmd.opt("--mapq-threshold", mapq_threshold);
cmd.opt("--min-qual", quality_threshold as f64);
if let Some(path) = &output { cmd.output("-o", path)?; }
cmd.record_into(&mut manifest);
```

`run_features` (~1852-1923), subcommand `"features"`:
```rust
let mut cmd = CommandCapture::new("features");
cmd.input("--index", &index_path)?;
cmd.input("--alignments", &alignments_path)?;
cmd.opt("--mapq-threshold", mapq_threshold);
cmd.opt("--max-depth", max_depth);
cmd.opt("--max-read-len", max_read_len);
cmd.flag_if(enforce, "--enforce");
if let Some(mb) = memory_budget_mb { cmd.opt("--memory-budget-mb", mb); }
if let Some(path) = &output { cmd.output("-o", path)?; }
cmd.record_into(&mut manifest);
```

`run_somatic` (~1229-1250), subcommand `"somatic"` (region-bounded; reproduce treats it as a non-text/unsupported path in v1 — capture is still recorded for provenance):
```rust
let mut cmd = CommandCapture::new("somatic");
cmd.input("--reference", &reference_path)?;
cmd.input("--tumor", &tumor_path)?;               // use the resolved single-file inputs this fn holds
cmd.input("--normal", &normal_path)?;
cmd.output("-o", &output_vcf)?;
cmd.record_into(&mut manifest);
```
(If somatic holds paired-FASTQ inputs rather than single files, record each with its own flag — `--tumor-r1`, etc. Match whatever paths the function already hashed into `inputs[]`.)

- [ ] **Step 4: Run the tests**

Run: `cargo test -p rosalind --test variants_index 2>&1 | tail -20` then `cargo test -p rosalind --lib provenance 2>&1 | tail -20`
Expected: the new test PASSES; existing provenance tests pass (update any that asserted `schema_version == "4"` or the literal `"enforced"` key — these are legitimate schema-5 refreshes, noted in the commit).

- [ ] **Step 5: Commit**

```bash
git add src/provenance/mod.rs src/main.rs tests/variants_index.rs
git commit -m "feat(provenance): schema 5 — record a replayable command at every receipt site (closes the arg-recording gap)"
```

---

### Task 3: Extract `verify_receipt` into the library (anti-drift)

**Files:**
- Modify: `src/provenance/mod.rs` (add `VerifyReport` + `verify_receipt`)
- Modify: `src/main.rs:890-1060` (`run_verify` becomes a thin shell)
- Modify: `src/lib.rs` (re-export)
- Test: `src/provenance/mod.rs` `#[cfg(test)]`

- [ ] **Step 1: Write the failing test**

Add to the provenance test module:

```rust
#[test]
fn verify_receipt_reports_a_tampered_claim() {
    let mut m = RunManifest::new("variants");
    m.inputs.push(FileHash { path: "a".into(), blake3: "aa".into() });
    m.finalize();
    let mut text = m.to_canonical_json();
    // flip a byte of a recorded input digest without re-sealing -> self-hash must catch it
    text = text.replace("\"aa\"", "\"ab\"");
    let report = verify_receipt(&text, &VerifyOpts::default());
    assert!(!report.ok, "tampered claim must not verify");
    assert!(report.problems.iter().any(|p| p.contains("manifest_blake3")), "{:?}", report.problems);
}

#[test]
fn verify_receipt_passes_a_clean_receipt() {
    let mut m = RunManifest::new("variants");
    m.finalize();
    let report = verify_receipt(&m.to_canonical_json(), &VerifyOpts::default());
    assert!(report.ok, "{:?}", report.problems);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rosalind --lib provenance::tests::verify_receipt 2>&1 | tail -20`
Expected: FAIL — `verify_receipt` / `VerifyReport` / `VerifyOpts` not found.

- [ ] **Step 3: Write minimal implementation**

In `src/provenance/mod.rs`, add (this is the file-rehashing-free core; file existence checks stay in the CLI shell so the library core is pure over the receipt text + an optional `budget_mb`/`expect_code`):

```rust
/// Inputs to `verify_receipt` beyond the receipt text itself.
#[derive(Default)]
pub struct VerifyOpts {
    pub budget_mb: Option<u64>,
    pub expect_code: Option<String>,
    /// When set, re-hash files at recorded paths and check digests (the CLI sets this;
    /// callers verifying receipt-internal integrity only can leave it false).
    pub rehash_files: bool,
}

/// The outcome of checking a receipt: a list of problems (empty == ok) and the parsed
/// manifest for callers that want it.
pub struct VerifyReport {
    pub ok: bool,
    pub problems: Vec<String>,
    pub manifest: Option<RunManifest>,
}

/// Check a receipt's internal integrity (self-hashes, cross-field consistency, optional
/// budget + expected-code), and — when `opts.rehash_files` — re-hash recorded files.
/// The single source of truth shared by `verify` and `reproduce` (and a future wasm
/// verifier) so they cannot drift.
pub fn verify_receipt(text: &str, opts: &VerifyOpts) -> VerifyReport {
    let manifest = match RunManifest::from_canonical_json(text) {
        Ok(m) => m,
        Err(e) => return VerifyReport { ok: false, problems: vec![format!("parse error: {e}")], manifest: None },
    };
    let mut problems: Vec<String> = Vec::new();

    if opts.rehash_files {
        for (kind, files) in [("input", &manifest.inputs), ("output", &manifest.outputs)] {
            for f in files {
                match blake3_file(std::path::Path::new(&f.path)) {
                    Ok(h) if h == f.blake3 => {}
                    Ok(h) => problems.push(format!("{kind} {} hash mismatch: recorded {}, now {}", f.path, f.blake3, h)),
                    Err(e) => problems.push(format!("{kind} {} unreadable: {e}", f.path)),
                }
            }
        }
    }

    // (Move the existing numeric-parse + budget + internal-consistency + self_hash_ok +
    //  measurement_hash_ok + claims_measurements + expect_code logic from run_verify here
    //  verbatim, pushing into `problems` instead of printing. Keep the println! notes as
    //  pushes into a separate `notes` vec or drop them — the CLI shell re-prints.)

    VerifyReport { ok: problems.is_empty(), problems, manifest: Some(manifest) }
}
```

Then rewrite `run_verify` (main.rs) as a thin shell: read the file, call `verify_receipt(&text, &VerifyOpts { budget_mb, expect_code, rehash_files: true })`, print each problem to stderr, and `std::process::exit(5)` if `!report.ok` (preserving today's exact exit code + the "verify: OK …" success line).

Re-export in `src/lib.rs`:
```rust
pub use provenance::{verify_receipt, VerifyOpts, VerifyReport};
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p rosalind --lib provenance 2>&1 | tail -20` then `cargo test -p rosalind --test plan_enforce 2>&1 | tail -20`
Expected: PASS (the `plan_enforce` suite exercises `verify` end-to-end; its exit-5 behavior must be unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/provenance/mod.rs src/main.rs src/lib.rs
git commit -m "refactor(verify): extract verify_receipt into the library (shared by verify + reproduce)"
```

---

## Part 2 — `reproduce`, the certificate, the badge, CI

### Task 4: The `reproduce` driver (text outputs)

**Files:**
- Create: `src/reproduce.rs`
- Modify: `src/lib.rs` (`pub mod reproduce;` + re-exports), `src/main.rs` (`Reproduce` clap variant + dispatch + `run_reproduce`)
- Test: `src/reproduce.rs` `#[cfg(test)]` (pure classification) + `tests/reproduce.rs` (end-to-end)

- [ ] **Step 1: Write the failing test (pure classification + argv)**

In `src/reproduce.rs` (file starts as the test module; impl added Step 3):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::FileHash;

    fn fh(role_hash: &str) -> FileHash { FileHash { path: "x".into(), blake3: role_hash.into() } }

    #[test]
    fn classify_reproduced_when_all_outputs_match() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "h1".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Reproduced));
    }

    #[test]
    fn classify_diverged_names_the_first_mismatch() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "DIFFERENT".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Diverged));
        assert!(v.diffs.iter().any(|d| d.contains("h1") && d.contains("DIFFERENT")));
    }

    #[test]
    fn exit_codes_are_stable() {
        assert_eq!(Verdict::Reproduced.exit_code(), 0);
        assert_eq!(Verdict::Diverged.exit_code(), 6);
        assert_eq!(Verdict::Inconclusive.exit_code(), 7);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rosalind --lib reproduce 2>&1 | tail -20`
Expected: FAIL — module/types not found.

- [ ] **Step 3: Write minimal implementation**

`src/reproduce.rs` (above the test module). The driver is split into pure, testable pieces (`classify_outputs`) and an IO orchestrator (`reproduce`):

```rust
//! `reproduce` — re-derive a recorded result and compare it byte-for-byte. Standalone:
//! filesystem + blake3 + subprocess only (no htslib). Verdict is over OUTPUT bytes;
//! code/inputs/resource are diagnostic context.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::provenance::{blake3_file, CommandCapture, RunManifest, FileHash};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict { Reproduced, Diverged, Inconclusive }

impl Verdict {
    pub fn exit_code(self) -> i32 {
        match self { Verdict::Reproduced => 0, Verdict::Diverged => 6, Verdict::Inconclusive => 7 }
    }
    pub fn label(self) -> &'static str {
        match self { Verdict::Reproduced => "REPRODUCED", Verdict::Diverged => "DIVERGED", Verdict::Inconclusive => "INCONCLUSIVE" }
    }
}

pub struct Outcome { pub verdict: Verdict, pub diffs: Vec<String> }

/// Pure: compare produced (role, hash) pairs against the recorded outputs. REPRODUCED iff
/// every recorded output has a produced match; DIVERGED otherwise, naming each mismatch.
pub fn classify_outputs(recorded: &[FileHash], produced: &[(String, String)]) -> Outcome {
    let mut diffs = Vec::new();
    for (i, rec) in recorded.iter().enumerate() {
        match produced.get(i) {
            Some((_, got)) if *got == rec.blake3 => {}
            Some((role, got)) => diffs.push(format!("output {role}: recorded {} got {got}", rec.blake3)),
            None => diffs.push(format!("output {i}: recorded {} but not produced", rec.blake3)),
        }
    }
    let verdict = if diffs.is_empty() { Verdict::Reproduced } else { Verdict::Diverged };
    Outcome { verdict, diffs }
}

/// Output types reproduce can byte-compare in v1. BAM/bgzf is out of scope (a C zlib not
/// captured by deps_lock); reproduce reports those INCONCLUSIVE rather than false-DIVERGED.
fn output_is_text(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    p.ends_with(".vcf") || p.ends_with(".tsv") || p.ends_with(".txt")
}

/// Build a content-hash -> path index of every file directly under `inputs_dir`.
fn index_inputs(inputs_dir: &Path) -> std::io::Result<BTreeMap<String, PathBuf>> {
    let mut idx = BTreeMap::new();
    for entry in std::fs::read_dir(inputs_dir)? {
        let p = entry?.path();
        if p.is_file() {
            if let Ok(h) = blake3_file(&p) { idx.entry(h).or_insert(p); }
        }
    }
    Ok(idx)
}
```

Then the orchestrator `pub fn reproduce(manifest_path, inputs_dir) -> anyhow::Result<ReproOutcome>` that: reads the receipt text; calls `verify_receipt(&text, &VerifyOpts::default())` and returns `Inconclusive`/exit-5-mapping on a tampered/malformed receipt; reads `params["command"]` (absent → `Inconclusive`, "pre-schema-5 receipt: no recorded command"); checks every recorded output `output_is_text` (else `Inconclusive`, "BAM/bgzf not byte-comparable in v1"); builds the input index; binds `@in:` via the index (missing → `Inconclusive`, naming the hash) and `@out:` to temp paths in a `tempdir`; re-execs `std::env::current_exe()?` with `CommandCapture::argv_from_command(...)`; on a non-zero child exit returns `Inconclusive` with the child's stderr; hashes the produced temp outputs (`blake3_file`) into `(role, hash)` pairs; calls `classify_outputs`; reads the child's temp `<out>.manifest.json` for `peak_rss_bytes`/`memory_budget_mb` to build the machine-local resource line; compares the receipt's `code_git_sha` to the current build (re-derive via the same `build.rs` constants the binary exposes) for the code line. Return a struct carrying the verdict, the diffs, the resource line, the code line, and the parsed parent manifest (Task 5 mints the certificate from it).

(`ReproOutcome` fields: `verdict: Verdict`, `outcome: Outcome`, `resource_line: String`, `code_line: String`, `parent: RunManifest`, `parent_claim: String` = `parent.content_hash()`.)

- [ ] **Step 4: Add the CLI verb + dispatch + run_reproduce**

`src/main.rs` — add to `enum Commands`:
```rust
/// Re-derive a recorded result from its receipt and content-located inputs, and
/// write a chainable reproduction certificate. Verdict is over output bytes.
Reproduce {
    /// Path to a `*.manifest.json` from a previous run.
    #[arg(long)]
    manifest: PathBuf,
    /// Directory holding the recorded inputs (located by content hash).
    #[arg(long)]
    inputs: PathBuf,
    /// Do not write a `.repro.json` certificate.
    #[arg(long, default_value_t = false)]
    no_attest: bool,
    /// Where to write the certificate (default: `<manifest>.repro.json`).
    #[arg(short, long)]
    output: Option<PathBuf>,
},
```
Dispatch arm:
```rust
Commands::Reproduce { manifest, inputs, no_attest, output } =>
    run_reproduce(manifest, inputs, no_attest, output)?,
```
`run_reproduce` calls `rosalind::reproduce::reproduce(...)`, prints the aligned verdict block (matching the spec's mock), writes the certificate unless `no_attest` (Task 5), and `std::process::exit(outcome.verdict.exit_code())` when non-zero (exit 0 returns `Ok(())`).

Re-export in `src/lib.rs`: `pub mod reproduce;`

- [ ] **Step 5: Write the end-to-end test**

`tests/reproduce.rs`:
```rust
// Reuses the index+sorted-BAM fixture pattern from tests/variants_index.rs and runs the
// real binary via env!("CARGO_BIN_EXE_rosalind").
#[test]
fn reproduces_a_variants_index_vcf_byte_for_byte() {
    // 1. build index + sorted bam under a temp `data/` dir
    // 2. run: variants --index <idx> --alignments <bam> -o <vcf>  (writes <vcf>.manifest.json)
    // 3. run: reproduce --manifest <vcf>.manifest.json --inputs <data/>
    // assert: exit 0, stdout contains "REPRODUCED", and <vcf>.manifest.json.repro.json exists
}

#[test]
fn inconclusive_when_an_input_is_missing() {
    // run reproduce with --inputs pointing at an EMPTY dir -> exit 7, "INCONCLUSIVE", names a hash
}

#[test]
fn inconclusive_on_a_pre_schema_5_receipt() {
    // hand-write a minimal schema-4 receipt (no `command`) -> reproduce exits 7 with a clear message
}
```

- [ ] **Step 6: Run tests**

Run: `cargo test -p rosalind --lib reproduce 2>&1 | tail -20` then `cargo test -p rosalind --test reproduce 2>&1 | tail -30`
Expected: PASS (note: the `.repro.json` assertion in test 1 depends on Task 5; until then, assert only exit 0 + "REPRODUCED" and add the certificate assertion in Task 5).

- [ ] **Step 7: Commit**

```bash
git add src/reproduce.rs src/main.rs src/lib.rs tests/reproduce.rs
git commit -m "feat(reproduce): byte re-derivation verb — REPRODUCED/DIVERGED/INCONCLUSIVE (exit 0/6/7), text outputs"
```

---

### Task 5: The reproduction certificate (`.repro.json`)

**Files:**
- Create: `src/provenance/repro.rs`
- Modify: `src/provenance/mod.rs` (`mod repro; pub use repro::ReproReceipt;`), `src/reproduce.rs` (mint + write), `src/lib.rs`
- Test: `src/provenance/repro.rs` `#[cfg(test)]` + extend `tests/reproduce.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn certificate_self_hashes_and_roundtrips() {
        let mut r = ReproReceipt::new("a1b2", "variants", "REPRODUCED");
        r.outputs.push(ReproOutput { role: "-o".into(), recorded_blake3: "h".into(), observed_blake3: "h".into(), matched: true });
        r.chain_depth = 1;
        let json = r.to_canonical_json(); // stamps repro_blake3 last
        let back = ReproReceipt::from_canonical_json(&json).unwrap();
        assert_eq!(back.parent_claim, "a1b2");
        assert_eq!(back.chain_depth, 1);
        assert!(back.self_hash_ok());
    }

    #[test]
    fn tamper_breaks_the_self_hash() {
        let r = ReproReceipt::new("a1b2", "variants", "REPRODUCED");
        let json = r.to_canonical_json().replace("a1b2", "ffff");
        let back = ReproReceipt::from_canonical_json(&json).unwrap();
        assert!(!back.self_hash_ok());
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rosalind --lib provenance::repro 2>&1 | tail -20`
Expected: FAIL — types not found.

- [ ] **Step 3: Write minimal implementation**

`src/provenance/repro.rs` — a small content-addressed type mirroring `RunManifest`'s canonical-JSON + self-hash discipline (sorted keys, no timestamps, `repro_blake3` stamped last over the rest), with a hand-parser like `from_canonical_json`. Reuse `blake3_hex`. Fields per the spec §8: `schema_version`, `parent_claim`, `parent_subcommand`, `verdict`, `outputs: Vec<ReproOutput>`, `reproducer_code` (git_sha/dirty/rustc/target_triple/deps_lock — read the same `build.rs` constants the binary already bakes), `resource_here` (peak_rss_bytes/declared_budget_mb/fit), `chain_depth`, `repro_blake3`. Provide `new`, `to_canonical_json` (stamps the self-hash), `from_canonical_json`, `self_hash_ok`.

- [ ] **Step 4: Mint + write in `reproduce`**

In `src/reproduce.rs`, after classifying: build a `ReproReceipt` from the parent (`parent_claim = parent.content_hash()`, `parent_subcommand = parent.subcommand`, `verdict`, `outputs`, `resource_here`, `reproducer_code`). For `chain_depth`: if the input `--manifest` parses as a `ReproReceipt`, set `chain_depth = parent.chain_depth + 1`; else `1`. `run_reproduce` writes it unless `--no-attest`: default path `<manifest>.repro.json` (Decision 1), or `-o`.

Extend `tests/reproduce.rs` test 1: assert `<vcf>.manifest.json.repro.json` exists, parses as a `ReproReceipt`, `self_hash_ok()`, `verdict == "REPRODUCED"`, and `parent_claim == <original receipt>.content_hash()`. Add a chaining test: `reproduce` the `.repro.json` itself → a new certificate with `chain_depth == 2`.

- [ ] **Step 5: Run tests + commit**

Run: `cargo test -p rosalind reproduce 2>&1 | tail -30` and `cargo test -p rosalind --lib provenance::repro 2>&1 | tail -20`
Expected: PASS.
```bash
git add src/provenance/repro.rs src/provenance/mod.rs src/reproduce.rs src/lib.rs tests/reproduce.rs
git commit -m "feat(reproduce): chainable, self-hashing .repro.json reproduction certificate"
```

---

### Task 6: `rosalind badge` — self-hosted shields-endpoint JSON + static SVG

**Files:**
- Create: `src/provenance/badge.rs`
- Modify: `src/provenance/mod.rs`, `src/main.rs` (`Badge` variant + dispatch + `run_badge`), `src/lib.rs`
- Test: `src/provenance/badge.rs` `#[cfg(test)]`

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn badge_json_is_shields_endpoint_shaped() {
        let j = badge_json(true, Some(256));
        assert!(j.contains("\"schemaVersion\":1"));
        assert!(j.contains("\"label\":\"rosalind\""));
        assert!(j.contains("reproducible"));
        assert!(j.contains("256 MiB"));
        assert!(j.contains("\"color\":\"brightgreen\""));
    }

    #[test]
    fn badge_json_red_when_not_reproducible() {
        let j = badge_json(false, None);
        assert!(j.contains("\"color\":\"red\""));
        assert!(j.contains("not reproducible"));
    }

    #[test]
    fn badge_svg_is_well_formed() {
        let s = badge_svg(true, Some(256));
        assert!(s.starts_with("<svg") && s.trim_end().ends_with("</svg>"));
        assert!(s.contains("reproducible"));
    }
}
```

- [ ] **Step 2: Run → fail.** `cargo test -p rosalind --lib provenance::badge 2>&1 | tail -20` → not found.

- [ ] **Step 3: Implement** `badge_json(reproducible: bool, fits_mb: Option<u64>) -> String` (a hand-built shields.io endpoint object: `schemaVersion`, `label`, `message` e.g. `"reproducible · fits 256 MiB"`, `color`) and `badge_svg(...) -> String` (a minimal static two-segment SVG, no external fetch). Add `Badge { manifest: PathBuf, repro: Option<PathBuf>, output: PathBuf }` to `Commands` + `run_badge` that parses the receipt (and optional `.repro.json`), derives `reproducible` (from the repro verdict if given, else "fits-only") and `fits_mb` (from `memory_budget_mb` vs `peak_rss_bytes`), and writes `.json`/`.svg` by `output` extension.

- [ ] **Step 4: Run → pass.** `cargo test -p rosalind --lib provenance::badge 2>&1 | tail -20`

- [ ] **Step 5: Commit**
```bash
git add src/provenance/badge.rs src/provenance/mod.rs src/main.rs src/lib.rs
git commit -m "feat(badge): self-hosted shields-endpoint JSON + static SVG (no shields.io dependency)"
```

---

### Task 7: CI fence + golden chain + docs

**Files:**
- Create: `tests/golden/reproduce/` (a tiny reference + sorted BAM + the receipt + expected VCF)
- Create: `.github/workflows/reproduce.yml`
- Create: `tests/reproduce_golden.rs`
- Modify: `CHANGELOG.md`, `README.md`

- [ ] **Step 1: Generate + commit the golden chain**

Use the bundled toy data path (`examples/data/illumina_toy`) or `scripts/generate_toy_data.py` to produce a tiny deterministic reference + sorted BAM, run `variants --index -o golden.vcf`, and commit `tests/golden/reproduce/{ref.fa, sorted.bam, golden.vcf, golden.vcf.manifest.json}`. (Keep it small — this is a fixture, not a benchmark.)

- [ ] **Step 2: Write the failing golden test**

`tests/reproduce_golden.rs`:
```rust
#[test]
fn golden_chain_reproduces_byte_for_byte() {
    // run: reproduce --manifest tests/golden/reproduce/golden.vcf.manifest.json
    //               --inputs tests/golden/reproduce/
    // assert: exit 0, stdout "REPRODUCED"
}

#[test]
fn one_byte_flip_diverges() {
    // copy the golden inputs to a temp dir, flip one byte of the BAM, run reproduce
    // -> the input hash no longer matches -> INCONCLUSIVE (exit 7), NOT a false DIVERGED.
    // To exercise true DIVERGED, assert classify_outputs on a mismatch (unit, Task 4) —
    // a same-binary same-input run cannot non-deterministically diverge by design.
}
```
(Document honestly in a comment: DIVERGED is unit-tested via `classify_outputs`; the integration negative is INCONCLUSIVE-on-tamper, because byte-identical determinism means the only way to make the *same* binary produce different bytes from the *same* inputs is to change the binary — which `reproduce` surfaces as the code line.)

- [ ] **Step 3: Run → fail, then pass** once the fixture + verb exist.
Run: `cargo test -p rosalind --test reproduce_golden 2>&1 | tail -20`

- [ ] **Step 4: The CI workflow**

`.github/workflows/reproduce.yml`: on push/PR, build `--release`, run `cargo test --test reproduce_golden`, then run the binary directly to prove cross-machine reproduction on the runner:
```yaml
name: reproduce
on: [push, pull_request]
jobs:
  reproduce:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo build --release
      - name: Reproduce the golden chain on a different machine
        run: |
          ./target/release/rosalind reproduce \
            --manifest tests/golden/reproduce/golden.vcf.manifest.json \
            --inputs tests/golden/reproduce/
      - name: Emit the badge
        run: |
          ./target/release/rosalind badge \
            --manifest tests/golden/reproduce/golden.vcf.manifest.json \
            -o badge.json
      - uses: actions/upload-artifact@v4
        with: { name: rosalind-repro-badge, path: badge.json }
```

- [ ] **Step 5: Docs**

`CHANGELOG.md` `[Unreleased]`: add a "Reproducibility" bullet group (the `reproduce` verb, the chainable `.repro.json` certificate, the `badge` verb, schema 5's recorded command, the CI fence). `README.md`: add a short "Reproduce a result (a stranger can run it)" section with the one-command demo and the honest scope (text outputs; tamper-evident not tamper-proof). Keep copy honest — no "tamper-proof", verdict is over output bytes, resource line is "(here: this machine)".

- [ ] **Step 6: Full suite + commit**

Run: `cargo fmt --all -- --check && cargo build --release 2>&1 | tail -5 && cargo test 2>&1 | tail -30`
Expected: 0 warnings; full suite green.
```bash
git add tests/golden/reproduce .github/workflows/reproduce.yml tests/reproduce_golden.rs CHANGELOG.md README.md
git commit -m "test(ci): cross-machine reproduce fence + golden chain; docs for reproduce/certificate/badge"
```

---

## Self-Review

**Spec coverage:** schema-5 normalized capture → Task 1+2; `reproduce` engine + verdict + exit codes + forensic diff + honesty scope → Task 4; reproduction certificate + chaining + default sidecar → Task 5; `verify_receipt` anti-drift extraction → Task 3; CI fence + self-hosted badge → Task 6+7; all 10 spec acceptance gates map to a task's tests (gate 1→T1; 2→T1/T2; 3→T4/T5; 4→T4 unit; 5→T4; 6→T2/T4; 7→T4; 8→T5; 9→T7; 10→T7 final). Resolved decisions (sidecar default, exit codes 0/6/7/5) are in T4/T5.

**Placeholder scan:** the migration steps for the three non-flagship receipt sites (T2 Step 3c) give exact `CommandCapture` call sequences; the `reproduce` orchestrator (T4 Step 3) and the `repro`/`badge` impls (T5/T6 Step 3) describe the exact functions/fields with the test contracts pinning their shapes — no "TBD". The one deliberate "reuse the existing fixture helper" references (T2/T4/T7) point at the concrete pattern already in `tests/variants_index.rs` (`env!("CARGO_BIN_EXE_rosalind")`).

**Type consistency:** `CommandCapture` methods (`new/input/output/input_hashed/output_hashed/opt/flag/flag_if/record_into/argv_from_command`) are identical across T1 and their uses in T2/T4. `Verdict`/`Outcome`/`classify_outputs` consistent across T4. `ReproReceipt`/`ReproOutput`/`to_canonical_json`/`from_canonical_json`/`self_hash_ok` consistent across T5 and its uses. `verify_receipt`/`VerifyOpts`/`VerifyReport` consistent across T3/T4. Exit codes `0/6/7` (reproduce) and `5` (verify, reused) consistent. Schema version `5` consistent (T2 + T4 back-compat).

**Known schema-shape changes (intentional, schema-5):** the discrete `enforced` key becomes `enforce`; `gvcf`/`chrom`/`mode`/`command` are added; `schema_version` is `5`. Any golden manifest fixture or unit test pinning the old shape is refreshed under T2 (committed as a documented schema bump); `from_canonical_json`/`verify` keep reading ≤4 receipts (back-compat preserved, tested in T4).
