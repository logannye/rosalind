# Phase C3 — `rosalind verify` + receipt-on-stdout + CI contract gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the contract's trust loop: every `variants --index` run persists a self-describing, BLAKE3-stamped receipt (even on stdout), and `rosalind verify` re-checks it — proving the realized peak landed inside the budget and the outputs came from exactly these inputs — without re-running.

**Architecture:** (1) Restructure `run_variants_index` to ALWAYS write a manifest (file → `<vcf>.manifest.json`, stdout → cwd `rosalind.variants.manifest.json`, or an explicit `--manifest <path>`), with new self-describing params (`memory_budget_mb`, `contract_verdict`, `enforced`, `max_depth`, `max_read_len`). (2) Add a small hand-parser `RunManifest::from_canonical_json` for the fixed canonical shape (all values are strings), guarded by a serialize→parse→serialize round-trip property test — no `serde_json`. (3) Add `rosalind verify --manifest <path> [--budget-mb B]` that re-hashes the listed inputs/outputs and re-checks `peak_rss_bytes` ≤ the recorded/supplied budget. (4) Add a deterministic CI gate asserting the pure estimator's working-set bound ≥ the realized `max_working_set_bytes` from an actual run.

**Tech Stack:** Rust 1.72 (MSRV), `clap` derive, `cargo test`/`fmt`/`build`. No new dependencies. Builds on C1+C2 (branch `rosalind/phase-c-contract`).

**Spec:** [`docs/superpowers/specs/2026-06-01-phase-c-contract-design.md`](../specs/2026-06-01-phase-c-contract-design.md) §7.

**Gate coverage note (spec §7.5):** the five listed gates are spread across the phase — exit-3 refuse is already in C2 (`tests/plan_enforce.rs`); working-set-flat-as-input-grows is the C1 library test (`call::whole_genome::…working_set_is_bounded…not_read_count`); verify round-trip + tamper is Task 3 here; predicted-envelope ≥ realized is Task 4 here. So C3 adds the *new* deterministic gates and does not duplicate the ones already proven.

---

## File Structure

- **Modify** `src/main.rs` — `Variants` gains `--manifest`; `run_variants_index` always writes a self-describing receipt; new `Verify` subcommand + `run_verify`.
- **Modify** `src/provenance/mod.rs` — `ManifestError` + `RunManifest::from_canonical_json` + round-trip property test.
- **Modify** `tests/plan_enforce.rs` — stdout-receipt test, verify round-trip + tamper test, and the predicted-≥-realized gate.

---

## Task 1: Always write a self-describing receipt + `--manifest`

**Files:**
- Modify: `src/main.rs` (`Variants` variant, dispatch, `run_variants_index`)

- [ ] **Step 1: Add the `--manifest` flag.** In `enum Commands`, in `Variants { … }`, after the `enforce: bool,` field, add:

```rust
        /// Where to write the reproducibility receipt. Default: `<output>.manifest.json`
        /// for file output, or `./rosalind.variants.manifest.json` for stdout output.
        #[arg(long)]
        manifest: Option<PathBuf>,
```

- [ ] **Step 2: Thread it through the dispatch.** In `main()`, in the `Commands::Variants { … }` destructure, add `manifest,` after `enforce,`; and add `manifest,` as the final argument to the `run_variants_index(…)` call (after `enforce,`).

- [ ] **Step 3: Extend `run_variants_index`'s signature.** Add the parameter (after `enforce: bool,`):

```rust
    manifest_out: Option<PathBuf>,
```

- [ ] **Step 4: Replace the file-only manifest block with an always-write, self-describing one.** Replace the entire `if let Some(path) = &output { … }` manifest block (the one starting `let mut manifest = RunManifest::new("variants");`, currently ~lines 1258–1288) with:

```rust
    // Compute the contract verdict before writing the receipt (so it records it).
    let verdict = match memory_budget_mb.map(|mb| MemoryBudget::from_mb(mb).admits(peak_rss)) {
        None => "unset",
        Some(true) => "within",
        Some(false) => "over",
    };

    // Reproducibility + memory receipt — ALWAYS written (every run is verifiable):
    // an explicit --manifest path wins; else a sidecar next to the VCF; else cwd.
    let mut manifest = RunManifest::new("variants");
    manifest.inputs.push(FileHash {
        path: index_path.display().to_string(),
        blake3: blake3_file(&index_path)?,
    });
    manifest.inputs.push(FileHash {
        path: alignments_path.display().to_string(),
        blake3: blake3_file(&alignments_path)?,
    });
    if let Some(path) = &output {
        manifest.outputs.push(FileHash {
            path: path.display().to_string(),
            blake3: blake3_file(path)?,
        });
    }
    manifest
        .params
        .insert("mapq_threshold".to_string(), mapq_threshold.to_string());
    manifest.params.insert(
        "min_qual".to_string(),
        (quality_threshold as f64).to_string(),
    );
    manifest
        .params
        .insert("max_depth".to_string(), max_depth.to_string());
    manifest
        .params
        .insert("max_read_len".to_string(), max_read_len.to_string());
    manifest
        .params
        .insert("enforced".to_string(), enforce.to_string());
    manifest
        .params
        .insert("peak_rss_bytes".to_string(), peak_rss.to_string());
    manifest.params.insert(
        "max_working_set_bytes".to_string(),
        max_ws.bytes.to_string(),
    );
    if let Some(mb) = memory_budget_mb {
        manifest
            .params
            .insert("memory_budget_mb".to_string(), mb.to_string());
    }
    manifest
        .params
        .insert("contract_verdict".to_string(), verdict.to_string());

    let manifest_path: PathBuf = match (&manifest_out, &output) {
        (Some(m), _) => {
            std::fs::write(m, manifest.to_canonical_json())
                .with_context(|| format!("failed to write manifest {}", m.display()))?;
            m.clone()
        }
        (None, Some(path)) => write_manifest(path, &manifest)?,
        (None, None) => {
            let p = PathBuf::from("rosalind.variants.manifest.json");
            std::fs::write(&p, manifest.to_canonical_json())
                .with_context(|| format!("failed to write manifest {}", p.display()))?;
            p
        }
    };
    eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
```

(The `// Memory receipt:` eprintln + the `if let Some(mb) = memory_budget_mb { … }` budget/exit-4 block that follow are **unchanged** and remain after this block.)

- [ ] **Step 5: Build**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -8`
Expected: success, 0 warnings.

- [ ] **Step 6: Add a stdout-receipt test.** Append to `tests/plan_enforce.rs`:

```rust
#[test]
fn stdout_run_persists_a_self_describing_receipt() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let manifest = dir.join("run.manifest.json");
    // stdout output (no -o), explicit --manifest so we know where to look.
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");
    let json = std::fs::read_to_string(&manifest).expect("manifest written");
    for needle in [
        "\"contract_verdict\":\"within\"",
        "\"enforced\":\"true\"",
        "\"max_depth\":\"1000\"",
        "\"memory_budget_mb\":\"4096\"",
        "\"peak_rss_bytes\":",
        "\"max_working_set_bytes\":",
    ] {
        assert!(json.contains(needle), "manifest missing {needle}: {json}");
    }
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 7: Run it**

Run: `cd ~/rosalind && cargo test --test plan_enforce stdout_run_persists 2>&1 | tail -12`
Expected: PASS.

- [ ] **Step 8: Commit**

```bash
cd ~/rosalind && git add src/main.rs tests/plan_enforce.rs && git commit -m "feat(cli): variants --index always writes a self-describing receipt (+ --manifest) (C3)"
```

---

## Task 2: Canonical-manifest parser (`from_canonical_json`)

**Files:**
- Modify: `src/provenance/mod.rs`

- [ ] **Step 1: Write the failing round-trip property test.** Append to the `tests` module in `src/provenance/mod.rs`:

```rust
    #[test]
    fn parse_round_trips_canonical_json_including_escapes() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "9.9.9".to_string();
        m.inputs.push(FileHash {
            path: "weird \"path\"\twith\\escapes/和.fa".to_string(),
            blake3: "aa".to_string(),
        });
        m.inputs.push(FileHash {
            path: "a.idx".to_string(),
            blake3: "bb".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.vcf".to_string(),
            blake3: "cc".to_string(),
        });
        m.params.insert("contract_verdict".to_string(), "within".to_string());
        m.params.insert("peak_rss_bytes".to_string(), "12345".to_string());
        m.params.insert("note".to_string(), "line1\nline2".to_string());

        let json = m.to_canonical_json();
        let parsed = RunManifest::from_canonical_json(&json).expect("parse");
        // Structural equality: the parse recovers exactly what was serialized
        // (inputs are stored sorted-by-path in the canonical form, so build the
        // expected by re-parsing rather than comparing to `m`'s push order).
        assert_eq!(parsed.to_canonical_json(), json, "serialize→parse→serialize identity");
        assert_eq!(parsed.tool_version, "9.9.9");
        assert_eq!(parsed.subcommand, "variants");
        assert_eq!(parsed.params.get("note").unwrap(), "line1\nline2");
        assert_eq!(parsed.params.get("contract_verdict").unwrap(), "within");
    }

    #[test]
    fn parse_rejects_malformed() {
        assert!(RunManifest::from_canonical_json("not json").is_err());
        assert!(RunManifest::from_canonical_json("{\"inputs\":[}").is_err());
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd ~/rosalind && cargo test -p rosalind --lib provenance::tests::parse_ 2>&1 | tail -10`
Expected: FAIL — `from_canonical_json` not found.

- [ ] **Step 3: Implement the parser.** In `src/provenance/mod.rs`, add the error type (after the `use` lines, before `FileHash`):

```rust
/// Failure parsing a canonical run manifest.
#[derive(Debug)]
pub struct ManifestError(pub String);

impl std::fmt::Display for ManifestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "malformed manifest: {}", self.0)
    }
}

impl std::error::Error for ManifestError {}
```

Add the parse entry point in `impl RunManifest` (after `to_canonical_json`):

```rust
    /// Parse a manifest from its canonical JSON form (the exact shape
    /// `to_canonical_json` emits; all values are strings). A small hand-parser —
    /// no general JSON dependency. Round-trips with `to_canonical_json`.
    pub fn from_canonical_json(s: &str) -> Result<RunManifest, ManifestError> {
        let mut p = Parser { b: s.as_bytes(), i: 0 };
        let m = p.parse_manifest()?;
        Ok(m)
    }
```

Add the parser implementation (after the `RunManifest` impl block, before `push_file_hashes`):

```rust
/// Minimal recursive parser for the fixed canonical-manifest shape. Every value
/// is a JSON string (inputs/outputs are arrays of `{blake3, path}` objects;
/// params is an object of string→string), so the parser only needs strings,
/// arrays, and objects — no numbers/bools/null.
struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn err(&self, m: &str) -> ManifestError {
        ManifestError(format!("{m} at byte {}", self.i))
    }

    fn expect(&mut self, c: u8) -> Result<(), ManifestError> {
        if self.i < self.b.len() && self.b[self.i] == c {
            self.i += 1;
            Ok(())
        } else {
            Err(self.err(&format!("expected '{}'", c as char)))
        }
    }

    fn parse_string(&mut self) -> Result<String, ManifestError> {
        self.expect(b'"')?;
        let mut buf: Vec<u8> = Vec::new();
        while self.i < self.b.len() {
            let c = self.b[self.i];
            self.i += 1;
            match c {
                b'"' => {
                    return String::from_utf8(buf).map_err(|_| self.err("invalid utf-8"));
                }
                b'\\' => {
                    let e = *self.b.get(self.i).ok_or_else(|| self.err("trailing escape"))?;
                    self.i += 1;
                    match e {
                        b'"' => buf.push(b'"'),
                        b'\\' => buf.push(b'\\'),
                        b'n' => buf.push(b'\n'),
                        b'r' => buf.push(b'\r'),
                        b't' => buf.push(b'\t'),
                        b'u' => {
                            let hex = self
                                .b
                                .get(self.i..self.i + 4)
                                .ok_or_else(|| self.err("short \\u"))?;
                            let cp = u32::from_str_radix(
                                std::str::from_utf8(hex).map_err(|_| self.err("bad \\u"))?,
                                16,
                            )
                            .map_err(|_| self.err("bad \\u"))?;
                            let ch = char::from_u32(cp).ok_or_else(|| self.err("bad codepoint"))?;
                            let mut tmp = [0u8; 4];
                            buf.extend_from_slice(ch.encode_utf8(&mut tmp).as_bytes());
                            self.i += 4;
                        }
                        _ => return Err(self.err("bad escape")),
                    }
                }
                _ => buf.push(c),
            }
        }
        Err(self.err("unterminated string"))
    }

    fn expect_key(&mut self, key: &str) -> Result<(), ManifestError> {
        let k = self.parse_string()?;
        if k != key {
            return Err(self.err(&format!("expected key \"{key}\", got \"{k}\"")));
        }
        self.expect(b':')
    }

    fn parse_file_array(&mut self) -> Result<Vec<FileHash>, ManifestError> {
        self.expect(b'[')?;
        let mut out = Vec::new();
        if self.i < self.b.len() && self.b[self.i] == b']' {
            self.i += 1;
            return Ok(out);
        }
        loop {
            self.expect(b'{')?;
            self.expect_key("blake3")?;
            let blake3 = self.parse_string()?;
            self.expect(b',')?;
            self.expect_key("path")?;
            let path = self.parse_string()?;
            self.expect(b'}')?;
            out.push(FileHash { path, blake3 });
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b']') => {
                    self.i += 1;
                    break;
                }
                _ => return Err(self.err("expected ',' or ']' in array")),
            }
        }
        Ok(out)
    }

    fn parse_params(&mut self) -> Result<std::collections::BTreeMap<String, String>, ManifestError> {
        self.expect(b'{')?;
        let mut map = std::collections::BTreeMap::new();
        if self.i < self.b.len() && self.b[self.i] == b'}' {
            self.i += 1;
            return Ok(map);
        }
        loop {
            let k = self.parse_string()?;
            self.expect(b':')?;
            let v = self.parse_string()?;
            map.insert(k, v);
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b'}') => {
                    self.i += 1;
                    break;
                }
                _ => return Err(self.err("expected ',' or '}' in object")),
            }
        }
        Ok(map)
    }

    fn parse_manifest(&mut self) -> Result<RunManifest, ManifestError> {
        self.expect(b'{')?;
        self.expect_key("inputs")?;
        let inputs = self.parse_file_array()?;
        self.expect(b',')?;
        self.expect_key("outputs")?;
        let outputs = self.parse_file_array()?;
        self.expect(b',')?;
        self.expect_key("params")?;
        let params = self.parse_params()?;
        self.expect(b',')?;
        self.expect_key("subcommand")?;
        let subcommand = self.parse_string()?;
        self.expect(b',')?;
        self.expect_key("tool_version")?;
        let tool_version = self.parse_string()?;
        self.expect(b'}')?;
        Ok(RunManifest {
            tool_version,
            subcommand,
            inputs,
            params,
            outputs,
        })
    }
}
```

- [ ] **Step 4: Run the parser tests**

Run: `cd ~/rosalind && cargo test -p rosalind --lib provenance 2>&1 | tail -12`
Expected: PASS (round-trip + reject-malformed + the existing provenance tests).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/provenance/mod.rs && git commit -m "feat(provenance): RunManifest::from_canonical_json hand-parser + round-trip test (C3)"
```

---

## Task 3: `rosalind verify` subcommand

**Files:**
- Modify: `src/main.rs` (`Verify` variant, dispatch, `run_verify`)
- Modify: `tests/plan_enforce.rs` (verify tests)

- [ ] **Step 1: Add the `Verify` variant.** In `enum Commands`, after the `Plan { … }` variant (before the enum's closing `}`), add:

```rust
    /// Re-check a reproducibility receipt without re-running: re-hash its inputs
    /// and outputs and confirm the realized peak landed within the budget.
    Verify {
        /// Path to a `*.manifest.json` written by a previous run.
        #[arg(long)]
        manifest: PathBuf,
        /// Budget (MiB) to check the recorded peak against (overrides the
        /// `memory_budget_mb` recorded in the manifest, if any).
        #[arg(long)]
        budget_mb: Option<u64>,
    },
```

- [ ] **Step 2: Add the dispatch arm.** In `main()`, after the `Commands::Plan { … } => run_plan(…)?,` arm, add:

```rust
        Commands::Verify {
            manifest,
            budget_mb,
        } => run_verify(manifest, budget_mb)?,
```

- [ ] **Step 3: Implement `run_verify`.** Add this function in `src/main.rs` after `run_plan` (before `run_locate`):

```rust
/// Re-check a reproducibility receipt without re-running: parse it, re-hash each
/// listed input/output and confirm the digests match, and confirm the recorded
/// realized peak RSS landed within the budget (supplied, or recorded in the
/// manifest). Exits non-zero with a per-check report on any mismatch.
fn run_verify(manifest_path: PathBuf, budget_mb: Option<u64>) -> Result<()> {
    use rosalind::provenance::{blake3_file, RunManifest};

    let text = std::fs::read_to_string(&manifest_path)
        .with_context(|| format!("failed to read manifest {}", manifest_path.display()))?;
    let manifest = RunManifest::from_canonical_json(&text)
        .map_err(|e| anyhow!("failed to parse manifest {}: {e}", manifest_path.display()))?;

    let mut problems: Vec<String> = Vec::new();

    // Re-hash inputs + outputs against the recorded digests.
    for (kind, files) in [("input", &manifest.inputs), ("output", &manifest.outputs)] {
        for f in files {
            match blake3_file(std::path::Path::new(&f.path)) {
                Ok(h) if h == f.blake3 => {}
                Ok(h) => problems.push(format!(
                    "{kind} {} hash mismatch: recorded {}, now {}",
                    f.path, f.blake3, h
                )),
                Err(e) => problems.push(format!("{kind} {} unreadable: {e}", f.path)),
            }
        }
    }

    // Re-check the recorded realized peak against the budget (CLI overrides manifest).
    let budget_mb = budget_mb.or_else(|| {
        manifest
            .params
            .get("memory_budget_mb")
            .and_then(|v| v.parse::<u64>().ok())
    });
    match (
        budget_mb,
        manifest
            .params
            .get("peak_rss_bytes")
            .and_then(|v| v.parse::<u64>().ok()),
    ) {
        (Some(mb), Some(peak)) => {
            let budget = rosalind::core::MemoryBudget::from_mb(mb);
            if budget.admits(peak) {
                println!(
                    "verify: peak {} MiB within budget {mb} MiB",
                    peak / (1 << 20)
                );
            } else {
                problems.push(format!(
                    "recorded peak {} MiB exceeded budget {mb} MiB",
                    peak / (1 << 20)
                ));
            }
        }
        (None, _) => println!("verify: no budget to check (none supplied or recorded)"),
        (Some(_), None) => problems.push("manifest has no recorded peak_rss_bytes".to_string()),
    }

    if problems.is_empty() {
        println!("verify: OK — {} input(s), {} output(s) match", manifest.inputs.len(), manifest.outputs.len());
        Ok(())
    } else {
        for p in &problems {
            eprintln!("verify: FAIL — {p}");
        }
        std::process::exit(5);
    }
}
```

- [ ] **Step 4: Build**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -8`
Expected: success, 0 warnings.

- [ ] **Step 5: Add verify tests.** Append to `tests/plan_enforce.rs`:

```rust
#[test]
fn verify_passes_on_an_untampered_run_and_fails_on_a_tampered_output() {
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

    // Untampered → verify OK (exit 0).
    let ok = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(ok.status.success(), "verify should pass: {}", String::from_utf8_lossy(&ok.stderr));
    assert!(String::from_utf8_lossy(&ok.stdout).contains("verify: OK"));

    // Tamper with the output VCF → verify FAILS (exit 5).
    std::fs::write(&vcf, b"##tampered\n").unwrap();
    let bad = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert_eq!(bad.status.code(), Some(5), "tampered output must fail verify: {bad:?}");
    assert!(String::from_utf8_lossy(&bad.stderr).contains("hash mismatch"));
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 6: Run the verify test**

Run: `cd ~/rosalind && cargo test --test plan_enforce verify_ 2>&1 | tail -12`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
cd ~/rosalind && git add src/main.rs tests/plan_enforce.rs && git commit -m "feat(cli): rosalind verify — re-check a receipt's hashes + budget without re-running (C3)"
```

---

## Task 4: CI contract gate — pure estimate bounds the realized working set

**Files:**
- Modify: `tests/plan_enforce.rs`

- [ ] **Step 1: Write the gate.** Append to `tests/plan_enforce.rs`:

```rust
#[test]
fn estimator_upper_bounds_the_realized_working_set() {
    // Run the real pipeline, read max_working_set_bytes from the receipt, and
    // assert the pure estimator (same shared constants) is a true upper bound for
    // the declared --max-depth / --max-read-len. Deterministic (working-set
    // numbers, not process RSS).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--max-depth", "1000", "--max-read-len", "250", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    let text = std::fs::read_to_string(&manifest).unwrap();
    let m = rosalind::provenance::RunManifest::from_canonical_json(&text).unwrap();
    let realized: u64 = m.params.get("max_working_set_bytes").unwrap().parse().unwrap();

    // The fixture's single contig is 32 bp; the estimator's bound at the declared
    // cap must dominate the realized working set.
    let predicted = rosalind::call::estimate_variants_working_set(32, 1000, 250).bytes;
    assert!(
        predicted >= realized,
        "estimator bound {predicted} must be >= realized working set {realized}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Confirm the estimator is reachable at the crate root for the test.** It is re-exported in `src/call/mod.rs` as `pub use plan::{estimate_variants_working_set, …}`, so `rosalind::call::estimate_variants_working_set` resolves. If the test fails to compile on the path, fall back to `rosalind::call::plan::estimate_variants_working_set` (the module is `pub`).

- [ ] **Step 3: Run the gate**

Run: `cd ~/rosalind && cargo test --test plan_enforce estimator_upper_bounds 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
cd ~/rosalind && git add tests/plan_enforce.rs && git commit -m "test(contract): pure estimator upper-bounds the realized working set (C3 CI gate)"
```

---

## Task 5: Full verification

**Files:** none (verification only)

- [ ] **Step 1: Format**

Run: `cd ~/rosalind && cargo fmt --all && cargo fmt --all -- --check && echo FMT_CLEAN`
Expected: `FMT_CLEAN`.

- [ ] **Step 2: Zero-warning builds**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -4 && cargo build --release 2>&1 | tail -4`
Expected: both 0 warnings.

- [ ] **Step 3: Full suite**

Run: `cd ~/rosalind && cargo test 2>&1 | grep -E "test result: FAILED|panicked|[1-9][0-9]* failed|^error" | head; cargo test 2>&1 | grep -cE "test result: ok\."`
Expected: no failures; `ok.` section count ≥ C2's 26.

- [ ] **Step 4: Commit any fmt fixups**

```bash
cd ~/rosalind && git add -A && git commit -m "style: rustfmt fixups (C3)" || true
```

---

## Self-Review notes

- **Spec §7 coverage:** §7.1 receipt-on-stdout + `--manifest` → Task 1; §7.2 self-describing params → Task 1 Step 4; §7.3 parser → Task 2; §7.4 `verify` → Task 3; §7.5 CI gate → Task 4 (others mapped in the gate-coverage note up top).
- **Type consistency:** `from_canonical_json` returns `Result<RunManifest, ManifestError>` (Task 2), consumed by `run_verify` (Task 3) and the Task-4 test. `manifest_out` param added to `run_variants_index` at its definition (Task 1 Step 3) and call site (Task 1 Step 2). `RunManifest`/`FileHash`/`MemoryBudget` are existing types. `verify` exit code `5` is distinct from `--enforce`'s 3/4.
- **Determinism:** the round-trip test (Task 2) includes escaped/multi-byte content; the canonical writer already sorts inputs/outputs by path and params by key, so parse→serialize is identity. The Task-4 gate compares working-set bytes (not RSS), so it is not flaky.
