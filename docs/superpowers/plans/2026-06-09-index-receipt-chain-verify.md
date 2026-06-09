# Index Receipt + Chain Verify Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `rosalind index` a content-addressed receipt so the content-hash edges downstream receipts already record resolve into a real, offline-walkable provenance DAG, then add `rosalind chain verify <dir>` to walk and verify it.

**Architecture:** PR1 makes `run_index` mirror the proven `variants`/`somatic` receipt path (`RunManifest` + `CommandCapture` + `write_manifest`), recording the reference FASTA file as the chain root and the `.idx` file as the edge — so an `index` receipt's `outputs[0].blake3` is bit-identical to a `variants` receipt's `inputs[--index].blake3` by construction. PR2 adds a pure, `std`-only `walk_chain` in the `rosalind-receipt` leaf crate (wasm-portable for a future `chain confirm` board) that classifies every edge (resolved / external / broken) and checks every node's self-hash, plus a `chain verify` CLI noun that drives it offline.

**Tech Stack:** Rust 1.83 (edition 2021), `clap` derive CLI, `blake3`, the in-repo `rosalind-receipt` crate (hand-rolled canonical JSON, no serde, no htslib). Tests: stdlib `std::process::Command` against `env!("CARGO_BIN_EXE_rosalind")` + `#[cfg(test)]` unit tests in the receipt crate.

**Spec:** `docs/superpowers/specs/2026-06-09-index-receipt-chain-verify-design.md`.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/main.rs` | Modify `run_index` (≈ line 690–718); add a `Chain` command + `ChainAction` enum (near the `Index` variant, ≈ line 236) + a dispatch arm (≈ line 540) + a new `run_chain_verify` fn | Emit the index receipt; the `chain verify` CLI surface + receipt-dir loading + rendering |
| `crates/receipt/src/chain.rs` | Create | Pure DAG walk: `EdgeStatus`, `ChainEdge`, `ChainNode`, `ChainReport`, `walk_chain`, `ChainReport::to_json`; unit tests |
| `crates/receipt/src/lib.rs` | Modify (add `pub mod chain;` + a `pub use`) | Export the chain API as `rosalind::provenance::*` |
| `tests/chain.rs` | Create | Integration gates: index receipt self-verifies; the edge resolves by construction; `chain verify` INTACT / BROKEN / external-integrity-only |

**Key API already in the codebase (do not re-derive):**
- `rosalind::provenance::{RunManifest, CommandCapture, write_manifest, blake3_hex, FileHash, walk_chain (new)}`.
- `CommandCapture::input(flag, &Path)` / `output(flag, &Path)` hash the file via `blake3_file` and record a `{path, blake3}` entry + an `@in:`/`@out:` token in the `command` string.
- `RunManifest`: public fields `inputs: Vec<FileHash>`, `outputs: Vec<FileHash>`, `params: BTreeMap<String,String>`, `measurements: BTreeMap<String,String>`; methods `new`, `record_measurement`, `finalize`, `content_hash`, `self_hash_ok() -> Option<bool>`, `from_canonical_json`, `to_canonical_json`.
- `FileHash { pub path: String, pub blake3: String }`.
- `write_manifest(output_path, &manifest)` writes `<output_path>.manifest.json` and returns its path.
- The `somatic` receipt block at `src/main.rs:1221`/`1264` is the exact template to mirror.

---

## PR1 — `index` emits a content-addressed receipt

### Task 1: `index` writes a self-verifying receipt sidecar

**Files:**
- Create: `tests/chain.rs`
- Modify: `src/main.rs` (`run_index`, after the build telemetry at ≈ line 717, before `Ok(())`)

- [ ] **Step 1: Write the failing integration test (creates `tests/chain.rs` with the shared harness)**

```rust
//! Provenance-DAG gates: the `index` receipt (PR1) and `chain verify` (PR2).
//! Driven through the real CLI binary (no rust-htslib dev-dependency), reusing the
//! `index -> align -> sort -> variants` toy pipeline from `tests/reproduce.rs`.

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
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = env::temp_dir().join(format!("rosalind-chain-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

fn write_fasta(dir: &Path, name: &str, seq: &str) -> PathBuf {
    let p = dir.join("ref.fa");
    std::fs::write(&p, format!(">{name}\n{seq}\n")).unwrap();
    p
}

fn write_fastq(dir: &Path, seq: &str, starts: &[usize], len: usize) -> PathBuf {
    let p = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in starts.iter().enumerate() {
        let read = &seq[start..start + len];
        let qual: String = std::iter::repeat_n('I', len).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&p, s).unwrap();
    p
}

/// `index` -> `align --format bam` -> `sort` -> `(idx, sorted.bam)`, all in `dir`.
fn build_index_and_sorted_bam(dir: &Path, fa: &Path, fq: &Path) -> (PathBuf, PathBuf) {
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
    assert!(run(&[
        "index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()
    ])
    .status
    .success());
    assert!(run(&[
        "align", "--reference", fa.to_str().unwrap(), "--reads", fq.to_str().unwrap(),
        "--format", "bam", "--output", raw.to_str().unwrap(),
    ])
    .status
    .success());
    assert!(run(&[
        "sort", "--input", raw.to_str().unwrap(), "--output", sorted.to_str().unwrap()
    ])
    .status
    .success());
    (idx, sorted)
}

#[test]
fn index_writes_a_self_verifying_receipt() {
    let d = tmpdir();
    let fa = write_fasta(&d, "chr1", "ACGTACGTACGTACGTACGTACGTACGTACGT");
    let idx = d.join("ref.idx");

    let out = run(&[
        "index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()
    ]);
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let manifest = d.join("ref.idx.manifest.json");
    assert!(
        manifest.exists(),
        "index must write a receipt sidecar at {}",
        manifest.display()
    );

    let v = run(&["verify", "--manifest", manifest.to_str().unwrap()]);
    assert!(
        v.status.success(),
        "verify must pass on the index receipt. stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&v.stdout),
        String::from_utf8_lossy(&v.stderr)
    );

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --test chain index_writes_a_self_verifying_receipt`
Expected: FAIL — the sidecar `ref.idx.manifest.json` does not exist (the assertion `index must write a receipt sidecar` panics), because `run_index` writes no manifest today.

- [ ] **Step 3: Implement the index receipt in `run_index`**

In `src/main.rs`, in `run_index`, replace the final `Ok(())` (≈ line 718) by inserting the receipt block immediately before it. The values `output`, `reference`, `total_bp`, `reference_blake3`, and `peak` are already in scope.

```rust
    // Content-addressed receipt: makes the index a chainable, verifiable node — the
    // root of every downstream provenance chain. Mirrors the `variants`/`somatic` path.
    // The `--output` operand records blake3(.idx file) into outputs[]; that digest is
    // bit-identical to what `variants` records as its `--index` input, so the chain edge
    // resolves by construction. `--reference` records blake3(FASTA file) as the root.
    {
        use rosalind::provenance::{blake3_hex, write_manifest, CommandCapture, RunManifest};

        let mut manifest = RunManifest::new("index");
        let mut cmd = CommandCapture::new("index");
        cmd.input("--reference", &reference)?;
        if let Some(mb) = memory_budget_mb {
            cmd.opt("--memory-budget-mb", mb);
        }
        cmd.output("--output", &output)?;
        cmd.record_into(&mut manifest);
        // `reference_blake3` is the in-memory NORMALIZED sequence hash (a "what genome"
        // id, stable across FASTA reformatting) — informational, NOT the chain edge.
        manifest
            .params
            .insert("reference_blake3".to_string(), blake3_hex(&reference_blake3));
        manifest
            .params
            .insert("total_bp".to_string(), total_bp.to_string());
        // Realized build peak is machine-dependent → a MEASUREMENT (relocated out of the
        // claim by finalize). The build is still O(reference) RAM; this records the cost,
        // it does NOT claim the budget was honored.
        manifest.record_measurement("peak_rss_bytes", peak.to_string());
        manifest.finalize();
        let dest = write_manifest(&output, &manifest)
            .with_context(|| format!("failed to write index receipt for {}", output.display()))?;
        eprintln!("wrote reproducibility receipt: {}", dest.display());
    }
    Ok(())
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test --test chain index_writes_a_self_verifying_receipt`
Expected: PASS — the sidecar exists and `verify` exits 0.

- [ ] **Step 5: Run the existing suite to confirm no regression**

Run: `cargo test --test index_cli && cargo test --test reproduce`
Expected: PASS. (`run_index` now also emits one extra stderr line + a sidecar file. If any `index_cli` assertion checks an exact directory listing or snapshots `index` stderr, update it to include the new `wrote reproducibility receipt:` line — the line is deterministic.)

- [ ] **Step 6: Commit**

```bash
git add src/main.rs tests/chain.rs
git commit -m "feat(index): emit a content-addressed receipt (chain root)

rosalind index now writes a <output>.manifest.json sidecar mirroring the
variants/somatic path: --reference records the FASTA file hash (chain root),
--output records the .idx file hash (the edge a downstream variants --index
input resolves to by construction). Realized build peak is a measurement, not
a budget claim. Additive; schema unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 2: the index→variants edge resolves by construction

**Files:**
- Modify: `tests/chain.rs` (append a test)

- [ ] **Step 1: Write the edge-resolution guard test**

Append to `tests/chain.rs`:

```rust
#[test]
fn index_output_hash_equals_variants_index_input_hash() {
    use rosalind::provenance::RunManifest;

    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");

    let made = run(&[
        "variants", "--index", idx.to_str().unwrap(),
        "--alignments", bam.to_str().unwrap(), "-o", vcf.to_str().unwrap(),
    ]);
    assert!(made.status.success(), "variants: {}", String::from_utf8_lossy(&made.stderr));

    let index_m = RunManifest::from_canonical_json(
        &std::fs::read_to_string(d.join("ref.idx.manifest.json")).unwrap(),
    )
    .unwrap();
    let variants_m = RunManifest::from_canonical_json(
        &std::fs::read_to_string(format!("{}.manifest.json", vcf.display())).unwrap(),
    )
    .unwrap();

    // The index's recorded .idx output hash == the variants' recorded --index input hash.
    let index_out = &index_m.outputs[0].blake3;
    let variants_index_in = variants_m
        .inputs
        .iter()
        .find(|f| f.path.ends_with("ref.idx"))
        .expect("variants receipt records the index as an input")
        .blake3
        .clone();
    assert_eq!(
        *index_out, variants_index_in,
        "the provenance edge must resolve: index outputs[0] == variants inputs[--index]"
    );

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run the test**

Run: `cargo test --test chain index_output_hash_equals_variants_index_input_hash`
Expected: PASS — both hashes are `blake3_file(ref.idx)`, equal by construction now that `index` emits `outputs[]`. (This guard would have been impossible before Task 1, when `index` wrote no receipt.)

- [ ] **Step 3: Commit**

```bash
git add tests/chain.rs
git commit -m "test(chain): the index->variants edge resolves by content hash

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## PR2 — `chain verify <dir>`: the offline DAG walker

### Task 3: pure `walk_chain` in the receipt crate

**Files:**
- Create: `crates/receipt/src/chain.rs`
- Modify: `crates/receipt/src/lib.rs` (add `pub mod chain;` + a `pub use`, near the existing `mod command; pub use command::CommandCapture;` at ≈ line 24)
- Test: `crates/receipt/src/chain.rs` (`#[cfg(test)] mod tests`)

- [ ] **Step 1: Wire the module into the crate**

In `crates/receipt/src/lib.rs`, after `pub use command::CommandCapture;` (≈ line 25), add:

```rust
mod chain;
pub use chain::{walk_chain, ChainEdge, ChainNode, ChainReport, EdgeStatus};
```

- [ ] **Step 2: Write `crates/receipt/src/chain.rs` with the types, an unimplemented walk, and failing unit tests**

```rust
//! Pure, `std`-only provenance-DAG walk over a set of run receipts. No filesystem, no
//! htslib — so it is unit-testable and wasm-portable (the future `chain confirm` board /
//! browser verifier reuse it). A node is a receipt (id = its claim `content_hash`); an
//! edge is one input operand, classified by whether its content hash is produced by
//! another node in the set.

use std::collections::HashMap;

use crate::{FileHash, RunManifest};

/// Input operand flags that MUST resolve to a producing node. An unresolved one is a
/// broken chain (a missing/mismatched upstream receipt). Everything else (reads,
/// alignments, reference FASTA) is an external source: integrity-verified by its
/// recorded hash, never a chain failure.
const EXPECTED_INTERNAL: &[&str] = &["--index"];

/// How one input operand resolved against the receipt set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EdgeStatus {
    /// The input hash equals some node's output hash — an internal provenance edge.
    Resolved { parent_id: String },
    /// Unresolved, but an expected-external operand (reference/alignments/reads).
    External,
    /// Unresolved AND an expected-internal operand (`--index`) — breaks the chain.
    Broken,
}

/// One classified input operand of one node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainEdge {
    pub child_id: String,
    pub child_subcommand: String,
    pub flag: String,
    pub input_blake3: String,
    pub status: EdgeStatus,
}

/// A receipt as a DAG node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainNode {
    pub id: String,
    pub subcommand: String,
    /// `Some(true/false)` if a claim self-hash is recorded and matches/mismatches;
    /// `None` for a pre-self-hash receipt.
    pub self_hash: Option<bool>,
}

/// The result of walking the set: nodes, classified edges, and the overall verdict.
#[derive(Debug, Clone)]
pub struct ChainReport {
    pub nodes: Vec<ChainNode>,
    pub edges: Vec<ChainEdge>,
    /// `true` iff no node fails its self-hash AND no edge is `Broken`.
    pub intact: bool,
}

impl ChainReport {
    /// A compact, dependency-free JSON summary for `--json` / a scheduler.
    pub fn to_json(&self) -> String {
        let resolved = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        let external = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::External))
            .count();
        let broken = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Broken))
            .count();
        format!(
            "{{\"intact\":{},\"nodes\":{},\"edges_resolved\":{},\"edges_external\":{},\"edges_broken\":{}}}",
            self.intact,
            self.nodes.len(),
            resolved,
            external,
            broken
        )
    }
}

/// Recover `(flag, input_blake3)` pairs from a recorded `command` string: each
/// `@in:<hash>` token is preceded by its operand flag.
fn input_operands(command: &str) -> Vec<(String, String)> {
    let toks: Vec<&str> = command.split(' ').collect();
    let mut out = Vec::new();
    for (i, t) in toks.iter().enumerate() {
        if let Some(h) = t.strip_prefix("@in:") {
            let flag = if i > 0 { toks[i - 1] } else { "?" };
            out.push((flag.to_string(), h.to_string()));
        }
    }
    out
}

/// Operands for a node: from its recorded `command` when present (carries the flags),
/// else a flag-less fallback over `inputs[]` (every input treated as external).
fn operands_for(m: &RunManifest) -> Vec<(String, String)> {
    match m.params.get("command") {
        Some(c) if !c.is_empty() => input_operands(c),
        _ => m
            .inputs
            .iter()
            .map(|f: &FileHash| ("?".to_string(), f.blake3.clone()))
            .collect(),
    }
}

/// Walk the receipt set as a provenance DAG. Pure: no I/O.
pub fn walk_chain(receipts: &[RunManifest]) -> ChainReport {
    let ids: Vec<String> = receipts.iter().map(|m| m.content_hash()).collect();

    // output content hash -> producing node id (first producer wins).
    let mut producer: HashMap<String, String> = HashMap::new();
    for (m, id) in receipts.iter().zip(&ids) {
        for o in &m.outputs {
            producer.entry(o.blake3.clone()).or_insert_with(|| id.clone());
        }
    }

    let nodes: Vec<ChainNode> = receipts
        .iter()
        .zip(&ids)
        .map(|(m, id)| ChainNode {
            id: id.clone(),
            subcommand: m.subcommand.clone(),
            self_hash: m.self_hash_ok(),
        })
        .collect();

    let mut edges = Vec::new();
    for (m, id) in receipts.iter().zip(&ids) {
        for (flag, hash) in operands_for(m) {
            let status = if let Some(parent) = producer.get(&hash) {
                EdgeStatus::Resolved {
                    parent_id: parent.clone(),
                }
            } else if EXPECTED_INTERNAL.contains(&flag.as_str()) {
                EdgeStatus::Broken
            } else {
                EdgeStatus::External
            };
            edges.push(ChainEdge {
                child_id: id.clone(),
                child_subcommand: m.subcommand.clone(),
                flag,
                input_blake3: hash,
                status,
            });
        }
    }

    let tampered = nodes.iter().any(|n| n.self_hash == Some(false));
    let broken = edges.iter().any(|e| matches!(e.status, EdgeStatus::Broken));
    ChainReport {
        nodes,
        edges,
        intact: !tampered && !broken,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FileHash, RunManifest};

    /// Build a finalized manifest with a recorded `command`, inputs, and outputs.
    fn mk(
        sub: &str,
        command: &str,
        inputs: &[(&str, &str)],
        outputs: &[(&str, &str)],
    ) -> RunManifest {
        let mut m = RunManifest::new(sub);
        m.params.insert("command".to_string(), command.to_string());
        m.inputs = inputs
            .iter()
            .map(|(p, h)| FileHash {
                path: p.to_string(),
                blake3: h.to_string(),
            })
            .collect();
        m.outputs = outputs
            .iter()
            .map(|(p, h)| FileHash {
                path: p.to_string(),
                blake3: h.to_string(),
            })
            .collect();
        m.finalize();
        m
    }

    fn index_node() -> RunManifest {
        mk(
            "index",
            "index --reference @in:rh --output @out:ih",
            &[("ref.fa", "rh")],
            &[("ref.idx", "ih")],
        )
    }

    fn variants_node(index_hash: &str) -> RunManifest {
        mk(
            "variants",
            &format!("variants --index @in:{index_hash} --alignments @in:bh -o @out:vh"),
            &[("ref.idx", index_hash), ("s.bam", "bh")],
            &[("calls.vcf", "vh")],
        )
    }

    #[test]
    fn resolves_the_index_edge_and_marks_external_inputs() {
        let report = walk_chain(&[index_node(), variants_node("ih")]);
        assert!(report.intact, "a complete chain is intact");
        // The --index edge resolves; --alignments (bh) and --reference (rh) are external.
        let resolved = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        let external = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::External))
            .count();
        assert_eq!(resolved, 1, "exactly the --index edge resolves");
        assert_eq!(external, 2, "--alignments and --reference are external sources");
        assert_eq!(report.to_json(),
            "{\"intact\":true,\"nodes\":2,\"edges_resolved\":1,\"edges_external\":2,\"edges_broken\":0}");
    }

    #[test]
    fn an_unresolved_index_input_is_broken() {
        // variants references an --index hash no node produces.
        let report = walk_chain(&[index_node(), variants_node("MISSING")]);
        assert!(!report.intact, "a dangling --index edge breaks the chain");
        assert!(report
            .edges
            .iter()
            .any(|e| e.flag == "--index" && e.status == EdgeStatus::Broken));
    }

    #[test]
    fn an_unresolved_external_input_does_not_break_the_chain() {
        // The index node alone: its --reference (rh) resolves to nothing, but --reference
        // is an external source, so the chain is still intact.
        let report = walk_chain(&[index_node()]);
        assert!(report.intact, "an unresolved external source is integrity-only, not broken");
        assert!(report
            .edges
            .iter()
            .any(|e| e.flag == "--reference" && e.status == EdgeStatus::External));
    }

    #[test]
    fn a_tampered_node_breaks_the_chain() {
        let mut tampered = index_node();
        // Edit a claim field AFTER finalize → the recorded manifest_blake3 no longer matches.
        tampered.params.insert("total_bp".to_string(), "999999".to_string());
        let report = walk_chain(&[tampered, variants_node("ih")]);
        assert!(!report.intact, "a self-hash mismatch breaks the chain");
        assert!(report.nodes.iter().any(|n| n.self_hash == Some(false)));
    }
}
```

- [ ] **Step 3: Run the unit tests to verify they pass**

Run: `cargo test -p rosalind-receipt chain::`
Expected: PASS — all four `chain::tests::*` pass. (If you staged this as red-first by stubbing `walk_chain` with `unimplemented!()`, restore the implementation above.)

- [ ] **Step 4: Confirm the crate still builds clean (wasm-portability guard: no new deps)**

Run: `cargo build -p rosalind-receipt`
Expected: PASS — `chain.rs` uses only `std::collections::HashMap` + the crate's own types.

- [ ] **Step 5: Commit**

```bash
git add crates/receipt/src/chain.rs crates/receipt/src/lib.rs
git commit -m "feat(receipt): pure walk_chain provenance-DAG walker

std-only, wasm-portable: classify each input operand resolved/external/broken
against the receipt set, check every node's self-hash, report an intact verdict.
EXPECTED_INTERNAL={--index}; reads/alignments/reference are external sources.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 4: the `chain verify` CLI noun

**Files:**
- Modify: `src/main.rs` (add `Chain` to `Commands` ≈ line 236; add a `ChainAction` enum next to it; add a dispatch arm ≈ line 540; add `run_chain_verify`)
- Test: `tests/chain.rs` (append the INTACT integration test)

- [ ] **Step 1: Write the failing integration test**

Append to `tests/chain.rs`:

```rust
#[test]
fn chain_verify_is_intact_on_a_real_index_variants_chain() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants", "--index", idx.to_str().unwrap(),
        "--alignments", bam.to_str().unwrap(), "-o", vcf.to_str().unwrap(),
    ])
    .status
    .success());

    // d now holds ref.idx.manifest.json + calls.vcf.manifest.json.
    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success(),
        "chain verify must exit 0 (INTACT). stdout:\n{stdout}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(stdout.contains("CHAIN INTACT"), "verdict line: {stdout}");

    let j = run(&["chain", "verify", d.to_str().unwrap(), "--json"]);
    let js = String::from_utf8_lossy(&j.stdout);
    assert!(js.contains("\"intact\":true"), "json: {js}");
    assert!(js.contains("\"edges_resolved\":1"), "json: {js}");

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --test chain chain_verify_is_intact_on_a_real_index_variants_chain`
Expected: FAIL — `clap` rejects the unknown `chain` subcommand (non-zero exit, stdout has no `CHAIN INTACT`).

- [ ] **Step 3: Add the `Chain` command + `ChainAction` enum**

In `src/main.rs`, in the `Commands` enum, after the `Index { … }` variant (≈ line 247), add:

```rust
    /// Walk a directory of receipts as a provenance DAG and verify it offline.
    Chain {
        #[command(subcommand)]
        action: ChainAction,
    },
```

Immediately after the `Commands` enum's closing brace, add the action enum:

```rust
#[derive(Subcommand, Debug)]
enum ChainAction {
    /// Verify every node self-hashes and every internal edge resolves by content hash.
    Verify {
        /// Directory of `*.manifest.json` receipts to walk.
        dir: PathBuf,
        /// Emit a compact JSON report instead of human-readable lines.
        #[arg(long)]
        json: bool,
    },
}
```

- [ ] **Step 4: Add the dispatch arm**

In `src/main.rs`, in the `match cli.command { … }` block, after the `Commands::Index { … } => …` arm (≈ line 540), add:

```rust
        Commands::Chain { action } => match action {
            ChainAction::Verify { dir, json } => run_chain_verify(dir, json)?,
        },
```

- [ ] **Step 5: Implement `run_chain_verify`**

In `src/main.rs`, add the function (e.g. just after `run_verify`, ≈ line 979):

```rust
/// Walk a directory of receipts as a provenance DAG: every node self-hashes and every
/// expected-internal edge (`--index`) resolves by content hash. External inputs
/// (reads/alignments/reference FASTA) are integrity-verified, not byte-reproduced.
/// Exit 0 = CHAIN INTACT, 5 = CHAIN BROKEN (matching `verify`).
fn run_chain_verify(dir: PathBuf, json: bool) -> Result<()> {
    use rosalind::provenance::{walk_chain, EdgeStatus, RunManifest};
    use std::collections::HashMap;

    let mut receipts: Vec<RunManifest> = Vec::new();
    for entry in std::fs::read_dir(&dir)
        .with_context(|| format!("failed to read directory {}", dir.display()))?
    {
        let path = entry?.path();
        let is_manifest = path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.ends_with(".manifest.json"))
            .unwrap_or(false);
        if !is_manifest {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        // Skip anything that is not a parseable run manifest (e.g. a stray file).
        if let Ok(m) = RunManifest::from_canonical_json(&text) {
            receipts.push(m);
        }
    }
    if receipts.is_empty() {
        bail!("no receipts (*.manifest.json) found in {}", dir.display());
    }

    let report = walk_chain(&receipts);

    if json {
        println!("{}", report.to_json());
    } else {
        let sub: HashMap<&str, &str> = report
            .nodes
            .iter()
            .map(|n| (n.id.as_str(), n.subcommand.as_str()))
            .collect();
        let tampered = report
            .nodes
            .iter()
            .filter(|n| n.self_hash == Some(false))
            .count();
        if tampered == 0 {
            println!("chain: {} nodes, all self-hash OK", report.nodes.len());
        } else {
            println!(
                "chain: {} nodes, {tampered} TAMPERED (self-hash mismatch)",
                report.nodes.len()
            );
        }
        for e in &report.edges {
            match &e.status {
                EdgeStatus::Resolved { parent_id } => {
                    let p = sub.get(parent_id.as_str()).copied().unwrap_or("?");
                    println!("edge: {}  {}-->  {p}  [resolved]", e.child_subcommand, e.flag);
                }
                EdgeStatus::External => {
                    println!(
                        "edge: {}  {}-->  (external)  [integrity-only]",
                        e.child_subcommand, e.flag
                    );
                }
                EdgeStatus::Broken => {
                    println!(
                        "edge: {}  {}-->  (unresolved)  [BROKEN]",
                        e.child_subcommand, e.flag
                    );
                }
            }
        }
        let resolved = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        if report.intact {
            println!(
                "VERDICT: CHAIN INTACT ({} nodes, {resolved} internal edges resolve)",
                report.nodes.len()
            );
        } else {
            println!("VERDICT: CHAIN BROKEN");
        }
    }

    if report.intact {
        Ok(())
    } else {
        std::process::exit(5);
    }
}
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `cargo test --test chain chain_verify_is_intact_on_a_real_index_variants_chain`
Expected: PASS — `chain verify` exits 0, prints `CHAIN INTACT`, and `--json` reports `"intact":true,"edges_resolved":1`.

- [ ] **Step 7: Commit**

```bash
git add src/main.rs tests/chain.rs
git commit -m "feat(cli): rosalind chain verify <dir> — offline provenance-DAG walk

A new 'chain' noun whose 'verify' action loads *.manifest.json from a dir and
walks it via walk_chain: CHAIN INTACT (exit 0) when every node self-hashes and
every --index edge resolves by content hash; CHAIN BROKEN (exit 5) otherwise.
--json emits a compact summary. Leaves room for chain confirm/show.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 5: failure-mode integration gates

**Files:**
- Modify: `tests/chain.rs` (append two tests)

- [ ] **Step 1: Write the broken-edge and tamper integration tests**

Append to `tests/chain.rs`:

```rust
#[test]
fn chain_verify_breaks_when_the_index_receipt_is_tampered() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants", "--index", idx.to_str().unwrap(),
        "--alignments", bam.to_str().unwrap(), "-o", vcf.to_str().unwrap(),
    ])
    .status
    .success());

    // Tamper with the index receipt's claim: flip a digit of total_bp without
    // recomputing manifest_blake3 → the node self-hash must fail.
    let mpath = d.join("ref.idx.manifest.json");
    let text = std::fs::read_to_string(&mpath).unwrap();
    let tampered = text.replacen("\"total_bp\":\"", "\"total_bp\":\"9", 1);
    assert_ne!(text, tampered, "the receipt must contain a total_bp field to tamper");
    std::fs::write(&mpath, tampered).unwrap();

    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(5), "a tampered node must exit 5");
    assert!(
        String::from_utf8_lossy(&out.stdout).contains("CHAIN BROKEN"),
        "stdout: {}",
        String::from_utf8_lossy(&out.stdout)
    );

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn chain_verify_reports_the_bam_input_as_external_not_a_failure() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants", "--index", idx.to_str().unwrap(),
        "--alignments", bam.to_str().unwrap(), "-o", vcf.to_str().unwrap(),
    ])
    .status
    .success());

    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    // The BAM alignments input has no producing receipt → external/integrity-only,
    // and it must NOT fail the chain.
    assert!(out.status.success(), "external inputs must not break the chain: {stdout}");
    assert!(
        stdout.contains("--alignments-->  (external)  [integrity-only]"),
        "the BAM edge must be reported external: {stdout}"
    );

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run the new tests**

Run: `cargo test --test chain chain_verify_breaks_when_the_index_receipt_is_tampered chain_verify_reports_the_bam_input_as_external_not_a_failure`
Expected: PASS — the tampered node yields exit 5 + `CHAIN BROKEN`; the BAM input is reported external and the chain stays INTACT.

- [ ] **Step 3: Full suite + lint + format gates**

Run:
```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```
Expected: PASS on all three. (`cargo fmt --check` must be clean; if it reports diffs, run `cargo fmt` and re-commit.)

- [ ] **Step 4: Commit**

```bash
git add tests/chain.rs
git commit -m "test(chain): tamper -> CHAIN BROKEN (exit 5); BAM input -> external

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**1. Spec coverage** (each spec section → a task):
- §4 index receipt (claim params, `--reference`/`--output` operands, peak as measurement, `write_manifest` sidecar, additive/schema-unchanged) → **Task 1**. The two-hash distinction (`outputs[0]` = file edge; `reference_blake3` = informational) → Task 1 code + **Task 2** guard.
- §5.1 CLI `chain verify <dir>` + `--json` → **Task 4**.
- §5.2 pure walker in `crates/receipt` → **Task 3**.
- §5.3 edge classification (resolved / external / broken; `EXPECTED_INTERNAL={--index}`; flag recovery from `command`) → **Task 3** (`operands_for`/`input_operands`/`walk_chain`).
- §5.4 verdicts/exit codes (0 / 5) → `run_chain_verify` (**Task 4**), gated by **Task 5**.
- §5.5 output shape + `--json` → **Task 4** (`run_chain_verify` + `ChainReport::to_json`).
- §5.6 honesty contract (external = integrity-only, never a failure) → **Task 5** (`…bam_input_as_external…`).
- §6 tests: unit resolved/external/broken/tampered → **Task 3**; integration index-self-verifies (1), edge-equality (2), CHAIN INTACT (3), tamper→BROKEN (4), broken edge (covered by Task 3 unit + the tamper integration), external integrity-only (6) → **Tasks 1,2,4,5**.
- §8 guardrails (no schema bump, keep stdout `IndexBuildReport`, std-only walker) → honored in Task 1 (additive, `finalize` keeps schema 5) + Task 3 (Step 4 build guard).

**2. Placeholder scan:** No TBD/TODO/"handle edge cases"/"similar to". Every step has the literal code or the exact command + expected result.

**3. Type consistency:** `EdgeStatus`/`ChainEdge`/`ChainNode`/`ChainReport`/`walk_chain`/`ChainReport::to_json` are defined in Task 3 and used identically in Task 4 (`run_chain_verify`) and exported from `lib.rs` in Task 3 Step 1. `self_hash: Option<bool>` is consistent across `ChainNode` definition, `walk_chain`, and `run_chain_verify`'s `== Some(false)` checks. The `ChainAction` enum (Task 4 Step 3) matches the dispatch arm (Step 4). `write_manifest`/`blake3_hex`/`CommandCapture`/`RunManifest` import paths match the existing `somatic` block.

**Note for the implementer:** PR1 = Tasks 1–2 (land/PR first; immediately bankable). PR2 = Tasks 3–5. The `index` receipt may surface a stderr/dir-listing assertion in `tests/index_cli.rs` (Task 1 Step 5) — update that snapshot to include the deterministic `wrote reproducibility receipt:` line if so.
