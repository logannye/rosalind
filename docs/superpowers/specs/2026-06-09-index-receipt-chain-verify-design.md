# Design: `index` receipt + `chain verify` (the rooted provenance DAG)

**Status:** Approved design — 2026-06-09. Companion to [`docs/ROADMAP.md`](../../ROADMAP.md)
(this realizes the under-named **P3.1** keystone) and [`CONTRACT.md`](../../../CONTRACT.md).
Audience: the implementer of the index-receipt + chain-walker work.

---

## 1. The goal, in one sentence

Give the `rosalind index` build a **content-addressed receipt** so the content-hash edges
that downstream receipts *already record* resolve into a real, offline-walkable **provenance
DAG** — making the "reproducibility web, no server" the brand markets **true at its root**.

## 2. Why this, why now

`rosalind index` is the **root of every provenance chain** yet is the **only producing
subcommand that emits no receipt**: `run_index` (`src/main.rs:639`) `print!`s a human
`IndexBuildReport` to stdout and writes nothing machine-verifiable. Meanwhile
`variants`/`somatic`/`features` all go through `RunManifest::new + CommandCapture +
finalize + write` and already record the index they consume **by content hash**
(`variants` records `cmd.input("--index", &index_path)` → `blake3_file(.idx)`). So the
edges exist; the **root node is missing**, and the chain dangles.

Every primitive needed to compose a verifiable provenance graph already ships and is
adversarially tested:

- `content_hash()` (`crates/receipt/src/lib.rs:246`) — the cross-machine claim address (node id).
- `FileHash { path, blake3 }` entries in `inputs[]`/`outputs[]` — the **edges**.
- `verify_receipt` (`lib.rs:679`) — the shared per-node tamper check (also used by CLI `verify`,
  `reproduce`, and the deployed WASM verifier).
- `CommandCapture` (`crates/receipt/src/command.rs`) — `cmd.input/​output` call `blake3_file(path)`,
  so an `index` output digest and a `variants` `--index` input digest over the **same `.idx`
  file are bit-identical by construction**.

What is missing is **one sidecar writer** (PR1) and **one graph walker** (PR2).

**The differentiation:** nf-core / WDL / Snakemake / Nextflow track provenance by *path +
timestamp* over *non-deterministic* callers — they structurally cannot offer "this whole graph
re-derives by content hash, verifiable offline, no server." Rosalind can, because determinism +
the content-addressed claim are properties of every run.

## 3. Scope

**In scope (this design):** two sequenced PRs, sandboxed entirely within `~/rosalind`.

- **PR1 — `index` receipt** *(effort: S)*: `rosalind index` writes a content-addressed sidecar
  manifest next to the `.idx`, mirroring the `variants` path.
- **PR2 — `chain verify <dir>`** *(effort: M)*: a new `chain` subcommand whose `verify` action
  auto-discovers and walks the receipt DAG in a directory, checking node self-hash + edge
  resolution, entirely offline.

**Out of scope (named, deferred):**

- Whole-pipeline **byte-reproduction** of the chain (the chain proves *integrity + edge
  resolution*, not reproduction).
- BAM/bgzf byte-reproduction (`reproduce` already reports INCONCLUSIVE for bgzf; the chain
  surfaces BAM inputs as integrity-only).
- Memory-bounded index **construction** (still `O(reference)` RAM — Phase D). The receipt records
  realized build peak as a *measurement*; it never implies the build budget was honored.
- `chain confirm` (multi-party `.repro.json` cert board) and `chain show --dot` (DAG render) —
  future work the new `chain` noun deliberately leaves room for.

## 4. PR1 — `index` emits a content-addressed receipt

**Where:** `run_index` (`src/main.rs:639`). **Additive** — the existing stdout `IndexBuildReport`
and the build-memory telemetry (stderr) stay unchanged. We add a receipt after the index is
written, reusing values already in scope (`output`, `total_bp`, `reference_blake3`, `peak`).

**The receipt construction** (mirrors `variants`, `main.rs:2236`):

```rust
let mut manifest = RunManifest::new("index");
let mut cmd = CommandCapture::new("index");
cmd.input("--reference", &reference)?;          // blake3 of the FASTA *file* → the chain ROOT source
if let Some(mb) = memory_budget_mb {
    cmd.opt("--memory-budget-mb", mb);
}
cmd.output("--output", &output)?;               // blake3 of the .idx *file* → THE chain edge
cmd.record_into(&mut manifest);                 // sets inputs/outputs; mode = "reference"
manifest.params.insert("reference_blake3".into(), blake3_hex(&reference_blake3)); // "what genome" id
manifest.params.insert("total_bp".into(), total_bp.to_string());
manifest.record_measurement("peak_rss_bytes", peak.to_string());  // machine-dependent → out of claim
manifest.finalize();
// sidecar destination = <output>.manifest.json (mirror the variants `-o` → `-o.manifest.json` rule)
std::fs::write(&dest, manifest.to_canonical_json())?;
eprintln!("wrote reproducibility receipt: {}", dest.display());
```

**The two-hash distinction (must be honored exactly):**

| Hash | Source | Role |
|---|---|---|
| `outputs[0].blake3` (via `cmd.output("--output", …)`) | `blake3_file(.idx)` — the index *file* bytes | **The chain edge.** Identical by construction to the `variants` `inputs[--index].blake3`. |
| `inputs[0].blake3` (via `cmd.input("--reference", …)`) | `blake3_file(reference.fa)` — the FASTA *file* bytes | The chain's **root source** (external; no producing receipt). |
| `reference_blake3` param | `blake3::hash(index.reference())` — the in-memory normalized sequence | **Informational** "what genome" identity, stable across FASTA reformatting. **NOT** an edge. |

**Recorded command flags:** the `index` subcommand's real flags are `--reference` and
`--output` (`-o`), so the recorded operands replay correctly.

**Claim vs measurement:** `reference_blake3` and `total_bp` are deterministic → **claim** params.
`peak_rss_bytes` is machine-dependent → **measurement** (relocated out of the claim by
`finalize()`, keeping `content_hash()` cross-machine stable).

**Back-compat:** `MANIFEST_SCHEMA_VERSION` stays **5** — `index` is simply a new subcommand value.
Pre-existing indexes have no sidecar; that is fine (they pre-date the feature).

**Honest framing:** the index now carries the same verifiable receipt the call paths do, and its
bytes are reproducible — **not** any memory or √t claim. The build is still `O(reference)` RAM.

## 5. PR2 — `chain verify <dir>`: the offline DAG walker

### 5.1 CLI surface

A new top-level **`chain`** noun (clap), with a `verify` action:

```
rosalind chain verify <dir>          # auto-discover + walk the whole receipt DAG in <dir>
rosalind chain verify <dir> --json   # structured ChainReport for a scheduler/CI
```

The `chain` namespace deliberately leaves room for future `chain confirm <certs-dir>`
(the multi-party `.repro.json` board) and `chain show <dir> --dot` (DAG render).

### 5.2 Architecture — library-first, wasm-ready

The **pure walk** lives in the `receipt` leaf crate (`crates/receipt/src/chain.rs`) — `std` +
the manifest model only, **no htslib** — so the future `chain confirm` board and the WASM
verifier can reuse it. `main.rs` does file I/O (load `*.manifest.json` from `<dir>`) and CLI
rendering.

```rust
// crates/receipt/src/chain.rs — pure, fs-free, unit-testable, wasm-portable
pub struct ChainReport { /* nodes, edges, verdict, counts */ }
pub enum EdgeStatus { Resolved { parent_id: String }, External, Broken }
pub fn walk_chain(receipts: &[RunManifest]) -> ChainReport;
```

### 5.3 The walk

1. **Load** every `*.manifest.json` in `<dir>`; parse via `RunManifest::from_canonical_json`.
   Skip `.repro.json` and unparseable files (they are not chain nodes).
2. Each parsed manifest = a **node**, id = `content_hash()`. Build an `output_blake3 → node` map.
3. **Per-node self-hash** via the existing tamper check (offline; needs no data files).
4. **Edge classification** (operand-aware — the honest, strong guarantee). For each node, for each
   `inputs[i]` (and its operand flag, recoverable from the recorded `command` token order):
   - input `blake3` **matches** some node's `outputs[*].blake3` → **`Resolved`** internal edge
     (e.g. `variants --index → index`).
   - **unresolved `--index`** input → **`Broken`** (the index receipt is missing or its `.idx`
     hash does not match) → fails the chain.
   - **unresolved `--reference` / `--alignments` / `--reads`** → **`External`** source —
     integrity-verified by recorded hash, **never** fails the chain.

   The expected-internal flag set is currently `{ --index }`; it is a small, explicit, extensible
   classification (a future `align --index` output would extend it).

### 5.4 Verdicts + exit codes (match `verify`'s convention)

- **`CHAIN INTACT`** — exit **0**: every node self-hashes **and** every expected-internal edge
  is `Resolved`.
- **`CHAIN BROKEN`** — exit **5**: any node fails self-hash (tampered) **or** any expected-internal
  (`--index`) edge is unresolved. (Exit 5 matches `verify`'s "receipt problem" code.)

### 5.5 Output

Human (the approved shape) and `--json` (the structured `ChainReport`):

```
chain: 3 nodes, all self-hash OK
edge: sample.vcf  --index-->  ref.idx         [resolved]
edge: sample.vcf  --alignments-->  (external) [integrity-only]
root: ref.idx  --reference-->  genome.fa      [external source]
VERDICT: CHAIN INTACT (3 nodes, 2 internal edges resolve, root = genome.fa)
```

### 5.6 What `CHAIN INTACT` asserts — precisely (the honesty contract)

> Every node is **tamper-evident** and every **internal provenance edge resolves by content
> hash** — offline, no re-run.

It does **not** claim the pipeline byte-reproduces; external reads/alignments are
**integrity-verified, not reproduced**; bgzf/BAM remains `reproduce`-INCONCLUSIVE by design.
Reserve the phrase "reproducibility web" for the future multi-party `.repro.json` layer; this is
a **verifiable provenance graph**.

## 6. Testing

**Unit (receipt crate, pure, no fs)** — `walk_chain` over synthetic `RunManifest`s with
hand-set `inputs`/`outputs`:
- a resolved `--index` edge → `Resolved`, `CHAIN INTACT`.
- an unresolved `--index` edge → `Broken`, `CHAIN BROKEN`.
- an unresolved `--alignments` / `--reference` edge → `External`, does not fail.
- a tampered node (claim mutated so the self-hash mismatches) → `CHAIN BROKEN`.

**Integration (`tests/chain.rs`, bundled toy data)** — a real reference → index → variants chain:
1. `index` writes a sidecar that self-hashes (`verify: OK`).
2. The `variants` receipt's `inputs[--index].blake3` **==** the `index` receipt's
   `outputs[0].blake3` (the edge resolves by construction).
3. `chain verify <dir>` → `CHAIN INTACT`, exit 0, expected node/edge counts.
4. **Tamper:** edit a byte of the index receipt's claim → node self-hash fails → `CHAIN BROKEN`
   (exit 5).
5. **Broken edge:** a `variants` receipt whose `--index` hash matches no node → `CHAIN BROKEN`
   (exit 5).
6. **External integrity-only:** the BAM `--alignments` input is reported `External` and does
   **not** fail the chain.

**Determinism/golden:** the existing byte-identity gates already cover the `.idx` and VCF; the new
receipt is canonical JSON (sorted keys, no timestamps) so it is byte-stable run-to-run by
construction.

## 7. Build sequence

1. **PR1** — `index` receipt + `tests/chain.rs` cases (1)+(2). Immediately bankable: an index
   build becomes `verify`-able on its own, and the edge is provably constructible.
2. **PR2** — `crates/receipt/src/chain.rs` `walk_chain` + unit tests; `chain verify` CLI in
   `main.rs`; `tests/chain.rs` cases (3)–(6). The differentiation lands here.

## 8. Non-goals / guardrails (do not drift)

- Do **not** imply the chain byte-reproduces the pipeline; bound the claim to edge-resolution +
  tamper-evidence, offline.
- Do **not** imply the index **build** is memory-bounded; the receipt records realized peak as a
  measurement only.
- Do **not** bump the schema version (additive new subcommand value).
- Do **not** remove the existing stdout `IndexBuildReport` (the receipt is additive).
- Keep the walker `std`-only (no htslib) so it stays wasm-portable for the future board.
