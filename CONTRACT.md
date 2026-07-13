# The memory contract

**Memory is a contract, not a hope.**

Schedulers can impose a ceiling, but a reusable analyzer still needs a declared
working-set model and evidence about what happened inside that ceiling. Rosalind
treats that relationship as a contract:

> **Rosalind predicts whether a declared job fits, refuses unsound enforcement,
> records the realized peak and assurance in a verifiable receipt, and can require
> an existing Linux cgroup ceiling when OS-backed enforcement is mandatory.**

(It does *not* claim to "never refuse": when a budget is genuinely too small the run declines cleanly
rather than crashing. Cooperative RSS polling is not a kernel memory sandbox; use
`--require-os-limit` for a detected cgroup-v2 ceiling. Graceful degrade-don't-die — sliding down a space/time curve to finish anyway — is
the Phase-D research direction, not a present claim.)

## The four verbs

The contract applies to whole-genome and selected `variants`, `features`, and
registered `analyze` paths. Preferred `.rref` and compatible `.idx` references use
the same `plan`/`--enforce`/`verify` envelope.

### 1. Declare

State the RAM you have. `--memory-budget-mb N` on `variants`, `--budget-mb N` on `plan`/`verify`.

### 2. Predict — `rosalind plan`

Ask *before committing a byte* whether the job fits. `plan` reads only reference metadata (plus your
declared depth/read-length assumptions) — it never opens the BAM:

```bash
rosalind plan --reference-pack genome.rref --max-depth 1000 --max-read-len 250 --budget-mb 2048
```

It prints a breakdown — reference decode + active set @ max-depth + engine overhead, atop a measured
process baseline — and a verdict: `[FITS]` or `[REFUSE]`.

### 3. Honor — `rosalind variants … --enforce`

```bash
rosalind variants --reference-pack genome.rref --alignments sample.sorted.bam \
  --memory-budget-mb 2048 --enforce -o sample.vcf
```

With `--enforce` (cooperative assurance):

- **predicted peak > budget → refuse up front** (exit **3**), before doing any work, with an actionable
  message (raise the budget, lower `--max-depth`, or drop `--enforce`);
- **realized peak > budget → fail loud** (exit **4**) after flushing the governed
  artifact to `<output>.partial` and sealing a breach receipt that names that partial
  artifact; the successful output name is never exposed;
- a **runtime governor** polls process RSS *during* the run and fails loud the moment the realized peak
  crosses the budget — exit **4** with the partial output + a `governor=tripped`, `contract_verdict=over`
  receipt — so a misprediction is caught mid-run instead of by a silent kernel OOM. The receipt also
  records the realized RSS residual (`rss_residual_bytes`) against the assumed margin
  (`io_rss_overhead_assumed_bytes`);
- otherwise the run completes within budget.

Without `--enforce`, the budget is **record-only**: the run always completes and the verdict is recorded in
the receipt. The active read set is capped at `--max-depth` (default 1000; `0` = uncapped) by an
**unbiased** content-hash downsampling — it bounds the working set without biasing allele balance, so a deep
variant is *not* silently dropped. Output changes only at sites deeper than the cap, and the dropped-read
count is surfaced (stderr + the receipt's `over_max_depth`).

`--require-os-limit` (used together with `--enforce`) raises the assurance level on
Linux: Rosalind requires the process to already be inside cgroup v2 with a finite
`memory.max` no greater than the declared budget. It never attempts privileged
cgroup creation. Missing, unlimited, or overly broad limits refuse before the output
is opened and include container/systemd remediation. macOS supports observed and
cooperative assurance; it cannot claim a Linux cgroup limit.

The runner's prediction includes both the pileup/reference model and the analyzer's
declared maximum retained memory. An unknown analyzer is valid in record-only mode,
but cannot make an enforcement claim and is therefore refused before output under
either enforcement mode.

### 4. Verify — `rosalind verify`

Re-check a receipt *without re-running*:

```bash
rosalind verify --manifest sample.vcf.manifest.json
```

It re-hashes the recorded inputs and outputs (BLAKE3) and re-checks the recorded peak against the budget
(supplied via `--budget-mb`, or read from the manifest). Exit **0** if everything matches and fits;
non-zero (exit **5**) with a per-check report on any drift, missing file, or over-budget peak. This is the
auditability story containers can't give you for a non-deterministic caller.

The receipt has two integrity layers. `manifest_blake3` protects the deterministic claim, including
content hashes, output-affecting parameters, replay recipe, and producer/analyzer/build identity.
`measurement_blake3` protects machine-local resource measurements. In schema 3 and newer, recorded
paths are intentionally portable metadata excluded from the claim hash: relocating a file preserves the
claim, while `verify` still re-hashes the artifact at its recorded path. Claim or measurement edits are
caught (exit **5**); a path-only edit is correctly reported as a non-claim change. This is
tamper-*evident*, not proof of authorship. See the separate
[receipt trust levels](docs/receipts-and-trust.md).

## What's bounded (honest scope)

- **Germline `variants`**, **`features`**, and registered analyzers are bounded:
  peak ≈ the largest selected reference interval + the depth-capped active set,
  independent of BAM size. Whole-genome execution streams the BAM; region/BED/shard
  execution uses indexed interval fetch and decodes only the selected reference span.
- **Somatic** (`somatic`) is **region-bounded**, not whole-genome-bounded (it collects both pileup streams
  for the region).
- **Analysis-reference build** (`reference build`) uses bounded buffering. Legacy
  search-index build (`index`) remains O(reference) RAM and is not required by an analyzer.
- The engine is **single-threaded**. Deterministic reference-span shards provide
  external parallelism and canonical first-party merge.

## Transactional output policy

File-producing commands use sibling temporary files and atomic rename. The default
is create-new: an existing destination or reserved `<output>.partial` causes exit 2
before computation. Pass `--force` to request atomic replacement. Successful runs
leave only the requested final artifact; input, analyzer, and ordinary I/O failures
remove temporary files. Receipts follow the same write discipline. This
non-overwrite default is the intentional v0.4 CLI compatibility change.

## Extend — build on the bounded substrate

Rosalind's kernel is a bounded, deterministic **`PileupColumn` stream**. Compute your own per-locus
analytics over it (coverage, QC, methylation, ML features) and inherit bounded memory + determinism for
free — no variant calling required:

```rust
use rosalind::{PileupEngine, PileupParams};
// PileupEngine<S: ReadSource> is an Iterator<Item = Result<PileupColumn, _>>.
// Each PileupColumn carries the locus, ref_base, depth(), allele_counts(), strand_counts().
for column in PileupEngine::new(source, reference, contig, region, PileupParams::default()) {
    let col = column?;
    // your bounded per-locus metric here
}
```

A complete, runnable example: [`examples/custom_pileup_analytics.rs`](examples/custom_pileup_analytics.rs)
(`cargo run --example custom_pileup_analytics`). For a first-class SDK that inherits the bounded contract,
implement the `ColumnAnalyzer` trait and call the public
`rosalind::contract::run_column_analysis` runner, or
`run_column_analysis_selected` for normalized intervals/shards. It owns planning, refusal, the governor, output,
and receipt sealing without ever terminating the host process. The CLI maps typed refusal and breach
outcomes to exits 3 and 4. Generate a complete standalone example with:

```bash
rosalind new analyzer my-analyzer --output ./my-analyzer
```

The RSS governor is process-wide, so only one enforced runner may be active per process; concurrent
attempts return a typed error. [`ARCHITECTURE.md`](ARCHITECTURE.md) describes the boundary.

## Reproducibility

Every receipt is canonical JSON (sorted keys, no timestamps) with BLAKE3 content hashes of inputs and
outputs, deterministic contract parameters, tokenized replay, producer/analyzer identity, and a separate
measurement block. Identical inputs produce a byte-identical primary artifact and the same portable claim;
machine-local measurements may differ without changing that claim. The current schema and historical
compatibility rules are published in [`docs/receipt-schema.md`](docs/receipt-schema.md).
