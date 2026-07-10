# The memory contract

**Memory is a contract, not a hope.**

Every other variant caller treats RAM as an emergent property you guess at (`-Xmx…`,
`--target-mem` "heuristics may not work well") and then crash on. Rosalind treats it as a contract:

> **Rosalind never silently OOM-kills you — it fits, or it tells you up front, and it records the
> realized peak in a receipt you can verify.**

(It does *not* claim to "never refuse": when a budget is genuinely too small the run declines cleanly
rather than crashing. Graceful degrade-don't-die — sliding down a space/time curve to finish anyway — is
the Phase-D research direction, not a present claim.)

## The four verbs

The contract applies to the bounded whole-genome paths — `rosalind variants --index` (germline calling) and `rosalind features --index` (per-locus feature egress), which share the same streaming engine and therefore the same `plan`/`--enforce`/`verify` envelope.

### 1. Declare

State the RAM you have. `--memory-budget-mb N` on `variants`, `--budget-mb N` on `plan`/`verify`.

### 2. Predict — `rosalind plan`

Ask *before committing a byte* whether the job fits. `plan` reads only the index header (plus your
declared depth/read-length assumptions) — it never opens the BAM:

```bash
rosalind plan --index genome.idx --max-depth 1000 --max-read-len 250 --budget-mb 2048
```

It prints a breakdown — reference decode + active set @ max-depth + engine overhead, atop a measured
process baseline — and a verdict: `[FITS]` or `[REFUSE]`.

### 3. Honor — `rosalind variants … --enforce`

```bash
rosalind variants --index genome.idx --alignments sample.sorted.bam \
  --memory-budget-mb 2048 --enforce -o sample.vcf
```

With `--enforce`:

- **predicted peak > budget → refuse up front** (exit **3**), before doing any work, with an actionable
  message (raise the budget, lower `--max-depth`, or drop `--enforce`);
- **realized peak > budget → fail loud** (exit **4**) *after* writing the VCF + receipt (you keep the data
  and the proof it overran) — never a silent overrun;
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
[receipt trust levels](docs/receipt-trust.md).

## What's bounded (honest scope)

- **Germline `variants --index`** and **`features --index`** are the bounded paths: peak ≈ the largest
  contig's reference + the depth-capped active set, **independent of BAM size**. Reads stream one record at
  a time; output rows (VCF calls or feature rows) stream straight to disk with no genome-wide buffer.
- **Somatic** (`somatic`) is **region-bounded**, not whole-genome-bounded (it collects both pileup streams
  for the region).
- **Index *build*** (`rosalind index`) is **O(reference)** in RAM today; `plan --reference` reports an
  advisory estimate. Sublinear-space construction is the Phase-D research direction (see
  `docs/OPEN_PROBLEMS.md`).
- The engine is **single-threaded** — outputs are deterministic, but there is no thread-invariance claim
  yet.

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
`rosalind::contract::run_column_analysis` runner. It owns planning, refusal, the governor, output,
and receipt sealing without ever terminating the host process. The CLI maps typed refusal and breach
outcomes to exits 3 and 4. Generate a complete standalone example with:

```bash
rosalind new analyzer my-analyzer --output ./my-analyzer
```

The RSS governor is process-wide, so only one enforced runner may be active per process; concurrent
attempts return a typed error. [`docs/architecture.md`](docs/architecture.md) describes the boundary.

## Reproducibility

Every receipt is canonical JSON (sorted keys, no timestamps) with BLAKE3 content hashes of inputs and
outputs, deterministic contract parameters, tokenized replay, producer/analyzer identity, and a separate
measurement block. Identical inputs produce a byte-identical primary artifact and the same portable claim;
machine-local measurements may differ without changing that claim. The current schema and historical
compatibility rules are published in [`docs/receipt-schema.md`](docs/receipt-schema.md).
