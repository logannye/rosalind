# The memory contract

**Memory is a contract, not a hope.**

Every other variant caller treats RAM as an emergent property you guess at (`-Xmx…`,
`--target-mem` "heuristics may not work well") and then crash on. Rosalind treats it as a contract:

> **Rosalind never silently OOM-kills you — it fits, or it tells you up front, and it proves the
> realized peak with a receipt.**

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
- otherwise the run completes within budget.

Without `--enforce`, the budget is **record-only**: the run always completes and the verdict is recorded in
the receipt. The active read set is capped at `--max-depth` (default 1000; `0` = uncapped) by an
**unbiased** content-hash reservoir — it bounds the working set without biasing allele balance, so a deep
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
(`cargo run --example custom_pileup_analytics`).

> **Legacy / non-bounded.** The `GenomicPlugin` trait (`src/plugin/`), the `framework/` evaluator, and the
> Python `run_rna_seq_plugin` demo still work but do **not** inherit the memory contract. Prefer the
> `PileupColumn` substrate above for bounded work.

## Reproducibility

Every receipt is canonical JSON (sorted keys, no timestamps) with BLAKE3 content hashes of the index, the
alignments, and the output VCF, plus the realized `peak_rss_bytes` / `max_working_set_bytes` and the
contract params (`memory_budget_mb`, `contract_verdict`, `enforced`, `max_depth`, `max_read_len`). Identical
inputs produce a byte-identical VCF, and a manifest identical except for the realized `peak_rss_bytes` (a
machine-dependent measurement) — `rosalind verify` re-checks the recorded hashes and peak against the budget.
