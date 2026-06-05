# Changelog

All notable changes to Rosalind are recorded here. Versions follow [Semantic Versioning](https://semver.org).

## [Unreleased]

A contract-hardening pass followed by moat-compounding capabilities, all on the
canonical bounded-memory substrate.

### Reproducibility — a stranger can re-derive your result (Track D)
- **`rosalind reproduce`** re-derives a recorded result byte-for-byte from its receipt and
  content-located inputs, reporting REPRODUCED / DIVERGED / INCONCLUSIVE (exit 0 / 6 / 7; a
  tampered receipt is exit 5). No incumbent caller can offer this — a non-deterministic caller
  reports DIVERGED on a *correct* run. Honest scope: the deterministic text outputs
  (`variants → VCF`, `features → TSV`); BAM/bgzf is reported INCONCLUSIVE, never a false DIVERGED.
- **Chainable reproduction certificate** (`<receipt>.repro.json`): a content-addressed,
  self-hashing attestation that names the original receipt's claim hash. N certificates over the
  same parent are N independent confirmations — a serverless reproducibility web. Signing-ready,
  and itself `verify`-able.
- **Schema 5 — a replayable command.** Every receipt now records a normalized, machine-independent
  `command` recipe (one `CommandCapture` chokepoint), closing the prior gap where
  `gvcf`/`chrom`/index-vs-reference mode went unrecorded. Pre-v5 receipts still parse and verify.
- **`rosalind badge`** emits a self-hosted shields.io endpoint JSON + a static SVG
  ("reproducible · fits N MiB") with no shields.io runtime dependency (works offline).
- **CI reproduce fence:** the `cli-e2e` job re-derives a result on the GitHub runner (a different
  machine than the author's) and confirms a byte-changed input is reported INCONCLUSIVE.
- Internals: `verify_receipt` is now a shared library function — `verify`, `reproduce`, and a
  future WASM verifier consume one source of truth, so they cannot drift.

### Contract hardening (true & trusted on real genomes)
- **Predicted peak is a true upper bound.** The per-contig reference decode no longer holds a
  transient second copy (`decode_window_arc`); the prediction is recorded in the receipt and carries an
  honest I/O margin. The contract's core inequality (predicted peak ≥ realized peak RSS) is now tested.
- **Real-genome correctness.** `eval-germline`/`eval-somatic` normalize each variant against its own
  contig (multi-contig benchmarks no longer miscompare or crash); a contig-naming mismatch (UCSC `chr1`
  vs Ensembl `1`) is refused up front instead of silently writing an empty VCF; `index` maps IUPAC
  ambiguity codes to `N` (stock references ingest).
- **Trust on-ramp.** The GitHub Action snippet resolves (`logannye/rosalind@v0.1.0`); `install.sh`
  verifies the `.sha256` it advertises; the CI fixture check is a real pinned hash, not a tautology;
  `verify` cross-checks the receipt's internal consistency; `Cargo.lock` is tracked.

### Fleet scheduling — prediction → placement
- **`rosalind pack`** packs many bounded `variants` jobs onto fixed-size nodes by their predicted peaks
  (read from each index header, additive) and *proves* a co-location fits before launching a byte, or
  refuses (exit 3). `plan --index --json` emits the predicted peak for a scheduler to read.

### ColumnKit SDK — implement one trait, inherit the contract
- **`ColumnAnalyzer` trait + `run_bounded_whole_genome` driver.** A builder's own per-locus analyzer
  inherits the bounded whole-genome walk, the working-set bound, and the verifiable receipt. The shipped
  `features` egress is the first impl (the SDK is the production path, not a parallel one).

### Cohort-ready gVCF
- **`variants --index --gvcf`** emits a banded gVCF (every callable locus → a variant or a `<NON_REF>`
  reference block with an `END=` span), so per-sample output joins into GLnexus/GATK — bounded (O(1)
  banding state) and byte-reproducible.

## [0.1.0] — 2026-06-02

First tagged release: a deterministic, low-memory genomics engine where **memory is a verifiable
contract** — predict it before you commit, honor it during the run, and verify it after.

### The memory contract
- `rosalind plan` — predict a job's peak memory against a declared budget *before* committing a byte.
- `rosalind variants … --enforce` — honor the budget: refuse up front (exit 3) or fail loud (exit 4),
  never a silent OOM-kill. Record-only without `--enforce`.
- `rosalind verify` — re-check a run's BLAKE3 receipt without re-running.
- The **Rosalind budget GitHub Action** (`action.yml`, used as `logannye/rosalind@v0.1.0`) — enforce the
  contract in *your* CI (fail the build on breach).

### Variant calling
- **Bounded whole-genome germline SNV calling** (`variants --index`) over a coordinate-sorted BAM and a
  persisted index: peak memory tracks coverage, not BAM size. Calibrated, abstention-aware.
- **Unbiased depth-cap downsampling** — a content-hash reservoir that bounds the working set without
  biasing allele balance (no silent variant drops); dropped-read counts surfaced in the receipt.
- **Tumor/normal somatic** SNV + simple-indel calling (`somatic`).
- **Measured detection accuracy** on simulated diploid truth (precision/recall 1.00/1.00 on clean data);
  `eval-germline` / `eval-somatic` truth-set comparison (the `eval-germline` path is GIAB-ready).

### Index, alignment, I/O
- Build-once, memory-mapped, byte-reproducible FM-index (`index` / `locate`).
- Single-contig FM-index aligner (`align`); deterministic external-merge coordinate sort (`sort`);
  streaming gzip/bgzf FASTA/FASTQ input.

### Reproducible ML feature substrate
- `rosalind features --index` — the same bounded stream as a per-locus feature **TSV**, byte-identical
  run-to-run, with a hash receipt: **bit-reproducible training inputs**.
- `python/rosalind.py` (stdlib + numpy) loads it; `examples/reproducible_features_demo.py` proves the
  bit-reproducibility end-to-end.

### Reproducibility & distribution
- Byte-identical primary outputs and a canonical-JSON BLAKE3 receipt per run.
- Prebuilt static/native binaries via a tagged release + `install.sh`; a 60-second quickstart.

### Research direction
- The `~√t` (square-root-space) framework is retained as honest framing for **sublinear-space index
  construction** (Phase D); today's index build is `O(reference)` and the contract covers call/query.
  See [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).

[Unreleased]: https://github.com/logannye/rosalind/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/logannye/rosalind/releases/tag/v0.1.0
