# Changelog

All notable changes to Rosalind are recorded here. Versions follow [Semantic Versioning](https://semver.org).

## [0.1.0] — 2026-06-02

First tagged release: a deterministic, low-memory genomics engine where **memory is a verifiable
contract** — predict it before you commit, honor it during the run, and verify it after.

### The memory contract
- `rosalind plan` — predict a job's peak memory against a declared budget *before* committing a byte.
- `rosalind variants … --enforce` — honor the budget: refuse up front (exit 3) or fail loud (exit 4),
  never a silent OOM-kill. Record-only without `--enforce`.
- `rosalind verify` — re-check a run's BLAKE3 receipt without re-running.
- The **`rosalind-budget` GitHub Action** — enforce the contract in *your* CI (fail the build on breach).

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

[0.1.0]: https://github.com/logannye/rosalind/releases/tag/v0.1.0
