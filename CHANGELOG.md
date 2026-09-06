# Changelog

All notable changes to Rosalind are recorded here. Versions follow Semantic Versioning.

## [Unreleased]

### Evidence-engine candidate track (0.5.0)

The feature-bearing evidence engine is being prepared as `0.5.0-rc.1`. No 0.4
release was published. The existing seven-day RC soak, independent-user gates,
protected publication review, and change-based caller evidence requirements apply.

- Indexed exact SNV evidence and panel QC, bounded Rust/Python batches, verified
  partition cache/resume, canonical workers, and dataset comparison are implemented.
- Panel callability now defaults to 10x consistently through Rust, CLI, and Python.
- Named, unknown, and explicitly pooled sample scope is recorded in evidence
  receipts and cache identities. Ambiguous sample headers require explicit intent;
  named runs refuse reads that cannot be assigned to a declared sample.
- Development status and the ordered next implementation sequence are tracked in
  `docs/ADOPTION_ROADMAP.md`; implementation does not imply publication or adoption.

### Analyzer-platform foundation

- Repositioned Rosalind around deterministic per-locus analysis contracts; the
  built-in caller is documented as a reference workload rather than a competitor
  to production callers.
- Added the `.rref` analysis-reference format, `ReferenceProvider`, streaming
  bounded-buffer construction, integrity validation, and legacy `.idx` conversion.
- Added `rosalind reference build|inspect|convert` and `--reference-pack` support
  for analysis, preflight, and planning while preserving legacy `.idx` replay.
- Added normalized region/BED selection, indexed BAM fetch, deterministic
  `reference-span-v1` shards, schema-5 partition claims, and receipt-validated
  canonical merge for first-party TSV, Arrow, sites VCF, and gVCF outputs.
- Added canonical Arrow IPC feature streams with a fixed schema and 65,536-row
  batches, plus a mixed Maturin/PyArrow package with lazy batch iteration.
- Added digest-pinned OCI build automation, generalized Action inputs, versioned
  Nextflow/Snakemake modules and local/Slurm examples, and an executable
  Rosalind/pysam/bcftools platform evidence harness.
- Changed the CLI variant quality threshold default to 30; the historical
  permissive behavior remains available with `--quality-threshold 10`.
- Version-gated the design-partner and seven-day RC soak requirements from 0.5.0;
  the subsequent feature-bearing candidate retains those gates.

## [0.4.0] — unpublished development checkpoint

This was the planned consolidated stabilization checkpoint after `0.1.0`, now
included in the feature-bearing 0.5 candidate track. The repository's
former `0.1.1`–`0.3.1` labels were internal development checkpoints, not public
GitHub or crates.io releases; their changes are consolidated here.

### Analyzer SDK and contract runner

- Public `ColumnAnalyzer` and `run_column_analysis` surfaces with declared analyzer
  memory models, typed refusal/breach outcomes, process-wide cooperative governance,
  and distinct observation-only/cooperative/cgroup-v2 assurance.
- Standalone analyzer scaffolding, conformance checks, contract testkit, and
  explicit external-binary replay.
- Maintained feature and coverage analyzers over the bounded whole-genome pileup.

### Transactional and reproducible execution

- Create-new atomic outputs by default, explicit atomic replacement with `--force`,
  and `<output>.partial` preservation on governed breach.
- Canonical schema-5 receipts with portable content claims, separate protected
  machine measurements, build identity, tokenized replay, offline verification,
  reproduction certificates, causal diff, and provenance-chain verification.
- Historical schema 1–5 parsing and capability-aware verification.
- Loopback-only Receipt Studio and client-side WASM verification with no genomic
  data upload or telemetry.

### Genomics workloads and evidence

- Bounded, coordinate-sorted BAM processing; deterministic feature TSV, coverage,
  germline sites/gVCF reference workloads, and external-merge BAM sort.
- Shared base-quality and MAPQ germline likelihood, genotype-aware evaluation, and
  a checksum-pinned opt-in HG002 GIAB workflow. No real GIAB score is claimed yet.
- A claims harness covering prediction, honor-or-refuse behavior, byte identity,
  replay/tamper evidence, and deterministic fleet packing.

### Release engineering

- Maintainer-only authenticated release plans, contract fingerprints, protected
  RC/stable workflows, resumable byte-verified crates.io publication, package
  checksums, and post-publish installation tests.

## [0.1.0] — 2026-06-02

Initial public release of the deterministic bounded pileup/calling engine, persisted
FM-index, memory plan/enforce/verify loop, reproducibility receipt, feature TSV,
binary installer, and budget GitHub Action.

[Unreleased]: https://github.com/logannye/rosalind/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/logannye/rosalind/releases/tag/v0.1.0
