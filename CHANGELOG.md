# Changelog

All notable changes to Rosalind are recorded here. Versions follow [Semantic Versioning](https://semver.org).

## [Unreleased]

### v0.4 trustworthy builder train

- Sound public analyzer contracts now distinguish unknown and fixed memory models,
  cooperative enforcement, observed-only evidence, and Linux cgroup-v2 assurance.
- File outputs and receipts are create-new and atomic by default; `--force` requests
  atomic replacement, while governed breaches preserve only `<output>.partial`.
- Replay schema 3 validates explicit execution plans before running, allowlists
  built-ins, requires `--binary` for external analyzers, isolates the working
  directory, and exposes non-executing `--dry-run` plans.
- One shared structured trust report powers `verify`, Receipt Studio, WASM, and
  badges. Missing reproduction evidence is neutral; only linked certificates earn
  `reproduced`.
- Added `doctor`, `receipt inspect`, `receipt sanitize`, unsigned in-toto export,
  local-only `studio`, complete demo provenance, analyzer conformance, and embedded
  maintained scaffold templates.
- `align` and `sort` now emit replayable schema-5 receipts with artifact roles. The
  opt-in HG002 workflow adds checksum-pinned hap.py 0.3.15/RTG vcfeval evaluation
  and GIAB v3.1 genome contexts without changing the caller.
- Added a maintainer-only `cargo xtask` release control plane with authenticated
  plans, contract fingerprints, protected RC/stable workflows, resumable
  byte-verified crates.io publication, pinned GHCR evaluator images, attested GIAB
  baseline PRs, and anonymized design-partner release gates.

Release candidates must soak for one week. After the v0.4.0 RC begins, public API,
receipt-field, or CLI changes require restarting the soak.

## [0.3.1] — 2026-07-09

### Scientific credibility and release hardening

- Germline sites and gVCF now share the `baseq-mapq-v1` likelihood: mapping
  uncertainty is marginalized into every observation; MAPQ 255 retains the historical
  base-quality-only model. Somatic calling is unchanged.
- `eval-germline --calls-filter all|pass --json` publishes both emitted and
  PASS-only metrics, including genotype concordance and call counts.
- The opt-in HG002 GIAB v5.0q GRCh38 chr20 workflow pins exact source URLs and
  SHA-256 hashes, verifies downloads, writes a local data manifest, and guards the
  first honest baseline. Pull requests run only a synthetic chr20-shaped smoke test.
- Packaged assets, workspace tests, clippy, MSRV, and release packaging remain gates.

## [0.3.0] — 2026-07-09

### Receipt Studio

- The existing `/verify/` page is now a framework-free, client-only studio with
  multi-file drag/drop, streaming BLAKE3 artifact hashing, content-based matching,
  causal receipt diff, provenance-chain traversal, mobile layout, keyboard access,
  and no third-party network requests.
- Trust is no longer collapsed into “verified”: receipt integrity, artifacts,
  resource contract, reproduction evidence, independent attestation, and signatures
  are distinct levels. Schema-3+ paths are visibly relocatable metadata.

## [0.2.1] — 2026-07-09

### Fork-builder onboarding

- `rosalind new analyzer NAME --output DIR` scaffolds a standalone analyzer binary,
  build identity, contract tests, CI, and a full external replay check without
  overwriting non-empty destinations.
- The `contract-testkit` feature exposes repeat-run, receipt-integrity, path-
  relocation, completion, and refusal assertions.
- `rosalind demo` is an embedded, offline index → align → sort → plan → enforce →
  verify → reproduce → chain walkthrough with human and JSON output.

## [0.2.0] — 2026-07-09

### Public contract runner and external replay

- `rosalind::contract::run_column_analysis` exposes typed configuration, refusal,
  breach, and successful outcomes. It owns output lifetime, the process-wide governor,
  and receipt sealing and never exits the host process. `features` and `analyze`
  delegate to it.
- Receipts bind producer and analyzer identity and add replay schema 2 with canonical
  `command_argv`; the legacy display command remains readable and replayable.
- `reproduce --binary PATH --json` explicitly selects third-party code, never trusts
  a recorded executable path, and preserves both original and rerun build identities.
- Added the publishable `rosalind-build-info` build dependency.
- Bumped the leaf `rosalind-receipt` crate to 0.2.0 for the argv replay and dual-identity API.

## [0.1.1] — 2026-07-09

### Trust and compatibility baseline

- Corrected receipt promises: claim fields and measurements are protected, while
  schema-3+ paths are portable metadata excluded from the claim hash.
- Published the schema-5 JSON envelope, extension namespaces, trust-level model, and
  canonical schema 1–5 fixtures. Historical receipts remain parseable and verify at
  their original capability level.
- `verify --json` and `reproduce --json` expose structured reports.
- Preserved the existing memory governor, replay, badge, causal diff, chain,
  deterministic gVCF, fleet packer, and build-identity foundations.

## [0.1.0] — 2026-06-02

First tagged release: a deterministic, low-memory genomics engine where **memory is a verifiable
contract** — predict it before you commit, honor it during the run, and verify it after.

### The memory contract
- `rosalind plan` — predict a job's peak memory against a declared budget *before* committing a byte.
- `rosalind variants … --enforce` — honor the budget: refuse up front (exit 3) or fail loud (exit 4),
  with an explicit cooperative verdict. Record-only without `--enforce`.
- `rosalind verify` — re-check a run's BLAKE3 receipt without re-running.
- The **Rosalind budget GitHub Action** (`action.yml`, used as `logannye/rosalind@v0.1.0`) — enforce the
  contract in *your* CI (fail the build on breach).

### Variant calling
- **Bounded whole-genome germline SNV calling** (`variants --index`) over a coordinate-sorted BAM and a
  persisted index: peak memory tracks coverage, not BAM size. Genotype-likelihood, abstention-aware.
- **Unbiased depth-cap downsampling** — deterministic content-hash selection that bounds the working set without
  biasing allele balance (no silent variant drops); dropped-read counts surfaced in the receipt.
- **Tumor/normal somatic** SNV calling (`somatic`).
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

[Unreleased]: https://github.com/logannye/rosalind/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/logannye/rosalind/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/logannye/rosalind/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/logannye/rosalind/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/logannye/rosalind/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/logannye/rosalind/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/logannye/rosalind/releases/tag/v0.1.0
