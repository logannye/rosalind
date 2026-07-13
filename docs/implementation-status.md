# Overhaul implementation status

Updated 2026-07-12. This file separates repository work from evidence and adoption
milestones that cannot be fabricated by implementation alone.

## Implemented on `codex/analyzer-platform-overhaul`

- Analyzer-platform positioning, intended-user/non-user boundaries, research-use
  warning, scaffold-first onboarding, current limitations, and claim discipline.
- Consolidated architecture, changelog, roadmap, SDK/reference/trust/benchmark docs,
  security policy, and citation metadata.
- CLI quality threshold default 30 with explicit historical threshold 10.
- 0.4 stabilization policy: partner and seven-day soak gates begin at 0.5.
- `.rref` v1, bounded-line/two-pass FASTA builder, deterministic packing, checksums,
  mmap reader, `ReferenceSequence`, `ReferenceProvider`, legacy `.idx` provider and
  converter, and atomic publication.
- `reference build|inspect|convert`; `--reference-pack` on features, analyze,
  variants, doctor, and plan; guided `--index` compatibility and replay flag/role.
- Public normalized selection types; samtools-style region and BED parsing;
  deterministic `reference-span-v1` shards; BAI discovery and indexed interval
  fetch; selected execution in features, analyze, variants, doctor, and plan.
- Schema-5 partition claims and BED input hashing, with whole-genome source
  compatibility and legacy `.idx` replay preserved.
- Receipt-driven canonical merge with stable exit codes and first-party codecs for
  feature/coverage TSV, Arrow IPC, sites VCF, and boundary-coalesced gVCF.
- Native uncompressed Arrow IPC v5 streams, a fixed 19-field schema, deterministic
  65,536-row batches, empty streams, and canonical shard rebatching.
- Mixed Maturin package (`rosalind-bio` / import `rosalind`) with the matching
  binary, lazy PyArrow batches, typed run/receipt objects, explicit collection,
  version enforcement, and early-child termination.
- macOS arm64/x86_64 and manylinux x86_64 wheel workflow with RC TestPyPI and
  stable PyPI trusted publishing through the protected release environment.
- Reference build/convert receipts; deterministic `.rref` and `.arrow` replay.
- Digest-pinned non-root OCI definition and provenance workflow; versioned
  Nextflow processes, a plan-derived memory request, Snakemake local/Slurm rules,
  deterministic fan-out/fan-in examples, and a generalized GitHub Action.
- GIAB harness converted to `.rref` with reference receipt verification, actual
  CLI defaults, exact argv, difficult-region evaluator output, and retained evidence.
- Pinned platform benchmark comparing Rosalind/pysam/bcftools from raw time,
  environment, semantic-diff, repeat-hash, memory-probe, and Arrow-merge evidence.
- macOS reference-pack smoke CI, RustSec audit job, scheduled receipt/replay parser
  fuzzing, and pinned local actionlint/shellcheck validation.
- Protected GitHub `rc` and `release` environments with required owner review.
- Claims-harness issue closed; crates.io issue updated with the remaining blocker.

## Verified

- Complete workspace test suite passes, including schema 1–5 fixtures, external
  analyzer conformance, selected analysis, canonical merge, and release policy.
- `.rref` deterministic rebuild, ambiguity round-trip, corruption/truncation
  rejection, Unicode/space-safe paths, legacy identity preservation, and feature
  byte equivalence against `.idx`.
- `cargo clippy --workspace --all-targets -- -D warnings`.
- Full workspace/all-target check passes under the declared Rust 1.83.0 MSRV with
  Arrow Rust crates locked to 55.2.0.
- Python syntax compilation; a native macOS arm64 `py3-none` wheel built and
  installed from scratch; its bundled binary version matched; lazy PyArrow
  iteration produced 3,972 toy rows and an inspectable receipt.
- Property-based randomized shard ownership; BED comments/unusual names/order,
  zero-length/overflow/unknown-contig cases; missing-BAI pre-publication refusal;
  empty TSV/Arrow/VCF artifacts; missing/tampered merge refusal; feature TSV,
  coverage TSV, sites VCF, gVCF, and Arrow shard-merge byte identity; Arrow schema
  and 0/1/65,535/65,536/65,537/131,072 batch-boundary tests.
- Claims harness: 7/7.
- Snakemake 9.8.1 lint/parser and dry-run succeeded; a real two-shard local toy
  workflow planned 87 MiB per shard, ran both governed Arrow jobs, merged them,
  and verified the merge receipt.
- Nextflow 25.04.8 with Java 21 parsed the versioned config and completed a real
  local doctor, two plans, two governed Arrow shards, canonical merge, and receipt
  verification with plan-derived scheduler memory.
- `rosalind-build-info` and `rosalind-receipt` package and verify from their
  archives; the root mixed Python/Rust wheel builds and installs successfully.
- actionlint 1.7.12 and shellcheck 0.10.0 pass.

## External blockers and evidence gates

`CARGO_REGISTRY_TOKEN` is not configured. No crates, tags, assets, release, or
fresh public install can be published until the owner supplies that secret. The
release doctor reports no other prerequisite blocker. A direct root `cargo
package` correctly stops until `rosalind-receipt 0.3.0` is published in the
declared dependency order; crates.io currently exposes only 0.1.0.

The local Docker daemon returned HTTP 500 on its socket, so the pinned platform
container could not be built or executed in this workspace. The harness passed
shell/Python/static validation, but no benchmark number is claimed. The GIAB data
is opt-in and its attested hap.py image lock still has no published digest, so the
roughly 2 GiB scientific run and human baseline review remain external work.

## Not yet implemented or evidenced

- Full coverage/QC metric set beyond the maintained coverage track.
- Actual HG002 GIAB execution, human scientific/provenance review, and merged
  baseline candidate.
- Successful pinned-container platform execution on toy and HG002 inputs, a clean
  second-machine reproduction, and a published report.
- CI execution of macOS x86_64/manylinux wheels, TestPyPI/PyPI publication,
  container publication/lock update, nf-test/container execution, and Snakemake
  container/Slurm execution.
- Published small microbial/human reference packs and their immutable source locks.
- Three-persona 0.5 RC validation, external repeat adoption, prevented-OOM evidence,
  independent reproduction, or an external analyzer repository.
- Every conditional Wave 4 item; none has met its user-evidence trigger.

These remaining items are not described as stable features in the README.
