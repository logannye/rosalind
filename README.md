# Rosalind

[![CI](https://github.com/logannye/rosalind/actions/workflows/ci.yml/badge.svg)](https://github.com/logannye/rosalind/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/rosalind-bio?logo=rust&label=crates.io&color=orange)](https://crates.io/crates/rosalind-bio)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

> **Build deterministic genomics analyses that fit the machine—and prove exactly what ran.**

Rosalind is a Rust SDK and CLI for turning coordinate-sorted short-read alignments
into reproducible per-locus artifacts under a declared memory budget. It predicts
resource needs before execution, can refuse or govern a run, publishes outputs
transactionally, and records content-addressed evidence for offline verification,
replay, and causal comparison.

Rosalind is for genomics research software engineers, bioinformatics platform
engineers, ML/data engineers building feature pipelines, and workflow engineers
operating fixed-resource HPC or CI jobs. It is not a clinical diagnostic system,
a turnkey interface for wet-lab users, a long-read platform, or a replacement for
GATK, DeepVariant, Clair3, or bcftools. The built-in SNV caller is a reference
workload used to test the platform contract, not the product's scientific moat.

**Research use only. Not for diagnostic use.**

## Why use it?

Most genomics tools let a scheduler impose a memory limit. Rosalind makes the
analysis itself participate in the decision and preserve evidence about it:

- **Resource confidence:** a versioned model predicts peak memory, a declared
  budget can refuse the run before destination creation, and cooperative or
  cgroup-backed enforcement is recorded distinctly from observation-only runs.
- **Reproducible artifacts:** deterministic ordering, create-new transactional
  output, canonical receipts, content hashes, tokenized replay, and causal diff.
- **Inherited infrastructure:** a custom `ColumnAnalyzer` receives the bounded
  pileup walk, planning, refusal/breach semantics, output safety, receipts,
  replay, and conformance checks without rebuilding them.
- **Offline privacy:** receipt and artifact inspection runs locally. Rosalind has
  no telemetry and does not upload genomic data.

This is differentiated from workflow engines, which can request resources but do
not know an analyzer's working set, and from mature callers, which optimize
scientific accuracy and throughput rather than providing a reusable deterministic
analysis contract. Rosalind complements both.

## Five-minute analyzer quickstart

The current source tree is the authoritative development surface. Until the
`0.4.0` stabilization release is published, build it from source:

```sh
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --locked --release
export PATH="$PWD/target/release:$PATH"

rosalind new analyzer locus-qc --output /tmp/locus-qc
cd /tmp/locus-qc
cargo test
```

The generated project is standalone: it implements `ColumnAnalyzer`, calls the
public contract runner, includes contract tests, records its own producer identity,
and can replay through an explicitly selected binary. Start with
[`docs/analyzer-sdk.md`](docs/analyzer-sdk.md).

To exercise the complete built-in contract offline:

```sh
rosalind demo --output-dir /tmp/rosalind-demo
rosalind verify --manifest /tmp/rosalind-demo/calls.vcf.manifest.json
rosalind reproduce \
  --manifest /tmp/rosalind-demo/calls.vcf.manifest.json \
  --inputs /tmp/rosalind-demo
```

File outputs are create-new and atomic by default. `--force` requests atomic
replacement. A governed breach preserves only `<output>.partial` plus evidence;
it never publishes the successful destination name.

## Resource contract

Prepare a deterministic analysis reference and ask whether the run fits:

```sh
rosalind reference build --fasta genome.fa --output genome.rref
rosalind reference inspect --reference-pack genome.rref --json

rosalind plan \
  --reference-pack genome.rref \
  --max-depth 1000 --max-read-len 250 \
  --budget-mb 4096 --json

rosalind analyze coverage \
  --reference-pack genome.rref \
  --alignments sample.sorted.bam \
  --memory-budget-mb 4096 --enforce \
  --output coverage.tsv

rosalind verify --manifest coverage.tsv.manifest.json
```

The `.rref` build uses bounded line buffering and two streaming FASTA passes;
construction memory scales with contig metadata rather than total reference
length. It stores contig metadata, packed sequence, ambiguity masks, normalized
source identity, versioning, and internal integrity checks. Analyzers no longer
need an FM-index.

Legacy `.idx` inputs and receipts remain supported. Analysis commands accept
`--index` with migration guidance, and conversion does not rebuild search data:

```sh
rosalind reference convert --index legacy.idx --output genome.rref
```

The FM-index remains the correct artifact for search:

```sh
rosalind index --reference genome.fa --output genome.idx
rosalind locate --index genome.idx --pattern GATTACA
```

See [`CONTRACT.md`](CONTRACT.md) for the estimator, exit codes, enforcement
assurances, and breach behavior.

## Receipts, replay, and trust

Receipts protect portable claims separately from machine-local measurements.
Modern paths are relocatable metadata; inputs and outputs are matched by content.
Receipts are **tamper-evident, not proof of authorship**: an unsigned receipt can
be resealed by someone who changes its contents.

```sh
rosalind receipt inspect --manifest coverage.tsv.manifest.json --json
rosalind verify --manifest coverage.tsv.manifest.json
rosalind reproduce --manifest coverage.tsv.manifest.json --inputs ./data
rosalind diff --a old.manifest.json --b new.manifest.json --json
rosalind chain verify --dir ./artifacts --json
```

Replay uses validated tokenized arguments, never a shell string. Historical
receipt schemas 1–5 remain parseable at their original capability level. The
loopback-only Receipt Studio and browser verifier inspect artifacts locally.
See [`docs/receipts-and-trust.md`](docs/receipts-and-trust.md).

## What is supported today

This table describes the current development tree. Region/shard/merge, Arrow,
Python wheels, and workflow integrations are unreleased v0.5 work until their
release and evidence gates pass; no published 0.4 artifact is implied.

| Surface | Current contract |
|---|---|
| Input | Coordinate-sorted short-read BAM; contig names and lengths must match the reference |
| Analysis reference | Preferred deterministic `.rref`; compatible legacy `.idx` |
| SDK | Streaming, single-process `ColumnAnalyzer` over depth-capped pileup columns |
| Built-in analyzers | Per-locus feature TSV/Arrow IPC and coverage track |
| Selection | Whole genome without an index; indexed region/BED and deterministic reference-span shards with BAM+BAI |
| Parallelism | External shard fan-out plus receipt-validated canonical first-party merge; no internal threads |
| Resource modes | Observation-only, cooperative enforcement, optional existing Linux cgroup-v2 limit |
| Outputs | Deterministic bytes within a Rosalind release; atomic create-new/replace |
| Evidence | Canonical schema-5 receipt, verify, replay, reproduction certificate, diff, chain |
| Caller | Diploid, single-sample, short-read SNVs; simulated evidence only; not a production caller |
| Search | Separate `.idx` FM-index and exact `locate` |

Current limitations are explicit:

- Execution is single-threaded; deterministic shards are the parallelism model.
  Canonical merge is intentionally limited to first-party TSV, Arrow IPC, sites
  VCF, and gVCF codecs.
- Region/BED/shard execution requires a `.bai`; CSI and CRAM/CRAI are deferred.
- Linux ARM, long reads, indels, configurable
  ploidy, and production-caller validation are not currently supported.
- A real HG002 GIAB baseline has not yet been established. The pinned workflow
  exists, and its first result will be published even if poor.
- Resource predictability is not a speed or accuracy claim.

## Built-in workloads

The two maintained analyzer examples share one bounded kernel:

```sh
rosalind analyze coverage \
  --reference-pack genome.rref --alignments sample.sorted.bam \
  --output coverage.tsv

rosalind features \
  --reference-pack genome.rref --alignments sample.sorted.bam \
  --format arrow-ipc --output features.arrow
```

Sparse selections use samtools-style 1-based inclusive regions or BED's zero-based
half-open coordinates. Shards own a complete, non-overlapping span of the ordered
reference, independent of read content:

```sh
rosalind features --reference-pack genome.rref --alignments sample.sorted.bam \
  --format arrow-ipc --shard-count 4 --shard-index 0 --output shard-0.arrow

rosalind merge \
  --manifest shard-0.arrow.manifest.json \
  --manifest shard-1.arrow.manifest.json \
  --manifest shard-2.arrow.manifest.json \
  --manifest shard-3.arrow.manifest.json \
  --inputs . --output features.arrow
```

The merge refuses missing, duplicate, overlapping, tampered, breached, or
contract-incompatible parents. A complete Arrow shard merge is rebatched into
canonical 65,536-row batches and is byte-identical to the unsharded stream.

Variant calling remains available as a reference workload. The CLI reporting
threshold defaults to 30 so it matches the published simulated configuration:

```sh
rosalind variants \
  --reference-pack genome.rref --alignments sample.sorted.bam \
  --quality-threshold 30 --output calls.vcf
```

Use `--quality-threshold 10` explicitly for the historical permissive behavior.
Do not infer clinical validity or competitiveness from the bundled simulated test.

## Evidence and benchmarks

The claims harness re-derives contract properties on bundled data:

```sh
bash benchmarks/run.sh
```

It checks conservative prediction, honor-or-refuse behavior, repeat-run byte
identity, replay/tamper detection, and deterministic packing. Measured findings
live under [`docs/findings/`](docs/findings/); the opt-in GIAB workflow lives under
[`benchmarks/giab/`](benchmarks/giab/). Claims without an executable evidence path
are intentionally excluded from this README.

The platform harness compares equivalent Rosalind, pysam, and bcftools extraction
on the same input and derives its report from retained raw measurements. Published
HG002 results remain pending; no result is inferred from the existence of the
harness. See [`benchmarks/platform/`](benchmarks/platform/) and
[`docs/benchmarks-and-limitations.md`](docs/benchmarks-and-limitations.md).

The unreleased Python package is a Maturin mixed project: distribution
`rosalind-bio`, import `rosalind`. `iter_features()` yields native PyArrow record
batches without accumulating the genome; `collect_features()` is the explicit
materialization path. Platform wheels and publishing remain release-gated.

## Rust API

The reference abstraction excludes search operations by design:

```rust,no_run
use rosalind::{ReferencePackReader, ReferenceProvider};

fn inspect(path: &std::path::Path) -> anyhow::Result<()> {
    let reference = ReferencePackReader::open(path)?;
    println!("{} bases", reference.len());
    for contig in reference.contigs().iter() {
        println!("{}\t{}", contig.name, contig.length);
    }
    Ok(())
}
```

The SDK front door re-exports `ColumnAnalyzer`, `PileupColumn`,
`ContractRunSpec`, `run_column_analysis`, typed refusal/breach outcomes, and the
contract testkit. API details are in [`docs/analyzer-sdk.md`](docs/analyzer-sdk.md)
and the generated rustdoc.

## Roadmap and contribution

The delivery sequence, evidence gates, and non-goals are maintained in
[`docs/ROADMAP.md`](docs/ROADMAP.md). The current architecture is
[`ARCHITECTURE.md`](ARCHITECTURE.md), and the exact implemented/deferred split is
[`docs/implementation-status.md`](docs/implementation-status.md). Historical implementation exploration under
`docs/superpowers/`, `docs/OPEN_PROBLEMS.md`, and `docs/GROWTH.md` is non-primary archive material; external-memory search-index work
does not return to the product roadmap without the documented user-evidence gate.

Contributions should preserve determinism, bounded analysis state, receipt
compatibility, and transactional output behavior. See
[`CONTRIBUTING.md`](CONTRIBUTING.md), [`SECURITY.md`](SECURITY.md), and
[`CITATION.cff`](CITATION.cff).

## License

Dual-licensed under Apache-2.0 and MIT.
