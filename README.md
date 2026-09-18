# Rosalind

[![CI](https://github.com/logannye/rosalind/actions/workflows/ci.yml/badge.svg)](https://github.com/logannye/rosalind/actions/workflows/ci.yml)
[![Latest published release](https://img.shields.io/github/v/release/logannye/rosalind)](https://github.com/logannye/rosalind/releases/latest)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

**Reusable genomic evidence for research and analysis tools.**

Rosalind is a CLI, Rust library, and Python interface for researchers inspecting
candidate SNVs and target regions, and builders developing analyses on that evidence.
Supply indexed BAM/CRAM plus candidate sites or BED intervals to obtain exact
read counts, quality summaries, and panel QC.

Its core idea: **memory settings control how work is scheduled, without changing
successful scientific results.** Smaller tiles can trade memory for extra indexed
I/O. Counts remain exact; insufficient resources cause refusal or an explicit failure.

**Development preview:** these features are in the 0.5.0 source tree. The latest
public stable release is 0.1.0 and does not provide this evidence workflow.
Use the [source installation guide](docs/installation.md) for the features below;
a source version or candidate build does not establish package publication. See
[implementation status](docs/implementation-status.md) for verified capabilities
and [release gates](docs/ROADMAP.md) for what remains.

## What you can do

| Task | Rosalind provides |
|---|---|
| Inspect candidate SNVs | Allele/strand counts, base and mapping quality summaries, read-position and filtering counts |
| Check a BED panel | Per-target coverage and configurable callability, with optional position evidence from the same traversal |
| Annotate variants | Evidence added to VCF, VCF.gz, or BCF while preserving records, genotypes, and allele order |
| Reuse an analysis | Verified Arrow datasets, cache/resume, subset queries, and Parquet export for Python, R, or SQL |
| Build a new analyzer | A Rust batch SDK with managed resource admission, cancellation, output publication, receipts, and replay |

## Get started

Choose a starting point:

- **[Analyze my data](docs/researcher-quickstart.md)**: inspect supplied SNVs, interpret evidence, verify, and replay.
- **[Build an analyzer](docs/builder-quickstart.md)**: implement a bounded analyzer over live or saved evidence.

Then [reuse saved evidence](docs/reuse-quickstart.md) to answer a second question,
or browse [all documentation](docs/index.md) for concepts, workflows, and troubleshooting.

For the current source preview, build with Rust 1.83 or newer and the
[native build prerequisites](docs/analyzer-sdk.md#build-prerequisites):

```sh
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --locked --bin rosalind
export PATH="$PWD/target/debug:$PATH"
```

**No data ready?** The [small real-data tutorial](examples/research-filter/README.md)
downloads 118 KB of pinned public data and walks through extraction, a candidate
filter, verification, and replay.

With your own indexed alignments, indexed FASTA, and SNV candidates:

```sh
rosalind analyze evidence \
  --reference genome.fa --alignments sample.sorted.bam --sites candidates.vcf \
  --memory-budget-mb 256 --format arrow-ipc --output evidence.arrow

rosalind verify --manifest evidence.arrow.manifest.json
```

Use `--regions targets.bed` instead of `--sites` for evidence at every selected
position, including zero coverage. Omit `--format arrow-ipc` for TSV. For target
summaries, use `rosalind analyze panel-qc` with the same reference, alignment, and
BED arguments. Add `--plan` to either analysis command to inspect resource planning.

Inputs need BAM+BAI/CSI or CRAM+CRAI; CRAM requires an explicit local FASTA.
The current CRAM path supports checked CRAM 3.0 layouts and validates the whole
file once per run; see [supported layouts and costs](docs/SEMANTICS.md#cram-decoder-admission).
References can be indexed FASTA, `.rref`, or compatible legacy `.idx` files.
The default evidence profile requires MAPQ 20 and base quality 20, and excludes
unmapped, secondary, supplementary, duplicate, and QC-failed reads. It counts
**reads**, including overlapping mates separately. For multiple samples, select
`--sample NAME` or explicitly use `--pool-samples`.
[Scientific semantics](docs/SEMANTICS.md) defines the full contract.

## Use evidence in your own tools

The [Python interface](python/README.md) yields bounded Arrow batches and can
materialize verified artifacts. [Saved-dataset examples](examples/persisted-evidence/README.md)
show Python, R, and SQL queries without reopening the original alignments.
Downstream code is responsible for memory it retains.

For a Rust analyzer, generate a working starting point:

```sh
rosalind new analyzer candidate-qc --api evidence --output ./candidate-qc
```

Follow the [SDK guide](docs/analyzer-sdk.md) to build it against the current source,
implement your reducer, and run conformance checks. The guide includes the local
dependency setup required until the candidate crates are published.

## Resource and verification guarantees

Scientific settings and execution settings are recorded separately. Successful
supported runs must preserve their evidence across admitted budgets, tile layouts,
and worker counts. Caching is opt-in; reused inputs and partitions are content-verified.

Memory planning and cooperative RSS checks are distinct from an OS-enforced cap.
Native validation can allocate before a cooperative check. CRAM planning now
includes a checked container envelope and a complete record-validation pass. Use the
[resource findings](docs/benchmarks-and-limitations.md) to assess measured behavior,
not a universal hard-memory guarantee.

Outputs are create-new and atomic by default. Receipts support byte verification,
comparison, and replay with the matching producer. They establish content identity,
not biological validity or authorship. See [receipts and trust](docs/receipts-and-trust.md).

## Guides and project status

| Start here | Guide |
|---|---|
| Install the right version | [Stable release versus source preview](docs/installation.md) |
| Research workflow | [Researcher quickstart](docs/researcher-quickstart.md) · [Pinned real-data tutorial](examples/research-filter/README.md) |
| Understand the results | [Core concepts](docs/concepts.md) · [Troubleshooting](docs/troubleshooting.md) |
| Variant interoperability | [VCF/BCF annotation](docs/variant-annotation.md) |
| Repeated analysis | [Executable reuse quickstart](docs/reuse-quickstart.md) · [Dataset reference](docs/reusable-evidence.md) |
| Workflow integration | [Integration guide](docs/workflow-integration.md) · [Nextflow](integrations/nextflow/) · [Snakemake](integrations/snakemake/) |
| Builder extension | [Builder quickstart](docs/builder-quickstart.md) · [Rust SDK](docs/analyzer-sdk.md) · [Python](python/README.md) |
| Validation | [Implementation status](docs/implementation-status.md) · [Benchmarks and limits](docs/benchmarks-and-limitations.md) |
| Try a real task and report results | [Researcher and builder validation kit](docs/adoption-validation.md) |
| Contribute | [Contributor guide](CONTRIBUTING.md) · [Roadmap](docs/ROADMAP.md) · [Security](SECURITY.md) · [Citation](CITATION.cff) |

Current scope is short-read DNA evidence and panel QC for research. Indel/haplotype
inference, UMI consensus, and production or clinical calling are outside this
contract. Legacy feature extraction, reference/search utilities, and the column
analyzer API remain available with their own documented semantics.

## License

Dual-licensed under Apache-2.0 and MIT.
