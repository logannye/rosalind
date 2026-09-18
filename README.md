# Rosalind

[![CI](https://github.com/logannye/rosalind/actions/workflows/ci.yml/badge.svg)](https://github.com/logannye/rosalind/actions/workflows/ci.yml)
[![Latest published release](https://img.shields.io/github/v/release/logannye/rosalind)](https://github.com/logannye/rosalind/releases/latest)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

**Reusable genomic evidence for research and analysis tools.**

Rosalind turns indexed DNA alignments into exact read evidence that you can inspect,
save, and reuse. Researchers use it to review supplied candidate SNVs and measure
coverage across target regions. Builders use its CLI, Python interface, or Rust SDK
to create their own reports and analyzers, with evidence extraction, resource
planning, and verifiable artifacts already provided.

```text
Indexed BAM/CRAM + reference + candidate sites or target regions
    → exact evidence → candidate review, panel QC, or your analyzer
    → portable dataset → another analysis without reopening the alignments
```

**Availability:** the workflows below are in the **0.5.0 source preview**. The latest
public stable release is **v0.1.0**, which predates this evidence engine. Start with
[installation](docs/installation.md); current examples require a source build.
[Implementation status](docs/implementation-status.md) distinguishes implemented,
validated, published, and independently used capabilities. Cohort CLI/Python APIs
are implemented in the separate, **unpublished 0.6 preview** [PR #145](https://github.com/logannye/rosalind/pull/145),
with a [source-bound validation report](https://github.com/logannye/rosalind/blob/codex/cohort-preview/docs/findings/cohort-preview-2026-09-18/README.md).
Current `main` remains the 0.5 line and does not include those cohort commands.

## Choose your starting point

- **[Analyze my data](docs/researcher-quickstart.md)** — extract candidate evidence, understand the counts, verify, and replay.
- **[Build an analyzer](docs/builder-quickstart.md)** — customize a working Rust reducer and run it on native or saved evidence.

No data ready? The [small real-data tutorial](examples/research-filter/README.md)
prepares four supplied candidate SNVs from 118 KB of pinned public inputs.

## What you can build

| Use case | Who it helps | What Rosalind supplies |
|---|---|---|
| Candidate-review report | Researchers revisiting supplied SNVs; developers building review tools | Exact allele and strand counts, quality summaries, and optional record-preserving VCF/BCF annotation |
| Reusable panel analysis | Workflow maintainers and developers answering new questions over the same target regions | Saved evidence for another candidate list, panel QC, bounded Python batches, or a custom Rust reducer |

The useful boundary is **save evidence once, then reuse the measured loci and fields**.
A later query can select different stored positions or SNV alleles and compute
compatible summaries. It cannot recover excluded reads, unmeasured positions, or
molecular information that was never stored.

## Builder quickstart: a working analyzer

Install Rust 1.83+ and the [native build prerequisites](docs/analyzer-sdk.md#build-prerequisites).
Clone the repository and build the source CLI:

```sh
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --locked --bin rosalind --target-dir target
export PATH="$PWD/target/debug:$PATH"
rosalind --version
```

Then build and check the included candidate-summary analyzer. These commands run
from the repository root and need no user data. Its local dependencies already
select this source checkout:

<!-- smoke:readme-builder -->
```sh
cargo test --locked --manifest-path examples/evidence-analyzer/Cargo.toml --target-dir target
cargo build --locked --manifest-path examples/evidence-analyzer/Cargo.toml --target-dir target
rosalind conformance analyzer --api evidence \
  --binary ./target/debug/rosalind-example-evidence-analyzer --json
```

The analyzer reports selected loci, callable read observations, and supplied ALT
observations. Its tests check the statistic; the 19 conformance checks exercise
execution behavior, including native/saved evidence agreement. The
[builder quickstart](docs/builder-quickstart.md) runs it on the four-candidate
fixture, adds a zero-coverage statistic, and verifies and replays its output.

To start a separate project, use `rosalind new analyzer candidate-qc --api evidence
--output ./candidate-qc`, then apply the [source dependency patches](docs/analyzer-sdk.md#candidate-dependency-setup)
before its first Cargo build. The generated project pins the SDK version; the
current evidence SDK is not yet available from the public registry.

Your reducer declares the fields it needs, its retained memory, and its scientific
parameters. Rosalind manages input validation, resource admission, cancellation,
output publication, and receipts. Your tests establish the scientific correctness
of the new statistic. See the [SDK reference](docs/analyzer-sdk.md) for that contract.

Prefer Python? Follow the [Python source installation](python/README.md) to install
`rosalind-bio` (import `rosalind`) with its matching native binary. You can stream
bounded Arrow batches or use the saved-evidence example below.

## Example 1: evidence for a candidate-review tool

With an indexed reference, indexed alignments, and a supplied SNV VCF:

<!-- smoke:readme-candidates -->
```sh
rosalind analyze evidence \
  --reference genome.fa --alignments sample.bam --sites candidates.vcf \
  --memory-budget-mb 256 --output candidate-evidence.tsv

rosalind verify --manifest candidate-evidence.tsv.manifest.json
```

Join the resulting counts to your candidate list to build a review table. The
[complete tutorial](examples/research-filter/README.md) includes that join,
interpretation of callable depth and ALT strand support, and an illustrative
research screen. [VCF/BCF annotation](docs/variant-annotation.md) can also add
evidence while preserving input records, genotypes, and allele order.

## Example 2: a second question from saved evidence

Extract all positions in your target BED once. Saving only the original candidate
sites would leave the other target positions unmeasured:

<!-- smoke:readme-save -->
```sh
rosalind analyze evidence \
  --reference genome.fa --alignments sample.bam --regions targets.bed \
  --cache-dir evidence-cache --memory-budget-mb 256 \
  --format arrow-ipc --output panel-evidence.arrow
```

With the matching Python package installed, use the portable manifest recorded in
the receipt to query another candidate list and produce a panel summary:

<!-- smoke:readme-reuse -->
```python
import json
from pathlib import Path
from rosalind import open_dataset

receipt = json.loads(Path("panel-evidence.arrow.manifest.json").read_text())
dataset = open_dataset(receipt["measurements"]["execution.evidence_dataset_manifest"])
dataset.verify()
dataset.materialize("second-candidates.tsv", format="tsv",
                    sites="second-candidates.vcf", fields=["depths", "alleles"],
                    memory_budget_mb=256)
dataset.panel_qc("targets.bed", "panel-qc.tsv", min_callable_depth=10,
                 memory_budget_mb=256)
```

Both queries use saved evidence without opening the original alignment or reference.
All requested loci and fields must be present; absent evidence is an error, while
stored zero depth remains an observed zero. The depth threshold is a technical
screen, not biological confidence. Both outputs receive verification/replay receipts.
Use fresh output names for each run. The [reuse tutorial](docs/reuse-quickstart.md)
checks equality with fresh extraction and relocates the dataset; copy the entire
portable dataset directory, not just its manifest. [Python, R, and SQL examples](examples/persisted-evidence/README.md)
show additional consumers and Parquet export.

## What makes the engine useful

- **Exact results across admitted execution settings.** Memory budgets change work
  scheduling and indexed I/O, while successful runs preserve the scientific result.
  The engine never silently downsamples; insufficient resources cause refusal or
  an explicit failure. Cooperative accounting is distinct from an OS-enforced cap.
- **Reusable local evidence.** Verified datasets preserve the measured counts and
  their scientific settings, so compatible analyses can run without the original
  alignments. Arrow and TSV support native byte replay; Parquet supports downstream
  analysis and integrity verification.
- **Managed analyzer execution.** The Rust SDK supplies canonical batches and a
  shared lifecycle for native and saved evidence. Python exposes the same native
  engine; memory retained by Python consumers remains their responsibility.

## Inputs, outputs, and current limits

Inputs are indexed BAM (BAI/CSI) or CRAM (CRAI), a reference, and candidate SNVs
(VCF/VCF.gz/BCF) or BED target intervals. CRAM requires a local FASTA and supports
[checked CRAM 3.0 layouts](docs/SEMANTICS.md#cram-decoder-admission), with a complete
file-validation pass per native run. Outputs include Arrow IPC, TSV, annotated
variants, portable evidence datasets, and dataset-derived Parquet.

The default evidence profile requires MAPQ 20 and base quality 20 and excludes
unmapped, secondary, supplementary, duplicate, and QC-failed reads. It counts
**reads**, including overlapping mates separately. Select `--sample NAME` for
multi-sample alignments or explicitly use `--pool-samples`. This is research
short-read DNA evidence; it does not provide UMI consensus, haplotype inference,
or production/clinical variant calling. Legacy APIs retain separate semantics.

Receipts establish content identity and reproducibility with the matching producer;
they do not establish biological validity or authorship. Resource behavior and
supported scope are documented in [scientific semantics](docs/SEMANTICS.md),
[benchmarks and limitations](docs/benchmarks-and-limitations.md), and
[receipts and trust](docs/receipts-and-trust.md).

## Continue

| Need | Guide |
|---|---|
| Install or troubleshoot | [Installation](docs/installation.md) · [Concepts](docs/concepts.md) · [Troubleshooting](docs/troubleshooting.md) |
| Build and integrate | [Builder quickstart](docs/builder-quickstart.md) · [Rust SDK](docs/analyzer-sdk.md) · [Python](python/README.md) · [Nextflow/Snakemake](docs/workflow-integration.md) |
| Check maturity or contribute | [Implementation status](docs/implementation-status.md) · [Roadmap](docs/ROADMAP.md) · [Contributor guide](CONTRIBUTING.md) |
| Try your own task and report results | [Researcher and builder validation kit](docs/adoption-validation.md) |
| Watch or rerun a complete workflow | [Evidence-reuse demonstration](docs/evidence-reuse-demo.md) |
| Find everything else | [Documentation index](docs/index.md) · [Security](SECURITY.md) · [Citation](CITATION.cff) |

## License

Dual-licensed under Apache-2.0 and MIT.
