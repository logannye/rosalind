# Rosalind documentation

Rosalind turns indexed short-read alignments into exact read evidence and target
QC that you can verify and reuse. The guides here describe the **0.5.0 source
preview**. The public stable release is **0.1.0**; it does not include these
evidence commands. Start with [installation](installation.md) to choose the right
version.

The separate **0.6 development branch** also has a [cohort CLI preview](cohort-cli-preview.md)
and [complete synthetic reanalysis example](../examples/cohort-reanalysis/README.md).
These commands are not included in the 0.5 release candidate.

## Choose your first task

| I want to… | Start here | Result |
|---|---|---|
| Inspect candidate SNVs or target coverage | [Researcher quickstart](researcher-quickstart.md) | Evidence, a panel summary, and receipts |
| Try Rosalind without preparing my own data | [Pinned NA18507 tutorial](../examples/research-filter/README.md) | Four supplied candidate SNVs and an illustrative research screen |
| Query evidence after the original inputs are unavailable | [Reuse quickstart](reuse-quickstart.md) | Relocated dataset, verified candidate evidence, and panel QC |
| Write an analyzer | [Builder quickstart](builder-quickstart.md) | A standalone reducer with lifecycle handling and conformance checks |
| Work in Python | [Python interface](../python/README.md) | Bounded Arrow batches or materialized artifacts |
| Resolve an installation or analysis failure | [Troubleshooting](troubleshooting.md) | Checks matched to common symptoms |

## Understand the contract

- [Core concepts](concepts.md): reads, selection, fields, zero versus missing, resources, and receipts.
- [Scientific semantics](SEMANTICS.md): coordinates, exact filters, sample scope, panel denominators, and supported CRAM layouts.
- [Receipts and trust](receipts-and-trust.md): content verification, replay, and what neither establishes.
- [Benchmarks and limitations](benchmarks-and-limitations.md): measured results and claims that remain unestablished.

## Extend a working analysis

- [VCF/BCF annotation](variant-annotation.md) preserves supplied records while adding evidence.
- [Portable dataset reference](reusable-evidence.md) covers field compatibility, missing-locus reuse, and export.
- [Python, R, and SQL examples](../examples/persisted-evidence/README.md) consume saved evidence.
- [Workflow integration](workflow-integration.md) explains the maintained Nextflow and Snakemake routes and the legacy GitHub Action.
- [Rust SDK reference](analyzer-sdk.md) defines fields, memory, lifecycle, conformance, and replay contracts.
- [Researcher and builder validation](adoption-validation.md) records whether a real task was completed and what work was avoided.

## Follow development

[Implementation status](implementation-status.md) distinguishes implemented,
verified, released, and adopted capabilities. The [delivery roadmap](ROADMAP.md)
and [adoption sequence](ADOPTION_ROADMAP.md) track remaining gates. For changes,
see the [contributor guide](../CONTRIBUTING.md), [changelog](../CHANGELOG.md), and
[security policy](../SECURITY.md).

Maintainers can use the [release procedures](MAINTAINER_RELEASES.md) and
[historical context](history.md). Archived specifications are background reading,
not prerequisites for either quickstart.
