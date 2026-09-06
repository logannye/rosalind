# Rosalind

[![CI](https://github.com/logannye/rosalind/actions/workflows/ci.yml/badge.svg)](https://github.com/logannye/rosalind/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/rosalind-bio?logo=rust&label=crates.io&color=orange)](https://crates.io/crates/rosalind-bio)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

**Reusable genomic evidence, with explicit semantics and verifiable results.**

Rosalind is a Rust library, CLI, and Python interface for turning indexed
short-read alignments into per-locus evidence and target QC. Researchers supply
candidate SNVs or target intervals. Builders supply bounded reducers. The engine
chooses an execution tile that fits the declared resources, counts every eligible
observation, and preserves content identity for verification and replay.

Its technical advantage is separating **scientific selection and reduction from
resource scheduling**: less memory can mean smaller tiles and more indexed I/O,
but cannot silently change a successful scientific result. Fixed output ordering,
integer summaries, canonical Arrow batches, and verified cache partitions make
the same evidence reusable across consumers.

**Development status:** exact evidence, panel QC, verified cache/resume, and dataset
comparison are implemented, with local and CI validation. The source package still has
version 0.4.0; this development work is unpublished, with staged 0.5/0.6/0.7
releases planned. The [retained validation report](docs/findings/evidence-engine-2026-09-05/README.md)
records 467 passing Rust tests, eight Python tests on both Python 3.9 and 3.11,
and fresh macOS arm64 wheel installs. Subsequent [wheel CI](https://github.com/logannye/rosalind/actions/runs/34009965756)
also passed fresh Python 3.9/3.11 installations on Linux x86_64 and both macOS
architectures. See [implementation status](docs/implementation-status.md)
for platform limits and [release gates](docs/ROADMAP.md) for publication and
independent adoption requirements.

Research use only. The evidence profile counts reads, not molecules; the built-in
experimental SNV caller is not a production or clinical caller.

## Try a small real-origin dataset

Build the current source checkout to use these development APIs (Rust 1.83 or
newer and a native C build toolchain).

```sh
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --locked --bin rosalind
export PATH="$PWD/target/debug:$PATH"
```

Follow the [small NA18507 research tutorial](examples/research-filter/README.md)
to download 118KB of content-locked upstream data, generate candidate SNVs with an
independent tool, extract exact evidence, join a research screen, and verify/replay
the artifact. The source provenance is documented; the example is not an accuracy
benchmark.

With existing indexed alignments and an indexed FASTA:

```sh
rosalind analyze evidence \
  --reference genome.fa --alignments sample.sorted.bam --sites candidates.vcf \
  --memory-budget-mb 256 --format arrow-ipc --output evidence.arrow
rosalind verify --manifest evidence.arrow.manifest.json
rosalind reproduce --manifest evidence.arrow.manifest.json --inputs .
```

Use `--regions targets.bed` instead of `--sites` to emit every selected locus,
including zero-depth positions. See [SEMANTICS.md](docs/SEMANTICS.md) for coordinate,
flag, quality, depth, and denominator definitions before comparing another tool.

## Panel QC and shared extraction

```sh
rosalind analyze panel-qc \
  --reference genome.fa --alignments sample.sorted.bam --regions targets.bed \
  --min-callable-depth 10 --memory-budget-mb 256 \
  --position-output positions.arrow --output panel.tsv
```

Each original BED target keeps its full denominator and uncovered positions;
overlapping targets remain distinct. Optional position evidence comes from the
same traversal. Reference-free coverage is available by omitting `--reference`
for panel QC, except that a quality-qualified stored SAM sequence symbol `=`
requires a reference to resolve its base. CRAM also needs a local decoding FASTA,
supplied with `--cram-reference` when the analysis reference is omitted.
Callability thresholds are configurable technical screens and have no intrinsic
clinical interpretation.

## Build on the evidence

The Rust `EvidenceRequest`/`EvidenceEngine`/`EvidenceAnalyzer` API accepts explicit
field, reference, context, and memory requirements. Schema v1 emits complete rows;
field requirements do not yet reduce physical row size. The
[standalone consumer crate](examples/evidence-analyzer/) demonstrates a bounded
integer reducer outside the workspace. The [SDK guide](docs/analyzer-sdk.md)
distinguishes extraction from the artifact/receipt lifecycle.

Python's mixed package is distribution `rosalind-bio`, import `rosalind`. A wheel
bundles its matching CLI. Until published, build a local wheel as described in
[python/README.md](python/README.md).

```python
from rosalind import iter_evidence, materialize_evidence

with iter_evidence("genome.fa", "sample.sorted.bam", sites="candidates.vcf",
                   memory_budget_mb=256) as run:
    rows = 0
    for batch in run:
        rows += batch.num_rows  # at most 1,024 rows per native Arrow batch
    result = run.result
print(f"Extracted {rows} loci")

artifact = materialize_evidence("genome.fa", "sample.sorted.bam", "python-evidence.arrow",
                                sites="candidates.vcf", memory_budget_mb=256)
```

Consumers own any retained Python memory. Exhaustion finalizes a streamed run;
early cancellation is not a successful artifact. Materialize when a persisted,
byte-verifiable output is required.

## Budget, reuse, and trust

Use the actual evidence command with `--plan` for its admitted model:

```sh
rosalind analyze evidence --reference genome.fa --alignments sample.sorted.bam \
  --sites candidates.vcf --memory-budget-mb 256 --format arrow-ipc --plan

rosalind analyze evidence --reference genome.fa --alignments sample.sorted.bam \
  --sites candidates.vcf --memory-budget-mb 256 --format arrow-ipc \
  --cache-dir ./evidence-cache --resume --workers 2 --output evidence-cached.arrow
```

Caching is opt-in. Later runs rehash inputs before verified partition reuse.
Workers extract first-party partitions and reducers consume them in canonical
order. Small budgets may increase repeated reads; there is no universal speedup
claim. Prepared alignment indexes and reference creation remain separate setup.
Keep alignment, reference, selection, and index files immutable for the whole run;
the [input contract](docs/SEMANTICS.md) describes mutation detection and its limits.

File outputs are create-new and atomic by default; `--force` permits replacement.
Failed runs do not publish a successful destination. A cooperative memory model
and sampled RSS are distinct from an OS cap; `--require-os-limit` checks an existing
Linux cgroup-v2 limit. htslib decodes records before checking their declared
envelope, so the software model is not a universal hard allocation guarantee.

Receipts bind content, scientific settings, producer identity, and outcomes.
Offline `verify`, `reproduce`, `diff`, and receipt inspection distinguish portable
claims from machine-local measurements. Unsigned receipts are tamper-evident,
not proof of authorship or biological validity. Historical schemas 1–5 still verify
at their original capabilities. See [receipts and trust](docs/receipts-and-trust.md).

## Compatibility and current limits

| Surface | Current development contract |
|---|---|
| Exact evidence input | Local indexed BAM+BAI/CSI or CRAM+CRAI; CRAM needs explicit local FASTA |
| Reference | Uncompressed FASTA+FAI, `.rref`, compatible `.idx`; optional for panel coverage |
| Evidence schema | Read-based exact SNV counts, allele strands, quality sums/histograms, stored-SEQ offset sums, filter counts |
| Legacy analyzers | `features`, `analyze coverage`, `ColumnAnalyzer`, and scaffold remain; separate filter defaults |
| Legacy capacity | New runs use `pileup.semantics=exact-or-fail-v1`; old semantic replay needs its producer |
| Outputs | TSV, canonical Arrow; legacy first-party shard merge remains supported |
| Platforms | Linux x86_64 and macOS arm64/x86_64 package workflows; publication/CI gates still apply |
| Deferred | Indels, UMI/fragment consensus, methylation, long reads, production calling, Linux ARM wheels |

Read-position summaries use offsets in the stored SEQ, adjusted for reverse
orientation. They do not reconstruct sequencing cycles removed before alignment.

The legacy extension entry point remains:

```sh
rosalind new analyzer locus-qc --output /tmp/locus-qc
```

Its contract runner adds output safety, receipt, and replay handling. Custom state
needs its own conservative memory declaration; an unknown model is observable but
cannot inherit a proven bound. Existing BAM region/shard work retains its original
BAI-only constraints; CSI/CRAM support belongs to the new evidence path.

Analysis does not require a search index. `rosalind reference build --fasta
genome.fa --output genome.rref` builds a reusable packed reference. The legacy
FM-index and `locate` remain separate search facilities.

## Evidence and contribution

The [local validation report](docs/findings/evidence-engine-2026-09-05/README.md)
retains raw test, package, workflow, and benchmark evidence. On the small NA18507
example, all 34 evidence fields matched an independent pysam oracle exactly across
three budgets and three repetitions; cache resume preserved output while avoiding
partition extraction. These checks establish the tested behavior on this dataset.
Larger-data efficiency, Linux/macOS Intel packages, container enforcement, and
independent adoption still have separate gates.

The [exact-evidence curve](benchmarks/evidence/) compares a semantically matched
streaming pysam oracle with three repeated declared-budget runs, retaining raw
time/RSS/I/O, hashes, and verification cost. The [platform harness](benchmarks/platform/)
covers legacy extraction and real cgroup probes when Docker is available.
bcftools default mpileup differences are reported explicitly rather than treated
as equivalence. A harness is not a published performance result.

The [GIAB caller baseline](benchmarks/giab/) remains unestablished. Its scientific
gate follows caller/shared-processing changes; evidence-engine correctness is
tested by exact observation agreement and invariance. See
[benchmarks and limitations](docs/benchmarks-and-limitations.md).

Maintained [Nextflow](integrations/nextflow/) and [Snakemake](integrations/snakemake/)
examples connect analysis to workflow verification. The [roadmap](docs/ROADMAP.md)
tracks independent researcher/builder adoption separately from authored examples.
Historical exploratory plans in `docs/superpowers/`, `docs/OPEN_PROBLEMS.md`, and
`docs/GROWTH.md` are archives, not the current delivery sequence.

Contributions preserve scientific invariance, receipt compatibility, and atomic
publication. See [CONTRIBUTING.md](CONTRIBUTING.md), [SECURITY.md](SECURITY.md), and
[CITATION.cff](CITATION.cff).

## License

Dual-licensed under Apache-2.0 and MIT.
