# Rosalind

**Deterministic, streaming genomics with O(√t) memory. Whole-genome runs in <100 MB RAM.**

**Rosalind** is a Rust engine for genome alignment, streaming variant calling, and custom bioinformatics analytics that runs on commodity or edge hardware. It achieves **O(√t)** working memory, deterministic replay, and drop-in extensibility for new pipelines (Rust plugins or Python bindings). Traditional pipelines often assume 50-100+ gigabytes of RAM, well-provisioned data centers, and uninterrupted connectivity; Rosalind is designed for the opposite: hospital workstations, clinic laptops, field kits, and classrooms.

---

## Quick Overview
- **Core problem**: standard tools such as BWA, GATK, or cloud-centric workflows frequently require >50 GB RAM, full copies of intermediate files, and high-bandwidth storage, placing them out of reach in many hospitals, public-health labs, and teaching environments.
- **Rosalind’s answer**: split workloads into √t blocks, reuse a rolling boundary between blocks, and evaluate a height-compressed tree so memory stays in L1/L2 cache while preserving deterministic results. The entire pipeline fits in well under 100 MB even for whole genomes.
- **How you use it**: run the CLI, embed the Rust APIs, or extend via plugins/Python to build bespoke genomics workflows—ideal for quick-turnaround clinical diagnostics, outbreak monitoring, or courses where students explore real data on laptops.

See [At a Glance](#at-a-glance), [How It Compares](#how-it-compares), and [What O(√t) memory means](#what-o√t-memory-means-and-what-t-is) for deeper context.

---

## At a Glance
- **O(√t) working memory** – whole-genome runs stay under ~100 MB without lossy approximations.
- **End-to-end deterministic** – outputs are bit-for-bit identical across runs and partition choices.
- **Full-history equivalent** – recomputation keeps results identical to unbounded-memory evaluations.
- **Streaming SAM/BAM/VCF** – standards-compliant outputs without materializing huge intermediates.
- **Edge-ready deployment** – runs on 8–16 GB laptops/desktops so PHI stays on-site.
- **Composable extensions** – plugins/Python bindings inherit the same memory and determinism guarantees.

---

## Why Rosalind Matters
1. **Clinical genomics on laptops** – In many hospitals outside major research centers, the fastest available hardware is a shared desktop with 8–16 GB RAM. Rosalind lets clinicians align patient genomes and call variants in a single shift without renting cloud machines or shipping sensitive data off-site.
2. **Outbreak monitoring at the edge** – During field responses (Ebola, Zika, SARS-CoV-2), portable sequencers and laptops are common, but high-memory servers are not. Rosalind streams reads, calls variants on the fly, and tolerates intermittent connectivity, enabling decisions while still on-site.
3. **Population-scale research** – Universities or small labs may have access to mid-tier clusters but still need to process cohorts of thousands. Rosalind keeps per-sample memory flat, so more genomes can be processed in parallel on cost-efficient instances.
4. **Education and outreach** – Students can experiment with alignment, variant calling, and plugin development on personal machines, making computational genomics coursework more hands-on and equitable.
5. **Custom analytics** – Developers can plug in coverage dashboards, quality metrics, or single-cell style aggregations while benefiting from the same space guarantees, simplifying the path from prototype to production.

---

## Key Guarantees
- **Space bound**: total memory `space_used ≤ block_size + num_blocks + O(log num_blocks) = O(√t)`; only the most recent block boundary is retained.
- **Deterministic replay**: every block is re-simulated from the previous boundary, producing the same results as a full history, even on resource-constrained devices.
- **Composable design**: block processors, plugins, and bindings use the same compressed evaluator, so new analyses inherit the guarantees—perfect for bespoke QC or epidemiology dashboards.
- **Guardrails included**: regression tests (`tests/space_bounds.rs`) and `scripts/run_scale_test.sh` fail if the O(√t) or sublinear scaling properties regress.
- **Partition invariance**: outputs are unchanged across valid choices of block size and chunking; merges are deterministic and order independent.
- **Full-history equivalence**: results match an unbounded-history evaluation; the space savings come from recomputation, not information loss.

### What O(√t) memory means (and what 't' is)
- `t` ≈ total bases processed. For 30× human whole-genome sequencing: coverage `C ≈ 30`, genome size `G ≈ 3.1×10⁹`, so `t ≈ C × G ≈ 9.3×10¹⁰`.
- `√t ≈ 3.0×10⁵`, which sets the block buffer size. The height-compressed merge stack is `log₂(t) ≈ 36` levels—negligible compared with the block.
- Working set ≈ `(α + β) · √t + γ`. With α ≈ 64–128 B per active position, whole-genome runs sit around 30–80 MB; even conservative assumptions keep the bound <100 MB.
- The bound holds because only the current block, rolling boundary, and compact merge stack reside in memory; older state is recomputed on demand.
- The FM-index/reference can be memory-mapped, so the O(√t) claim concerns Rosalind’s dynamic working set relative to dataset size.

---

## Why This Is Different
- **Partition-invariant determinism** – bit-for-bit identical outputs across runs, regardless of block size or partitioning; ideal for clinical audits, SOP lock-down, and incident investigation.
- **Strict, test-enforced O(√t) memory** – whole-genome runs fit in well under 100 MB; CI gates fail if the bound regresses.
- **Full-history equivalence** – recomputation (not truncation) guarantees identical results to an unbounded-memory evaluation.
- **True streaming reads → variants** – no need to materialize large intermediates; reduces IO/storage pressure and time-to-first-result.
- **Standards without heavyweight infra** – interoperable SAM/BAM/VCF emission while retaining streaming and memory advantages.
- **Cache-resident execution** – keeping state in L1/L2 minimizes cache misses and paging on modest hardware, improving real-world throughput outside of large servers.
- **Composable extensions that inherit guarantees** – plugins and Python bindings share the same compressed evaluator and workspace pool, preserving memory bounds and determinism.

### Clinical Relevance
- Keeps PHI on-site—runs comfortably on 8–16 GB hospital desktops and field laptops without cloud transfer.
- Bit-for-bit reproducibility simplifies CAP/CLIA audits, SOP lock-down, and incident review.
- Predictable <100 MB working set avoids OOMs and keeps shared schedulers/workstations stable.
- Streaming-friendly operation tolerates intermittent connectivity and minimizes temp-file churn in the field.

---

## How It Compares
| Capability | Rosalind | Typical Stack |
| --- | --- | --- |
| Peak RAM (WGS) | `<100 MB` working set; no multi-GB temp files | `1–16+ GB` RAM plus large intermediates |
| Determinism | Bit-for-bit identical outputs; partition invariant | Often varies with thread ordering or sharding |
| Partition invariance | Validated in CI across block sizes | Repartitioning can alter outputs |
| Streaming outputs | Reads → SAM/BAM → VCF without materializing huge files | Batch stages typically require full intermediates |
| Standards | SAM/BAM/VCF with streaming-friendly pipeline | Standards supported, but streaming often limited |
| On-prem edge viability | Runs on 8–16 GB laptops; PHI stays on-site | Assumes high-RAM servers or cloud resources |
| Guardrails/tests | CI enforces O(√t), determinism, property tests | Unit tests common; resource/determinism guards rare |

---

## Core Capabilities
- **FM-index Alignment** – Blocked Burrows–Wheeler/FM-index search with per-block rank/select checkpoints uses only O(√t) working state, making it practical to align entire reference genomes on a mid-spec laptop.
- **Streaming Variant Calling** – On-the-fly pileups with Bayesian scoring keep memory bounded while emitting variants live; ideal for remote surveillance, bedside genomics, or interactive notebooks.
- **Standards-compliant outputs** – Interoperable SAM/BAM/VCF with streaming-friendly IO, minimizing large intermediate files and sort/index overhead.
- **Plugin & Python Ecosystem** – Implement the `GenomicPlugin` trait or call into the PyO3 bindings to add custom analyses without duplicating memory—e.g., RNA expression summaries or coverage drop-out checks for diagnostics.
- **Rolling Boundary** – During DFS evaluation only the latest block summary is retained; older summaries are discarded, guaranteeing `O(b + T + log T) = O(√t)`.

---

## How Rosalind Achieves O(√t) Space
1. **Block decomposition** – Workloads are partitioned into √t blocks with deterministic summaries; only one block’s buffer is in memory at a time, avoiding whole-read caches.
2. **Height-compressed trees** – Child summaries merge in an implicit tree of height O(log√t); pointerless DFS stores just 2 bits per level and recomputes endpoints on the fly.
3. **Streaming ledger** – Two bits per block track completion, preventing duplicate merges or stored intermediate trees.
4. **Workspace pooling** – A single reusable allocation is shared across components (alignment, pileups, plugins), avoiding allocator churn and preserving the bound.
5. **Execution flow** – *(reads)* → block alignment → rolling boundary update → tree merge → streaming outputs (variants, metrics, analytics). This pipeline mirrors typical aligner + variant-caller workflows, but with drastically lower memory.

---

## Repository Layout
- `src/framework/` – Generic compressed evaluator and configuration helpers.
- `src/genomics/` – FM-index, alignment summaries, pileups, variant calling, and shared types.
- `src/plugin/` – Plugin trait, executor, registry, and example RNA-seq quantification plugin.
- `src/python_bindings/` – PyO3 bridge exposing Rosalind to Python.
- `src/main.rs` – Command-line interface (align/variants subcommands).
- `examples/` – CLI examples (including `verify_installation.rs`) and the scale performance benchmark; includes small synthetic datasets for experimentation.
- `scripts/` – Utility scripts (toy dataset generation, SAM vs BAM benchmarking).
- `tests/` – Unit, property, and integration tests (alignment, variant calling, space bounds).
- `tests/common/`, `tests/snapshots/` – Snapshot harness and golden fixtures used by CLI/VCF tests.

---

## Install & Build
### Prerequisites
- Rust 1.72+ (`rustup` recommended)
- Python 3.9+ (for PyO3 bindings; set `PYO3_PYTHON=/path/to/python` if the default interpreter is unsuitable)
- Native compression headers for BAM output (`libbz2-dev` & `liblzma-dev` on Debian/Ubuntu, `brew install bzip2 xz` on macOS)

### Build
```bash
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo test           # run the full suite
cargo build --release
```

Verify the CLI is available:

```bash
cargo run --release -- --help
```

Optional: run the smoke test example to confirm the full pipeline:

```bash
cargo run --example verify_installation
```

To use Rosalind in another crate:

```toml
[dependencies]
rosalind = { path = "./rosalind" }
```

Sample data: `examples/data/` contains small FASTA/FASTQ snippets and alignment inputs mirroring the “Quick Start” commands, so you can reproduce the workflows without hunting for external datasets.

Need something bigger? Generate a deterministic ~10× toy genome with:

```bash
python scripts/generate_toy_data.py examples/data/illumina_toy
cat examples/data/illumina_toy/SHA256SUMS
```

The script emits `reference.fa`, paired `reads_R1.fastq`/`reads_R2.fastq`, and reproducible SHA256 checksums so larger demos can be cached or shared safely.

---

## Quick Start
### 1. Align reads (Rust API)
Ideal for embedding Rosalind in your own pipeline or unit tests.
```rust
use rosalind::genomics::{BWTAligner, AlignmentResult};

fn align_reads(reads: &[Vec<u8>], reference: &[u8]) -> anyhow::Result<Vec<AlignmentResult>> {
    let mut aligner = BWTAligner::new(reference)?;
    aligner.align_batch(reads.iter().map(|r| r.as_slice()))
}
```
Each `AlignmentResult` includes FM intervals (candidate positions), mismatch counts, and heuristic scores—sufficient to drive downstream filtering or assembly logic.

### 2. Call variants (Rust API)
```rust
use rosalind::genomics::{AlignedRead, StreamingVariantCaller};

fn call_variants(reads: Vec<AlignedRead>, reference: &[u8]) -> anyhow::Result<Vec<rosalind::genomics::Variant>> {
    let chrom = std::sync::Arc::from("chr1");
    let reference = std::sync::Arc::from(reference.to_vec().into_boxed_slice());
    let mut caller = StreamingVariantCaller::new(chrom, reference, 0, 1024, 10.0, 1e-6)?;
    caller.call_variants(reads)
}
```
Each `Variant` reports position, reference/alternate alleles, allele fraction, and quality—enough to feed clinical reporting, surveillance dashboards, or QC scripts.

### 3. Command line
Perfect for pilots, demos, or quick analyses with sample data.
```bash
# Align FASTQ reads against the bundled reference and capture SAM output
cargo run --release -- align \
  --reference examples/data/ref.fa \
  --reads examples/data/reads.fastq \
  --format sam \
  --max-mismatches 2 \
  --reference-offset 0 > examples/data/alignments.sam

# Emit coordinate-sorted BAM directly to disk (recommended for downstream tooling)
cargo run --release -- align \
  --reference examples/data/ref.fa \
  --reads examples/data/reads.fastq \
  --format bam \
  --output examples/data/alignments.bam

# Call variants from the SAM/BAM alignments (VCF to stdout by default)
cargo run --release -- variants \
  --reference examples/data/ref.fa \
  --alignments examples/data/alignments.sam \
  --mapq-threshold 10 \
  --region-start 0

# Or write the VCF to disk
cargo run --release -- variants \
  --reference examples/data/ref.fa \
  --alignments examples/data/alignments.sam \
  --output examples/data/variants.vcf
```

Key knobs:
- `--max-mismatches` bounds per-read Hamming distance during seeding.
- `--format {sam|bam}` toggles between plain-text SAM and BGZF-compressed BAM (requires `--output` for BAM).
- `--mapq-threshold` filters low-confidence alignments before variant calling.
- `--region-start` allows offsetting reported genomic coordinates for tiled analyses.
- `--output/-o` writes results to disk instead of stdout for both subcommands.

### 4. Python bindings
Great for exploratory analysis, rapid prototyping, or teaching.

Install the bindings with [maturin](python/README.md):

```bash
pip install maturin
maturin develop --release
```

```python
from rosalind_py import PyGenomicEngine
engine = PyGenomicEngine()
print(engine.list_plugins())

# region_start, region_end, reads[(pos, seq)], block_size
depth = engine.run_rna_seq_plugin(
    region_start=100_000,
    region_end=101_000,
    reads=[(100_020, "ACGTACGT"), (100_050, "TTTACGT")],
    block_size=512,
)
```
The example plugin emits per-base coverage suitable for expression quantification or QC charts.

---

## Verify the Guarantees
| Command | Purpose |
| --- | --- |
| `cargo test --test space_bounds` | Verifies O(√t) bound, component limits, and √t scaling ratios; includes assert-based checks. |
| `./scripts/run_scale_test.sh` | Runs the long-form benchmark; exits non-zero if the bound or sublinear scaling checks are violated (use `--csv` to capture output). |
| `cargo run --example scale_performance_test` | Prints per-scale metrics and component breakdowns for manual inspection (leaf buffer, stack depth, ledger). |
| `cargo test` | Runs all unit/integration tests covering alignment, variant calling, plugins, and supporting modules. |
| `cargo test --test determinism` | Confirms bit-for-bit identical outputs across repeated runs given the same inputs/params. |
| `cargo test --test fm_index_props` | Property tests that validate FM-index rank/total invariants against naive counting. |
| `cargo test --test golden_vcf` | Snapshot test for stable VCF rendering; refresh with `ROSALIND_UPDATE_SNAPSHOTS=1`. |

### Reproducibility & Validation
- `cargo test --test determinism` — locks down partition-invariant, bit-for-bit outputs for auditability and SOPs.
- `cargo test --test fm_index_props` — property-based checks for FM-index rank/total invariants.
- `ROSALIND_UPDATE_SNAPSHOTS=1 cargo test` — refreshes golden VCF/SAM outputs when expected results change; the default run verifies no drift.
- Together with CI enforcement of the O(√t) bound, these tests provide a defensible validation story for regulated or clinical environments.

Optional: enable an RSS regression check with `--features rss` if you want to monitor process RSS in addition to logical counters.

### Benchmarks & Snapshots
- `python scripts/benchmark_formats.py --release` – compare SAM vs BAM throughput using the bundled toy dataset (it will be generated on first run).
- Refresh golden outputs with `ROSALIND_UPDATE_SNAPSHOTS=1 cargo test`.

### Notes on Guarantees and Trade-offs
- Variant scoring remains statistical (Bayesian), but execution is deterministic—no sampling or thread-order nondeterminism—simplifying validation.
- FM-index seeding is exact over the indexed sequence; MAPQ/filters are explicit and covered by tests.
- Recomputing block boundaries to preserve the space bound can add CPU vs server-optimized pipelines; in exchange, you get predictable, cache-friendly performance on commodity hardware.

---

## Extend Rosalind
- **Rust plugins**: implement `GenomicPlugin` to process blocks with custom summaries and merges (see `src/plugin/examples.rs`). Handy for adding expression metrics, QC counts, or domain-specific analytics.
- **CLI extensions**: add subcommands in `src/main.rs` to orchestrate new workflows, such as `rosalind coverage` or `rosalind qc`.
- **Python orchestration**: use `rosalind_py.PyGenomicEngine` to blend Rosalind with pandas, scikit-learn, or visualization libraries.
- **Sample datasets**: the `examples/data/` directory shows how to package synthetic or anonymized data for reproducible demos.

---

## Troubleshooting
- `cargo run` fails with “bin not found” → rerun `cargo build --release` and ensure you invoke `cargo run --release -- …`.
- CLI complains about missing files → verify the sample data exists in `examples/data/` or provide your own reference/reads.
- `ImportError: rosalind_py` → run `maturin develop --release` (see [python/README.md](python/README.md)) or export `PYO3_PYTHON=/path/to/python3.9+`.
- Build errors on install → confirm `rustc --version` ≥ 1.72 and `cargo clean` before rebuilding.
- Writing BAM without `--output` → specify a filename via `-o`/`--output`; BAM is binary and cannot be streamed to stdout.

---

## Contributing
Rosalind welcomes pull requests that:
- Add domain pipelines (RNA-seq, metagenomics, QC dashboards).
- Improve FM-index performance (SIMD rank/select, multi-threaded DFS paths).
- Expand Python bindings or add CLI subcommands for new workflows.

Before submitting, ensure `cargo fmt`, `cargo clippy`, and `cargo test` pass.

---

## License & Support
- Licensed under Apache-2.0 + MIT dual license.
- Use GitHub Issues for bugs and feature requests.
- Community Discord coming soon (ideal for sharing plugins, datasets, or course material).

---

**Built for researchers, clinicians, educators, and developers who need genome-scale analysis without genome-scale infrastructure.**
