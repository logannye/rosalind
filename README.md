# Rosalind: Efficient, accessible genomics for every device

**Rosalind** is a Rust engine for genome alignment, streaming variant calling, and custom bioinformatics analytics that runs on commodity or edge hardware. It achieves **O(√t)** working memory, deterministic replay, and drop-in extensibility for new pipelines (Rust plugins or Python bindings). Traditional pipelines often assume tens of gigabytes of RAM, well-provisioned data centers, and uninterrupted connectivity; Rosalind is designed for the opposite: hospital workstations, clinic laptops, field kits, and classrooms.

---

## Quick Overview
- **Core problem**: standard tools such as BWA, GATK, or cloud-centric workflows frequently require >50 GB RAM, full copies of intermediate files, and high-bandwidth storage, placing them out of reach in many hospitals, public-health labs, and teaching environments.
- **Rosalind’s answer**: split workloads into √t blocks, reuse a rolling boundary between blocks, and evaluate a height-compressed tree so memory stays in L1/L2 cache while preserving deterministic results. The entire pipeline fits in well under 100 MB even for whole genomes.
- **How you use it**: run the CLI, embed the Rust APIs, or extend via plugins/Python to build bespoke genomics workflows—ideal for quick-turnaround clinical diagnostics, outbreak monitoring, or courses where students explore real data on laptops.

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

---

## Core Capabilities
- **FM-index Alignment** – Blocked Burrows–Wheeler/FM-index search with per-block rank/select checkpoints uses only O(√t) working state, making it practical to align entire reference genomes on a mid-spec laptop.
- **Streaming Variant Calling** – On-the-fly pileups with Bayesian scoring keep memory bounded while emitting variants live; ideal for remote surveillance, bedside genomics, or interactive notebooks.
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
- `tests/` – Unit and integration tests (alignment, variant calling, space bounds).
- `examples/` – CLI examples and the scale performance benchmark; includes small synthetic datasets for experimentation.

---

## Install & Build
Requires Rust 1.72+.

```bash
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo test           # run the full suite
cargo build --release
```

To use Rosalind in another crate:

```toml
[dependencies]
rosalind = { path = "./rosalind" }
```

Sample data: `examples/data/` contains small FASTA/FASTQ snippets and alignment inputs mirroring the “Quick Start” commands, so you can reproduce the workflows without hunting for external datasets.

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
# Align a FASTA reference and read list (one sequence per line)
cargo run --release -- align --reference examples/data/ref.fa --reads examples/data/reads.txt
# Call variants from an alignments file (tab-separated position and sequence)
cargo run --release -- variants --reference examples/data/ref.fa --alignments examples/data/aln.tsv --chrom chr7
```
Results are printed to stdout in simple text formats that can be redirected into downstream tools or notebooks.

### 4. Python bindings
Great for exploratory analysis, rapid prototyping, or teaching.
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

Optional: enable an RSS regression check with `--features rss` if you want to monitor process RSS in addition to logical counters.

---

## Extend Rosalind
- **Rust plugins**: implement `GenomicPlugin` to process blocks with custom summaries and merges (see `src/plugin/examples.rs`). Handy for adding expression metrics, QC counts, or domain-specific analytics.
- **CLI extensions**: add subcommands in `src/main.rs` to orchestrate new workflows, such as `rosalind coverage` or `rosalind qc`.
- **Python orchestration**: use `rosalind_py.PyGenomicEngine` to blend Rosalind with pandas, scikit-learn, or visualization libraries.
- **Sample datasets**: the `examples/data/` directory shows how to package synthetic or anonymized data for reproducible demos.

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
