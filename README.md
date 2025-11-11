# Rosalind: Accessible genomics for every device.

## Rosalind is a Rust engine that performs genome alignment, streaming variant calling, and custom bioinformatics workloads using **O(√t)** working memory, letting you run billion-step analyses on laptops, tablets, and edge devices without sacrificing deterministic accuracy.

---

## Why Rosalind Matters

- **Democratized genomics** – Align a human genome or call variants from a full sequencing run in under 100 MB of RAM, keeping the entire pipeline within L1/L2 cache on commodity CPUs.
- **Edge-first architecture** – Run clinical diagnostics, infectious-disease surveillance, and population-scale screening on field hardware with limited memory or sporadic connectivity.
- **Composable analytics** – Extend Rosalind with bespoke pipelines (e.g., RNA expression, quality control, phylogenetics) by implementing a single trait or loading from Python.

---

## Core Capabilities

- **FM-index Alignment**  
  A space-efficient Burrows–Wheeler/FM-index implementation partitions references into √t blocks, keeps per-block rank/select structures, and runs backward search while storing only O(√t) state. Billion-base references can be aligned using tens of kilobytes.

- **Streaming Variant Calling**  
  Rosalind constructs pileups on-the-fly, performs Bayesian scoring using O(√t) working space, and emits variant calls as data arrives. Streaming lets you watch coverage, allele fractions, and quality metrics in real time on the edge.

- **Plugin & Python Ecosystem**  
  The `GenomicPlugin` trait lets you implement new analyses (expression quantification, quality control, coverage audits) as block processors. The PyO3 bindings expose the entire engine to Python without duplicating memory.

---

## How Rosalind Achieves O(√t) Space

- **Block decomposition**  
  Every workload is split into √t blocks with deterministic summaries (FM intervals, pileup stats, etc.). Only one block’s replay buffer (O(√t) cells) lives in memory at a time.

- **Height-compressed trees**  
  Summaries from child blocks are merged in an implicit binary tree whose height is O(log√t). A pointerless DFS stores just 2 bits per level; endpoints are recomputed on demand.

- **Streaming ledger**  
  Two bits per block track whether the left/right children have completed. This ensures merges happen once, without storing intermediate trees or full tape state.

- **Workspace pooling**  
  A single reusable allocation is shared across blocks and algorithms (alignment, pileup, plugins), preventing allocator churn while preserving the space bound.

The result: `space_used ≤ block_size + num_blocks + O(log num_blocks) = O(√t)` with deterministic correctness.

---

## Repository Layout

- `src/framework/` – Generic compressed evaluator and configuration helpers.  
- `src/genomics/`
  - `compressed_dna.rs` – 2-bit encoding for DNA plus ambiguity masks.  
  - `fm_index.rs` – Blocked BWT/FM-index + rank/select checkpoints.  
  - `block_alignment.rs` / `bwt_aligner.rs` – Block-wise backward search and orchestrator.  
  - `pileup.rs` / `variant_caller.rs` – Streaming pileups, Bayesian scoring, variant extraction.  
  - `statistics.rs`, `types.rs` – Shared data models.  
- `src/plugin/` – Plugin trait, executor, registry, and RNA-seq example plugin.  
- `src/python_bindings/` – PyO3 bridge for Python users.  
- `src/main.rs` – Command-line interface for alignment and variant-calling tasks.  
- `tests/` – Unit and integration tests for FM-index behaviour, variant calling, and block merges.  
- `examples/` – CLI snippets and benchmarking harnesses.

---

## Install & Build

Requires Rust 1.72+ (edition 2021).

```bash
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo test           # run the full suite
cargo build --release
```

Add Rosalind to your project:

```toml
[dependencies]
rosalind = { path = "./rosalind" }
```

---

## Quick Start

### Align reads against a reference

```rust
use rosalind::genomics::{BWTAligner, AlignmentResult};

fn align_reads(reads: &[Vec<u8>], reference: &[u8]) -> anyhow::Result<Vec<AlignmentResult>> {
    let mut aligner = BWTAligner::new(reference)?;
    let results = aligner.align_batch(reads.iter().map(|r| r.as_slice()))?;
    Ok(results)
}
```

### Streaming variant calling

```rust
use rosalind::genomics::{AlignedRead, StreamingVariantCaller};

fn call_variants(
    reads: Vec<AlignedRead>,
    reference: &[u8],
) -> anyhow::Result<Vec<rosalind::genomics::Variant>> {
    let chrom = std::sync::Arc::from("chr1");
    let reference = std::sync::Arc::from(reference.to_vec().into_boxed_slice());

    let mut caller = StreamingVariantCaller::new(
        chrom,
        reference,
        0,       // region start
        1024,    // bases per block
        10.0,    // minimum quality
        1e-6,    // prior
    )?;
    caller.call_variants(reads).map_err(Into::into)
}
```

### Command-line usage

```bash
# Align a FASTA reference and read list (one read per line)
cargo run --release -- align --reference data/ref.fa --reads data/reads.txt

# Call variants from an alignments file (position \t sequence)
cargo run --release -- variants --reference data/ref.fa --alignments data/aln.tsv --chrom chr7
```

### Python bindings

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

---

## Testing & Benchmarks

- `cargo test` – runs unit + integration tests (alignment, FM-index, variant calling, plugins).  
- `cargo run --example align_genome` – reference alignment demo.  
- `cargo run --example call_variants` – streaming variant calling demo.  
- `cargo test -- --ignored` – optional long-running scale experiments.

---

## Contributing

Rosalind welcomes pull requests that:

- Add domain pipelines (RNA-seq, metagenomics, quality control).  
- Improve FM-index performance (SIMD rank/select, multi-threaded DFS paths).  
- Expand Python bindings or add CLI subcommands for new workflows.

Before submitting, ensure `cargo fmt`, `cargo clippy`, and `cargo test` pass.

---

## License & Support

- Licensed under Apache-2.0 + MIT dual license.  
- File issues or feature requests in GitHub Issues.  
- Join the discussion on the community Discord (link coming soon).

---

**Built for researchers, clinicians, and developers who need genome-scale analysis without genome-scale infrastructure.**


