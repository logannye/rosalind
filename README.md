# Rosalind

A deterministic, low-memory genomics engine in Rust for read alignment and variant calling on commodity hardware.

Rosalind streams variant calling over coordinate-sorted alignments with a working set bounded by local read coverage rather than by file size, and emits BAM/VCF that are designed to be bit-for-bit reproducible across runs. It is built to be embedded and extended: a Rust library and CLI you can call directly, add plugins to, or drive from Python.

---

## What it does today

- **Alignment** — Builds a Burrows–Wheeler / FM-index (SA-IS suffix array with blocked rank/select) over a reference contig and aligns reads via exact-match seeding, deterministic diagonal chaining, and banded affine-gap refinement. Emits SAM or BGZF-compressed BAM.
- **Streaming I/O** — Reads plain or gzip/bgzf-compressed FASTA/FASTQ, auto-detected from the magic bytes; FASTQ can stream from stdin (`-`). Parsing lives in the `io::` layer (`io::fasta`, `io::fastq`).
- **Coordinate sort** — A deterministic external merge sort (spills to disk) that orders a BAM by position within a configurable memory budget.
- **Germline variant calling** — Streams a pileup over a coordinate-sorted BAM and calls SNVs to VCF, keeping the in-memory working set proportional to coverage, not to the size of the input.
- **Somatic (tumor/normal) calling** — Calls somatic SNVs and simple indels from a paired tumor/normal BAM set using a deterministic binomial log-likelihood-ratio model with explicit depth and allele-fraction filters.
- **Truth-set evaluation** — Compares a call set against a truth VCF over confident regions (BED), with variant normalization (left-align + trim) and precision / recall / F1.
- **Extensibility** — Implement the `GenomicPlugin` trait to run custom per-block analyses on the same bounded-memory evaluator (an example RNA-seq coverage plugin is included), or call the PyO3 bindings from Python.
- **Determinism** — Primary artifacts are emitted in a canonical, stable order and are designed to be byte-for-byte identical across repeated runs given identical inputs and configuration. See [`docs/determinism.md`](docs/determinism.md).

## Current scope

At the **command line**, Rosalind currently operates on a **single reference contig per run**, reads plain or **gzip/bgzf-compressed** FASTQ/FASTA (auto-detected, including from stdin), and runs **single-threaded**. Variant calling is **single-sample** (germline) or a **tumor/normal pair** (somatic); calling is SNV-focused, with simple indels in the somatic path. Alignment uses exact-match seeding. The FM-index is built in memory at the start of each run (memory proportional to the reference); the bounded-memory property applies to the streaming pileup and variant-calling stages. These boundaries define what the engine targets well today — small-to-moderate references, targeted regions, and per-sample streaming workloads.

The library also provides a multi-contig FM-index over the concatenated genome (`genomics::GenomeIndex`) that resolves matches to `(contig, position)`; it is not yet used by the CLI. See the roadmap below.

## Who it's for

- **Edge, field, and low-resource settings** — sequencing on a laptop or portable device where large servers aren't available and predictable memory matters more than peak throughput.
- **Reproducibility-sensitive work** — pipelines where byte-identical, auditable outputs are a first-class requirement.
- **Teaching and learning** — a readable, end-to-end Rust implementation of FM-index alignment, streaming pileup, and variant calling to study, modify, and extend.
- **Builders** — anyone who wants a hackable Rust genomics engine to embed, extend with plugins, or drive from Python, rather than a black-box pipeline.
- **Somatic SNV / simple-indel exploration** on small references and targeted regions.

## Roadmap

The core primitive is a streaming, CIGAR-aware pileup column stream; variant calling and custom plugins consume it.

- **Phase A (done):** the streaming pileup engine; calibrated, abstention-aware germline SNV calling; tumor/normal somatic SNV calling; spec-valid VCF output; a BLAKE3 reproducibility receipt per run.
- **Phase B (in progress):** streaming gzip/bgzf input and a multi-contig FM-index over the concatenated genome (`genomics::GenomeIndex`, with `(contig, position)` resolution and boundary-aware exact-match lookup) have landed. Next: wiring multi-contig through the CLI (whole-genome alignment and calling), a build-once memory-mapped index (`rosalind index`), and pipe-native composition across subcommands.
- **Later:** germline indel calling and richer read QC; deterministic multithreading with an enforced memory budget; a Python binding exposing the pileup stream.

Target architecture and per-phase plans: [`docs/superpowers/specs/`](docs/superpowers/specs/), [`docs/superpowers/plans/`](docs/superpowers/plans/).

---

## Install & build

### Prerequisites
- Rust 1.72+ (`rustup` recommended)
- Native compression headers for BAM output: `libbz2-dev` & `liblzma-dev` on Debian/Ubuntu, `brew install bzip2 xz` on macOS
- Python 3.9+ (only for the PyO3 bindings; set `PYO3_PYTHON=/path/to/python` if the default interpreter is unsuitable)

### Build
```bash
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --release
cargo test            # run the full suite
cargo run --release -- --help
```

Use Rosalind as a library in another crate:
```toml
[dependencies]
rosalind = { path = "./rosalind" }
```

Bundled sample data lives in `examples/data/` (small FASTA/FASTQ + alignments) so the commands below run without external downloads. For a larger, deterministic toy dataset:
```bash
python scripts/generate_toy_data.py examples/data/illumina_toy
```

---

## Use it

### Command line

```bash
# 1. Align FASTQ reads to a reference contig → SAM (stdout) or BAM (to disk)
cargo run --release -- align \
  --reference examples/data/ref.fa \
  --reads examples/data/reads.fastq \
  --format sam \
  --max-mismatches 2 > examples/data/alignments.sam

cargo run --release -- align \
  --reference examples/data/ref.fa \
  --reads examples/data/reads.fastq \
  --format bam \
  --output examples/data/alignments.bam

# 2. Call germline SNVs from the alignments → VCF (stdout, or --output FILE)
cargo run --release -- variants \
  --reference examples/data/ref.fa \
  --alignments examples/data/alignments.sam \
  --mapq-threshold 10
```

Other subcommands (run `rosalind <subcommand> --help` for exact flags):
- `rosalind sort` — deterministic coordinate sort of a BAM within a memory budget.
- `rosalind somatic` — tumor/normal somatic SNV + simple-indel calling from a paired BAM set over a region.
- `rosalind eval-somatic` — compare a call set to a truth VCF over confident regions.

`align` indexes the first FASTA record; additional records are ignored with a warning (single-contig scope). `variants` reads coordinate-sorted SAM/BAM alignments. Inputs may be plain or gzip/bgzf-compressed (auto-detected); pass `-` to read FASTQ from stdin, e.g. `gzip -dc reads.fastq.gz | rosalind align --reads - --reference ref.fa --format sam`.

### Rust API

```rust
use rosalind::genomics::{AlignerError, AlignmentResult, BWTAligner};

// Align reads to a single reference contig via exact-match FM-index seeding.
fn align_reads(reads: &[Vec<u8>], reference: &[u8]) -> Result<Vec<AlignmentResult>, AlignerError> {
    let mut aligner = BWTAligner::new(reference)?;
    aligner.align_batch(reads.iter().map(|r| r.as_slice()))
}
```

```rust
use std::sync::Arc;

use rosalind::call::{call_germline_region, GermlineCall, GermlineParams};
use rosalind::core::{AlignedRead, CoreError, Locus};
use rosalind::pileup::{PileupParams, SliceSource};

// Calibrated, abstention-aware germline SNV calls over one contig, from
// coordinate-sorted reads. Each emitted site is (locus, ref_base, call).
fn call_germline(
    reads: Vec<AlignedRead>,
    reference: Arc<[u8]>,
) -> Result<Vec<(Locus, u8, GermlineCall)>, CoreError> {
    let region = 0..reference.len() as u32;
    call_germline_region(
        SliceSource::new(reads),
        reference,
        0, // contig id
        region,
        PileupParams::default(),
        &GermlineParams::default(),
    )
}
```

```rust
use rosalind::core::Locus;
use rosalind::genomics::{GenomeIndex, GenomeIndexError};

// Multi-contig exact match: build an index over the concatenated genome; hits
// resolve to (contig, position), and matches that would straddle a contig
// boundary are rejected.
fn locate(query: &[u8]) -> Result<Vec<Locus>, GenomeIndexError> {
    let index = GenomeIndex::from_named_sequences(&[
        ("chr1".to_string(), b"ACGTACGT".to_vec()),
        ("chr2".to_string(), b"TTTTGGGG".to_vec()),
    ])?;
    Ok(index.locate_exact(query, 16))
}
```

### Python

```bash
pip install maturin
maturin develop --release
```
```python
from rosalind_py import PyGenomicEngine

engine = PyGenomicEngine()
print(engine.list_plugins())

# Per-base coverage over a region via the example plugin.
depth = engine.run_rna_seq_plugin(
    region_start=100_000,
    region_end=101_000,
    reads=[(100_020, "ACGTACGT"), (100_050, "TTTACGT")],
    block_size=512,
)
```

---

## Extend

- **Rust plugins** — implement `GenomicPlugin` (see `src/plugin/examples.rs`) to run custom per-block analyses (coverage, QC counts, domain-specific summaries) on the same bounded-memory evaluator.
- **CLI subcommands** — add workflows in `src/main.rs`.
- **Python** — drive the engine from `rosalind_py.PyGenomicEngine` alongside pandas / NumPy / scikit-learn.

---

## Tests

```bash
cargo test                          # full unit + integration suite
cargo test --test determinism       # byte-identical outputs across repeated runs
cargo test --test space_bounds      # working-set scaling checks for the streaming evaluator
cargo test --test fm_index_props    # property tests: FM-index rank/total invariants vs. naive counts
cargo test --test golden_vcf        # snapshot test for stable VCF rendering
```
Refresh golden snapshots with `ROSALIND_UPDATE_SNAPSHOTS=1 cargo test`. See [`docs/determinism.md`](docs/determinism.md) for the determinism contract and canonicalization rules.

---

## License

Dual-licensed under Apache-2.0 and MIT. Use GitHub Issues for bugs and feature requests.
