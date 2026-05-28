# Rosalind

**A deterministic, low-memory genomics engine in Rust — call variants across a whole genome on a laptop, with memory you can predict and verify, and results that are byte-for-byte reproducible.**

Most variant callers use memory that grows with your data, so "will this finish on my machine?" is something you find out the hard way. Rosalind is built around a different promise: **memory as a contract you can see before you commit and verify after the run.** It streams a coordinate-sorted BAM one read at a time, reads the reference from a compact portable index (no second copy of the genome in RAM), keeps its working set proportional to *local read depth* rather than file size, and prints a receipt showing the memory it actually used. Run it twice on the same inputs and you get identical output, bit for bit. And where the evidence is too thin to be sure, it **abstains** instead of guessing.

It's a Rust **library and CLI** you can call directly, extend with plugins, or drive from Python — not a black-box pipeline.

---

## The headline: bounded whole-genome calling from a portable index

You already have an aligner you trust (bwa-mem2, minimap2, or Rosalind's own). Bring Rosalind a **coordinate-sorted BAM** and a **prebuilt index**, and it calls germline variants across **every contig** in bounded, predictable memory:

```bash
# 1. Build a portable, memory-mappable index of your reference — once.
rosalind index --reference genome.fa --output genome.idx

# 2. (Align reads with your favorite aligner and coordinate-sort the BAM.)
#    `rosalind sort` will do the sort deterministically within a memory budget.

# 3. Call germline variants across the WHOLE genome, streaming, bounded.
rosalind variants \
  --index genome.idx \
  --alignments sample.sorted.bam \
  --memory-budget-mb 4096 \
  -o sample.vcf
# ...writes a multi-contig VCF, plus to stderr:
#   memory: peak RSS 412 MiB; max pileup working set 18 KiB
#   wrote reproducibility receipt: sample.vcf.manifest.json
```

What makes this different:

- **Bounded memory, independent of BAM size.** Reads stream one record at a time; peak memory is roughly *the largest contig's reference + the local pileup working set* — not the size of your alignments. A human genome calls comfortably on a laptop.
- **Self-contained.** The reference comes from the `.idx`; you don't need the original FASTA at call time.
- **A memory receipt.** Every run reports its realized peak RSS and max pileup working set — to stderr and into a reproducibility manifest. `--memory-budget-mb` flags a run that exceeds your declared budget *(it records the verdict; it does not yet abort — enforcement is on the roadmap)*.
- **Reproducible + auditable.** Identical inputs produce a byte-identical VCF; a BLAKE3 manifest records the index, the BAM, the output, and the memory used.

---

## What it does today

- **Bounded whole-genome germline calling** — `rosalind variants --index` streams a coordinate-sorted BAM over all contigs of a persisted index, calling SNVs to a multi-contig VCF with a working set bounded by coverage. Calls are calibrated and **abstention-aware** (no confident call → no row, rather than a guess).
- **Build-once, query-many index** — `rosalind index` builds a portable, memory-mapped FM-index over a (multi-contig) reference; `rosalind locate` answers exact-match queries against it in milliseconds. The index is **never rebuilt on load** and is **byte-identical across builds** of the same reference.
- **Alignment** — `rosalind align` builds a Burrows–Wheeler / FM-index over a reference contig and aligns reads via exact-match seeding, deterministic diagonal chaining, and banded affine-gap refinement. Emits SAM or BGZF-compressed BAM.
- **Streaming I/O** — Reads plain or gzip/bgzf-compressed FASTA/FASTQ, auto-detected from the magic bytes; FASTQ can stream from stdin (`-`).
- **Deterministic coordinate sort** — `rosalind sort`: an external merge sort (spills to disk) that orders a BAM by position within a configurable memory budget.
- **Somatic (tumor/normal) calling** — `rosalind somatic` calls somatic SNVs and simple indels from a paired tumor/normal BAM set using a deterministic binomial log-likelihood-ratio model with explicit depth and allele-fraction filters.
- **Truth-set evaluation** — `rosalind eval-somatic` compares a call set against a truth VCF over confident regions (BED), with variant normalization (left-align + trim) and precision / recall / F1.
- **Extensibility** — Implement the `GenomicPlugin` trait to run custom per-block analyses on the same bounded-memory evaluator, or call the PyO3 bindings from Python.
- **Determinism by design** — Primary artifacts are emitted in a canonical, stable order, byte-for-byte identical across repeated runs given identical inputs. See [`docs/determinism.md`](docs/determinism.md).

## Why it matters

Three properties, treated as first-class guarantees rather than nice-to-haves:

1. **Predictable memory.** The bet is that for a growing set of users — sequencing in the field, in the clinic, on a laptop, or at genome scale on modest hardware — *"it fits, and I knew it would"* matters more than raw throughput. Rosalind makes the streaming working set bounded by local coverage and surfaces the realized peak so the bound is **verifiable, not just claimed**.
2. **Reproducibility.** Byte-identical outputs and a per-run BLAKE3 manifest make results auditable — a hard requirement for clinical and regulated pipelines, and a sanity-saver for everyone else.
3. **Honest uncertainty.** Calibrated, abstention-aware calling refuses to emit a call where the evidence is insufficient, instead of papering over it.

Under the hood, Rosalind is also a research vehicle for **space-bounded genomics**: a `~√t` (square-root-space) evaluation framework as a continuous space/time knob — trade time for memory along a curve a declared budget selects. That direction (sublinear-space index *construction*, budget *enforcement*, `rosalind plan`/`verify`) is on the roadmap below; the bounded streaming engine you can use today is the practical foundation it builds on.

## Who it's for

- **Edge, field, and low-resource settings** — sequencing on a laptop or portable device where predictable memory matters more than peak throughput.
- **Reproducibility-sensitive work** — pipelines where byte-identical, auditable outputs are a first-class requirement.
- **Builders** — anyone who wants a hackable Rust genomics engine to embed, extend with plugins, or drive from Python.
- **Teaching and learning** — a readable, end-to-end Rust implementation of FM-index alignment, streaming pileup, and variant calling to study, modify, and extend.

## Current scope (what's single-contig vs. whole-genome today)

- **Whole-genome:** germline variant calling via `rosalind variants --index` (all contigs, streaming, bounded) and exact-match lookup via `rosalind index` / `rosalind locate`.
- **Single-contig:** Rosalind's own **aligner** (`rosalind align`) and the FASTA-based `variants --reference` path operate on one reference contig per run. For whole-genome calling, align with any standard aligner and bring the coordinate-sorted BAM to `variants --index`. (Wiring the *aligner* onto the persisted multi-contig index is a later phase — see the roadmap.)
- Variant calling is **single-sample** (germline) or a **tumor/normal pair** (somatic); calling is SNV-focused, with simple indels in the somatic path.
- The engine runs **single-threaded** today. `--memory-budget-mb` is **record-only** (it reports a verdict but does not yet enforce).

## Roadmap

The core primitive is a streaming, CIGAR-aware pileup column stream; variant calling and custom plugins consume it. Performance work deliberately *follows* the unique capability — the target user needs "it fits and is predictable" before "it's fastest."

- **Phase A (done):** the streaming pileup engine; calibrated, abstention-aware germline SNV calling; tumor/normal somatic calling; spec-valid VCF; a BLAKE3 reproducibility receipt per run.
- **Phase B (done):** streaming gzip/bgzf input; a multi-contig FM-index over the concatenated genome with `(contig, position)` resolution; a build-once, memory-mapped, byte-reproducible persisted index (`rosalind index`/`locate`); zero-copy reference access from the index; and **bounded whole-genome germline calling over a sorted BAM** (`rosalind variants --index`) with a realized-memory receipt.
- **Phase C (next):** memory as an *enforceable* contract — `rosalind plan` (a checkable memory envelope before you commit), budget **enforcement** with graceful degradation (never OOM on a real device), and `rosalind verify`.
- **Later:** sublinear-space index construction (the `~√t` space/time knob across the full curve); the aligner over the persisted multi-contig index (`align --index`, whole-genome alignment); germline indels and richer read QC; deterministic multithreading; a Python binding over the pileup stream.

Target architecture and per-phase specs/plans live in [`docs/superpowers/specs/`](docs/superpowers/specs/) and [`docs/superpowers/plans/`](docs/superpowers/plans/); the guiding thesis is in [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).

---

## Install & build

### Prerequisites
- Rust 1.72+ (`rustup` recommended)
- Native compression headers for BAM I/O: `libbz2-dev` & `liblzma-dev` on Debian/Ubuntu, `brew install bzip2 xz` on macOS
- Python 3.9+ (only for the PyO3 bindings; set `PYO3_PYTHON=/path/to/python` if the default interpreter is unsuitable)

### Build
```bash
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --release
cargo test              # run the full suite
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

### Whole-genome germline calling (the flagship path)

```bash
# Build the index once.
rosalind index --reference genome.fa --output genome.idx

# Call across all contigs from a coordinate-sorted BAM, in bounded memory.
rosalind variants \
  --index genome.idx \
  --alignments sample.sorted.bam \
  --mapq-threshold 20 \
  --memory-budget-mb 4096 \
  -o sample.vcf
```

`variants --index` requires a **coordinate-sorted BAM** (use `rosalind sort` or `samtools sort`). It reads the reference from the index — no `--reference` FASTA needed — and writes a multi-contig VCF plus a memory + reproducibility receipt. `--memory-budget-mb` records (does not yet enforce) a verdict against the realized peak.

### Single-contig alignment + calling

```bash
# Align FASTQ reads to a single reference contig → SAM (stdout) or BAM (to disk).
rosalind align --reference examples/data/ref.fa --reads examples/data/reads.fastq \
  --format bam --output examples/data/alignments.bam

# Call germline SNVs from a single-contig reference + sorted alignments → VCF.
rosalind variants --reference examples/data/ref.fa \
  --alignments examples/data/alignments.sam --mapq-threshold 10
```

Inputs may be plain or gzip/bgzf-compressed (auto-detected); pass `-` to read FASTQ from stdin, e.g. `gzip -dc reads.fastq.gz | rosalind align --reads - --reference ref.fa --format sam`. `align` indexes the first FASTA record (single-contig scope).

### Other subcommands

Run `rosalind <subcommand> --help` for exact flags.

- `rosalind locate --index genome.idx --pattern GATTACA` — exact-match positions in a prebuilt index (memory-mapped, never rebuilt). Exact-match only; seed/chain/extend alignment against the persisted index is a later phase.
- `rosalind sort` — deterministic coordinate sort of a BAM within a memory budget.
- `rosalind somatic` — tumor/normal somatic SNV + simple-indel calling from a paired BAM set over a region.
- `rosalind eval-somatic` — compare a call set to a truth VCF over confident regions.

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

The bounded whole-genome drive (`rosalind::call::call_germline_whole_genome`) wraps the per-contig caller above: it streams a sorted read source over a persisted index's `ReferenceView`, calling every contig in a single pass and returning the calls plus the max working set observed.

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
cargo test --test variants_index    # bounded whole-genome `variants --index` gates
cargo test --test determinism       # byte-identical outputs across repeated runs
cargo test --test space_bounds      # working-set scaling checks for the streaming evaluator
cargo test --test fm_index_props    # property tests: FM-index rank/total invariants vs. naive counts
cargo test --test golden_vcf        # snapshot test for stable VCF rendering
```
Refresh golden snapshots with `ROSALIND_UPDATE_SNAPSHOTS=1 cargo test`. See [`docs/determinism.md`](docs/determinism.md) for the determinism contract and canonicalization rules.

---

## License

Dual-licensed under Apache-2.0 and MIT. Use GitHub Issues for bugs and feature requests.
