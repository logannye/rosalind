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

# 3. Will it fit in 4 GB? Ask before committing a byte.
rosalind plan --index genome.idx --max-depth 1000 --budget-mb 4096

# 4. Call germline variants across the WHOLE genome, streaming, honoring the budget.
rosalind variants \
  --index genome.idx \
  --alignments sample.sorted.bam \
  --memory-budget-mb 4096 --enforce \
  -o sample.vcf
# ...writes a multi-contig VCF, plus to stderr:
#   memory: peak RSS 412 MiB; max pileup working set 18 KiB
#   contract: OK — realized peak 412 MiB within declared 4096 MiB
#   wrote reproducibility receipt: sample.vcf.manifest.json

# 5. Re-check the receipt later — no re-run — to prove it fit and is reproducible.
rosalind verify --manifest sample.vcf.manifest.json
```

What makes this different:

- **Bounded memory, independent of BAM size.** Reads stream one record at a time; peak memory is roughly *the largest contig's reference + the local pileup working set* — not the size of your alignments. A human genome calls comfortably on a laptop.
- **Self-contained.** The reference comes from the `.idx`; you don't need the original FASTA at call time.
- **A contract, honored.** `rosalind plan` predicts the peak *before you commit a byte*; `--enforce` honors the budget — refusing up front (exit 3) or failing loud (exit 4) rather than silently OOM-killing you; `rosalind verify` re-checks the receipt without re-running. Without `--enforce`, the budget is record-only. The full story: [the memory contract](CONTRACT.md).
- **Reproducible + auditable.** Identical inputs produce a byte-identical VCF; a BLAKE3 manifest records the index, the BAM, the output, and the memory used.

**Proof — the contract on a real genome.** On the real *E. coli* K-12 MG1655 chromosome (4,641,652 bp, 30× simulated reads), a declared **256 MiB** budget *fits* — `plan` → `variants --enforce` → `verify: OK`, **realized peak 22 MiB** — while an **8 MiB** budget is *refused up front* (exit 3, no work). The contract honored both ways; this demo's claim is *memory* (the reads are simulated) — calling accuracy is measured separately ([Accuracy](#accuracy)). Full numbers + one-command reproduction (`bash scripts/flagship_ecoli_demo.sh`): [`docs/findings/2026-06-01-flagship-ecoli-contract.md`](docs/findings/2026-06-01-flagship-ecoli-contract.md).

---

## Quickstart (60 seconds)

Grab a prebuilt binary and watch the contract fire on the bundled data — no toolchain, no build:

```sh
curl -fsSL https://raw.githubusercontent.com/logannye/rosalind/main/install.sh | sh
cd rosalind-*/

# Build a portable index, sort the bundled BAM, then declare a budget and honor it.
./rosalind index --reference examples/data/illumina_toy/reference.fa --output ref.idx
./rosalind sort  --input examples/data/illumina_toy/alignments.bam --output sorted.bam
./rosalind plan  --index ref.idx --budget-mb 512                       # FITS?  predicted peak
./rosalind variants --index ref.idx --alignments sorted.bam \
    --memory-budget-mb 512 --enforce -o calls.vcf                      # honors it (exit 3/4)
./rosalind verify --manifest calls.vcf.manifest.json                   # re-checks the receipt
```

You'll see `plan` predict `[FITS]`, `variants --enforce` print `contract: OK — realized peak … within`, and `verify: OK`. Tighten `--budget-mb` to `1` and `variants --enforce` *refuses up front* (exit 3, no VCF). That is the whole differentiator, in one minute. (Releases are cut from tags; if none is published yet, build from source below.)

---

## What it does today

- **Bounded whole-genome germline calling** — `rosalind variants --index` streams a coordinate-sorted BAM over all contigs of a persisted index, calling SNVs to a multi-contig VCF with a working set bounded by coverage. Calls are calibrated and **abstention-aware** (no confident call → no row, rather than a guess). Detection accuracy is measured against simulated diploid truth (precision/recall 1.00/1.00 on clean data — see [Accuracy](#accuracy)).
- **A reproducible ML feature substrate** — `rosalind features --index` emits the same bounded pileup stream as a per-locus feature **table** (TSV) instead of variant calls — byte-identical run-to-run, with a hash receipt, for **bit-reproducible** training inputs (see below).
- **Build-once, query-many index** — `rosalind index` builds a portable, memory-mapped FM-index over a (multi-contig) reference; `rosalind locate` answers exact-match queries against it in milliseconds. The index is **never rebuilt on load** and is **byte-identical across builds** of the same reference.
- **Alignment** — `rosalind align` builds a Burrows–Wheeler / FM-index over a reference contig and aligns reads via exact-match seeding, deterministic diagonal chaining, and banded affine-gap refinement. Emits SAM or BGZF-compressed BAM.
- **Streaming I/O** — Reads plain or gzip/bgzf-compressed FASTA/FASTQ, auto-detected from the magic bytes; FASTQ can stream from stdin (`-`).
- **Deterministic coordinate sort** — `rosalind sort`: an external merge sort (spills to disk) that orders a BAM by position within a configurable memory budget.
- **Somatic (tumor/normal) calling** — `rosalind somatic` calls somatic SNVs and simple indels from a paired tumor/normal BAM set using a deterministic binomial log-likelihood-ratio model with explicit depth and allele-fraction filters.
- **Truth-set evaluation** — `rosalind eval-germline` / `rosalind eval-somatic` compare a call set against a truth VCF over confident regions (BED), with variant normalization (left-align + trim) and precision / recall / F1. `eval-germline` is the drop-in interface for a GIAB benchmark.
- **Extensibility** — Build custom bounded per-locus analytics over the `PileupColumn` iterator substrate (see [`examples/custom_pileup_analytics.rs`](examples/custom_pileup_analytics.rs)), inheriting bounded memory + determinism for free. *(The legacy `GenomicPlugin` trait + PyO3 RNA-seq demo still work but are **not** memory-bounded — see [CONTRACT.md](CONTRACT.md).)*
- **Determinism by design** — Primary artifacts are emitted in a canonical, stable order, byte-for-byte identical across repeated runs given identical inputs. See [`docs/determinism.md`](docs/determinism.md).

## Why it matters

Three properties, treated as first-class guarantees rather than nice-to-haves:

1. **Predictable memory.** The bet is that for a growing set of users — sequencing in the field, in the clinic, on a laptop, or at genome scale on modest hardware — *"it fits, and I knew it would"* matters more than raw throughput. Rosalind makes the streaming working set bounded by local coverage and surfaces the realized peak so the bound is **verifiable, not just claimed**.
2. **Reproducibility.** Byte-identical outputs and a per-run BLAKE3 manifest make results auditable — a hard requirement for clinical and regulated pipelines, and a sanity-saver for everyone else.
3. **Honest uncertainty.** Calibrated, abstention-aware calling refuses to emit a call where the evidence is insufficient, instead of papering over it.

The contract is real today: `rosalind plan` predicts before you commit, `--enforce` honors the budget, and `rosalind verify` re-checks the receipt (see [CONTRACT.md](CONTRACT.md)). Under the hood, Rosalind is *also* a research vehicle for **space-bounded genomics**: a `~√t` (square-root-space) evaluation framework as a continuous space/time knob, aimed at **sublinear-space index *construction*** — the future Phase-D direction that would extend the contract to the index build step. That layer is not yet load-bearing; the bounded streaming engine you use today is the practical foundation it builds on.

## Accuracy

The contract wraps a real caller, and its detection accuracy is **measured, not assumed.** Against simulated diploid truth (known het/hom SNVs, reads sampled from both haplotypes, the BAM built directly so this isolates the *caller*, not alignment):

- **40× / 0.5% error (clean):** precision **1.00**, recall **1.00**, F1 **1.00** — every het/hom SNV recovered, zero false positives.
- **12× / 1.5% error (stress):** the caller still finds *every* true variant (unfiltered recall 1.00); the default `min_qual` / `min_depth` PASS filter then trades a little recall for precision under noise (0.90 / 0.68).

Scope is honest: simulated (not yet GIAB), detection-only (position + ref + alt), SNV-only. A real **GIAB HG002** benchmark plugs into the same comparator via `rosalind eval-germline --reference … --calls … --truth … --regions highconf.bed`. Full numbers + the gated harness: [`docs/findings/2026-06-02-germline-accuracy.md`](docs/findings/2026-06-02-germline-accuracy.md).

## Who it's for

- **Edge, field, and low-resource settings** — sequencing on a laptop or portable device where predictable memory matters more than peak throughput.
- **Reproducibility-sensitive work** — pipelines where byte-identical, auditable outputs are a first-class requirement.
- **Builders** — anyone who wants a hackable Rust genomics engine to embed, extend with plugins, or drive from Python.
- **Teaching and learning** — a readable, end-to-end Rust implementation of FM-index alignment, streaming pileup, and variant calling to study, modify, and extend.

## Current scope (what's single-contig vs. whole-genome today)

- **Whole-genome:** germline variant calling via `rosalind variants --index` (all contigs, streaming, bounded) and exact-match lookup via `rosalind index` / `rosalind locate`.
- **Single-contig:** Rosalind's own **aligner** (`rosalind align`) and the FASTA-based `variants --reference` path operate on one reference contig per run. For whole-genome calling, align with any standard aligner and bring the coordinate-sorted BAM to `variants --index`. (Wiring the *aligner* onto the persisted multi-contig index is a later phase — see the roadmap.)
- Variant calling is **single-sample** (germline) or a **tumor/normal pair** (somatic); calling is SNV-focused, with simple indels in the somatic path.
- The engine runs **single-threaded** today. `--memory-budget-mb` is record-only by default; add `--enforce` to honor it (refuse up front / fail loud — see [CONTRACT.md](CONTRACT.md)).

## The memory contract in *your* CI

Drop the Rosalind budget Action (this repo's [`action.yml`](action.yml)) into any pipeline to make a declared memory budget a **gate** — the build fails if a whole-genome calling step would breach it. No toolchain on your runner; the Action fetches a prebuilt binary.

```yaml
- uses: logannye/rosalind@v0.1.0
  with:
    index: ref.idx          # built by `rosalind index`
    alignments: sorted.bam  # coordinate-sorted
    budget-mb: 4096         # refuse up front (exit 3) / fail after (exit 4) on breach
    max-depth: 1000         # optional (default 1000)
    max-read-len: 250       # optional (default 250)
```

It runs `plan` (predicts the peak), then `variants --index --enforce` (honors the budget), and uploads the BLAKE3 receipt as a build artifact. This is the one thing a `--max-mem` flag on another caller can't give you: a portable, declarative, **verifiable** memory budget that fails a stranger's build loudly — the contract, enforced where your pipeline already lives. (Available once a release is published; see Quickstart.)

## Pack a fleet: prediction → placement

Because the predicted peak is computed from the **index header alone** — no run, no read I/O, milliseconds — and is a conservative, *additive* upper bound, you can turn the prediction into a scheduling decision: **prove** that N calling jobs co-located on one node will all fit, before launching a single read.

```sh
# jobs.tsv: one job per line — <index>[\t<max-depth>[\t<max-read-len>]]
rosalind pack --jobs jobs.tsv --node-mb 64000          # pack onto 64 GB nodes
rosalind pack --jobs jobs.tsv --node-mb 64000 --nodes 4 # …or refuse (exit 3) if it needs more than 4
```

```
pack: 37 job(s) → 3 node(s) of 64000 MiB — every node proven within capacity
  node 0: 61920 / 64000 MiB  [sampleA.idx, sampleB.idx, …]
  ...
```

`rosalind plan --index ref.idx --budget-mb 64000 --json` emits the same predicted peak as one line of JSON for a workflow engine to read. This is what an emergent-peak caller (GATK, DeepVariant) structurally cannot do: their peak is only known *after* a possible OOM-kill, so every co-location is a gamble. Here, `Packed` is a proof — each node's summed predicted peak is `≤` its capacity, established up front, and the schedule is deterministic.

## A reproducible feature substrate for ML

The same bounded streaming engine that calls variants can emit **per-locus features** instead — one tabular row per callable position, ready for a model:

```sh
rosalind features --index ref.idx --alignments sorted.bam -o features.tsv
# columns: contig, pos, ref, depth, raw_depth, A/C/G/T counts,
#          per-allele fwd/rev strand counts, mean base-qual, mean mapq
```

```python
import pandas as pd
df = pd.read_csv("features.tsv", sep="\t")   # one line; ready for sklearn/PyTorch/JAX
```

Two properties no other pileup gives you together: it is **bounded** (the whole-genome table streams to disk; peak memory tracks coverage, not genome size — a 1 Mbp toy genome's ~983k-row table is produced in ~6 MiB), and it is **byte-identical run-to-run**, with a BLAKE3 receipt over the output. That means **bit-reproducible training inputs**: hash your feature file, and you can prove this quarter's model saw exactly the same data as last quarter's. `features` honors the same `plan`/`--enforce`/`verify` memory contract as `variants`.

A dependency-light Python boundary ([`python/rosalind.py`](python/rosalind.py), stdlib + numpy) loads the table directly, and [`examples/reproducible_features_demo.py`](examples/reproducible_features_demo.py) trains a small model on it and *proves* the inputs are bit-reproducible (two independent extractions → matching receipt hashes → bit-identical trained weights; see [`docs/findings/2026-06-02-reproducible-features-demo.md`](docs/findings/2026-06-02-reproducible-features-demo.md)). *(TSV today; an Arrow/Parquet egress and a zero-copy `pyarrow` in-process binding are the next step.)*

## ColumnKit: implement one trait, inherit the contract

Want your *own* per-locus metric — methylation, a coverage/QC track, custom ML features, a star-allele genotyper? Implement one trait and run it through the driver; you inherit the **same** bounded whole-genome walk, the **same** working-set bound that `plan`/`--enforce` admit, and the **same** verifiable receipt the shipped subcommands enjoy — without re-deriving any of it.

```rust
use rosalind::{ColumnAnalyzer, run_bounded_whole_genome, PileupColumn};
use std::io::{self, Write};

struct CoverageTrack;
impl ColumnAnalyzer for CoverageTrack {
    fn header(&self) -> Option<String> { Some("#contig\tpos\tdepth\n".into()) }
    fn on_column(&mut self, col: &PileupColumn, contig: &str, out: &mut dyn Write) -> io::Result<()> {
        writeln!(out, "{contig}\t{}\t{}", col.locus.pos.0 + 1, col.depth())
    }
}
// run_bounded_whole_genome(&mut CoverageTrack, source, &ref_view, contigs, params, &mut out)
//   → returns the bounded WorkingSet `plan`/`--enforce` reason about.
```

The trait is *welded* to the bounded kernel — `run_bounded_whole_genome` drives the exact same column stream as the shipped `features` egress (which is itself just the first `ColumnAnalyzer`), so the estimator that admits a run provably upper-bounds the realized working set of *your* analyzer too. No other genomics library can offer "implement this trait, inherit a machine-checkable memory budget + a hash-verifiable receipt," because no other library has a memory contract to inherit. See [`examples/columnkit_coverage.rs`](examples/columnkit_coverage.rs).

## Roadmap

The core primitive is a streaming, CIGAR-aware pileup column stream; variant calling and custom plugins consume it. Performance work deliberately *follows* the unique capability — the target user needs "it fits and is predictable" before "it's fastest."

- **Phase A (done):** the streaming pileup engine; calibrated, abstention-aware germline SNV calling; tumor/normal somatic calling; spec-valid VCF; a BLAKE3 reproducibility receipt per run.
- **Phase B (done):** streaming gzip/bgzf input; a multi-contig FM-index over the concatenated genome with `(contig, position)` resolution; a build-once, memory-mapped, byte-reproducible persisted index (`rosalind index`/`locate`); zero-copy reference access from the index; and **bounded whole-genome germline calling over a sorted BAM** (`rosalind variants --index`) with a realized-memory receipt.
- **Phase C (done):** memory as a *verifiable contract* — `rosalind plan` (a checkable envelope before you commit), `--enforce` (honor-or-refuse: refuse up front / fail loud, never a silent OOM-kill), and `rosalind verify`. See [CONTRACT.md](CONTRACT.md).
- **Hardening & reach (done):** unbiased depth-cap downsampling (no silent variant drops) and a CI-enforced memory gate; **measured** germline detection accuracy ([Accuracy](#accuracy)); the **`rosalind features`** reproducible ML feature substrate; and a one-command adoption on-ramp — prebuilt binaries (`install.sh`, with checksum verification) plus the **Rosalind budget GitHub Action** (`action.yml`, used as `logannye/rosalind@v0.1.0`) that enforces the contract in *your* CI.
- **Fleet scheduling (done):** [prediction → placement](#pack-a-fleet-prediction--placement) — `rosalind pack` proves a co-location of N calling jobs fits a node before launching a byte (predicted peaks are additive and read from the index header); `plan --index --json` for a scheduler to read.
- **ColumnKit SDK (done):** [implement one trait, inherit the contract](#columnkit-implement-one-trait-inherit-the-contract) — a `ColumnAnalyzer` trait + `run_bounded_whole_genome` driver so a builder's own per-locus analyzer inherits bounded memory, determinism, and a verifiable receipt. The shipped `features` egress is the first impl.
- **Phase D (research):** sublinear-space index construction — the `~√t` space/time knob across the full curve — extending the contract to the index *build* step (today's build is O(reference)). The headline space-complexity bet; see [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).
- **Later:** the aligner over the persisted multi-contig index (`align --index`, whole-genome alignment); germline indels and richer read QC; deterministic multithreading; a Python/tensor binding over the pileup stream.

Target architecture and per-phase specs/plans live in [`docs/superpowers/specs/`](docs/superpowers/specs/) and [`docs/superpowers/plans/`](docs/superpowers/plans/); the guiding thesis is in [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).

---

## Install & build

### Prerequisites
- Rust 1.72+ (`rustup` recommended)
- Native compression headers for BAM I/O: `libbz2-dev` & `liblzma-dev` on Debian/Ubuntu, `brew install bzip2 xz` on macOS
- Python 3.9+ with `numpy` for the feature-substrate boundary (`python/rosalind.py`); the optional legacy PyO3 bindings additionally need `maturin` (set `PYO3_PYTHON=/path/to/python` if the default interpreter is unsuitable)

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

# Predict the peak before committing; then call all contigs, honoring the budget.
rosalind plan --index genome.idx --max-depth 1000 --budget-mb 4096
rosalind variants \
  --index genome.idx \
  --alignments sample.sorted.bam \
  --mapq-threshold 20 \
  --memory-budget-mb 4096 --enforce \
  -o sample.vcf
rosalind verify --manifest sample.vcf.manifest.json
```

`variants --index` requires a **coordinate-sorted BAM** (use `rosalind sort` or `samtools sort`). It reads the reference from the index — no `--reference` FASTA needed — and writes a multi-contig VCF plus a memory + reproducibility receipt. With `--enforce` the declared budget is honored (refuse up front / fail loud); without it, it is record-only. The full contract: [CONTRACT.md](CONTRACT.md).

### Try the contract end-to-end (bundled data, in-house tools only)

```bash
D=examples/data/illumina_toy
rosalind index  --reference $D/reference.fa --output /tmp/toy.idx
rosalind sort   --input $D/alignments.bam   --output /tmp/toy.sorted.bam
rosalind plan   --index /tmp/toy.idx --budget-mb 512
rosalind variants --index /tmp/toy.idx --alignments /tmp/toy.sorted.bam \
  --memory-budget-mb 512 --enforce -o /tmp/toy.vcf
rosalind verify --manifest /tmp/toy.vcf.manifest.json
```

This bundled demo is **single-contig** because Rosalind's own aligner is single-contig. For **whole-genome** calling, align with bwa-mem2/minimap2, coordinate-sort, and bring the BAM to `variants --index` — which calls *every* contig in bounded memory.

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

- `rosalind plan` — predict a job's peak memory vs a declared budget *before* committing (`--index` for the variants peak, `--reference` for the index build).
- `rosalind verify` — re-check a reproducibility receipt without re-running: re-hash its inputs/outputs and confirm the realized peak landed within budget.
- `rosalind features --index genome.idx --alignments sorted.bam -o features.tsv` — stream a bounded, byte-identical per-locus feature table (TSV) under the same memory contract (see [the ML substrate](#a-reproducible-feature-substrate-for-ml)).
- `rosalind locate --index genome.idx --pattern GATTACA` — exact-match positions in a prebuilt index (memory-mapped, never rebuilt). Exact-match only; seed/chain/extend alignment against the persisted index is a later phase.
- `rosalind sort` — deterministic coordinate sort of a BAM within a memory budget.
- `rosalind somatic` — tumor/normal somatic SNV + simple-indel calling from a paired BAM set over a region.
- `rosalind eval-germline` / `rosalind eval-somatic` — compare a germline / somatic call set to a truth VCF over confident regions (BED), reporting precision / recall / F1.

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

The dependency-light boundary ([`python/rosalind.py`](python/rosalind.py), stdlib + numpy, no build) runs `rosalind features` and loads the per-locus table into numpy:

```python
import sys; sys.path.insert(0, "python")
from rosalind import features

ft = features("genome.idx", "sample.sorted.bam", binary="rosalind")
X = ft.data            # (n_loci, n_numeric_features) float64, ready for a model
print(len(ft), "loci,", X.shape[1], "features; receipt:", ft.manifest_path)
```

The binary streams the whole-genome table in bounded memory; this just loads the result. [`examples/reproducible_features_demo.py`](examples/reproducible_features_demo.py) is a runnable demo that trains a model and *proves* the inputs are bit-reproducible. A zero-copy in-process `pyarrow` binding is the next step.

> **Legacy.** A PyO3 extension (`maturin develop --release` → `rosalind_py.PyGenomicEngine`) exposes an older plugin demo; it is **not** memory-bounded and does not use the feature substrate. Prefer `python/rosalind.py` above.

---

## Extend

Rosalind's kernel is a **bounded, deterministic `PileupColumn` stream** — build your own per-locus analytics (coverage, QC, methylation, ML features) over it and inherit bounded memory + determinism for free:

- **Rust (recommended)** — consume the `PileupEngine` iterator over any `ReadSource`. See [`examples/custom_pileup_analytics.rs`](examples/custom_pileup_analytics.rs) (`cargo run --example custom_pileup_analytics`) for a non-caller consumer computing per-locus coverage. The contract verbs and the substrate are re-exported at the crate root (`use rosalind::{PileupEngine, PileupColumn, ReadSource, …}`).
- **CLI subcommands** — add workflows in `src/main.rs`; compose subcommands over pipes.

> **Legacy / non-bounded.** The `GenomicPlugin` trait (`src/plugin/`), the `framework/` evaluator, and the Python `rosalind_py.PyGenomicEngine` RNA-seq demo still work but do **not** inherit the memory contract. Prefer the `PileupColumn` substrate above for bounded work. See [CONTRACT.md](CONTRACT.md).

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
