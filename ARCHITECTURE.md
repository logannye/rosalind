# Rosalind Architecture

A map of the codebase, so you know where a change belongs. Companion to
[`docs/ROADMAP.md`](docs/ROADMAP.md) (where we're going) and
[`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md) (the research thesis).

## The kernel

Rosalind's spine is a **bounded, deterministic `PileupColumn` stream**. Variant calling, gVCF, the
feature substrate, and any ColumnKit analyzer all consume that one stream, so they inherit the same
memory contract and byte-reproducibility.

## Modules (`src/`)

| Module | Responsibility |
|---|---|
| `core/` | Shared types: `Locus`/`ContigSet`/`AlignedRead`, the memory `budget` + `WorkingSet`, the runtime `governor`, the typed `CoreError`. |
| `io/` | Streaming FASTA/FASTQ (gzip/bgzf auto-detect), BAM/VCF read+write (via rust-htslib), decompression. |
| `genomics/` | The FM-index: `suffix_array` (SA-IS), `fm_index`/`fm_backing`, `rank_select`, the persisted `index/` (format, mmap `view`), `compressed_dna`, the `bwt_aligner`, deterministic `sort`, and `eval/` (truth-set comparison). |
| `pileup/` | The streaming pileup `engine` + `column` + read `source` — the bounded kernel. |
| `call/` | Consumers of the kernel: `germline`/`somatic`/`gvcf` calling, the `features` egress, the `columnkit` SDK, the `plan` estimator, fleet `pack`, the `pipeline`/`whole_genome` drivers. |
| `provenance/` | The canonical-JSON, self-hashing BLAKE3 run receipt (`RunManifest`). |
| `util/` | Process RSS (`getrusage`), mmap helpers. |

## Roadmap phases → code

- **A** (streaming pileup + calling): `pileup/`, `call/{germline,somatic,gvcf,pipeline}`.
- **B** (genome-scale index): `genomics/{suffix_array,fm_index,index}`, `io/bam`.
- **C** (the memory contract): `core/{budget,governor}`, `call/plan`, `provenance/`, `verify` in `main.rs`.
- **D** (sublinear-space index *build* — research): `genomics/suffix_array` + a future external-memory constructor. See `docs/OPEN_PROBLEMS.md`.
- **E** (reach + ML substrate): `call/{features,columnkit}`, `python/`, `genomics/eval`.

## CLI

`src/main.rs` is the CLI surface (`index`, `align`, `sort`, `variants`, `somatic`, `features`, `plan`,
`pack`, `verify`, `locate`, `eval-*`), thin over the library.
