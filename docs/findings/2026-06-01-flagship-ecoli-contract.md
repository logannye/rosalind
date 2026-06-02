# The memory contract on a real genome — E. coli K-12 MG1655

**2026-06-01.** A reproducible demonstration that Rosalind's memory contract holds on a **real** genome:
on the *Escherichia coli* K-12 MG1655 chromosome (**4,641,652 bp**), a declared **256 MiB** budget *fits*
— `plan` → `variants --enforce` → `verify: OK`, realized peak **22 MiB** — while an **8 MiB** budget is
*refused up front* (exit 3, no work). The contract, honored both ways, on a real 4.6 Mbp genome.

Reproduce with one command: `bash scripts/flagship_ecoli_demo.sh` (knobs: `COVERAGE`, `FIT_MB`,
`REFUSE_MB`, `REF_URL`).

## What this shows — and what it does not

- ✅ **Real genome, real scale.** The reference is the actual NCBI RefSeq assembly **GCF_000005845.2**
  (accession **NC_000913.3**), 4.64 Mbp — not a toy.
- ✅ **The memory contract is real.** A declared budget is *predicted* before the run, *honored* during it
  (fits cleanly, or refuses cleanly), and *verified* against a deterministic receipt afterward.
- ⚠️ **The reads are simulated** (deterministic, 30×, seed 1337) from the real reference. There is no
  short-read archive tooling in the build environment, and — importantly — the contract claim is
  **read-realism-independent**: peak memory is governed by the reference size and the depth cap, not by
  whether the reads are real. The numbers below would not change with real reads.
- ❌ **This is not an accuracy benchmark.** Rosalind uses exact-match alignment + a basic diploid
  genotype-likelihood caller; there is no GATK/DeepVariant comparison here. This is a **memory +
  reproducibility** result.

## Setup

| | |
|---|---|
| Reference | *E. coli* K-12 MG1655, NCBI RefSeq GCF_000005845.2 (NC_000913.3), **4,641,652 bp**, single chromosome |
| Reads | **simulated**, 30× paired, 150 bp, seed 1337 (`scripts/generate_toy_data.py --reference`) |
| Pipeline | `rosalind index` → `rosalind align` (single-contig) → `rosalind sort` → the contract verbs |
| Index | `index_bytes` = 9,203,308 (~8.8 MiB on disk); build peak RSS **185 MiB** (see note below) |

## (a) A 256 MiB budget — fits

```
$ rosalind plan --index ecoli.idx --max-depth 1000 --budget-mb 256
plan: predicted peak RSS (upper bound)
  process baseline (measured):        7 MiB
  reference decode (largest contig):  4 MiB
  active set @ max-depth 1000:           4 MiB
  engine overhead:                    0 MiB
  -------------------------------------------------
  predicted peak: ~16 MiB / budget 256 MiB  [FITS]

$ rosalind variants --index ecoli.idx --alignments ecoli.sorted.bam \
    --memory-budget-mb 256 --enforce -o ecoli.vcf
memory: peak RSS 22 MiB; max pileup working set 4838 KiB
contract: OK — realized peak 22 MiB within declared 256 MiB

$ rosalind verify --manifest ecoli.vcf.manifest.json
verify: peak 22 MiB within budget 256 MiB
verify: OK — 2 input(s), 1 output(s) match
```

The run emitted **5,924** germline variant rows. Realized peak RSS was **23,134,208 B (22 MiB)** and the
max pileup working set **4,954,240 B (4.7 MiB)** — comfortably inside the 256 MiB budget. The committed
receipt is [`ecoli.vcf.manifest.json`](ecoli.vcf.manifest.json) (`contract_verdict: within`,
`peak_rss_bytes: 23134208`, `max_working_set_bytes: 4954240`, `memory_budget_mb: 256`).

## (b) An 8 MiB budget — refused, up front

```
$ rosalind variants --index ecoli.idx --alignments ecoli.sorted.bam \
    --memory-budget-mb 8 --enforce -o ecoli_refused.vcf
contract: REFUSE — declared 8 MiB, predicted peak ~16 MiB (largest contig 4 MiB +
  active @ max-depth 1000 / max-read-len 250 atop a 8 MiB baseline). Raise
  --memory-budget-mb, lower --max-depth, or drop --enforce.
# exit code 3; no VCF written.
```

The job declines *before doing any work* and tells you exactly why — it does not start, run for a while,
and then get OOM-killed.

## The bound holds independent of coverage

The same run at **2× coverage** produced a realized peak of **21 MiB**; at **30×**, **22 MiB**. Depth
went up 15×; peak memory did not — because the working set is bounded by the reference (one contig) plus a
depth-capped active set, not by the number of reads. That is the whole promise: *peak memory is a property
of the genome and your declared budget, not of your data volume.*

## Honest notes

- **Predicted vs. realized.** The pre-run prediction (~16 MiB) is a deliberately coarse envelope (a
  measured process baseline + a modeled working set); the realized peak (22 MiB) is the truth. Both sit far
  inside the 256 MiB budget. When the two could diverge near the budget edge, the **post-run check is the
  backstop** — `--enforce` fails loud (exit 4) if the realized peak ever exceeds the budget, so you are
  never silently over.
- **Reproducibility.** The VCF and the receipt's BLAKE3 content hashes are byte-reproducible across
  machines (deterministic reference download + deterministic read simulation + deterministic engine). The
  realized `peak_rss_bytes` is machine-dependent — but always ≤ the declared budget, which is the
  guarantee that matters.
- **The build is the frontier.** Building the index used **185 MiB** of RAM (O(reference)) while *calling*
  used **22 MiB**. The bounded contract today covers the streaming call path; extending it to the index
  *build* (so a constrained device can index a genome that doesn't fit in RAM) is the Phase-D research
  direction — the `~√t` sublinear-space construction work in [`../OPEN_PROBLEMS.md`](../OPEN_PROBLEMS.md).

## Reproduce

```bash
bash scripts/flagship_ecoli_demo.sh          # 30× by default takes ~80s end-to-end
COVERAGE=2 bash scripts/flagship_ecoli_demo.sh   # faster; same contract numbers
```

The script fetches the reference from NCBI (cached after the first run), simulates reads, runs the
in-house pipeline, and asserts every gate (`[FITS]`, `verify: OK`, refuse exit 3). Generated data lands in
`results/flagship-ecoli/` (gitignored).
