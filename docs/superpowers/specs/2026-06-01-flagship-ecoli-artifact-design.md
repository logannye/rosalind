# Move #5 — Flagship artifact: the contract on a real genome (design)

**Status:** Spec for review — 2026-06-01. Strategy Move #5 (the 2026-06-01 synthesis): a reproducible,
committed proof that the memory contract holds on a **real genome** — the repo's legible "see it actually
work" headline. Built on a fresh branch off the merged `main` (Phase C + the front door are live).

## 1. Goal + the precise claim

Produce a reproducible artifact showing Rosalind's memory contract — **`plan` → `variants --enforce` →
`verify`** — on the **real E. coli K-12 MG1655 genome** (NCBI RefSeq GCF_000005845.2, a single ~4.64 Mbp
chromosome), demonstrating **both sides** of the contract (it *fits* a sane budget and *refuses* an
impossible one). The artifact is a committed **receipt + writeup**, reproducible via a committed script.

**The claim is the memory contract, stated honestly:**
- ✅ on a **real genome** at real bacterial scale (4.64 Mbp), the run **honors a declared budget** — fits
  cleanly under a sane budget, refuses cleanly under an impossible one — and `verify` confirms the realized
  peak landed inside the budget.
- ⚠️ the **reads are simulated** (deterministic) from the real reference — there is no SRA tooling in this
  environment and the contract claim is read-realism-independent (memory behavior does not depend on read
  realism). The writeup states this plainly.
- ❌ **NOT** a calling-accuracy claim. No GATK/DeepVariant comparison (Rosalind is exact-match alignment +
  a basic diploid GL caller). The artifact is a *memory + reproducibility* proof, not an accuracy benchmark.

## 2. The reproducible script — `scripts/flagship_ecoli_demo.sh`

Bash (curl + python3 + the built `rosalind` binary; no external aligner/samtools/SRA needed). Idempotent;
caches downloads; writes everything under a `results/` dir (gitignored). Steps:

1. **Fetch the real reference** (cache if present):
   `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`
   → gunzip → `results/ecoli.fa`. (GCF_000005845.2 is the chromosome only — single record. If extra
   records are present, the script keeps only the first/chromosome so the single-contig in-house aligner
   applies.)
2. **Simulate deterministic reads** from `results/ecoli.fa` via the existing
   `scripts/generate_toy_data.py`, **extended with a `--reference <fa>` option** that reads a given FASTA's
   sequence (instead of inventing a random one) and emits `reads_R1/R2.fastq` from it. Default coverage
   `10×` via a `COVERAGE` env knob (dial down if alignment is slow).
3. **Build the pipeline (in-house):** `rosalind index --reference results/ecoli.fa --output results/ecoli.idx`
   → `rosalind align --reference results/ecoli.fa --reads <reads> --format bam --output results/raw.bam`
   → `rosalind sort --input results/raw.bam --output results/ecoli.sorted.bam`.
4. **Contract — side (a), fits** (the headline loop), capturing stdout/stderr to `results/`:
   - `rosalind plan --index results/ecoli.idx --max-depth 1000 --budget-mb 256` → expect `[FITS]`.
   - `rosalind variants --index results/ecoli.idx --alignments results/ecoli.sorted.bam --memory-budget-mb 256 --enforce -o results/ecoli.vcf` → completes; receipt at `results/ecoli.vcf.manifest.json`.
   - `rosalind verify --manifest results/ecoli.vcf.manifest.json` → expect `verify: OK`.
5. **Contract — side (b), refuses:**
   - `rosalind variants --index results/ecoli.idx --alignments results/ecoli.sorted.bam --memory-budget-mb 8 --enforce -o results/ecoli_refused.vcf` → expect **exit 3 (REFUSE)**, no VCF written, with the
     actionable message. (Script tolerates the non-zero exit and records it.)
6. **Capture** the plan breakdown, the realized memory line, the verify output, the refuse message, and
   wall-clock timings into `results/SUMMARY.txt` for transcription into the writeup.

Risk handling: the in-house aligner is not optimized; if step 3 is impractical at 4.64 Mbp / 10×, lower
`COVERAGE` (the working-set / contract behavior is unaffected by coverage given the depth cap), or fall
back to a smaller real genome (e.g., phiX174 / a plasmid) — the script's `REF_URL`/`COVERAGE` knobs make
this a one-line change. The **contract step (4–5) is fast regardless** (bounded streaming).

## 3. The committed artifact

- **`docs/findings/2026-06-01-flagship-ecoli-contract.md`** — the narrative + the real numbers (genome +
  accession, read sim params, the predicted-vs-realized peak for the 256 MiB fit, the verify confirmation,
  the 8 MiB refuse with its predicted peak), the exact reproduction commands (`scripts/flagship_ecoli_demo.sh`),
  and the §1 honest caveats. This is the headline doc.
- **`docs/findings/ecoli.vcf.manifest.json`** — the committed receipt from the fit run (small; the proof
  object: `contract_verdict=within`, `peak_rss_bytes`, `max_working_set_bytes`, `memory_budget_mb=256`,
  BLAKE3 input/output hashes). Note in the writeup: the VCF + content hashes are byte-reproducible across
  machines; the realized `peak_rss` is machine-dependent but always ≤ the budget.
- **`.gitignore`** — add `results/` and the downloaded genome / reads / BAM / VCF (do not commit the
  ~4.6 Mbp FASTA or the reads/BAM; the script regenerates them).

## 4. README pin

A short **"Proof: the contract on a real genome"** block (near the headline) linking the findings doc with
the one-line headline numbers (e.g., "real 4.64 Mbp E. coli, declared 256 MiB, realized peak NNN MiB,
`verify: OK`; an 8 MiB budget is refused up front"). Honest one-liner that reads are simulated.

## 5. Verification

The run **is** the verification — the script asserts each gate (plan `[FITS]`, variants exit 0 + receipt
written, `verify: OK`, refuse exit 3) and fails loudly otherwise. No new unit tests (the contract is
already covered by `tests/plan_enforce.rs` + the C1 library tests); this is a real-genome *demonstration*,
not a test. The artifact's numbers are transcribed from the actual captured run.

## 6. Non-goals

No accuracy benchmark / GATK comparison; no real reads (no SRA tooling); no multi-contig real genome (the
in-house aligner is single-contig; the multi-contig flagship is documented as a bring-your-own-aligner
command in the README). No new engine code — this is a script + a writeup + a README pin. The genome and
reads are not committed (reproduced by the script). Phase D (√t construction) is unrelated.

## 7. Branch

`rosalind/flagship-ecoli`, off the merged `main`. Its own PR when done (the user decides merge).
