# Move #5 — Flagship Artifact (the contract on a real genome) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A reproducible, committed proof that Rosalind's memory contract holds on the real E. coli K-12 MG1655 genome — fits a sane budget and refuses an impossible one — as the repo's headline "see it work" artifact.

**Architecture:** Extend `scripts/generate_toy_data.py` to simulate reads from a *given* FASTA; add a `scripts/flagship_ecoli_demo.sh` orchestrator that fetches the real reference from NCBI, simulates reads, runs `index → align → sort`, then demonstrates the contract both ways (`plan`/`variants --enforce`/`verify` fits at 256 MiB; refuses at 8 MiB). Commit a `docs/findings/` writeup + the fit receipt + a README pin; gitignore the big regenerated data. No engine code changes.

**Tech Stack:** bash, `curl`, `python3` (all present); the release `rosalind` binary. No external aligner/samtools/SRA. MSRV unaffected; no new Rust deps.

**Spec:** [`docs/superpowers/specs/2026-06-01-flagship-ecoli-artifact-design.md`](../specs/2026-06-01-flagship-ecoli-artifact-design.md).

**Note on numbers:** the writeup (Task 4) transcribes the **actual captured numbers** from the Task-3 run — they are not known until the run completes (that is the nature of a findings artifact), and are not "placeholders" in the plan-failure sense.

---

## File Structure

- **Modify** `scripts/generate_toy_data.py` — add a `--reference <fa>` option (simulate reads from a given FASTA instead of inventing one).
- **Create** `scripts/flagship_ecoli_demo.sh` — the reproducible orchestrator.
- **Modify** `.gitignore` — ignore the regenerated `results/` data.
- **Create** `docs/findings/2026-06-01-flagship-ecoli-contract.md` — the writeup (with real captured numbers).
- **Create** `docs/findings/ecoli.vcf.manifest.json` — the committed fit receipt.
- **Modify** `README.md` — a short "Proof: the contract on a real genome" pin.

---

## Task 1: Extend `generate_toy_data.py` to simulate reads from a given FASTA

**Files:**
- Modify: `scripts/generate_toy_data.py`

- [ ] **Step 1: Add the import + a FASTA reader.** Change the typing import and add a helper. Replace:

```python
from typing import Tuple
```
with:
```python
from typing import Optional, Tuple
```

Add this helper (after `revcomp`):

```python
def read_fasta_sequence(path: Path) -> str:
    """Concatenate the FIRST FASTA record's sequence, uppercased (single-contig use)."""
    parts: list[str] = []
    started = False
    with path.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if started:
                    break  # only the first record (the chromosome)
                started = True
                continue
            parts.append(line.strip().upper())
    return "".join(parts)
```

- [ ] **Step 2: Use the given reference in `generate_dataset`.** Replace the signature + the reference-construction block:

```python
def generate_dataset(output_dir: Path, ref_length: int, coverage: int, seed: int) -> Tuple[Path, Path, Path]:
    rng = random.Random(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference = random_dna(ref_length, rng)
    reference_path = output_dir / "reference.fa"
    with reference_path.open("w", encoding="ascii") as handle:
        handle.write(">chrToy\n")
        for i in range(0, len(reference), 80):
            handle.write(reference[i : i + 80] + "\n")

    total_reads = max(1, int((ref_length * coverage) / READ_LENGTH))
```
with:
```python
def generate_dataset(
    output_dir: Path,
    ref_length: int,
    coverage: int,
    seed: int,
    reference_fa: Optional[Path] = None,
) -> Tuple[Path, Path, Path]:
    rng = random.Random(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    if reference_fa is not None:
        # Simulate reads FROM a given (e.g. real, downloaded) reference; do not
        # invent or overwrite one.
        reference = read_fasta_sequence(reference_fa)
        reference_path = reference_fa
    else:
        reference = random_dna(ref_length, rng)
        reference_path = output_dir / "reference.fa"
        with reference_path.open("w", encoding="ascii") as handle:
            handle.write(">chrToy\n")
            for i in range(0, len(reference), 80):
                handle.write(reference[i : i + 80] + "\n")

    total_reads = max(1, int((len(reference) * coverage) / READ_LENGTH))
```

(The only other change is `ref_length` → `len(reference)` in the `total_reads` line, so coverage is computed from the actual reference whether random or provided.)

- [ ] **Step 3: Wire the CLI arg.** In `main()`, after the `--seed` argument, add:

```python
    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="Simulate reads from this FASTA's sequence (do not invent/overwrite a reference)",
    )
```

And change the call:
```python
    reference_path, r1_path, r2_path = generate_dataset(args.output, args.length, args.coverage, args.seed)
```
to:
```python
    reference_path, r1_path, r2_path = generate_dataset(
        args.output, args.length, args.coverage, args.seed, args.reference
    )
```

- [ ] **Step 4: Test it on a tiny FASTA.**

Run:
```bash
cd ~/rosalind && T=$(mktemp -d) && printf '>tinychr\n%s\n' "$(python3 -c "import random;r=random.Random(1);print(''.join(r.choice('ACGT') for _ in range(600)))")" > "$T/ref.fa" && python3 scripts/generate_toy_data.py "$T/out" --reference "$T/ref.fa" --coverage 4 --seed 7 && echo "R1 reads: $(($(wc -l < "$T/out/reads_R1.fastq")/4))" && head -2 "$T/out/reads_R1.fastq" && rm -rf "$T"
```
Expected: a non-zero R1 read count and a 150-or-less-bp read line; no traceback. (The reference is 600 bp so reads fit.)

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add scripts/generate_toy_data.py && git commit -m "feat(scripts): generate_toy_data --reference simulates reads from a given FASTA (Move #5)"
```

---

## Task 2: The reproducible orchestrator + gitignore

**Files:**
- Create: `scripts/flagship_ecoli_demo.sh`
- Modify: `.gitignore`

- [ ] **Step 1: Write the script.** Create `scripts/flagship_ecoli_demo.sh`:

```bash
#!/usr/bin/env bash
# Flagship demo: Rosalind's memory contract on the REAL E. coli K-12 MG1655 genome.
#
# Reads are SIMULATED from the real reference (deterministic). The claim is the
# memory CONTRACT (plan -> enforce -> verify; bounded peak), NOT calling accuracy.
#
# Knobs (env): REF_URL, COVERAGE (default 5), FIT_MB (256), REFUSE_MB (8), SEED (1337).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS="$ROOT/results/flagship-ecoli"
mkdir -p "$RESULTS"

REF_URL="${REF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz}"
COVERAGE="${COVERAGE:-5}"
FIT_MB="${FIT_MB:-256}"
REFUSE_MB="${REFUSE_MB:-8}"
SEED="${SEED:-1337}"

BIN="$ROOT/target/release/rosalind"
if [ ! -x "$BIN" ]; then
  echo ">> building release binary..."
  (cd "$ROOT" && cargo build --release)
fi

REF="$RESULTS/ecoli.fa"
if [ ! -s "$REF" ]; then
  echo ">> fetching real E. coli K-12 MG1655 reference from NCBI"
  curl -fsSL "$REF_URL" -o "$RESULTS/ecoli.fa.gz"
  gunzip -f "$RESULTS/ecoli.fa.gz"
  # Keep only the first record (the chromosome) so the single-contig aligner applies.
  python3 - "$RESULTS/ecoli.fa" > "$RESULTS/ecoli.first.fa" <<'PY'
import sys
recs = 0
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith(">"):
            recs += 1
            if recs > 1:
                break
        print(line, end="")
PY
  mv "$RESULTS/ecoli.first.fa" "$REF"
fi
BP=$(grep -v '^>' "$REF" | tr -d '\n' | wc -c | tr -d ' ')
echo ">> reference: $REF ($BP bp)"

echo ">> simulating ${COVERAGE}x reads from the real reference (seed $SEED)"
python3 "$ROOT/scripts/generate_toy_data.py" "$RESULTS" --reference "$REF" --coverage "$COVERAGE" --seed "$SEED"

echo ">> index -> align -> sort (in-house, single-contig)"
"$BIN" index --reference "$REF" --output "$RESULTS/ecoli.idx"
"$BIN" align --reference "$REF" --reads "$RESULTS/reads_R1.fastq" --format bam --output "$RESULTS/raw.bam"
"$BIN" sort --input "$RESULTS/raw.bam" --output "$RESULTS/ecoli.sorted.bam"

SUMMARY="$RESULTS/SUMMARY.txt"
: > "$SUMMARY"
echo "== Rosalind flagship: memory contract on real E. coli K-12 MG1655 ($BP bp) ==" | tee -a "$SUMMARY"
echo "(real reference; reads SIMULATED ${COVERAGE}x seed $SEED; claim = memory contract, not accuracy)" | tee -a "$SUMMARY"
echo | tee -a "$SUMMARY"

echo ">> (a) FITS: plan -> variants --enforce (--memory-budget-mb $FIT_MB) -> verify"
echo "--- plan (budget ${FIT_MB} MiB) ---" | tee -a "$SUMMARY"
"$BIN" plan --index "$RESULTS/ecoli.idx" --max-depth 1000 --budget-mb "$FIT_MB" 2>&1 | tee -a "$SUMMARY"
echo "--- variants --enforce (budget ${FIT_MB} MiB) ---" | tee -a "$SUMMARY"
"$BIN" variants --index "$RESULTS/ecoli.idx" --alignments "$RESULTS/ecoli.sorted.bam" \
  --memory-budget-mb "$FIT_MB" --enforce -o "$RESULTS/ecoli.vcf" 2>&1 | tee -a "$SUMMARY"
echo "--- verify ---" | tee -a "$SUMMARY"
"$BIN" verify --manifest "$RESULTS/ecoli.vcf.manifest.json" 2>&1 | tee -a "$SUMMARY"
echo | tee -a "$SUMMARY"

echo ">> (b) REFUSES: variants --enforce (--memory-budget-mb $REFUSE_MB) — expect exit 3"
echo "--- variants --enforce (budget ${REFUSE_MB} MiB) ---" | tee -a "$SUMMARY"
set +e
"$BIN" variants --index "$RESULTS/ecoli.idx" --alignments "$RESULTS/ecoli.sorted.bam" \
  --memory-budget-mb "$REFUSE_MB" --enforce -o "$RESULTS/ecoli_refused.vcf" 2>&1 | tee -a "$SUMMARY"
RC=${PIPESTATUS[0]}
set -e
echo "refuse exit code: $RC (expected 3)" | tee -a "$SUMMARY"
if [ "$RC" != "3" ]; then echo "FAIL: expected refuse exit 3, got $RC"; exit 1; fi
if [ -f "$RESULTS/ecoli_refused.vcf" ]; then echo "FAIL: refuse must not write a VCF"; exit 1; fi

echo ">> DONE. Summary: $SUMMARY ; receipt: $RESULTS/ecoli.vcf.manifest.json"
```

- [ ] **Step 2: Make it executable + gitignore the regenerated data.**

```bash
cd ~/rosalind && chmod +x scripts/flagship_ecoli_demo.sh
```

Append to `.gitignore`:
```
# Flagship demo: regenerated by scripts/flagship_ecoli_demo.sh (not committed).
/results/
```

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add scripts/flagship_ecoli_demo.sh .gitignore && git commit -m "feat(scripts): flagship_ecoli_demo.sh — the contract on the real E. coli genome (Move #5)"
```

---

## Task 3: Run the demo (produce the real numbers)

**Files:** none (produces `results/flagship-ecoli/`, gitignored)

- [ ] **Step 1: Run it.** (Use a generous timeout — the in-house aligner is unoptimized; the build + download + align dominate. If alignment is impractically slow, re-run with a lower coverage, e.g. `COVERAGE=2`, then `COVERAGE=1`; the contract behavior is coverage-independent given the depth cap. If the NCBI URL fails, the spec's fallback is a smaller real genome via `REF_URL=`.)

Run: `cd ~/rosalind && COVERAGE=5 bash scripts/flagship_ecoli_demo.sh`
Expected: the script completes (exit 0): `plan … [FITS]`, `variants` completes + writes a receipt, `verify: OK`, and the refuse run exits 3 with no VCF. The asserts inside the script gate all of this.

- [ ] **Step 2: Read the captured numbers** for the writeup.

Run: `cd ~/rosalind && cat results/flagship-ecoli/SUMMARY.txt`
Expected: the full transcript — note the reference bp, the plan predicted peak + `[FITS]`, the realized `memory: peak RSS … MiB; max pileup working set … KiB`, the `contract: OK …`, the `verify: OK …`, and the refuse message + `predicted peak ~… MiB`. These are the numbers for Task 4.

(No commit — `results/` is gitignored. The captured numbers flow into Task 4.)

---

## Task 4: The findings writeup + the committed receipt

**Files:**
- Create: `docs/findings/2026-06-01-flagship-ecoli-contract.md`
- Create: `docs/findings/ecoli.vcf.manifest.json`

- [ ] **Step 1: Copy the fit receipt into `docs/findings/`** (the committed proof object).

Run: `cd ~/rosalind && mkdir -p docs/findings && cp results/flagship-ecoli/ecoli.vcf.manifest.json docs/findings/ecoli.vcf.manifest.json`

- [ ] **Step 2: Write `docs/findings/2026-06-01-flagship-ecoli-contract.md`** with these sections, filling the **actual captured numbers** from Task 3's `SUMMARY.txt` (no invented numbers):

  1. **Title + one-line result:** "Rosalind's memory contract on the real E. coli K-12 MG1655 genome (4,641,652 bp): declared 256 MiB, realized peak `<X>` MiB, `verify: OK`; an 8 MiB budget is refused up front (exit 3)." (Use the real bp from `SUMMARY.txt`.)
  2. **What this shows / what it does not** — the §1 claim from the spec, verbatim in spirit: real reference; **reads simulated** (deterministic, `<COVERAGE>×`, seed 1337); the claim is the **memory contract**, not calling accuracy; no GATK comparison.
  3. **Side (a) — fits:** the `plan` breakdown (predicted peak components), the realized `peak RSS` + `max pileup working set`, the `contract: OK` line, the `verify: OK` line. Paste the relevant `SUMMARY.txt` lines in a fenced block.
  4. **Side (b) — refuses:** the 8 MiB `variants --enforce` refuse message (predicted peak vs budget) + `exit 3` + "no VCF written."
  5. **Reproduce:** `bash scripts/flagship_ecoli_demo.sh` (knobs: `COVERAGE`, `FIT_MB`, `REFUSE_MB`, `REF_URL`). Note: the VCF + the receipt's BLAKE3 content hashes are byte-reproducible across machines; the realized `peak_rss_bytes` is machine-dependent but always ≤ the budget.
  6. **The committed receipt:** link `ecoli.vcf.manifest.json` and note the key fields (`contract_verdict=within`, `peak_rss_bytes`, `max_working_set_bytes`, `memory_budget_mb=256`).

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add docs/findings/2026-06-01-flagship-ecoli-contract.md docs/findings/ecoli.vcf.manifest.json && git commit -m "docs(findings): flagship — the memory contract on the real E. coli genome (Move #5)"
```

---

## Task 5: README pin

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add a "Proof" block.** Immediately after the headline command block's closing (the section that ends with the `rosalind verify` example, before the `---` separator), insert:

```markdown

**Proof — the contract on a real genome.** On the real *E. coli* K-12 MG1655 chromosome (4,641,652 bp), a declared **256 MiB** budget *fits* (`plan` → `variants --enforce` → `verify: OK`, realized peak well inside the budget), while an **8 MiB** budget is *refused up front* (exit 3, no work) — the contract honored both ways. Reads are simulated from the real reference; the claim is the memory contract, not calling accuracy. Full numbers + one-command reproduction: [`docs/findings/2026-06-01-flagship-ecoli-contract.md`](docs/findings/2026-06-01-flagship-ecoli-contract.md).
```

(Verify the exact realized-peak phrasing against the writeup; keep it qualitative in the README — "well inside the budget" — and let the findings doc carry the exact figure.)

- [ ] **Step 2: Commit**

```bash
cd ~/rosalind && git add README.md && git commit -m "docs(readme): pin the real-genome contract proof (Move #5)"
```

---

## Task 6: Final verification

**Files:** none

- [ ] **Step 1: Confirm the artifact is coherent + the tree is clean of regenerated data.**

Run: `cd ~/rosalind && git status --short && echo "---" && ls docs/findings/ && echo "--- results gitignored? ---" && git check-ignore results/flagship-ecoli/ecoli.fa && echo "(ignored ✓)"`
Expected: clean tree (no `results/` tracked); `docs/findings/` has the writeup + the receipt; `results/...` is ignored.

- [ ] **Step 2: Confirm the repo still builds + tests pass** (no engine code changed, but verify nothing regressed).

Run: `cd ~/rosalind && cargo build 2>&1 | tail -2 && cargo test 2>&1 | grep -E "FAILED|[1-9][0-9]* failed" || echo "no failures"`
Expected: builds; no failures.

---

## Self-Review notes

- **Spec coverage:** §2 script → Tasks 1 (`generate_toy_data --reference`) + 2 (`flagship_ecoli_demo.sh`); §2 run → Task 3; §3 artifact (writeup + receipt + gitignore) → Tasks 2 (gitignore) + 4; §4 README pin → Task 5; §5 verification (script self-asserts) → Task 3 Step 1 + Task 6.
- **Placeholder scan:** the only `<X>`/`<COVERAGE>` tokens are real numbers captured at run time for the findings doc (Task 4) — explicitly sourced from `SUMMARY.txt`, not plan placeholders. All code steps (the Python change, the bash script) are complete.
- **Consistency:** the script's paths (`results/flagship-ecoli/ecoli.{fa,idx,sorted.bam,vcf,vcf.manifest.json}`) are used identically in Tasks 2/3/4; `--reference` (Task 1) is the flag the script passes (Task 2); the 256/8 MiB budgets and the exit-3 refuse match the spec.
- **Risk:** the aligner-speed unknown is handled by the `COVERAGE` knob (Task 3 Step 1) and the `REF_URL` fallback — the contract step is coverage-independent.
