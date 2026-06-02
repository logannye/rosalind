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
