#!/usr/bin/env bash
# D0 measure-first probe: build the index for a small and a ~3x-larger REAL genome,
# capture each build receipt (realized peak RSS + modeled n-scale memory + bytes/base),
# and surface the pre-registered gate verdict. Build-only — no reads/alignment.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS="$ROOT/results/build-probe"
mkdir -p "$RESULTS"
BIN="$ROOT/target/release/rosalind"
[ -x "$BIN" ] || (cd "$ROOT" && cargo build --release)

# E. coli (gate genome) — reuse the Move-#5 cached reference if present.
ECOLI="$ROOT/results/flagship-ecoli/ecoli.fa"
if [ ! -s "$ECOLI" ]; then
  ECOLI="$RESULTS/ecoli.fa"
  if [ ! -s "$ECOLI" ]; then
    curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" -o "$RESULTS/ecoli.fa.gz"
    gunzip -f "$RESULTS/ecoli.fa.gz"
  fi
fi

# S. cerevisiae R64 (~12.1 Mbp, ~3x; multi-contig — build handles it).
YEAST="$RESULTS/yeast.fa"
if [ ! -s "$YEAST" ]; then
  curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz" -o "$RESULTS/yeast.fa.gz"
  gunzip -f "$RESULTS/yeast.fa.gz"
fi

SUMMARY="$RESULTS/SUMMARY.txt"
: > "$SUMMARY"
probe() {
  local name="$1" ref="$2"
  echo "================ $name ================" | tee -a "$SUMMARY"
  "$BIN" index --reference "$ref" --output "$RESULTS/$name.idx" 2>&1 | tee -a "$SUMMARY"
  echo | tee -a "$SUMMARY"
}
probe "ecoli" "$ECOLI"
probe "yeast" "$YEAST"

echo ">> D0 gate (assess in the findings doc): realized peak within +/-25% of ~185 MiB on E. coli" | tee -a "$SUMMARY"
echo ">> AND >=70% model/realized attribution; bytes/base ~constant across E. coli & yeast = n-scale." | tee -a "$SUMMARY"
echo ">> Summary: $SUMMARY"
