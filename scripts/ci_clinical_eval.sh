#!/bin/bash
#
# CI-friendly clinical evaluation smoke run.
#
# Generates a tiny paired-end toy dataset, runs end-to-end somatic twice, and checks determinism.
#
set -euo pipefail

WORKDIR="${WORKDIR:-/tmp/rosalind_ci_eval}"
DATA="${WORKDIR}/data"
BUNDLE="${WORKDIR}/bundle"
OUT="${WORKDIR}/out.vcf"

mkdir -p "${WORKDIR}"
rm -rf "${DATA}" "${BUNDLE}"
mkdir -p "${DATA}" "${BUNDLE}"

echo "[1/3] Generate toy paired-end data"
python scripts/generate_toy_data.py "${DATA}" --length 20000 --coverage 2 --seed 1337

echo "[2/3] Build a minimal bundle (chrToy)"
python scripts/prepare_eval_bundle.py \
  --reference "${DATA}/reference.fa" \
  --contig chrToy \
  --tumor-r1 "${DATA}/reads_R1.fastq" --tumor-r2 "${DATA}/reads_R2.fastq" \
  --normal-r1 "${DATA}/reads_R1.fastq" --normal-r2 "${DATA}/reads_R2.fastq" \
  --out "${BUNDLE}" \
  --seed 1337 --max-pairs 2000

echo "[3/3] Run end-to-end twice and hash compare"
bash scripts/run_clinical_eval.sh --bundle "${BUNDLE}" --out "${OUT}" --workdir "${WORKDIR}/runs" --memory-mb 256

echo "OK: CI clinical eval smoke passed"


