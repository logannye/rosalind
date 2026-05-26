#!/bin/bash
#
# End-to-end clinical evaluation runner (tumor/normal WGS slice).
#
# This script is intentionally deterministic:
# - Runs the same pipeline twice with the same inputs
# - Compares output hashes (VCF + manifest)
#
# Usage (after preparing a bundle with scripts/prepare_eval_bundle.py):
#   bash scripts/run_clinical_eval.sh --bundle /path/to/bundle --out /tmp/out.vcf
#
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_clinical_eval.sh --bundle DIR --out OUT.vcf [--workdir DIR] [--memory-mb N]

Bundle must contain:
  reference.fa
  tumor_R1.fastq
  tumor_R2.fastq
  normal_R1.fastq
  normal_R2.fastq

Optional (for truth comparison once implemented):
  truth.vcf(.gz)
  regions.bed
EOF
}

BUNDLE=""
OUT=""
WORKDIR=""
MEMORY_MB="1024"
TRUTH=""
REGIONS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bundle) BUNDLE="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --memory-mb) MEMORY_MB="$2"; shift 2;;
    --truth) TRUTH="$2"; shift 2;;
    --regions) REGIONS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

if [[ -z "${BUNDLE}" || -z "${OUT}" ]]; then
  usage
  exit 2
fi

if [[ -z "${WORKDIR}" ]]; then
  WORKDIR="$(dirname "${OUT}")/rosalind_eval_work"
fi

REF="${BUNDLE}/reference.fa"
T_R1="${BUNDLE}/tumor_R1.fastq"
T_R2="${BUNDLE}/tumor_R2.fastq"
N_R1="${BUNDLE}/normal_R1.fastq"
N_R2="${BUNDLE}/normal_R2.fastq"

for f in "${REF}" "${T_R1}" "${T_R2}" "${N_R1}" "${N_R2}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing required bundle file: ${f}" >&2
    exit 2
  fi
done

mkdir -p "${WORKDIR}"

echo "[1/4] Build (release)"
cargo build --release

VCF1="${OUT%.vcf}.run1.vcf"
VCF2="${OUT%.vcf}.run2.vcf"
W1="${WORKDIR}/run1"
W2="${WORKDIR}/run2"
mkdir -p "${W1}" "${W2}"

echo "[2/4] Run 1"
time cargo run --release -- somatic \
  --reference "${REF}" \
  --tumor-r1 "${T_R1}" --tumor-r2 "${T_R2}" \
  --normal-r1 "${N_R1}" --normal-r2 "${N_R2}" \
  --output "${VCF1}" \
  --workdir "${W1}" \
  --memory-mb "${MEMORY_MB}"

echo "[3/4] Run 2"
time cargo run --release -- somatic \
  --reference "${REF}" \
  --tumor-r1 "${T_R1}" --tumor-r2 "${T_R2}" \
  --normal-r1 "${N_R1}" --normal-r2 "${N_R2}" \
  --output "${VCF2}" \
  --workdir "${W2}" \
  --memory-mb "${MEMORY_MB}"

echo "[4/4] Determinism check (hash compare)"
shasum -a 256 "${VCF1}" "${VCF2}"
if ! cmp -s "${VCF1}" "${VCF2}"; then
  echo "ERROR: VCF outputs differ across runs" >&2
  exit 1
fi

if [[ -f "${W1}/somatic.manifest.txt" && -f "${W2}/somatic.manifest.txt" ]]; then
  shasum -a 256 "${W1}/somatic.manifest.txt" "${W2}/somatic.manifest.txt"
  if ! cmp -s "${W1}/somatic.manifest.txt" "${W2}/somatic.manifest.txt"; then
    echo "ERROR: manifest outputs differ across runs" >&2
    exit 1
  fi
fi

cp "${VCF1}" "${OUT}"
echo "OK: deterministic run. Final VCF: ${OUT}"

if [[ -n "${TRUTH}" ]]; then
  echo "[optional] Truth comparison"
  eval_args=(--reference "${REF}" --calls "${OUT}" --truth "${TRUTH}")
  if [[ -n "${REGIONS}" ]]; then
    eval_args+=(--regions "${REGIONS}")
  fi
  cargo run --release -- eval-somatic "${eval_args[@]}"
fi


