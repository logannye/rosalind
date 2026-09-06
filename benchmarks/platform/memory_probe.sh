#!/usr/bin/env bash
set -euo pipefail
implementation="${1:?implementation required}"
budget_mb="${2:?matched container/declaration budget required}"
case "$implementation" in
  rosalind)
    exec rosalind features --reference-pack /data/reference.rref --alignments /data/sorted.bam \
      --memory-budget-mb "$budget_mb" --max-depth 1000 --enforce --require-os-limit --output /outputs/features.tsv
    ;;
  pysam)
    exec python3 /opt/rosalind-platform/pysam_features.py --reference /data/reference.fa \
      --bam /data/sorted.bam --max-depth 1000 --output /outputs/features.tsv
    ;;
  bcftools)
    exec python3 /opt/rosalind-platform/bcftools_features.py --reference /data/reference.fa \
      --bam /data/sorted.bam --max-depth 1000 --output /outputs/features.tsv
    ;;
  *) exit 2 ;;
esac
