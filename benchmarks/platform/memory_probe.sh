#!/usr/bin/env bash
set -euo pipefail
implementation="${1:?implementation required}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cp /data/reference.fa "$work/reference.fa"
cp /data/input.bam "$work/input.bam"
cp "$work/input.bam" "$work/sorted.bam"
if ! samtools index "$work/sorted.bam" 2>/dev/null; then
  rm -f "$work/sorted.bam.bai"
  samtools view -b -F 4 "$work/input.bam" |
    samtools sort -m 8M -o "$work/sorted.bam"
  samtools index "$work/sorted.bam"
fi
samtools faidx "$work/reference.fa"
case "$implementation" in
  rosalind)
    rosalind reference build --fasta "$work/reference.fa" --output "$work/reference.rref" >/dev/null
    # Must refuse before creating the destination under an intentionally tiny declaration.
    set +e
    rosalind features --reference-pack "$work/reference.rref" --alignments "$work/sorted.bam" \
      --memory-budget-mb 1 --max-depth 1000 --enforce --output "$work/refused.tsv"
    code=$?
    set -e
    test "$code" -eq 3 && test ! -e "$work/refused.tsv"
    exit "$code"
    ;;
  pysam) python3 /opt/rosalind-platform/pysam_features.py --reference "$work/reference.fa" --bam "$work/sorted.bam" --output /dev/stdout ;;
  bcftools) python3 /opt/rosalind-platform/bcftools_features.py --reference "$work/reference.fa" --bam "$work/sorted.bam" --output /dev/stdout ;;
  *) exit 2 ;;
esac
