#!/usr/bin/env bash
# The "reproduce duel": Rosalind re-derives a result byte-for-byte from its receipt —
# offline, on a different machine — something a non-deterministic caller cannot do.
#
#   cargo build --release && ./scripts/reproduce_demo.sh
#
# To record it as an asciinema cast (for the README / a "Show HN"):
#   asciinema rec -c ./scripts/reproduce_demo.sh reproduce.cast
set -euo pipefail

ROSALIND="${ROSALIND:-./target/release/rosalind}"
[ -x "$ROSALIND" ] || ROSALIND="./target/debug/rosalind"
[ -x "$ROSALIND" ] || { echo "build first: cargo build --release"; exit 1; }
D="examples/data/illumina_toy"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT

echo "# 1. Call variants. Every run writes a content-addressed BLAKE3 receipt."
"$ROSALIND" index    --reference "$D/reference.fa" --output "$work/ref.idx"           >/dev/null
"$ROSALIND" sort     --input "$D/alignments.bam"   --output "$work/sorted.bam"        >/dev/null
"$ROSALIND" variants --index "$work/ref.idx" --alignments "$work/sorted.bam" \
                     -o "$work/calls.vcf"                                             >/dev/null
echo "  -> $work/calls.vcf  (+ calls.vcf.manifest.json)"
echo

echo "# 2. Hand a stranger ONLY the receipt + the inputs. They re-derive it, offline,"
echo "#    with no GATK / Docker / Nextflow / re-aligning:"
echo
"$ROSALIND" reproduce --manifest "$work/calls.vcf.manifest.json" --inputs "$work"
echo

echo "# A non-deterministic caller (GATK, DeepVariant) would report DIVERGED on a"
echo "# CORRECT re-run — its output is not byte-deterministic. Rosalind's is, so"
echo "# \"reproduce this exact result\" is a command a stranger can run."
