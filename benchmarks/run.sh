#!/usr/bin/env bash
# Rosalind claims harness — one command, zero downloads, runs on the bundled toy data.
#
#   bash benchmarks/run.sh
#
# Each claim is a property you can re-derive; any failure exits non-zero (so this is also
# a CI gate). It is about the VERIFIABLE MEMORY CONTRACT + BYTE-REPRODUCIBILITY — not
# speed, not real-world accuracy. See benchmarks/README.md.
set -euo pipefail
cd "$(dirname "$0")/.."

BIN="${ROSALIND_BIN:-}"
if [ -z "$BIN" ]; then
  echo "==> building the release binary (cargo build --release)…"
  cargo build --release -q
  BIN="target/release/rosalind"
fi

ROSALIND_BIN="$BIN" python3 benchmarks/claims.py
