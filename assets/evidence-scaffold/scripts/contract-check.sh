#!/usr/bin/env bash
set -euo pipefail
ROSALIND_BIN="${ROSALIND_BIN:-rosalind}"
cargo build --release --locked
"$ROSALIND_BIN" conformance analyzer --api evidence \
  --binary target/release/__PACKAGE_NAME__ --json
