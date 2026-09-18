#!/usr/bin/env bash
# Evidence reuse demonstration. Select an exact binary/version; no silent fallback.
# See docs/evidence-reuse-demo.md for setup and recording instructions.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec "${PYTHON:-python3}" "$SCRIPT_DIR/reproduce_demo.py" "$@"
