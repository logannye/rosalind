#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [[ "$#" -lt 2 ]]; then
  echo "usage: $0 bundle.tar.gz (--candidate-source PATH | --registry-sdk)" >&2
  exit 2
fi
ARCHIVE="$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve(strict=True))' "$1")"
shift
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$ARCHIVE" -C "$WORK"
BUNDLES=("$WORK"/rosalind-*)
if [[ "${#BUNDLES[@]}" -ne 1 || ! -x "${BUNDLES[0]}/rosalind" ]]; then
  echo "expected one native Rosalind bundle" >&2
  exit 2
fi
python3 "$ROOT/scripts/onboarding.py" smoke --bundle "${BUNDLES[0]}" \
  --binary "${BUNDLES[0]}/rosalind" "$@"
