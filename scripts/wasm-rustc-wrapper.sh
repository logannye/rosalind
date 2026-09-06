#!/usr/bin/env bash
# Confined to the standalone receipt verifier build. Cargo unifies features for
# a package/version on one target; keep different package versions disjoint while
# removing the checkout path from the wasm crate's symbol namespace.
set -euo pipefail
compiler="$1"
shift
is_wasm=false
for argument in "$@"; do
  if [[ "$argument" == "wasm32-unknown-unknown" ]]; then
    is_wasm=true
  fi
done
if [[ "$is_wasm" != true || -z "${CARGO_PKG_NAME:-}" || -z "${CARGO_PKG_VERSION:-}" ]]; then
  exec "$compiler" "$@"
fi
arguments=()
while (( $# )); do
  case "$1" in
    -C)
      if [[ "${2:-}" == metadata=* ]]; then
        shift 2
        continue
      fi
      ;;
    -Cmetadata=*)
      shift
      continue
      ;;
  esac
  arguments+=("$1")
  shift
done
arguments+=("-Cmetadata=rosalind-receipt-verifier-v1:${CARGO_PKG_NAME}:${CARGO_PKG_VERSION}")
exec "$compiler" "${arguments[@]}"
