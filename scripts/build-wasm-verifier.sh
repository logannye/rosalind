#!/usr/bin/env sh
# Build the in-browser "caught-you" receipt verifier (wasm) into web/verify/pkg/.
#
# Requires:
#   rustup target add wasm32-unknown-unknown
#   cargo install wasm-pack --version 0.14.0
# CI sets RUSTUP_TOOLCHAIN=1.83.0 so committed output is reproducible at the MSRV.
#
# NOTE: RUSTFLAGS is cleared on purpose. A global `target-cpu` (e.g. the common
# `~/.cargo/config.toml` with rustflags = ["-C","target-cpu=native"]) resolves to a
# host CPU that is invalid for wasm32 and makes the wasm-bindgen step fail with
# "failed to find intrinsics to enable `clone_ref`". Clearing it avoids that.
set -eu
cd "$(dirname "$0")/.."
RUSTFLAGS="" wasm-pack build crates/receipt-wasm \
  --target web \
  --out-dir ../../web/verify/pkg \
  --out-name rosalind_verify
echo
echo "Built -> web/verify/pkg/. Serve it (ES modules need http, not file://):"
echo "  (cd web/verify && python3 -m http.server 8000)   # then open http://localhost:8000/"
