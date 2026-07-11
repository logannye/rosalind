#!/usr/bin/env sh
# Build the in-browser "caught-you" receipt verifier (wasm) into web/verify/pkg/.
#
# Requires:
#   rustup target add wasm32-unknown-unknown
#   cargo install wasm-pack --version 0.14.0
# CI sets RUSTUP_TOOLCHAIN=1.83.0. The committed bytes are canonical on the
# Linux/amd64 release platform; other host architectures may emit valid but
# byte-different optimized Wasm.
#
# NOTE: RUSTFLAGS is replaced on purpose. A global `target-cpu` (e.g. the common
# `~/.cargo/config.toml` with rustflags = ["-C","target-cpu=native"]) resolves to a
# host CPU that is invalid for wasm32. Remapping the checkout and Cargo registry
# also prevents macOS/Linux absolute paths from leaking into the committed Wasm.
set -eu
cd "$(dirname "$0")/.."
repo_root=$(pwd)
cargo_home=${CARGO_HOME:-"$HOME/.cargo"}
# Populate the registry source directory before deriving canonical remaps. This
# matters in a fresh CI home where the loop below would otherwise see no crates.
RUSTFLAGS="" cargo fetch --manifest-path crates/receipt-wasm/Cargo.toml --locked
separator=$(printf '\037')
encoded_flags="--remap-path-prefix=$repo_root=/workspace${separator}--remap-path-prefix=$cargo_home=/cargo"
for registry_source in "$cargo_home"/registry/src/*; do
  [ -d "$registry_source" ] || continue
  encoded_flags="${encoded_flags}${separator}--remap-path-prefix=$registry_source=/cargo/registry/src/index"
done
RUSTFLAGS="" CARGO_ENCODED_RUSTFLAGS="$encoded_flags" \
  wasm-pack build crates/receipt-wasm \
  --target web \
  --out-dir ../../web/verify/pkg \
  --out-name rosalind_verify
echo
echo "Built -> web/verify/pkg/. Serve it (ES modules need http, not file://):"
echo "  (cd web/verify && python3 -m http.server 8000)   # then open http://localhost:8000/"
