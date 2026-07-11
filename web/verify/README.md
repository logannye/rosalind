# Receipt Studio (entirely in your browser)

Receipt Studio is a static, framework-free client for Rosalind receipts. Drop run
receipts, reproduction certificates, and artifacts together to inspect trust evidence
without uploading data or contacting a third party.

It uses the same `rosalind-receipt` Rust crate as the CLI, compiled to WebAssembly.
The browser can:

- check the claim and measurement self-hashes;
- stream artifacts through BLAKE3 in fixed-size chunks and match by content, not name;
- distinguish receipt integrity, artifact completeness, resource-contract result,
  reproduction evidence, and signature status;
- localize causal parameter/output differences between two run receipts; and
- render a compact provenance chain for multiple receipts and certificates.

## Build

```sh
rustup target add wasm32-unknown-unknown --toolchain 1.83.0
cargo install wasm-pack --version 0.14.0
RUSTUP_TOOLCHAIN=1.83.0 ./scripts/build-wasm-verifier.sh
```

The committed package is generated and checked byte-for-byte on the Linux/amd64
release platform. Rust's optimized Wasm code generation may differ on other host
architectures even with the same target and tool versions; local builds remain
valid but are not the canonical release artifact.

## Run locally

ES modules and the WASM fetch require HTTP rather than `file://`:

```sh
cd web/verify
python3 -m http.server 8000
# open http://localhost:8000/
```

## Precise integrity guarantee

The deterministic claim protects content hashes, output-affecting parameters, replay
recipe, and producer/analyzer/build identity. The independent measurement hash protects
recorded resource telemetry. In schema 3 and newer, recorded paths are intentionally
portable metadata excluded from the claim hash, so a path-only edit is displayed as a
non-claim change rather than tampering.

“Receipt intact” does not mean artifacts were supplied, the computation was reproduced,
or an author was authenticated. Those are separate trust levels documented in
[receipt trust levels](../../docs/receipt-trust.md). Receipts are tamper-evident, not
signed; signature status is explicitly unavailable until signing ships.

## Native equivalent

```sh
rosalind verify --manifest path/to/sample.manifest.json --json
```
