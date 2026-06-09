# The "caught-you" receipt verifier (in your browser)

A static page that checks a Rosalind reproducibility receipt's tamper-evident BLAKE3
self-hash **entirely client-side** — drag in a `*.manifest.json`, then edit one byte and
watch it flip to **TAMPERED**. It runs the *same* Rust verification that ships in
`rosalind verify`, compiled to WebAssembly from the [`rosalind-receipt`](../../crates/receipt)
crate (std + blake3, no htslib).

## Build

```sh
# one-time: rustup target add wasm32-unknown-unknown && cargo install wasm-pack
./scripts/build-wasm-verifier.sh        # -> web/verify/pkg/  (71 KB wasm)
```

## Run locally

ES modules + the `.wasm` fetch need http (not `file://`):

```sh
cd web/verify && python3 -m http.server 8000
# open http://localhost:8000/
```

The page loads with a real sample receipt; edit any character to break the hash. Drag a
`*.manifest.json` from a `rosalind variants` / `features` run onto the box to check your own.

## What it proves (and doesn't)

It re-derives the receipt's own integrity — `manifest_blake3` over the canonical claim,
plus the independent `measurement_blake3` — exactly the self-hash check `rosalind verify`
runs, minus re-hashing the input/output *files* (which aren't in the browser). **VERIFIED**
means intact + untampered; a single edited byte breaks the hash. It does **not** re-run the
analysis, and it is tamper-*evident*, not tamper-*proof* (cryptographic signing is a
separate, planned step).

## Native equivalent

```sh
cargo run -p rosalind-receipt --example verify_file -- path/to/sample.manifest.json
```

## Deploy

`pkg/` is a build artifact (git-ignored). To publish this as a shareable link, build it in
CI and serve `web/verify/` from GitHub Pages — a follow-up.
