# Managed evidence analyzer validation — 2026-09-06

These are checks of local development candidates, not published packages or
independent external users. Native and wheel checks use explicit candidate source
patches for unpublished Rust dependencies. The examples run outside the workspace.

## Reproduction caught by fresh installation

The first local wheel passed its Python/dataset checks and built both standalone
scaffolds, but its generated evidence analyzer intermittently failed receipt
finalization. [The retained report](initial-wheel-failure.json) records the failed
128MiB run and resource partial. The sealing loop incorrectly required the observed
RSS high-water mark to reach a fixed point within eight serialization attempts.
The parent output/receipt remained absent. The correction and its validation are recorded below.

The pre-correction source archive `/tmp/rosalind-priority6-source.tar.gz` was
SHA256 `a277f9dd84d444340fe17485db5e0c40170f9cf632aa9f948f5989f9146bfc43`.
This archive predates the final receipt sealing correction and is not the final SDK.

## Shared resource and cancellation checks

- `receipt.log`: exact non-MiB boundaries, conflicting/malformed budgets, legacy
  MiB fallback, and unchanged historical schema 1–5 fixtures pass.
- `cancellation.log`, `signal-generations.log`, `signal-lifetimes.log`: scoped
  cancellation, signal restoration, old-token isolation, delayed signal races,
  and compile-time scope lifetime protection pass.
- `python39.log` and `python311.log`: all 23 source-package tests pass on both
  supported Python versions. Unindexed VCF/BCF diagnostics come from deliberately
  sequential fixtures, not failed assertions.

## Canonical browser verifier

`wasm-validation.json` records Linux/amd64 Rust1.83.0 and wasm-pack0.14.0,
73 native receipt tests and four native verifier tests. `wasm-first.log` and
`wasm-second.log` retain two clean generations with identical package bytes.
Actual Node/WebAssembly checks cover exact byte boundaries, contradictory claims,
and historical MiB fallback. Generated WASM SHA256:
`beb699f879f36f1fbc46893e9b41bd1b5d9163b22ae2e7ee3c08098dfd192de1`.

Subsequent PR #103 CI exposed a checkout-path dependency in Cargo's metadata
namespace: same-path rebuilds matched, but a different checkout changed the bytes.
Commit `9db09f7` gives the standalone verifier a stable package/version namespace
and adds an alternate-checkout regression. Both paths now produce SHA256
`56bfc6734c15e539cb52c8b5ac34f9d07cbdcdf575dc34eefbbfb65982cfeb73`;
the actual WASM boundary/legacy checks and GitHub verifier gate pass. The earlier
same-path result above is preserved with its narrower scope.

## Corrected SDK checks

The runner now patches a fixed-width decimal peak and the measurement digest in
the already allocated JSON buffer. It observes after encoding, file sync and
source revalidation; the receipt names this sampling phase. A final budget breach
changes the transaction to explicitly identified partials. Cancellation publishes
neither a success nor a partial result.

- `artifact.log`: eight lifecycle integration tests pass, including three-budget
  and tile equality, native/projected saved-source equality, a 25,000-site
  pre-admission case, failure cleanup, input mutation and resource partials.
- `receipt-sealing.log`: a publication-time resource failure seals partials; the
  reseal preserves claim/measurement integrity, canonical roundtrip and allocation
  capacity/pointer for within/over/unset and zero/uint64-maximum measurements.
- `conformance.json`: the independently compiled example passes all 19 executable
  checks. Native originals are unavailable during relocated replay; saved-dataset
  analysis succeeds with original alignments/reference deleted. Custom suffixes
  replay physically and exact non-MiB budgets verify.
- `fresh-evidence-msrv.log` / `fresh-column-msrv.log`: new external generated
  projects resolve and test on Rust1.83 without borrowing the repository lockfile.
- `workspace.log`: the full workspace suite passed before the final narrow sealing
  correction; `artifact.log`, `receipt-sealing.log` and `conformance.json` cover
  that correction. Final workspace Clippy (`clippy.log`), Rust1.83 (`msrv.log`),
  28-guide closure (`docs.log`) and nine onboarding helper tests pass.

The frozen Linux validation archive is
`/tmp/rosalind-priority6-final-source.tar.gz`, SHA256
`ad5b01229f4b383b8efabf4cc9d304965ec5b3961e7e00cfcf8131deb3b017b3`.
It includes the corrected production runner; it precedes only the additional
allocation-preservation unit test and final documentation updates.

## Corrected fresh wheel

`wheel.log` and `wheel-validation.json` record a fresh Python3.11 installation,
matching CLI0.5.0, evidence/panel/dataset/replay/Parquet and README examples, both
actual CLI-generated external crates, legacy conformance12/12, evidence
conformance19/19, and all four packaged standalone analyzer tests. Cargo resolution
started with an empty cache and explicit candidate source patches; this is not a
registry-only SDK installation. `wheel-conformance.json` retains the corrected
all-19 result.

Wheel SHA256: `eb7ede881ff5a9e517d6b19cd8473cbeb049b1ef20f7431a39b31f7f96a8568c`.
Bundled native SHA256: `a284df54a5d5da4f0053b4520d5402f7005c34b40df3e347ee8d73d91271baba`.

## Linux native SDK

`linux-validation.json` and the three `linux-*.log` files retain the frozen
candidate's eight artifact and three cancellation integration checks, publication
resource regression, four standalone example tests, and all 19 executable
conformance checks. This is a Linux/amd64 container on the same Apple Silicon host
via Rosetta; it is platform validation, not an independent user's machine.
Native CLI SHA256: `b3ed85fec8c7de920d2983e8ef72e8c5c16860057a88e36f9f37ca47fbc12755`.
