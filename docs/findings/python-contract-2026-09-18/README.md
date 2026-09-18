# Python contract fingerprint validation — 2026-09-18

The exact committed candidate `222b8e1e2b0182bb3b19ba3545e7d18ff6cb30ed` passed `cargo xtask contract snapshot`: **63 components**, status `ready`, no blockers, exit 0. Its detached worktree was clean before and after the run. The previous exact-commit baseline `bc6e9e6e051d233da98573dc40f894f56ccd58b5` also passed its older 62-component snapshot, which did not fingerprint Python.

| Snapshot | Aggregate BLAKE3 | Components |
| --- | --- | ---: |
| [Baseline bc6e9e6](baseline-bc6e9e6/contract-snapshot.json) | `cc493862f498187d12ece7c8a73b8bc024a5659694deaddca3069b427feaa5a1` | 62 |
| [Candidate 222b8e1](candidate-222b8e1/contract-snapshot.json) | `8724eed303c014c2b52d74337eab5f935599a8fd2a841743fd2ae41ef2cafd91` | 63 |

The [comparison](comparison.json) shows exactly one added component, `python.contract` (`546cd430a14f97ca0b25ee0b38e58c5dd31ca90c19db4ddddd8b384decf0eeea`), and one changed component, `package.rosalind-bio.files`. The retained package inventories grew from 332 to 335 paths: `scripts/python-api-contract.py`, `scripts/test_python_api_contract.py`, and `scripts/test_installer.py`. No paths were removed. The other 61 component fingerprints, all 10 primary artifact fingerprints, receipt-field inventories (8 measurements, 74 parameters, 6 top-level fields), and package versions were identical.

## Coverage and execution

The final snapshot covers 49 CLI help entries, 2 Rust public APIs, 3 schema files, 5 legacy fixtures, 3 Cargo package file inventories, and the new Python component. That component parses all 5 Python package modules into a normalized AST, including implementation and defaults, and includes packaging/dependency metadata. It ignores comments, docstrings, source locations, empty optional AST fields, and the explicitly excluded descriptive project metadata. The Python package is never imported or executed by this fingerprint step. Cargo package components fingerprint file inventories, not every file's contents.

Both full snapshots ran on macOS arm64 in separate detached worktrees, using an explicitly borrowed private build cache and cached Cargo dependencies with network access disabled for Cargo. The workspace's shared build target was not used. The candidate ran from `17:00:31` to `17:01:26` UTC on 2026-09-18. The native fixture checks and generated Rust analyzer conformance are part of the snapshot command. Exact commands, paths, timing, environment, exit status, compiler configuration hash, and tool versions are retained in each snapshot directory's `invocation.json` and `runtime-tools.json`.

Pinned tools: Rust/Cargo toolchain `1.95.0` (rustc `59807616e`), public-API toolchain `nightly-2026-04-15` (rustc `1.97.0-nightly`, `a5c825cd8`), and `cargo-public-api 0.50.1`. The final snapshot's AST helper used Python `3.14.6`.

```sh
cargo +1.95.0 xtask contract snapshot \
  --ref 222b8e1e2b0182bb3b19ba3545e7d18ff6cb30ed \
  --output /absolute/path/contract-snapshot.json --json
```

## Supplemental implementation checks

The following previously completed checks are retained unchanged; they were not repeated as part of the exact-commit snapshot run:

- [Rust tests](rust-tests.log): 32 unit tests and 5 integration tests passed; no failures.
- [Script tests](script-tests.log): 41 tests passed. These include mocked onboarding/installer tests and should not be read as an installed-package runtime smoke.
- [Clippy](clippy.log) and [MSRV](msrv.log): successful xtask checks.
- [Cross-version AST comparison](cross-version.json): Python 3.9.6, 3.11.15, and 3.14.6 produced identical canonical bytes for all 5 modules: 179,800 bytes, SHA256 `b667369a352a735a163589c22c09bfbd87fd682b503ea06db9022926792e400d`. This validates source canonicalization across those interpreters; it does not import the package or exercise Python runtime behavior.

This evidence validates the local source contract and its added Python coverage. It does not constitute a published release-candidate attestation, all-platform wheel validation, or proof that every possible behavioral change is captured. Earlier snapshots remain readable but lack the new Python component, so they cannot establish equivalence to this candidate's complete fingerprint. No release was published by this validation.

[SHA256SUMS](SHA256SUMS) records the bytes of the retained evidence and this README.
