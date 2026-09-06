# Portable evidence validation — 2026-09-06

These checks validate the development tree, not a published candidate or independent
adoption. Fixtures are synthetic correctness/resource probes, not representative
biological performance benchmarks. Source tests reconstruct the inputs.

- [Full workspace](workspace.log), [Clippy](clippy.log), and [Rust 1.83](msrv.log) pass.
- [Python 3.11](python311.log) and [Python 3.9](python39.log): 23 tests pass each,
  including 5 portable dataset tests, three-budget byte equality, partial-overlap reuse,
  offline relocated replay, early termination, corruption, and Parquet cleanup.
- [Linux dataset Rust](linux-dataset-rust.log): 30 tests pass, with public library
  budget checks, source-role and symlink guards, malformed ownership rejection,
  source-mask memory costs and exact uint64 export.
- [Linux dataset Python](linux-dataset-python.log): 5 tests pass.
- [R 4.2.2](linux-r-adapter.log): a relocated 55-position query works without the
  original BAM/reference, paths containing spaces work, failure propagates, replay
  matches, and character parsing preserves uint64 maximum and 2^53+1.
- [DuckDB 1.4.4](sql-results.json) executes the documented SQL against the native
  Parquet export: 17,020 positions, 16,940 zero-depth, 960 callable observations,
  80 positions at 10x, 320 C observations and C-quality sum 12,640. Export peak is
  21,676,032 bytes under 128 MiB; this single small run is not a resource curve.
- [Fresh wheel smoke](wheel-smoke.log) installs outside the checkout, exercises
  portable Python queries/export and byte replay, and builds/runs packaged SDK
  examples with explicit candidate source patches. Registry-only SDK publication
  remains a separate gate.

Frozen final Linux source archive SHA256:
`dbd1cda670318c18d3bc1c82b2a57fc1a810c4ec461606885859b89c2332b3c3`.
Linux binary SHA256:
`953852e620a48ee2e92f06060d68aef57c5d9543473e1d8b7687ef09df8b4226`.
That snapshot precedes two focused reporting fixes: preserve invalid query REF as
an input error, and refresh the query receipt's finalization memory observation.
Those changes have focused Rust/Python validation in the development tree.

The Rust schema 1 and pre-extension projection golden bytes remain frozen. New
Parquet files preserve physical uint64 scalar and fixed-list types, with a direct
uint64 maximum roundtrip regression. Parquet directory byte replay is explicitly
unsupported; Arrow/TSV materialization supplies native byte replay.
