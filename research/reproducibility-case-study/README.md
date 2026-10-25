# A bounded real-read reproducibility case study

This command checks whether exact HG002 read measurements survive changed queries,
relocation, and a clean saved-only replay. It compares fresh extraction (A),
existing native saved queries (B), and checked zstd Parquet (D). D shares B's
native extraction and pays the same setup plus conversion. This does not compare
Rosalind with an independent raw-to-Parquet pipeline.

**Execution status:** implementation prepared; exact-source hosted results are
pending. No new real-data measurements or successful clean replay are claimed
by this document. The findings and claim ledger must identify the tested commit,
run and artifacts before an application result can be quoted.

The [frozen protocol](PROTOCOL.md) specifies one public HG002 BAM window,
chr20:[10,000,000,11,000,000), 1,000 reference-derived probes, legacy mask63,
MAPQ/BQ20, sample HG002, one worker, 256 MiB cooperative budget, and 16,384-base
tiles. These are evidence probes, not truth variants. Q1 uses the first 250
probes; Q2 changes ALT at 500 even-index probes; Q3 retains all 1,000 positions
and depth≥5/10/20 flags. All zero-depth observations remain. General-CIGAR pysam
independently checks every legacy field before downstream answers are compared.

## Run the exact source on hosted CI

The [two-job workflow](../../.github/workflows/reproducibility-case-study.yml)
builds/tests an exact requested commit, runs the unchanged NA18507 preflight,
then executes this case study. A second fresh job receives only the complete
compact bundle and the exact binary; it has no raw-data download or cache step.
The runner does not build, install packages, or find an executable on PATH.
Public stable v0.1.0 does not contain the required evidence engine.

From an authenticated clone, choose the complete commit SHA being evaluated:

```sh
SHA="$(git rev-parse HEAD)"
gh workflow run reproducibility-case-study.yml \
  --ref codex/reproducibility-case-study -f source_sha="$SHA"
```

The initial implementation push automatically runs the same workflow; do not
also dispatch it and create a duplicate attempt. Preserve failure artifacts;
allow only one explained correction. Repository CI remains a separate gate.
The workflow uses pinned Rust 1.99.1, Python 3.11.17, wheel-only dependencies from
[the existing lock](../cohort-reanalysis/workloads/requirements.lock), one native
target per job, and disabled incremental/dev/test debug data. It checks 20 GiB
free scratch and 8 GiB available RAM before work. Acquisition is about 1.98 GB;
preparation expands the full reference (about 3.14 GB) and creates unused CRAM,
chromosome and panel files. Only the selected BAM case is measured.

In that prepared hosted environment the complete case study is one command:

```sh
PYTHONPATH=python python research/reproducibility-case-study/run.py \
  --binary "$ROSALIND_BINARY" \
  --build-identity "$ROSALIND_BUILD_IDENTITY" \
  --expected-source "$ROSALIND_SOURCE_SHA" \
  --protocol research/reproducibility-case-study/protocol-v1.json \
  --download-cache "$RUNNER_TEMP/rosalind-public-downloads" \
  --output "$RUNNER_TEMP/rosalind-case-study"
```

Outputs must be new paths. The explicit build record binds binary SHA256, native
compiled commit/tree/compiler/lock identities, the protocol, dependencies and
relevant source files. A version string is only an assertion. Authentication
checks consistency with the trusted hosted record; unsigned hashes do not prove
authorship or protect against coordinated replacement of all artifacts.

## Verify only the saved bundle

Download `replay-bundle-RUN-ATTEMPT` and `exact-binary-RUN-ATTEMPT` from the same
identified run into separate directories. Hosted retention is 30 days; rerun
the exact commit and pinned public sources if artifacts expire. The binary
artifact contains `rosalind` and `build-identity.json`; no native target is needed.
Checkout that exact source SHA and install the pinned wheels without installing
the repository. Obtain the expected inventory hash from the producer job output
or independently retained `artifact.json`, outside the bundle being checked.

```sh
PYTHONPATH=python python research/reproducibility-case-study/verify_bundle.py \
  --bundle "$ROSALIND_REPLAY_BUNDLE" \
  --binary "$ROSALIND_BINARY" \
  --build-identity "$ROSALIND_BUILD_IDENTITY" \
  --expected-inventory-sha256 "$ROSALIND_INVENTORY_SHA256" \
  --output "$RUNNER_TEMP/rosalind-independent-replay"
```

Use a clean environment with the producer's original input/cache paths absent.
The verifier reconstructs requests and expected answers from the retained
independent oracle, checks every inventory item and full legacy field, reruns
B/D, and requires matching-code native `REPRODUCED` for the original Q1 TSV and
unaltered receipt. The reduced JSON answer receives a separate outer check.
It repeats missing/REF, synthetic zero/missing and consumed-payload corruption
controls. A producer's `passed` field cannot satisfy these gates.

`report.json` retains all 27 query observations, parent-clock prefixes, individual
verification times, actual exits, output hashes, source rechecks, process-tree RSS
and setup/conversion costs. `trials.csv` provides the compact trial table.
The replay bundle includes every native partition, compressed Parquet and
inventory, requests, expected results, unmodified native receipts, protocol,
identity, oracle and producer report. Raw data stays outside it. A 128 MiB
ceiling is enforced, with a 64 MiB target. Raw receipt paths are retained as
provenance; other than public sources, no personal data or credentials belong in
these artifacts.

## Interpret the result

There is no speed-win requirement. Shared setup is measured once; query medians
and full ranges have three observations each. Prefixes come from the parent
clock at verified publication, include startup once, and are never sums of
inner medians. Native integrated hashing remains inside total query time;
separately observable verification is reported alongside it. D's explicit
inventory/payload verification is timed separately but already included in total.
RSS is a 10 ms sampled simultaneous process-tree sum, not a hard memory cap.
No cold-cache, p99, extrapolated crossover, clinical outcome, release acceptance
or outside-user adoption follows.

Both earlier stops remain visible: the [bounded cohort analyzer](../../docs/findings/bounded-cohort-reanalysis-2026-10-22/README.md)
was correct but slower than existing saved queries and checked Parquet; the
[frozen HG001 threshold transfer](../threshold-transfer/CONTINUATION_RESULTS.md)
removed two true candidates and failed practical-benefit gates. Neither is
reopened. See [claim ledger](CLAIMS.md) and [interview notes](INTERVIEW.md).
