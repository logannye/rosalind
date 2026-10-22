# Bounded saved-evidence reanalysis: correctness passed, optimization stopped

The research analyzer now produces exact candidate-SNV summaries over explicit
portable specimen leaves, preserves missing observations, refuses incompatible
measurements, and supports verified saved-data replay. The diagnostic pilot did
**not** show a useful performance advantage over the existing saved-query reducer
or a competent Parquet cache. Keep the simple correctness implementation; stop
optional coalescing and further performance expansion under the specification's
G2 rule. This completes the working vertical slice and its lifecycle/value check,
not the proposed representative G4/G5 study or public cohort roadmap.

## Exact evidence identity

[Hosted run 39816575369](https://github.com/logannye/rosalind/actions/runs/39816575369)
tested PR merge `a7f3b40daff9aa319a0e089c5af1ef7f37fd5ec4`, incorporating branch
head `9627c936c5e84121c97b5167cc694df3f78d3bce` into baseline
`8874508844e106484fa290b04e3fc6349686beaf`. This is PR validation, not an actual-main
or release acceptance. Subsequent findings/documentation commits do not change
which source produced these measurements.

The native SHA256 is
`b2638ccc3b5e734f2c15dacbdf1b4e223cd2c71d778c877aabc2acfeb4e320ce`.
The clean build used x86_64 Linux, Rust 1.99.1 and Python 3.11.17, with PyArrow
20.0.0, pysam 0.23.3, psutil 6.1.1 and blake3 1.0.5. Full platform and dependency
identity is in [identity.json](identity.json). The pilot records four logical CPUs;
a specific CPU model was not retained, limiting cross-host interpretation. Generated validation outputs explain
the later untracked directory in the benchmark's source-status record; recorded
Python source hashes were checked against the local frozen implementation.
All native compilation and executable native tests ran on hosted CI.

## Correctness and lifecycle results

The focused run passed `cargo test --locked --workspace`, including native saved
data, coverage, replay and CRAM regression gates, then built the release CLI. It
also passed 29 existing Python tests, **58 research tests** (47 pure checks and
11 native integration cases), and 38 independent-oracle regression tests. Logs
are retained in [validation/](validation/). These are test cases, not biological
specimens or clinical validations.

The unchanged three-specimen fixture and independent validator passed. Its four
summary rows reproduce across windows 2 and 1 with the generated raw paths hidden;
all portable leaves verify. At position 40, the result retains three requested
specimens but two observed specimens: positive counts in A, missing B, measured
zero in C. [The retained summary](demo/summary.jsonl), [plan](demo/plan.json),
[demo report](demo/report.json), and [replay certificate](demo/replay.json) make
that distinction inspectable. Tests additionally cover windows 1/2/4 with native
budgets 128/256/512 MiB, parent admission 96/128 MiB, thresholds 5/10/20, duplicate
and reordered requests, new ALTs, empty requests, wholly absent candidates,
integers above 2^53, overflow, corrupt payloads, earlier-input mutation, changed
coverage masks, cancellation, failed consumers and create-new publication.

Current saved-only reanalysis of a historical **HG002 leaf with 1,000 loci across
62 partitions** matched retained independent CIGAR-oracle outputs exactly:
3,000 summary rows at thresholds 5/10/20, identical at windows 128 and 512.
Both summaries have SHA256
`f8c045cca2ca312b7a9963554ad6b4e1763f2616fcf93682e815b7b6e487666a`.
See [the report](validation/archived-hg002.json) and
[input provenance](validation/archived-inputs.json). The raw extraction and oracle
were originally run at `9b8e12f`; they were **not rerun** here. This is current
saved-data correctness evidence for one real specimen, not fresh extraction,
amortization, or representative multi-specimen performance.

## Locked synthetic diagnostic

[protocol.json](protocol.json) preserves the workload committed before the pilot:
three distinct synthetic named BAM sources, 1,024 loci spread over four ownership
partitions per member, explicit missing masks, and three sequentially revealed
requests. Prefixes K=1/2/3 include overlapping requests, new ALTs, depth screens
and revisited partitions. Every path sees the same observations and masks; no
future request knowledge is used. There are three rotating repetitions of each
path/prefix, **45 successful trials with exact agreement and zero trial failures**.
All trials, commands, output hashes, ranges and failures are in
[results.json](results.json); per-request work is in
[trial-metrics.json](trial-metrics.json).

The table leads with **one measured setup plus median whole-worker query time**,
in seconds. Setup includes extraction/materialization and validation of the
respective cache; it is measured once, so these totals have no repeated-setup
uncertainty estimate. Synthetic input generation is outside both setup timers.

| Path | K=1 total | K=2 total | K=3 total |
|---|---:|---:|---:|
| B: existing verified saved queries + simple exact reducer | 0.689 | 0.738 | 0.755 |
| C: analyzer, W=64 | 0.872 | 1.151 | 1.384 |
| C: analyzer, W=256 | 0.856 | 1.071 | 1.251 |
| D: bounded Parquet, performance only | 0.452 | 0.475 | 0.469 |
| D: bounded Parquet, content checked | 0.449 | 0.468 | 0.470 |

Native setup took 0.497 seconds; columnar setup took 0.247 seconds. The Parquet
baseline independently counts the same generated BAM observations, preserves
stored zero rows and masks, uses uint64 columns, zstd, partition/predicate pruning,
requested-ALT projection, 1,024-row batches and controlled readahead. Its
content-checked mode hashes consumed files before and after querying and checks
inventory identity. That is additional integrity work, not native receipt parity.
B has native per-invocation verification; C additionally performs comparison,
planning, consumed-input union checks, source/binary checks and outer publication.

| Path at K=3 | Query median [range], seconds | Whole-workflow sampled peak median [range], MiB |
|---|---:|---:|
| B | 0.258 [0.256, 0.259] | 103.83 [103.64, 104.29] |
| C, W=64 | 0.887 [0.865, 0.897] | 110.88 [110.64, 110.91] |
| C, W=256 | 0.754 [0.703, 0.757] | 109.73 [109.70, 110.19] |
| D, performance only | 0.223 [0.216, 0.230] | 115.05 [114.72, 115.58] |
| D, content checked | 0.223 [0.216, 0.245] | 115.58 [115.16, 115.95] |

Memory uses the maximum of the setup and query simultaneous process-tree RSS
samples, never their sum. Sampling is every 10 ms and may miss short peaks.
Budgets are cooperative, not OS-enforced caps. Cache state was not evicted or
controlled; three rotations are not a fully balanced order design. Small timing
differences between D modes are not evidence that hashing accelerates queries.
No p99 or cold-cache claim is made.

## Mechanism and stop decision

At K=3, W=64 launched 27 native query children and made 30 partition visits, with
15 repeated decodes within requests. W=256 launched nine query children, made
15 visits and repeated none within requests. Including metadata commands and
version probes, these paths launched 114 and 78 native processes respectively;
B launched 18. Increasing the window reduced query time but did not yield a
material sampled-memory tradeoff on this small workload. The recorded work also
includes roughly 94–96 MB of parent hashing before finalization and a further
45 MB final union check per three-request prefix, including the native executable
and provenance files. These counters identify extra work; they are not an isolated
time attribution to hashing or decoding.

C at W=256 cost about 1.66 times B and 2.66 times content-checked D in total at K=3.
It used about 5% less sampled memory than D and more than B. C did not beat either
baseline's total time at any measured prefix. **No C-over-baseline crossover was
observed.** No raw A comparison was performed, so raw-versus-saved crossover is
unestablished. Do not fit or extrapolate a future crossover from three prefixes.

The added contract and replay behavior is useful correctness work, but this pilot
does not establish a performance advantage worth extending for the application.
Retain the simple analyzer; do not pursue optional coalescing, select a post hoc
resource target, or run the reserved held-back workload to search for a win.
The held-back configuration remains unexecuted. The broader real-read A/B/C/D
study and final scale target are deliberately not claimed complete. Mechanisms,
falsifying checks and retained failures are detailed in [decisions.md](decisions.md).

## Storage, retention and reproduction

Logical retained input sizes were 355,948 bytes of synthetic raw inputs,
762,765 bytes of native portable cache, and 37,134 bytes of columnar cache/metadata.
Native setup also retained 469,600 bytes of aggregate Arrow exports, 18,488 bytes
of their receipts and 3,258 other bytes: total 1,254,111 bytes. Those exports are
not the minimum needed for later leaf queries. All 45 retained query outputs sum
to 44,884,634 bytes. These are logical file sizes, not measured peak temporary
disk usage; the raw and cache formats provide different provenance envelopes.

[artifact.json](artifact.json) identifies the full hosted bundle and its archive
digest. Its hosted retention is 14 days; the reports, all trials, metrics,
validation logs, identities and selected output records here are durable in git.
Full generated datasets remain outside git in the local downloaded artifact.
`demo/run-manifest-record.json` is an audit copy referring to the full run's
inventory, not a standalone complete run directory. Use the downloaded
`retained-demo/run` for wrapper verification; do not point native receipt tools at
the outer Python aggregation or assume staged child paths relocate automatically.

Reproduce using the exact-source hosted
[workflow](../../../.github/workflows/cohort-reanalysis.yml) and the commands in
[walkthrough.md](walkthrough.md). The saved HG002 input is selected from the already
tracked, hash-pinned historical archive; no bulk data acquisition is required.
The full artifact's four completed run directories were also SHA-verified locally
without executing native code. Earlier diagnostic runs and the telemetry failure
remain in [development/index.json](development/index.json).

Historical evidence remains separate: the
[3,159-position reuse demo](../evidence-reuse-demo-2026-09-18/README.md),
[108-run representative study](../representative-evidence-2026-09-06/README.md),
and [corrected CRAM rerun](../cram-decoder-validation-2026-09-06/verified-9b8e12f-2026-09-18/README.md)
retain their own revisions, failures and memory-plan findings. They do not certify
this extension. [Interview notes](interview-notes.md) connect the representation,
missingness and evaluation lessons to ML engineering without claiming clinical
or brain-interface performance.
