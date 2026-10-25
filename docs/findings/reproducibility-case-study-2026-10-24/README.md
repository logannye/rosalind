# Exact HG002 evidence, saved queries, and clean replay

The bounded case study passed at source `401b06b12b875e13e4c101305f41d4bb0ed252fc` in
[hosted run 40151881458](https://github.com/logannye/rosalind/actions/runs/40151881458) on October 24 Pacific / October 25 UTC, 2026.
An independent general-CIGAR pysam traversal matched all 34 legacy columns at
1,000 HG002 probes. All 27 timed query answers matched expectations. A separate
clean hosted job verified the complete saved evidence, reran all three native
and Parquet requests, and returned matching-code native `REPRODUCED` for the
original Q1 TSV and unmodified receipt without the original alignments/reference.

This result establishes exact measurements and reproducible downstream analysis
for one frozen workload. These reference-derived probes are not truth variants.
It does not establish clinical accuracy, production adoption, BCI benefit,
independent raw-to-Parquet superiority, or a general speed or memory guarantee.
The application wording and interview explanation are in the
[claim ledger](../../../research/reproducibility-case-study/CLAIMS.md) and
[interview notes](../../../research/reproducibility-case-study/INTERVIEW.md).

## What was fixed before the run

The [protocol](../../../research/reproducibility-case-study/PROTOCOL.md) uses public
HG002 chr20:[10,000,000,11,000,000), explicit sample HG002, MAPQ/BQ20, legacy
all/mask63, one native worker, 256 MiB cooperative budget, and 16,384-base tiles.
Preparation preserved the full 195-contig matching reference dictionary and all
261,761 source records overlapping the window; observed maximum read length was
151. Native sparse extraction reported 248,188 record visits and 62 microtiles;
record visits and prepared source-record counts have different denominators.

Q1 is the first 250 probes with original ALT. Q2 is the 500 even-index probes with
the second sorted non-REF base as ALT. Q3 is all 1,000 probes with per-locus
depth≥5/10/20 flags. A/B project depths, alleles and strands. The full source
representation retains uint64 counts, sums and histograms. Zero observations stay
in the answers. Query definitions and expected outputs are bound in the bundle.

A performs fresh native extraction for each request. B queries the existing native
saved dataset. D queries checked zstd Parquet converted from the same verified
native Arrow evidence. B and D both pay native materialization/verification setup;
D additionally pays conversion and readback. D is a downstream representation
comparison, not an independently populated raw-to-Parquet implementation.

## Costs with explicit boundaries

All durations below are seconds. The acquisition cache was initially empty.
Common stages were measured once. RSS is the largest simultaneous worker-plus-
descendants sum sampled at 10 ms, so short peaks can be missed.

| Stage | Wall seconds | Sampled peak RSS, bytes |
|---|---:|---:|
| Unchanged small NA18507 preflight | 0.798649 | 102,367,232 |
| Shared public acquisition/preparation | 33.937795 | 297,918,464 |
| Independent raw pysam oracle | 1.357455 | 48,762,880 |
| Common B/D materialization and verification | 1.325328 | 81,440,768 |
| D conversion, validation and readback | 0.405434 | 117,485,568 |

Common B/D setup includes extraction, receipt verification, portable copy, full
dataset verification, full TSV export/receipt verification, and comparisons of
native TSV/Arrow with the oracle. Oracle counting is the separate evaluation stage
above. Conversion includes zstd Parquet creation, schema/value validation and
readback. Build, dependency setup, controls and final artifact auditing are not
charged as query time. Hosted tool/dependency setup through the binary identity
record took 212 seconds by retained whole-second timestamps; the one hosted native
target occupied 1,549,678,853 bytes before job cleanup.

Three balanced rounds ran A/B/D, B/D/A, D/A/B, with Q1/Q2/Q3 in each worker.
Each cell below is median [minimum, maximum] across the three parent-clock
prefixes, plus the single measured setup. Each prefix starts at worker launch and
ends when Qk is verified and published, including startup once. It is not a sum
of per-query medians.

| Path, including its setup | Through Q1 | Through Q2 | Through Q3 |
|---|---:|---:|---:|
| A fresh | 0.938530 [0.933609, 0.956589] | 1.908173 [1.902705, 1.909947] | 2.880463 [2.855964, 2.890056] |
| B native saved | 1.370024 [1.368487, 1.403491] | 1.399018 [1.394912, 1.434242] | 1.435907 [1.431209, 1.471755] |
| D checked Parquet | 1.862467 [1.861837, 1.873144] | 1.874178 [1.873874, 1.884770] | 1.895494 [1.893129, 1.906313] |

Both saved paths first have a lower median than A at the observed K=2 prefix.
There is no claim of a general crossover or statistical significance. Native B's
K=3 median is 1.435907 s versus A's 2.880463 s and D's 1.895494 s for this protocol.
If the common 33.937795 s acquisition/preparation is added, these arithmetic totals
become A 36.818258 s, B 35.373702 s, D 35.833289 s; oracle/controls/audit costs remain
separate. The result is not a twofold acquisition-inclusive pipeline speedup.

The inner query durations below already include their observable verification.
They exclude initial worker startup that the parent prefixes above include.
The verification column is a subset of query time, not an additional charge.

| Path / query | Median [min, max], seconds | Explicit verification median, seconds |
|---|---:|---:|
| A / Q1 | 0.907398 [0.893659, 0.927213] | 0.430271 |
| A / Q2 | 0.968790 [0.953064, 0.969352] | 0.429209 |
| A / Q3 | 0.971703 [0.952729, 0.979560] | 0.421255 |
| B / Q1 | 0.014505 [0.014040, 0.027285] | 0.002460 |
| B / Q2 | 0.028735 [0.026179, 0.030499] | 0.004219 |
| B / Q3 | 0.036321 [0.035794, 0.036962] | 0.006259 |
| D / Q1 | 0.101938 [0.101270, 0.112768] | 0.006761 |
| D / Q2 | 0.011966 [0.011427, 0.012004] | 0.004334 |
| D / Q3 | 0.018333 [0.017632, 0.021072] | 0.005665 |

D was faster than B on Q2 and Q3. D's first query includes initial PyArrow work;
B was faster for the complete setup-inclusive sequence. A repeatedly hashes the
full reference, including integrated native hashing and a separate receipt
verification. Native integrated hashing remains inside query time and cannot be
added again. A/B's explicit column includes the external receipt verifier; D's includes
internal inventory/payload/schema checks. Both also include canonical result
loading/readback/equality checks, so these subcolumns are not identical
instrumentation scopes. This experiment does not isolate faster read counting.

| Path | Post-preparation setup/query sampled peak RSS, bytes |
|---|---:|
| A | 30,724,096 |
| B | 81,440,768 |
| D | 117,485,568 |

These path maxima are maxima across setup/conversion/query phases, never sums of
nonconcurrent peaks. The original report key calls these `workflow_max...`; its
scope excludes shared preparation. The largest sampled phase was shared
acquisition/preparation at **297,918,464 bytes**. The cooperative 256 MiB native
budget excludes Python consumers and is not a whole-workflow hard cap.

| Retained item | Logical bytes |
|---|---:|
| Downloaded and authenticated public source files | 1,980,479,678 |
| All prepared files, including unused CRAM/panel/chromosome outputs | 3,162,748,946 |
| Minimum native saved leaf (payloads and metadata) | 3,732,895 |
| Minimum checked Parquet (payload and inventory) | 197,471 |
| Materialization exports and receipts | 3,388,782 |
| All query outputs and logs | 5,190,873 |
| Complete replay bundle, 161 files including its inventory | 5,840,123 |
| Exact binary, supplied separately | 13,917,792 |

The native leaf comprises 3,347,992 payload bytes and 384,903 metadata bytes;
Parquet comprises 34,649 payload bytes and 162,822 inventory bytes. These are
logical file sizes, not disk allocation or peak scratch requirements. The
complete replay bundle is below the 64 MiB target and 128 MiB refusal ceiling.

## Correctness and refusal evidence

The direct-CIGAR oracle is `benchmarks/evidence/streaming_pysam.py`. Native full TSV
was byte-equal over 1,001 lines including its header; full Arrow had zero
mismatches over 1,000 loci × 34 columns. Complete native and Parquet values were
also compared to the retained oracle in the clean job. All A/B/D query answers
matched independently derived expectations; all six saved B/D answers were
recomputed after relocation. Parquet agreement checks format preservation rather
than a third independent read-counting implementation. Native and pysam still
share input files, htslib decoding and the written semantics, so a shared
specification error remains possible.

The clean replay checked 39 absent original input/cache paths, fully verified
124 partition files in the 62-partition native dataset, resolved native replay
inputs by content, and retained explicit `execution.alignment_record_visits=0`
and `original_sources_rehashed=false` in every saved-query receipt. It reproduced Q1's native TSV bytes with the original receipt and
matching compiled source identity. The custom reduced answer was separately
rederived and checked. Absence of counters was never interpreted as zero work.
The isolation boundary is a separate clean hosted job with only the bundle and
exact binary downloaded; this is not evidence of outside-user adoption.

Corrupting one byte in a copied native Arrow payload consumed by Q1, while leaving
manifests unchanged, caused dataset verification and querying to refuse with exit
5 and no successful final output/receipt. The intact original remained unchanged.
Absent-locus requests refused with exit 3 and REF-mismatch requests with exit 2,
also without successful publication. Checked Parquet refused missing and REF-mismatched requests. A
separate labeled synthetic fixture returned a stored zero-depth position while
refusing an adjacent unmeasured position; the real probes had no zero-callable
rows. Both hosted jobs repeated these
controls; the verifier recomputed them instead of trusting producer flags.

## Exact-source validation and retained attempts

The producer used Ubuntu 24.04 x86_64 (image 20261023.334.3), Rust 1.99.1,
Python 3.11.17, pysam 0.23.3, pyarrow 20.0.0, psutil 6.1.1, blake3 1.0.5 and
numpy 2.2.6. Its exposed CPU was AMD EPYC 9V45, four logical CPUs. Replay ran on
a separate Ubuntu job exposing AMD EPYC 7763, four logical CPUs, with the same
binary and Python dependency versions. The replay environment's available Rust
compiler was not used to build anything. Build identity is checked against the
compiled binary, not a version label alone.

`cargo fmt --check`, locked workspace native tests and release build passed.
Python suites executed 3 demo, 29 native integration, 38 oracle, and 38 new
case-study tests: 108 total, no skips. The unchanged public NA18507 preflight
passed. The [validation records](validation/) retain exact test outputs and
resource checks. No native compilation or executable native tests ran locally.
The later findings/documentation commit is not the measured scientific source.

The [first attempt](attempt-1/README.md), run 40151648804 at `385e65b`, failed
before acquisition because an existing generated-analyzer acceptance test looked
under a private target directory while Cargo used the configured shared target.
The one allowed correction made that test resolve and use the same target path.
No measurement-engine Rust code changed. This successful second attempt exhausted
the authorized scientific attempts; no query/configuration search or third run
was performed. Both the initial failure and correction remain visible.

## Artifact audit and retrieval

[Combined checked disposition](checks.json), [producer report](report.json), [all 27 trials](trials.csv),
[clean replay report](replay-report.json), [build identity](build-identity.json),
[artifact anchor](artifact.json), [independent artifact audit](independent-artifact-audit.json),
and [independent cost audit](cost-audit.json) are retained here. The independent
local audits rechecked hashes, retained answers, ordering, prefix arithmetic and
storage; they did not rerun native binaries or Arrow decoding locally. The
hosted clean verifier performed those executable checks.

The producer's original `complete_case_study_status` says `awaiting_clean_replay`,
and the earlier copy inside the bundle still has bundle creation marked running.
Those immutable records predate the second job. The separate completed replay
report and provider job records close the result; no original report or receipt
bytes were rewritten to manufacture completion.

Download the four artifacts from the identified hosted run. Hosted retention is
30 days. Original ZIPs were also preserved locally under
`research/reproducibility-case-study/artifacts/40151881458/` (ignored by Git).
[Archive verification](archive-verification.json) recomputes each ZIP's SHA256 and
matches the provider digest; [provider records](provider-artifacts.json) identify
the archive IDs. Full native payloads are in the replay bundle, not the producer's
compact diagnostics. The exact binary is a separate archive. Reproduce with the
[documented saved-only command](../../../research/reproducibility-case-study/README.md#verify-only-the-saved-bundle)
from the exact tested checkout; the runner never builds or installs dependencies.
Artifact retention is finite, so retain the ZIPs before expiry rather than treating
the run URL as permanent storage. Any new scientific rerun needs fresh scope.

Binary SHA256:
`9bf4ce1f897b3a299bf8c17bac500ff9df01227f6ec9ff900e732840e17a5705`.
Externally retained bundle inventory SHA256:
`2aaf4c42a4a916d5351d38647496dee4b89cdcc0e0c80ac0309e6e09babe5398`.
These checks bind artifacts to the trusted hosted record. Unsigned hashes do not
prove authorship or protect against coordinated replacement of all provenance.

## Earlier decisions remain stopped

The [bounded cohort analyzer](https://github.com/logannye/rosalind/blob/70af512ecb803baa61b79abcf68ccb6d42cb990f/docs/findings/bounded-cohort-reanalysis-2026-10-22/README.md) was correct but slower: at K=3,
analyzer W256 took 1.251 s, existing native saved queries 0.755 s, and independently
populated checked Parquet 0.470 s (merge `a7f3b40daff9aa319a0e089c5af1ef7f37fd5ec4`,
run 39816575369). Optional optimization stopped. That stronger independent
Parquet comparison is separate from this study's shared-native-input comparator.

In [frozen HG001 threshold transfer](https://github.com/logannye/rosalind/blob/70af512ecb803baa61b79abcf68ccb6d42cb990f/research/threshold-transfer/CONTINUATION_RESULTS.md), false retained candidates fell
69→53 and all retained 5,697→5,672, while true retained fell 4,983→4,981 (source
`10096419d00c18a03f64098c8c5f35701560bb59`, run 39997188782). Practical-benefit gates
failed and the simple baseline remained. Unknown labels remain distinct. This
case study does not reopen acquisition, refitting, threshold choice, a held-back
cohort run, or optimization.
