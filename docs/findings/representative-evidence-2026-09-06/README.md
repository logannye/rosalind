# Representative evidence validation — 2026-09-06

**Baseline result:** all 108 requested runs completed and all emitted fields agreed
with an independent streaming pysam traversal. Equivalent BAM and CRAM results
matched. The post-run audit verified unchanged inputs, binary and harness; every
Rosalind run stayed below its declared budget. It also found **eight CRAM memory
planning underestimates**. Scientific invariance passed; the baseline decoder
memory model needs correction. These findings do not establish publication,
independent adoption, a high-depth assay benchmark, or a universal memory bound.

## Workloads and identities

The pinned HG002 NovaSeq PCR-free chr20 BAM contains 19,430,746 mapped records
and is 1,093,938,924 bytes. Source URLs and SHA256 values come from the existing
[GIAB resource lock](../../../benchmarks/giab/resources.tsv); [preparation.json](preparation.json)
records the actual source verification and derivation. The full matching 195-contig
reference dictionary is retained, including the 3,144,230,986-byte FASTA.

| Workload | Selection | Reads available |
|---|---|---|
| Chromosome BAM, sparse | 999 reference-derived probes across chr20 | Original chr20 BAM |
| Window BAM and CRAM, sparse | 1,000 probes across chr20:10,000,000–11,000,000 (0-based) | 261,761 overlapping reads |
| Window BAM and CRAM, panel positions | Twenty separate 500-base BED targets; 10,000 positions | Same 261,761 reads |

No reads were simulated or replicated. The matched CRAM3.0 was encoded locally
from the same unfiltered source records, using pysam 0.23.3 / htslib 1.21.
Read lengths are at most 151 bases. The panel is an analysis request on WGS data,
not a capture assay or a high-depth cohort. Probe ALT alleles are arbitrary;
these are observation tests, not validated candidate variants or caller scores.
One all-ambiguous chromosome bin has no probe, giving 999 rather than 1,000 loci.

The macOS arm64 runs used the installed candidate wheel from the
[SDK validation](../evidence-sdk-2026-09-06/README.md), not an unverified debug binary:

- Wheel SHA256: `eb7ede881ff5a9e517d6b19cd8473cbeb049b1ef20f7431a39b31f7f96a8568c`.
- Installed CLI SHA256: `a284df54a5d5da4f0053b4520d5402f7005c34b40df3e347ee8d73d91271baba`.
- Manifest SHA256: `6dc5cc1b60213fda93a1da2bcafb5dfff85a7195fb5cf10eaa65fe9ff03a8c82`.

The source-candidate build identity and environment are retained in every native
measurement. This wheel predates the later CRAM planner correction. Do not label
these results as measurements of that correction or as registry installation.

## Design and results

[report.json](report.json) retains every argv, three repetitions per case, raw
resource observations, phase measurements, byte hashes, field comparisons, and
verification outcome. Jobs were shuffled with seed 17; cold/resume/persisted jobs
kept their required ordering. All runs use the same read-count profile, MAPQ/BQ20,
explicit HG002 sample, and the original full field set (mask 63).

Serial cases pair 96/128/256 MiB with tile limits of 128/4,096/16,384 bases.
This is a **budget and tile configuration sweep**, not an isolated experiment on
budget alone. Window workloads also use 2 workers at 512 MiB and 8 workers at
1,024 MiB, followed by cache resume and serial saved-dataset extraction. The
chromosome workload uses one worker. The streaming baseline runs three times per
workload and uses bounded 1,024-base windows.

![Baseline resource curves](resource-curve.png)

Points show medians and bars show the three-run range. Dotted RSS lines show the
baseline plan, making its CRAM underprediction visible. Full values, including
parallel, cache, chromosome and oracle runs, are in [medians.csv](medians.csv).
The [plot script](plot.py) regenerates the PNG and [SVG](resource-curve.svg).

| Serial case at 256 MiB / 16,384-base limit | Native wall, median | Native RSS, median | Separate verification, median | pysam wall, median |
|---|---:|---:|---:|---:|
| Chromosome BAM, 999 loci | 3.004 s | 10.42 MiB | 1.985 s | 0.796 s |
| Window BAM, 1,000 loci | 1.689 s | 9.94 MiB | 1.558 s | 0.609 s |
| Window CRAM, 1,000 loci | 2.029 s | 31.05 MiB | 1.556 s | 4.268 s |
| Window BAM, 10,000 positions | 1.606 s | 10.72 MiB | 1.551 s | 0.328 s |
| Window CRAM, 10,000 positions | 1.565 s | 22.72 MiB | 1.454 s | 0.464 s |

Native wall time includes setup, input hashing, extraction, encoding and output
finalization. Verification is measured separately and rechecks source content.
The pysam baseline extracts equivalent observations but does not perform Rosalind's
source hashing, receipt publication or verification. Its wall time is therefore
not an equal-trust end-to-end comparison. Recorded phase timings can overlap;
setup and hashing must not be added as if all phases were disjoint.

Smaller windows changed actual work. For window BAM sparse loci, the three serial
cases used 1,000 / 245 / 62 microtiles and 39,233 / 207,924 / 248,188 record visits.
Coalescing reduces seeks but may decode reads in unselected gaps. For the panel,
128-base tiles used 80 microtiles and 5,769 visits, versus 20 and 3,417 for larger
tiles. The receipt's tile width is an admitted limit; it is not a measurement of
the maximum actual fetch span. [audit.json](audit.json) corroborates execution
changes using microtile and record-visit counters.

## Repeated work and its costs

All 12 cache-resumed runs had zero alignment-record visits and preserved evidence.
Resume still rehashes the original inputs; on this small window workload that cost
dominates, so resume alone is not a general speedup. Saved-dataset extraction also
has zero alignment visits and verifies the stored dataset instead of reopening
original alignments. Across the four window workloads its median wall time was
0.041–0.060 seconds, with another 0.011–0.024 seconds for separate verification.
These are repeat-query costs after the initial materialization.

The complete three-repeat study retained 208,433,016 logical cache bytes across
1,020 files and 60,454,330 artifact bytes across 216 files. These totals include
independent repeated caches and receipts, not the storage required for a single
sample. The audit inventories logical and allocated storage separately. Initial
creation, cache storage and verification must be included in any amortized savings.

## Memory finding and actual Linux limits

The independent post-run audit reports eight predictions below measured CRAM RSS,
all within the requested budget. The largest observed shortfall was approximately
40 MiB RSS against a 21 MiB plan, at a 96 MiB budget. CRAM decodes slices and codec
workspaces whose sizes are not controlled by the genomic microtile. A per-record
allowance and generic decoder slack did not account for that coexistence.
The original observations remain retained while a conservative preflight model
and explicit supported decoder envelopes are implemented and tested separately.

The same source candidate was built for Linux/amd64 and tested in a separate
Colima VM on this Apple Silicon host. This is a second execution environment, not
an independent physical machine. Baseline binary SHA256:
`b3ed85fec8c7de920d2983e8ef72e8c5c16860057a88e36f9f37ca47fbc12755`.
[BAM](cgroup-bam.json) and [CRAM](cgroup-cram.json) each passed seven scenarios:

- Both nominal 128 MiB modes completed, verified, and produced equal bytes. One
  recorded cooperative assurance; the other checked and recorded the existing OS cap.
- An oversized declared record envelope caused typed preflight refusal (exit 3).
- A 1 MiB startup budget and a record-capacity case produced resource failures
  (exit 4), without a successful primary output. Partial receipts are classified
  by their content, not only their filenames.
- A separate allocation control and an 8 MiB native Rosalind run were killed
  (exit 137), with both Docker `OOMKilled` and increasing `memory.events.oom_kill`.

The native OOM run is separate from the control. Raw `memory.max`, `memory.swap.max`,
`memory.events`, `memory.peak`, exit codes, output hashes and verification results
are retained. Cgroup memory includes file-cache charges and reached its 128 MiB cap;
it must not be called process RSS. No host limit or protected environment was changed.

## Retained evidence and reproduction

[baseline-raw.tar.xz](baseline-raw.tar.xz) contains the complete 108-run directory:
TSVs, Arrow cache partitions, manifests, command logs, raw timing, verification
outputs, report, audit and the exact executed harness snapshot. It expands to
approximately 263 MiB. [BAM raw cgroup evidence](cgroup-bam-raw.tar.xz) and
[CRAM raw cgroup evidence](cgroup-cram-raw.tar.xz) include their executed harness
copies. [retained-identities.json](retained-identities.json) gives archive hashes.
Original genomic inputs are downloaded from the locked public sources, not
embedded here. The actual [preparation script](preparation-script.py) and final
[workload configuration](workloads.json) are preserved separately: the final
configuration tightened the small tile to 128 bases and enabled persisted queries.

To run a new study, follow the [representative harness guide](../../../benchmarks/evidence/README.md)
and [Linux cgroup guide](../../../benchmarks/evidence/cgroup-probe.md). New
preparations have their own path/tool-derived identities and provenance; do not
rewrite this baseline to make them match. The harness has since added direct
budget/final-identity gates and more robust failed-run retention. Its later source
hash is not the source hash used for these measurements.

No filesystem-cache eviction was performed. “Cold” means no Rosalind cache existed
for that case; it does not mean a cold OS cache. Raw filesystem counters are
operations, not decoded bytes or total I/O; their zero values on this cached run do
not establish zero I/O. Timings were collected on a development workstation with
other validation activity. The oracle establishes the tested valid-input
observations, not equivalent rejection of every malformed record. High-depth
assays, other CRAM codecs/layouts, whole-genome throughput, aligned samtools or
bam-readcount comparisons, and independent users remain additional validation.

## Adoption readiness

The [researcher and builder kit](../../adoption-validation.md) supplies task
protocols and blank measurement forms. Its researcher commands were executed
against the pinned tutorial: candidate review matched after saved-evidence reuse,
receipts verified, and replay reproduced. [Smoke output](adoption-researcher-smoke.log),
[schema checks](adoption-schema-validation.log), and [kit tests](adoption-all-tests.log)
are retained. These are authored checks. No non-author sessions or 30-day returns
have been collected, and no outreach or follow-up has been scheduled by this work.
