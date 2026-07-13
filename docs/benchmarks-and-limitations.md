# Benchmarks and limitations

Benchmark results are evidence artifacts, not product adjectives. Each published
result must identify exact inputs and hashes, command argv, code identity, platform,
resource controls, output hashes, and a reproduction path.

The executable platform comparison runs equivalent feature extraction through Rosalind,
pysam, and bcftools where possible using one BAM/reference. It reports setup steps,
peak RSS, wall time, output size, repeat byte identity, insufficient-cgroup behavior,
and second-machine reproduction effort. Resource predictability is reported
separately from speed and scientific accuracy.

The pinned implementation is in `benchmarks/platform/`. The report is generated
from retained GNU-time files, outputs, environment locks, receipts, and hashes.
bcftools-only semantic gaps are explicit fields in the report. The harness also
asserts repeated TSV/Arrow identity and sharded Arrow merge identity. No real-data
number is published until the HG002 run completes and its artifacts are reviewed.

The reference-construction comparison will measure legacy `.idx`, `.rref`, BWA
indexing, and a relevant external-memory constructor where outputs are comparable.
Its purpose is to demonstrate that analyzers no longer depend on FM-index
construction, not to market a new indexing algorithm.

Current scientific boundaries: the built-in caller is diploid, short-read, and
SNV-focused; the committed accuracy evidence is simulated; the pinned HG002 GIAB
baseline remains pending. Current platform boundaries: coordinate-sorted BAM,
BAI for sparse access, single-threaded execution, first-party-only merge codecs,
no CRAM/CSI, and no published Linux ARM wheel.
