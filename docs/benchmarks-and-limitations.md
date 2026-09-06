# Benchmarks and limitations

Benchmark results are evidence artifacts, not product adjectives. Each published
result must identify exact inputs and hashes, command argv, code identity, platform,
resource controls, output hashes, and a reproduction path.

The legacy platform comparison runs related feature extraction through Rosalind,
pysam, and bcftools using one BAM/reference. Only explicitly matched fields and
filter rules support equivalence claims. It reports setup steps,
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

The new `benchmarks/evidence` harness uses direct independent pysam CIGAR traversal
under `shortread-dna-readcount-v1` and compares every emitted evidence field. Its
default three budgets and three repeats retain input/package identity, raw timing
and RSS, output hashes, and verification time; optional cache runs measure reuse.
The tiny real-origin NA18507 tutorial is a semantic/repeatability check, not a
whole-genome performance comparison. Hashing, setup, and combined analysis/encoding
costs are retained separately where the implementation measures them. It imposes
no OS cap and does not attribute hard-budget assurance to the baseline.

The caller remains diploid, short-read, and SNV-focused; committed accuracy evidence
is simulated and the HG002 GIAB baseline remains pending. Legacy feature selection
retains BAM/BAI constraints. New evidence supports BAM/BAI/CSI and CRAM/CRAI with
explicit local reference, and bounded first-party workers; custom reducers remain
serial. Linux ARM wheels, physical field projection, UMI/fragment counting, and
production calling remain outside the current contract.
