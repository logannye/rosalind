# Evidence engine delivery roadmap

Rosalind turns existing indexed alignments into reusable research evidence. A
researcher chooses the loci and scientific filters; the engine chooses a bounded
execution strategy and records inputs, semantics, and results for verification.

The invariant across every phase is: **identical inputs and scientific settings
produce identical successful results across memory budgets, tile sizes, and
worker counts**. Insufficient resources cause refusal or failure, never hidden
downsampling. Speed, measured memory, biological accuracy, and adoption require
different evidence and are reported separately.

These are staged release targets. The development package remains **0.4.0** until
a release is authorized and published. See [implementation status](implementation-status.md)
for implemented, verified, released, and adopted states.

## M0 — 0.4.0 stabilization

- Preserve existing filter defaults; make pileup capacity limits exact-or-fail and
  record `pileup.semantics=exact-or-fail-v1` in new recipes.
- Preserve historical receipt schemas 1–5 verification. Historical scientific
  replay requires the original producer selected with an explicit binary.
- Stream canonical merge, separating scientific recipe compatibility from run
  outcomes. Include Arrow transients and the full lifecycle in memory accounting.
- Repair scheduled validation and genuinely streaming comparison baselines.
  Require fresh GIAB evidence for changed caller/default/shared-input behavior,
  without making unrelated evidence-engine work depend on caller accuracy.
- Validate one candidate SHA across Linux x86_64 and macOS arm64/x86_64 binaries
  and wheels, Linux OCI, minimum Python, CLI inventory, and unique RC versions.
  Invoke publication jobs explicitly from the authorized release workflow.

Exit: regression and historical-receipt tests, release-policy tests, fresh package
smoke, and candidate platform builds pass. Publication also needs registry/trusted
publisher configuration, protected environment approval, and a reviewed caller
baseline when watched source changed. Those external prerequisites do not prevent
completing or reviewing the code.

## M1 — 0.5 exact SNV evidence

- Indexed BAM with BAI/CSI, CRAM with CRAI and explicit local FASTA; analysis
  reference `.fa` plus `.fai`, `.rref`, or compatible `.idx`.
- Supplied plain-text SNV VCF or every locus of a BED union, including zero depth.
  `shortread-dna-readcount-v1`: MAPQ20, BQ20, read counts, excluded unmapped,
  secondary, supplementary, duplicate, and QC-failed reads; unavailable qualities
  never pass a threshold.
- Prefilter/aligned/callable depths, allele/strand counts, quality sums and
  histograms, sequencing-cycle sums, and per-locus filter counts. Checked integer
  aggregation over bounded indexed tiles and canonical 1,024-row Arrow batches.
- Rust evidence API, Python `iter_evidence`/`materialize_evidence`, and a small
  real-origin alignment tutorial joining supplied candidates to a research screen.

Exit: independent observation oracle, adversarial CIGAR/filter fixtures,
budget/tile invariance, empty/zero-depth selections, relocated replay, and fresh
wheel execution. A real researcher exercises the RC; seven-day soak and external
release gates begin here. Caller precision is not an extractor validation metric.

## M2 — 0.6 panel QC and builder API

- Original BED targets keep full-length denominators, uncovered positions,
  overlaps, and configurable callability. Deletions/reference skips contribute
  no aligned-base coverage. No threshold has clinical meaning.
- Optional reference for coverage-only work; fused position evidence and target
  summaries from one extraction pass.
- `EvidenceRequest`, `EvidenceBatch`, `EvidenceAnalyzer`, and explicit field,
  reference/context, and retained-memory requirements. Schema v1 can validate
  requirements while emitting full rows; physical projection is separate work.
- Keep `ColumnAnalyzer` compatibility with unknown memory models observable but
  ineligible for a declared bound. Maintain a standalone analyzer and equivalent
  Nextflow/Snakemake workflows.

Exit: an external crate builds and consumes the SDK; target summaries match a
semantically aligned oracle; fused output equals separate runs; a real workflow
completes from prepared inputs to verified artifacts.

## M3 — 0.7 reusable datasets and measured execution

- Opt-in cache/resume with verified atomic canonical partitions keyed by input
  content, reference, normalized selection, profile, fields, and schema. Hash each
  distinct input once per run; rehash on later runs before reuse.
- Bounded parallel extraction with canonical output; custom reducers remain
  serial. Incomplete/corrupt cache entries cannot become successful evidence.
- Dataset comparison identifies changed loci/metrics and changed inputs, profile,
  and analyzer. Measure hashing, extraction, encoding, and verification costs.

Exit: interrupted/resumed, cold/warm-cache, and single/multiple-worker results
match byte-for-byte; changed content/profile invalidates reuse; three-repeat
memory/time/I/O curves expose small-tile costs; an external consumer reuses a
dataset and an independent user replays it.

## Validation and adoption

| Dimension | Minimum retained evidence |
|---|---|
| Semantics | Independent CIGAR oracle; all filters/unavailable qualities; deletions/skips/clips; overlapping mates; duplicate/overlapping/empty selections |
| Invariance | Three admitted budgets, several tile sizes, one/multiple workers; equal successful bytes; separate refusal/runtime-limit tests |
| Memory | Whole-process peak including hash/decode/analyzer/encoder/verification; actual cgroup limit where claimed; depth/length/target-count curves |
| Time/I/O | Three repeats, warm/cold cache distinction, exact versions/argv/hashes, setup separate, end-to-end cost retained |
| Reuse | Fresh/resumed, corrupt/missing partitions, changed profile, relocated inputs, independent replay |
| Installation | Fresh binaries/wheels on promised platforms, minimum Python, outside-checkout import, exact candidate identity |

Adoption targets: three non-authors, one external analyzer, one integrated
workflow, and two teams returning after 30 days. Record install/integration time,
semantic rework, and measured resource effects. Stars, authored examples, and CI
are not substitutes for these gates.

## Deferred

Production calling, indels, UMI/fragment models, methylation, long reads,
single-cell analysis, pangenomes, plugin registries, and a new search-index
algorithm are outside this sequence. Revisit when concrete users demonstrate a
need that the existing evidence contract cannot meet.
