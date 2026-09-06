# Short-read evidence semantics

`analyze evidence` and `analyze panel-qc` use evidence schema **1** and profile
**shortread-dna-readcount-v1**. Existing `features`, `analyze coverage`, and
`variants` keep their legacy defaults. Do not mix their similarly named columns
without explicitly matching filtering rules.

## Inputs, loci, and coordinates

Input is local, coordinate-sorted indexed BAM (BAI/CSI) or CRAM (CRAI). CRAM
requires an explicit local FASTA and adjacent FAI; no network reference lookup is
needed. Analysis references may be uncompressed FASTA plus FAI, `.rref`, or
compatible `.idx`. Contig names and lengths must match the alignment dictionary.

Supply exactly one plain-text SNV VCF (`--sites`) or BED (`--regions`). VCF REF
and each ALT must be a single A/C/G/T base. REF is checked against the reference;
conflicting duplicates, indels, symbolic alleles, unknown contigs, and invalid
coordinates fail. Duplicate loci combine ALT alleles and emit one row. VCF
genotypes, FILTER labels, and INFO annotations do not change read filtering.

BED uses zero-based half-open coordinates. Its normalized union defines evidence
loci, including uncovered positions. Overlap does not duplicate evidence rows.
Output `pos` is one-based; Rust `EvidenceRow.position` is zero-based. Order follows
the reference dictionary and then numerical position. Empty selections produce a
valid empty artifact. Selection is scientific; execution tiles are not selection.

## Observation and filtering rules

Sample scope is resolved before locus filtering. By default, one unambiguous
`@RG SM` sample is selected automatically; a header without sample names has
explicitly unknown sample identity. Multiple declared samples, or a mixture of
named and unnamed read groups, require `--sample NAME` or `--pool-samples`.
Named selection counts only the selected sample's read groups and refuses reads
whose group is missing, undeclared, or has no sample assignment. Multiple read
groups for the same sample are combined. Explicit pooling includes all records,
including unassigned records, and is recorded as pooled rather than named.

Receipts record the resolved scope as versioned `evidence.sample_scope` JSON;
scientific and cache identities include it. `--sample NAME` and automatic
selection of the same named sample have the same scientific scope. Reads from
other named samples contribute no locus counters. Skipped indexed record visits
are execution diagnostics, not unique biological read counts.

The unit is **one aligned read base**, not a fragment or UMI consensus. Overlapping
mates count separately, and orphaned/improper pairs are not automatically removed.
Only CIGAR M, =, and X contribute matched observations. D and N consume reference
coordinates without contributing depth; I and S consume read bases without a
reference observation. H/P contribute neither. There is no BAQ adjustment,
overlap correction, duplicate discovery, local realignment, or allele calling.

Unmapped records have no counted locus. At each matched locus:

1. `prefilter_depth` counts matched observations before profile filtering.
2. Exclude secondary, supplementary, QC-failed, duplicate-flagged reads in that
   order, then unavailable MAPQ255 and MAPQ below the threshold (default 20).
3. `aligned_depth` counts the remaining matched observations.
4. Exclude unavailable BQ255, BQ below the threshold (default 20), then non-ACGT
   observations. `callable_depth` and allele/strand statistics count the remainder.

Filter counters are **exclusive first-failure** counts at a locus. A read with
multiple excluded flags increments the first applicable counter. Missing quality
never passes, even at threshold 0. Supported BQ values are 0–93, reported MAPQ0–254.
Profile overrides are explicit recipe values and alter scientific identity.

The SAM sequence symbol `=` means equality to the reference and resolves to the
supplied reference base ([SAM specification, section 1.4](https://samtools.github.io/hts-specs/SAMv1.pdf)).
An ambiguous reference remains noncallable. Coverage without reference sequence
refuses a quality-qualified `=` observation with guidance to provide a reference;
it does not silently omit that observation.

Thus `prefilter_depth >= aligned_depth >= callable_depth`, and A+C+G+T equals
callable depth. `callable` here means passing the declared technical filters, not
that a genotype is reliable or clinically interpretable.

## Counts, quality, and cycles

All reducers use checked integers. Allele counts are A/C/G/T; strand counts use
alignment strand. Quality sums and histograms include only callable observations:
94 BQ bins and 255 MAPQ bins, with no bin for missing255. Histograms permit exact
downstream quantiles without retained reads. TSV histograms encode ordered
`quality:count` pairs; `.` means all bins zero.

`read_position_sum` uses the zero-based offset in stored SEQ: reverse-strand reads
map the aligned query position back through `read_length - 1 - query_position`.
Stored sequence length includes soft clipping and excludes hard-clipped bases.
It does not reconstruct cycles removed before alignment. Derive a mean by dividing a sum by
callable depth; zero depth has no defined quality/cycle mean.

## Panel QC

The BED union is extracted once, but each original nonempty target receives its
own summary in input order. Overlapping targets each receive the overlapping
loci; equal target identifiers do not merge targets. BED column4 names a target;
otherwise a line-based identifier is supplied.

Mean callable depth divides by the **complete target length**. Uncovered bases
remain in the denominator and minimum depth. Breadth counts positions reaching
1x/10x/20x/30x callable depth. `--min-callable-depth` (default 10) supplies a separate
callable-position threshold. Without a reference, BAM dictionary coordinates are
used and reference bases are N; no reference-dependent statistic is inferred.
`--position-output` writes Arrow evidence from the same pass as the target summary.
The Rust, CLI, and Python panel interfaces share the 10x default; Python delegates
the default to its matching native executable and passes only explicit overrides.

## Execution and memory

Evidence extraction never downsamples to satisfy a budget. The engine scans
indexed intervals into bounded aggregate tiles, rather than retaining every active
read. Budget changes may reduce microtile width and repeat indexed reads. A single
record/read beyond the declared envelope or integer overflow fails the run.
Successful scientific output is invariant to budget, microtile, and worker count.

Canonical ownership tiles span 16,384 reference bases. Arrow encoding uses fixed
**1,024-row** batches independently of computation width. This differs from the
legacy feature encoder's 65,536-row batches. Schema v1 always emits full rows;
field requirements are validated capabilities, not physical projection savings.

Plans include engine, consumer/encoder, and process overhead. The htslib record
envelope is checked after decode, so it is a cooperative check rather than an
allocation cap inside htslib. A declared budget and sampled RSS cannot prove a
universal hard OS bound. An existing Linux cgroup-v2 limit is an additional,
separately recorded assurance. Python memory retained by a consumer is outside the
native child-process budget.

Legacy pileup instead retains active observations within declared read/depth
capacities. New runs record `pileup.semantics=exact-or-fail-v1` and fail at capacity
instead of silently choosing a subset. Old receipts still verify at their original
schema capability; replay of old behavior requires the matching producer.

## Identity, replay, and reuse

Scientific identity includes content, normalized selection, profile, fields, and
schema. Execution settings and outcomes are recorded separately. Input files are
hashed once per distinct path identity in a run; later runs rehash before cache
reuse. Verified atomic partitions may be reused with `--cache-dir --resume`.
Workers extract first-party partitions; custom reducers consume serially in
canonical order.

All alignment, reference, selection, and index files must remain immutable for
the entire run. Initial content hashes establish identity; before/after metadata
and inode guards detect ordinary in-place edits or file replacement during
execution. These metadata checks are not cryptographic re-verification and do not
defend against an intentional modification that restores the observed metadata.
Workers aborted by an input-change guard must not publish changed-source cache
partitions. Keep mutable upstream jobs and evidence extraction in separate stages.

Materialized outputs and their receipts support byte verification and replay.
A stdout stream has no persisted output artifact unless the consumer explicitly
stores one. Receipts are tamper-evident records, not signatures or proof of
biological accuracy. See [receipts and trust](receipts-and-trust.md).
