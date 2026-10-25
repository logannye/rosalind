# Interview notes

A measurement engine turns reads into a specified set of counts. Before choosing
a storage format, define the counting unit, reference, sample, filters, requested
positions and available fields. Here the unit is an aligned read base; overlapping
mates count separately. Changing an ALT query is possible because all A/C/G/T
counts were retained. Lowering MAPQ/BQ later is impossible because excluded
observations were discarded. Unstored positions, unretained fields, fragment/UMI
identity and broader biological interpretation cannot be reconstructed.

Native extraction matching native replay is useful for portability, but both can
share the same bug. The direct-CIGAR pysam traversal applies the scientific rules
independently and checks complete legacy fields. It still shares source files,
htslib decoding, and the written semantics; agreement cannot establish reference
truth or rule out a shared specification error. Arrow-to-Parquet readback checks
format preservation, not a third independent extraction implementation.

Observed zero says a stored locus had no qualifying read observations. Missing
means the locus was never measured in this selection. Unknown labels refer to
biological adjudication in a different evaluation; they are neither false nor
zero. A rejected query is an execution outcome, not a scientific observation.
Keeping these states distinct prevents silent denominator and accuracy errors.

The experiment reports acquisition/preparation, evidence materialization and
verification, columnar conversion/readback, query execution, explicit verification,
and retained bytes. Whole-worker and parent-clock prefix measurements include
startup. Native fresh extraction repeatedly hashes the full reference, so a
saved-query benefit would not isolate faster counting. The experiment does not
separately identify every cost of decoding, counting, startup and hashing. A
10 ms process-tree RSS maximum is an observation that can miss short peaks; the
native budget is cooperative and excludes Python consumers.

Parquet can win by using a small compressed sorted table with column projection
and predicate pushdown. Native provenance and replay have additional work and
semantics; those do not imply universal performance superiority. The prior cohort
analyzer was slower than simpler saved queries and independent checked Parquet,
so optional optimization stopped. The frozen HG001 richer policy removed 16 false
candidates but also two true candidates and failed practical gates, so the simple
baseline remained. Those decisions are evidence of choosing a release boundary,
not evidence that greater complexity is always bad.

A changed measurement engine would require exact-source native tests, independent
raw-field checks, stable semantics, unchanged answers where promised, verified
missingness and failure behavior, and an evaluated resource envelope on affected
workloads. This one-window case cannot establish all-platform/all-workload safety,
variant-calling accuracy, clinical usefulness, production adoption or neural
interface performance. Each of those needs its own evaluation.
