# Core concepts

These concepts apply to the current exact evidence interface. The complete
definitions live in [scientific semantics](SEMANTICS.md); legacy analyzers have
their own profiles.

## Evidence is an observation summary

Rosalind answers what eligible aligned read bases show at selected positions.
The current unit is a **read base**, not an independent DNA molecule. Overlapping
mates count separately. Default filtering excludes duplicate-flagged reads, but
Rosalind does not discover duplicates, infer UMI families, correct sequencing
errors, or turn counts into a validated genotype call.

Three depths describe successive filters:

| Depth | Meaning |
|---|---|
| `prefilter_depth` | Matched aligned observations before profile filtering |
| `aligned_depth` | Observations passing the read-flag and mapping-quality filters |
| `callable_depth` | Observations also passing base-quality and A/C/G/T filters |

`prefilter_depth >= aligned_depth >= callable_depth`, and A+C+G+T equals
callable depth. Filter counters identify the first failing rule for each excluded
observation. A read may have other reasons for exclusion that do not receive an
additional counter.

## Selection and scheduling are different

A VCF selects SNV loci; its genotypes and FILTER labels do not change the read
filters. BED selects the union of its zero-based half-open intervals. Evidence
positions in TSV are one-based. The scientific selection fixes which loci must
appear, including positions with no eligible observations.

Execution tiles control how the engine visits those loci. A smaller memory budget
may cause smaller tiles and more indexed I/O. Successful outputs must remain
identical for the same inputs and scientific settings across admitted budgets
and worker counts. If resources do not fit, the engine refuses or fails explicitly;
it does not silently lower depth by downsampling.

## Zero is different from missing

| Situation | Meaning |
|---|---|
| Stored locus with `callable_depth = 0` | The position was processed, but no observation passed the profile |
| Positive depth below a chosen threshold | Some evidence exists, but it does not meet that technical screen |
| Requested locus absent from a saved dataset | The dataset cannot answer; the query refuses |
| Requested field group was not stored | Its values are unavailable; projection cannot reconstruct them |
| No ALT reads at a covered position | No eligible ALT observations under this profile; not proof of biological absence |

Panel summaries use the full original target length, including uncovered bases.
This prevents an apparently good mean or breadth from silently discarding the
hard-to-measure parts of a target.

## A saved dataset has a scientific profile

Profiles include reference identity, sample scope, quality thresholds, flag
filters, counting unit, and schema. Stored fields are selected at extraction.
Queries can request a subset of stored positions and fields, or a different ALT
at a stored SNV locus using its A/C/G/T counts. They cannot lower an extraction
threshold and recover excluded reads.

For example, storing `depths,alleles,quality-sums` supports candidate counts and
panel QC. It does not support later strand analysis unless `strands` was also
stored. [Dataset reuse](reusable-evidence.md) can fill missing positions from
verified source inputs when the scientific profile is compatible.

Sample identity comes from read-group metadata. Ambiguous mixtures need explicit
`--sample NAME` or deliberate `--pool-samples`. Pooling is recorded as pooling;
it does not become a named individual's evidence.

## Resource admission is not an operating-system cap

The planner accounts for the native engine, input decoding, consumer, encoding,
and metadata. Cooperative checks and sampled RSS cannot prevent every native
allocation before it occurs. Linux cgroup checks provide a separate assurance
when an appropriate limit is already configured. CRAM has additional envelope
admission and whole-file validation costs.

The Python iterator yields bounded batches, but Python, R, and SQL consumers own
the memory they retain. Collecting every batch defeats bounded consumer memory.
See [benchmarks and limitations](benchmarks-and-limitations.md) for measurements.

## Receipts connect results to their inputs

A receipt records content identities, scientific settings, producer identity,
and outcomes. Verification checks bytes and recorded dependencies. Replay reruns
the recipe using the matching producer. Portable dataset verification works
without the original alignments, so their recorded identities remain provenance
claims rather than a fresh inspection of those files.

Receipts do not establish biological validity or authorship. A valid hash cannot
make an inappropriate scientific profile appropriate. The
[trust guide](receipts-and-trust.md) explains these boundaries.

Try these ideas in the [researcher](researcher-quickstart.md) and
[reuse](reuse-quickstart.md) quickstarts.
