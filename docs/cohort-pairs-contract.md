# Explicit paired candidate comparisons (0.6 source preview)

`cohort compare-pairs` and Python `EvidenceCohort.compare_pairs()` compare read
observations for user-supplied SNVs across explicitly ordered pairs of saved
samples. This interface belongs to the separate 0.6 development line; it is not
part of the current 0.5 release. Engineering fixtures do not establish independent
partner use or clinical validity.

The input table has exactly these tab-separated columns, in this order:

```tsv
id	left	right
baseline-to-followup	baseline	followup
control-to-treated	control	treated
```

Pair IDs are unique, nonempty and at most 256 UTF-8 bytes, without controls. Both
member IDs must exist and differ. Members may participate in multiple pairs;
distinct pair IDs can request both directions. Pair table order is preserved, then
candidates are ordered by reference dictionary, position and ALT. Names, subjects,
groups and timestamps never imply pairing or verified biological identity.
Header-only tables explicitly request no comparisons. Optional `--member` flags
limit the permitted pair scope; only members actually referenced by pairs are
queried. An unknown member or pair outside this scope is an error.

```sh
rosalind cohort compare-pairs --cohort saved-cohort --snapshot "$SNAPSHOT" \
  --pairs pairs.tsv --sites candidates.vcf --plan
rosalind cohort compare-pairs --cohort saved-cohort --snapshot "$SNAPSHOT" \
  --pairs pairs.tsv --sites candidates.vcf --missing partial \
  --memory-budget-mb 512 --enforce --format tsv --output pairs.tsv.result
```

```python
from rosalind import open_cohort
cohort = open_cohort("saved-cohort", snapshot_id)
plan = cohort.plan(sites="candidates.vcf", operation="compare-pairs", pairs="pairs.tsv")
result = cohort.compare_pairs("comparison.arrow", pairs="pairs.tsv",
                              sites="candidates.vcf", missing="partial",
                              memory_budget_mb=512)
```

To use the exact difference from an Arrow row, parse its decimal strings with
Python's arbitrary-precision `int` and `Fraction`:

```python
from fractions import Fraction
import pyarrow as pa

with pa.ipc.open_stream(result.artifact_path) as reader:
    for batch in reader:  # canonical batches contain at most 1,024 rows
        for row in batch.to_pylist():
            if row["difference_numerator"] is None:
                difference = None  # an unmeasured side or an observed zero denominator
            else:
                difference = Fraction(int(row["difference_numerator"]),
                                      int(row["difference_denominator"]))
                if row["difference_negative"]:
                    difference = -difference
            print(row["pair_id"], row["pos"], difference, row["both_depth_eligible"])
```

The [authored paired example](../examples/cohort-reanalysis/README.md#explicit-paired-follow-up-paired-development-branch)
uses A→B with these results:

| Position | Left ALT/depth | Right ALT/depth | Right-minus-left | Both depth-eligible |
|---|---|---|---|---|
| 10 | 4/10 | 1/9 | `negative=true`, `26/90` (−13/45) | false: right depth is below 10 |
| 20 | 0/0, observed | 4/12 | null: left denominator is zero | false |
| 40 | 2/12 | unmeasured | null: right was not measured | null |

The reverse B→A comparison flips the nonzero sign. A low-depth fraction still
exists mathematically; its technical eligibility flag remains visible. These
asserted specimen IDs and an explicit table do not establish biological pairing.
The engineering preview is implemented separately from 0.5; independently
observed paired use and user value remain unestablished.

The scientific fields are `depths,alleles`. Each row preserves each side's
`status`, exact uint64 `callable_depth` and `alt_count`, `depth_eligible`,
`alt_supported` and observed ALT fraction numerator/denominator. A prefix of
`left_` or `right_` identifies the side. Default `min_callable_depth=10` is a
technical screen, not confidence. Strict coverage refuses absent loci before
publication. Explicit partial coverage emits null measurements for an unmeasured
side. Stored zero-depth observations retain their zero counts and observed
status. Missing required fields and incompatible comparisons are errors.

The direction is **right observed ALT fraction minus left observed ALT fraction**:

`(right_ALT × left_depth − left_ALT × right_depth) / (right_depth × left_depth)`

`difference_negative` is a nullable Boolean. `difference_numerator` is the
unsigned magnitude; `difference_denominator` is the original depth product,
without implicit reduction. Both are exact decimal strings in Arrow and TSV,
since products of uint64 counts can exceed signed int128. No floating-point
rounding is used. Zero has a nonnegative sign. With either side unmeasured or
zero depth, all three difference fields are null (`.` in TSV). Positive low-depth
observations still have defined arithmetic differences. `both_depth_eligible` is
null when either side is unmeasured; otherwise it is the conjunction of the two
eligibility flags.

For example, right 3/4 minus left 1/2 yields positive 2/8, even though both sides
fail the default depth screen. If M is uint64 maximum, right (M−2)/(M−1) minus
left (M−1)/M is exactly −1/[M(M−1)]. Swapping sides reverses the sign. These are
observed read fractions, not genotype changes, population frequencies or clinical
conclusions.

Execution traverses one pair at a time through bounded genomic windows, using
saved evidence only. Canonical 1,024-row encoding batches are independent of the
execution window or admitted budget. Pair metadata, current-window observations,
decoder state, output and receipt costs participate in the existing single managed
resource/cancellation lifetime. Default table/count envelopes are 8 MiB and
65,536 pairs, with a separate 64 MiB conservative pair/scope metadata envelope.
Envelopes and cooperative budgets do not imply an OS cap.

Receipts bind the snapshot, explicit ordered pair table and normalized pair
identities, direction, normalized candidate query, required fields, depth screen,
missingness policy, consumed datasets/partitions and result bytes. The pair table
is a guarded input throughout publication. Relocated Arrow/TSV results use the
existing `reproduce --inputs` workflow, retaining the pair table beside the
candidate file and copied cohort. Changing pair order, direction or content is a
changed scientific request. Parquet replay is not added.

See [cohort CLI and Python onboarding](cohort-cli-preview.md) for importing
samples, selecting candidates, verification, extension and replay. This preview
still needs independently observed paired research usage; no pairing, outcome or
adoption is inferred from authored tests.
