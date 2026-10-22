# Bounded diagnostic workload

`pilot-v1.json` fixes a synthetic three-member, 1,024-locus experiment before
measurement. Its selected coordinates span four native ownership partitions.
Distinct BAM/read-group/sample identities, deterministic read counts, excluded
reads, stored zeros and missing masks are generated locally on the selected
validation host. These members are synthetic observations, not participants.

Run after an exact-revision hosted build and correctness gates:

```sh
python research/cohort-reanalysis/benchmark.py \
  --binary /absolute/path/to/rosalind \
  --config research/cohort-reanalysis/workloads/pilot-v1.json \
  --output /absolute/path/to/new-pilot
```

The command never downloads data or builds a binary. It requires the wheels in
`requirements.lock`; the Arrow and pysam pins match the existing platform harness.
The focused GitHub workflow installs wheels only and builds Rust on the hosted
runner. Its identity report records the actual checked-out SHA, platform, compiler,
Cargo lock, binary digest/version and Python dependency versions. A prior local
binary cannot validate the changed native adapter.

B uses existing verified saved-query streams and a straightforward exact reducer.
C runs the research analyzer at both predefined windows. D independently counts
the same generated BAMs using explicit one-base CIGAR/filter logic, materializes
uint64 Parquet with row-group statistics, then uses file pruning, predicate and
projection pushdown, 1,024-row batches, one-batch/fragment readahead, and no threaded
readahead. Its two reported modes distinguish performance-only queries from
consumed-file SHA256 checks before and after reduction. These checks are useful
integrity costs, but do not supply native receipts or every analyzer guarantee.
No native or analyzer output is used to populate D's cache.

Every worker's outer wall time includes interpreter startup, planning/imports,
verification where applicable, decoding, reduction and output publication. Setup
includes extraction, materialization and initial verification, and is charged to
all total-time prefixes. The pilot measures setup once per representation and
labels this limitation. Three rotations counterbalance the principal query paths;
all trials and failures are retained. Requests are revealed sequentially within
each K prefix; no implementation receives a coalesced future query. The prefixes
are K=1,2,3. No raw-input A path is measured here, so this pilot cannot establish a
raw-versus-saved amortization crossover.

Process-tree RSS is the sampled simultaneous sum for a worker and its active
children at 10 ms intervals. It can miss brief peaks, excludes the measurement
driver, and is not an OS-enforced cap. Native cooperative budgets are reported
separately. Storage records logical raw/cache/query-output bytes; it is not a
sampled peak temporary-disk metric. The entire retained pilot is refused above
512 MiB. There is no OS cache eviction and no cold-cache or p99 claim.

Exact canonical answers must agree across B, both C windows, both D modes and all
repetitions before the pilot passes. `results.json` records trial order, all
commands, failures, setup costs, medians/ranges, memory observations and identities.
Retained native query receipts and analyzer manifests give their actual narrower
scopes. The harness deliberately does not turn an unrun configuration into a
measured advantage or select an optimization automatically.

`held-back-v1.json` is reserved for an engineering generalization check after the
algorithm and protocol are frozen. Do not tune against it. Real-read HG002
measurement remains a separate Tier B run requiring the locked local resources;
this entirely synthetic pilot cannot replace it.
