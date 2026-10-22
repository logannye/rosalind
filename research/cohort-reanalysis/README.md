# Exact saved-evidence candidate reanalysis

This research analyzer summarizes candidate SNVs across explicitly named specimens,
using one existing portable evidence leaf per specimen. It answers how many
specimens have a stored observation, meet a depth threshold, and have eligible ALT
support. Missing observations remain distinct from stored zero-depth observations.
It never reopens original alignments or references to fill missing coverage.

**Status:** implementation and pure Python tests are available. The new native
coverage adapter, end-to-end analyzer, failure/replay paths, and diagnostic pilot
require a hosted build of the exact changed revision before integration or
performance claims. No new benchmark result or representative cohort performance
is established by this README. See [IMPLEMENTATION.md](IMPLEMENTATION.md) for gates.
This directory does not establish a public `cohort` command or `open_cohort` API,
and does not complete roadmap C02–C11.

## Supported questions and information boundaries

All successful comparisons retain exact read observations. Overlapping mates
count separately. No UMI/fragment consensus, genotype, haplotype, disease
classification, or clinical interpretation is inferred.

| Question | Retained fields required | Fixed information and limitations |
|---|---|---|
| Stored observation and callable-depth screen | `depths`; this analyzer always requires `depths,alleles` | Same stored locus mask, read unit, sample scope and extraction filters |
| Different positive depth threshold | `depths,alleles` | Changes the downstream screen; never changes extraction MAPQ/BQ settings |
| New A/C/G/T ALT at a stored position | `depths,alleles` | Uses retained A/C/G/T counts; observed REF must match saved evidence |
| Several ALT alleles at the same position | `depths,alleles` | Each ALT row repeats the locus denominator; repeated denominators are not combined |
| Query positions absent from one or every leaf | Stored-coordinate coverage metadata | `strict` refuses; `partial` retains observation denominators and reports missing loci |
| Lower MAPQ/BQ, different flag exclusions, fragment or UMI counts | Not recoverable by this analyzer | Discarded reads and required relationships cannot be reconstructed from counts |
| Indels, symbolic alleles, haplotypes, linkage or genotype inference | Unsupported | Native SNV normalization rejects unsupported candidate alleles |

Saved evidence is sufficient only for the declared questions under its immutable
extraction contract. It is not a universal substitute for raw reads.

## Inputs and comparison contract

`members.json` has a versioned, explicitly enumerated shape:

```json
{
  "version": 1,
  "members": [
    {"specimen_id": "specimen-a", "dataset": "portable/a/evidence-dataset.manifest.json"},
    {"specimen_id": "specimen-b", "dataset": "portable/b/evidence-dataset.manifest.json"}
  ]
}
```

Relative paths resolve against the manifest's directory. Each path must identify
an existing leaf manifest file; directories are not searched. Unique nonempty IDs
sort lexically and are limited to 256 UTF-8 bytes without control characters.
Optional `group`, `subject`, and `timepoint` strings obey the same bound and remain
user assertions. They establish no pairing or biological independence.

Native `inspect` validates each descriptor. Comparison then requires named sample
scope, complete resolved read-group metadata, equal explicit read/filter semantics,
equal ordered contig dictionaries, and effective reference source role, byte hash,
and byte size. `reference` takes precedence over `cram-reference`. Descriptor v1
does not record a format enum: exact original-file byte identity conservatively
preserves representation without guessing from a suffix or equating FASTA and
`.rref` sequence content. Original source paths are retained as provenance only.

Physical masks may differ if every member supplies `depths,alleles`, with supported
native schema/mask versions. Legacy `all` is mask 63; `all-supported` is mask 127.
Different alignment bytes and sample labels are allowed. The same alignment
content and complete named scope entered under aliases is rejected, even if its
index or projection differs. The existing native `compatibility_blake3` reuse key
is unchanged. The separate research `comparison_sha256` binds shared measurement
semantics rather than declaring that different specimens have the same source.

`--sites` is a byte-bounded uncompressed SNV VCF handled by the native parser.
The research CLI refuses compressed VCF/BCF to bound expanded selection input;
the underlying native command retains its broader format support. Native coverage
owns coordinate conversion and canonicalization. Internally normalized sites are
zero-based; output `pos` is one-based. Identical candidates collapse, ALT alleles
are unioned and sorted, and contradictory REF annotations are rejected. A VCF with
no candidate records is a valid empty request. Python never implements a second
permissive VCF/BED input parser.

## Commands

Use an already built executable for the exact changed revision. Required Python
runtime dependencies are `pyarrow` and `blake3`; `psutil` enables sampled process
tree memory. The synthetic native fixture and pilot also require `pysam`. Do not
install the project from source as a shortcut: that can compile Rust. Local
compilation remains opt-in; run native build and validation on hosted CI.

```sh
python research/cohort-reanalysis/analyze.py \
  --binary /absolute/path/to/built/rosalind \
  --members /absolute/path/to/members.json \
  --sites /absolute/path/to/candidates.vcf \
  --min-callable-depth 5 --min-callable-depth 10 --min-callable-depth 20 \
  --missingness partial --window-size 256 \
  --native-memory-budget-mb 512 --parent-memory-budget-mb 128 \
  --output /absolute/path/to/new-run
```

Default missingness is `strict`; default depth threshold is 10. Repeated positive
thresholds are deduplicated and sorted. The parent of `--output` must exist, and
the output path must be new. The selected executable must match the installed
Python wrapper version and implement `dataset coverage`; an older binary cannot
validate the changed Rust source.

Add `--plan` to publish metadata, normalized selections, coverage gaps and admitted
resource assumptions without reading partition bodies or publishing a completed
aggregation manifest. A plan does not validate REF against rows or certify payload
integrity. `strict` refuses coverage gaps in preflight. `partial` queries only
each member's covered intersection through the existing strict native reader.

```sh
python research/cohort-reanalysis/analyze.py verify /absolute/path/to/new-run
python research/cohort-reanalysis/analyze.py replay /absolute/path/to/new-run \
  --members /absolute/path/to/relocated-members.json \
  --binary /absolute/path/to/built/rosalind \
  --window-size 128 --output /absolute/path/to/new-replay
```

`verify` checks the research manifest, retained output bytes, and consistency of
the scientific document's identity. It does not reopen original sources, rerun
the reducer, authenticate authorship, or certify untouched partition bodies.
`replay` first verifies the retained run, then reruns with an explicit member map
and saved normalized candidates. It compares canonical summary bytes and
scientific metadata identity, excluding timing and execution settings. The
relocated leaves must preserve their native manifests/descriptors. A separate
`<new-replay>.replay.json` records agreement or divergence; an independently
completed divergent run is retained. Replay cannot certify biological truth.

## Scientific outputs

`summary.jsonl` is canonical long-form JSON: reference-dictionary order, position,
ALT, then threshold. Count values are JSON integers with explicit uint64 bounds;
consumers must preserve those integers rather than converting them to binary64.

| Columns | Meaning |
|---|---|
| `contig,pos,ref,alt,min_callable_depth` | Candidate and positive downstream threshold |
| `n_requested,n_observed` | Requested specimens and specimens with a stored row |
| `n_depth_eligible,n_alt_supported` | Observed depth at threshold; eligible ALT count at least one |
| `total_callable_depth,total_alt_count` | Additive observations across the observed set |
| `eligible_callable_depth,eligible_alt_count` | Additive observations restricted to eligible specimens |
| `alt_fraction_numerator,alt_fraction_denominator` | Exact ALT/callable pair across observed specimens |
| `eligible_alt_fraction_numerator,eligible_alt_fraction_denominator` | Exact ALT/callable pair restricted to eligible specimens |
| `reference_validation` | `checked_on_observed_rows` or `unavailable` when no member stores the locus |

Both members of a fraction pair are null when its denominator is zero. A stored
zero-depth row contributes one to `n_observed`, zeros to totals and no eligibility.
A wholly absent candidate contributes zero to additive totals over the empty set,
has no defined fractions, and retains REF as an unchecked query assertion. These
summary zeros do not turn absent member measurements into observed zeros. A zero
ALT count is not a reference genotype. Any uint64 input/sum overflow fails the run.

`scientific.json` binds comparison/query identities, members, missingness,
thresholds and fraction policy. `plan.json`, `coverage.jsonl`, per-member metadata
and coverage files retain the input and admission explanation. `children/` retains
the completed native stdout-query receipts and exact window VCFs.
`input-identities.json`, `metrics.json` and `run-manifest.json` bind input checks,
measurements, source/binary/dependency identities and output bytes. The wrapper
manifest is a separate research format, not a native receipt for the outer sum.

## Bounded execution and lifecycle

Execution is serial: one admitted locus window, one member and one active native
decoder, with at most 1,024 projected rows per batch. Only window counters and a
current-member observation mask are accumulated; the cohort is not materialized
as a dense member-by-locus matrix. Input/member/descriptor metadata and receipts
also occupy memory and are included in admission or explicitly refused.

| Current limit | Bound |
|---|---|
| Specimens / normalized ALT rows | 128 / 100,000 |
| Distinct depth thresholds | 8 |
| Distinct loci per requested window | 1–100,000; default 256 |
| Native query children / consumed input identities | 4,096 / 16,384 |
| Summary rows | 1,000,000 ceiling, also constrained by ALT and threshold limits |
| Retained run files | 256 MiB before the bounded manifest; explicit file inventory cap |
| Individual input/metadata response | At most 4 MiB and tightened by parent budget |
| Parent fixed reserve | 64 MiB, plus modeled metadata, candidates, counters and receipts |

The model charges JSON/Python amplification, physical decoding plans, projected
batches, output inventories and finalization. Requested budgets are cooperative
admission limits, not OS-enforced RSS caps. Native physical partitions may be wider
than the requested columns or windows. Small windows can repeat partition hashing,
decoding and process launches; no performance benefit is assumed.

The union of every consumed child's dependencies is checked again before
publication, including files used by earlier children. Inputs are expected to
remain immutable; this is not snapshot isolation or an adversarial race guarantee.
Create-new atomic publication is supported on Linux and macOS. Refusal, corruption,
overflow, cancellation or native failure leaves the final output unpublished;
hidden staging may retain a clearly failed diagnostic record without a completed
run manifest. Existing destinations are never overwritten.

`metrics.json` labels process-tree RSS as sampled, with its sample interval and
coverage. Internal phase times may overlap and must not be added. Internal sampled
RSS excludes final publication; the outer benchmark measures the full worker
lifetime including publication. Bytes hashed by the parent, final union rechecks,
native invocations, partitions and repeated decodes are reported in their actual
scope. Unavailable metrics are not silently treated as measured zeros.

## Validation and diagnostic study

Pure Python checks require no native executable or build:

```sh
python -m unittest discover -s research/cohort-reanalysis/tests -p 'test_contract.py' -v
python -m unittest discover -s research/cohort-reanalysis/tests -p 'test_reduce.py' -v
python -m unittest discover -s research/cohort-reanalysis/tests -p 'test_lifecycle.py' -v
```

The reducer test reads the unchanged authored oracle at
[`examples/cohort-reanalysis/expected.json`](../../examples/cohort-reanalysis/expected.json).
Native integration tests require `ROSALIND_BINARY` explicitly and are skipped
without it. An exact-revision hosted build must pass integration and relevant
existing conformance/replay/CRAM gates before new end-to-end claims.

`benchmark.py --binary /absolute/path/to/built/rosalind --output /absolute/path/to/new-pilot`
uses the locked [`workloads/pilot-v1.json`](workloads/pilot-v1.json) diagnostic
protocol: three synthetic members, 1,024 spaced loci spanning several partitions,
three sequential requests, prefixes 1/2/3, two windows, and three rotated trial
orders. It compares existing saved queries (B), this analyzer (C), and bounded
PyArrow/Parquet scans (D). D has separately labeled performance-only and content-
checked modes; neither claims native receipt semantics. Initial native/cache setup,
all trials/failures, exact agreement, sampled memory and logical storage are retained.
No cache eviction, cold-cache claim, p99 estimate or general speedup is implied.

This pilot does not implement the final real-read fresh-extraction comparison (A),
prove an amortization crossover, establish a final resource target, or select an
optimization automatically. A useful next decision depends on observed full costs.
Preserve losses and stop extending the wrapper if the strongest fair baseline
reveals no meaningful operational advantage.
