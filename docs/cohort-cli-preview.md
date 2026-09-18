# Cohort candidate reanalysis: CLI preview

This guide describes the **unpublished 0.6 development preview**. Use a binary built
from the cohort development branch; the public stable release and 0.5 candidate do
not include these commands. Check `rosalind --version` and `rosalind cohort --help`.
The hidden binary bridge is not a supported Rust cohort SDK.

With the [source build prerequisites](analyzer-sdk.md#build-prerequisites) installed:

```sh
git clone --branch codex/cohort-preview https://github.com/logannye/rosalind.git rosalind-cohort
cd rosalind-cohort
cargo build --locked --release --bin rosalind
export PATH="$PWD/target/release:$PATH"
rosalind --version
rosalind cohort --help
```

The [three-sample tutorial](../examples/cohort-reanalysis/README.md) supplies
redistributable data, an independent arithmetic oracle, a readable report, a
Python example and measurements of import, storage, verification and repeated queries.

A cohort stores intact copies of verified evidence datasets and immutable sample
snapshots. It answers a new candidate-SNV question across selected samples, reports
which observations were actually measured, and can explicitly extract missing loci
when original local sources are supplied. It counts read observations; it does not
call genotypes, infer biological pairing, or provide clinical interpretation.

## Import saved evidence

First create portable datasets with the [researcher workflow](researcher-quickstart.md)
and [evidence reuse guide](reuse-quickstart.md). Each member needs a named, single-sample
dataset. A member ID and optional group, subject, or timepoint are your assertions;
Rosalind does not verify biological identity from those labels.

Create a UTF-8 **tab-separated** table such as `members.tsv`:

```text
id	manifest	group
sample-A	datasets/A/evidence-dataset.manifest.json	study
sample-B	datasets/B/evidence-dataset.manifest.json	study
```

Paths are relative to the table's actual directory, or absolute. Required columns
are `id` and `manifest`; optional columns are `group`, `subject`, and `timepoint`,
in any order. Empty optional cells mean absent metadata. To import multiple disjoint
leaf datasets for one member, repeat its ID with the same metadata and a different
manifest. Duplicate member/manifest rows, conflicting metadata, overlapping genomic
ownership, and one source/sample imported under multiple member IDs are refused.
A header-only table is refused. IDs must be nonempty, at most 256 UTF-8 bytes, and
contain no control characters. Each table row is bounded to 8,192 bytes; the default
table envelope is 8 MiB.

```sh
rosalind cohort create --cohort research-cohort --members members.tsv --plan
rosalind cohort create --cohort research-cohort --members members.tsv > created.json
```

The plan reads metadata; it does not verify payload bytes or publish anything.
Creation copies and fully verifies each declared inventory before publishing the
snapshot. It prints JSON with `snapshot_id`, a 64-character content hash. Use this
exact ID in the following commands. Reimporting identical content is idempotent.
`--parent ID` records ancestry within the same cohort; it does not automatically
inherit members omitted from the supplied table.

```sh
SNAPSHOT=$(python3 -c 'import json; print(json.load(open("created.json"))["snapshot_id"])')
rosalind cohort inspect --cohort research-cohort --snapshot "$SNAPSHOT"
rosalind cohort verify --cohort research-cohort --snapshot "$SNAPSHOT"
```

`inspect` validates snapshot ancestry, object metadata, and ownership without
reading evidence rows. `verify` additionally scans all current snapshot payloads.
It does not scan unrelated ancestor payloads or reopen original BAM/CRAM/reference
files. Keep every cohort file immutable. Copying the complete cohort directory
preserves its usability; always identify the intended snapshot explicitly.

## Ask a candidate question

Supply A/C/G/T SNVs in local VCF, VCF.gz, or BCF. Candidate files need no index.
Repeated loci and ALT alleles are normalized into canonical reference/position/ALT
order. A new ALT at a stored locus can use the stored A/C/G/T counts.

```sh
rosalind cohort extract --cohort research-cohort --snapshot "$SNAPSHOT" \
  --sites candidates.vcf --plan
```

The JSON plan reports `status: ready` or `blocked`, selected members, covered and
missing loci, required field groups, structured incompatibilities, and resource
reservations. A successfully produced blocked plan exits successfully so tools can
inspect its issues. Actual materialization refuses those issues. Planning does not
prove stored REF values or payload integrity; readers verify those when consumed.

The default is strict coverage. If any requested sample/locus is absent, no result
or successful receipt is published. An explicit partial query preserves absence:

```sh
rosalind cohort extract --cohort research-cohort --snapshot "$SNAPSHOT" \
  --sites candidates.vcf --missing partial --format tsv --output evidence.tsv
rosalind cohort summarize --cohort research-cohort --snapshot "$SNAPSHOT" \
  --sites candidates.vcf --missing partial --min-callable-depth 10 \
  --format tsv --output summary.tsv
```

Use `--member sample-A --member sample-B` for a subset; omitted selection means all
members. The default projection is `depths,alleles`. Additional stored groups can be
requested with `--fields depths,alleles,strands`, or `all-supported`. Missing groups,
unequal scientific profiles, incompatible reference content, or unsupported sample
scope are errors even in partial mode. A different stored mask is acceptable when
it supplies every requested group.

Extracted keys are `sample_id`, `contig`, **one-based** `pos`, `ref`, `alt`, and
`status`. Scalar evidence uses exact unsigned 64-bit integers. Arrow uses actual
nulls for unmeasured evidence; TSV uses `.`. Optional fixed-size metric lists in TSV
are dense comma-separated integer values; a missing whole list is `.`.

| Observation | Status | Callable depth / ALT | Depth eligible | ALT fraction |
|---|---|---|---|---|
| Locus absent from stored evidence | `unmeasured` | null / null | null | undefined |
| Stored row with no callable reads | `observed` | 0 / 0 | false | undefined |
| Stored low-depth row | `observed` | exact measured counts | false | ALT/depth when depth > 0 |
| Stored row passing the depth screen | `observed` | exact measured counts | true | ALT/depth |

Fractions are separate exact numerator and denominator columns. A zero denominator
produces nulls, never 0/0 or a fabricated zero fraction. Summaries report
`n_requested`, `n_observed`, `n_depth_eligible`, and `n_alt_supported`, plus observed
and eligible read-count totals. `n_alt_supported` counts depth-eligible members with
at least one ALT read. The default depth threshold 10 is a technical screen, not
confidence; sample support proportions are not population allele frequencies.

Run another candidate file against the same snapshot using the same commands.
Saved-only extraction and summaries never implicitly reopen alignments or reference
files. Arrow IPC uses canonical 1,024-row output batches; TSV and Arrow bytes do not
depend on admitted execution window width. Output order is member then candidate
for extraction, candidate for summaries. Empty selections produce schema/header-only
outputs. Existing output/receipt paths are refused unless `--force` is explicit.

## Explicitly fill missing loci

Extension creates missing-only datasets and a **new snapshot**. The original
snapshot remains unchanged. Create a source table with the exact header below,
one recorded source role per row. Include exactly the affected selected members.
All scientific source roles recorded in their existing datasets must be mapped;
this commonly includes the following BAM roles:

```text
id	role	path
sample-B	alignments	raw/sample-B.bam
sample-B	alignment-index	raw/sample-B.bam.bai
sample-B	reference	raw/reference.fa
sample-B	reference-fai	raw/reference.fa.fai
```

Paths are table-relative or absolute. CRAM or reference-pack datasets may record
additional/different reference roles; use those exact roles from the portable
dataset descriptor. Paths may relocate, but source byte identities must match.
There is no source discovery and no silent field backfill at already stored loci.

```sh
rosalind cohort extend --cohort research-cohort --snapshot "$SNAPSHOT" \
  --sites second-candidates.vcf --sources sources.tsv --plan
mkdir extension-work
rosalind cohort extend --cohort research-cohort --snapshot "$SNAPSHOT" \
  --sites second-candidates.vcf --sources sources.tsv --work-dir extension-work \
  > extended.json
```

The extension plan lists affected members and whether the supplied member mapping
set matches, without opening original sources. It does not certify source identities;
execution hashes and validates them. The scratch directory must already exist
outside the cohort. Temporary extracted data are cleaned up. If nothing is missing,
use a source table containing only its header; the operation returns the existing
snapshot without opening raw sources.

The result JSON reports the new `snapshot_id`, `changed`, and per-member retained
and computed loci. Native indexed read visits, CRAM full-file validation records,
source hashing, and existing-object byte verification are separate measurements.
Already stored evidence rows are not decoded by extension, although existing object
bytes are rehashed for publication integrity. CRAM can still require a complete raw
source validation pass; missing-only extraction is not a promise of zero other I/O.
Any failure preserves existing snapshots. Query the returned snapshot ID with strict
coverage to confirm the newly covered question.

## Resources, verification, and replay

Add `--memory-budget-mb 512 --enforce` to admit and monitor the complete process.
Without `--enforce`, a declared budget is recorded rather than guaranteed.
`--require-os-limit` additionally requires an existing compatible Linux cgroup-v2
limit. Scheduler resource requests alone do not establish OS enforcement.
`--tile-bases` changes execution windows; it cannot change scientific settings.
Variant header/record limits are cooperative decoder envelopes, not protection
against every transient native allocation.

Successful extracted and summarized files have `<output>.manifest.json` receipts,
or the explicit `--manifest` path. Receipts bind the snapshot, normalized query,
member selection, consumed dataset files, reducer settings, and output hashes.

```sh
rosalind verify --manifest evidence.tsv.manifest.json
rosalind reproduce --manifest evidence.tsv.manifest.json \
  --inputs /path/containing/relocated-cohort-and-candidate-file \
  --binary /absolute/path/to/the-matching/rosalind
```

Retain the complete relocated cohort directory structure and the original candidate
file for file-based queries. Replay supports cohort extraction and summaries in
Arrow/TSV. It does not claim Parquet replay, cohort-extension replay, or independent
biological validation. Receipts establish integrity and reproducibility within their
recorded scope; they do not establish clinical correctness or external adoption.

The [synthetic fixture](../examples/cohort-reanalysis/README.md) supplies authored
multi-sample expectations. Binary integration coverage in
[`tests/cohort_cli.rs`](../tests/cohort_cli.rs) exercises exact missingness/counts,
compressed variants, explicit extension, portability, and relocated replay. These
are engineering checks, not design-partner sessions or performance benchmarks.
