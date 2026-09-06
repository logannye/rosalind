# Portable and reusable evidence

The current development API can persist exact evidence independently of the
original alignments. It supports offline queries, target QC, serial partial-overlap
reuse, and bounded Parquet export. These additions are not yet a published release.

## Create and move a dataset

```sh
rosalind analyze evidence --reference reference.fa --alignments reads.bam \
  --regions panel.bed --fields depths,alleles,quality-sums,allele-quality \
  --cache-dir evidence-cache --format arrow-ipc -o evidence.arrow
```

The primary receipt records the portable receipt path in
`measurements.execution.evidence_dataset_manifest`. Its containing directory has
`evidence-dataset.manifest.json`, `dataset.descriptor.json`, and canonical partition
subdirectories. Copy the entire directory together. The older
`dataset.manifest.json` remains available for existing cache/lookup consumers.

New portable datasets use a versioned namespace minted by `VerifiedInputSession`
from actual source hashes. The older API accepting a caller-supplied digest remains
available but does not itself establish verified source identity. Old cache files
are never silently relabeled as portable datasets.

```sh
rosalind dataset inspect --dataset portable/evidence-dataset.manifest.json
rosalind dataset verify --dataset portable/evidence-dataset.manifest.json
```

`inspect` verifies bounded metadata. `verify` also hashes and validates every
partition's evidence, field mask, row order, ownership, and expected positions.
Neither command opens the original BAM/CRAM/reference: those recorded identities
remain provenance claims. Locally verified dataset files and original-source claims
are distinct. Receipts detect changed bytes; they do not establish authorship or
independently prove biological accuracy.

## Query without alignment decoding

```sh
rosalind dataset extract --dataset portable/evidence-dataset.manifest.json \
  --regions subset.bed --fields depths,alleles --format arrow-ipc -o subset.arrow
rosalind dataset panel-qc --dataset portable/evidence-dataset.manifest.json \
  --regions subset.bed --min-callable-depth 10 -o panel.tsv
rosalind dataset extract --dataset portable/evidence-dataset.manifest.json \
  --fields depths --memory-budget-mb 256 --plan
```

Omitting selection means all stored positions. `--sites` also accepts SNV VCF,
VCF.gz, and BCF. Queries validate requested reference alleles and preserve their
normalized ALT annotations; a new ALT can use the already stored A/C/G/T counts.
A position absent from the stored selection causes refusal. A stored position with
zero callable reads remains a genuine zero-depth row. Both cases are explicit.

Stored fields must cover the requested groups. Quality thresholds, flag filters,
sample scope, reference mode, counting unit, and scientific semantics belong to the
original extraction. An aggregate cannot reconstruct observations already excluded
by a different filter. Panel summaries retain the same target denominators and
quality definitions as live extraction.

Output batches have at most 1,024 rows. A projection reserves the physical source
decoder as well as its smaller output arrays. Query metadata grows with selected
sites, intervals, and partition inventory; memory is not independent of these
cardinalities. Plans include actual metadata, decoder, projection, consumer,
encoding, and finalization costs. Library callers receive local budget checks even
without a global governor. OS limits remain separate from cooperative checks.

The default maximum descriptor/receipt size is 32 MiB. For a large inventory,
`--max-dataset-metadata-bytes BYTES` raises the accepted envelope; it does not waive
memory admission. This option is available during creation, reuse, and offline
queries. Python `open_dataset(..., max_dataset_metadata_bytes=...)` records the
same limit. An oversized envelope or insufficient working set causes refusal.

## Reuse a subset while filling missing loci

```sh
rosalind analyze evidence --reference reference.fa --alignments reads.bam \
  --regions expanded-panel.bed --fields depths,alleles \
  --reuse-dataset portable/evidence-dataset.manifest.json \
  --workers 1 --memory-budget-mb 256 --format arrow-ipc -o expanded.arrow
```

This rehashes the current alignment, index, and reference inputs once, verifies
exact compatibility with the persisted source, and computes only missing loci.
The receipt separates `execution.reused_loci` and `execution.computed_loci` from
biological evidence counters. Successful output bytes equal fresh extraction with
the same scientific settings. A fully covered request performs no alignment-record
decoding, although input verification and preflight still run.

The first reuse implementation keeps one projected canonical partition plus a
source decoder and native extraction state. It supports one worker and refuses
combination with `--cache-dir`, `--resume`, or `--annotated-variants`. Existing exact
cache resume and parallel extraction remain separate paths. An expanded result is
a materialized Arrow/TSV artifact; this command does not publish an expanded cache.

## Python, R, and SQL

```python
from rosalind import open_dataset

dataset = open_dataset("portable/evidence-dataset.manifest.json")
with dataset.batches(regions="subset.bed", fields=["depths", "alleles"]) as run:
    for batch in run:
        consume(batch)
    print(run.result.manifest_path)

result = dataset.materialize("subset.arrow", regions="subset.bed")
dataset.panel_qc("subset.bed", "panel.tsv")
export = dataset.export_parquet("evidence_export", fields=["depths", "alleles"])
```

A handle is lazy. Every native invocation verifies metadata and the partitions it
consumes. Stream exhaustion finalizes a run result; early cancellation does not.
Python-retained batches and downstream SQL/R memory are outside the native budget.
The [executable language examples](../examples/persisted-evidence/README.md) show
R character-column TSV import and DuckDB querying preserved unsigned integers.

Parquet exports atomically publish a new directory after data, receipts, input
mutation checks, and resource checks succeed. They refuse existing directories.
Default files contain at most 16,384 rows, row groups at most 1,024. Compression and
dictionaries are disabled; a bounded footer avoids retaining every row group's
metadata for an arbitrarily large single file. Empty output still contains a
schema-bearing Parquet file. Integer scalar and allele-quality list columns retain
uint64 values, including values above floating-point's exact integer range.

## Verification and replay

Materialized Arrow/TSV query receipts bind the parent, descriptor, consumed
partitions, selection inputs, scientific settings, and actual output bytes:

```sh
rosalind verify --manifest subset.arrow.manifest.json
rosalind reproduce --manifest subset.arrow.manifest.json --inputs relocated-inputs
```

`--inputs` searches the supplied tree for expected hashes using a bounded directory
stack; it does not follow directory symlinks. Keep portable dataset directories
intact so their relative partition paths remain usable. Include external selection
files in that input tree. Replay compares physical output bytes.

Parquet receipts bind physical output files, normalized query, verified lineage,
and export settings. Native byte replay for the directory export is explicitly
unsupported in this version; use materialized Arrow/TSV for replay. Historical
receipt schemas 1–5 and pre-projection evidence goldens remain unchanged.
