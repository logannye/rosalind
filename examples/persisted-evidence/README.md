# Reuse the same evidence from Python, R, and SQL

Requires a build containing the portable dataset commands. These examples use
short-read DNA read counts and keep native integer counters exact. They do not
infer fragments, molecules, genotypes, or clinical callability.

First create evidence with the fields your later analyses need:

```sh
rosalind analyze evidence --reference reference.fa --alignments reads.bam \
  --regions panel.bed --fields depths,alleles,quality-sums,allele-quality \
  --cache-dir evidence_cache -o evidence.arrow --format arrow-ipc
```

The output receipt's `measurements.execution.evidence_dataset_manifest` names
`evidence-dataset.manifest.json` inside the cache namespace. Copy its entire
containing directory to move the dataset. The original BAM/reference are needed
for extraction and filling missing loci, but are not needed for offline queries.

```sh
rosalind dataset verify --dataset portable/evidence-dataset.manifest.json
rosalind dataset extract --dataset portable/evidence-dataset.manifest.json \
  --regions subset.bed --fields depths,alleles --format arrow-ipc -o subset.arrow
rosalind dataset panel-qc --dataset portable/evidence-dataset.manifest.json \
  --regions subset.bed -o panel.tsv
rosalind dataset export --dataset portable/evidence-dataset.manifest.json \
  --fields depths,alleles,allele-quality -o evidence_export
```

The query includes every selected position, including genuine zero-depth rows.
A position absent from the stored selection causes refusal. Field projection
cannot recover omitted groups or change the extraction's quality filters.

Python (`rosalind-bio`, `pyarrow`):

```sh
python query.py portable/evidence-dataset.manifest.json subset.bed evidence_export
```

Each native iterator batch contains at most 1,024 rows. Python may retain arbitrary
numbers of batches; that retained memory is outside the native budget. Use the
iterator as a context manager when stopping early.

R uses base R and the native executable:

```sh
Rscript query.R portable/evidence-dataset.manifest.json subset.bed subset.tsv
```

The adapter reads TSV columns as character vectors to preserve the full uint64
range. Do not silently convert counters or sums to R doubles. Use an integer or
decimal representation appropriate to the values when doing arithmetic.

Run `query.sql` in DuckDB from the directory containing `evidence_export`.
Parquet files preserve unsigned 64-bit columns and fixed-list allele arrays;
list entries are A,C,G,T. DuckDB's list indexing starts at 1. SQL engines may
materialize data; their memory settings are separate from the native export plan.

Parquet exports use new directories, at most 16,384 rows per file and 1,024 rows
per row group. The manifest verifies physical Parquet files and lineage. Native
byte replay currently supports the materialized Arrow/TSV queries; the Parquet
receipt explicitly reports unsupported replay.
