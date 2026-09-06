# Nextflow integration

The maintained exact-evidence workflow is `examples/evidence/main.nf`. It plans
the actual panel analyzer, derives scheduler memory from predicted bytes plus a
margin, extracts panel TSV and position Arrow in one pass, and verifies their
shared receipt with all inputs staged. Use Nextflow 25.04.8 and Java 21:

```sh
nextflow run integrations/nextflow/examples/evidence/main.nf \
  -c integrations/nextflow/examples/evidence/nextflow.config \
  --reference reference.rref --bam sample.bam --index sample.bam.bai \
  --targets targets.bed --outdir results/evidence \
  --image 'ghcr.io/logannye/rosalind@sha256:REPLACE_WITH_PUBLISHED_DIGEST'
```

Replace the placeholder only with a published immutable digest. For local source
validation, put `rosalind` on PATH and use `nextflow.test.config` instead; it
disables Docker. The workflow requests scheduler memory but does not claim a
verified OS cap. The equivalent Snakemake workflow is `../snakemake/Evidence.smk`.
No container/Slurm execution is implied by a local source run.

## Legacy feature shards

The versioned module exposes doctor, plan, deterministic shard analysis, canonical
merge, verify, and optional reproduction processes. Set `image` to the immutable
`ghcr.io/logannye/rosalind@sha256:…` value from `integrations/container-lock.json`;
mutable tags are intentionally unsupported in production examples. Fan out shard
indices with `Channel.of(0..<shard_count).flatten()` and collect all artifact/receipt
tuples before `ROSALIND_MERGE`.

The local example joins each emitted `plan.json` to its shard and derives the
scheduler request from `predicted_peak_rss_bytes`, rounded up to MiB plus
`scheduler_margin_mb`. Run it after replacing the placeholder with the published
digest:

```sh
nextflow run integrations/nextflow/examples/local/main.nf \
  -c integrations/nextflow/examples/local/nextflow.config \
  --reference reference.rref --bam sample.bam --bai sample.bam.bai
```

For Slurm, add
`-c integrations/nextflow/examples/slurm/nextflow.config`. Reproduction remains a
separate optional process because it reruns the complete analysis. Verification
tuples include the evidence inputs named by the receipt; staging only an artifact
and its JSON sidecar is insufficient when `verify` must re-hash parent shards.
