# Integrate evidence into a workflow

The maintained Nextflow and Snakemake examples run exact panel QC and position
evidence in one extraction, plan the actual consumer, and verify the resulting
artifact pair. They use current source/candidate capabilities, which are not in
public stable v0.1.0. Start with [installation and source status](../README.md) and
[the evidence semantics](SEMANTICS.md).

## Local source validation

Build the current CLI and put it on PATH from the repository root:

```sh
cargo build --locked --bin rosalind
export PATH="$PWD/target/debug:$PATH"
```

Prepare a `.rref` reference, coordinate-sorted indexed BAM and target BED. The
[Nextflow guide](../integrations/nextflow/README.md) specifies its tested tool versions.
Its local test configuration disables Docker:

```sh
nextflow run integrations/nextflow/examples/evidence/main.nf \
  -c integrations/nextflow/examples/evidence/nextflow.test.config \
  --reference reference.rref --bam sample.bam --index sample.bam.bai \
  --targets targets.bed --outdir results/evidence
```

For Snakemake, copy `integrations/snakemake/evidence.example.yaml` to
`my-evidence.yaml` and set absolute input paths. The local profile runs the CLI on
PATH; its image placeholder is not a published container:

```sh
snakemake --snakefile integrations/snakemake/Evidence.smk \
  --configfile my-evidence.yaml --profile integrations/snakemake/profiles/local
```

See the [Snakemake guide](../integrations/snakemake/README.md) for container/Slurm
configuration. Local execution is evidence for that local workflow only.

## Plan, schedule, extract and verify

Both workflows derive scheduler memory from `predicted_peak_rss_bytes`, rounded up
to MiB with an explicit scheduler margin. The analysis plan includes panel reducer
and encoding state. A scheduler request is not a checked OS cap; Linux cgroup
assurance must be established separately if required.

Stage the reference, alignment, index and selection files needed by the receipt,
as well as artifacts and sidecars. Keep them immutable throughout the run.
Verification rehashes dependencies; staging only the JSON receipt is insufficient.
Replay is a separate operation that repeats work. For supported CRAM input in a
custom workflow, account for whole-file validation even during `--plan`; native
allocation can precede its checkpoints. The maintained examples above use BAM.

Pin a matching tested binary and workflow version. Container deployments require
an actually published immutable image digest; placeholder digests in source are
not usable releases. Current evidence candidates are unpublished. See
[implementation status](implementation-status.md) before selecting an artifact.

## Legacy wrappers

The repository-root [GitHub action](../action.yml) supports legacy `variants`,
`features` and `analyze coverage`, not `analyze evidence` or portable dataset
queries. Its existing inputs/behavior are preserved. Likewise, the legacy feature
shard/merge examples have different scientific defaults from the evidence profile.
Use the maintained evidence workflows for new panel integrations; do not substitute
legacy top-level `plan` for an evidence-consumer plan.

Custom reducers can follow the [builder quickstart](builder-quickstart.md) and
reuse the [managed SDK](analyzer-sdk.md) from their workflow process. Retain their
explicit executable identity for replay and reserve their declared consumer state.
