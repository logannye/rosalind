# Snakemake integration

`Evidence.smk` is the maintained exact-evidence/panel workflow. It plans with the
actual analyzer, derives scheduler memory from predicted bytes plus an explicit
margin, writes panel TSV and position Arrow in one extraction, and verifies their
shared receipt. It is equivalent to the [Nextflow evidence example](../nextflow/).

Prepare a `.rref`, indexed BAM, and BED; copy `evidence.example.yaml` with absolute
input paths and the immutable published image digest. No production image digest
is claimed until `integrations/container-lock.json` is populated by publication.

```sh
snakemake --snakefile integrations/snakemake/Evidence.smk \
  --configfile my-evidence.yaml --profile integrations/snakemake/profiles/local
```

The local command uses `rosalind` on PATH. Add Snakemake container execution options
for a real pinned image, or use the Slurm profile with the matching executor plugin
and site configuration. A scheduler request is distinct from a checked cgroup
limit; the workflow does not claim OS enforcement. Keep inputs and receipts for
later replay. Existing `Snakefile` preserves legacy feature shard/merge behavior
for compatibility; its scientific profile is different.
