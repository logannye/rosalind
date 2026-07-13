# Product and delivery roadmap

Rosalind's product promise is:

> Build deterministic genomics analyses that fit the machine—and prove exactly what ran.

The primary user is a genomics research software or bioinformatics platform engineer
building short-read, per-locus analysis. ML/data, reproducibility, and workflow/HPC
engineers are secondary users. Clinical diagnostics, wet-lab turnkey calling,
long-read field work, and production-caller competition are not current targets.

## Wave 0 — truthful 0.4 stabilization

- Analyzer-platform positioning and scaffold-first onboarding.
- CLI default quality threshold 30, matching the documented simulated workload.
- One consolidated 0.4.0 changelog; internal checkpoints are not presented as releases.
- Internal release, clean-fingerprint, package, asset, and fresh-install checks.
- Partner and seven-day soak gates start at 0.5.0, not the stabilization release.

Publication still requires protected `rc`/`release` environments,
`CARGO_REGISTRY_TOKEN`, pinned actionlint/shellcheck, ordered crates publication,
platform assets, and fresh-install execution of every documented stable command.

## Wave 1 — measured truth

- Publish the pinned HG002 GIAB chr20 result even if poor.
- Compare equivalent Rosalind, pysam, and bcftools-based coverage/QC work on setup,
  RSS, time, size, repeated bytes, insufficient-cgroup behavior, and reproduction.
- Compare `.idx` and `.rref` construction with relevant established index builders.
- Link every performance or scientific claim to an executable evidence artifact.

## Wave 2 — 0.5 analyzer foundation

- `.rref`, `ReferenceProvider`, build/inspect/convert, and guided `.idx` compatibility.
- Region/BED selection, deterministic complete non-overlapping shards, and canonical
  merge with strict analyzer/reference/schema/parameter compatibility.
- Canonical Arrow IPC with a versioned schema and fixed 65,536-row batches.
- `rosalind-bio` Python wheels importing as `rosalind`, bundling the matching binary,
  and exposing bounded lazy PyArrow batch iteration.
- Coverage/QC and ML feature analyzers, modular CLI, updated external scaffold.

Stable 0.5 requires installable RCs exercised by an analyzer builder, workflow/HPC
user, and constrained/offline user. Contract-changing feedback creates another RC.

## Wave 3 — 0.6 workflow adoption

- Versioned Nextflow module, Snakemake wrapper, digest-pinned OCI image, and Action.
- Scheduler examples using `doctor` and `plan --json`.
- Local, Slurm, and Nextflow deterministic shard fan-out/fan-in.
- Hashed small microbial and human reference packs.
- CRAM/CRAI unless design partners make it a 0.5 blocker.

Evidence gates: one external real workflow, one prevented OOM or reduced allocation,
two teams repeating after 30 days, and one external analyzer repository.

## Conditional work

- **Long reads/field:** only after three teams identify it as their primary unmet need;
  then ARM, ploidy, indels, technology-specific models, and device measurements.
- **Generic wrapper:** only after three teams prefer governance around mature tools;
  execution remains tokenized, allowlisted, and explicit about unverifiable bytes.
- **External-memory search index:** only after three users need fixed-RAM Rosalind
  search-index construction, `.rref` is insufficient, established tools cannot
  meet the need, and a bounded spike finds a competitive point.
- **Production caller:** only with scientific ownership, full GIAB evaluation,
  indels/complex variants, ploidy, calibration, difficult regions, multiple
  technologies, and independent review.

## Success measures

- Median install-to-verified-analyzer receipt under 15 minutes; 80% partner
  completion without maintainer help.
- Three non-author teams, one external analyzer, one production-like workflow, and
  two 30-day repeat users.
- A documented prevented OOM/reduced allocation, an independent reproduction, and
  a causal diff identifying real drift.
- No command/version mismatch and no unsupported clinical, field, scientific, or
  competitive claim.
