# HG002 GIAB v5.0q credibility benchmark

This opt-in workflow uses the current HG002 v5.0q small-variant benchmark on
GRCh38. NIST identifies v5.0q as current and deprecates v4.2.1 for HG002. The
truth VCF and benchmark BED come directly from NCBI's GIAB release; aligned
35× NovaSeq chr20 reads are the public DeepVariant case-study data. Every source
URL and SHA-256 is pinned in [`resources.tsv`](resources.tsv).

No large genomic artifact is committed. Preparation downloads roughly 2 GiB,
verifies every source before use, requires `samtools`, performs indexed chr20
extraction with a version-neutral header, and writes a local SHA-256 data manifest:

```sh
benchmarks/giab/prepare.sh
cargo xtask giab image plan --output image-plan.json --json
cargo xtask giab image dispatch --plan image-plan.json --confirm PLAN_ID
# Merge the generated digest-lock PR, then:
benchmarks/giab/run.sh
```

The report includes both Rosalind's internal regression evaluator and an external
Illumina hap.py v0.3.15 evaluation using RTG vcfeval 3.12.1. Both run all-emitted
and PASS-only calls inside the high-confidence BED; hap.py additionally reports
SNV-specific metrics across pinned GIAB v3.1 low-complexity, low-mappability,
segmental-duplication, and all-difficult contexts. The container runs with no
network and must be selected by immutable digest. Source archives, RTG, and the
linux/amd64 base image are checksum-pinned in [`happy/lock.json`](happy/lock.json).

The combined report also carries genotype concordance and call counts, memory
prediction/realization, receipt claim, tool versions, and exact argv. Routine runs
never edit the committed baseline. The first run emits a pending candidate; later
divergence emits an honest candidate and marks the workflow failed only after its
evidence is attested and uploaded. `cargo xtask giab baseline propose` requires a
reason-specific `GIAB baseline update` changelog entry and opens a pull request.

This is a credibility baseline, not a competitive threshold. Rosalind remains
SNV-focused, and the v5.0q truth includes difficult regions and variant classes
the caller does not yet model. Pull requests run only [`smoke.sh`](smoke.sh);
the full benchmark is manual/scheduled.

Primary guidance: [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle).
