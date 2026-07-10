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
benchmarks/giab/run.sh
```

The report includes both all-emitted and PASS-only precision, recall, F1,
genotype concordance and call counts, plus memory prediction/realization, the
receipt claim, and exact replay argv. The first successful run replaces the
explicit pending baseline. Later metric changes fail unless invoked with
`--update-baseline` and accompanied by a changelog entry containing
`GIAB baseline update`.

This is a credibility baseline, not a competitive threshold. Rosalind remains
SNV-focused, and the v5.0q truth includes difficult regions and variant classes
the caller does not yet model. Pull requests run only [`smoke.sh`](smoke.sh);
the full benchmark is manual/scheduled.

Primary guidance: [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle).
