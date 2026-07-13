# Platform comparison

This harness compares Rosalind, a maintained pysam implementation, and the
explicitly equivalent subset of bcftools mpileup/query on the same reference and
BAM. It retains raw GNU-time measurements, three output hashes, normalized TSVs,
package locks, receipts, verification, Arrow output, and the sharded merge.

Run the bundled smoke dataset without network access inside the benchmark:

```sh
benchmarks/platform/run.sh
```

Run the real prepared HG002 chr20 data with:

```sh
benchmarks/platform/run.sh \
  benchmarks/giab/data/prepared/GRCh38.chr20.fa \
  benchmarks/giab/data/prepared/HG002.chr20.bam \
  benchmarks/platform/results-hg002
```

The environment image is built from immutable base-image digests, exact Python
requirements, and the current Cargo lock. Its report records resolved Debian and
Python packages. The narrative deliberately separates resource prediction and
refusal from runtime and scientific semantics. bcftools does not expose all of
Rosalind's canonical feature columns through the selected interface; those fields
are recorded as unmatched rather than silently treated as equivalent.
