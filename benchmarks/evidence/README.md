# Exact-evidence budget curve

Prepare the content-locked [NA18507 tutorial](../../examples/research-filter/README.md).
Then run with the same Python environment containing pysam 0.23.3:

```sh
python benchmarks/evidence/run.py --binary target/debug/rosalind \
  --data /tmp/research-filter --output /tmp/evidence-curve --cache
```

Default: three repeats at 96/128/256MiB, three independent pysam extractions, and
optional cold/resumed cache runs. Results directories are create-new. Any differing
row, repeat, failed verification, or nonzero extraction exit makes the run fail;
the report and raw commands/timing/stderr remain available for diagnosis.

The baseline directly traverses CIGAR and applies the documented flag, MAPQ, BQ,
missing-quality, read-cycle, allele, strand, and histogram semantics. It streams
one locus at a time and retains only supplied candidate metadata. It does not use
the Rosalind implementation or sampled pileup defaults.

`report.json` retains input/binary hashes, source preparation provenance, versions,
exact argv, process time/RSS, raw timing I/O counters, artifact sizes/hashes, every
semantic comparison, receipt setup/hash/combined-analysis-encoding measurements,
and independently timed verification. The baseline does not implement provenance;
its extraction cost is labeled separately from Rosalind's end-to-end contract.

No cgroup is imposed by this portable harness. A declared budget is not a measured
hard limit for the baseline. The tiny four-locus example establishes repeatability
and agreement only; scale the supplied candidate set and aligned data, retaining
the same evidence, before making performance claims. Kernel-only encoding timing,
large-data curves, and independent-machine reproduction remain separate gates.
