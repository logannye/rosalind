# Python boundary for Rosalind

The dependency-light way to use Rosalind from Python is [`rosalind.py`](rosalind.py) — a
stdlib + numpy module that runs the `rosalind` binary and loads its **bit-reproducible**
per-locus feature table.

## Prerequisites
- Python 3.8+ with `numpy`
- A `rosalind` binary (`cargo build --release` → `target/release/rosalind`, or `install.sh`)

## Usage

```python
from rosalind import features

# Streams the whole-genome feature table to disk in bounded memory, then loads it.
ft = features("ref.idx", "sorted.bam", binary="target/release/rosalind")
print(ft.columns)            # contig, pos, ref, depth, A/C/G/T counts, strand, mean bq/mapq, …
print(ft.data.shape)         # one row per callable locus
```

Because `rosalind features` is byte-identical run-to-run with a BLAKE3 receipt, the loaded
table is a **bit-reproducible** ML training input: hash it, and you can prove this quarter's
model saw exactly the same data as last quarter's. The runnable proof is
[`../examples/reproducible_features_demo.py`](../examples/reproducible_features_demo.py)
(`python3 examples/reproducible_features_demo.py target/release/rosalind`): it extracts
features twice and shows the trained weights are bit-identical.

> A zero-copy in-process `pyarrow` binding is a future step; the subprocess boundary above is
> the dependency-light entry point today.
