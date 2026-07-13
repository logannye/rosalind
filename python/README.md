# rosalind-bio for Python

The wheel bundles the matching `rosalind` executable and exposes a typed Python
package imported as `rosalind`.

```python
from rosalind import iter_features, inspect_receipt

run = iter_features("GRCh38.rref", "sample.bam")
for batch in run.batches():
    train(batch)

receipt = inspect_receipt(run.result.manifest_path)
```

`iter_features` consumes Rosalind's native canonical Arrow IPC stream and keeps
Python memory bounded by one record batch. `collect_features` is the explicit
whole-result materialization path. Genomic data remains local.
