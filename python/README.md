# rosalind-bio for Python

The mixed Maturin wheel bundles the matching `rosalind` executable. The distribution
is `rosalind-bio`; import `rosalind`. Python 3.9+ is declared and tested by candidate
wheel CI. The source version remains 0.4.0; new evidence APIs are under development,
and no future-version public wheel is implied.

Build and install from the repository root:

```sh
python3 -m venv /tmp/rosalind-python
/tmp/rosalind-python/bin/pip install maturin==1.9.4
/tmp/rosalind-python/bin/maturin build --locked --release --out /tmp/rosalind-wheels
/tmp/rosalind-python/bin/pip install /tmp/rosalind-wheels/*.whl
```

Run Python outside the repository to verify the installed package is being used.
RC versions are normalized between native `0.4.0-rc.N` and Python `0.4.0rcN`;
different RC numbers and stable versions do not count as a match.

```python
from rosalind import iter_evidence, materialize_evidence, panel_qc

with iter_evidence("genome.fa", "sample.bam", sites="candidates.vcf",
                   memory_budget_mb=256) as run:
    for batch in run:
        consume(batch)
    assert run.result is not None

result = materialize_evidence("genome.fa", "sample.bam", "evidence.arrow",
                              regions="targets.bed", memory_budget_mb=256)
summary = panel_qc("sample.bam", "targets.bed", "panel.tsv",
                   reference="genome.fa", min_callable_depth=10)
```

`iter_evidence` yields canonical native Arrow batches of at most 1,024 rows and
requires exactly one `sites` or `regions`. The native budget excludes memory
retained by Python consumers. Do not accumulate batches when bounded Python
memory matters. Stream exhaustion sets `run.result`; a context manager cancels
the child when iteration ends early. A stream has no persisted output artifact.
`materialize_evidence` writes an artifact and receipt for byte verification/replay.

The versioned profile defaults to MAPQ20/BQ20, counts overlapping mates separately,
excludes unavailable qualities, and never downsamples. See
[SEMANTICS.md](../docs/SEMANTICS.md). Errors expose native exit codes through
`EvidenceProcessError`. Cache/resume and workers are explicit execution options;
all genomic processing stays local.

Legacy `iter_features` and `collect_features` remain available with their original
feature schema. `collect_features` explicitly materializes the whole result in
Python memory; its semantics differ from the new exact evidence profile.
