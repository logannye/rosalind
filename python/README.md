# rosalind-bio for Python

The mixed Maturin wheel bundles the matching `rosalind` executable. The distribution
is `rosalind-bio`; import `rosalind`. Python 3.9+ is declared and tested by candidate
wheel CI. The source version is 0.5.0; these changes are unpublished until the
release gates complete. A source version does not establish registry availability.

Installing a prebuilt wheel needs Python, pip, and its declared PyArrow dependency;
it does not need a compiler. Building a wheel from source needs Rust/Cargo and the
[native build prerequisites](../docs/analyzer-sdk.md#build-prerequisites), including
CMake and the C/C++ toolchain. The commands below also fetch build dependencies;
the source build is not an offline installation procedure.

Build and install from the repository root:

```sh
python3 -m venv /tmp/rosalind-python
/tmp/rosalind-python/bin/pip install maturin==1.9.4
/tmp/rosalind-python/bin/maturin build --locked --release --out /tmp/rosalind-wheels
/tmp/rosalind-python/bin/pip install /tmp/rosalind-wheels/*.whl
```

Run Python outside the repository to verify the installed package is being used.
RC versions are normalized between native `0.5.0-rc.N` and Python `0.5.0rcN`;
different RC numbers and stable versions do not count as a match.

Supply an indexed `genome.fa`, indexed `sample.bam`, SNV-only `candidates.vcf`,
and `targets.bed`. Start in a fresh output directory so create-new artifacts do
not collide with a previous run. This example counts rows without retaining
whole-genome Python state:

<!-- smoke:python-evidence -->
```python
from rosalind import iter_evidence, materialize_evidence, panel_qc

rows = 0
with iter_evidence("genome.fa", "sample.bam", sites="candidates.vcf",
                   memory_budget_mb=256) as run:
    for batch in run:
        rows += batch.num_rows
    assert run.result is not None
print(f"Extracted {rows} candidate loci")

result = materialize_evidence("genome.fa", "sample.bam", "evidence.arrow",
                              regions="targets.bed", memory_budget_mb=256)
summary = panel_qc("sample.bam", "targets.bed", "panel.tsv",
                   reference="genome.fa", min_callable_depth=10)
assert summary.returncode == 0
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

Maintainers run `scripts/smoke-wheel.sh path/to/candidate.whl --candidate-source
/absolute/path/to/rosalind` to install into a fresh environment, check exact
native/Python version identity, and execute this README example outside the
checkout. The generated analyzer is explicitly patched to that candidate source.
After the matching SDK is published, `--registry-sdk` replaces the source option
and requires a fresh registry-only analyzer build; a source-patched smoke is not
evidence of registry installability. The smoke fetches dependencies before its
locked offline SDK build and then performs local receipt verification/replay.

Legacy `iter_features` and `collect_features` remain available with their original
feature schema. `collect_features` explicitly materializes the whole result in
Python memory; its semantics differ from the new exact evidence profile.
