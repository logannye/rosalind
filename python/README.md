# rosalind-bio for Python

For version selection and the current source workflow, start with
[installation](../docs/installation.md). The public stable release is 0.1.0;
an unqualified registry install is not evidence that this preview is installed.
The [documentation index](../docs/index.md) connects this interface to the
researcher and saved-dataset quickstarts.

The mixed Maturin wheel bundles the matching `rosalind` executable. The distribution
is `rosalind-bio`; import `rosalind`. Python 3.9+ is declared and tested by candidate
wheel CI. This development line is 0.6.0 with experimental cohort interfaces;
the 0.5 evidence release remains a separate release effort. These changes are
unpublished. A source version does not establish registry availability.

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
With the environment created above:

```sh
cd /tmp
/tmp/rosalind-python/bin/python -c 'import rosalind; print(rosalind.__version__); print(rosalind.__file__)'
/tmp/rosalind-python/bin/rosalind --version
/tmp/rosalind-python/bin/rosalind analyze evidence --help
```

The module path should be inside that environment's `site-packages`; compare the
printed Python and native versions before following a tutorial.
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
CRAM follows the native [supported decoder profile](../docs/SEMANTICS.md#cram-decoder-admission)
and performs one complete file-validation pass before yielding evidence, including
on native cache resume. This work is included in run setup and recorded separately
from indexed record visits. `open_dataset()` avoids the original alignment inputs.

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

For physical projection, pass `fields=["depths", "alleles"]` to `iter_evidence`
or `materialize_evidence`. Native accumulation, Arrow batches, and persisted
artifacts contain only those groups and the locus identity columns. Missing
columns are absent. `panel_qc` defaults to depths and quality sums; an explicit
field list must include both groups. Use `fields="all"` for complete evidence.

The next-minor cohort preview adds immutable local sample collections above
portable datasets. Create a snapshot with `rosalind cohort create --cohort cohort
--members members.tsv`; the table contains `id` and `manifest` columns. Each
member must have named single-sample evidence. Use the returned `snapshot_id`:

```python
from rosalind import open_cohort

cohort = open_cohort("cohort", snapshot="<snapshot_id from cohort create>")
plan = cohort.plan(sites="candidates.vcf")
result = cohort.materialize("candidate-review.tsv", sites="candidates.vcf", format="tsv")
summary = cohort.summarize("candidate-summary.arrow", sites="candidates.vcf",
                           min_callable_depth=10)
# Ask a second question using only the cohort directory and a new candidate file.
with cohort.batches(sites="second-candidates.vcf", missing="partial") as run:
    for batch in run:
        print(batch.num_rows)
    assert run.result is not None
```

`open_cohort()` is lazy. Each native call verifies the immutable snapshot and
consumed evidence. Strict coverage refuses unmeasured loci; explicit partial
coverage emits null evidence and status while retaining measured zero depth.
The depth threshold is a technical screen. ALT support proportions are not
genotypes or population allele frequencies. `members=None` selects all members;
`members=[]` selects none. Sample labels and metadata are user assertions.

`batches()` currently materializes a complete native Arrow artifact on disk
before yielding its first batch, then reads at most 1,024 rows per batch. A
context manager owns the child process and reader. Its `workdir` (a new temporary
directory by default) retains the artifact and receipt; remove it when no longer
needed. For an explicitly named artifact use `materialize()`. Native budgets
exclude Python-retained state; `memory_budget_mb` enables cooperative admission,
and `require_os_limit=True` additionally requires a verified Linux cgroup limit.

To fill missing loci, supply a TSV with `id`, `role`, and `path` columns containing
exactly the affected samples and all original scientific source roles. Paths are
relative to the table; their content identities must match the original sources.

```python
cohort.plan(operation="extend", sites="second-candidates.vcf", sources="sources.tsv")
extension = cohort.extend("sources.tsv", sites="second-candidates.vcf", workdir="scratch")
extended = extension.cohort
print(extension.report["members"])
extended.verify()
```

The scratch directory must exist outside the cohort. Extension computes missing
loci only, publishes a new snapshot after all additions succeed, and preserves
the original. It cannot add missing fields at already stored loci. A source-free
no-op uses a header-only source table and returns the same snapshot. Plan reports
do not open raw sources; execution verifies their hashes and reports native
record visits, CRAM full-validation work and existing-object hashing separately.
Large explicit member selections for extraction, summaries and their plans use
a bounded temporary native argument file, removed when the invocation closes.
Extension currently supports at most 64 KiB of command arguments; select all
members (the default) and use the source table for large extension operations.
