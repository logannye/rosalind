# Reproduce the evidence-reuse demonstration

This is a runnable **source-preview demonstration**, prepared for R07. It does not
establish publication of a release, independent adoption, caller accuracy, or a
speedup. The [recorded source-preview run](findings/evidence-reuse-demo-2026-09-18/README.md)
identifies its exact binary, source identity, pinned inputs, commands, and results.
A release-bound demonstration must be rerun after its candidate or stable artifacts
are published; the label alone is not evidence of publication.

## Run it

Use a supported macOS/Linux environment and an evidence-capable Rosalind binary.
For a source build, follow [installation](installation.md) and build the CLI.
The demonstration uses only the [public NA18507 tutorial inputs](../examples/research-filter/README.md).
Preparation downloads about 118 KB and verifies the pinned source hashes; the
subsequent saved-evidence queries need no original alignment/reference files.

Prepare a Python 3.9+ environment with the exact preparation dependency:

```sh
python3 -m venv /tmp/rosalind-demo-env
/tmp/rosalind-demo-env/bin/pip install pysam==0.23.3
```

From the repository root, select the binary and its expected version explicitly.
The output directory must not exist; everything generated is retained there:

```sh
PYTHON=/tmp/rosalind-demo-env/bin/python scripts/reproduce_demo.sh \
  --binary "$PWD/target/debug/rosalind" --expected-version 0.5.0 \
  --label source-preview --output /tmp/rosalind-evidence-demo
```

The [shell entry point](../scripts/reproduce_demo.sh) invokes the
[Python driver](../scripts/reproduce_demo.py). There is no fallback to a different
installed binary. For a published candidate, pass its actual executable and exact
native version, such as `0.5.0-rc.N`, and use `--label release-candidate`. For stable
publication, use the actual stable version and `--label stable`. These labels are
operator assertions; the report independently records the executable SHA256,
reported version, and code identity returned by replay.

## What the demonstration checks

1. Prepare the pinned reference, alignments, indexes, targets, and four supplied
   candidate SNVs. Extract evidence and verify its receipt.
2. Save depth, allele, and quality summaries at all 3,159 target positions.
3. Select a two-candidate shortlist and retain fresh candidate/panel answers for
   comparison. The shortlist selects supplied calls; it does not generate new calls.
4. Move the complete portable dataset, retain the VCF/BED selections, and delete
   **only this run's generated input directory**. Original user files are never inputs.
5. Query the relocated evidence, compare both result files byte-for-byte with the
   fresh answers, verify their receipts, and replay the candidate output with the
   explicitly selected executable.

The four candidate depths are 36, 34, 44, and 27; the corresponding supplied-ALT
counts are 17, 16, 22, and 12. Panel denominators are 1,575 and 1,584 positions;
1,493 and 1,517 respectively meet the illustrative depth-10 screen. These expected
values check the pinned fixture, not the correctness of arbitrary biological
interpretation. Counts are read observations; overlapping mates are counted
separately. Stored zero depth remains different from an unmeasured locus.

## Inspect and record it

The output contains results, receipts, the relocated portable dataset,
`transcript.txt`, and `demo-report.json`. The transcript records executed commands
and native output. It replaces machine-local paths with `$ROSALIND`, `$PYTHON`,
`$REPOSITORY`, and `$DEMO`; the report records hashes and observed pass/fail results.
Original receipt files retain their local execution paths and should be reviewed
before public sharing. Raw receipts are not included in the checked-in recording.

To create a terminal recording, run the same command through your terminal
recorder. The generated text transcript is available without an extra recording
dependency. Before linking a recording from a release, rerun against its published
binary, verify `passed: true`, and link the release plus the transcript and report.
Retain a failed report separately when correcting a failed demonstration.

No comparison with GATK, DeepVariant, or any other caller is established here.
This replaces the old caller-centered demonstration and its unsupported blanket
claims about other tools' determinism.
