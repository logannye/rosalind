# Evidence sufficiency: HG002 P0

This experiment asks whether compact saved read evidence improves candidate-SNV
triage enough to justify retaining it. The practical target is fewer false
positive candidates at a fixed recall, with reusable evidence and explicit costs.
The current [protocol](protocol-p0.json) tests feasibility within HG002. It cannot
establish independent-person performance or clinical validity.

The full [implementation plan](PLAN.md) defines the subsequent locked study.
Previous cohort performance stops remain in force; this experiment does not
establish a speed advantage for Rosalind over Parquet.

## Execute on hosted compute

Use `.github/workflows/evidence-sufficiency.yml`. Native compilation, native
tests, acquisition, extraction and model fitting run on hosted Ubuntu with
Python 3.11 and the pinned binary wheels in `requirements.lock`. Pushes to the
implementation branch run P0 when its code/protocol changes. A manual dispatch
can run validation only or repeat the exact frozen pilot in a new run directory.
No P1 job is implicitly launched by a successful Actions job.

The initial hosted run `39878127360` completed correctly but had no assessed
calibration/threshold labels in the 30–31 Mb window. No models were fitted.
The workflow now defaults to the sole permitted count-driven expansion,
[`protocol-p0-expanded.json`](protocol-p0-expanded.json), frozen before its run.
It covers 15.99 Mb; the former descriptive test becomes development, and
55–59 Mb supplies a fresh test. Candidate rules, models and numerical gates are
unchanged. Manual dispatch can still select the original protocol.

The workflow runs the workspace native tests, existing Python integration tests,
independent evidence oracle suite and new scientific boundary/estimator tests.
Only then does it invoke:

```bash
python -u research/evidence-sufficiency/run_pilot.py \
  --binary target/release/rosalind \
  --output validation-results/hg002-p0
```

Native builds and scientific execution belong on hosted compute. Locally, the
standard-library tests can run without installing Rosalind or native dependencies:

```bash
python3 -m unittest discover -s research/evidence-sufficiency/tests -v
```

Estimator tests explicitly skip if sklearn/SciPy are absent. Hosted validation
requires those dependencies, so its estimator tests must run.

## Frozen experiment and limits

The 3.99 Mb of disjoint intervals use separate training, probability calibration,
threshold-selection and descriptive-test roles. Original-depth candidates are
generated before truth is loaded; both read-depth conditions use that same list.
The reduced BAM hashes the source, seed, read group and template name while
preserving physical records and mate decisions. Counts remain read counts.

Four fixed models cross counts/rich evidence with logistic regression/gradient
boosting. Preprocessing fits on training only. Sigmoid calibration and the 99%
recall threshold use separate roles. Simple depth/fraction thresholds and original
depth bcftools confidence provide baselines. Reduced-depth bcftools performance
is explicitly unavailable because that caller has not been rerun there.

The run admits at most 4 GiB downloaded, 30 GiB analysis scratch, 2 GiB retained
evidence, 100,000 candidate loci and 300,000 candidate ALTs per condition. Analysis
has a two-hour deadline, a sampled 14 GiB process-tree RSS stop within the 16 GiB
host allocation, and a workflow timeout. Sampled limits can miss brief peaks;
they are not a formal whole-process memory proof. Input acquisition and build
costs are separate from model timing. No GPU is used.

## Evidence and decision

The artifact retains exact source, executable, protocol, input and output
identities; raw downloaded BAM/reference/truth files are excluded from upload.
`run-manifest.json` binds the retained result files. `failure.json` means the
scientific run did not complete; an Actions success alone is not a P1 gate.

| Artifact | Purpose |
|---|---|
| `report.json`, `study/gate.json` | Combined correctness and development-signal decision |
| `identity.json`, `prepared.json`, `input-recheck.json` | Exact execution and input identities |
| `candidate-scan.json`, `candidates.jsonl` | Truth-blind candidate provenance and original depth |
| `label-domain.json`, `candidate-recall.json` | Assessment scope, exclusions and candidate-generation recall |
| `truth-recovery.json` | Accepted positives divided by all assessed truth SNVs, including missed candidates |
| `records.jsonl`, `study/predictions.jsonl` | Per-candidate features, labels, missingness and frozen decisions |
| `study/frozen-models.json`, `.pkl` | Model/threshold freeze before descriptive test evaluation; load only trusted bundles |
| `count-audit-*.json`, `replay.json`, `prediction-equivalence.json` | Independent counts, saved-only replay and native/Parquet equality |
| `telemetry.json`, `study/metrics.json` | End-to-end phase costs, sampled memory/storage and model metrics |

P1 requires all correctness checks, both classes in each fitting/calibration role
and training fold, at least 50 original-development negatives, and mean blocked-CV
AP improvement of at least 0.005 for richer evidence within at least one learner.
The descriptive test cannot select advancement. The runner accepts only the
recognized development gate with every required external check true.

If counts are inadequate, the plan permits one newly frozen count-driven
expansion up to 16 Mb. Adequate counts with no qualifying signal stop automatic
feature/model expansion. Failed correctness invalidates the run; infrastructure
repairs must preserve the scientific protocol. Passing P0 authorizes preparing
and running P1 only after its reference, inputs, splits and analysis are frozen.

Saved evidence and Parquet contain identical observations. This pilot verifies
their equivalence; its timings do not constitute a matched backend performance
benchmark. Model recall is conditional on assessed candidates; total truth-SNV
recovery is reported separately, and both retain the conservative GIAB label mask.
