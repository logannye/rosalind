# Saved HG002 depth calibration

This separate exploratory experiment diagnoses the reduced-depth threshold
failure found in [P0](../evidence-sufficiency/RESULTS.md). It compares depth-aware
counts with rich evidence while retaining the original `stop_no_signal` decision
and P1 lock. All data, including the saved `test` region, have previously been
inspected; this is not fresh validation.

The [protocol](PROTOCOL.md) and [executable configuration](protocol.json) freeze
the arms, fitting roles, matched-target-recall comparisons, paired genomic-window
uncertainty, costs and decision screen before execution. Original-depth and
10/35-template-retained rows remain paired within their existing roles. The
original candidate list and conservative GIAB label scope are unchanged.

## Hosted execution

`.github/workflows/depth-calibration.yml` runs on pushes of experiment source to
`codex/hg002-depth-calibration`. It installs the exact original binary wheels,
runs synthetic scientific and provenance tests, and retrieves the pinned artifact
from successful P0 run `39879349288`. It verifies provider run identity, ZIP SHA256,
every consumed member, the original inventory, model and dependency identities
before deserializing models or fitting new ones.

No native source is changed or compiled and no raw reads or new truth data are
downloaded. All scientific fitting, prediction, evaluation and timing run on
hosted Ubuntu with one numerical thread, a 15-minute analysis limit, sampled
8 GiB RSS limit and 1 GiB retained-output ceiling. Infrastructure and wheel
installation are outside the scientific timer; the workflow has a 25-minute cap.

The result artifact contains the exact protocol/source/input identities, a model
freeze written before new test scoring, predictions, metrics, paired uncertainty,
cost repetitions and an output SHA256 inventory. Scientific and infrastructure
failures are distinct; a successful Actions job does not authorize P1 or selective
escalation. Source raw alignments are unavailable to this workflow.

## Cost and decision boundary

The saved artifact already contains rich evidence. This experiment measures
saved feature size and repeated parsing/matrix/inference costs. It cannot measure
the cost of initially extracting richer fields or falling back to reads/callers.
Even a positive exploratory predictive screen can only nominate fresh-data and
complete-cost validation. Selective escalation remains unauthorized until its
practical advantage and full workflow cost are established.

Local checks are limited to static inspection and provenance unit tests:

```sh
python3 -m unittest discover -s research/depth-calibration/tests -p test_provenance.py -v
git diff --check
```

Run the complete test suite and `run.py` only on hosted compute. The runner
explicitly requires the hosted Linux Actions environment.
