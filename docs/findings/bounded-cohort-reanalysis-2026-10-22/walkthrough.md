# Two-minute saved-evidence walkthrough

**Evidence status:** the [dated report](README.md) records the passing hosted
execution, exact source/binary and observed performance losses. Prepare or open
the retained artifacts before presenting; the walkthrough duration is narration
time, not a promise about setup time.

## Prepare on the validated host

Run from the repository root with `ROSALIND_BINARY` set to the explicit hosted-built
executable for the revision identified in the final report. The output directories
must be new. No command below builds or installs Rosalind.

```sh
mkdir -p validation-results
python examples/cohort-reanalysis/prepare.py validation-results/fixture
python examples/cohort-reanalysis/validate.py validation-results/fixture \
  --binary "$ROSALIND_BINARY"
python research/cohort-reanalysis/demo.py \
  --prepared validation-results/fixture \
  --binary "$ROSALIND_BINARY" \
  --output validation-results/retained-demo
```

The focused hosted workflow already runs this sequence and retains its outputs.
Use those artifacts directly when available. Confirm
`validation-results/retained-demo/report.json` says `status: passed`, then check its
oracle hash, native binary hash and recorded command exit codes. The demo requires
the same binary as the unchanged independent fixture validator.

## 0:00–0:35: Pose a precise question

“I have saved evidence from three named specimens. For this A-to-C candidate at
position 40, how many specimens were measured, how many reach callable depth 10,
and how many of those have ALT support?”

Show `validation-results/retained-demo/members.json` and
`validation-results/retained-demo/plan/plan.json`. The catalog names exactly one
portable leaf per specimen. The plan states `metadata_only`, records the shared
comparison contract, and makes missing coverage explicit. It has not checked
candidate REF against stored row payloads yet.

## 0:35–1:10: Show missing versus measured zero

The unchanged authored fixture specifies:

| Specimen at position 40 | Stored observation | Callable depth | C count |
|---|---|---|---|
| specimen-a | Yes | 12 | 2 |
| specimen-b | No | Unmeasured | Unmeasured |
| specimen-c | Yes | 0 | 0 |

Show the corresponding row in
`validation-results/retained-demo/run/summary.jsonl`. A passing run must report
`n_requested=3`, `n_observed=2`, `n_depth_eligible=1`, `n_alt_supported=1`,
`total_callable_depth=12`, and `total_alt_count=2`. Its ALT fraction is the exact
integer pair `2/12`. The cohort summary has no floating-point conversion.

“Specimen C was measured and had zero callable observations. Specimen B's position
was not saved. Treating both as zero would destroy the observation denominator.
ALT support is a read-count screen, not a genotype or clinical conclusion.”

The fixture's raw inputs happen to contain reads at B:40. Partial saved-data
reanalysis still reports B:40 as unmeasured; silently reopening those reads would
answer a different question. The independent validator is unchanged at
[`examples/cohort-reanalysis/validate.py`](../../../examples/cohort-reanalysis/validate.py),
with authored values in
[`expected.json`](../../../examples/cohort-reanalysis/expected.json).

## 1:10–1:50: Show completed verification and replay

```sh
python research/cohort-reanalysis/analyze.py verify \
  validation-results/retained-demo/run
```

Open `validation-results/retained-demo/replay.replay.json` and the demo's
`report.json`. The retained demonstration temporarily makes only the generated
fixture's original `raw` paths unavailable, verifies every portable leaf through
the native reader, and repeats the query with window size 1 instead of 2. A passing
certificate states `reproduced`; the summary bytes and scientific identity agree.
The script restores the generated raw directory in `finally`.

“The wrapper verifies the retained output bytes and replays the actual scientific
calculation. Native children verify the saved data they consume. Their receipts
do not independently certify the outer Python sum. The result bundle records the
tested executable, source, comparison contract and output identities.”

## 1:50–2:30: Explain the engineering decision

“A narrower window reduces retained counters, but may reread and hash the same
physical Arrow partition and launch more native processes. Projection also does
not eliminate the physical decoder's cost. We compare total setup-plus-query cost
and sampled process-tree memory against existing saved queries and a bounded
Parquet baseline before choosing any optimization.”

Show the final report's measured table only after reviewing its passing run and
data scope. Retain observed losses. Do not present the diagnostic synthetic pilot
as a universal speedup, a hard memory cap, or fresh-versus-saved amortization.

If useful, show `validation-results/archived-hg002/results.json`. The hosted
`research/cohort-reanalysis/archived_real.py` check selects a hash-pinned historical
HG002 saved leaf with 1,000 loci across 62 physical partitions and compares current
summaries at thresholds 5, 10 and 20 against retained independent oracle outputs,
using windows 128 and 512. Those are declared input/configuration counts; the
report's passing status establishes whether the current checks succeeded.
This is one real specimen and a saved-only correctness check. The raw extraction
and independent oracle were not rerun, and its single-run timings are not a
performance study.

The earlier 108-run HG002 engine study remains historical evidence tied to its
own source and binary. It is not 108 new runs of this research extension.
