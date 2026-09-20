# Paired comparison preview checkpoint — September 20, 2026

Progress is committed and pushed to [PR #148](https://github.com/logannye/rosalind/pull/148),
stacked above the isolated 0.6 cohort preview. The tested native implementation is
`3e75ca0101daff70c744fe407297a333051e1b16`. The follow-up `9c485c4` includes the
paired script and input table in offline tutorial bundles; it does not change
native or Python runtime code.

## Implemented and exercised

- Explicit ordered pairs; exact counts and right-minus-left fractions;
  missing/zero/low-depth distinctions; bounded canonical Arrow/TSV output;
  CLI/Python planning and reports; guarded inputs and relocated saved-only replay.
- Reviewed build identity correction #147 integrated locally. The clean local
  binary records the correct source SHA and `code_dirty=false`.
- Readable paired reports with exactly simplified fractions and explicit
  coverage labels. All 12 synthetic comparisons match authored counts and
  independent `Fraction` arithmetic, including side fraction/support fields.
- Complete local Rust workspace tests and all 31 Python tests passed. The
  obsolete replay-error assertion was corrected; other failed local attempts
  remain separately retained rather than relabeled as passes.
- Linux CI's paired demo and all five actual cgroup cases passed, including
  paired completion under a 512 MiB OS cap and refusal with a 1 MiB budget.
  Linux CI tests GitHub merge commit `4cc38df`; its tree is byte-identical to
  the tested head. This small fixture is not a scalability claim.
- Six cgroup harness tests, 12 onboarding tests, seven roadmap tests, actionlint
  and formatting passed. The complete staged cohort/paired example passed from
  its copied tutorial assets using the clean tested native binary.
- The CI-built Apple Silicon wheel passed all eight cohort Python tests after
  fresh installation outside the checkout with `PYTHONPATH` unset.

The [retained evidence](https://github.com/logannye/rosalind/blob/codex/cohort-pairs-preview/docs/findings/cohort-pairs-2026-09-20/README.md)
contains exact source identities, reports, failures, log hashes and remote run
links. Use its validation inventory for the final platform/wheel outcomes.
Subsequent evidence-only commits use `[skip ci]` to avoid rebuilding identical
native code; they are not represented as independently CI-tested source.

## Next decisions

1. Inspect the existing [wheel run](https://github.com/logannye/rosalind/actions/runs/35545019389):
   Linux and Apple Silicon passed both Python versions; Intel Mac installation
   checks remain running at this checkpoint. All ten CI jobs and both workflow
   integrations passed. Record the remaining result without launching duplicate
   builds. Then review the paired contract and retained platform results and integrate the
   0.6 previews in their intended release sequence. Keep them outside 0.5.
2. Merge the reviewed 0.5 fixes #146/#147 after explicit approval. Automatic
   approval review blocked those main-branch merges; both remain open. Their CI
   passed. #147's local integration here does not mean it is merged on main.
3. Follow the release continuation in `docs/EXECUTION_CHECKPOINT.md`: recover the
   already-published evaluator with #149, prepare immutable 0.5.0-rc.3 from
   corrected main, verify actual public artifacts, complete the technical
   benchmark/account prerequisites and observe the seven-day soak.
4. Measure recurring research use and total costs when users are available.
   Independent participation remains advisory. Do not claim biological response,
   molecular evidence, clinical validity or real-cohort economic savings from
   the authored demonstration.

The [paired contract](cohort-pairs-contract.md) defines interpretation and limits.
No source discovery, inferred pairing, production caller or molecule model is
introduced by this work.
