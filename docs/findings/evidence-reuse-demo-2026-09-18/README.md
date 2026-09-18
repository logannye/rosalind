# Evidence-reuse demonstration — source preview, 2026-09-18

The [demonstration driver](../../evidence-reuse-demo.md) was executed against the
clean macOS arm64 binary built from
[`9b8e12f3b6e802d2103c5ebeb5ff89a08c8d85ce`](https://github.com/logannye/rosalind/commit/9b8e12f3b6e802d2103c5ebeb5ff89a08c8d85ce),
installed in an isolated Python 3.11 environment. The report records the executable
SHA256 and the source identity verified by replay. This is source-preview evidence;
0.5 publication, independent user acceptance, and the release-bound R07 recording
remain separate gates.

The corrected run passed:

- Four supplied candidates with callable depths 36, 34, 44, 27 and ALT observations
  17, 16, 22, 12, matching the pinned tutorial's expected evidence.
- 3,159 positions verified in the saved dataset.
- Relocation of the complete portable dataset and removal of generated BAM/FASTA
  inputs before saved-only queries.
- Exact TSV byte agreement with fresh extraction for a two-candidate shortlist and
  panel QC; both receipts verified.
- Candidate-output replay with the explicit binary, matching code identity, and
  verdict `REPRODUCED`.

The public NA18507 slices and the independent candidate-generation steps are
identified in the report's preparation record. Depth 10 is an illustrative
technical screen. This small-data run is neither a whole-genome performance study
nor a caller-accuracy comparison. Local execution paths in these copies are
replaced by descriptive variables; original artifact receipts remain in the local
run and are not published here. Source hashes and result-file hashes are retained.

| Run | Evidence | Result |
|---|---|---|
| Initial script development | [Transcript](initial-format-error/transcript.txt) · [Report](initial-format-error/demo-report.json) | Failed the byte check: the script omitted `--format tsv` on `dataset extract`, whose default is Arrow despite a `.tsv` filename |
| Corrected source-preview script | [Transcript](source-preview/transcript.txt) · [Report](source-preview/demo-report.json) | Passed with explicit TSV format and matching-code replay |

The initial failure is a demonstration-script error, not an evidence-counting
failure. Its report remains `passed: false`. The corrected script explicitly
selects TSV and checks byte equality, fixture values, replay verdict, and code
identity before reporting success. The complete corrected demonstration was run
twice; both produced identical candidate and panel output hashes.

Use the [rerun instructions](../../evidence-reuse-demo.md) with an exact published
candidate/stable executable to create its eventual release-bound recording.

## README checks

The recurring onboarding runner now extracts and executes the README's candidate,
saved-evidence, and Python reuse snippets. It reuses the researcher tutorial's
prepared fixture and the existing SDK steps; it does not build another analyzer.
Native bundles run the shell snippets; matching installed-wheel runs additionally
exercise Python after deleting only the copied original input files. Native-only
bundles explicitly report that the Python package check belongs to wheel smoke.

The updated staged workflows passed outside the checkout on this macOS arm64
machine. The demonstration also passed from the staged bundle, with output hashes
matching the recorded run. This targeted validation repeats the affected scientific workflows, not
the unchanged full SDK build or a registry installation. The separate source-only
README analyzer commands passed four fixture tests and 19 conformance checks when
the README was introduced in commit `92ac1cc`.

- [Native workflow log](onboarding-native.log)
- [Installed-wheel workflow log](onboarding-wheel.log)
- [44 passing script tests](script-tests.log), including wrong-version refusal and
  preservation of an existing demonstration output directory
- [Validation summary](onboarding-validation.json)
