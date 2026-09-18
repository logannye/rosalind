# Release bootstrap correction, 2026-09-18

This is engineering evidence, not a published release or independent adoption.

## Retained failure

[RC2 run 35378970879](https://github.com/logannye/rosalind/actions/runs/35378970879)
ran workflow source `15934253caba6a26f9e00824a9515cd8fc46d283`. Its `gates`
job failed during checkout; all build and publication jobs were skipped. At
2026-09-18T18:14:13Z, checkout reported:

```text
A branch or tag with the name 'HEAD' could not be found
```

The local plan had resolved a commit but retained the literal input `HEAD` in
`metadata.ref`. Dispatch passed that local alias to the remote checkout action.
This failure occurred before building or publishing RC assets. It supplied no
release soak, installation, or scientific-validation result.

## Corrected behavior

- Release preflight records the full resolved commit as `metadata.ref`. RC and
  stable snapshots consume that frozen commit, not a second resolution of the
  original alias. Evaluator-image planning uses the same rule.
- Dispatch passes `report.commit` for these three workflow families. A legacy
  plan whose recorded `ref` still differs from its commit refuses before any
  GitHub call and asks the maintainer to regenerate the plan. Plans made from
  `HEAD`, its branch name and its SHA now have identical authenticated intent.
- RC and stable annotated-tag steps declare `GIT_COMMITTER_NAME` and
  `GIT_COMMITTER_EMAIL` in the step environment. No global Git configuration,
  required reviewers, workflow permissions or scientific release gate changed.

## Focused acceptance evidence

- A real local Git fixture exercises release preflight with `HEAD`, `main` and
  the full SHA, checks equal plan IDs, then moves the branch. The dispatcher still
  submits the original commit. Legacy mutable-ref plans for all three workflows
  are rejected before remote operations.
- The actual checked-in RC and stable publication shell blocks run against local
  bare repositories with ambient Git configuration disabled and
  `user.useConfigOnly=true`. They create and push annotated tags at the expected
  commit with the declared bot identity. A local `gh` stub records the release
  request; these tests perform no external publication.
- The full xtask test suite passed: 38 unit tests and 5 publication integration
  tests. The release-helper suite passed 13 tests. Clippy with warnings denied,
  Rust 1.83 all-target checking and Actionlint validate the changed source and
  workflows.

A successful replacement RC run must be linked separately. These local checks
do not claim that TestPyPI accepted an upload, GHCR is publicly accessible or the
prerelease exists. Its actual GitHub publication timestamp starts the soak.
