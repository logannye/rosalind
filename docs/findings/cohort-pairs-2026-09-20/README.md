# Paired research preview — September 20, 2026

The paired preview compares explicitly ordered saved samples for supplied SNVs.
It preserves each side's counts, missingness and technical depth screen, and
reports exact right-minus-left ALT/read-depth fractions. No pairing is inferred
from sample metadata. This is separate from the unpublished 0.5 release.

Tested implementation: `3e75ca0101daff70c744fe407297a333051e1b16`, in
[PR #148](https://github.com/logannye/rosalind/pull/148). It integrates the reviewed
worktree producer-identity fix from #147 and the latest #145 cohort checkpoint.
The local development binary identifies this exact commit with
`code_dirty=false`; its SHA-256, target and dependency-lock identity are retained
in the demonstration reports. The follow-up `9c485c4` adds the paired tutorial assets to offline bundles;
its staged cohort/paired example and 12 onboarding tests passed locally. Later
evidence/documentation commits do not change the tested native implementation.

## What a researcher can inspect

Open the [readable comparison report](REPORT.md) and [native rows](pairs.tsv).
Twelve comparisons cover A→B, B→A and A→C at four candidates. The example includes
measured zero, low depth, absent evidence and swapped direction. Core regressions
also cover equal fractions and unsigned-64-bit maximum counts.
Each side's integer counts, fraction components, support and eligibility are
checked against authored expectations; differences are independently calculated
with Python `Fraction`. The native receipt is verified after the comparisons.
The saved-only query records zero original alignment records decoded.

- [Paired demonstration and source identity](report.json)
- [Complete cohort workflow, extension, relocation and total costs](cohort-demo.json)
- [Native paired receipt](pairs.tsv.manifest.json)
- [Validation inventory and log hashes](validation.json)

The small synthetic demonstration is author-run engineering evidence. It does
not establish real-cohort economics, independent use, biological response,
molecule counts, rare-variant sensitivity or clinical validity.

## Validation and retained failures

The initial integrated suite found an obsolete replay-error substring assertion.
The native allowlist correctly refused the forbidden requests; the regression
now checks the stable refusal prefix. One subsequent local attempt failed during
reference-pack setup in `partition_arrow_merge`; its isolated rerun passed.
Its root cause was not established and the failed log remains retained. A later
sandboxed run reached Receipt Studio but could not start its loopback server;
the final local run with localhost socket access passed the full workspace suite. No failed run is recorded
as successful.

The Linux probe now accepts an explicit read-only pair table and adds successful
512 MiB OS-limited paired execution plus tiny-budget refusal. It compares output
with the outside-probe baseline, checks receipts and observed cgroup counters,
and mounts no original alignments or references. All five actual Linux cases passed; the [retained cgroup report](linux-cgroup.json)
and [raw CI artifact archive](linux-ci-evidence.tar.gz) preserve their evidence.
The paired run peaked at 4,853,760 cgroup bytes on this tiny fixture; that is not
a representative cohort memory measurement. GitHub built merge commit
`4cc38df76ff2f6996af4cd41630aa1a72f2b0630`, whose tree matches `3e75ca0` exactly
(`fa48c5946d3c01f081c0f09620ac9166df51c909`). The [commit metadata](linux-merge-commit.json)
and [Linux paired report](linux-pairs.json) retain these distinct identities.

Remote runs for the tested implementation:

- [Linux CI, cohort demonstration and cgroup cases](https://github.com/logannye/rosalind/actions/runs/35545019393)
- [Linux x86_64 and macOS arm64/x86_64 wheels](https://github.com/logannye/rosalind/actions/runs/35545019389)
- [Nextflow and Snakemake integrations](https://github.com/logannye/rosalind/actions/runs/35545019399)

Commands, completed outcomes and any outstanding checks are recorded in the
validation inventory. The local binary is a development build, not a release
bundle. There is no local Docker daemon, so Linux execution evidence comes from
GitHub's Linux runner. The fixture does not establish cohort-scale performance
or a universal memory guarantee.

## Release continuation

Release fixes #146 and #147 passed their checks but remain unmerged pending
explicit approval. Automatic approval review rejected the attempted main-branch
merges because the request to resume work did not explicitly authorize that
shared publication step. Local preview integration and draft-PR updates continued.

The next release task remains an immutable **0.5.0-rc.3** from corrected main,
followed by actual public-artifact verification and the existing seven-day soak.
The prior evaluator publication should be recovered via #149 without rebuilding
the already-published image. Stable registry configuration and the technical
benchmark gate remain separate tasks. Independent participants are advisory and
must not become a new shipping gate.
