# Verified CRAM decoder and reader-lifetime candidate

Candidate `9b8e12f3b6e802d2103c5ebeb5ff89a08c8d85ce` passed the fresh 108-run
HG002 representative matrix on 2026-09-18. Every completed scientific output
matches the original baseline byte for byte, and the post-run audit found no
prediction underestimates or changed final input, binary or harness identities.
The prior [failed attempt](../attempt-9862692/README.md) remains intact.

## Correction and regression

The earlier correction accounts for complete CRAM decoder state and records its
whole-file validation work. The rerun exposed a separate ownership error in the
pinned rust-htslib 0.44.1 indexed reader: file closure preceded destruction of a
CRAM index whose destructor still needed that file. The retained crash points to
`cram_index_free_recurse`, consistent with the
[upstream lifetime fix](https://github.com/rust-bio/rust-htslib/pull/518).

Rosalind's private indexed-reader owner now destroys the active iterator and index
before its ordinary reader closes the file. This keeps the pinned HTSlib ABI and
the Rust 1.83 minimum. Exclusive ownership preserves the existing `Send` contract;
the wrapper is not shared between workers. The regression repeatedly opens,
reads and drops eight parallel CRAM workers, compares with BAM, exercises failed
initialization and checks worker-factory survival after the original engine drops.

The new regression was also run against `9862692` in a separate temporary
worktree and reproduced the SIGSEGV ([negative-control log](lifetime-before.log)).
The original user checkout was not modified.

## Exact candidate evidence

| Check | Observed result | Evidence |
| --- | --- | --- |
| Representative workloads | 108 requested and completed; all schedule, equality, verification and budget gates pass | [Report](report.json), [raw archive](matrix-raw.tar.xz) |
| Historical output comparison | No differing output hashes | [Comparison](baseline-comparison.json) |
| Independent post-run audit | Pass; zero issues and zero prediction underestimates | [Audit](audit.json) |
| Linux focused tests | 31 tests pass, including the reader-lifetime regression | [Tests](linux-tests.log), [platform](linux-platform.txt) |
| Real cgroup cases | Seven scenarios pass for each of BAM and CRAM | [BAM](cgroup-bam.json), [CRAM](cgroup-cram.json) |
| Installed wheel Python API | 23 tests pass on each of Python 3.9 and 3.11 outside the checkout | [3.9 log](python39.log), [3.11 log](python311.log) |
| Packaged onboarding | Both analyzer APIs scaffold and pass conformance; standalone evidence analyzer and installed-wheel README pass | [Log](onboarding.log) |

[Candidate identities](candidate-identities.json) records the exact commit, clean
source archive, wheel, installed macOS executable, Linux executable and workload
manifest hashes. The macOS wheel was built from a clean temporary worktree at
that commit. Linux was built from a separate clean worktree at the same commit.
Matrix receipts record this commit and `code_dirty: false`.

The Linux executable is a debug build running as Linux/amd64 through the existing
local `colima-rosalind-ci` VM on an aarch64 Mac. Its configuration remained four
CPUs, approximately 4 GiB RAM and the existing disk allocation. This is a real
Linux cgroup-v2 test on the same physical development workstation, not an
independent-machine replication. [Build](linux-build.log),
[dependencies](linux-dependencies.log), [BAM raw evidence](cgroup-bam-raw.tar.xz)
and [CRAM raw evidence](cgroup-cram-raw.tar.xz) preserve the executed environment.

Both encodings completed and verified at 128 MiB in cooperative and OS-limit
modes, with equal output bytes. Preflight refusal exited 3; startup-budget and
record-capacity failures exited 4. The separate allocation control and native
8 MiB run exited 137 with Docker OOM classification and increased kernel
`oom_kill` counters. Cgroup memory includes file-cache charges and must not be
described as process RSS or a universal memory bound.

## Inputs, retained artifacts and limitations

Public source downloads were re-created from the locked HG002 source list with
pysam 0.23.3. Reference, BAM and non-CRAM prepared operands match the historical
hashes. Locally re-encoded CRAM/CRAI hashes differ because the reference location
changed; preparation and final inputs have their own recorded identities.
[Preparation](preparation.json), [preparation script](prepare_representative.py),
[source lock](source-lock.tsv) and [workloads](workloads.json) preserve that
derivation. The dataset is a chromosome-wide sparse query plus a 1-Mb window,
not a whole-genome performance study or validated variant-calling benchmark.

The matrix archive retains all TSVs, Arrow cache partitions, receipts, command
arguments, verification outputs, timing logs and the frozen harness. It expands
to approximately 263 MiB. The cgroup archives retain the actual harnesses and
kernel measurements. [Source archive](source-9b8e12f.tar.gz) and
[retained identities](retained-identities.json) permit independent inspection.
Genomic inputs, executables and the wheel are identified by hashes rather than
embedded in this report. Reproduction uses the frozen source and the
[representative harness guide](../../../../benchmarks/evidence/README.md) and
[cgroup guide](../../../../benchmarks/evidence/cgroup-probe.md), with new output
directories; paths in retained commands describe the original run.

Timings overlap other local validation and are not an isolated performance
ranking. Filesystem caches were not evicted. A resumed application cache with no
indexed extraction can still perform CRAM validation. These results do not prove
all CRAM codecs, high-depth assays, independent adoption, 30-day return usage or
clinical validity. The local artifacts are not published release packages.

## Supplementary foundation and release checks

[Full workspace tests](workspace.log) and the [Rust 1.83 check](msrv.log) passed
in the shared checkout before the final `Send` compile assertion; they include
release-inventory changes and are supplementary to the clean-candidate runs.
Do not attribute them to an unchanged clean `9b8e12f` tree. The
[CLI inventory](cli-inventory.json) matches all 49 policy entrypoints, including
dataset commands and both analyzer API modes; that policy change is separately
committed as `48c63f3`.

Read-only [release prerequisites](release-prerequisites.json) and
[doctor with the explicit Colima context](release-doctor-colima.json) show that
the existing Linux daemon was healthy. Those historical snapshots reported
missing allowed refs in the protected `rc` and `release` environments, an absent
`CARGO_REGISTRY_TOKEN`, and local actionlint 1.7.7 rather than pinned 1.7.12.
The [default-context doctor](release-doctor.json) additionally reports the unused
`desktop-linux` daemon unavailable. Its failure does not describe the explicitly
selected Colima validation environment. PyPI/TestPyPI publisher configuration
was not verified in authenticated index settings. No credentials, approval
rules or deployment settings were changed.

The later [private-tool check](actionlint-1.7.12.log) passes workflow validation
with official actionlint 1.7.12. Its macOS/arm64 archive was verified against both
the official release checksum file and GitHub's asset digest before extraction
into a task-private `/tmp` directory. [Provenance](actionlint-1.7.12-provenance.json)
records the download, executable and workflow hashes. The global 1.7.7 binary
remained unchanged. A new [doctor snapshot](release-doctor-private-tools.json)
using that private PATH and the explicit Colima context passes every local tool
and daemon check; only the protected-environment ref policies and registry token
remain reported blockers. Historical doctor snapshots are unchanged.

After these checks, the existing Colima VM was restored to its earlier stopped
state, and the default Docker context remained `desktop-linux`
([cleanup record](local-environment-restoration.json)). A future Linux or release
check must start the VM and select its context again. The private tool directory
is temporary; future runs must install or select the pinned version explicitly.

These publication prerequisites are distinct from the later stable-release
adoption and soak gates. Neither a successful local check nor this report
authorizes bypassing protected release approval.
