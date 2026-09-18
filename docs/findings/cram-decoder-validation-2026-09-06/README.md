# CRAM decoder correction and native reader lifetime

The first decoder-envelope correction exposed a native index-lifetime failure.
The corrected candidate now completes the same representative workload and Linux
failure probes. Both attempts are retained; the passing rerun does not rewrite
the historical failure.

| Candidate | Representative matrix | Status |
| --- | --- | --- |
| [`9862692`](attempt-9862692/README.md) | 106 of 108 completed; parallel CRAM SIGSEGV and blocked dependent query | Failed |
| [`9b8e12f`](verified-9b8e12f-2026-09-18/README.md) | 108 of 108 completed; all output hashes match the original baseline | Passed within the tested scope |

The [retention manifest](retained-attempt-provenance.json) records the original
failed attempt's source location and per-file hashes. It was copied byte for byte
from the existing checkout without modifying those source files.

The verified candidate also passed 14 real Linux cgroup scenarios, installed
Python tests on 3.9 and 3.11, and packaged analyzer onboarding. These are local,
authored checks. They do not establish public release, independent adoption,
high-depth cohort behavior or clinical validity. The verified report separates
the exact candidate evidence from later documentation and release-inventory work.
