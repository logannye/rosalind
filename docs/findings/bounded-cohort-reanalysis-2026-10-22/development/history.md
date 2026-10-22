# Development failures and superseded runs

These records are development evidence, separate from the final measured run
39816575369. A superseded or cancelled job is not counted as a passing gate.

- Initial source `1185419` had a workflow YAML error (a colon in a plain run scalar).
  Run [39813234945](https://github.com/logannye/rosalind/actions/runs/39813234945)
  could not execute that focused workflow. Commit `e9fa4c2` uses a block scalar.
- The first broad CI on `e9fa4c2` encountered a newly deprecated Rust
  `AtomicU64::fetch_update` call under `-D warnings`. Commit `e86bd88` adds a
  narrow deprecation allowance at the existing call; it preserves the Rust 1.83
  MSRV and cancellation semantics. Final-source Clippy and MSRV checks pass.
- Early passing focused runs [39813358266](https://github.com/logannye/rosalind/actions/runs/39813358266)
  and [39813787043](https://github.com/logannye/rosalind/actions/runs/39813787043)
  preceded benchmark-worker fairness corrections. The later development run
  [39814541120](https://github.com/logannye/rosalind/actions/runs/39814541120)
  also passed the offline demo and saved HG002 check, but preceded the final
  coverage-copy pin and requested-ALT-only columnar projection. Retained early
  measurements are diagnostic and are not pooled with the final trials.
- Source `6e7d9e1`, measured merge `2e98785`, passed native build/tests and the
  existing Python suite, then failed the retained demo in run
  [39815498177](https://github.com/logannye/rosalind/actions/runs/39815498177).
  Its memory-metric guard expected an integer although native receipt measurements
  are decimal strings. The demo refused publication. Commit `9627c93` fixes the
  parser and adds a pure regression; the final full rerun passes. Its report and
  stderr are retained beside this file.
- Superseded CI and wheel matrices were cancelled where still running so the
  latest implementation could use hosted resources. Cancellation is not evidence
  of a software failure or a successful platform gate.

[index.json](index.json) binds the retained development artifacts to their run,
branch head, actual tested merge and archive digest. No scientific or throughput
claim is transferred from an earlier binary to the final source.
