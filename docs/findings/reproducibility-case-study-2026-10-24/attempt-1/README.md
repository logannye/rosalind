# First attempt: existing nested-build path assumption

[Run 40151648804](https://github.com/logannye/rosalind/actions/runs/40151648804)
tested source `385e65b024a4cf83158290c50aa1f3989109765f` and failed in the existing
Rust generated-analyzer conformance test before release build, Python gates,
NA18507 preflight, HG002 acquisition or scientific measurements. Replay was
unstarted. This is an infrastructure/test-location failure, not scientific
agreement or a performance result.

`tests/conformance.rs` invoked Cargo with inherited `CARGO_TARGET_DIR`, but
looked for the generated analyzer in the scaffold's private `target/debug`.
The correction explicitly passes the one shared target directory to nested
Cargo and resolves the analyzer from that same directory. No measurement code,
query definition, interval, filtering policy, native engine, or gate is weakened.
The entire locked workspace test suite and exact-source build must run again.

The first producer artifact is `12583196065`, 20,965 compressed bytes, provider
SHA256 `b51e013fd158d2d9a3874c0c6e3f86f8240163faf15e3f53fe2752fcc8fc3da4`.
Provider run/jobs/artifact records, resource preflight and full failing native
test log are preserved alongside this note. The job's resource preflight passed;
no local native build or executable test was used to repair it.

This consumes the first attempt. The next hosted run is the single permitted
correction attempt under the frozen case-study protocol; no performance-selection
rerun is authorized.
