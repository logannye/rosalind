# RC2 container prerequisite failure and prepublication regression check

The 0.5.0-rc.2 container build failed before an image was produced from candidate
[`8874508`](https://github.com/logannye/rosalind/commit/8874508844e106484fa290b04e3fc6349686beaf).
The [original job](https://github.com/logannye/rosalind/actions/runs/35379805820/job/105726566760)
and [retained error excerpt](rc2-failure.txt) remain separate from corrected-run
evidence. The failure was `libz-sys v1.1.23` invoking the missing `cmake` executable
to build zlib-ng. The Rust `cmake` crate was present; its native tool was absent.
The [failure record](rc2-failure.json) identifies the candidate, tag and full
downloaded job-log hash.

The builder now explicitly installs the same Linux compiler and compression
prerequisites as the SDK/native release routes: build-essential, clang, CMake,
pkg-config, bzip2, lzma and zlib headers. Runtime compression libraries are
explicitly installed too. The pinned base images and Debian snapshot remain
unchanged.

Ordinary pull-request CI now builds the actual release Dockerfile for linux/amd64
without registry login, push, a protected environment or write-token permissions.
It runs the image with networking disabled, capabilities dropped and privilege
escalation disabled, verifies the default UID10001 and exact package version,
checks current evidence/dataset commands, and completes the embedded demo's
native processing, artifact verification and byte replay. The job has a 30-minute
timeout and retains build/smoke output for both success and failure.

Local shell and workflow lint can run without Docker. This workstation's Docker
daemon was unavailable, so an actual successful build must be established by the
new CI job before publication proceeds. The demo is a packaging check; it does
not establish independent use, representative resource performance or new
scientific agreement. Publication workflows and protected release settings are
unchanged by this fix.
