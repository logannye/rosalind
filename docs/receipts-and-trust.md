# Receipts and trust

Rosalind receipts are portable, tamper-evident execution records. Claim-protected
fields include input/output content, deterministic parameters, producer/analyzer
identity, output-affecting format choices, the execution plan, and run status.
Machine-local RSS and enforcement observations live in a separately protected
measurement block so they do not change the portable claim.

Receipt integrity does not establish authorship. Someone who can alter an unsigned
receipt can recompute its hashes. Independent trust comes from matching artifacts,
reproducing bytes with independently obtained code and inputs, or an external
attestation system. Rosalind does not currently sign receipts.

Replay validates a tokenized argv plan, content-locates inputs, uses an allowlist for
built-ins, requires explicit `--binary` selection for external analyzers, and never
executes a recorded shell string. Modern recorded paths are relocatable metadata.
Historical schemas 1–5 parse and verify at their original capability level.

Receipt Studio and the browser verifier process data locally and make no automatic
third-party requests. Rosalind has no automatic telemetry.
