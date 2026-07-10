# Receipt compatibility and extension contract

Rosalind writes one stable top-level envelope: `inputs`, optional `measurements`,
`outputs`, `params`, `subcommand`, and `tool_version`. The published schema for the
current envelope is [`schema/receipt-v5.schema.json`](schema/receipt-v5.schema.json).
Additive claim parameters do not require a top-level schema bump.

## What each hash protects

- `manifest_blake3` protects the deterministic claim: content hashes, parameters,
  producer/build identity, analyzer identity, replay recipe, and output hashes.
- `measurement_blake3` independently protects machine-local measurements such as
  peak RSS and the resource-contract verdict.
- For schema 3 and newer, `inputs[].path` and `outputs[].path` are portable lookup
  metadata. Paths remain in the file for humans and native artifact checks, but are
  deliberately excluded from the claim hash. Relocating data therefore preserves
  the claim; changing a content hash does not.

Receipt integrity is tamper-evidence, not authorship. A party able to rewrite and
reseal an unsigned receipt can create a new internally consistent claim. Signatures
are a separate trust layer.

## Parameter namespaces

- Unprefixed keys and Rosalind-owned namespaces such as `model.*` are reserved for
  Rosalind.
- `producer.name`, `producer.version`, `producer.repository`, and
  `producer.binary` identify the binary that created a run.
- `analyzer.id` and `analyzer.version` identify the analyzer implementation.
  Analyzer-owned parameters should use `analyzer.*`.
- Third-party extensions must use `x.<reverse-dns>.*`, for example
  `x.org.example.filter.window`.

Unknown additive parameters must be preserved by general-purpose receipt tools and
are part of the deterministic claim.

## Replay compatibility

New receipts set `replay_schema=3`, record `replay.kind=rosalind|external-analyzer`,
and store `command_argv` as a canonical JSON array encoded in a string. Tokens are
never interpreted by a shell, so spaces and Unicode in paths and option values are
unambiguous. Before execution, Rosalind requires an intact claim and measurements,
then cross-checks the subcommand, argv, parameters, input/output markers, and hashes.
Built-in replay is allowlisted; external-analyzer replay requires schema 3 and an
explicit `--binary PATH`. A dry run returns the validated execution plan without
starting a child process.

The legacy human-readable `command` remains present. Applicable historical Rosalind
schema-5 recipes continue through strict legacy validation. Legacy third-party
recipes are reported `INCONCLUSIVE` with upgrade guidance because they do not carry
enough information for safe automatic execution.

`rosalind reproduce --binary PATH` is the only way to select a third-party replay
executable. A receipt's recorded `producer.binary` is informational and is never
executed automatically. Replay runs in a fresh working directory with redirected
`HOME` and `TMPDIR`, a minimal locale/timezone, and only the caller's existing
`PATH`. This contains ordinary path side effects; trusting an explicitly selected
binary is still required, and the working directory is not a full OS sandbox.

## Historical capabilities

| Schema | Capability added | Compatibility rule |
|---|---|---|
| 1 | Claim self-hash | Paths and measurements are inside the claim. |
| 2 | Claim/measurement split | Measurements have an independent hash. |
| 3 | Content-only portable claim | Recorded paths no longer affect claim identity. |
| 4 | Build identity | Git, dirty state, rustc, target, and lockfile hash enter the claim. |
| 5 | Replay recipe | `reproduce`, causal diff, and chain analysis can reconstruct a run. |

Canonical fixtures for all five versions live in
[`crates/receipt/fixtures`](../crates/receipt/fixtures) and are regression-tested.
