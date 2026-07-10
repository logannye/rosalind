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

New receipts set `replay_schema=2` and record `command_argv` as a canonical JSON
array encoded in a string. Tokens are never interpreted by a shell, so spaces in
paths and option values are unambiguous. The legacy human-readable `command` remains
present. Schema-5 receipts without `command_argv` replay through the historical
`command` parser.

`rosalind reproduce --binary PATH` is the only way to select a third-party replay
executable. A receipt's recorded `producer.binary` is informational and is never
executed automatically.

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
