# Receipt trust levels

Receipt Studio and the CLI report separate evidence instead of collapsing every
check into “verified.” These levels accumulate; none silently implies the next.

1. **Receipt intact** — the deterministic claim and any measurement block match
   their recorded self-hashes.
2. **Artifacts match** — every supplied input/output artifact matches a recorded
   content hash. Matching is by BLAKE3 content, not filename or path.
3. **Resource contract honored** — the receipt records a declared budget and a
   within-budget result. With enforcement, upfront refusal and live breach are also
   explicit outcomes.
4. **Reproduced** — a linked, intact reproduction certificate reports
   `REPRODUCED`, names the original `parent_claim`, and records byte-matching output.
5. **Independently attested** — reproduction evidence was produced by an identified
   independent party or controlled environment. Identity is visible evidence, not
   cryptographic proof by itself.
6. **Signed** — a future cryptographic signature authenticates an identified signer.
   Unsigned receipts correctly display this level as unavailable.

Paths are relocatable metadata in schema 3 and newer. A path-only edit can leave
“receipt intact” unchanged; artifact matching still requires the bytes supplied at
verification time. Claim fields, content hashes, replay tokens, identities, and
measurements remain protected by their appropriate hash.

`verify --json` emits trust-report schema 2, published as
[`schema/trust-report-v2.schema.json`](schema/trust-report-v2.schema.json). The
schema-1 consumer fields (`ok`, `claim`, input/output counts, notes, and problems)
remain at the top level with their original meanings; schema 2 adds the independent
`trust` facets. Consumers that do not understand schema 2 may continue reading
those stable fields and ignoring the additive object.
