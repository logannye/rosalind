# Historical context

Rosalind's original motivation was to exchange additional computation and I/O for
lower memory use. Today's evidence engine applies that idea concretely through
resource-driven genomic tiles, indexed rereads, bounded canonical batches and
admitted worker counts. Successful evidence retains the same scientific meaning;
an unsupported layout or insufficient budget can still cause refusal. This does
not establish a universal square-root-memory algorithm or an internal hard
allocation sandbox.

The current commitments live in the [capability inventory](implementation-status.md),
[resource contract](../CONTRACT.md) and [delivery roadmap](ROADMAP.md).

The following material is retained as historical context. Its archive banners
remain authoritative about its historical status:

- [Original research thesis and open problems](https://github.com/logannye/rosalind/blob/main/docs/OPEN_PROBLEMS.md)
- [Production and growth snapshot](https://github.com/logannye/rosalind/blob/main/docs/GROWTH.md)
- [Original target architecture](https://github.com/logannye/rosalind/blob/main/docs/superpowers/specs/2026-05-26-rosalind-target-architecture.md)
- [Dated design specifications](https://github.com/logannye/rosalind/tree/main/docs/superpowers/specs)

Historical clinical, field, broad-modality and never-refusing scenarios are
research possibilities, not supported workflows. New work needs an explicit
scientific contract, measured validation and a useful user task before it becomes
a product commitment. Dated implementation specifications may describe a shipped
subsystem, but current API and semantics references take precedence.
