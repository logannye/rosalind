## Determinism Contract

Rosalind’s primary product guarantee is **bit-for-bit reproducibility**: given identical inputs and configuration, Rosalind produces identical outputs across runs.

This document defines what is (and is not) covered by that guarantee, and the engineering constraints required to preserve it.

### Definitions

- **Deterministic output**: the emitted bytes of primary artifacts (BAM/VCF/JSON manifests) are identical across repeated runs with the same inputs and configuration, on the same Rosalind version and target platform.
- **Canonicalization**: where file formats permit multiple equivalent encodings (e.g., header ordering), Rosalind defines and emits a canonical form.
- **Stable ordering**: when results are produced from unordered collections or parallel execution, Rosalind explicitly defines a stable sort order and uses it everywhere.

### Covered inputs

Determinism is guaranteed with respect to:

- **Reference**: the exact FASTA bytes (including contig names and sequences) used for indexing and calling.
- **Reads**: the exact FASTQ bytes (including read names and base qualities).
- **Configuration**: all CLI flags / config file values, including thread count and memory limits.
- **Rosalind version**: the exact crate version and build profile.

### Covered outputs (v1 target)

Rosalind aims to guarantee deterministic bytes for:

- **BAM/CRAM**: alignment records and headers emitted by Rosalind.
- **VCF**: variant records and headers emitted by Rosalind.
- **Run manifest**: a machine-readable file recording inputs, parameters, and versions.

### Not covered (unless explicitly stated)

- **Wall-clock timestamps** in logs or headers (Rosalind should avoid emitting them into primary artifacts).
- **OS-dependent path prefixes** in provenance unless normalized.
- **Nondeterministic hardware behavior** outside our control (e.g., bit flips, kernel bugs).

### Engineering rules (hard requirements)

To preserve determinism, core pipeline code must follow these rules:

1. **No unordered iteration in output-critical paths**:
   - Do not rely on iteration order of `HashMap`/`HashSet`.
   - Prefer `BTreeMap`/`BTreeSet` or sort keys before emitting results.
2. **Explicit tie-breakers**:
   - Whenever two candidates have the same primary score, define secondary keys (position, strand, read name, etc.).
3. **Deterministic sorting**:
   - Use stable sorts where key collisions are possible.
   - Define a single canonical ordering for records (e.g., coordinate sort with full tie-break key).
4. **Deterministic parallelism**:
   - Parallel execution is allowed only when results are combined via an explicitly ordered reduction.
   - Otherwise, default to single-threaded execution for output-critical stages.
5. **Numerics**:
   - Avoid `f32` for ranking/tie-breaks in persisted outputs.
   - When floating point is necessary, use `f64` and ensure reductions are performed in a stable order.
6. **Canonical output emission**:
   - Canonicalize headers and record ordering.
   - Do not emit run-dependent data (timestamps, randomized IDs) into primary artifacts.

### Determinism testing

Rosalind should maintain:

- **Repeat-run tests**: run the same pipeline twice and compare output bytes.
- **Thread-count invariance tests** (when multithreading is added): `--threads 1` vs `--threads N` must match.

### The memory governor and determinism

Under `--enforce`, a background governor thread reads process RSS and, on a breach,
cooperatively aborts the run (exit 4). This does **not** weaken determinism: the guard
never touches output bytes, and a run that fits never fires it — identical inputs still
produce a byte-identical VCF/feature table. A *breach* is a non-deterministic abort
(which locus trips depends on timing), but a breach exits non-zero and is a failure, not
a reproducible artifact, so no determinism guarantee is affected. (The receipt's
`peak_rss_bytes` / `baseline_rss_bytes` / `rss_residual_bytes` are machine-dependent
measurements, like any realized-memory field — the primary output stays byte-identical.)


