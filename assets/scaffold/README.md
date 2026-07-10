# __PACKAGE_NAME__

A standalone Rosalind `ColumnAnalyzer` that explicitly declares its retained-memory
model and inherits bounded execution, transactional output, safe replay, and a
content-addressed receipt.

```sh
cargo build --release
rosalind demo --output-dir /tmp/rosalind-demo
target/release/__PACKAGE_NAME__ run \
  --index /tmp/rosalind-demo/ref.idx \
  --alignments /tmp/rosalind-demo/sorted.bam \
  --memory-budget-mb 128 --enforce -o output.tsv
rosalind receipt inspect --manifest output.tsv.manifest.json --artifact output.tsv
rosalind reproduce --manifest output.tsv.manifest.json \
  --inputs /tmp/rosalind-demo --binary target/release/__PACKAGE_NAME__ --dry-run
rosalind conformance analyzer --binary target/release/__PACKAGE_NAME__ --json
```

`--force` opts into atomic replacement. `--require-os-limit` additionally requires
an already-active Linux cgroup-v2 `memory.max` at or below the declared budget.

Run `bash scripts/contract-check.sh` for repeat determinism, verification,
sanitization, safe replay, reproduction, refusal, breach, diff, and collision checks.
