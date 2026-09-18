# Native evidence under Linux cgroup limits

`cgroup_probe.py` runs the real evidence CLI in fresh Docker containers with
`memory.max` set explicitly and swap disabled. It uses an explicitly named
Docker context, copies and rehashes inputs in a temporary volume, disables container network
access, and removes only its own containers and volume after retaining results.
It never changes Docker, Colima, host cgroups, or daemon settings.

Requirements: Python3.9+, Docker API1.49+ (`image inspect --platform`), a cgroup-v2
Linux daemon, and an already available Linux/amd64 image containing `sh`, `awk`
and the native binary's runtime libraries. Supply a Linux/amd64 Rosalind binary,
an uncompressed indexed FASTA, indexed BAM or CRAM, and exactly one VCF or BED.
The files must describe a real nonempty workload with reads longer than one base,
because one scenario intentionally declares a one-base read limit.

```sh
python3 benchmarks/evidence/cgroup_probe.py \
  --docker-context YOUR_ISOLATED_CONTEXT --image rust:1.95-bookworm \
  --binary /path/to/linux-amd64/rosalind \
  --reference /path/to/slice.fa --alignments /path/to/slice.bam \
  --sites /path/to/candidates.vcf --output /tmp/cgroup-probe \
  --memory-mib 128 --native-oom-mib 8
```

The image tag is resolved and execution uses its immutable repository digest.
Use `--alignment-index` for a nonstandard BAI, CSI or CRAI path. The same command
supports CRAM with the explicit local FASTA. Output directories are create-new.
The default evidence projection is `depths,alleles`; `--fields` selects another
recorded projection. This probe is separate from the representative timing curve.

The default cases are cooperative completion, completion with
`--require-os-limit`, preflight refusal of a declared record envelope larger than
the budget, an observed one-MiB startup-budget breach, a declared read-capacity failure,
and a separately labeled allocation control that triggers the kernel OOM killer
in a32MiB container. `--native-oom-mib 8` additionally attempts an intentionally
undersized native invocation and requires kernel OOM evidence. This is a test of
failure handling, not an admitted execution strategy. If a selected native limit
does not trigger OOM, the probe fails rather than claiming it did.

`report.json` and each case directory retain actual `memory.max`,
`memory.swap.max`, before/after `memory.events`, `memory.peak` when available,
the native and container exit statuses, complete Docker inspection including
`OOMKilled`, exact commands, stdout/stderr, input/binary/image identities, and
every output/partial/receipt file's size and SHA256. Completed artifacts are
independently verified using the native CLI and compared by their physical bytes.
Capacity failures must have an integrity-checked receipt explicitly identifying
the partial. First-party CLI failure receipts retain the requested receipt path;
their failed status and partial output identity distinguish them from completed
artifacts. An exit137 alone is insufficient evidence of kernel OOM.

An OS SIGKILL cannot execute cooperative cleanup or promise a partial receipt;
the probe requires failed invocations to leave no successful destination and
retains any unpublished staging files. Missing controller observations remain
unknown; the report does not invent `memory.events` after a killed controller.
Docker memory charges include the small controller and page cache, whereas the
Rosalind receipt reports its process RSS. The allocation control is not a
Rosalind scientific workload, and this harness makes no throughput claims.

Run the outcome-classification regressions without Docker:

```sh
python3 -m unittest discover -s benchmarks/evidence -p test_cgroup_probe.py -v
```
