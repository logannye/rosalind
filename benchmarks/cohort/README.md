# Cohort CI evidence

The Python CI job prepares the authored three-sample fixture and executes
[`run_cohort.py`](../../examples/cohort-reanalysis/run_cohort.py). It retains the
whole generated directory on success or failure, including commands, native
outputs, receipts, oracle comparisons, and measured import/reuse/extension costs.
These are recurring engineering checks, not independent adoption evidence.

`cgroup_probe.py` then exercises saved-only `cohort extract` and `cohort summarize`
in fresh Linux/amd64 Docker containers. The probe uses the immutable Rust Bookworm
runtime already pinned in the repository, sets `memory.max` to 512 MiB and swap to
zero, disables networking, drops capabilities, and mounts only the binary, cohort,
candidate VCF, and a new output directory. Original alignments and references are
not mounted. The Python CI runner is Ubuntu 22.04 with a separate Cargo cache so the host-built binary's glibc
baseline remains compatible with the pinned Bookworm runtime.

For each operation, an outside-probe-container run supplies exact TSV comparison
bytes. The constrained runs require `--enforce --require-os-limit`, native receipt
verification, matching output bytes, and a recorded zero original-alignment decode
count. A separate 1 MiB declared-budget invocation must exit 3 without successful
or partial artifacts. All containers retain actual before/after `memory.events`,
`memory.max`, `memory.swap.max`, `memory.peak`, native/container exit statuses,
Docker inspection including `OOMKilled`, image identity, commands and full logs.
Missing controller evidence fails the probe. No OS-enforcement result is claimed
until an actual Linux CI run passes and its retained report is linked.

```sh
python3 -m unittest discover -s benchmarks/cohort -p 'test_*.py'
python3 benchmarks/cohort/cgroup_probe.py \
  --binary /absolute/path/to/linux-amd64/rosalind \
  --demo-report /path/to/prepared/cohort-validation/report.json \
  --output /new/path/to/cohort-cgroup --docker-context default
```

The output directory must be new. Pulling the immutable image uses the Docker
client; analytical containers have no network access. The harness removes only
containers it created and leaves all logs/reports intact. It does not reconfigure
Docker or host cgroups. This tiny synthetic probe is a failure-contract check;
it establishes neither representative cohort performance nor universal low-memory
execution. Cgroup memory includes controller and page-cache charges, whereas
Rosalind receipts measure native process RSS.
