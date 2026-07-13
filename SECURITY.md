# Security policy

## Supported versions

Security fixes are applied to the latest published release and the default branch.
Pre-release branches may change without compatibility guarantees.

## Reporting a vulnerability

Please use GitHub's private security-advisory reporting flow for this repository.
Do not open a public issue for a vulnerability that could expose genomic data,
execute untrusted receipt content, bypass output/refusal guarantees, or accept a
corrupt reference/index/receipt. Include affected versions, reproduction steps,
impact, and any proposed mitigation. Maintainers will acknowledge a complete report
as soon as practical and coordinate disclosure after a fix is available.

## Security boundaries

- Rosalind is research software and is not validated for clinical diagnosis.
- Receipts are tamper-evident, not proof of authorship; signing is not implemented.
- Replay is tokenized and allowlisted, but users must independently trust any
  external binary supplied with `--binary`.
- Receipt Studio is loopback-only and browser inspection is local; Rosalind has no
  telemetry or automatic genomic-data transmission.
- `--enforce` is cooperative unless an existing cgroup-v2 limit is required. It is
  not a process sandbox.
