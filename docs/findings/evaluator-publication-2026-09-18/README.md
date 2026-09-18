# Published evaluator recovery — September 18, 2026

[Run 35378796613](https://github.com/logannye/rosalind/actions/runs/35378796613)
built and published the pinned linux/amd64 hap.py evaluator, passed the offline
hap.py 0.3.15 / RTG 3.12.1 version checks and identical-SNV vcfeval smoke, and
published its SBOM and GitHub provenance attestation. The run concluded **failed**
because its final automatic PR step required a base branch for detached HEAD;
the failure is retained in [automation-failure.txt](automation-failure.txt).
Repository Actions PR creation is also disabled and remains unchanged.

The exact generated lock was recovered from the uploaded `happy-image-evidence`
artifact. Its source/dependency lock is unchanged. The resulting image is:

```text
ghcr.io/logannye/rosalind-happy@sha256:e2a6bed76ec740eec691a9383790b130f0882705b0aca4f8563e74cd233cc6cf
```

Source commit: `92ac1ccd7bb898060a7e6bb7297411cde895a3fd`.
Source lock: `fdcf635037992264eb95dae5b3d1c6a7c36ebb702e586b67f82a64bad4d55a7c`.
Anonymous registry manifest retrieval returned HTTP 200 (maintainer observation).
Fresh cryptographic verification succeeded with the exact source and workflow:

```sh
gh attestation verify oci://ghcr.io/logannye/rosalind-happy@sha256:e2a6bed76ec740eec691a9383790b130f0882705b0aca4f8563e74cd233cc6cf \
  --repo logannye/rosalind \
  --source-digest 92ac1ccd7bb898060a7e6bb7297411cde895a3fd \
  --signer-workflow logannye/rosalind/.github/workflows/happy-image.yml
```

[Verified provenance](verified-provenance.json) retains the verified certificate
identity, signed subject digest, source commit, and workflow invocation. The full
attestation remains retrievable from GitHub/OCI; this excerpt alone is not a
signature-verification substitute. [Smoke output](smoke.json) is one synthetic
identical SNV with one true positive and no false positives or negatives.
**This is evaluator readiness, not scientific HG002/GIAB validation.**

The workflow now uploads the same evidence and writes a manual PR handoff summary.
It does not need repository write/PR permissions or change the repository's
Actions settings. Review the artifact, verify provenance and source identity,
then submit its exact lock through an ordinary maintainer-created PR.

After this lock is reviewed and merged, use a clean, current main checkout:

```sh
cargo xtask giab image status --json
cargo xtask giab benchmark plan --output /tmp/rosalind-giab-plan.json --json
# Review the plan and use its exact plan_id only when authorizing the real run:
cargo xtask giab benchmark dispatch --plan /tmp/rosalind-giab-plan.json --confirm PLAN_ID
cargo xtask giab benchmark status --json
```

The real benchmark prepares pinned HG002 chr20 data and runs the external
evaluator. Download and verify the resulting attested candidate, retain failures,
and review any baseline proposal separately. No image build, publication, or
scientific benchmark was launched for this recovery PR.

Full downloaded job-log SHA-256: `5eda486f5481980c6f955b1ef6dcee0a861245b40fa5aba6b2719ca707ca1d8d`.
