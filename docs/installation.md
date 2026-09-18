# Install Rosalind

Choose the installation that matches the workflow you need.

| Version | Availability | What these guides can assume |
|---|---|---|
| 0.1.0 stable | Public historical release | Its own versioned documentation; not the new evidence workflow |
| 0.5.0 source preview | Build from this repository | Evidence, panel QC, saved datasets, and the current SDK |
| Candidate artifacts | Only an explicitly provided and identified artifact | The exact candidate's commands and validation record; not stable publication |

The [stable release documentation](https://github.com/logannye/rosalind/tree/v0.1.0)
describes the older interface. For the workflows linked from the current README,
use the source instructions below. Neither `cargo install rosalind-bio` nor
`pip install rosalind-bio` should be assumed to install the unpublished evidence
preview.

## Build the current CLI

You need Git, Rust/Cargo **1.83 or newer**, a C/C++ toolchain, CMake, pkg-config,
and native compression headers. Follow the platform commands in
[build prerequisites](analyzer-sdk.md#build-prerequisites) first. Source builds
download dependencies and require network access on a fresh machine.

```sh
git clone https://github.com/logannye/rosalind.git
cd rosalind
cargo build --locked --bin rosalind
export PATH="$PWD/target/debug:$PATH"
rosalind --version
rosalind analyze evidence --help
rosalind dataset --help
```

The version should identify the source preview you built. The two help commands
confirm the evidence and saved-dataset interfaces are present. Keep this shell's
PATH setting when following the quickstarts; a new shell needs it again.
Record `git rev-parse HEAD` when sharing a source-build result.

Candidate platform workflows cover Linux x86_64 and macOS arm64/x86_64. This is
not a claim that a current candidate has been published for each platform. See
[implementation status](implementation-status.md) for the retained platform
evidence. Linux ARM wheels and Windows installation are outside these guides.

## Install the Python interface

The distribution is `rosalind-bio`; the import is `rosalind`. Its wheel includes
the matching CLI. Follow [Python source installation](../python/README.md) in a
fresh virtual environment using Python 3.9 or newer. Source wheel builds require
the native prerequisites above and Maturin; an already built wheel requires no
compiler.

Run Python from outside the checkout to check the installed package rather than
accidentally importing local source. Keep the wheel and native executable from
the same build. The Python guide includes the maintained example used by wheel
onboarding checks.

## Use an explicitly supplied candidate

For a native bundle, use its included executable and guides. Its
`ONBOARDING-BUNDLE.json` records the candidate source and SDK dependency rendering.
For a wheel, install the exact local wheel path in a fresh virtual environment;
do not replace it with an unqualified registry install. A candidate SDK may still
need the source dependency patches documented in the [SDK guide](analyzer-sdk.md).

Package installation, external SDK installation, scientific validation, and
stable release are separate gates. Consult the artifact's recorded validation
before treating a candidate as interchangeable with another build.

Continue with the [researcher quickstart](researcher-quickstart.md),
[reuse quickstart](reuse-quickstart.md), or [builder quickstart](builder-quickstart.md).
For missing commands or build errors, see [troubleshooting](troubleshooting.md).
