FROM rust:1.95-bookworm@sha256:6258907abe69656e41cd992e0b705cdcfabcbbe3db374f92ed2d47121282d4a1 AS build
RUN apt-get update && apt-get install -y --no-install-recommends clang libbz2-dev liblzma-dev pkg-config && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN cargo build --release --locked --bin rosalind && strip target/release/rosalind

FROM python:3.12.10-slim-bookworm@sha256:fd95fa221297a88e1cf49c55ec1828edd7c5a428187e67b5d1805692d11588db
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools=1.16-1 samtools=1.16.1-1 time=1.9-0.2 jq=1.6-2.1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY benchmarks/platform/requirements.lock /tmp/requirements.lock
RUN pip install --no-cache-dir --requirement /tmp/requirements.lock
COPY --from=build /src/target/release/rosalind /usr/local/bin/rosalind
COPY benchmarks/platform /opt/rosalind-platform
ENV PYTHONHASHSEED=0 LC_ALL=C LANG=C
WORKDIR /work
ENTRYPOINT ["/opt/rosalind-platform/benchmark.sh"]
