FROM rust:1.95-bookworm@sha256:6258907abe69656e41cd992e0b705cdcfabcbbe3db374f92ed2d47121282d4a1 AS build
RUN apt-get update && apt-get install -y --no-install-recommends clang libbz2-dev liblzma-dev pkg-config && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN cargo build --release --locked --bin rosalind && strip target/release/rosalind

FROM debian:bookworm-slim@sha256:60eac759739651111db372c07be67863818726f754804b8707c90979bda511df
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates libbz2-1.0 liblzma5 && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 10001 rosalind
COPY --from=build /src/target/release/rosalind /usr/local/bin/rosalind
USER 10001:10001
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/rosalind"]
CMD ["--help"]
