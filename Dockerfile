# Mei-rs container image: the binary on a slim Debian with `ps` (workflow
# engines such as Nextflow need it). Built against glibc: a musl build is
# much slower with several threads (musl's allocator, #84).
#
#   docker build -t mei-rs .
#   docker run --rm -v "$PWD":/data -w /data mei-rs scan --bam sample.bam ...

FROM rust:1-trixie AS build
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY benches ./benches
COPY data ./data
COPY README.md LICENSE ./
RUN cargo build --release --locked

FROM debian:trixie-slim
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/mei-rs /usr/local/bin/mei-rs
LABEL org.opencontainers.image.source="https://github.com/VictorYt/Mei-rs" \
      org.opencontainers.image.licenses="AGPL-3.0-or-later"
ENTRYPOINT ["mei-rs"]
