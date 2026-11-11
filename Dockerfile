# Mei-rs container image: the binary on a slim Debian with `ps` (workflow
# engines such as Nextflow need it). Built against glibc: a musl build is
# much slower with several threads (musl's allocator, #84).
#
#   docker build -t mei-rs .
#   docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data mei-rs scan --bam sample.bam ...
#
# Base images are pinned by digest (#132); the image runs as a non-root user.

FROM rust:1-trixie@sha256:8a4e90c7946b74626488b1446f712c0fe0ed1a596cbd3d871684bc9c5433ae07 AS build
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY benches ./benches
COPY data ./data
COPY README.md LICENSE ./
RUN cargo build --release --locked

FROM debian:trixie-slim@sha256:a72b4433de1cf7d6bae0c4cee74d0547f58601a643b3375f070c51600d63b6b2
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --uid 10001 --no-create-home mei
COPY --from=build /src/target/release/mei-rs /usr/local/bin/mei-rs
LABEL org.opencontainers.image.source="https://github.com/VictorYt/Mei-rs" \
      org.opencontainers.image.licenses="AGPL-3.0-or-later"
USER mei
ENTRYPOINT ["mei-rs"]
