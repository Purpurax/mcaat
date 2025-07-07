# Stage 1: Build
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    zlib1g-dev \
    libomp-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Clone your repo (or mount it during build)
RUN git clone --recurse-submodules https://github.com/RNABioInfo/mcaat.git .

# If needed:
# RUN git submodule update --init --recursive

RUN mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc)

# Stage 2: Runtime
FROM debian:bookworm-slim

# Install only runtime dependencies
RUN apt-get update && apt-get install -y \
    libomp5 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /app/build/mcaat /usr/local/bin/mcaat

# Set default command
ENTRYPOINT ["mcaat"]
