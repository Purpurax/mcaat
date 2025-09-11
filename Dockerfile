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

# Clone your repo
RUN git clone --recurse-submodules https://github.com/RNABioInfo/mcaat.git .

# Build with default release (no option needed)
RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

# For debug build (uncomment if needed):
# RUN mkdir build && cd build && \
#     cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_DEBUG=ON && \
#     make -j$(nproc)

# Stage 2: Runtime
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libomp5 \
    zlib1g \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /app/build/mcaat /usr/local/bin/mcaat

# Set default command
ENTRYPOINT ["mcaat"]