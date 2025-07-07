#!/bin/bash

set -e

INSTALL=false

# Parse arguments
for arg in "$@"; do
    case $arg in
        --install)
            INSTALL=true
            shift
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: ./install_and_build.sh [--install]"
            exit 1
            ;;
    esac
done


echo "🔧 Installing system dependencies..."
sudo apt update
sudo apt install -y \
  build-essential \
  cmake \
  git \
  zlib1g-dev \
  libomp-dev \
  libpthread-stubs0-dev

echo "📦 Cloning submodules..."
git submodule update --init --recursive

echo "🏗️  Creating build directory..."
mkdir -p build
cd build

echo "⚙️  Running CMake..."
cmake ..

echo "🚀 Building the project..."
make -j$(nproc)

echo "✅ Build complete. Binary is located at: build/mcaat"

if $INSTALL; then
    echo "📥 Installing binary to /usr/local/bin..."
    sudo cp mcaat /usr/local/bin/
    echo "✅ Installed. You can now run 'mcaat' from anywhere."
fi
