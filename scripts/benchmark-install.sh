#!/bin/bash
set -e

echo "=== Installing Mapping Tools ==="

if [[ $EUID -eq 0 ]]; then
    SUDO=""
else
    SUDO="sudo"
fi

echo "Updating package lists..."
$SUDO apt-get update

echo "Installing dependencies..."
$SUDO apt-get install -y wget curl git build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libncurses-dev

echo "Installing samtools..."
if ! command -v samtools &> /dev/null; then
    $SUDO apt-get install -y samtools
    echo "samtools installed successfully"
else
    echo "samtools already installed"
fi

echo "Installing bcftools..."
if ! command -v bcftools &> /dev/null; then
    $SUDO apt-get install -y bcftools
    echo "bcftools installed successfully"
else
    echo "bcftools already installed"
fi

echo "Checking for rustup (Rust toolchain)..."
if ! command -v rustup &> /dev/null; then
    echo "rustup not found. Installing rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    export PATH="$HOME/.cargo/bin:$PATH"
    source "$HOME/.cargo/env" || true
    echo "rustup installed successfully"
else
    echo "rustup already installed"
fi

echo "Installing bowtie2..."
if ! command -v bowtie2 &> /dev/null; then
    $SUDO apt-get install -y bowtie2
    echo "bowtie2 installed successfully"
else
    echo "bowtie2 already installed"
fi

echo "Installing dwgsim..."
if ! command -v dwgsim &> /dev/null; then
    cd /tmp
    if [ ! -d "DWGSIM" ]; then
        git clone --recursive https://github.com/nh13/DWGSIM.git
    fi
    cd DWGSIM
    make
    $SUDO cp dwgsim /usr/local/bin/
    cd ..
    echo "dwgsim installed successfully"
else
    echo "dwgsim already installed"
fi

echo "Installing minimap2..."
if ! command -v minimap2 &> /dev/null; then
    cd /tmp
    if [ ! -d "minimap2-2.30_x64-linux" ]; then
        curl -L https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 | tar -jxvf -
    fi
    $SUDO cp /tmp/minimap2-2.30_x64-linux/minimap2 /usr/local/bin/
    $SUDO chmod +x /usr/local/bin/minimap2
    echo "minimap2 installed successfully"
else
    echo "minimap2 already installed"
fi

echo "Installing bwa-mem2..."
if ! command -v bwa-mem2 &> /dev/null || [ ! -f "/usr/local/bin/bwa-mem2.avx2" ]; then
    cd /tmp
    if [ ! -d "bwa-mem2-2.2.1_x64-linux" ]; then
        curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
    fi
    $SUDO cp -f /tmp/bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/
    $SUDO chmod +x /usr/local/bin/bwa-mem2*
    echo "bwa-mem2 installed successfully"
else
    echo "bwa-mem2 already installed"
fi

echo ""
echo "=== Installation Complete ==="
echo "Installed tools:"
command -v samtools && samtools --version | head -n1
command -v bcftools && bcftools --version | head -n1
command -v dwgsim && echo "dwgsim installed" || echo "dwgsim not found"
command -v minimap2 && minimap2 --version
command -v bwa-mem2 && bwa-mem2 version
echo ""
echo "All mapping tools are ready for benchmarking!"
