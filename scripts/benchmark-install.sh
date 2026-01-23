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

echo "Installing Columba..."
if ! command -v columba &> /dev/null; then
    cd /tmp
    if [ ! -d "columba" ]; then
        git clone https://github.com/biointec/columba.git
    fi
    cd columba
    
    bash build_script.sh Vanilla
    
    $SUDO cp build_Vanilla/columba /usr/local/bin/
    $SUDO cp build_Vanilla/columba_build /usr/local/bin/
    cd ..
    echo "Columba installed successfully"
else
    echo "Columba already installed"
fi

echo ""
echo "=== Installation Complete ==="
echo "Installed tools:"
command -v samtools && samtools --version | head -n1
command -v dwgsim && echo "dwgsim installed" || echo "dwgsim not found"
command -v minimap2 && minimap2 --version
command -v bwa-mem2 && bwa-mem2 version
command -v columba && echo "Columba installed" || echo "Columba not found"
command -v columba_build && echo "Columba_build installed" || echo "Columba_build not found"

echo ""
echo "All mapping tools are ready for benchmarking!"
