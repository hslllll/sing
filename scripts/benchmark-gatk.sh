#!/usr/bin/env bash
set -e

THREADS=$1
MODE=$2
if [ -z "$MODE" ]; then
    echo "Usage: $0 <threads> [y|a|m|h|b]"
    exit 1
fi

python3 benchmark.py "$THREADS" gatk "$MODE"
