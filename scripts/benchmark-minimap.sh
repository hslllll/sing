#!/bin/bash
set -e

THREADS=$1
MODE=$2
if [ -z "$MODE" ]; then
    echo "Usage: $0 <threads> [h|y|a|m|b]"
    exit 1
fi

python3 benchmark.py "$THREADS" minimap "$MODE"
