#!/bin/bash
set -e

cd "$(dirname "$0")/../sing-browser"

echo "Building sing-browser with wasm-pack..."
wasm-pack build --release --target web --out-dir pkg

echo "Build complete! Output is in the 'sing-browser/pkg' directory."
