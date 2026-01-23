#!/bin/bash
set -e

cd "$(dirname "$0")/../sing-web"

echo "Building sing-web with wasm-pack..."
wasm-pack build --release --target web --out-dir pkg

echo "Build complete! Output is in the 'sing-web/pkg' directory."
