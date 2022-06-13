#!/bin/bash
WASM=./pkg/nbezier_wasm.js
if [ -f "$WASM" ]; then
    echo "Visit http://0.0.0.0:8000/playground/"
    /usr/bin/env python3 -m http.server
else 
    echo "Please run from inside the wasm crate root"
    echo "Make sure you have compiled it using wasm-pack --target web"
fi
