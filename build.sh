#!/bin/sh
# Assemble the single-file page from its three parts.
cd "$(dirname "$0")"
cat src/head.html src/sim.js src/ui.js > plates.html
echo "wrote plates.html ($(wc -c < plates.html) bytes)"
