#!/bin/sh
# Assemble the single-file page from its four parts.
cd "$(dirname "$0")"
cat src/head.html src/sim.js src/l1.js src/ui.js > plates.html
echo "wrote plates.html ($(wc -c < plates.html) bytes)"
