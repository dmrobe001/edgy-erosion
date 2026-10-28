#!/bin/sh
# Assemble the single-file page from its parts; BUILD_ID is a hash of the sources, recorded in session exports.
cd "$(dirname "$0")"
id=$(cat src/head.html src/sim.js src/mesh.js src/l1.js src/l2.js src/ui.js | sha1sum | cut -c1-12)
{ cat src/head.html; echo "const BUILD_ID=\"$id\";"; cat src/sim.js src/mesh.js src/l1.js src/l2.js src/ui.js; } > plates.html
echo "wrote plates.html ($(wc -c < plates.html) bytes, build $id)"
