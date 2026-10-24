// Level 0 in a Web Worker, for the browser gate (scripts/browser-gate.js).
import { generatePlanet } from '../../src/world/l0.js';
import { fieldHash } from '../../src/world/hash.js';

const hashes = (t) => Object.fromEntries(Object.entries(t.fields).map(([k, v]) => [k, fieldHash(v)]));

self.onmessage = (e) => {
  const { seed, options, golden } = e.data;
  const out = {};
  // the timed run goes first, in a fresh worker, as a user's first planet would
  const t0 = performance.now();
  const tile = generatePlanet(seed, options || {});
  out.seconds = (performance.now() - t0) / 1000;
  if (golden) out.golden = hashes(generatePlanet(golden.seed, golden.options));
  out.hashes = hashes(tile);
  out.cells = tile.stats.cells;
  out.events = tile.events.length;
  self.postMessage(out);
};
