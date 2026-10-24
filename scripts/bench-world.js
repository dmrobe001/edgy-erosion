#!/usr/bin/env node
// Time the hierarchical generator's stages against the design's budgets:
//   node scripts/bench-world.js            level 0 (full, ~10 s) and the level-1 solver
//   node scripts/bench-world.js --quick    level-1 solver only
import { generatePlanet, L0_DEFAULTS } from '../src/world/l0.js';
import { RasterSolver } from '../src/world/fastscape.js';
import { worldKey, deriveKey, uniform } from '../src/world/hash.js';

const quick = process.argv.includes('--quick');
if (!quick) {
  const t0 = performance.now();
  const t = generatePlanet(1);
  const s = (performance.now() - t0) / 1000;
  console.log(`level 0: ${t.stats.cells} cells, ${(t.stats.spacing / 1000).toFixed(0)} km apart, ` +
    `${(L0_DEFAULTS.duration / 1e6).toFixed(0)} Myr in ${L0_DEFAULTS.duration / L0_DEFAULTS.dt} steps: ${s.toFixed(1)} s (budget 10 s); ` +
    `${t.events.length} events, ${t.rivers.nodes.length} river nodes`);
}

// level-1 solver: a block rising at 1 mm/yr, walls on three sides, draining south
const key = deriveKey(worldKey(1), 7);
for (const N of [128, 256, 512, 1024]) {
  const s = new RasterSolver(N, N, 1000);
  const z = new Float64Array(N * N), fixed = new Uint8Array(N * N), outlet = new Uint8Array(N * N);
  for (let i = 0; i < N * N; i++) {
    const c = i % N, r = (i - c) / N;
    z[i] = 10 * uniform(key, i) + 3 * r;
    if (c === 0 || r === 0 || c === N - 1 || r === N - 1) fixed[i] = 1;
    if (r === 0) { outlet[i] = 1; z[i] = 0; }
  }
  const steps = N >= 1024 ? 4 : 12;
  const t = { route: 0, accumulate: 0, erode: 0, diffuse: 0 };
  let flooded = 0;
  const time = (k, f) => { const t0 = performance.now(); f(); t[k] += performance.now() - t0; };
  for (let k = 0; k < steps; k++) {
    time('route', () => s.route(z, outlet, fixed));
    if (s.flooded) flooded++;
    time('accumulate', () => s.accumulate());
    time('erode', () => s.erode(z, { dt: 1e5, K: 2e-6, m: 0.5, n: 1, U: 1e-3, fixed }));
    time('diffuse', () => s.diffuse(z, { dt: 1e5, D: 0.5, fixed }));
  }
  // the fast routing path alone, on the surface just after erosion (no depressions)
  s.route(z, outlet, fixed); s.accumulate(); s.erode(z, { dt: 1e5, K: 2e-6, U: 1e-3, fixed });
  const t0 = performance.now();
  const ok = s.routeDescent(z, outlet, fixed);
  const fast = performance.now() - t0;
  const total = Object.values(t).reduce((a, b) => a + b) / steps;
  console.log(`level 1, ${N}x${N}: ${total.toFixed(1)} ms/step (` +
    Object.entries(t).map(([k, v]) => `${k} ${(v / steps).toFixed(1)}`).join(', ') +
    `; flood used ${flooded}/${steps}); fast routing alone ${fast.toFixed(1)} ms${ok ? '' : ' (incomplete)'}`);
}
