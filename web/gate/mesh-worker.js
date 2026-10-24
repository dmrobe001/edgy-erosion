// Level-1 tiles on the irregular mesh, on WebGPU, for the mesh gate (scripts/browser-gate.js --only=mesh).
import { generatePlanet } from '../../src/world/l0.js';
import { meshRegionContext, meshRegionTile } from '../../src/world/regionmesh.js';
import { gpuGraphSolver } from '../../src/world/gpugraph.js';
import { getDevice } from '../../src/world/gpu.js';
import { tileAt, tileKey, WORLD_LAYOUT } from '../../src/world/cubesphere.js';

function landTile(planet, layout) {
  const e = planet.fields.elevation, P = planet.mesh.pos, M = planet.mesh;
  // land cells with land all around (an island cell alone smooths into the sea)
  const inland = (i) => { for (let k = 0; k < M.valence[i]; k++) { const j = M.neighbors[6 * i + k]; if (j >= 0 && e[j] <= 200) return false; } return true; };
  const land = [];
  for (let i = 0; i < e.length; i++) if (e[i] > 200 && inland(i)) land.push(i);
  land.sort((a, b) => e[a] - e[b] || a - b);
  const i = land[Math.floor(0.6 * land.length)];
  return { addr: tileAt(1, P[3 * i], P[3 * i + 1], P[3 * i + 2], layout), planetElevation: e[i] };
}

const summary = (r) => ({ hack: r.stats.hack.exponent, hackR2: r.stats.hack.r2, concavity: r.stats.slopeArea.concavity, concavityR2: r.stats.slopeArea.r2, consistency: r.stats.consistency });

self.onmessage = async (e) => {
  try {
    const { seed, full } = e.data;
    const device = await getDevice();
    if (!device) { self.postMessage({ error: 'WebGPU is not available in this worker' }); return; }
    const info = device.adapterInfo || {};
    const planet = generatePlanet(seed);
    const gpu = (win) => gpuGraphSolver(device, win);
    // agreement with the CPU reference, on a small layout (the CPU is slow)
    const small = { ...WORLD_LAYOUT, tileCells: 64 };
    const s = landTile(planet, small);
    const opts = { iterations: 200 };
    let t = performance.now();
    const g = await meshRegionTile(meshRegionContext(planet, opts, small), s.addr, gpu);
    const gpuSmall = (performance.now() - t) / 1000;
    t = performance.now();
    const c = await meshRegionTile(meshRegionContext(planet, opts, small), s.addr);
    const cpuSmall = (performance.now() - t) / 1000;
    let ss = 0, lo = Infinity, hi = -Infinity;
    const a = g.fields.elevation, b = c.fields.elevation;
    for (let k = 0; k < a.length; k++) { ss += (a[k] - b[k]) ** 2; lo = Math.min(lo, b[k]); hi = Math.max(hi, b[k]); }
    const out = {
      adapter: { vendor: info.vendor, architecture: info.architecture },
      small: { tile: tileKey(s.addr), cells: small.tileCells, gpuSeconds: gpuSmall, cpuSeconds: cpuSmall, rms: Math.sqrt(ss / a.length), relief: hi - lo, gpu: summary(g), cpu: summary(c) },
    };
    if (full) {
      // a full-size tile with its strips and corners (shader compilation already done)
      const f = landTile(planet, WORLD_LAYOUT);
      t = performance.now();
      const r = await meshRegionTile(meshRegionContext(planet, {}), f.addr, gpu);
      out.full = { tile: tileKey(f.addr), planetElevation: f.planetElevation, seconds: (performance.now() - t) / 1000, ...summary(r) };
    }
    self.postMessage(out);
  } catch (err) {
    self.postMessage({ error: String((err && err.stack) || err) });
  }
};
