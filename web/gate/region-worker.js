// Level-1 tiles on WebGPU, for the region gate (scripts/browser-gate.js).
import { generatePlanet } from '../../src/world/l0.js';
import { regionInput, runRegion, cpuSolver } from '../../src/world/region.js';
import { gpuSolver, getDevice } from '../../src/world/gpu.js';
import { tileAt, tileKey } from '../../src/world/cubesphere.js';

/** Land tiles at quantiles of the planet's land elevation (deterministic). */
function pickTiles(planet, count) {
  const e = planet.fields.elevation, P = planet.mesh.pos;
  const land = [];
  for (let i = 0; i < e.length; i++) if (e[i] > 200) land.push(i);
  land.sort((a, b) => e[a] - e[b] || a - b);
  const seen = new Set(), out = [];
  for (let k = 0; k < count; k++) {
    const i = land[Math.floor(((k + 0.5) / count) * land.length)];
    const t = tileAt(1, P[3 * i], P[3 * i + 1], P[3 * i + 2]);
    if (!seen.has(tileKey(t))) { seen.add(tileKey(t)); out.push({ addr: t, planetElevation: e[i] }); }
  }
  return out;
}

const summary = (r) => ({
  hack: r.stats.hack.exponent, hackR2: r.stats.hack.r2,
  concavity: r.stats.slopeArea.concavity, concavityR2: r.stats.slopeArea.r2,
  consistency: r.stats.consistency, riverNodes: r.rivers.nodes.length, maxOrder: Math.max(0, ...r.rivers.edges.map((e) => e.order)),
});

self.onmessage = async (e) => {
  try {
    const { seed, tiles: count, compareCPU } = e.data;
    const device = await getDevice();
    if (!device) { self.postMessage({ error: 'WebGPU is not available in this worker' }); return; }
    const info = device.adapterInfo || {};
    const t0 = performance.now();
    const planet = generatePlanet(seed);
    const planetSeconds = (performance.now() - t0) / 1000;
    const results = [];
    // one warm-up tile (shader compilation), not timed
    await runRegion(regionInput(planet, pickTiles(planet, 1)[0].addr), gpuSolver(device, regionInput(planet, pickTiles(planet, 1)[0].addr)));
    for (const { addr, planetElevation } of pickTiles(planet, count)) {
      const t1 = performance.now();
      const input = regionInput(planet, addr);
      const solver = gpuSolver(device, input);
      const r = await runRegion(input, solver);
      solver.destroy();
      const seconds = (performance.now() - t1) / 1000;
      const row = { tile: tileKey(addr), planetElevation, seconds, ...summary(r) };
      if (compareCPU && results.length === 0) {
        const t2 = performance.now();
        const input2 = regionInput(planet, addr);
        const c = await runRegion(input2, cpuSolver(input2));
        row.cpuSeconds = (performance.now() - t2) / 1000;
        row.cpu = summary(c);
        let ss = 0, sr = 0;
        const a = r.fields.elevation, b = c.fields.elevation;
        for (let k = 0; k < a.length; k++) { ss += (a[k] - b[k]) ** 2; sr += b[k] ** 2; }
        row.gpuCpuRms = Math.sqrt(ss / a.length);
        let lo = Infinity, hi = -Infinity;
        for (const v of b) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
        row.relief = hi - lo;
      }
      results.push(row);
    }
    self.postMessage({ adapter: { vendor: info.vendor, architecture: info.architecture }, planetSeconds, results });
  } catch (err) {
    self.postMessage({ error: String(err && err.stack || err) });
  }
};
