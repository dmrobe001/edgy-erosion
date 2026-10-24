// Level 0: the planet, built once per world.
//
// The kinematic plate model of src/planet.js, run without river routing (land is worn
// down by a cheap elevation-dependent denudation instead), for `duration` years. The
// result is packed into an immutable tile record: raster fields on the icosahedral mesh,
// a coarse river graph, and the dated event log. Every finer level reads this record and
// never the Planet object itself, so the record is the whole contract.

import { Planet, MATERIAL_NAMES } from '../planet.js';
import { MinHeap } from '../flow.js';
import { worldKey, deriveKey } from './hash.js';

export const GENERATOR_VERSION = 'hg-0.1';

/**
 * Level-0 settings. The defaults run in about 11 s in Node on one core (measured on the
 * development container): 23 042 cells about 160 km apart, 250 steps of 2 Myr. The
 * design's 30-50 km cells would need about 60x the work; see docs/PLAN.md.
 */
export const L0_DEFAULTS = {
  subdivisions: 48,
  duration: 5e8,                 // yr of plate history
  dt: 2e6,                       // yr per step
  upliftWindow: 1e7,             // yr: uplift is averaged over the last part of the history
  riverArea: 2e11,               // m^2: drainage area at which a level-0 river is recorded
  lapseRate: 0.0065,             // K per m
  planet: {
    surfaceProcesses: false,
    crust: { maxFlowSubsteps: 8, flowThreshold: 55000, maxThickness: 70000 },
    continents: { fraction: 0.4 },
    plates: { reorganizeEvery: 1e8 },
    denudation: { timescale: 6e7 },
  },
};

/** Mean annual sea-level temperature (deg C) at latitude lat (rad): a smooth caricature of Earth's. */
export function seaLevelTemperature(lat) {
  const s = Math.sin(lat);
  return 27 - 52 * s * s;
}

/**
 * Build the planet tile.
 * @param {number|string|bigint} worldSeed
 * @param {object} [options] overrides of L0_DEFAULTS (planet: Planet parameter overrides)
 * @param {(fraction:number) => void} [onProgress]
 */
export function generatePlanet(worldSeed, options = {}, onProgress) {
  const o = { ...L0_DEFAULTS, ...options, planet: mergePlanet(L0_DEFAULTS.planet, options.planet || {}) };
  const root = worldKey(worldSeed);
  const key = deriveKey(root, 0, -1, 0, 0);
  const planet = new Planet({ ...o.planet, seed: key[0] % 2147483647 || 1, subdivisions: o.subdivisions, time: { dt: o.dt } });
  const steps = Math.round(o.duration / o.dt);
  const windowSteps = Math.max(1, Math.min(steps, Math.round(o.upliftWindow / o.dt)));
  const n = planet.mesh.size;
  const upliftSum = new Float64Array(n);
  for (let k = 0; k < steps; k++) {
    planet.step();
    if (k >= steps - windowSteps) for (let i = 0; i < n; i++) upliftSum[i] += planet.upliftRate[i];
    if (onProgress && (k % 10 === 9 || k === steps - 1)) onProgress((k + 1) / steps);
  }
  return packPlanet(planet, { o, key, worldSeed, upliftSum, windowSteps });
}

function mergePlanet(a, b) {
  const out = { ...a };
  for (const [k, v] of Object.entries(b)) out[k] = v && typeof v === 'object' && !Array.isArray(v) ? { ...(a[k] || {}), ...v } : v;
  return out;
}

function packPlanet(planet, { o, key, worldSeed, upliftSum, windowSteps }) {
  const mesh = planet.mesh;
  const n = mesh.size;
  const sea = planet.seaLevel;
  const elevation = new Float32Array(n);
  const uplift = new Float32Array(n);
  const precip = new Float32Array(n);
  const temperature = new Float32Array(n);
  const regolith = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    // elevations are stored relative to sea level: sea level is 0 at every finer level
    elevation[i] = planet.sim.z[i] - sea;
    uplift[i] = upliftSum[i] / windowSteps;
    precip[i] = 1000 * planet.sim.precipWeight[i];
    temperature[i] = seaLevelTemperature(mesh.lat[i]) - o.lapseRate * Math.max(0, elevation[i]);
  }
  const material = planet.material();
  const { area, receiver, order } = drainage(mesh, elevation);
  for (let i = 0; i < n; i++) {
    // [conjecture] a placeholder until a weathering model runs at level 1: deeper where
    // warm and wet and low, nothing under the sea or on the highest ground
    if (elevation[i] <= 0) continue;
    const warm = Math.max(0, Math.min(1, (temperature[i] + 5) / 30));
    const wet = Math.min(1, precip[i] / 1500);
    regolith[i] = 10 * warm * wet * Math.exp(-elevation[i] / 2000);
  }
  const T = planet.time;
  const events = planet.eventLog().map((e) => ({
    kind: e.kind,
    time: T - e.end,             // years before present that it ended (0: still active)
    duration: e.end - e.start,
    region: { type: 'cap', center: e.center, radius: e.radius },
    params: { plateA: e.plates[0], plateB: e.plates[1], cells: e.cells },
  }));
  return {
    version: GENERATOR_VERSION,
    worldSeed: String(worldSeed),
    addr: { level: 0, face: -1, i: 0, j: 0 },
    key,
    mesh,
    seaLevel: 0,
    fields: { elevation, material, uplift, precip, temperature, regolith, drainageArea: area },
    plates: { owner: Int32Array.from(planet.owner), boundary: Uint8Array.from(planet.boundary) },
    rivers: riverGraph(mesh, elevation, area, receiver, order, o.riverArea, precip),
    events,
    materialNames: MATERIAL_NAMES,
    stats: { ...planet.stats(), duration: T, cells: n, spacing: mesh.spacing },
  };
}

/**
 * Drainage on the sphere: priority flood from the sea (cells at or below 0), so every
 * land cell drains, through any depressions, to the sea. Each land cell's receiver is
 * the cell it was reached from. Area includes the cell itself (m^2).
 */
function drainage(mesh, z) {
  const n = mesh.size;
  const receiver = new Int32Array(n).fill(-1);
  const level = new Float64Array(n);
  const done = new Uint8Array(n);
  const heap = new MinHeap(7 * n);
  const order = [];
  for (let i = 0; i < n; i++) if (z[i] <= 0) { heap.push(i, z[i]); level[i] = z[i]; done[i] = 1; }
  while (heap.n) {
    const i = heap.pop();
    order.push(i);
    for (let k = 0; k < mesh.valence[i]; k++) {
      const j = mesh.neighbors[6 * i + k];
      if (j < 0 || done[j]) continue;
      done[j] = 1;
      receiver[j] = i;
      level[j] = Math.max(z[j], level[i]);
      heap.push(j, level[j]);
    }
  }
  const area = new Float32Array(n);
  const acc = Float64Array.from(mesh.area);
  for (let k = order.length - 1; k >= 0; k--) {
    const i = order[k];
    if (receiver[i] >= 0) acc[receiver[i]] += acc[i];
  }
  for (let i = 0; i < n; i++) area[i] = z[i] > 0 ? acc[i] : mesh.area[i];
  return { area, receiver, order };
}

/** Rivers draining more than `minArea`, as a graph of cell centres with Strahler orders. */
function riverGraph(mesh, z, area, receiver, order, minArea, precip) {
  const n = mesh.size;
  const id = new Int32Array(n).fill(-1);
  const nodes = [];
  for (const i of order) {
    if (!(z[i] > 0 && area[i] >= minArea)) continue;
    id[i] = nodes.length;
    nodes.push({ cell: i, x: mesh.pos[3 * i], y: mesh.pos[3 * i + 1], z: mesh.pos[3 * i + 2], elevation: z[i] });
  }
  // Strahler order, upstream to downstream
  const strahler = new Int32Array(n);
  const top = new Int32Array(n), topCount = new Int32Array(n);
  for (let k = order.length - 1; k >= 0; k--) {
    const i = order[k];
    if (id[i] < 0) continue;
    strahler[i] = top[i] === 0 ? 1 : topCount[i] > 1 ? top[i] + 1 : top[i];
    const r = receiver[i];
    if (r < 0) continue;
    if (strahler[i] > top[r]) { top[r] = strahler[i]; topCount[r] = 1; } else if (strahler[i] === top[r]) topCount[r]++;
  }
  const edges = [];
  for (const nd of nodes) {
    const i = nd.cell, r = receiver[i];
    if (r < 0) continue;
    // a river ends in the first sea cell, which is not a node: record it as the mouth
    const to = id[r] >= 0 ? id[r] : -1;
    edges.push({ from: id[i], to, toCell: r, discharge: (area[i] * precip[i] * 0.5) / 1000 / 3.15576e7, order: strahler[i] });
  }
  return { nodes, edges, crossings: [] };
}

