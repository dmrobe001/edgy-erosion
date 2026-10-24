import { test } from 'node:test';
import assert from 'node:assert/strict';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU, regionSeams } from '../src/world/region.js';
import { landscapeSeams, landscapeInput, runLandscape, landscapeSummary } from '../src/world/landscape.js';
import { tileAt, tileKey, neighborOf, nodePosition } from '../src/world/cubesphere.js';
import { seamContext, SIDES } from '../src/world/seams.js';
import { RasterSolver } from '../src/world/fastscape.js';
import { fieldHash } from '../src/world/hash.js';

const small = { rootTiles: 32, branching: 4, tileCells: 32 };
const planet = generatePlanet(31, { subdivisions: 16, duration: 1e8 });
const hashes = (r) => Object.fromEntries(Object.entries(r.fields).map(([k, v]) => [k, fieldHash(v)]));
function borderOf(r, side) {
  const m = r.nodes, e = r.fields.elevation, out = [];
  for (let t = 0; t < m; t++) {
    const [ci, cj] = side === 'W' ? [0, t] : side === 'E' ? [m - 1, t] : side === 'S' ? [t, 0] : [t, m - 1];
    out.push(e[cj * m + ci]);
  }
  return out;
}

async function twoParents() {
  const e = planet.fields.elevation, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > e[c]) c = i;
  const A1 = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], small);
  const B1 = neighborOf(A1, 'E', small);
  const s1 = regionSeams(planet, { iterations: 40 }, small);
  const rA = await runRegionCPU(regionInput(planet, A1, { iterations: 40 }, small, s1));
  const rB = await runRegionCPU(regionInput(planet, B1, { iterations: 40 }, small, s1));
  return { A1, B1, s1, rA, rB };
}

test('level-2 seams on a level-1 seam: tiles under different parents agree bitwise, in any order', async () => {
  const { A1, s1, rA, rB } = await twoParents();
  const b = small.branching;
  const T = { level: 2, face: A1.face, i: A1.i * b + b - 1, j: A1.j * b + 1 }; // on A1's east edge
  const Ut = neighborOf(T, 'E', small);                                        // under B1
  const V = neighborOf(T, 'N', small);                                         // under A1
  const opts = { steps: 6 };
  const run = (parent, ctx, t) => runLandscape(landscapeInput(parent, ctx, t, planet.events, opts, small));
  const ctxA = landscapeSeams(rA, s1, planet.worldSeed, opts, small), ctxB = landscapeSeams(rB, s1, planet.worldSeed, opts, small);
  const first = { T: run(rA, ctxA, T), U: run(rB, ctxB, Ut), V: run(rA, ctxA, V) };
  const ctxA2 = landscapeSeams(rA, s1, planet.worldSeed, opts, small), ctxB2 = landscapeSeams(rB, s1, planet.worldSeed, opts, small);
  const second = { V: run(rA, ctxA2, V), U: run(rB, ctxB2, Ut), T: run(rA, ctxA2, T) };
  for (const k of ['T', 'U', 'V']) assert.deepEqual(hashes(second[k]), hashes(first[k]), `${k} depends on order`);
  assert.deepEqual(borderOf(first.T, 'E'), borderOf(first.U, 'W'), 'border across the level-1 seam');
  assert.deepEqual(borderOf(first.T, 'N'), borderOf(first.V, 'S'), 'border inside one parent');
  // and the level-1 seam's crossings reach level 2 at the same place, from both sides
  const xT = first.T.border.E.crossings, xU = first.U.border.W.crossings;
  assert.equal(xT.length, xU.length);
  for (const x of xT) assert.ok(xU.some((y) => y.node === x.node && y.inflow === !x.inflow && y.area === x.area));
});

test('consistency at every level: block means within tolerance of the parent, planet to level 2', async () => {
  const layout = { rootTiles: 32, branching: 8, tileCells: 64 };
  const e = planet.fields.elevation, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > e[c]) c = i;
  const A1 = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
  const s1 = regionSeams(planet, { iterations: 100 }, layout);
  const rA = await runRegionCPU(regionInput(planet, A1, { iterations: 100 }, layout, s1));
  assert.ok(rA.stats.consistency < 2, `level 1 off its parent by ${rA.stats.consistency} m`);
  const T = { level: 2, face: A1.face, i: A1.i * 8 + 3, j: A1.j * 8 + 4 };
  const ctx = landscapeSeams(rA, s1, planet.worldSeed, {}, layout);
  const r2 = runLandscape(landscapeInput(rA, ctx, T, planet.events, { steps: 8 }, layout));
  assert.ok(r2.stats.consistency < 1, `level 2 off its parent by ${r2.stats.consistency} m`);
});

// --- variety: one parent class, different histories, against noise alone ---------------

/** A synthetic level-1 parent: the same upland everywhere, with chosen rock, climate and uplift. */
function syntheticParent({ material, temperature, uplift, precip }, layout) {
  const n1 = layout.tileCells, m = n1 + 1;
  const addr = { level: 1, face: 0, i: 5, j: 7 };
  const R = 6.371e6;
  const hill = (v) => 1500 + 500 * Math.sin(40 * v[1]) * Math.cos(37 * v[2]) + 250 * Math.sin(140 * v[1] + 90 * v[2]);
  const elev = new Float32Array(m * m), pos = [];
  for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) { const v = nodePosition(addr, n1, i, j, undefined, layout); pos.push(v); elev[j * m + i] = hill(v); }
  const a = pos[Math.floor(m / 2) * m + Math.floor(m / 2)], b2 = pos[Math.floor(m / 2) * m + Math.floor(m / 2) + 1];
  const dx = R * Math.acos(a[0] * b2[0] + a[1] * b2[1] + a[2] * b2[2]);
  const rs = new RasterSolver(m, m, dx);
  const outlet = new Uint8Array(m * m);
  for (let i = 0; i < m; i++) outlet[i] = 1; // drains south
  rs.route(Float64Array.from(elev), outlet);
  const area = Float32Array.from(rs.accumulate());
  const receiver = Int32Array.from(rs.receiver, (r, c) => (r === c ? -1 : r));
  const fill = (v) => new Float32Array(m * m).fill(v);
  const parent = {
    addr, nodes: m, dx,
    fields: { elevation: elev, uplift: fill(uplift), precip: fill(precip), material: new Uint8Array(m * m).fill(material), temperature: Float32Array.from(elev, (z) => temperature - 0.0065 * z), drainageArea: area },
    drainage: { receiver, area },
  };
  const parentSeams = seamContext({
    level: 1, worldSeed: 'variety', layout, radius: R, linkSpacing: 1e9, links: [],
    elevation: hill, base: hill, notchDepth: () => 0,
  });
  return { parent, parentSeams, addr };
}

test('variety: landscapes under one parent class differ by rock and climate far more than by noise', () => {
  const layout = { rootTiles: 32, branching: 32, tileCells: 96 };
  const histories = [
    { name: 'granite, temperate', material: 1, temperature: 12, uplift: 5e-4, precip: 900 },
    { name: 'schist, temperate', material: 2, temperature: 12, uplift: 5e-4, precip: 900 },
    { name: 'lava and tuff, wet', material: 3, temperature: 15, uplift: 5e-4, precip: 2000 },
    { name: 'sandstone and shale, dry', material: 4, temperature: 15, uplift: 2e-4, precip: 300 },
    { name: 'granite, glaciated', material: 1, temperature: 3, uplift: 5e-4, precip: 900 },
  ];
  const events = [{ kind: 'glaciation', time: 2e4, duration: 1e5, region: { type: 'global' }, params: { cooling: 6 } }];
  const run = (h, salt) => {
    const { parent, parentSeams, addr } = syntheticParent(h, layout);
    const T = { level: 2, face: addr.face, i: addr.i * 32 + 13, j: addr.j * 32 + 17 };
    const opts = { steps: 12, noiseSalt: salt };
    const ctx = landscapeSeams(parent, parentSeams, 'variety', opts, layout);
    return landscapeSummary(runLandscape(landscapeInput(parent, ctx, T, events, opts, layout)));
  };
  const byHistory = histories.map((h) => run(h, 0));
  const byNoise = [0, 1, 2, 3].map((salt) => run(histories[1], salt));
  const keys = Object.keys(byHistory[0]);
  const spread = (rows, k) => {
    const v = rows.map((r) => r[k]), mean = v.reduce((a, b) => a + b, 0) / v.length;
    return Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / (v.length - 1));
  };
  const table = keys.map((k) => ({ k, history: spread(byHistory, k), noise: spread(byNoise, k) }));
  const beyond = table.filter((r) => r.history > 2 * r.noise + 1e-9);
  // eslint-disable-next-line no-console
  if (process.env.VERBOSE) console.table(table);
  assert.ok(beyond.length >= Math.ceil(keys.length * 0.6), `only ${beyond.map((r) => r.k).join(', ')} vary beyond noise`);
});
