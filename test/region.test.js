import { test } from 'node:test';
import assert from 'node:assert/strict';
import { JacobiSolver, FLAG_OUTLET } from '../src/world/jacobi.js';
import { RasterSolver } from '../src/world/fastscape.js';
import { upstreamOrder, mainStreamLength, extractRiverGraph, drainageStatistics } from '../src/world/rivers.js';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU } from '../src/world/region.js';
import { tileAt } from '../src/world/cubesphere.js';
import { uniform, worldKey, deriveKey } from '../src/world/hash.js';

function block(N, dx, seed) {
  const n = N * N, k = deriveKey(worldKey(seed), 1);
  const z = new Float32Array(n), U = new Float32Array(n).fill(1e-3), K = new Float32Array(n).fill(1e-5), P = new Float32Array(n).fill(1), flags = new Uint32Array(n);
  for (let i = 0; i < n; i++) {
    const r = Math.floor(i / N);
    z[i] = 1 + uniform(k, i) + 0.01 * r;
    if (r === 0) { flags[i] = FLAG_OUTLET; z[i] = 0; }
  }
  return { z, U, K, P, flags };
}

test('parallel stream power: steady state is exactly the implicit one, K A^m S = U', () => {
  const N = 48, dx = 500;
  const s = new JacobiSolver(N, N, dx, block(N, dx, 5));
  s.run(1500, { dt: 5e3, m: 0.5 });
  let worst = 0, cnt = 0;
  for (let i = 0; i < N * N; i++) {
    const r = s.receiverIndex(i);
    if (r < 0 || s.A[i] < 20 * dx * dx) continue;
    const S = (s.z[i] - s.z[r]) / s.receiverDistance(i);
    worst = Math.max(worst, Math.abs(S / (100 * s.A[i] ** -0.5) - 1));
    cnt++;
  }
  assert.ok(cnt > 100);
  assert.ok(worst < 1e-2, `worst relative misfit ${worst}`);
});

test('parallel and implicit solvers reach the same steady landscape', () => {
  const N = 40, dx = 500;
  const inp = block(N, dx, 6);
  const s = new JacobiSolver(N, N, dx, inp);
  s.run(1500, { dt: 5e3, m: 0.5 });
  const r = new RasterSolver(N, N, dx);
  const z = Float64Array.from(inp.z), outlet = Uint8Array.from(inp.flags, (v) => v & FLAG_OUTLET);
  for (let it = 0; it < 750; it++) { r.route(z, outlet); r.accumulate(); r.erode(z, { dt: 1e4, K: 1e-5, m: 0.5, U: 1e-3 }); }
  // the networks differ in detail (each picks its own receivers on ties), so compare relief statistics
  const stats = (a) => { let s = 0, mx = 0; for (const v of a) { s += v; mx = Math.max(mx, v); } return [s / a.length, mx]; };
  const [mj, xj] = stats(s.z), [mi, xi] = stats(z);
  assert.ok(Math.abs(mj / mi - 1) < 0.1 && Math.abs(xj / xi - 1) < 0.1, `mean ${mj} vs ${mi}, max ${xj} vs ${xi}`);
});

test('river graph: heads, confluences and mouths, Strahler orders, and Hack lengths on a hand-made tree', () => {
  // a 5 x 3 grid; two streams join at cell 7 and leave through cell 12 (bottom middle)
  //  0  1  2  3  4      rows top to bottom; 0 -> 6 -> 7, 4 -> 8 -> 7, 7 -> 12
  //  5  6  7  8  9
  // 10 11 12 13 14
  const receiver = new Int32Array(15).fill(-1);
  receiver[0] = 6; receiver[6] = 7; receiver[4] = 8; receiver[8] = 7; receiver[7] = 12;
  const dist = new Float32Array(15).fill(1);
  const area = new Float32Array(15).fill(1);
  for (const i of upstreamOrder(receiver)) if (receiver[i] >= 0) area[receiver[i]] += area[i];
  const z = new Float32Array(15).map((_, i) => 10 - Math.floor(i / 5));
  const L = mainStreamLength(receiver, dist);
  assert.equal(L[7], 2); assert.equal(L[12], 3);
  const g = extractRiverGraph({ receiver, dist, area, z, nx: 5 }, { channelArea: 1 });
  const kinds = g.nodes.map((n) => `${n.cell}:${n.kind}`).sort();
  assert.deepEqual(kinds, ['0:head', '12:mouth', '4:head', '7:confluence'].sort());
  const toMouth = g.edges.find((e) => g.nodes[e.from].cell === 7);
  assert.equal(toMouth.order, 2);
  assert.deepEqual(toMouth.path, [7, 12]);
  assert.equal(g.edges.length, 3);
});

test('a region tile from the planet drains everywhere, keeps the parent block means, and obeys Hack and slope-area laws', async () => {
  const planet = generatePlanet(11, { subdivisions: 16, duration: 1e8 });
  const e = planet.fields.elevation, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > e[c] && e[i] < 2500) c = i;
  const layout = { rootTiles: 32, branching: 32, tileCells: 96 };
  const addr = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
  const input = regionInput(planet, addr, { halo: 12, iterations: 160 }, layout);
  const r = await runRegionCPU(input);
  const d = r.drainage;
  let undrained = 0;
  for (let k = 0; k < d.receiver.length; k++) if (d.receiver[k] < 0 && !(input.flags[k] & FLAG_OUTLET)) undrained++;
  assert.equal(undrained, 0);
  assert.ok(r.stats.consistency < 1, `block means off by ${r.stats.consistency} m`);
  assert.ok(r.stats.hack.exponent > 0.45 && r.stats.hack.exponent < 0.7, `Hack exponent ${r.stats.hack.exponent}`);
  assert.ok(r.stats.slopeArea.concavity > 0.3 && r.stats.slopeArea.concavity < 0.7, `concavity ${r.stats.slopeArea.concavity}`);
  assert.ok(r.rivers.edges.length > 20);
});
