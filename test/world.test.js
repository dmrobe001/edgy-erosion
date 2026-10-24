import { test } from 'node:test';
import assert from 'node:assert/strict';
import { philox4x32, worldKey, deriveKey, uniform, normal, fieldHash } from '../src/world/hash.js';
import {
  faceToVector, vectorToFace, tileAt, tileKey, neighborOf, cornerTiles, edgeKey, cornerKey,
  parentOf, childrenOf, cellCenter, tilesPerSide,
} from '../src/world/cubesphere.js';
import { generatePlanet } from '../src/world/l0.js';
import { RasterSolver } from '../src/world/fastscape.js';

const hex = (a) => Array.from(a, (x) => x.toString(16).padStart(8, '0')).join(' ');

test('Philox4x32-10 reproduces the Random123 known-answer vectors', () => {
  assert.equal(hex(philox4x32(0, 0, 0, 0, 0, 0)), '6627e8d5 e169c58d bc57ac4c 9b00dbd8');
  const f = 0xffffffff;
  assert.equal(hex(philox4x32(f, f, f, f, f, f)), '408f276d 41c83b0e a20bc7c6 6d5451fd');
  assert.equal(hex(philox4x32(0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344, 0xa4093822, 0x299f31d0)), 'd16cfe09 94fdcceb 5001e420 24126ea1');
});

test('keys depend on every part and its position; draws are uniform and repeatable', () => {
  const w = worldKey(42);
  assert.deepEqual(deriveKey(w, 2, 3, 5, 7), deriveKey(worldKey(42), 2, 3, 5, 7));
  const variants = [[2, 3, 5, 7], [2, 3, 7, 5], [2, 3, 5], [2, 3, 5, 7, 0], [-2, 3, 5, 7]].map((p) => deriveKey(w, ...p).join(','));
  assert.equal(new Set(variants).size, variants.length);
  assert.notDeepEqual(deriveKey(worldKey(43), 2, 3, 5, 7), deriveKey(w, 2, 3, 5, 7));
  assert.notDeepEqual(worldKey('alpha'), worldKey('alphb'));
  // moments of 1e5 uniforms and normals
  const k = deriveKey(w, 1);
  let s = 0, s2 = 0, g = 0, g2 = 0;
  const N = 100000;
  for (let i = 0; i < N; i++) {
    const u = uniform(k, i);
    assert.ok(u >= 0 && u < 1);
    s += u; s2 += u * u;
    const z = normal(k, i);
    g += z; g2 += z * z;
  }
  assert.ok(Math.abs(s / N - 0.5) < 0.005 && Math.abs(s2 / N - 1 / 3) < 0.005);
  assert.ok(Math.abs(g / N) < 0.02 && Math.abs(g2 / N - 1) < 0.02);
  assert.equal(uniform(k, 12345, 2, 9), uniform(k, 12345, 2, 9));
});

test('cube-sphere: face coordinates round-trip, and cells are within 1.5x of each other in width', () => {
  for (let f = 0; f < 6; f++) {
    for (const [a, b] of [[0.01, 0.5], [0.5, 0.5], [0.99, 0.02], [0.3, 0.8]]) {
      const v = faceToVector(f, a, b);
      const r = vectorToFace(...v);
      assert.equal(r.face, f);
      assert.ok(Math.abs(r.a - a) < 1e-12 && Math.abs(r.b - b) < 1e-12);
    }
  }
  const ang = (p, q) => Math.acos(Math.min(1, p[0] * q[0] + p[1] * q[1] + p[2] * q[2]));
  const n = 64;
  const centre = ang(faceToVector(0, 0.5, 0.5), faceToVector(0, 0.5 + 1 / n, 0.5));
  const edge = ang(faceToVector(0, 0.5 / n, 0.5), faceToVector(0, 1.5 / n, 0.5));
  assert.ok(centre / edge < 1.5 && edge / centre < 1.5, `${centre} vs ${edge}`);
});

test('cube-sphere: neighbours are mutual, edges and corners have one identity from every side', () => {
  const level = 2;
  const N = tilesPerSide(level);
  const k = deriveKey(worldKey(1), 99);
  const samples = [];
  // every face's corner tiles (the cube's edges and corners) plus random tiles
  for (let f = 0; f < 6; f++) for (const [i, j] of [[0, 0], [N - 1, 0], [0, N - 1], [N - 1, N - 1], [0, 7], [N - 1, 9]]) samples.push({ level, face: f, i, j });
  for (let s = 0; s < 40; s++) samples.push({ level, face: Math.floor(6 * uniform(k, s, 0)), i: Math.floor(N * uniform(k, s, 1)), j: Math.floor(N * uniform(k, s, 2)) });
  for (const t of samples) {
    for (const side of ['W', 'E', 'S', 'N']) {
      const nb = neighborOf(t, side);
      const back = ['W', 'E', 'S', 'N'].map((sd) => tileKey(neighborOf(nb, sd)));
      assert.ok(back.includes(tileKey(t)), `${tileKey(t)} -> ${side} -> ${tileKey(nb)} does not lead back`);
      const backSide = ['W', 'E', 'S', 'N'][back.indexOf(tileKey(t))];
      assert.equal(edgeKey(nb, backSide), edgeKey(t, side));
    }
    for (const c of ['SW', 'SE', 'NW', 'NE']) {
      const tiles = cornerTiles(t, c);
      assert.ok(tiles.length === 4 || tiles.length === 3, `${tiles.length} tiles at a corner`);
      assert.ok(tiles.some((u) => tileKey(u) === tileKey(t)));
      // every tile sharing the corner names it the same way
      const key = cornerKey(t, c);
      for (const u of tiles) assert.ok(['SW', 'SE', 'NW', 'NE'].some((c2) => cornerKey(u, c2) === key), `${key} from ${tileKey(u)}`);
    }
  }
  // the 8 cube corners are where only three tiles meet
  let three = 0;
  for (let f = 0; f < 6; f++) for (const [c, i, j] of [['SW', 0, 0], ['SE', N - 1, 0], ['NW', 0, N - 1], ['NE', N - 1, N - 1]]) if (cornerTiles({ level, face: f, i, j }, c).length === 3) three++;
  assert.equal(three, 24); // 8 corners x 3 faces
});

test('cube-sphere: children tile their parent', () => {
  const t = { level: 2, face: 3, i: 17, j: 200 };
  const kids = childrenOf(t);
  assert.equal(kids.length, 32 * 32);
  for (const c of [kids[0], kids[33], kids[1023]]) {
    assert.deepEqual(parentOf(c), t);
    const v = cellCenter(c, 256, 128, 128);
    assert.deepEqual(tileAt(2, ...v), t);
  }
  assert.deepEqual(parentOf({ level: 1, face: 0, i: 3, j: 4 }), { level: 0, face: -1, i: 0, j: 0 });
});

const SMALL_L0 = { subdivisions: 12, duration: 2e7, dt: 2e6, upliftWindow: 4e6 };

test('level 0 is a pure function of the world seed (golden hashes)', () => {
  const a = generatePlanet(7, SMALL_L0);
  const b = generatePlanet(7, SMALL_L0);
  const c = generatePlanet(8, SMALL_L0);
  const hashes = (t) => Object.fromEntries(Object.entries(t.fields).map(([k, v]) => [k, fieldHash(v)]));
  assert.deepEqual(hashes(a), hashes(b));
  assert.notDeepEqual(hashes(a).elevation, hashes(c).elevation);
  assert.deepEqual(a.events, b.events);
  // Stored on the development machine (Node 22, V8). A change here is either a
  // deliberate, versioned model change (bump GENERATOR_VERSION and update these) or a bug.
  assert.equal(a.version, 'hg-0.1');
  assert.deepEqual(hashes(a), GOLDEN_L0);
});

const GOLDEN_L0 = {
  elevation: '4dba2976', material: '37106f2d', uplift: '61763e85', precip: 'b46f0d05',
  temperature: '22ad8421', regolith: 'f5956342', drainageArea: 'c73b4d90',
};

test('level 0: sea level is 0, materials are consistent with the crust, rivers run downhill to the sea', () => {
  const t = generatePlanet(3, SMALL_L0);
  const f = t.fields;
  const n = f.elevation.length;
  let land = 0;
  for (let i = 0; i < n; i++) {
    if (f.elevation[i] > 0) land++;
    assert.ok(f.material[i] <= 4);
    assert.ok(f.temperature[i] > -90 && f.temperature[i] < 35);
    assert.ok(f.drainageArea[i] >= t.mesh.area[i] * 0.999);
  }
  assert.ok(land > 0.1 * n && land < 0.6 * n, `land ${land / n}`);
  // ocean floor is oceanic material
  let deep = 0, deepOceanic = 0;
  for (let i = 0; i < n; i++) if (f.elevation[i] < -3000) { deep++; if (f.material[i] === 0) deepOceanic++; }
  assert.ok(deepOceanic > 0.95 * deep);
  for (const e of t.rivers.edges) {
    if (e.to < 0) continue;
    assert.ok(f.drainageArea[t.rivers.nodes[e.to].cell] > f.drainageArea[t.rivers.nodes[e.from].cell]);
  }
});

test('stream power: steady state recovers S = (U / K)^(1/n) A^(-m/n), the slope-area law', () => {
  // a 64 x 64 block rising at U, draining to its south edge
  const N = 64, dx = 500;
  const s = new RasterSolver(N, N, dx);
  const z = new Float64Array(N * N), fixed = new Uint8Array(N * N), outlet = new Uint8Array(N * N);
  const k = deriveKey(worldKey(5), 1);
  for (let i = 0; i < N * N; i++) {
    const c = i % N, r = (i - c) / N;
    z[i] = 1 + uniform(k, i) + 0.01 * r;
    if (r === 0) { outlet[i] = 1; z[i] = 0; }
  }
  const U = 1e-3, K = 1e-5, m = 0.5;
  for (let it = 0; it < 400; it++) {
    s.route(z, outlet);
    s.accumulate();
    s.erode(z, { dt: 2e4, K, m, n: 1, U });
  }
  // compare every channel cell (A > 20 cells) with the law
  s.route(z, outlet);
  const A = s.accumulate();
  let worst = 0, cnt = 0;
  for (let i = 0; i < N * N; i++) {
    if (outlet[i] || A[i] < 20 * dx * dx) continue;
    const r = s.receiver[i];
    const S = (z[i] - z[r]) / s.recDist[i];
    worst = Math.max(worst, Math.abs(S / ((U / K) * A[i] ** -m) - 1));
    cnt++;
  }
  assert.ok(cnt > 100);
  assert.ok(worst < 1e-3, `worst relative misfit ${worst}`);
});

test('routing: the fast path and the priority flood agree on a surface without depressions', () => {
  const N = 48;
  const s = new RasterSolver(N, N, 100);
  const z = new Float64Array(N * N), outlet = new Uint8Array(N * N);
  const k = deriveKey(worldKey(9), 2);
  for (let i = 0; i < N * N; i++) {
    const c = i % N, r = (i - c) / N;
    z[i] = r + 0.5 * Math.sin(c / 3) + 0.01 * uniform(k, i) + 0.2 * c; // tilted: every cell has a lower neighbour
    if (r === 0) outlet[i] = 1;
  }
  assert.ok(s.routeDescent(z, outlet));
  const fast = Int32Array.from(s.receiver);
  const areaFast = Float64Array.from(s.accumulate());
  s.routeFlood(z, outlet);
  assert.deepEqual(Array.from(s.receiver), Array.from(fast));
  assert.deepEqual(Array.from(s.accumulate()), Array.from(areaFast));
  // and a depression makes the fast path decline, while the flood still drains everything
  z[20 * N + 20] = -5;
  assert.equal(s.routeDescent(z, outlet), false);
  s.route(z, outlet);
  assert.ok(s.flooded);
  assert.equal(s.orderCount, N * N);
});
