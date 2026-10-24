import { test } from 'node:test';
import assert from 'node:assert/strict';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU, regionSeams } from '../src/world/region.js';
import { tileAt, tileKey, neighborOf, tilesPerSide } from '../src/world/cubesphere.js';
import { tileBorder, SIDES } from '../src/world/seams.js';
import { fieldHash } from '../src/world/hash.js';

const layout = { rootTiles: 32, branching: 32, tileCells: 40 };
const opts = { iterations: 60 };
const planet = generatePlanet(21, { subdivisions: 16, duration: 1e8 });

/** A land tile, its east, north and north-east neighbours, and a pair across a cube edge. */
function tiles() {
  const e = planet.fields.elevation, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > e[c]) c = i;
  const A = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
  const B = neighborOf(A, 'E', layout), C = neighborOf(A, 'N', layout), D = neighborOf(B, 'N', layout);
  const N = tilesPerSide(1, layout);
  const X = { level: 1, face: 2, i: N - 1, j: Math.floor(N / 2) };
  const Y = neighborOf(X, 'E', layout);
  assert.notEqual(X.face, Y.face, 'X and Y straddle a cube edge');
  return { A, B, C, D, X, Y };
}

async function generate(order, fresh) {
  const ctx = regionSeams(planet, opts, layout, { fresh });
  const out = new Map();
  for (const t of order) out.set(tileKey(t), await runRegionCPU(regionInput(planet, t, opts, layout, ctx)));
  return out;
}

const hashes = (r) => Object.fromEntries(Object.entries(r.fields).map(([k, v]) => [k, fieldHash(v)]));

/** Border values of side `side` of a tile result, as an array in the tile's node order. */
function borderOf(r, side) {
  const m = r.nodes, e = r.fields.elevation, out = [];
  for (let t = 0; t < m; t++) {
    const [ci, cj] = side === 'W' ? [0, t] : side === 'E' ? [m - 1, t] : side === 'S' ? [t, 0] : [t, m - 1];
    out.push(e[cj * m + ci]);
  }
  return out;
}

test('seams: generation order and cache state change nothing, and shared borders are bitwise identical', async () => {
  const T = tiles();
  const list = Object.values(T);
  const first = await generate(list, false);
  const second = await generate([...list].reverse(), true);
  for (const t of list) assert.deepEqual(hashes(second.get(tileKey(t))), hashes(first.get(tileKey(t))), `tile ${tileKey(t)} depends on order`);
  // every pair of neighbours in the set agrees on its shared border, node for node
  let pairs = 0;
  for (const a of list) {
    for (const side of SIDES) {
      const b = neighborOf(a, side, layout);
      const rb = first.get(tileKey(b));
      if (!rb) continue;
      const back = SIDES.find((s) => tileKey(neighborOf(b, s, layout)) === tileKey(a));
      const va = borderOf(first.get(tileKey(a)), side), vb = borderOf(rb, back);
      const same = va.every((v, k) => v === vb[k]) || va.every((v, k) => v === vb[vb.length - 1 - k]);
      assert.ok(same, `${tileKey(a)} ${side} vs ${tileKey(b)} ${back}`);
      pairs++;
    }
  }
  assert.ok(pairs >= 10, `${pairs} shared borders checked`); // 2x2 block: 8 (both ways) + the cube-edge pair: 2
});

test('seams: a crossing seen from both sides is the same river, flowing out of one tile and into the other', () => {
  const ctx = regionSeams(planet, opts, layout, { fresh: true });
  const { A, B, C, D, X, Y } = tiles();
  let matched = 0;
  for (const a of [A, B, C, D, X, Y]) {
    for (const side of SIDES) {
      const b = neighborOf(a, side, layout);
      const back = SIDES.find((s) => tileKey(neighborOf(b, s, layout)) === tileKey(a));
      const xa = tileBorder(ctx, a)[side].crossings, xb = tileBorder(ctx, b)[back].crossings;
      assert.equal(xa.length, xb.length);
      const n = layout.tileCells;
      for (const x of xa) {
        const y = xb.find((q) => q.node === x.node || q.node === n - x.node);
        assert.ok(y, `crossing at node ${x.node} of ${tileKey(a)} ${side} not found from ${tileKey(b)}`);
        assert.equal(y.inflow, !x.inflow);
        assert.equal(y.area, x.area);
        assert.equal(y.z, x.z);
        matched++;
      }
    }
  }
  assert.ok(matched > 0, 'some crossings to compare');
});
