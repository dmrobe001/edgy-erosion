// Level 1 on the irregular mesh: simulated seams (corner patches, edge strips, tiles).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { generatePlanet } from '../src/world/l0.js';
import { tileAt, neighborOf, tileKey, childrenOf, tilesPerSide } from '../src/world/cubesphere.js';
import { SIDES, sideNode } from '../src/world/seams.js';
import { meshRegionContext, meshRegionTile } from '../src/world/regionmesh.js';
import { World } from '../src/world/pipeline.js';

const planet = generatePlanet(21, { subdivisions: 16, duration: 1e8 });
const layout = { rootTiles: 32, branching: 32, tileCells: 32 };
const opts = { iterations: 60 };

function landTile() {
  const F = planet.fields, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < F.elevation.length; i++) if (Math.abs(F.elevation[i] - 1500) < Math.abs(F.elevation[c] - 1500)) c = i;
  return tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
}

async function run(tiles) {
  const ctx = meshRegionContext(planet, opts, layout);
  const out = new Map();
  for (const t of tiles) out.set(tileKey(t), await meshRegionTile(ctx, t));
  return out;
}

function seamMismatches(out, tiles) {
  const n = layout.tileCells;
  let bad = 0, checked = 0;
  for (const t of tiles) for (const s of SIDES) {
    const nb = neighborOf(t, s, layout);
    if (!out.has(tileKey(nb))) continue;
    const back = SIDES.find((q) => tileKey(neighborOf(nb, q, layout)) === tileKey(t));
    const ea = out.get(tileKey(t)).fields.elevation, eb = out.get(tileKey(nb)).fields.elevation;
    const ba = out.get(tileKey(t)).borders[s], bb = out.get(tileKey(nb)).borders[back];
    // node order along the shared edge may run the other way on the other tile
    const rev = ba.z[0] !== bb.z[0];
    for (let u = 0; u <= n; u++) {
      const [i, j] = sideNode(s, u, n), [i2, j2] = sideNode(back, rev ? n - u : u, n);
      if (ea[j * (n + 1) + i] !== eb[j2 * (n + 1) + i2]) bad++;
      checked++;
    }
    // the same rivers cross, flowing out of one and into the other
    assert.equal(ba.crossings.length, bb.crossings.length);
    for (const x of ba.crossings) {
      const y = bb.crossings.find((c) => c.node === (rev ? n - x.node : x.node));
      assert.ok(y, `crossing at ${x.node}`);
      assert.equal(y.inflow, !x.inflow);
      assert.deepEqual(y.sed, x.sed);
    }
  }
  return { bad, checked };
}

test('mesh region: tiles are the same in any order, and shared borders agree bitwise', async () => {
  const A = landTile();
  const tiles = [A, neighborOf(A, 'E', layout), neighborOf(A, 'N', layout), neighborOf(neighborOf(A, 'E', layout), 'N', layout)];
  const a = await run(tiles), b = await run([...tiles].reverse());
  for (const k of a.keys()) assert.equal(Buffer.compare(Buffer.from(a.get(k).fields.elevation.buffer), Buffer.from(b.get(k).fields.elevation.buffer)), 0, k);
  const { bad, checked } = seamMismatches(a, tiles);
  assert.ok(checked > 100);
  assert.equal(bad, 0);
  // the planet's drainage carries on: rivers and sediment cross the shared edges
  const crossing = tiles.flatMap((t) => Object.values(a.get(tileKey(t)).borders).flatMap((x) => x.crossings));
  assert.ok(crossing.some((x) => x.sed && x.sed[0] > 0), 'some crossing carries sediment');
  const r = a.get(tileKey(A));
  assert.ok(r.stats.hack.exponent > 0.4 && r.stats.hack.exponent < 0.8, `Hack ${r.stats.hack.exponent}`);
});

test('mesh region: borders agree across a cube-face edge', async () => {
  const N = tilesPerSide(1, layout);
  const A = { level: 1, face: 2, i: N - 1, j: 7 };
  const B = neighborOf(A, 'E', layout);
  assert.notEqual(A.face, B.face);
  const out = await run([A, B]);
  const { bad, checked } = seamMismatches(out, [A, B]);
  assert.ok(checked > 30);
  assert.equal(bad, 0);
});

test('mesh region in the pipeline: level-2 tiles under two mesh parents agree on the parents\' seam', async () => {
  const lay = { rootTiles: 32, branching: 4, tileCells: 32 };
  const w = new World(21, { layout: lay, region: 'mesh', regionOptions: opts, planetOptions: { subdivisions: 16, duration: 1e8 } });
  const F = w.planet().fields, P = w.planet().mesh.pos;
  let c = 0;
  for (let i = 0; i < F.elevation.length; i++) if (Math.abs(F.elevation[i] - 1500) < Math.abs(F.elevation[c] - 1500)) c = i;
  const A1 = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], lay), B1 = neighborOf(A1, 'E', lay);
  // the child of A1 on its east edge, and the one next to it under B1
  const a2 = childrenOf(A1, lay).find((t) => tileKey(neighborOf(t, 'E', lay)).startsWith(`2/${B1.face}`) && childrenOf(B1, lay).some((u) => tileKey(u) === tileKey(neighborOf(t, 'E', lay))));
  const b2 = neighborOf(a2, 'E', lay);
  const [ra, rb] = [await w.tile(a2), await w.tile(b2)];
  const n = lay.tileCells, m = n + 1;
  for (let j = 0; j <= n; j++) assert.equal(ra.fields.elevation[j * m + n], rb.fields.elevation[j * m], `node ${j}`);
});
