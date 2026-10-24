import { test } from 'node:test';
import assert from 'node:assert/strict';
import { World, memoryStore } from '../src/world/pipeline.js';
import { tileAt, neighborOf, childrenOf, tileKey } from '../src/world/cubesphere.js';
import { fieldHash } from '../src/world/hash.js';

const layout = { rootTiles: 32, branching: 4, tileCells: 32 };
const planetOptions = { subdivisions: 16, duration: 1e8 };
const hashes = (r) => Object.fromEntries(Object.entries(r.fields).map(([k, v]) => [k, fieldHash(v)]));

function siteAddress(world) {
  const p = world.planet(), e = p.fields.elevation, P = p.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > 200 && e[i] < 2000 && (c === 0 || e[i] > e[c])) c = i;
  const l1 = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
  return childrenOf(childrenOf(l1, layout)[3], layout)[9];
}

test('pipeline: a stored tile is exactly the generated one, and prefetch builds neighbours', async () => {
  const store = memoryStore();
  const a = new World(51, { layout, store, planetOptions });
  const addr = siteAddress(a);
  const first = await a.tile(addr);
  assert.equal(a.stats.generated, 3, 'site, landscape and region tiles generated');
  // a second world on the same store: no generation, identical fields
  const b = new World(51, { layout, store, planetOptions });
  const again = await b.tile(addr);
  assert.equal(b.stats.generated, 0);
  assert.equal(b.stats.storeHits, 1);
  assert.deepEqual(hashes(again), hashes(first));
  // a third world with no store at all regenerates the same tile
  const c = new World(51, { layout, planetOptions });
  assert.deepEqual(hashes(await c.tile(addr)), hashes(first));
  // prefetch: the four neighbours get built in the background
  const nbs = ['W', 'E', 'S', 'N'].map((s) => neighborOf(addr, s, layout));
  a.prefetch(nbs);
  const deadline = Date.now() + 120000;
  while (nbs.some((t) => !a.lru.has(a.key(t))) && Date.now() < deadline) await new Promise((r) => setTimeout(r, 20));
  for (const t of nbs) assert.ok(a.lru.has(a.key(t)), `${tileKey(t)} prefetched`);
});
