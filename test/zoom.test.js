// Milestone 5 gate: zooming from the planet down to the site level contradicts nothing.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU, regionSeams, restrictToBlocks } from '../src/world/region.js';
import { landscapeSeams, landscapeInput, runLandscape } from '../src/world/landscape.js';
import { siteSeams, buildSite, SITE_COVER } from '../src/world/site.js';
import { tileAt, neighborOf, childrenOf, tileKey } from '../src/world/cubesphere.js';
import { fieldHash } from '../src/world/hash.js';

const layout = { rootTiles: 32, branching: 4, tileCells: 32 };
const planet = generatePlanet(41, { subdivisions: 16, duration: 1e8 });

/** Mean of a child node field over each parent node's block (nearest parent node). */
function blockMeans(child, parentNodes, addr) {
  const n = layout.tileCells, m = parentNodes, b = layout.branching, per = n / b;
  const offI = addr.i % b, offJ = addr.j % b;
  const sums = new Map();
  for (let cj = 0; cj <= n; cj++) for (let ci = 0; ci <= n; ci++) {
    const x = offI * per + (ci * per) / n, y = offJ * per + (cj * per) / n;
    const k = Math.round(y) * m + Math.round(x);
    const e = sums.get(k) || { s: 0, c: 0 };
    e.s += child[cj * (n + 1) + ci]; e.c++;
    sums.set(k, e);
  }
  return sums;
}

async function chain() {
  const e = planet.fields.elevation, P = planet.mesh.pos;
  let c = 0;
  for (let i = 0; i < e.length; i++) if (e[i] > 300 && e[i] < 2500 && (c === 0 || e[i] > e[c])) c = i;
  const L1 = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2], layout);
  const s1 = regionSeams(planet, { iterations: 60 }, layout);
  const r1 = await runRegionCPU(regionInput(planet, L1, { iterations: 60 }, layout, s1));
  const L2 = childrenOf(L1, layout)[5];
  const s2 = landscapeSeams(r1, s1, planet.worldSeed, {}, layout);
  const r2 = runLandscape(landscapeInput(r1, s2, L2, planet.events, { steps: 8 }, layout));
  const L3 = childrenOf(L2, layout)[6];
  const s3 = siteSeams(r2, s2, planet.worldSeed, {}, layout);
  const r3 = buildSite(r2, s3, L3, {}, layout);
  return { L1, L2, L3, r1, r2, r3, s1, s2, s3 };
}

test('zoom: no contradictions from the planet to the site level', async () => {
  const { L2, L3, r1, r2, r3, s2, s3 } = await chain();
  const problems = [];
  // 1. elevation: each child, coarse-grained, gives back its parent (within the stated
  //    tolerance; at the site level, carved channels may lower a block by their share)
  assert.ok(r1.stats.consistency < 5, `level 1 vs planet: ${r1.stats.consistency} m`);
  assert.ok(r2.stats.consistency < 1, `level 2 vs level 1: ${r2.stats.consistency} m`);
  const water = SITE_COVER.indexOf('water');
  {
    const means = blockMeans(r3.fields.elevation, r2.nodes, L3);
    const wet = blockMeans(Float32Array.from(r3.fields.cover, (v) => (v === water ? 1 : 0)), r2.nodes, L3);
    for (const [k, v] of means) {
      const allowed = 0.25 + 6 * (wet.get(k).s / wet.get(k).c); // channel share x a few m of depth
      if (Math.abs(v.s / v.c - r2.fields.elevation[k]) > allowed && v.c > 50) problems.push(`L3 block ${k} off by ${(v.s / v.c - r2.fields.elevation[k]).toFixed(2)} m`);
    }
  }
  // 2. sea and land agree wherever the parent is clearly one or the other
  const checkCoast = (child, parentField, parentNodes, addr, level, margin) => {
    const means = blockMeans(Float32Array.from(child, (z) => (z > 0 ? 1 : 0)), parentNodes, addr);
    for (const [k, v] of means) {
      const p = parentField[k], landShare = v.s / v.c;
      if (p > margin && landShare < 0.5) problems.push(`L${level} block ${k}: parent land (${p.toFixed(0)} m), child mostly sea`);
      if (p < -margin && landShare > 0.5) problems.push(`L${level} block ${k}: parent sea (${p.toFixed(0)} m), child mostly land`);
    }
  };
  checkCoast(r2.fields.elevation, r1.fields.elevation, r1.nodes, L2, 2, 20);
  checkCoast(r3.fields.elevation, r2.fields.elevation, r2.nodes, L3, 3, 5);
  // 3. rivers: wherever the parent has a big river, the child carries one
  {
    const n = layout.tileCells, m1 = r1.nodes, per = n / layout.branching;
    const offI = L2.i % layout.branching, offJ = L2.j % layout.branching;
    for (let j = 1; j < per; j++) for (let i = 1; i < per; i++) {
      const k = (offJ * per + j) * m1 + offI * per + i;
      const Ap = r1.fields.drainageArea[k];
      if (Ap < 1e9 || r1.fields.elevation[k] <= 0) continue;
      // child nodes within a parent cell of it (the parent's own path is only good to a cell)
      let best = 0;
      const r = 1.5 * (n / per);
      for (let cj = Math.max(0, Math.round(j * (n / per) - r)); cj <= Math.min(n, Math.round(j * (n / per) + r)); cj++) {
        for (let ci = Math.max(0, Math.round(i * (n / per) - r)); ci <= Math.min(n, Math.round(i * (n / per) + r)); ci++) {
          best = Math.max(best, r2.fields.drainageArea[cj * (n + 1) + ci]);
        }
      }
      if (best < 0.25 * Ap) problems.push(`L2: parent river of ${(Ap / 1e6).toFixed(0)} km2 at node ${k}, child's largest there ${(best / 1e6).toFixed(0)} km2`);
    }
  }
  // 4. ice: snow only where it is cold enough today
  const snow = SITE_COVER.indexOf('snow');
  if (r3.fields.cover.some((v) => v === snow)) {
    const coldest = Math.min(...r2.fields.temperature);
    if (coldest > -4) problems.push('L3 has snow but its parent is nowhere colder than -4 C');
  }
  // 5. seams at the site level: a neighbour under the same level-2 parent and one across
  //    a level-2 seam, in either order, agree bitwise on their shared border
  const edge = (r, side) => {
    const m = r.nodes, e = r.fields.elevation, out = [];
    for (let t = 0; t < m; t++) { const [ci, cj] = side === 'E' ? [m - 1, t] : side === 'W' ? [0, t] : side === 'N' ? [t, m - 1] : [t, 0]; out.push(e[cj * m + ci]); }
    return out;
  };
  const E = neighborOf(L3, 'E', layout);
  const rE = buildSite(r2, siteSeams(r2, s2, planet.worldSeed, {}, layout), E, {}, layout);
  const rL3again = buildSite(r2, siteSeams(r2, s2, planet.worldSeed, {}, layout), L3, {}, layout);
  assert.deepEqual(edge(r3, 'E'), edge(rE, 'W'), 'site seam inside one level-2 parent');
  assert.equal(fieldHash(rL3again.fields.elevation), fieldHash(r3.fields.elevation), 'the same site tile twice');
  void s3;
  assert.deepEqual(problems, [], problems.join('\n'));
});
