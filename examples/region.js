// Build a world's planet, then a 2 x 2 block of level-1 region tiles around its highest
// land, and save them as one picture (seams included):
//   node examples/region.js [seed] [out-dir]
import { mkdirSync, writeFileSync } from 'node:fs';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU, regionSeams } from '../src/world/region.js';
import { tileAt, tileKey, neighborOf } from '../src/world/cubesphere.js';
import { encodePng } from '../bin/png.js';

const [seed = 'example', out = 'out/region'] = process.argv.slice(2);
mkdirSync(out, { recursive: true });
let t0 = performance.now();
const planet = generatePlanet(seed);
console.log(`planet: ${((performance.now() - t0) / 1000).toFixed(1)} s`);
const e = planet.fields.elevation, P = planet.mesh.pos;
let c = 0;
for (let i = 0; i < e.length; i++) if (e[i] > e[c] && e[i] < 3500) c = i;
const A = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2]);
const B = neighborOf(A, 'E'), C = neighborOf(A, 'N'), D = neighborOf(B, 'N');
const ctx = regionSeams(planet);
const grid = [[C, D], [A, B]]; // rows top to bottom
const res = new Map();
for (const t of [A, B, C, D]) {
  t0 = performance.now();
  const r = await runRegionCPU(regionInput(planet, t, {}, undefined, ctx));
  res.set(tileKey(t), r);
  const s = r.stats;
  console.log(`tile ${tileKey(t)}: ${((performance.now() - t0) / 1000).toFixed(1)} s (CPU reference), Hack ${s.hack.exponent.toFixed(2)}, ` +
    `concavity ${s.slopeArea.concavity.toFixed(2)}, block means within ${s.consistency.toFixed(1)} m, ${r.rivers.crossings.length} crossings`);
}
const m = res.get(tileKey(A)).nodes;
const W = 2 * m, H = 2 * m, rgba = new Uint8ClampedArray(4 * W * H);
let zmax = 0;
for (const r of res.values()) for (const v of r.fields.elevation) zmax = Math.max(zmax, v);
grid.forEach((row, gy) => row.forEach((t, gx) => {
  const r = res.get(tileKey(t)), el = r.fields.elevation, ar = r.fields.drainageArea;
  for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) {
    const k = (m - 1 - j) * m + i, o = 4 * ((gy * m + j) * W + gx * m + i), z = el[k];
    let col = z <= 0 ? [45, 85, 150] : [70 + 175 * z / zmax, 115 + 120 * z / zmax, 65 + 160 * z / zmax];
    if (ar[k] > 5e7 && z > 0) { const tt = Math.min(1, Math.log10(ar[k] / 5e7) / 2.5); col = [40, 90 + 100 * tt, 205]; }
    const zl = el[(m - 1 - j) * m + Math.max(0, i - 1)], zu = el[Math.max(0, m - 2 - j) * m + i];
    const sh = Math.max(0.55, Math.min(1.3, 1 + ((z - zl) + (z - zu)) / (2 * r.dx) * 6));
    rgba[o] = col[0] * sh; rgba[o + 1] = col[1] * sh; rgba[o + 2] = col[2] * sh; rgba[o + 3] = 255;
  }
}));
writeFileSync(`${out}/region-2x2.png`, encodePng(W, H, rgba));
console.log(`wrote ${out}/region-2x2.png (${W} x ${H}, ${(res.get(tileKey(A)).dx / 1000).toFixed(2)} km per pixel)`);
