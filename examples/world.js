// Build a world's level 0 and save maps of it:   node examples/world.js [seed] [out-dir]
import { mkdirSync, writeFileSync } from 'node:fs';
import { generatePlanet, MATERIAL_NAMES, spherePixelMap, paint } from '../src/index.js';
import { encodePng } from '../bin/png.js';

const [seed = 'example', out = 'out/world'] = process.argv.slice(2);
mkdirSync(out, { recursive: true });
const t0 = performance.now();
const tile = generatePlanet(seed, {}, (f) => process.stdout.write(`\r  level 0: ${(100 * f).toFixed(0)}%`));
console.log(`\r  level 0 of world "${seed}": ${((performance.now() - t0) / 1000).toFixed(1)} s, ${tile.stats.cells} cells, ` +
  `${tile.events.length} events, ${tile.rivers.nodes.length} river nodes`);

const W = 1200, H = 600;
const pm = spherePixelMap(tile.mesh, W, H);
const n = tile.mesh.size;
const save = (name, colorOf) => {
  const colors = new Uint8ClampedArray(3 * n);
  for (let i = 0; i < n; i++) colors.set(colorOf(i), 3 * i);
  const rgba = new Uint8ClampedArray(4 * W * H);
  paint(pm, colors, rgba, [0, 0, 0]);
  writeFileSync(`${out}/${name}.png`, encodePng(W, H, rgba));
};
const lerp = (a, b, t) => a.map((v, k) => v + (b[k] - v) * Math.max(0, Math.min(1, t)));
const f = tile.fields;
save('elevation', (i) => {
  const z = f.elevation[i];
  if (z <= 0) return lerp([60, 110, 170], [10, 25, 70], -z / 6000);
  return z < 1500 ? lerp([90, 140, 70], [170, 150, 90], z / 1500) : lerp([170, 150, 90], [245, 245, 245], (z - 1500) / 3500);
});
const MATERIAL_COLORS = [[30, 50, 110], [190, 120, 80], [120, 160, 90], [200, 60, 60], [220, 200, 120]];
save('material', (i) => MATERIAL_COLORS[f.material[i]]);
save('drainage', (i) => (f.elevation[i] <= 0 ? [20, 20, 30] : lerp([60, 60, 60], [120, 200, 255], Math.log10(f.drainageArea[i] / 1e10) / 3)));
console.log(`  materials: ${MATERIAL_NAMES.map((m, k) => `${m} ${(100 * f.material.filter((v) => v === k).length / n).toFixed(0)}%`).join(', ')}`);
const kinds = {};
for (const e of tile.events) kinds[e.kind] = (kinds[e.kind] || 0) + 1;
console.log(`  events: ${Object.entries(kinds).map(([k, v]) => `${v} ${k}`).join(', ')}; still active: ${tile.events.filter((e) => e.time === 0).length}`);
console.log(`  maps in ${out}/ (material colours: ${MATERIAL_NAMES.join(', ')} = blue, brown, green, red, sand)`);
