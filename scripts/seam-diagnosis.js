#!/usr/bin/env node
// Diagnosis: is the incongruity between neighbouring region tiles caused by solving
// them separately (walls at the seams), or by the inputs?
// Solves a 2 x 2 block of level-1 tiles (a) as four tiles with seam records and (b) as
// one domain with the same node positions (a layout of 16 root tiles of 512 cells), then
// compares them by eye and with two measures:
//   seam jump     how abruptly the smoothed flow-direction field changes across a line:
//                 mean |v(x - 4) - v(x + 4)| along the internal seams, against the same
//                 lines in the single domain (and lines through the tiles' middles)
//   radial index  mean cosine between flow direction and the direction from each
//                 quadrant's centre: how much drainage organises around tile centres
//   slope step    rms change of slope across the internal seams (m per node)
//
// With --mesh, (c) the same four tiles on the irregular mesh with simulated seams
// (regionmesh.js) join the comparison, as a third panel.
//
//   node scripts/seam-diagnosis.js [seed] [out-dir] [--mesh]
import { mkdirSync, writeFileSync } from 'node:fs';
import { generatePlanet } from '../src/world/l0.js';
import { regionInput, runRegionCPU, regionSeams } from '../src/world/region.js';
import { tileAt, WORLD_LAYOUT } from '../src/world/cubesphere.js';
import { meshRegionContext, meshRegionTile } from '../src/world/regionmesh.js';
import { encodePng } from '../bin/png.js';

const withMesh = process.argv.includes('--mesh');
const [seed = 'example', out = 'out/seams'] = process.argv.slice(2).filter((a) => !a.startsWith('--'));
mkdirSync(out, { recursive: true });
const planet = generatePlanet(seed);
const e = planet.fields.elevation, P = planet.mesh.pos;
let c = 0;
for (let i = 0; i < e.length; i++) if (e[i] > e[c] && e[i] < 3500) c = i;
const t = tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2]);
const i0 = t.i - (t.i % 2), j0 = t.j - (t.j % 2);
const n = WORLD_LAYOUT.tileCells, M = 2 * n + 1; // mosaic nodes per side

// (a) four tiles (and (c) on the mesh)
async function mosaic(make) {
const quad = new Float32Array(M * M), quadRec = new Int32Array(M * M).fill(-1);
let t0 = performance.now();
for (const [di, dj] of [[0, 0], [1, 0], [0, 1], [1, 1]]) {
  const r = await make({ level: 1, face: t.face, i: i0 + di, j: j0 + dj });
  const m = r.nodes, rec = r.drainage.receiver;
  for (let cj = 0; cj < m; cj++) for (let ci = 0; ci < m; ci++) {
    const g = (dj * n + cj) * M + di * n + ci, k = cj * m + ci;
    quad[g] = r.fields.elevation[k];
    const q = rec[k];
    // store the receiver as a mosaic index (it lies inside the same tile)
    quadRec[g] = q < 0 ? -1 : (dj * n + Math.floor(q / m)) * M + di * n + (q % m);
  }
  console.log(`tile ${di},${dj}: Hack ${r.stats.hack.exponent.toFixed(2)}, concavity ${r.stats.slopeArea.concavity.toFixed(2)}`);
}
console.log(`four tiles: ${((performance.now() - t0) / 1000).toFixed(0)} s`);
return { quad, quadRec };
}
const s1 = regionSeams(planet);
const { quad, quadRec } = await mosaic((a) => runRegionCPU(regionInput(planet, a, {}, WORLD_LAYOUT, s1)));
let meshQ = null;
if (withMesh) {
  const mctx = meshRegionContext(planet, {});
  meshQ = await mosaic((a) => meshRegionTile(mctx, a));
}
let t0;
// (b) one domain
const big = { rootTiles: WORLD_LAYOUT.rootTiles / 2, branching: WORLD_LAYOUT.branching, tileCells: 2 * n };
t0 = performance.now();
const rb = await runRegionCPU(regionInput(planet, { level: 1, face: t.face, i: i0 / 2, j: j0 / 2 }, {}, big, regionSeams(planet, {}, big)));
console.log(`one domain: ${((performance.now() - t0) / 1000).toFixed(0)} s, Hack ${rb.stats.hack.exponent.toFixed(2)}, concavity ${rb.stats.slopeArea.concavity.toFixed(2)}`);
const one = rb.fields.elevation, oneRec = rb.drainage.receiver;

// flow direction (unit vector to receiver), smoothed over a (2w+1)^2 box
function flowField(rec, w = 4) {
  const vx = new Float32Array(M * M), vy = new Float32Array(M * M);
  for (let k = 0; k < M * M; k++) {
    const q = rec[k];
    if (q < 0) continue;
    const dx = (q % M) - (k % M), dy = Math.floor(q / M) - Math.floor(k / M), L = Math.hypot(dx, dy) || 1;
    vx[k] = dx / L; vy[k] = dy / L;
  }
  const sx = new Float32Array(M * M), sy = new Float32Array(M * M);
  for (let j = 0; j < M; j++) for (let i = 0; i < M; i++) {
    let ax = 0, ay = 0, cnt = 0;
    for (let b = Math.max(0, j - w); b <= Math.min(M - 1, j + w); b++) for (let a = Math.max(0, i - w); a <= Math.min(M - 1, i + w); a++) { ax += vx[b * M + a]; ay += vy[b * M + a]; cnt++; }
    sx[j * M + i] = ax / cnt; sy[j * M + i] = ay / cnt;
  }
  return { sx, sy };
}
const jump = (f, line, vertical, off = 5) => {
  let s = 0, cnt = 0;
  for (let u = 8; u < M - 8; u++) {
    const a = vertical ? u * M + line - off : (line - off) * M + u;
    const b = vertical ? u * M + line + off : (line + off) * M + u;
    s += Math.hypot(f.sx[a] - f.sx[b], f.sy[a] - f.sy[b]); cnt++;
  }
  return s / cnt;
};
const radial = (rec) => {
  let s = 0, cnt = 0;
  for (const [di, dj] of [[0, 0], [1, 0], [0, 1], [1, 1]]) {
    const cx = di * n + n / 2, cy = dj * n + n / 2;
    for (let j = dj * n + 8; j < dj * n + n - 8; j++) for (let i = di * n + 8; i < di * n + n - 8; i++) {
      const k = j * M + i, q = rec[k];
      if (q < 0) continue;
      const dx = (q % M) - i, dy = Math.floor(q / M) - j, rx = i - cx, ry = j - cy;
      const L = Math.hypot(dx, dy) * Math.hypot(rx, ry);
      if (L > 0) { s += (dx * rx + dy * ry) / L; cnt++; }
    }
  }
  return s / cnt;
};
const runs = [['four tiles', quad, quadRec], ['one domain', one, oneRec]];
if (meshQ) runs.push(['four mesh tiles', meshQ.quad, meshQ.quadRec]);
const rows = [];
for (const [name, , rec] of runs) {
  const f = flowField(rec);
  const seamsJ = (jump(f, n, true) + jump(f, n, false)) / 2;
  const middles = (jump(f, n / 2, true) + jump(f, n / 2, false) + jump(f, n + n / 2, true) + jump(f, n + n / 2, false)) / 4;
  rows.push({ run: name, 'seam jump': seamsJ.toFixed(3), 'jump across tile middles': middles.toFixed(3), 'radial index': radial(rec).toFixed(3) });
}
// elevation and slope steps across the seams (difference of one-sided gradients)
const steps = (z) => {
  let ds = 0, cnt = 0;
  for (let u = 2; u < M - 2; u++) for (const vertical of [true, false]) {
    const at = (k) => (vertical ? u * M + n + k : (n + k) * M + u);
    const left = z[at(0)] - z[at(-1)], right = z[at(1)] - z[at(0)];
    ds += (left - right) ** 2; cnt++;
  }
  return Math.sqrt(ds / cnt);
};
runs.forEach(([, z], k) => { rows[k]['slope step at seam (m per node, rms)'] = steps(z).toFixed(1); });
console.table(rows);

// pictures, side by side
const W = runs.length * (M + 8) - 8, H = M, rgba = new Uint8ClampedArray(4 * W * H).fill(255);
let zmax = 0;
for (const v of one) zmax = Math.max(zmax, v);
const dx = rb.dx;
runs.map(([, z], k) => [z, k * (M + 8)]).forEach(([z, x0]) => {
  for (let j = 0; j < M; j++) for (let i = 0; i < M; i++) {
    const k = (M - 1 - j) * M + i, o = 4 * (j * W + x0 + i), h = z[k];
    let col = h <= 0 ? [45, 85, 150] : [70 + 175 * h / zmax, 115 + 120 * h / zmax, 65 + 160 * h / zmax];
    const hl = z[(M - 1 - j) * M + Math.max(0, i - 1)], hu = z[Math.max(0, M - 2 - j) * M + i];
    const sh = Math.max(0.55, Math.min(1.3, 1 + ((h - hl) + (h - hu)) / (2 * dx) * 6));
    rgba[o] = col[0] * sh; rgba[o + 1] = col[1] * sh; rgba[o + 2] = col[2] * sh; rgba[o + 3] = 255;
  }
});
writeFileSync(`${out}/four-tiles-vs-one-domain.png`, encodePng(W, H, rgba));
console.log(`wrote ${out}/four-tiles-vs-one-domain.png (left to right: ${runs.map((r) => r[0]).join('; ')})`);
