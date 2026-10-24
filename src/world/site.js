// Level 3: a site tile (~305 m, 1.2-m nodes) refined from its level-2 parent. No time
// stepping: every feature is placed by rules conditioned on what the parent says.
//
//   terrain      the parent's surface (bicubic) plus fine relief noise, block means then
//                returned to the parent's
//   channels     the parent's drainage links with enough area become channels: a
//                smoothed centre line (Chaikin), width W = 4 Q^0.5 and depth D = 0.4 Q^0.35
//                from bankfull discharge (Leopold & Maddock 1953 hydraulic geometry), a
//                parabolic bed, and gravel bars along the banks
//   debris       talus on slopes below the parent's landslide scars and deposits, bare rock
//                where the exposed bed stands steeper than soil can; boulders as a point
//                process, denser on talus and rock, sizes from a power law
//   ice and snow where it is cold enough now
//   trees        a point process on a global jittered grid: acceptance and height from
//                temperature, rainfall, slope and ground, clustered by a noise field; each
//                tree a crown in a canopy-height raster
//
// Everything that decides whether a feature exists is keyed on global positions (face
// coordinates or the sphere), so a tree or boulder near a seam is the same tree from
// both tiles. The border nodes come from the shared edge records (seams.js), as at the
// levels above.

import { nodePosition, tileRect, WORLD_LAYOUT } from './cubesphere.js';
import { worldKey, deriveKey, mix32 } from './hash.js';
import { fbm3k } from './noise.js';
import { landscapeSampler, landscapeSeams } from './landscape.js';
import { tileBorder, sideNode, SIDES } from './seams.js';
import { correctBlockMeans } from './region.js';
import { chaikin } from '../meander.js';

const LAPSE = 0.0065;

export const SITE_COVER = ['sea', 'water', 'gravel', 'meadow', 'forest', 'shrub', 'rock', 'talus', 'snow', 'bare', 'beach'];
const C = Object.fromEntries(SITE_COVER.map((n, i) => [n, i]));

export const SITE_DEFAULTS = {
  relief: { base: 0.25, perSlope: 0.2, wavelength: 40, octaves: 4, gain: 0.5 },
  channelArea: 2e5,             // m^2 for a channel
  runoff: 0.5,                  // share of rain that runs off
  bankfull: 3,                  // bankfull discharge over mean discharge
  treeSpacing: 2.5,             // m between tree candidates
  boulderSpacing: 4,            // m between boulder candidates
};

function crWeights(t, w) {
  const t2 = t * t, t3 = t2 * t;
  w[0] = -0.5 * t3 + t2 - 0.5 * t; w[1] = 1.5 * t3 - 2.5 * t2 + 1; w[2] = -1.5 * t3 + 2 * t2 + 0.5 * t; w[3] = 0.5 * t3 - 0.5 * t2;
  return w;
}

function bilinear(f, m, x, y) {
  const xi = Math.min(m - 2, Math.max(0, Math.floor(x))), yi = Math.min(m - 2, Math.max(0, Math.floor(y)));
  const u = Math.min(1, Math.max(0, x - xi)), v = Math.min(1, Math.max(0, y - yi));
  return (1 - u) * (1 - v) * f[yi * m + xi] + u * (1 - v) * f[yi * m + xi + 1] + (1 - u) * v * f[(yi + 1) * m + xi] + u * v * f[(yi + 1) * m + xi + 1];
}

/** Fractional parent-node coordinates of node (ci, cj) of addr. */
function toParent(addr, n, ci, cj, layout) {
  const b = layout.branching, per = layout.tileCells / b;
  const offI = addr.i - b * Math.floor(addr.i / b), offJ = addr.j - b * Math.floor(addr.j / b);
  return [offI * per + (ci * per) / n, offJ * per + (cj * per) / n];
}

/** Uniform [0, 1) from a key and integer coordinates (fast, stateless). */
const u01 = (key, a, b, c = 0) => mix32(key, a, b, c) / 4294967296;

/**
 * The seams of the site tiles under one level-2 parent (built with the site's relief
 * settings, so the records and the interiors use the same noise).
 * @param {object} parent level-2 result
 * @param {object} parentSeams the level-2 seam context the parent was built with
 */
export function siteSeams(parent, parentSeams, worldSeed, options = {}, layout = WORLD_LAYOUT) {
  const relief = { ...SITE_DEFAULTS.relief, ...(options.relief || {}) };
  return landscapeSeams(parent, parentSeams, worldSeed, { relief }, layout);
}

/**
 * Build a site tile.
 * @param {object} parent level-2 result (runLandscape) of parentOf(addr)
 * @param {object} ctx landscapeSeams(parent, level-2 seam context, ...) — the seams of level 3
 */
export function buildSite(parent, ctx, addr, options = {}, layout = WORLD_LAYOUT) {
  const o = { ...SITE_DEFAULTS, ...options, relief: { ...SITE_DEFAULTS.relief, ...(options.relief || {}) } };
  const smp = landscapeSampler(parent, ctx.worldSeed, { relief: o.relief, radius: ctx.radius });
  const n = layout.tileCells, nx = n + 1, N = nx * nx, m = parent.nodes;
  const PF = parent.fields;
  const R = ctx.radius;
  const a = nodePosition(addr, n, n / 2, n / 2, undefined, layout), b2 = nodePosition(addr, n, n / 2 + 1, n / 2, undefined, layout);
  const dx = R * Math.acos(Math.min(1, a[0] * b2[0] + a[1] * b2[1] + a[2] * b2[2]));
  const z = new Float32Array(N), base = new Float32Array(N), cover = new Uint8Array(N), canopy = new Float32Array(N);
  const block = new Int32Array(N), locIds = new Int32Array(16 * N), locW = new Float32Array(16 * N);
  const keep = new Uint8Array(N);
  const v = [0, 0, 0], wx = [0, 0, 0, 0], wy = [0, 0, 0, 0];
  const x2s = new Float32Array(N), y2s = new Float32Array(N);
  for (let cj = 0; cj <= n; cj++) for (let ci = 0; ci <= n; ci++) {
    const c = cj * nx + ci;
    nodePosition(addr, n, ci, cj, v, layout);
    const [x2, y2] = toParent(addr, n, ci, cj, layout);
    x2s[c] = x2; y2s[c] = y2;
    const s = smp.sample(x2, y2, v);
    z[c] = s.z; base[c] = s.e;
    block[c] = Math.min(m - 1, Math.round(y2)) * m + Math.min(m - 1, Math.round(x2));
    const fx = Math.floor(x2), fy = Math.floor(y2);
    crWeights(x2 - fx, wx); crWeights(y2 - fy, wy);
    for (let bb = 0; bb < 4; bb++) for (let aa = 0; aa < 4; aa++) {
      locIds[16 * c + 4 * bb + aa] = Math.min(m - 1, Math.max(0, fy - 1 + bb)) * m + Math.min(m - 1, Math.max(0, fx - 1 + aa));
      locW[16 * c + 4 * bb + aa] = wx[aa] * wy[bb];
    }
  }
  // border from the shared records
  const border = tileBorder(ctx, addr);
  for (const side of SIDES) for (let t = 0; t <= n; t++) {
    const [ci, cj] = sideNode(side, t, n);
    z[cj * nx + ci] = border[side].z[t];
    keep[cj * nx + ci] = 1;
  }

  // block means back to the parent's before anything is carved: a channel is detail the
  // parent could not resolve, and pushing the carved volume into the banks (as a
  // correction afterwards would) raises walls along every river
  const consistency = correctBlockMeans(z, { block, locIds, locW, base, locK: 16 }, keep, 40, 64);

  // channels from the parent's drainage
  const rec = parent.drainage.receiver, area = parent.drainage.area;
  const per = layout.tileCells / layout.branching;
  const offI = addr.i - layout.branching * Math.floor(addr.i / layout.branching);
  const offJ = addr.j - layout.branching * Math.floor(addr.j / layout.branching);
  const toMe = (pc) => { const pi = pc % m, pj = (pc - pi) / m; return [((pi - offI * per) * n) / per, ((pj - offJ * per) * n) / per]; };
  const inView = (pc, margin) => { const [xi, yi] = toMe(pc); return xi >= -margin && yi >= -margin && xi <= n + margin && yi <= n + margin; };
  const isChan = (pc) => pc >= 0 && area[pc] >= o.channelArea;
  const visited = new Set();
  const paths = [];
  const margin = 3 * n / per;
  for (let pc = 0; pc < rec.length; pc++) {
    if (!isChan(pc) || !inView(pc, margin)) continue;
    // a head here: no channel donor among the nodes in view
    let hasDonor = false;
    const pi = pc % m, pj = (pc - pi) / m;
    for (let dj = -1; dj <= 1 && !hasDonor; dj++) for (let di = -1; di <= 1; di++) {
      const q = (pj + dj) * m + pi + di;
      if (q >= 0 && q < rec.length && (di || dj) && rec[q] === pc && isChan(q) && inView(q, margin)) { hasDonor = true; break; }
    }
    if (hasDonor) continue;
    const path = [];
    let q = pc;
    while (q >= 0 && isChan(q) && inView(q, margin + 2)) {
      const [xi, yi] = toMe(q);
      // bankfull discharge (m^3/s) from drainage area and rainfall
      const Q = (area[q] * (PF.precip[q] / 1000) * o.runoff * o.bankfull) / 3.15576e7;
      path.push({ x: xi, y: yi, w: Math.max(0.6, 4 * Math.sqrt(Q)) / dx, d: Math.max(0.2, 0.4 * Q ** 0.35), z: PF.elevation[q] });
      if (visited.has(q)) break;
      visited.add(q);
      q = rec[q];
    }
    if (path.length >= 2) paths.push(path);
  }
  for (const raw of paths) {
    const line = chaikin(raw, 3);
    // water surface falls monotonically downstream
    for (let k = 1; k < line.length; k++) line[k].z = Math.min(line[k].z, line[k - 1].z);
    for (let k = 0; k < line.length - 1; k++) stampSegment(line[k], line[k + 1]);
  }
  function stampSegment(p, q) {
    const hw = Math.max(p.w, q.w) / 2 + 2;
    const x0 = Math.max(1, Math.floor(Math.min(p.x, q.x) - hw)), x1 = Math.min(n - 1, Math.ceil(Math.max(p.x, q.x) + hw));
    const y0 = Math.max(1, Math.floor(Math.min(p.y, q.y) - hw)), y1 = Math.min(n - 1, Math.ceil(Math.max(p.y, q.y) + hw));
    const ex = q.x - p.x, ey = q.y - p.y, L2 = ex * ex + ey * ey || 1;
    for (let cj = y0; cj <= y1; cj++) for (let ci = x0; ci <= x1; ci++) {
      const t = Math.max(0, Math.min(1, ((ci - p.x) * ex + (cj - p.y) * ey) / L2));
      const d = Math.hypot(ci - (p.x + t * ex), cj - (p.y + t * ey));
      const w = (p.w + t * (q.w - p.w)) / 2;
      const zw = p.z + t * (q.z - p.z);
      const depth = p.d + t * (q.d - p.d);
      const c = cj * nx + ci;
      if (d < w) {
        const bed = zw - depth * (1 - (d / w) ** 2);
        if (bed < z[c]) z[c] = Math.fround(bed);
        cover[c] = C.water;
        keep[c] = 1;
      } else if (d < w + 1.5 && cover[c] !== C.water) {
        // banks and bars, graded down to the water's edge
        const bank = zw + 0.3 * (d - w);
        if (bank < z[c]) z[c] = Math.fround(bank);
        cover[c] = C.gravel;
      }
    }
  }


  // ground cover
  const wkey = deriveKey(worldKey(ctx.worldSeed), addr.level, 0x67726f75);
  const rr = tileRect(addr, layout);
  const hard = PF.hardRock, slideRec = parent.record?.landslide, depRec = parent.record?.deposit;
  const tempNow = new Float32Array(N);
  for (let cj = 0; cj <= n; cj++) for (let ci = 0; ci <= n; ci++) {
    const c = cj * nx + ci, x2 = x2s[c], y2 = y2s[c];
    const T = bilinear(PF.temperature, m, x2, y2) - LAPSE * (z[c] - bilinear(PF.elevation, m, x2, y2));
    tempNow[c] = T;
    if (cover[c] === C.water || cover[c] === C.gravel) continue;
    if (z[c] <= 0) { cover[c] = C.sea; continue; }
    const ip = Math.min(nx - 1, ci + 1), im = Math.max(0, ci - 1), jp = Math.min(nx - 1, cj + 1), jm = Math.max(0, cj - 1);
    const slope = Math.hypot(z[cj * nx + ip] - z[cj * nx + im], z[jp * nx + ci] - z[jm * nx + ci]) / ((ip - im) * dx || dx);
    const slid = slideRec ? bilinear(slideRec, m, x2, y2) : 0;
    const dep = depRec ? bilinear(depRec, m, x2, y2) : 0;
    const hardHere = hard ? hard[Math.round(y2) * m + Math.round(x2)] : 0;
    if (z[c] < 2.5 && nearSea(c)) cover[c] = C.beach;
    else if (T < -4) cover[c] = C.snow;
    else if (slope > (hardHere ? 0.9 : 1.2)) cover[c] = C.rock;
    else if ((slid > 0.3 || dep > 0.8) && slope > 0.2) cover[c] = C.talus;
    else {
      const P = bilinear(PF.precip, m, x2, y2);
      cover[c] = P < 300 ? C.bare : T < 0 ? C.shrub : C.meadow;
    }
  }
  function nearSea(c) {
    const ci = c % nx, cj = (c - ci) / nx;
    for (let dj = -3; dj <= 3; dj++) for (let di = -3; di <= 3; di++) {
      const q = (cj + dj) * nx + ci + di;
      if (ci + di >= 0 && cj + dj >= 0 && ci + di < nx && cj + dj < nx && z[q] <= 0) return true;
    }
    return false;
  }

  // trees and boulders on global jittered grids (face coordinates), so a seam splits nothing
  const objects = { trees: [], boulders: [] };
  const faceSpan = rr.size; // face coordinate width of the tile
  const metresPerFace = (dx * n) / faceSpan;
  const clusterKey = deriveKey(worldKey(ctx.worldSeed), 3, 0x74726565);
  const scatter = (spacing, salt, accept) => {
    const g = spacing / metresPerFace;
    const i0 = Math.floor(rr.a0 / g), i1 = Math.ceil((rr.a0 + faceSpan) / g);
    const j0 = Math.floor(rr.b0 / g), j1 = Math.ceil((rr.b0 + faceSpan) / g);
    for (let gj = j0; gj < j1; gj++) for (let gi = i0; gi < i1; gi++) {
      const fa = (gi + u01(wkey, gi, gj, salt)) * g, fb = (gj + u01(wkey, gi, gj, salt + 1)) * g;
      const ci = ((fa - rr.a0) / faceSpan) * n, cj = ((fb - rr.b0) / faceSpan) * n;
      if (ci < 0 || cj < 0 || ci >= n || cj >= n) continue; // belongs to a neighbour
      accept(ci, cj, u01(wkey, gi, gj, salt + 2), u01(wkey, gi, gj, salt + 3));
    }
  };
  scatter(o.treeSpacing, 10, (ci, cj, r1, r2) => {
    const c = Math.round(cj) * nx + Math.round(ci);
    const cv = cover[c];
    if (cv !== C.meadow && cv !== C.shrub) return;
    const T = tempNow[c], P = bilinear(PF.precip, m, x2s[c], y2s[c]);
    // forest where warm enough and wet enough; thinner towards the tree line and the dry edge
    const warm = Math.max(0, Math.min(1, (T + 1) / 8)), wet = Math.max(0, Math.min(1, (P - 350) / 600));
    nodePosition(addr, n, ci, cj, v, layout);
    const clump = 0.5 + 0.5 * fbm3k(clusterKey, v[0] * R, v[1] * R, v[2] * R, 60, 3);
    const p = warm * wet * Math.min(1, 1.6 * clump);
    if (r1 >= p) return;
    const h = (6 + 24 * warm * wet) * (0.6 + 0.4 * r2);
    const rad = 0.18 * h + 0.6;
    objects.trees.push({ x: ci * dx, y: cj * dx, height: +h.toFixed(2), crown: +rad.toFixed(2) });
    stampCrown(ci, cj, h, rad / dx);
  });
  scatter(o.boulderSpacing, 20, (ci, cj, r1, r2) => {
    const c = Math.round(cj) * nx + Math.round(ci);
    const cv = cover[c];
    const p = cv === C.talus ? 0.6 : cv === C.rock ? 0.25 : cv === C.gravel ? 0.1 : cv === C.meadow || cv === C.shrub ? 0.01 : 0;
    if (r1 >= p) return;
    const radius = 0.25 * (1 - r2) ** -0.6; // Pareto tail: many small, a few big
    objects.boulders.push({ x: ci * dx, y: cj * dx, radius: +Math.min(4, radius).toFixed(2) });
  });
  function stampCrown(ci, cj, h, rn) {
    const r0 = Math.ceil(rn);
    for (let dj = -r0; dj <= r0; dj++) for (let di = -r0; di <= r0; di++) {
      const a2 = Math.round(ci) + di, b3 = Math.round(cj) + dj;
      if (a2 < 0 || b3 < 0 || a2 >= nx || b3 >= nx) continue;
      const q = (di * di + dj * dj) / (rn * rn);
      if (q > 1) continue;
      const k = b3 * nx + a2;
      canopy[k] = Math.max(canopy[k], h * Math.sqrt(1 - q));
      if (cover[k] === C.meadow || cover[k] === C.shrub) cover[k] = C.forest;
    }
  }
  return {
    addr, n, nodes: nx, dx, border,
    fields: { elevation: z, cover, canopy },
    objects, stats: { consistency, channels: paths.length, trees: objects.trees.length, boulders: objects.boulders.length },
  };
}
