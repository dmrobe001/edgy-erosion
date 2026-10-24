// Level 1: a region tile (~310 km, 1.2-km cells) refined from the planet tile.
//
// The tile is a raster of tileCells x tileCells cells plus a halo of `halo` cells on
// every side; the halo is computed and then discarded, so the tile's processes never
// see a hard edge directly. Initial condition: the planet's fields interpolated to each
// cell (barycentric on the icosahedral mesh), plus subgrid relief as seeded noise
// evaluated at each cell's position on the sphere (so neighbouring tiles agree). The
// noise is the unresolved detail; erosion turns it into drainage structure.
//
// Erodibility comes from the planet's material class (blended across class borders)
// times a climate factor; uplift and rainfall come from the planet. The sea (cells
// below sea level) and, until the seam records exist, the halo's outer ring are the
// outlets.

import { nodePosition, WORLD_LAYOUT } from './cubesphere.js';
import { seamContext, tileBorder, sideNode, SIDES, locateInTile } from './seams.js';
import { worldKey, deriveKey, uniform } from './hash.js';
import { fbm3k } from './noise.js';
import { JacobiSolver, FLAG_OUTLET, FLAG_FIXED } from './jacobi.js';
import { drainageStatistics, extractRiverGraph } from './rivers.js';
import { RasterSolver } from './fastscape.js';

const LAPSE = 0.0065; // K per m, as at level 0

export const REGION_DEFAULTS = {
  halo: 16,                     // cells computed around the tile and discarded
  duration: 5e6,                // yr simulated: a few drainage response times
  iterations: 500,              // parallel iterations (dt = duration / iterations)
  m: 0.5,                       // area exponent (n = 1)
  K: 4e-6,                      // erodibility for 'terrane' rock, m^(1-2m)/yr
  // erodibility relative to terrane, by material class: oceanic, craton, terrane, arc, basin
  materialK: [1, 0.4, 1, 1.6, 3],
  D: 0.05,                      // creep diffusivity, m^2/yr (barely felt at 1.2-km cells)
  eps: 1e-5,                    // gradient across filled depressions (lowland trunks run at 1e-4 to 1e-5)
  relief: { base: 40, perMetre: 0.12, wavelength: 40000, octaves: 5, gain: 0.5 }, // subgrid noise amplitude (m) and spectrum
  channelArea: 2e7,             // m^2 drained for a channel in the river graph and statistics
  nudgeEvery: 25,               // iterations between consistency adjustments of erodibility
  nudgeTime: 0.1,               // integral timescale of that control, as a fraction of the duration
  kHeight: 100,                 // m of block-mean error that changes erodibility by a factor e
  kRange: 1000,                 // erodibility may be scaled by 1/kRange ... kRange
  minorSpacing: 0,              // nodes between minor stream crossings of tile edges (0: none; see docs/PLAN.md)
};

/**
 * What the planet says at a unit vector, for level 1: smooth fields interpolated
 * barycentrically, and the initial elevation with subgrid noise.
 */
export function regionSampler(planet, options = {}, { warmStart = true } = {}) {
  const o = { ...REGION_DEFAULTS, ...options, relief: { ...REGION_DEFAULTS.relief, ...(options.relief || {}) } };
  const mesh = planet.mesh, R = mesh.radius, F = planet.fields;
  const noiseKey = deriveKey(worldKey(planet.worldSeed), 1, 0x6e6f6973); // one field for the whole level
  let hint = -1;
  const amplitude = (e) => o.relief.base + o.relief.perMetre * Math.max(0, e);
  // smooth interpolation: a Gaussian-weighted mean over the nearest cell's two rings of
  // neighbours (sigma 0.5 cell spacings). Linear interpolation on the triangles would
  // leave straight creases along every triangle edge, 160 km long.
  const sig2 = 2 * (0.5 * mesh.spacing / R) ** 2;
  const rings = new Map();
  const ringsOf = (cell) => {
    let r = rings.get(cell);
    if (r) return r;
    const set = new Set([cell]);
    for (let pass = 0; pass < 2; pass++) for (const c of [...set]) for (let k = 0; k < mesh.valence[c]; k++) { const j = mesh.neighbors[6 * c + k]; if (j >= 0) set.add(j); }
    r = Int32Array.from([...set].sort((a, b) => a - b));
    rings.set(cell, r);
    return r;
  };
  const Pm = mesh.pos;
  // first-order correction of the smoothing, per planet cell (cached per planet): the
  // cell's bias (its value minus the plain smooth field's mean over it) is interpolated
  // with the same kernel and added, which restores most of what smoothing takes off
  // peaks and adds to basins
  if (!biasCache.has(planet)) biasCache.set(planet, new Map());
  const biases = biasCache.get(planet);
  const cellBias = (j) => {
    let b = biases.get(j);
    if (b === undefined) { b = F.elevation[j] - quadrature(j, (p) => plain(p)); biases.set(j, b); }
    return b;
  };
  function plain(v) {
    const loc = mesh.locate(v[0], v[1], v[2], -1);
    let e = 0, wsum = 0;
    for (const id of ringsOf(loc.cell)) {
      const wgt = Math.exp(-((Pm[3 * id] - v[0]) ** 2 + (Pm[3 * id + 1] - v[1]) ** 2 + (Pm[3 * id + 2] - v[2]) ** 2) / sig2);
      wsum += wgt; e += wgt * F.elevation[id];
    }
    return e / wsum;
  }
  function sample(v) {
    // warm-starting the mesh walk from the last point is only for speed; records shared
    // between tiles start cold, so nothing depends on which tile asked first
    const loc = mesh.locate(v[0], v[1], v[2], warmStart ? hint : -1);
    if (warmStart) hint = loc.cell;
    let e = 0, u = 0, p = 0, k = 0, tsl = 0, wsum = 0;
    for (const id of ringsOf(loc.cell)) {
      const dxv = Pm[3 * id] - v[0], dyv = Pm[3 * id + 1] - v[1], dzv = Pm[3 * id + 2] - v[2];
      const wgt = Math.exp(-(dxv * dxv + dyv * dyv + dzv * dzv) / sig2);
      wsum += wgt;
      e += wgt * (F.elevation[id] + cellBias(id));
      u += wgt * F.uplift[id];
      p += wgt * F.precip[id];
      k += wgt * o.materialK[F.material[id]];
      tsl += wgt * (F.temperature[id] + LAPSE * Math.max(0, F.elevation[id])); // at sea level
    }
    e /= wsum; u /= wsum; p /= wsum; k /= wsum; tsl /= wsum;
    // subgrid relief: rougher where the planet's surface stands higher. It is damped,
    // not removed, near sea level: enough to give the coast bays and headlands, while
    // the planet's coastline stays where it was on average
    const amp = amplitude(e);
    const noise = fbm3k(noiseKey, v[0] * R, v[1] * R, v[2] * R, o.relief.wavelength, o.relief.octaves, o.relief.gain);
    const coastFade = Math.min(1, 0.5 + Math.abs(e) / (2 * amp));
    const z = Math.fround(e + amp * noise * coastFade);
    return { loc, e, u, p, k, z, amp, tsl };
  }
  /**
   * How far the smooth field's mean over planet cell b falls from the cell's own value
   * (13-point quadrature: centre, and points 0.4 and 0.8 of the way to each edge
   * midpoint, weighted by the area they stand for). Smoothing flattens peaks and fills
   * basins; this bias is added to each block's target, so a child that matches its
   * targets gives back the parent exactly where it covers a whole cell.
   */
  function quadrature(b, f) {
    const c = [Pm[3 * b], Pm[3 * b + 1], Pm[3 * b + 2]];
    const along = (a, d, t) => { const x = a[0] + t * (d[0] - a[0]), y = a[1] + t * (d[1] - a[1]), z = a[2] + t * (d[2] - a[2]); const L = Math.hypot(x, y, z); return [x / L, y / L, z / L]; };
    let q = 0.04 * f(c), ring = 0;
    const mids = [];
    for (let k = 0; k < mesh.valence[b]; k++) { const j = mesh.neighbors[6 * b + k]; if (j >= 0) mids.push(along(c, [Pm[3 * j], Pm[3 * j + 1], Pm[3 * j + 2]], 0.5)); }
    for (const mp of mids) ring += 0.32 * f(along(c, mp, 0.4)) + 0.64 * f(along(c, mp, 0.8));
    return q + ring / mids.length;
  }
  function blockBias(b) {
    return F.elevation[b] - quadrature(b, (p) => sample(p).e);
  }
  return { sample, amplitude, blockBias, options: o };
}

const seamCache = new WeakMap();
const biasCache = new WeakMap();

/** The seam context of level-1 tiles under a planet (records are cached per planet: they are pure). */
export function regionSeams(planet, options = {}, layout = WORLD_LAYOUT, { fresh = false } = {}) {
  const k = JSON.stringify([options, layout]);
  if (!fresh && seamCache.get(planet)?.has(k)) return seamCache.get(planet).get(k);
  const smp = regionSampler(planet, options, { warmStart: false });
  const P = planet.mesh.pos, rec = planet.flow.receiver, area = planet.fields.drainageArea;
  const links = [];
  for (let i = 0; i < rec.length; i++) {
    const r = rec[i];
    if (r < 0) continue;
    links.push({ id: i, to_id: r, from: [P[3 * i], P[3 * i + 1], P[3 * i + 2]], to: [P[3 * r], P[3 * r + 1], P[3 * r + 2]], area: area[i] });
  }
  const ctx = seamContext({
    level: 1, worldSeed: planet.worldSeed, layout, radius: planet.mesh.radius, linkSpacing: planet.mesh.spacing, links,
    minorSpacing: smp.options.minorSpacing,
    elevation: (v) => smp.sample(v).z,
    base: (v) => smp.sample(v).e,
    notchDepth: (v) => 0.5 * smp.amplitude(smp.sample(v).e),
  });
  if (!fresh) {
    if (!seamCache.has(planet)) seamCache.set(planet, new Map());
    seamCache.get(planet).set(k, ctx);
  }
  return ctx;
}

/**
 * Inputs of a region tile from the planet tile.
 *
 * The tile is (n + 1) x (n + 1) nodes, its border nodes on the tile's edges. With
 * `seams` (the default) the border takes its elevations and river crossings from the
 * shared edge records (seams.js): the border is a fixed wall except where a parent river
 * crosses it. An outflow crossing is a sink; an inflow crossing injects the parent's
 * upstream drainage area at the node inside it, and a channel is pre-carved from there
 * along the parent's flow path to where that river leaves the tile (or reaches the sea).
 * Without seams (seams: false) the tile has a halo of `halo` nodes whose outer ring is
 * open: each tile then drains wherever its halo edge is lowest.
 *
 * @param {object} planet level-0 tile record (generatePlanet)
 * @param {{level:1, face:number, i:number, j:number}} addr
 */
export function regionInput(planet, addr, options = {}, layout = WORLD_LAYOUT, seamCtx) {
  const smp = regionSampler(planet, options);
  const o = smp.options;
  if (addr.level !== 1) throw new Error('regionInput builds level-1 tiles');
  const seams = o.seams !== false;
  const n = layout.tileCells, h = seams ? 0 : o.halo, nx = n + 1 + 2 * h, ny = nx, N = nx * ny;
  const R = planet.mesh.radius;
  const key = deriveKey(worldKey(planet.worldSeed), addr.level, addr.face, addr.i, addr.j);
  // node spacing at the tile's centre
  const a = nodePosition(addr, n, n / 2, n / 2, undefined, layout), b = nodePosition(addr, n, n / 2 + 1, n / 2, undefined, layout);
  const dx = R * Math.acos(Math.min(1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]));
  const z = new Float32Array(N), U = new Float32Array(N), K = new Float32Array(N), P = new Float32Array(N), inflow = new Float32Array(N);
  const flags = new Uint32Array(N), material = new Uint8Array(N), pos = new Float64Array(3 * N);
  // the parent's view of each node: its planet cell (block), barycentric interpolation
  // weights, and the interpolated smooth elevation (what the parent says, before noise)
  const block = new Int32Array(N), locIds = new Int32Array(3 * N), locW = new Float32Array(3 * N), base = new Float32Array(N);
  const tempSL = new Float32Array(N);
  const v = [0, 0, 0];
  for (let cj = 0; cj < ny; cj++) {
    for (let ci = 0; ci < nx; ci++) {
      const c = cj * nx + ci;
      nodePosition(addr, n, ci - h, cj - h, v, layout);
      pos[3 * c] = v[0]; pos[3 * c + 1] = v[1]; pos[3 * c + 2] = v[2];
      const sm = smp.sample(v);
      material[c] = planet.fields.material[sm.loc.cell];
      block[c] = sm.loc.cell;
      for (let t = 0; t < 3; t++) { locIds[3 * c + t] = sm.loc.ids[t]; locW[3 * c + t] = sm.loc.weights[t]; }
      base[c] = sm.e;
      tempSL[c] = sm.tsl;
      z[c] = sm.z;
      U[c] = sm.u;
      P[c] = Math.max(0.05, sm.p / 1000);
      K[c] = o.K * sm.k;
      if (z[c] <= 0) flags[c] |= FLAG_OUTLET;
      if (!seams && (ci === 0 || cj === 0 || ci === nx - 1 || cj === ny - 1)) flags[c] |= FLAG_OUTLET | FLAG_FIXED;
    }
  }
  let border = null;
  if (seams) {
    const ctx = seamCtx || regionSeams(planet, options, layout);
    border = tileBorder(ctx, addr);
    const idx = (ci, cj) => cj * nx + ci;
    for (const side of SIDES) {
      const bd = border[side];
      for (let t = 0; t <= n; t++) {
        const c = idx(...sideNode(side, t, n));
        z[c] = bd.z[t];
        flags[c] = FLAG_FIXED | (bd.z[t] <= 0 ? FLAG_OUTLET : 0);
      }
    }
    const inward = { W: [1, 0], E: [-1, 0], S: [0, 1], N: [0, -1] };
    // outflow crossings first (the carved channels end at them)
    const exits = new Map(); // parent link id -> border node index
    for (const side of SIDES) for (const x of border[side].crossings) {
      if (x.inflow) continue;
      const c = idx(...sideNode(side, x.node, n));
      flags[c] = FLAG_FIXED | FLAG_OUTLET;
      for (const id of x.links) exits.set(id, c);
    }
    const rec = planet.flow.receiver, Pm = planet.mesh.pos, E = planet.fields.elevation;
    for (const side of SIDES) for (const x of border[side].crossings) {
      if (!x.inflow) continue;
      const [bi, bj] = sideNode(side, x.node, n);
      const ii = bi + inward[side][0], jj = bj + inward[side][1];
      inflow[idx(ii, jj)] += x.area;
      if (x.minor) continue; // a minor stream finds its own way
      // the parent's path from here: follow its drainage links while inside the tile
      const path = [[ii, jj]];
      let cell = rec[x.links[0]], link = x.links[0], end = -1;
      for (let guard = 0; guard < 64 && cell >= 0; guard++) {
        const at = locateInTile(addr, [Pm[3 * cell], Pm[3 * cell + 1], Pm[3 * cell + 2]], n, layout);
        if (!at) { end = exits.get(link) ?? -1; break; }
        path.push([Math.min(n - 1, Math.max(1, at[0])), Math.min(n - 1, Math.max(1, at[1]))]);
        if (E[cell] <= 0) break;
        link = cell;
        cell = rec[cell];
      }
      if (end < 0 && exits.has(link)) end = exits.get(link);
      if (end >= 0) path.push([end % nx, Math.floor(end / nx)]);
      carve(z, flags, nx, path, x.z, dx);
    }
  }
  // sea is only sea if it reaches the open sea: below-sea-level nodes connected to a sea
  // node of the border (or to the halo's edge), or a body larger than 2% of the tile.
  // Anything else is a closed depression (it fills as a lake), not a sink in the land.
  {
    const label = new Int32Array(N).fill(-1);
    const stack = [];
    for (let c0 = 0; c0 < N; c0++) {
      if (label[c0] >= 0 || !(z[c0] <= 0)) continue;
      const members = [];
      let open = false;
      label[c0] = c0; stack.push(c0);
      while (stack.length) {
        const c = stack.pop();
        members.push(c);
        const ci = c % nx, cj = (c - ci) / nx;
        if (flags[c] & FLAG_FIXED) open = true;
        for (const [di, dj] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
          const a2 = ci + di, b2 = cj + dj;
          if (a2 < 0 || b2 < 0 || a2 >= nx || b2 >= ny) continue;
          const q = b2 * nx + a2;
          if (label[q] < 0 && z[q] <= 0) { label[q] = c0; stack.push(q); }
        }
      }
      const sea = open || members.length > 0.02 * N;
      for (const c of members) {
        if (flags[c] & FLAG_FIXED) continue;
        if (sea) flags[c] |= FLAG_OUTLET; else flags[c] &= ~FLAG_OUTLET;
      }
    }
  }
  // Rho8 (Fairfield & Leymarie 1991, Water Resour. Res. 27): a random weight per node
  // on diagonal slopes when choosing receivers, which removes the straight parallel
  // channels that plain steepest descent draws at 0 and 45 degrees
  const rho = new Float32Array(N);
  for (let c = 0; c < N; c++) rho[c] = Math.fround(2 - uniform(key, c, 0, 0x72686f38));
  // each block's target bias (see blockBias): targets are the smooth field's means plus this
  const bias = new Map();
  for (const bk of new Set(block)) bias.set(bk, smp.blockBias(bk));
  // the solvers' runoff weight per node: its own rain plus any river flowing in
  const runoff = Float32Array.from(P, (p, c) => p + inflow[c] / (dx * dx));
  return { addr, key, n, halo: h, nx, ny, dx, z, U, K, P, inflow, runoff, rho, tempSL, locK: 3, bias, flags, material, pos, block, locIds, locW, base, border, options: o };
}

/**
 * Lower the nodes along a polyline (fractional node coordinates) to a profile falling
 * evenly from zStart to the elevation of its last node. Only lowers; fixed nodes stay.
 */
function carve(z, flags, nx, path, zStart, dx) {
  const nodes = [];
  const seen = new Set();
  for (let s = 0; s < path.length - 1; s++) {
    const [x0, y0] = path[s], [x1, y1] = path[s + 1];
    const steps = Math.max(1, Math.ceil(2 * Math.hypot(x1 - x0, y1 - y0)));
    for (let k = 0; k <= steps; k++) {
      const c = Math.round(y0 + ((y1 - y0) * k) / steps) * nx + Math.round(x0 + ((x1 - x0) * k) / steps);
      if (!seen.has(c)) { seen.add(c); nodes.push(c); }
    }
  }
  if (nodes.length < 2) return;
  const zEnd = z[nodes[nodes.length - 1]];
  const top = Math.max(zStart, zEnd + 1e-4 * dx * nodes.length);
  for (let k = 0; k < nodes.length; k++) {
    const c = nodes[k];
    if (flags[c] & FLAG_FIXED) continue;
    const zz = Math.fround(top + ((zEnd - top) * k) / (nodes.length - 1));
    if (zz < z[c]) z[c] = zz;
  }
}

/**
 * Restriction R: the mean of a field over the cells of each planet cell (block).
 * @returns {Map<number, {sum:number, count:number}>}
 */
export function restrictToBlocks(field, block, mask) {
  const out = new Map();
  for (let c = 0; c < field.length; c++) {
    if (mask && !mask[c]) continue;
    let e = out.get(block[c]);
    if (!e) out.set(block[c], (e = { sum: 0, count: 0 }));
    e.sum += field[c]; e.count++;
  }
  return out;
}

/**
 * Bring a child field's block means back to the parent's (design principle 1): the
 * residual per block is spread back by the same barycentric interpolation that built
 * the child, a few times over (each pass removes most of what is left). Cells with
 * keep[c] set (the sea) are not moved. With `sampledOnly`, blocks with fewer than
 * minCells free cells, or less than half of them free, are left alone (a narrow window
 * sees only slivers of them, and a window with pinned borders can move only a few of
 * their cells: correcting those is ill-conditioned). `weight` (0-1 per cell) scales the
 * correction, e.g. to taper it off towards cells that are held. Returns the largest remaining |residual| (m) over
 * blocks with at least minCells cells.
 */
export function correctBlockMeans(field, input, keep, passes = 12, minCells = 50, sampledOnly = false, weight = null) {
  const { block, locIds, locW, base } = input;
  const L = input.locK || 3;
  const N = field.length;
  // blocks as dense indices (Map look-ups in the inner loop cost more than the arithmetic)
  const ids = [...new Set(block)];
  const dense = new Map(ids.map((b, k) => [b, k]));
  const B = ids.length;
  const blk = Int32Array.from(block, (b) => dense.get(b));
  const loc = Int32Array.from(locIds, (b) => (dense.has(b) ? dense.get(b) : -1));
  const target = new Float64Array(B), count = new Float64Array(B), sum = new Float64Array(B), resid = new Float64Array(B);
  const free = new Float64Array(B);
  for (let c = 0; c < N; c++) { target[blk[c]] += base[c]; count[blk[c]]++; if (!(keep && keep[c])) free[blk[c]]++; }
  // with sampledOnly, a block is corrected (and counted) only if enough of it is free to move
  const sampled = (k) => free[k] >= minCells && free[k] >= 0.3 * count[k];
  for (let k = 0; k < B; k++) target[k] = target[k] / count[k] + (input.bias?.get(ids[k]) || 0);
  let worst = 0;
  for (let pass = 0; pass <= passes; pass++) {
    sum.fill(0);
    for (let c = 0; c < N; c++) sum[blk[c]] += field[c];
    worst = 0;
    for (let k = 0; k < B; k++) {
      resid[k] = target[k] - sum[k] / count[k];
      if (sampledOnly && !sampled(k)) resid[k] = 0;
      else if (count[k] >= minCells) worst = Math.max(worst, Math.abs(resid[k]));
    }
    if (pass === passes) break;
    for (let c = 0; c < N; c++) {
      if (keep && keep[c]) continue;
      let d = 0;
      for (let t = 0; t < L; t++) { const q = loc[L * c + t]; if (q >= 0) d += locW[L * c + t] * resid[q]; }
      field[c] += weight ? weight[c] * d : d;
    }
  }
  return worst;
}

/** Barycentric prolongation of per-block values to every cell. */
export function prolongBlocks(values, input, out = new Float32Array(input.nx * input.ny)) {
  const { locIds, locW } = input;
  const L = input.locK || 3;
  for (let c = 0; c < out.length; c++) {
    let d = 0;
    for (let t = 0; t < L; t++) d += locW[L * c + t] * (values.get(locIds[L * c + t]) || 0);
    out[c] = d;
  }
  return out;
}

/**
 * Run a region tile on any solver with the interface
 *   run(count, {dt, m, D, eps}), readZ() -> Float32Array, setK(Float32Array), readResult()
 * (sync or async).
 *
 * Consistency is steered while the tile erodes, not imposed afterwards: every
 * `nudgeEvery` iterations the block means of the surface are compared with the
 * parent's, and each block's erodibility is scaled (proportional-integral control on
 * log K, spread across blocks by barycentric interpolation): a child standing too high
 * gets more erodible rock, one standing too low less. Uplift stays as the planet gave
 * it, and scaling K within a block leaves the slope-area exponent alone. The planet's
 * coarse denudation and the region's stream power need not agree on how high a
 * landscape must stand to erode at the uplift rate; this reconciles them per block.
 * Factors stay within [1/kRange, kRange]. A final small correction (correctBlockMeans)
 * removes what is left.
 */
export async function runRegion(input, solver) {
  const o = input.options;
  const dt = o.duration / o.iterations;
  const tau = o.nudgeTime * o.duration;
  // the sea does not erode: compare the land cells of each block with the parent's
  // values over those same cells
  const land = new Uint8Array(input.nx * input.ny);
  for (let c = 0; c < land.length; c++) if (!((input.flags[c] & FLAG_OUTLET) && !(input.flags[c] & FLAG_FIXED))) land[c] = 1;
  const target = restrictToBlocks(input.base, input.block, land);
  const K = Float32Array.from(input.K);
  const integral = new Map();
  const lnMax = Math.log(o.kRange);
  for (let done = 0; done < o.iterations;) {
    const count = Math.min(o.nudgeEvery, o.iterations - done);
    await solver.run(count, { dt, m: o.m, D: o.D, eps: o.eps });
    done += count;
    if (done >= o.iterations) break;
    const z = await solver.readZ();
    const now = restrictToBlocks(z, input.block, land);
    const lnK = new Map();
    for (const [b, e] of now) {
      // err > 0: the child stands too high, so its rock should erode faster
      const err = (e.sum / e.count - target.get(b).sum / target.get(b).count - (input.bias?.get(b) || 0)) / o.kHeight;
      const I = Math.max(-lnMax, Math.min(lnMax, (integral.get(b) || 0) + (err * count * dt) / tau));
      integral.set(b, I);
      lnK.set(b, Math.max(-lnMax, Math.min(lnMax, err + I)));
    }
    const f = prolongBlocks(lnK, input);
    for (let c = 0; c < K.length; c++) K[c] = input.K[c] * Math.exp(f[c]);
    await solver.setK(K);
  }
  const result = finishRegion(input, await solver.readResult());
  result.kFactor = [...integral.values()].map(Math.exp);
  return result;
}

/** The CPU reference solver behind the runRegion interface. */
export function cpuSolver(input) {
  const s = new JacobiSolver(input.nx, input.ny, input.dx, { ...input, P: input.runoff });
  return {
    run: (count, p) => { s.run(count, p); },
    readZ: () => s.z,
    setK: (K) => { s.K.set(K); },
    readResult: () => ({ z: s.z, w: s.z, rec: s.rec, A: s.A, iterations: s.iterations }),
    solver: s,
  };
}

/** Run a region tile on the CPU reference of the parallel solver. */
export function runRegionCPU(input) {
  return runRegion(input, cpuSolver(input));
}

const OFF = [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, 1], [1, -1], [-1, -1]];

/**
 * Package solver output: drainage raster for the whole halo-extended grid, the river
 * graph and statistics for the tile proper (halo excluded), and the tile's fields.
 * @param {object} raw {z, w, rec (direction codes, -1 none), A}
 */
export function finishRegion(input, raw) {
  const { nx, ny, dx, halo: h, n } = input;
  const N = nx * ny;
  // consistency with the parent: block means back to the parent's; the fixed border
  // (shared with the neighbours) is left alone
  const elevation = Float32Array.from(raw.w);
  const keep = new Uint8Array(N);
  for (let c = 0; c < N; c++) if (input.flags[c] & FLAG_FIXED) keep[c] = 1;
  for (let c = 0; c < N; c++) if (input.flags[c] & FLAG_FIXED) elevation[c] = input.z[c];
  const consistency = correctBlockMeans(elevation, input, keep, input.options.finalPasses ?? 40);
  // the sea floor takes part in the correction but stays under the sea
  for (let c = 0; c < N; c++) if ((input.flags[c] & FLAG_OUTLET) && !(input.flags[c] & FLAG_FIXED) && elevation[c] > -0.5) elevation[c] = -0.5;
  // the drainage of the final surface, routed exactly on the CPU (priority flood where
  // needed): the extracted network does not depend on the solver's last-iteration
  // transients or on GPU rounding, and every node drains
  const rs = new RasterSolver(nx, ny, dx);
  const outlet = Uint8Array.from(input.flags, (f) => (f & FLAG_OUTLET ? 1 : 0));
  const fixed = Uint8Array.from(input.flags, (f) => (f & FLAG_FIXED ? 1 : 0));
  rs.route(Float64Array.from(elevation), outlet, fixed, input.rho);
  const area = Float32Array.from(rs.accumulate(Float64Array.from(input.P), Float64Array.from(input.inflow)));
  const receiver = new Int32Array(N), dist = new Float32Array(N);
  for (let c = 0; c < N; c++) {
    receiver[c] = rs.receiver[c] === c ? -1 : rs.receiver[c];
    dist[c] = rs.recDist[c] || dx;
  }
  // statistics describe the network this tile made: the tile proper, without reaches
  // that mostly carry water from beyond its edges (their profiles and lengths are set
  // by the parent's boundary condition and by land outside the tile)
  const inherited = rs.accumulate(new Float64Array(N), Float64Array.from(input.inflow));
  const mask = new Uint8Array(N);
  for (let cj = h; cj <= h + n; cj++) for (let ci = h; ci <= h + n; ci++) {
    const c = cj * nx + ci;
    if (inherited[c] <= 0.5 * area[c]) mask[c] = 1;
  }
  const d = { receiver, dist, area, z: elevation, nx, dx, mask };
  const o = input.options;
  const stats = drainageStatistics(d, { channelArea: o.channelArea, flatSlope: 3 * o.eps });
  const rivers = extractRiverGraph(d, { channelArea: o.channelArea });
  // the tile proper ((n + 1) x (n + 1) nodes), without the halo
  const m = n + 1;
  const crop = (arr, Type = Float32Array) => {
    const out = new Type(m * m);
    for (let cj = 0; cj < m; cj++) for (let ci = 0; ci < m; ci++) out[cj * m + ci] = arr[(cj + h) * nx + ci + h];
    return out;
  };
  // the crossings, with the drainage area the child actually delivers or receives there
  if (input.border) {
    for (const side of SIDES) for (const x of input.border[side].crossings) {
      const [ci, cj] = sideNode(side, x.node, n);
      const c = (cj + h) * nx + ci + h;
      rivers.crossings.push({ edge: input.border[side].key, side, node: x.node, inflow: x.inflow, parentArea: x.area, childArea: x.inflow ? x.area : area[c] });
    }
  }
  return {
    addr: input.addr,
    n, nodes: m, dx,
    border: input.border,
    fields: {
      elevation: crop(elevation),
      uplift: crop(input.U),
      precip: crop(Float32Array.from(input.P, (p) => 1000 * p)),
      material: crop(input.material, Uint8Array),
      drainageArea: crop(area),
      temperature: crop(Float32Array.from(elevation, (zz, c) => input.tempSL[c] - LAPSE * Math.max(0, zz))),
    },
    rivers,
    stats: { ...stats, consistency },
    drainage: d,
    iterations: raw.iterations,
  };
}
