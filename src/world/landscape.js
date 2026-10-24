// Level 2: a landscape tile (~9.8 km, 38-m nodes) refined from its level-1 parent.
//
// Initial condition: the parent's elevation interpolated bicubically (Catmull-Rom) to
// each node, plus subgrid relief as seeded noise in sphere coordinates, scaled by how
// much relief the parent's local slope implies within one of its cells. Rock comes from
// the parent's material class as a stratigraphy (below): what is exposed at a node, and
// so how fast it erodes and how steep it can stand, depends on its current elevation.
//
// Processes over `duration` years, in implicit steps (CPU; fastscape.js):
//   rivers       stream power with the exposed rock's erodibility
//   creep        linear diffusion
//   landslides   threshold failure, a sandpile automaton: a slope steeper than the
//                exposed rock's critical slope topples half its excess onto the node
//                below until stable (Bak, Tang & Wiesenfeld 1987 in spirit; the critical
//                slopes are angles of repose and of rock-mass strength)
//   fans         transport-limited deposition: sediment eroded upstream travels down the
//                drainage tree and is laid down where the carrying capacity (Kc A S) drops,
//                at mountain fronts and valley mouths
//   glaciers     during a glaciation in the event log (if recent enough to show), nodes
//                colder than the threshold carry ice: erosion by sliding ice scales with
//                ice discharge and reaches across the whole ice-covered width, so valleys
//                become troughs and valley heads cirques
//
// Seams: edge and corner records as at level 1 (seams.js), from this tile's parent; an
// edge lying on the parent's own seam takes the parent seam's record, so tiles under
// different parents agree without reading each other's parents.

import { nodePosition, tileRect, parentOf, tileKey, vectorToFace, WORLD_LAYOUT } from './cubesphere.js';
import { worldKey, deriveKey, uniform } from './hash.js';
import { fbm3k } from './noise.js';
import { RasterSolver } from './fastscape.js';
import { seamContext, tileBorder, sideNode, SIDES } from './seams.js';
import { correctBlockMeans } from './region.js';
import { drainageStatistics, extractRiverGraph } from './rivers.js';

const LAPSE = 0.0065;

/**
 * Stratigraphy by material class (oceanic, craton, terrane, arc, basin): two rock types
 * alternating in beds, the first `hardShare` of each bed pair hard. Erodibility is
 * relative to the landscape's K; critical slope is tan(angle) a face of that rock holds.
 * `dip` tilts the beds (as a gradient), `thickness` is one bed pair (m).
 */
export const STRATA = [
  { name: 'sea floor', hard: { k: 1, sc: 0.6 }, soft: { k: 1, sc: 0.6 }, hardShare: 0.5, thickness: 50, dip: 0 },
  { name: 'granite', hard: { k: 0.25, sc: 1.6 }, soft: { k: 0.6, sc: 1.0 }, hardShare: 0.85, thickness: 120, dip: 0.05 },
  { name: 'schist and quartzite', hard: { k: 0.4, sc: 1.3 }, soft: { k: 1.8, sc: 0.7 }, hardShare: 0.4, thickness: 24, dip: 1.5 },
  { name: 'lava flows and tuff', hard: { k: 0.3, sc: 1.4 }, soft: { k: 3, sc: 0.55 }, hardShare: 0.45, thickness: 30, dip: 0.04 },
  { name: 'sandstone and shale', hard: { k: 0.45, sc: 1.1 }, soft: { k: 2.5, sc: 0.45 }, hardShare: 0.35, thickness: 40, dip: 0.03 },
];

export const LANDSCAPE_DEFAULTS = {
  duration: 1e5,                 // yr
  steps: 25,
  K: 1e-5,                       // stream-power erodibility of a K=1 rock, m^(1-2m)/yr (m = 0.5)
  m: 0.5,
  D: 0.005,                      // creep, m^2/yr
  relief: { base: 2, perSlope: 0.35, wavelength: 1200, octaves: 5, gain: 0.55 },
  landslideSweeps: 12,
  fans: { capacity: 2e-3, deposit: 0.5 },   // Kc (m^3/yr per m^2 of area per unit slope); share of the excess laid down per node
  glacier: { threshold: -2, K: 2e-4, widen: 1.5 }, // deg C mean annual during the glaciation; ice erodibility; ice spreads this many nodes per sqrt(km^2) of area
  channelArea: 2e5,              // m^2 for the river graph and statistics
  finalPasses: 40,
};

/** Catmull-Rom weights for fractional offset t in [0, 1). */
function crWeights(t, w) {
  const t2 = t * t, t3 = t2 * t;
  w[0] = -0.5 * t3 + t2 - 0.5 * t;
  w[1] = 1.5 * t3 - 2.5 * t2 + 1;
  w[2] = -1.5 * t3 + 2 * t2 + 0.5 * t;
  w[3] = 0.5 * t3 - 0.5 * t2;
  return w;
}

/** Bicubic (Catmull-Rom) interpolation of a node raster m x m at fractional (x, y), clamped. */
function bicubic(f, m, x, y) {
  const xi = Math.floor(x), yi = Math.floor(y);
  const wx = crWeights(x - xi, [0, 0, 0, 0]), wy = crWeights(y - yi, [0, 0, 0, 0]);
  let s = 0;
  for (let b = 0; b < 4; b++) {
    const yy = Math.min(m - 1, Math.max(0, yi - 1 + b));
    let r = 0;
    for (let a = 0; a < 4; a++) r += wx[a] * f[yy * m + Math.min(m - 1, Math.max(0, xi - 1 + a))];
    s += wy[b] * r;
  }
  return s;
}

function bilinear(f, m, x, y) {
  const xi = Math.min(m - 2, Math.max(0, Math.floor(x))), yi = Math.min(m - 2, Math.max(0, Math.floor(y)));
  const u = Math.min(1, Math.max(0, x - xi)), v = Math.min(1, Math.max(0, y - yi));
  return (1 - u) * (1 - v) * f[yi * m + xi] + u * (1 - v) * f[yi * m + xi + 1] + (1 - u) * v * f[(yi + 1) * m + xi] + u * v * f[(yi + 1) * m + xi + 1];
}

/** 1-D Catmull-Rom along an array at fractional index t, clamped. */
function cubic1(a, t) {
  const i = Math.floor(t), w = crWeights(t - i, [0, 0, 0, 0]);
  let s = 0;
  for (let k = 0; k < 4; k++) s += w[k] * a[Math.min(a.length - 1, Math.max(0, i - 1 + k))];
  return s;
}

/**
 * What the parent says at a level-2 position. Positions are given as fractional node
 * coordinates of the parent tile (exact: both tiles use the same face coordinates).
 */
export function landscapeSampler(parent, worldSeed, options = {}) {
  const o = { ...LANDSCAPE_DEFAULTS, ...options, relief: { ...LANDSCAPE_DEFAULTS.relief, ...(options.relief || {}) } };
  const m = parent.nodes, F = parent.fields, R = options.radius ?? 6.371e6;
  const noiseKey = deriveKey(worldKey(worldSeed), 2, 0x6e6f6973, o.noiseSalt || 0);
  const dx1 = parent.dx;
  const noiseAt = (v) => fbm3k(noiseKey, v[0] * R, v[1] * R, v[2] * R, o.relief.wavelength, o.relief.octaves, o.relief.gain);
  /** Elevation from a smooth value e, the parent's relief across one of its cells (m), and the noise at v. */
  const withNoise = (e, cellRelief, v) => {
    const amp = o.relief.base + o.relief.perSlope * cellRelief;
    const fade = Math.min(1, 0.5 + Math.abs(e) / (2 * amp));
    return { e, amp, z: Math.fround(e + amp * noiseAt(v) * fade) };
  };
  function sample(x1, y1, v) {
    const e = bicubic(F.elevation, m, x1, y1);
    // relief within one parent cell, from the parent's local slope
    // (bicubic, whose gradient is continuous: a bilinear one jumps at every parent cell
    // edge, and the texture would change in squares)
    const gx = (bicubic(F.elevation, m, x1 + 0.5, y1) - bicubic(F.elevation, m, x1 - 0.5, y1));
    const gy = (bicubic(F.elevation, m, x1, y1 + 0.5) - bicubic(F.elevation, m, x1, y1 - 0.5));
    return withNoise(e, Math.hypot(gx, gy), v);
  }
  return { sample, withNoise, options: o, dx1, m };
}

/** Fractional parent-node coordinates of level-2 node (ci, cj) of addr. */
function toParent(addr, n2, ci, cj, layout) {
  const b = layout.branching, n1 = layout.tileCells;
  const offI = addr.i - b * Math.floor(addr.i / b), offJ = addr.j - b * Math.floor(addr.j / b);
  const per = n1 / b; // parent cells per child tile
  return [offI * per + (ci * per) / n2, offJ * per + (cj * per) / n2];
}

/**
 * Seam context for the level-2 tiles under one level-1 parent.
 * @param {object} parent level-1 tile result (finishRegion)
 * @param {object} parentSeams the level-1 seam context (regionSeams)
 */
export function landscapeSeams(parent, parentSeams, worldSeed, options = {}, layout = WORLD_LAYOUT) {
  const smp = landscapeSampler(parent, worldSeed, { ...options, radius: parentSeams.radius });
  const n1 = layout.tileCells, m = parent.nodes, P1 = parent.addr;
  const b = layout.branching, per = n1 / b;
  // drainage links of the parent (node to receiver), as unit vectors
  const rec = parent.drainage.receiver, area = parent.drainage.area;
  const links = [];
  for (let c = 0; c < rec.length; c++) {
    const r = rec[c];
    if (r < 0 || area[c] < 4 * parent.dx * parent.dx) continue;
    const ci = c % m, cj = (c - ci) / m, ri = r % m, rj = (r - ri) / m;
    links.push({ id: c, from: nodePosition(P1, n1, ci, cj, undefined, layout), to: nodePosition(P1, n1, ri, rj, undefined, layout), area: area[c] });
  }
  // which side of the parent an edge of child C lies on, if any
  const onParentSide = (C, sC) => {
    const offI = C.i - b * Math.floor(C.i / b), offJ = C.j - b * Math.floor(C.j / b);
    if (sC === 'W' && offI === 0) return { side: 'W', k0: offJ * per };
    if (sC === 'E' && offI === b - 1) return { side: 'E', k0: offJ * per };
    if (sC === 'S' && offJ === 0) return { side: 'S', k0: offI * per };
    if (sC === 'N' && offJ === b - 1) return { side: 'N', k0: offI * per };
    return null;
  };
  const n = layout.tileCells;
  const ctx = seamContext({
    level: 2, worldSeed, layout, radius: parentSeams.radius, linkSpacing: parent.dx, links,
    // elevation / base are only asked for points inside the parent (seams on the
    // parent's border come from the overrides below)
    elevation: (v) => { throw new Error(`level-2 seam sampled at a parent border point ${v}`); },
    base: () => 0,
    notchDepth: () => 0,
    edgeOverride(C, sC, pos) {
      const on = onParentSide(C, sC);
      const z = new Float32Array(n + 1);
      if (on) {
        // the parent's own seam: only its shared record and this level's noise (never
        // either parent's interior, so tiles under both parents compute the same)
        const pb = tileBorder(parentSeams, parentOf(C, layout))[on.side];
        const xs = pb.crossings.filter((x) => x.node >= on.k0 && x.node <= on.k0 + per);
        for (let t = 0; t <= n; t++) {
          z[t] = borderValue(pb.z, on.k0 + (t * per) / n, pos[t]);
          const near = xs.some((x) => Math.abs((x.node - on.k0) * (n / per) - t) < n / per / 2);
          if (near) z[t] = Math.fround(cubic1(pb.z, on.k0 + (t * per) / n));
        }
        const crossings = xs.map((x) => ({ node: Math.round((x.node - on.k0) * (n / per)), area: x.area, inflowToCanonical: x.inflow, links: [], z: x.z }));
        for (const x of crossings) z[x.node] = x.z;
        return { z, crossings };
      }
      if (tileKey(parentOf(C, layout)) !== tileKey(P1)) throw new Error('level-2 seam context used for an edge under another parent');
      return null; // inside the parent: default path (parent surface, parent drainage links)
    },
    cornerOverride(tile, corner, pos) {
      const [ci, cj] = { SW: [0, 0], SE: [n, 0], NW: [0, n], NE: [n, n] }[corner];
      const [x1, y1] = toParent(tile, n, ci, cj, layout);
      // a corner on its parent's border takes the border record's value there
      const onB = x1 === 0 ? ['W', y1] : x1 === n1 ? ['E', y1] : y1 === 0 ? ['S', x1] : y1 === n1 ? ['N', x1] : null;
      if (onB) return borderValue(tileBorder(parentSeams, parentOf(tile, layout))[onB[0]].z, onB[1], pos);
      if (tileKey(parentOf(tile, layout)) !== tileKey(P1)) throw new Error(`level-2 corner under another parent: ${tileKey(tile)} ${corner} at parent node ${x1},${y1}; this parent ${tileKey(P1)}`);
      return smp.sample(x1, y1, pos).z;
    },
  });
  /** A point on a parent's border: the record interpolated along it, plus noise scaled by the record's own slope. */
  function borderValue(arr, k, v) {
    const e = cubic1(arr, k);
    const rel = Math.abs(cubic1(arr, k + 0.5) - cubic1(arr, k - 0.5));
    return smp.withNoise(e, rel, v).z;
  }
  // for crossings on edges inside the parent the default path needs base/notch at a
  // unit vector: invert the face coordinates of that point into the parent's node frame
  const r1 = tileRect(P1, layout);
  function sampleAtVector(v) {
    const g = faceCoords(P1.face, v);
    const x1 = ((g[0] - r1.a0) / r1.size) * n1, y1 = ((g[1] - r1.b0) / r1.size) * n1;
    return smp.sample(x1, y1, v);
  }
  ctx.base = (v) => sampleAtVector(v).e;
  ctx.notchDepth = (v) => 0.5 * sampleAtVector(v).amp;
  ctx.elevation = (v) => sampleAtVector(v).z;
  ctx.sampler = smp;
  return ctx;
}

function faceCoords(face, v) {
  const f = vectorToFace(v[0], v[1], v[2]);
  return [f.a, f.b];
}

/** Stratum at a node: 1 if the exposed rock is the hard bed, else 0. */
function hardAt(strata, zeta, z) {
  const ph = (z - zeta) / strata.thickness;
  return ph - Math.floor(ph) < strata.hardShare ? 1 : 0;
}

/**
 * Inputs of a level-2 tile.
 * @param {object} parent level-1 result (finishRegion) of parentOf(addr)
 * @param {object} ctx landscapeSeams(parent, ...)
 * @param {object[]} events the planet's event log
 */
export function landscapeInput(parent, ctx, addr, events, options = {}, layout = WORLD_LAYOUT) {
  const smp = landscapeSampler(parent, ctx.worldSeed, { ...options, radius: ctx.radius });
  const o = smp.options;
  const n = layout.tileCells, nx = n + 1, N = nx * nx, m = parent.nodes;
  const F = parent.fields;
  const key = deriveKey(worldKey(ctx.worldSeed), addr.level, addr.face, addr.i, addr.j);
  const a = nodePosition(addr, n, n / 2, n / 2, undefined, layout), b2 = nodePosition(addr, n, n / 2 + 1, n / 2, undefined, layout);
  const dx = ctx.radius * Math.acos(Math.min(1, a[0] * b2[0] + a[1] * b2[1] + a[2] * b2[2]));
  const z = new Float32Array(N), base = new Float32Array(N), U = new Float32Array(N), P = new Float32Array(N);
  const tempSL = new Float32Array(N), material = new Uint8Array(N), zeta = new Float32Array(N);
  // the parent's view of each node: its block (nearest parent node) and the weights that
  // spread a per-block correction back to it: bicubic (16 nodes), so corrections have a
  // continuous gradient (bilinear ones shade flats in squares, one per parent cell)
  const block = new Int32Array(N), locIds = new Int32Array(16 * N), locW = new Float32Array(16 * N);
  const wx = [0, 0, 0, 0], wy = [0, 0, 0, 0];
  const flags = new Uint8Array(N); // 1 outlet, 2 fixed
  const v = [0, 0, 0];
  // the beds' structural surface: a gentle global dip plane (direction from the world
  // seed) with long-wavelength undulation, the same for every tile
  const dipKey = deriveKey(worldKey(ctx.worldSeed), 2, 0x64697021);
  const contactKey = deriveKey(worldKey(ctx.worldSeed), 2, 0x636f6e74);
  const dd = [uniform(dipKey, 0) - 0.5, uniform(dipKey, 1) - 0.5, uniform(dipKey, 2) - 0.5];
  const R = ctx.radius;
  for (let cj = 0; cj <= n; cj++) {
    for (let ci = 0; ci <= n; ci++) {
      const c = cj * nx + ci;
      nodePosition(addr, n, ci, cj, v, layout);
      const [x1, y1] = toParent(addr, n, ci, cj, layout);
      const s = smp.sample(x1, y1, v);
      z[c] = s.z; base[c] = s.e;
      U[c] = bilinear(F.uplift, m, x1, y1);
      P[c] = Math.max(0.05, bilinear(F.precip, m, x1, y1) / 1000);
      tempSL[c] = bilinear(F.temperature, m, x1, y1) + LAPSE * Math.max(0, bilinear(F.elevation, m, x1, y1));
      const nxi = Math.min(m - 1, Math.max(0, Math.round(x1))), nyi = Math.min(m - 1, Math.max(0, Math.round(y1)));
      block[c] = nyi * m + nxi;
      const xi = Math.min(m - 2, Math.floor(x1)), yi = Math.min(m - 2, Math.floor(y1)), u = x1 - xi, w = y1 - yi;
      const ids = [yi * m + xi, yi * m + xi + 1, (yi + 1) * m + xi, (yi + 1) * m + xi + 1];
      const wts = [(1 - u) * (1 - w), u * (1 - w), (1 - u) * w, u * w];
      {
        const fx = Math.floor(x1), fy = Math.floor(y1);
        crWeights(x1 - fx, wx); crWeights(y1 - fy, wy);
        for (let bb = 0; bb < 4; bb++) for (let aa = 0; aa < 4; aa++) {
          const qx = Math.min(m - 1, Math.max(0, fx - 1 + aa)), qy = Math.min(m - 1, Math.max(0, fy - 1 + bb));
          locIds[16 * c + 4 * bb + aa] = qy * m + qx;
          locW[16 * c + 4 * bb + aa] = wx[aa] * wy[bb];
        }
      }
      // rock class: one of the four surrounding parent nodes', picked with probability
      // equal to its interpolation weight by a smooth noise field, so contacts between
      // rock classes wander instead of following the parent's grid
      const pick = 0.5 + 0.5 * Math.max(-0.999, Math.min(0.999, 1.6 * fbm3k(contactKey, v[0] * R, v[1] * R, v[2] * R, 600, 3)));
      let acc = 0, which = 3;
      for (let t = 0; t < 4; t++) { acc += wts[t]; if (pick < acc) { which = t; break; } }
      material[c] = F.material[ids[which]];
      const st = STRATA[material[c]];
      const p = [v[0] * R, v[1] * R, v[2] * R];
      zeta[c] = st.dip * (p[0] * dd[0] + p[1] * dd[1] + p[2] * dd[2]) + 25 * fbm3k(dipKey, p[0], p[1], p[2], 3000, 3);
      if (z[c] <= 0) flags[c] |= 1;
    }
  }
  // border from the seam records
  const border = tileBorder(ctx, addr);
  for (const side of SIDES) {
    const bd = border[side];
    for (let t = 0; t <= n; t++) {
      const c = sideNode(side, t, n)[1] * nx + sideNode(side, t, n)[0];
      z[c] = bd.z[t];
      flags[c] = 2 | (bd.z[t] <= 0 ? 1 : 0);
    }
  }
  const inflow = new Float32Array(N);
  const inward = { W: [1, 0], E: [-1, 0], S: [0, 1], N: [0, -1] };
  for (const side of SIDES) for (const x of border[side].crossings) {
    const [bi, bj] = sideNode(side, x.node, n);
    if (x.inflow) inflow[(bj + inward[side][1]) * nx + bi + inward[side][0]] += x.area;
    else flags[bj * nx + bi] = 3;
  }
  // was this place glaciated recently enough to still show it?
  const glaciation = events.find((e) => e.kind === 'glaciation' && e.time < (o.glacialMemory ?? 2e5));
  return {
    addr, key, n, nx, ny: nx, dx, z, base, U, P, inflow, tempSL, material, zeta, flags, border,
    block, locIds, locW, locK: 16, glaciation, options: o,
  };
}

/**
 * Run a level-2 tile (CPU).
 * @returns {object} {fields, rivers, stats, drainage, record}
 */
export function runLandscape(input) {
  const { nx, ny, dx, n, flags, U, P, inflow, tempSL, material, zeta } = input;
  const o = input.options;
  const N = nx * ny;
  const z = Float64Array.from(input.z);
  const outlet = Uint8Array.from(flags, (f) => f & 1), fixed = Uint8Array.from(flags, (f) => (f & 2) >> 1);
  const rs = new RasterSolver(nx, ny, dx);
  const dt = o.duration / o.steps;
  const K = new Float64Array(N), Sc = new Float64Array(N), D = new Float64Array(N).fill(o.D);
  const record = { landslide: new Float32Array(N), deposit: new Float32Array(N), ice: new Uint8Array(N), hard: new Uint8Array(N) };
  const gl = input.glaciation;
  const glacialSteps = gl ? Math.round(o.steps * Math.min(1, gl.duration / o.duration)) : 0;
  const before = new Float64Array(N);
  const rec = rs.receiver;
  const held = Uint8Array.from(flags, (f) => (f ? 1 : 0));
  const Qs = new Float64Array(N), reach = new Float32Array(N);
  for (let step = 0; step < o.steps; step++) {
    // exposed rock
    for (let c = 0; c < N; c++) {
      const st = STRATA[material[c]];
      const h = hardAt(st, zeta[c], z[c]);
      const r = h ? st.hard : st.soft;
      K[c] = o.K * r.k;
      Sc[c] = r.sc;
      record.hard[c] = h;
    }
    // glaciers in the first part of the run (the glaciation is in the recent past)
    const icy = step < glacialSteps;
    rs.route(z, outlet, fixed);
    const A = rs.accumulate(P, inflow);
    before.set(z);
    if (icy) {
      const cool = gl.params.cooling;
      const ice = record.ice;
      // ice where it was cold enough, then down the valleys and across their width
      for (let c = 0; c < N; c++) if (tempSL[c] - cool - LAPSE * Math.max(0, z[c]) < o.glacier.threshold) ice[c] = 1;
      for (let k = rs.orderCount - 1; k >= 0; k--) {
        const c = rs.order[k], r = rec[c];
        if (ice[c] && r !== c && A[c] > 1e6) ice[r] = 1; // glacier tongues reach below the snowline
      }
      // the ice spreads sideways by a reach that grows with ice discharge: a max-minus-
      // distance transform in two raster sweeps (chamfer weights 1 and sqrt 2), O(n)
      for (let c = 0; c < N; c++) reach[c] = ice[c] ? Math.min(8, o.glacier.widen * Math.sqrt(A[c] / 1e6)) : -1;
      const s2 = Math.SQRT2;
      for (let cj = 0; cj < ny; cj++) for (let ci = 0; ci < nx; ci++) {
        const c = cj * nx + ci;
        let r = reach[c];
        if (ci > 0) r = Math.max(r, reach[c - 1] - 1);
        if (cj > 0) {
          r = Math.max(r, reach[c - nx] - 1);
          if (ci > 0) r = Math.max(r, reach[c - nx - 1] - s2);
          if (ci < nx - 1) r = Math.max(r, reach[c - nx + 1] - s2);
        }
        reach[c] = r;
      }
      for (let cj = ny - 1; cj >= 0; cj--) for (let ci = nx - 1; ci >= 0; ci--) {
        const c = cj * nx + ci;
        let r = reach[c];
        if (ci < nx - 1) r = Math.max(r, reach[c + 1] - 1);
        if (cj < ny - 1) {
          r = Math.max(r, reach[c + nx] - 1);
          if (ci < nx - 1) r = Math.max(r, reach[c + nx + 1] - s2);
          if (ci > 0) r = Math.max(r, reach[c + nx - 1] - s2);
        }
        reach[c] = r;
      }
      for (let c = 0; c < N; c++) if (reach[c] >= 0) ice[c] = 1;
      for (let c = 0; c < N; c++) if (ice[c]) K[c] = Math.max(K[c], o.glacier.K);
    }
    rs.erode(z, { dt, K, m: o.m, n: 1, U, fixed: held });
    if (icy) {
      // ice reaches across the trough: erosion spreads sideways within the ice
      const ice = record.ice;
      for (let c = 0; c < N; c++) D[c] = ice[c] ? o.D * 400 : o.D;
    } else D.fill(o.D);
    rs.diffuse(z, { dt, D, fixed: held });
    // fans: sediment down the tree, deposited where capacity falls short
    Qs.fill(0);
    for (let k = rs.orderCount - 1; k >= 0; k--) {
      const c = rs.order[k], r = rec[c];
      if (flags[c]) continue;
      Qs[c] += Math.max(0, before[c] - z[c]) * dx * dx;
      const S = r !== c ? Math.max(0, (z[c] - z[r]) / rs.recDist[c]) : 0;
      const cap = o.fans.capacity * A[c] * S * dt;
      if (Qs[c] > cap) {
        const dep = Math.min((Qs[c] - cap) * o.fans.deposit, 0.5 * dx * dx * Math.max(0, (r !== c ? z[r] - z[c] : 0) + 2));
        if (dep > 0) { z[c] += dep / (dx * dx); Qs[c] -= dep; record.deposit[c] += dep / (dx * dx); }
      }
      if (r !== c) Qs[r] += Qs[c];
    }
    // landslides: topple slopes beyond the rock's critical slope. The first sweep visits
    // every node; later ones only the nodes a failure touched (a worklist)
    let work = null;
    for (let sw = 0; sw < o.landslideSweeps; sw++) {
      const next = [];
      const visit = (c) => {
        const ci = c % nx, cj = (c - ci) / nx;
        if (flags[c] || ci === 0 || cj === 0 || ci === nx - 1 || cj === ny - 1) return;
        let best = 0, bj = -1, bd = dx;
        for (let k = 0; k < 8; k++) {
          const q = c + OFFS[k] + OFFS[k + 8] * nx;
          const d = k < 4 ? dx : dx * Math.SQRT2;
          const sl = (z[c] - z[q]) / d;
          if (sl > best) { best = sl; bj = q; bd = d; }
        }
        if (bj < 0 || best <= Sc[c]) return;
        const ex = 0.5 * (z[c] - z[bj] - Sc[c] * bd);
        z[c] -= ex;
        record.landslide[c] += ex;
        if (!(flags[bj] & 1)) { z[bj] += ex; record.deposit[bj] += ex; }
        next.push(c, bj);
        for (let k = 0; k < 8; k++) next.push(c + OFFS[k] + OFFS[k + 8] * nx);
      };
      if (!work) for (let c = 0; c < N; c++) visit(c);
      else for (const c of work) visit(c);
      if (!next.length) break;
      work = [...new Set(next)].sort((p, q) => p - q);
    }
  }
  return finishLandscape(input, z, record);
}

// the 8 neighbour offsets: OFFS[k] in x, OFFS[k + 8] in y (cardinal first)
const OFFS = [1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, 1, -1, -1];
function finishLandscape(input, z64, record) {
  const { nx, ny, dx, flags } = input;
  const N = nx * ny;
  const elevation = Float32Array.from(z64);
  for (let c = 0; c < N; c++) if (flags[c] & 2) elevation[c] = input.z[c];
  const keep = Uint8Array.from(flags, (f) => (f ? 1 : 0));
  const consistency = correctBlockMeans(elevation, input, keep, input.options.finalPasses, 30);
  const rs = new RasterSolver(nx, ny, dx);
  rs.route(Float64Array.from(elevation), Uint8Array.from(flags, (f) => f & 1), Uint8Array.from(flags, (f) => (f & 2) >> 1));
  const area = Float32Array.from(rs.accumulate(Float64Array.from(input.P), Float64Array.from(input.inflow)));
  const receiver = new Int32Array(N), dist = new Float32Array(N);
  for (let c = 0; c < N; c++) { receiver[c] = rs.receiver[c] === c ? -1 : rs.receiver[c]; dist[c] = rs.recDist[c] || dx; }
  const d = { receiver, dist, area, z: elevation, nx, dx };
  const o = input.options;
  const stats = drainageStatistics(d, { channelArea: o.channelArea });
  const rivers = extractRiverGraph(d, { channelArea: o.channelArea });
  const temperature = Float32Array.from(elevation, (zz, c) => input.tempSL[c] - LAPSE * Math.max(0, zz));
  return {
    addr: input.addr, n: input.n, nodes: nx, dx, border: input.border,
    fields: { elevation, material: input.material, drainageArea: area, temperature, hardRock: record.hard },
    record, rivers, stats: { ...stats, consistency }, drainage: d,
  };
}

/** Summary statistics of a landscape tile, for the variety test. */
export function landscapeSummary(r) {
  const e = r.fields.elevation, N = e.length, nx = r.nodes, dx = r.dx;
  const land = [];
  for (let c = 0; c < N; c++) if (e[c] > 0) land.push(e[c]);
  land.sort((a, b) => a - b);
  const q = (p) => land[Math.floor(p * (land.length - 1))] ?? 0;
  let slopeSum = 0, steep = 0, cnt = 0;
  for (let j = 1; j < nx - 1; j++) for (let i = 1; i < nx - 1; i++) {
    const c = j * nx + i;
    const s = Math.hypot(e[c + 1] - e[c - 1], e[c + nx] - e[c - nx]) / (2 * dx);
    slopeSum += s; if (s > 0.7) steep++; cnt++;
  }
  const relief = q(0.95) - q(0.05);
  const hyps = relief > 0 ? (land.reduce((a, b) => a + b, 0) / land.length - q(0.05)) / relief : 0;
  const A = r.fields.drainageArea;
  let chan = 0;
  for (let c = 0; c < N; c++) if (A[c] > 1e6) chan++;
  const frac = (arr, thr) => arr.reduce((a, v) => a + (v > thr ? 1 : 0), 0) / N;
  return {
    relief, meanSlope: slopeSum / cnt, cliffFraction: steep / cnt, hypsometricIntegral: hyps,
    drainageDensity: chan / N, landslideFraction: frac(r.record.landslide, 0.5), depositFraction: frac(r.record.deposit, 0.5),
    iceFraction: frac(r.record.ice, 0), hardRockFraction: frac(r.fields.hardRock, 0),
  };
}

