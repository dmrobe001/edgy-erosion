// Level 1 on the irregular mesh, with simulated seams.
//
// Generation is dimension-ordered, every piece from the planet alone:
//
//   corner patch  a square window around a tile corner (shared by its 3 or 4 tiles),
//                 eroded like a tile. Its core (the middle) becomes fixed for everything
//                 that touches it.
//   edge strip    a window straddling a tile edge, a quarter of a tile deep on each side,
//                 its ends pinned to the corner cores. Its band (the nodes within `band`
//                 cells of the edge) becomes fixed for the two tiles.
//   tile          the tile and a margin, its bands and corner cores pinned to the strips'
//                 and corners' solutions; it erodes everything else.
//
// A pinned node is a sink for the window that pins it (its elevation is given), and
// wherever the shared solution sent water from a pinned node into the window's free
// interior, that water and its sediment enter there. So channels cross a seam inside a
// band that both tiles share, carrying the same water and sediment, instead of meeting
// a wall. Strip and corner windows have an open outer ring at the parent's elevation;
// the planet's rivers enter them where its flow links cross into them, with drainage
// area and sediment flux.
//
// All windows are views of one global point set (pointmesh.js); a node has the same
// id, position, initial elevation and rock in every window that holds it.

import { WORLD_LAYOUT, tileKey, tileRect, neighborOf, cornerTiles, nodePosition } from './cubesphere.js';
import { windowPoints, triangulateWindow, extendedFace, locateInWindow, cellsPerFace } from './pointmesh.js';
import { GraphSolver, FLAG_OUTLET, FLAG_FIXED, CLASSES, GRAPH_DEFAULTS } from './graphsolver.js';
import { regionSampler, restrictToBlocks, correctBlockMeans, prolongBlocks, REGION_DEFAULTS } from './region.js';
import { drainageStatistics, extractRiverGraph, upstreamOrder } from './rivers.js';
import { RasterSolver } from './fastscape.js';
import { SIDES, sideNode } from './seams.js';

export const MESH_REGION_DEFAULTS = {
  ...REGION_DEFAULTS,
  // window sizes in cells; left unset, they scale with the tile (values for 256 cells)
  band: undefined,  // cells either side of an edge pinned to the strip (8)
  strip: undefined, // cells a strip reaches either side of its edge (56)
  corner: undefined, // half-width of a corner patch (56)
  core: undefined,  // half-width of a corner patch's pinned core (20)
  margin: 8,        // cells of extra window for triangulation (connectivity settles within ~5)
  G: GRAPH_DEFAULTS.G,
  // share of eroded rock in each grain class (fine, sand, gravel) by material class:
  // oceanic, craton (granite), terrane (schist), arc (lavas), basin (sediments)
  grains: [[0.5, 0.3, 0.2], [0.3, 0.4, 0.3], [0.4, 0.35, 0.25], [0.45, 0.3, 0.25], [0.6, 0.35, 0.05]],
};

const LAPSE = 0.0065;
const CORNERS = ['SW', 'SE', 'NW', 'NE'];
const CORNER_OF = { W: ['SW', 'NW'], E: ['SE', 'NE'], S: ['SW', 'SE'], N: ['NW', 'NE'] };

/** Everything the mesh level shares under one planet: sampler, flow links, caches. */
export function meshRegionContext(planet, options = {}, layout = WORLD_LAYOUT) {
  const o = { ...MESH_REGION_DEFAULTS, ...options, relief: { ...MESH_REGION_DEFAULTS.relief, ...(options.relief || {}) } };
  const nc = layout.tileCells;
  o.band ??= Math.max(3, Math.round(nc / 32));
  o.strip ??= Math.max(10, Math.round((56 * nc) / 256));
  o.corner ??= o.strip;
  o.core ??= Math.max(o.band + 1, Math.round((20 * nc) / 256));
  const smp = regionSampler(planet, o);
  const P = planet.mesh.pos, F = planet.fields, rec = planet.flow.receiver;
  const links = [];
  for (let i = 0; i < rec.length; i++) {
    const r = rec[i];
    if (r < 0) continue;
    const sed = F.sedimentFlux ? F.sedimentFlux[i] : 0;
    links.push({ from: [P[3 * i], P[3 * i + 1], P[3 * i + 2]], to: [P[3 * r], P[3 * r + 1], P[3 * r + 2]], area: F.drainageArea[i], sed: o.grains[F.material[i]].map((g) => g * sed) });
  }
  return { planet, o, layout, smp, links, corners: new Map(), strips: new Map(), N: cellsPerFace(1, layout), radius: planet.mesh.radius };
}

// --- geometry in a frame tile's extended face coordinates ---------------------------

function cellSize(ctx) { return tileRect({ level: 1, face: 0, i: 0, j: 0 }, ctx.layout).size / ctx.layout.tileCells; }

function cornerAB(tile, name, layout) {
  const r = tileRect(tile, layout);
  return [r.a0 + (name[1] === 'E' ? r.size : 0), r.b0 + (name[0] === 'N' ? r.size : 0)];
}

/** Canonical corner: the lowest-keyed tile sharing it, and that tile's name for it. */
function canonicalCorner(tile, name, layout) {
  const tiles = cornerTiles(tile, name, layout);
  const c0 = tiles[0];
  const me = nodePosition(tile, layout.tileCells, name[1] === 'E' ? layout.tileCells : 0, name[0] === 'N' ? layout.tileCells : 0, undefined, layout);
  let best = null, bd = Infinity;
  for (const nm of CORNERS) {
    const p = nodePosition(c0, layout.tileCells, nm[1] === 'E' ? layout.tileCells : 0, nm[0] === 'N' ? layout.tileCells : 0, undefined, layout);
    const d = (p[0] - me[0]) ** 2 + (p[1] - me[1]) ** 2 + (p[2] - me[2]) ** 2;
    if (d < bd) { bd = d; best = nm; }
  }
  return { key: `C:${tiles.map(tileKey).join('|')}`, tile: c0, name: best };
}

function canonicalEdge(tile, side, layout) {
  const nb = neighborOf(tile, side, layout);
  const [ka, kb] = [tileKey(tile), tileKey(nb)];
  if (ka < kb) return { key: `E:${ka}|${kb}`, tile, side, other: nb };
  for (const s of SIDES) if (tileKey(neighborOf(nb, s, layout)) === ka) return { key: `E:${kb}|${ka}`, tile: nb, side: s, other: tile };
  throw new Error('edge without a way back');
}

// --- windows ---------------------------------------------------------------------

/**
 * Build a window: points, triangulation, the parent's fields at every node.
 * `classify(fa, fb)` returns 'open' (the outer ring: fixed, an outlet at the parent's
 * elevation), or 'free'.
 */
function buildWindow(ctx, frame, rect, classify) {
  const { o, smp } = ctx;
  const W = windowPoints(ctx.planet.worldSeed, 1, frame, rect, ctx.layout);
  const g = triangulateWindow(W.pos, ctx.radius);
  const n = g.n;
  const z = new Float32Array(n), base = new Float32Array(n), U = new Float32Array(n), K = new Float32Array(n), Pw = new Float32Array(n);
  const tempSL = new Float32Array(n), flags = new Uint32Array(n), material = new Uint8Array(n), block = new Int32Array(n);
  const locIds = new Int32Array(3 * n), locW = new Float32Array(3 * n), frac = new Float32Array(CLASSES * n), fab = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    const v = W.pos[i];
    const s = smp.sample(v);
    z[i] = s.z; base[i] = s.e; U[i] = s.u; Pw[i] = Math.max(0.05, s.p / 1000); K[i] = o.K * s.k; tempSL[i] = s.tsl;
    material[i] = ctx.planet.fields.material[s.loc.cell];
    block[i] = s.loc.cell;
    for (let t = 0; t < 3; t++) { locIds[3 * i + t] = s.loc.ids[t]; locW[3 * i + t] = s.loc.weights[t]; }
    const gr = o.grains[material[i]];
    for (let k = 0; k < CLASSES; k++) frac[k * n + i] = gr[k];
    const e = extendedFace(frame.face, v);
    fab[2 * i] = e[0]; fab[2 * i + 1] = e[1];
    if (z[i] <= 0) flags[i] |= FLAG_OUTLET;
    if (classify(e[0], e[1]) === 'open') flags[i] = FLAG_FIXED | FLAG_OUTLET;
  }
  const idIndex = new Map();
  for (let i = 0; i < n; i++) idIndex.set(W.ids[i], i);
  const bias = new Map();
  for (const b of new Set(block)) bias.set(b, smp.blockBias(b));
  return { frame, rect, ids: W.ids, pos: W.pos, g, n, idIndex, z, base, U, K, P: Pw, tempSL, flags, material, block, locIds, locW, locK: 3, bias, frac, fab,
    inflow: new Float32Array(n), sedIn: new Float32Array(CLASSES * n), pinned: new Uint8Array(n) };
}

/** Pin the nodes a shared solution holds in `ids`, and take in what it sends into the window. */
function pinFrom(win, sol, ids) {
  for (const id of ids) {
    const i = win.idIndex.get(id);
    if (i === undefined) continue;
    const s = sol.idIndex.get(id);
    win.z[i] = sol.z[s];
    win.flags[i] = FLAG_FIXED | FLAG_OUTLET;
    win.pinned[i] = 1;
  }
  // hand-off: a pinned node whose receiver (in the shared solution) is free here
  for (const id of ids) {
    const i = win.idIndex.get(id);
    if (i === undefined) continue;
    const s = sol.idIndex.get(id), r = sol.rec[s];
    if (r < 0) continue;
    const ri = win.idIndex.get(sol.ids[r]);
    if (ri === undefined || win.pinned[ri] || (win.flags[ri] & FLAG_FIXED)) continue;
    win.inflow[ri] += sol.A[s];
    for (let k = 0; k < CLASSES; k++) win.sedIn[k * win.n + ri] += sol.Qs[k * sol.n + s];
  }
}

/** The planet's rivers entering a window: where a flow link crosses into it, its area and sediment enter at the nearest free node. */
function parentInflows(ctx, win) {
  const { rect, frame } = win;
  const inside = (v) => { const e = extendedFace(frame.face, v); return e && e[0] >= rect.a0 && e[0] <= rect.a1 && e[1] >= rect.b0 && e[1] <= rect.b1; };
  const free = [];
  for (let i = 0; i < win.n; i++) if (!(win.flags[i] & FLAG_FIXED)) free.push(i);
  if (!free.length) return 0;
  let count = 0;
  for (const L of ctx.links) {
    if (inside(L.from)) continue;
    // first point along the link inside the window
    let hit = null;
    for (let s = 1; s <= 32 && !hit; s++) {
      const t = s / 32, v = [L.from[0] + t * (L.to[0] - L.from[0]), L.from[1] + t * (L.to[1] - L.from[1]), L.from[2] + t * (L.to[2] - L.from[2])];
      const l = Math.hypot(...v);
      const u = [v[0] / l, v[1] / l, v[2] / l];
      if (inside(u)) hit = u;
    }
    if (!hit) continue;
    let best = -1, bd = Infinity;
    for (const i of free) { const p = win.pos[i]; const d = (p[0] - hit[0]) ** 2 + (p[1] - hit[1]) ** 2 + (p[2] - hit[2]) ** 2; if (d < bd) { bd = d; best = i; } }
    win.inflow[best] += L.area;
    for (let k = 0; k < CLASSES; k++) win.sedIn[k * win.n + best] += L.sed[k];
    count++;
  }
  return count;
}

/** The CPU solver behind solveWindow's interface. */
export function cpuGraphSolver(win) {
  const s = new GraphSolver(win.g, { ...win, P: win.P });
  return {
    run: (count, p) => { s.run(count, p); },
    readZ: () => s.z,
    setK: (K) => { s.K.set(K); },
    readResult: () => ({ z: s.z, rec: s.rec, recDist: s.recDist, A: s.A, Qs: s.Qs, dep: s.dep, iterations: s.iterations }),
  };
}

/**
 * Erode a window, steering block means to the parent's (proportional-integral control
 * of erodibility, as in region.js), then correct what is left.
 */
async function solveWindow(ctx, win, makeSolver, { gentle = false } = {}) {
  const o = ctx.o;
  const dt = o.duration / o.iterations, tau = o.nudgeTime * o.duration;
  const n = win.n;
  const solver = await makeSolver(win);
  const land = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (!(win.flags[i] & FLAG_OUTLET)) land[i] = 1;
  const target = restrictToBlocks(win.base, win.block, land);
  const K = Float32Array.from(win.K), integral = new Map(), lnMax = Math.log(o.kRange);
  const params = { dt, m: o.m, D: o.D, eps: o.eps, G: o.G };
  for (let done = 0; done < o.iterations;) {
    const count = Math.min(o.nudgeEvery, o.iterations - done);
    await solver.run(count, params);
    done += count;
    if (done >= o.iterations) break;
    const z = await solver.readZ();
    const now = restrictToBlocks(z, win.block, land);
    const lnK = new Map();
    for (const [b, e] of now) {
      const err = (e.sum / e.count - target.get(b).sum / target.get(b).count - (win.bias.get(b) || 0)) / o.kHeight;
      const I = Math.max(-lnMax, Math.min(lnMax, (integral.get(b) || 0) + (err * count * dt) / tau));
      integral.set(b, I);
      lnK.set(b, Math.max(-lnMax, Math.min(lnMax, err + I)));
    }
    const fk = prolongBlocks(lnK, win, new Float32Array(n));
    for (let i = 0; i < n; i++) K[i] = win.K[i] * Math.exp(fk[i]);
    await solver.setK(K);
  }
  const raw = await solver.readResult();
  solver.destroy?.();
  const z = Float32Array.from(raw.z);
  const keep = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (win.flags[i] & FLAG_FIXED) { keep[i] = 1; z[i] = win.z[i]; }
  // strips and corners see slivers of many blocks: correct only the well-sampled ones, a few times
  const consistency = gentle ? correctBlockMeans(z, win, keep, 4, 400, true) : correctBlockMeans(z, win, keep, o.finalPasses ?? 40, 400, true);
  for (let i = 0; i < n; i++) if ((win.flags[i] & FLAG_OUTLET) && !(win.flags[i] & FLAG_FIXED) && z[i] > -0.5) z[i] = -0.5;
  return { ...win, z, rec: Int32Array.from(raw.rec), recDist: Float32Array.from(raw.recDist), A: Float32Array.from(raw.A), Qs: Float32Array.from(raw.Qs), dep: raw.dep ? Float32Array.from(raw.dep) : null, consistency, iterations: raw.iterations };
}

// --- corners, strips, tiles -------------------------------------------------------

/** The corner patch at corner `name` of `tile` (computed once, in its canonical frame). */
export async function cornerPatch(ctx, tile, name, makeSolver = cpuGraphSolver) {
  const cc = canonicalCorner(tile, name, ctx.layout);
  if (ctx.corners.has(cc.key)) return ctx.corners.get(cc.key);
  const job = (async () => {
    const h = cellSize(ctx), { corner: R, core, margin } = ctx.o;
    const [ac, bc] = cornerAB(cc.tile, cc.name, ctx.layout);
    const rect = { a0: ac - (R + margin) * h, a1: ac + (R + margin) * h, b0: bc - (R + margin) * h, b1: bc + (R + margin) * h };
    const win = buildWindow(ctx, cc.tile, rect, (fa, fb) => (Math.abs(fa - ac) > R * h || Math.abs(fb - bc) > R * h ? 'open' : 'free'));
    parentInflows(ctx, win);
    const sol = await solveWindow(ctx, win, makeSolver, { gentle: true });
    sol.key = cc.key;
    sol.coreIds = [];
    for (let i = 0; i < sol.n; i++) if (Math.abs(sol.fab[2 * i] - ac) <= core * h && Math.abs(sol.fab[2 * i + 1] - bc) <= core * h) sol.coreIds.push(sol.ids[i]);
    // the corner node's values, sampled once at the canonical tile's corner position
    const nn = ctx.layout.tileCells;
    sol.point = pointValues(sol, nodePosition(cc.tile, nn, cc.name[1] === 'E' ? nn : 0, cc.name[0] === 'N' ? nn : 0, undefined, ctx.layout));
    return sol;
  })();
  ctx.corners.set(cc.key, job);
  return job;
}

/** The edge strip on `side` of `tile` (computed once, in its canonical frame). */
export async function edgeStrip(ctx, tile, side, makeSolver = cpuGraphSolver) {
  const ce = canonicalEdge(tile, side, ctx.layout);
  if (ctx.strips.has(ce.key)) return ctx.strips.get(ce.key);
  const job = (async () => {
    const h = cellSize(ctx), { strip: Ws, core, band, margin } = ctx.o;
    const C = ce.tile, sC = ce.side;
    const corners = await Promise.all(CORNER_OF[sC].map((nm) => cornerPatch(ctx, C, nm, makeSolver)));
    const r = tileRect(C, ctx.layout);
    const vertical = sC === 'W' || sC === 'E';
    const line = sC === 'W' ? r.a0 : sC === 'E' ? r.a0 + r.size : sC === 'S' ? r.b0 : r.b0 + r.size;
    const lo = vertical ? r.b0 : r.a0, hi = lo + r.size;
    const across = (fa, fb) => (vertical ? fa : fb) - line, along = (fa, fb) => (vertical ? fb : fa);
    const ext = (core + margin) * h, deep = (Ws + margin) * h;
    const rect = vertical ? { a0: line - deep, a1: line + deep, b0: lo - ext, b1: hi + ext } : { a0: lo - ext, a1: hi + ext, b0: line - deep, b1: line + deep };
    const win = buildWindow(ctx, C, rect, (fa, fb) => (Math.abs(across(fa, fb)) > Ws * h || along(fa, fb) < lo - core * h || along(fa, fb) > hi + core * h ? 'open' : 'free'));
    for (const c of corners) pinFrom(win, c, c.coreIds);
    parentInflows(ctx, win);
    const sol = await solveWindow(ctx, win, makeSolver, { gentle: true });
    sol.key = ce.key;
    sol.canonical = C; sol.side = sC;
    const coreSet = new Set(corners.flatMap((c) => c.coreIds));
    sol.bandIds = [];
    for (let i = 0; i < sol.n; i++) {
      const fa = sol.fab[2 * i], fb = sol.fab[2 * i + 1];
      if (Math.abs(across(fa, fb)) <= band * h && along(fa, fb) >= lo && along(fa, fb) <= hi && !coreSet.has(sol.ids[i])) sol.bandIds.push(sol.ids[i]);
    }
    sol.geometry = { line, lo, hi, vertical };
    // the edge's node values, sampled once at the canonical tile's nodes (both tiles and
    // the level below read these, so they agree bitwise); its ends are the corners'
    const nn = ctx.layout.tileCells;
    sol.record = [corners[0].point];
    for (let t = 1; t < nn; t++) sol.record.push(pointValues(sol, nodePosition(C, nn, ...sideNode(sC, t, nn), undefined, ctx.layout)));
    sol.record.push(corners[1].point);
    return sol;
  })();
  ctx.strips.set(ce.key, job);
  return job;
}

/**
 * A level-1 tile on the mesh.
 * @param {object} ctx meshRegionContext
 * @param {{level:1, face, i, j}} addr
 * @param {(win:object) => object} [makeSolver] solver factory (CPU by default; gpu.js for WebGPU)
 */
export async function meshRegionTile(ctx, addr, makeSolver = cpuGraphSolver) {
  const o = ctx.o, h = cellSize(ctx), L = ctx.layout, n = L.tileCells;
  const strips = {};
  for (const s of SIDES) strips[s] = await edgeStrip(ctx, addr, s, makeSolver);
  const corners = await Promise.all(CORNERS.map((nm) => cornerPatch(ctx, addr, nm, makeSolver)));
  const r = tileRect(addr, L);
  const m = o.band + o.margin;
  const rect = { a0: r.a0 - m * h, a1: r.a0 + r.size + m * h, b0: r.b0 - m * h, b1: r.b0 + r.size + m * h };
  const bandOut = o.band * h;
  const win = buildWindow(ctx, addr, rect, (fa, fb) => (fa < r.a0 - bandOut || fa > r.a0 + r.size + bandOut || fb < r.b0 - bandOut || fb > r.b0 + r.size + bandOut ? 'open' : 'free'));
  for (const c of corners) pinFrom(win, c, c.coreIds);
  for (const s of SIDES) pinFrom(win, strips[s], strips[s].bandIds);
  const sol = await solveWindow(ctx, win, makeSolver);
  const out = finishMeshTile(ctx, addr, sol, strips, corners);
  out.borders = {};
  for (const s of SIDES) out.borders[s] = await meshSideBorder(ctx, addr, s, makeSolver);
  return out;
}

/**
 * The parent-border interface level 2 needs (`sideBorder(addr, side)`, `radius`), from a
 * mesh level-1 tile's own border records: its own sides, and its neighbours' sides
 * facing it (the same shared record, seen from the other tile).
 */
export function meshBorders(tile, radius, layout = WORLD_LAYOUT) {
  const n = layout.tileCells, me = tileKey(tile.addr);
  const d = (p, q) => (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2;
  const corners = [['W', 0], ['W', n], ['E', 0], ['E', n]].map(([sd, t]) => ({ p: nodePosition(tile.addr, n, ...sideNode(sd, t, n), undefined, layout), z: tile.borders[sd].z[t] }));
  return {
    radius,
    /** The value at one of this tile's corners (the one nearest p). */
    cornerZ(p) {
      let best = corners[0];
      for (const c of corners) if (d(c.p, p) < d(best.p, p)) best = c;
      return best.z;
    },
    sideBorder(addr, side) {
      if (tileKey(addr) === me) return tile.borders[side];
      const mine = SIDES.find((s) => tileKey(neighborOf(tile.addr, s, layout)) === tileKey(addr) && tileKey(neighborOf(addr, side, layout)) === me);
      if (!mine) throw new Error(`no border of ${tileKey(addr)} ${side} from ${me}`);
      const b = tile.borders[mine];
      const p0 = nodePosition(addr, n, ...sideNode(side, 0, n), undefined, layout);
      const reversed = d(p0, nodePosition(tile.addr, n, ...sideNode(mine, n, n), undefined, layout)) < d(p0, nodePosition(tile.addr, n, ...sideNode(mine, 0, n), undefined, layout));
      const z = new Float32Array(n + 1);
      for (let t = 0; t <= n; t++) z[t] = b.z[reversed ? n - t : t];
      const crossings = b.crossings.map((x) => ({ ...x, node: reversed ? n - x.node : x.node, inflow: !x.inflow }));
      return { key: b.key, z, crossings };
    },
  };
}

// --- products ----------------------------------------------------------------------

/** Value of a solution's field at a unit vector (barycentric in 3-D on its triangulation). */
function sampleSolution(sol, v, field, hint) {
  const loc = locateInWindow(sol.g, sol.pos, v, hint);
  if (!loc) return null;
  return { value: loc.w[0] * field[loc.ids[0]] + loc.w[1] * field[loc.ids[1]] + loc.w[2] * field[loc.ids[2]], loc };
}

function maxAt(sol, loc, field, stride = 0, k = 0) {
  return Math.max(field[k * stride + loc.ids[0]], field[k * stride + loc.ids[1]], field[k * stride + loc.ids[2]]);
}

/** What a tile node takes from a solution at unit vector v: elevation interpolated, flow and rock from the nodes around it. */
function pointValues(sol, v, hint = -1) {
  const s = sampleSolution(sol, v, sol.z, hint);
  if (!s) throw new Error(`point outside its window (${sol.key || 'tile'})`);
  const near = s.loc.ids[s.loc.w.indexOf(Math.max(...s.loc.w))];
  return {
    z: Math.fround(s.value), A: maxAt(sol, s.loc, sol.A),
    sed: Array.from({ length: CLASSES }, (_, k) => maxAt(sol, s.loc, sol.Qs, sol.n, k)),
    material: sol.material[near], U: sol.U[near], P: sol.P[near], tsl: sol.tempSL[near], tri: s.loc.tri,
  };
}

/** Does this tile's side run against the strip's record (node 0 at the record's node n)? */
function reversedOn(ctx, strip, addr, side) {
  const n = ctx.layout.tileCells, L = ctx.layout;
  const d = (p, q) => (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2;
  const p0 = nodePosition(addr, n, ...sideNode(side, 0, n), undefined, L);
  return d(p0, nodePosition(strip.canonical, n, ...sideNode(strip.side, n, n), undefined, L)) < d(p0, nodePosition(strip.canonical, n, ...sideNode(strip.side, 0, n), undefined, L));
}

function finishMeshTile(ctx, addr, sol, strips, corners) {
  const L = ctx.layout, n = L.tileCells, m = n + 1;
  const elevation = new Float32Array(m * m), area = new Float32Array(m * m), temperature = new Float32Array(m * m);
  const uplift = new Float32Array(m * m), precip = new Float32Array(m * m), material = new Uint8Array(m * m);
  const sed = new Float32Array(CLASSES * m * m);
  const v = [0, 0, 0];
  // border nodes come from the shared records (strip and corner values sampled once,
  // at the canonical tile's nodes), so neighbours agree bitwise; the rest from the tile
  const rev = Object.fromEntries(SIDES.map((sd) => [sd, reversedOn(ctx, strips[sd], addr, sd)]));
  const borderAt = new Map();
  for (const sd of SIDES) for (let t = 0; t <= n; t++) {
    const [ci, cj] = sideNode(sd, t, n);
    borderAt.set(cj * m + ci, strips[sd].record[rev[sd] ? n - t : t]);
  }
  let hint = -1;
  for (let cj = 0; cj <= n; cj++) for (let ci = 0; ci <= n; ci++) {
    const c = cj * m + ci;
    let pv = borderAt.get(c);
    if (!pv) {
      nodePosition(addr, n, ci, cj, v, L);
      pv = pointValues(sol, v, hint);
      hint = pv.tri;
    }
    elevation[c] = pv.z; area[c] = pv.A;
    for (let k = 0; k < CLASSES; k++) sed[k * m * m + c] = pv.sed[k];
    material[c] = pv.material; uplift[c] = pv.U; precip[c] = 1000 * pv.P;
    temperature[c] = pv.tsl - LAPSE * Math.max(0, elevation[c]);
  }
  // drainage on the raster for the levels below (they read links and areas from it)
  const rs = new RasterSolver(m, m, sol.g.dist[0] ? meanSpacing(sol) : 1000);
  const outlet = new Uint8Array(m * m);
  for (let c = 0; c < m * m; c++) { const ci = c % m, cj = (c - ci) / m; if (ci === 0 || cj === 0 || ci === n || cj === n || elevation[c] <= 0) outlet[c] = 1; }
  rs.route(Float64Array.from(elevation), outlet);
  const receiver = Int32Array.from(rs.receiver, (q, c) => (q === c ? -1 : q));
  // statistics and river graph on the mesh itself, over the tile's free nodes
  const N = sol.n, inTile = new Uint8Array(N), r = tileRect(addr, L);
  for (let i = 0; i < N; i++) {
    const fa = sol.fab[2 * i], fb = sol.fab[2 * i + 1];
    if (fa >= r.a0 && fa <= r.a0 + r.size && fb >= r.b0 && fb <= r.b0 + r.size) inTile[i] = 1;
  }
  // consistency over the tile itself: block means of its nodes against the parent's
  const blocks = new Map();
  for (let i = 0; i < N; i++) {
    if (!inTile[i] || (sol.flags[i] & FLAG_OUTLET && !(sol.flags[i] & FLAG_FIXED))) continue;
    const e = blocks.get(sol.block[i]) || { n: 0, d: 0 };
    e.n++; e.d += sol.z[i] - sol.base[i] - (sol.bias.get(sol.block[i]) || 0);
    blocks.set(sol.block[i], e);
  }
  let consistency = 0;
  for (const e of blocks.values()) if (e.n >= 400) consistency = Math.max(consistency, Math.abs(e.d / e.n));
  const recv = Int32Array.from(sol.rec);
  const order = upstreamOrder(recv);
  const inherited = Float64Array.from(sol.inflow);
  for (const i of order) if (recv[i] >= 0) inherited[recv[i]] += inherited[i];
  const mask = Uint8Array.from(inTile, (t, i) => (t && !sol.pinned[i] && inherited[i] <= 0.5 * sol.A[i] ? 1 : 0));
  const d = { receiver: recv, dist: Float32Array.from(recv, (q, i) => (q < 0 ? 1 : sol.recDist[i])), area: sol.A, z: sol.z, mask, coords: (i) => [sol.fab[2 * i], sol.fab[2 * i + 1]] };
  const stats = drainageStatistics(d, { channelArea: ctx.o.channelArea, flatSlope: 3 * ctx.o.eps });
  const rivers = extractRiverGraph(d, { channelArea: ctx.o.channelArea });
  return {
    addr, n, nodes: m, dx: meanSpacing(sol), mesh: true,
    fields: { elevation, uplift, precip, material, drainageArea: area, temperature, sediment: sed },
    drainage: { receiver, area, nx: m, dx: meanSpacing(sol) },
    stats: { ...stats, consistency },
    rivers,
    solution: { n: sol.n, ids: sol.ids, z: sol.z, rec: sol.rec, A: sol.A, pinned: sol.pinned, inTile },
    strips: Object.fromEntries(SIDES.map((s) => [s, strips[s].key])),
  };
}

function meanSpacing(sol) {
  let s = 0;
  for (let e = 0; e < sol.g.dist.length; e++) s += sol.g.dist[e];
  return s / sol.g.dist.length;
}

/**
 * The border record of a tile side, for the level below (same shape as seams.js
 * sideBorder): elevations at the side's n + 1 nodes and the rivers crossing it, both
 * from the shared strip (and corners), so the two tiles of a seam give the same record.
 */
export async function meshSideBorder(ctx, addr, side, makeSolver = cpuGraphSolver) {
  const L = ctx.layout, n = L.tileCells;
  const strip = await edgeStrip(ctx, addr, side, makeSolver);
  const reversed = reversedOn(ctx, strip, addr, side);
  const z = new Float32Array(n + 1);
  for (let t = 0; t <= n; t++) z[t] = strip.record[reversed ? n - t : t].z;
  // crossings: strip drainage links that cross the edge line, with enough area
  const { line, lo, hi, vertical } = strip.geometry;
  const across = (i) => (vertical ? strip.fab[2 * i] : strip.fab[2 * i + 1]) - line;
  const along = (i) => (vertical ? strip.fab[2 * i + 1] : strip.fab[2 * i]);
  const raw = [];
  for (let i = 0; i < strip.n; i++) {
    const q = strip.rec[i];
    if (q < 0 || strip.A[i] < ctx.o.channelArea) continue;
    const a0 = across(i), a1 = across(q);
    if ((a0 >= 0) === (a1 >= 0)) continue;
    const lam = a0 / (a0 - a1), pos = along(i) + lam * (along(q) - along(i));
    if (pos < lo || pos > hi) continue;
    raw.push({ t: ((pos - lo) / (hi - lo)) * n, area: strip.A[i], sed: Array.from({ length: CLASSES }, (_, k) => strip.Qs[k * strip.n + i]), towardPlus: a1 > a0 });
  }
  // keep the larger of crossings closer than 3 nodes
  raw.sort((a, b) => b.area - a.area);
  const kept = [];
  for (const x of raw) if (!kept.some((k) => Math.abs(k.t - x.t) < 3)) kept.push(x);
  // orientation: is this tile the strip's canonical tile, and is its node order the same?
  const mine = tileKey(strip.canonical) === tileKey(addr);
  // for side W (or S) of the canonical tile its interior has across > 0; for E (or N), across < 0
  const interiorPlus = strip.side === 'W' || strip.side === 'S';
  const crossings = kept.map((x) => {
    const intoC = x.towardPlus === interiorPlus;
    const node = Math.max(1, Math.min(n - 1, Math.round(reversed ? n - x.t : x.t)));
    return { node, area: x.area, sed: x.sed, inflow: mine ? intoC : !intoC, z: z[node], links: [] };
  });
  return { key: strip.key, z, crossings };
}
