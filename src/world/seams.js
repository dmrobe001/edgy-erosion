// Seams: corner and edge records shared by neighbouring tiles.
//
// Neighbouring tiles must agree exactly where they meet, whatever order they are made
// in, without either reading the other. So everything on a tile's border comes from
// records that depend only on the parent and the world seed (design: dimension-ordered
// generation):
//
//   corner (0D)  shared by three or four tiles. Its position is computed in the frame of
//                the lowest-keyed tile that shares it, and its elevation is the parent's,
//                plus the level's subgrid noise at that position.
//   edge (1D)    shared by two tiles: n + 1 nodes, computed in the frame of the
//                lower-keyed tile; the end nodes are the corners. Elevations as for
//                corners. Plus the river crossings: every link of the parent's drainage
//                tree that crosses the edge, placed (with seeded jitter) within the
//                stretch the coarse link spans, with its drainage area and direction.
//                A crossing notches the edge profile down to the river.
//
// Each tile reads the records for its four sides, oriented into its own node order by
// matching end points (so it works across cube edges, where neighbouring faces' axes are
// rotated or reversed). The border values are then bitwise identical from both sides.

import {
  WORLD_LAYOUT, tileKey, neighborOf, cornerTiles, nodePosition, tileRect, gnomonic, toGnomonic, vectorToFace,
} from './cubesphere.js';
import { worldKey, deriveKey, uniform } from './hash.js';

export const SIDES = ['W', 'E', 'S', 'N'];
const CORNER_OF = { W: ['SW', 'NW'], E: ['SE', 'NE'], S: ['SW', 'SE'], N: ['NW', 'NE'] };
const CORNER_NODE = { SW: [0, 0], SE: [1, 0], NW: [0, 1], NE: [1, 1] };

/** Node (ci, cj) of side `side` at position t = 0..n, in the tile's own order. */
export function sideNode(side, t, n) {
  if (side === 'W') return [0, t];
  if (side === 'E') return [n, t];
  if (side === 'S') return [t, 0];
  return [t, n];
}

/**
 * A seam context: what the records of one level are made from.
 * @param {object} o
 * @param {number} o.level the level of the tiles whose seams these are
 * @param {string} o.worldSeed
 * @param {(v:number[]) => number} o.elevation initial elevation at a unit vector (parent + noise)
 * @param {(v:number[]) => number} o.base the parent's smooth elevation at a unit vector
 * @param {(v:number[]) => number} o.notchDepth how far a crossing cuts below the parent's surface
 * @param {{from:number[], to:number[], area:number, id:number}[]} o.links the parent's drainage links (unit vectors)
 * @param {number} o.linkSpacing typical length of a parent link (m)
 * @param {number} o.radius planet radius (m)
 * @param {(C:object, sC:string, pos:number[][]) => ({z:Float32Array, crossings:object[]}|null)} [o.edgeOverride]
 *   for edges that must come from somewhere other than the parent's interior (an edge
 *   lying on the parent's own seam takes that seam's record); null to use the default
 * @param {(tile:object, corner:string, pos:number[]) => (number|null)} [o.cornerOverride] likewise for corners
 */
export function seamContext(o) {
  return { layout: WORLD_LAYOUT, cache: new Map(), ...o, n: (o.layout || WORLD_LAYOUT).tileCells };
}

function cornerRecord(ctx, addr, corner) {
  const tiles = cornerTiles(addr, corner, ctx.layout);
  const key = `C:${tiles.map(tileKey).join('|')}`;
  if (ctx.cache.has(key)) return ctx.cache.get(key);
  // the canonical tile's corner that coincides with this one
  const c0 = tiles[0];
  const n = ctx.n;
  const me = nodePosition(addr, n, CORNER_NODE[corner][0] * n, CORNER_NODE[corner][1] * n, undefined, ctx.layout);
  let pos = null, best = Infinity;
  for (const nm of Object.keys(CORNER_NODE)) {
    const p = nodePosition(c0, n, CORNER_NODE[nm][0] * n, CORNER_NODE[nm][1] * n, undefined, ctx.layout);
    const d = (p[0] - me[0]) ** 2 + (p[1] - me[1]) ** 2 + (p[2] - me[2]) ** 2;
    if (d < best) { best = d; pos = p; }
  }
  const over = ctx.cornerOverride ? ctx.cornerOverride(c0, nameAt(c0, pos, ctx), pos) : null;
  const rec = { key, pos, z: over ?? ctx.elevation(pos) };
  ctx.cache.set(key, rec);
  return rec;
}

function nameAt(tile, pos, ctx) {
  const n = ctx.n;
  let best = null, bd = Infinity;
  for (const nm of Object.keys(CORNER_NODE)) {
    const p = nodePosition(tile, n, CORNER_NODE[nm][0] * n, CORNER_NODE[nm][1] * n, undefined, ctx.layout);
    const d = dist2(p, pos);
    if (d < bd) { bd = d; best = nm; }
  }
  return best;
}

function sideToward(tile, other, layout) {
  const k = tileKey(other);
  for (const s of SIDES) if (tileKey(neighborOf(tile, s, layout)) === k) return s;
  throw new Error(`${tileKey(tile)} is not next to ${k}`);
}

/** The canonical record of the edge on `side` of `addr`. */
export function edgeRecord(ctx, addr, side) {
  const nb = neighborOf(addr, side, ctx.layout);
  const [ka, kb] = [tileKey(addr), tileKey(nb)];
  const key = `E:${ka < kb ? ka : kb}|${ka < kb ? kb : ka}`;
  if (ctx.cache.has(key)) return ctx.cache.get(key);
  const C = ka < kb ? addr : nb;
  const D = ka < kb ? nb : addr;
  const sC = ka < kb ? side : sideToward(nb, addr, ctx.layout);
  const n = ctx.n;
  const pos = [], z = new Float32Array(n + 1);
  const [cs, ce] = CORNER_OF[sC].map((c) => cornerRecord(ctx, C, c));
  for (let t = 0; t <= n; t++) {
    if (t === 0) { pos.push(cs.pos); z[t] = cs.z; continue; }
    if (t === n) { pos.push(ce.pos); z[t] = ce.z; continue; }
    const [ci, cj] = sideNode(sC, t, n);
    const p = nodePosition(C, n, ci, cj, undefined, ctx.layout);
    pos.push(p);
    z[t] = ctx.elevation(p);
  }
  const over = ctx.edgeOverride ? ctx.edgeOverride(C, sC, pos) : null;
  if (over) {
    // keep the shared corners exact
    over.z[0] = cs.z; over.z[n] = ce.z;
    const rec = { key, canonical: tileKey(C), other: tileKey(D), side: sC, pos, z: over.z, crossings: over.crossings };
    ctx.cache.set(key, rec);
    return rec;
  }
  const crossings = findCrossings(ctx, C, sC, key);
  if (ctx.minorSpacing > 0) crossings.push(...minorCrossings(ctx, C, sC, z, crossings));
  crossings.sort((x, y) => x.node - y.node);
  for (const x of crossings) {
    if (x.minor) continue; // a minor crossing sits at a low point of the profile already
    // the notch: the river's level, below the parent's smooth surface
    const notch = ctx.base(pos[x.node]) - ctx.notchDepth(pos[x.node]);
    x.z = Math.fround(Math.min(notch, ...Array.from(z.subarray(Math.max(1, x.node - 2), Math.min(n, x.node + 3)))));
    z[x.node] = x.z;
  }
  const rec = { key, canonical: tileKey(C), other: tileKey(D), side: sC, pos, z, crossings };
  ctx.cache.set(key, rec);
  return rec;
}

/**
 * Parent drainage links crossing side sC of tile C. Great circles and the tile's edges
 * are both straight lines in the gnomonic projection of C's face, so the crossing is a
 * line-segment intersection there.
 */
function findCrossings(ctx, C, sC, edgeKeyStr) {
  const n = ctx.n;
  const r = tileRect(C, ctx.layout);
  const horizontal = sC === 'S' || sC === 'N';
  // the edge as a line u = const (W/E) or v = const (S/N) in gnomonic coordinates
  const fixedA = sC === 'W' ? r.a0 : sC === 'E' ? r.a0 + r.size : sC === 'S' ? r.b0 : r.b0 + r.size;
  const g0 = gnomonic(fixedA);
  const lo = gnomonic(horizontal ? r.a0 : r.b0), hi = gnomonic((horizontal ? r.a0 : r.b0) + r.size);
  const outward = sC === 'E' || sC === 'N' ? 1 : -1;
  const edgeLen = (r.size * Math.PI / 2) * ctx.radius;
  const key = deriveKey(worldKey(ctx.worldSeed), ctx.level, 0x65646765, ...Array.from(worldKey(edgeKeyStr)));
  const found = [];
  for (const L of ctx.links) {
    const p = toGnomonic(C.face, ...L.from), q = toGnomonic(C.face, ...L.to);
    if (!p || !q) continue;
    const pa = horizontal ? p[1] : p[0], qa = horizontal ? q[1] : q[0]; // across the edge
    // an end point exactly on the line counts as being on its far (+) side, so a river
    // passing through a node that lies on the edge is seen crossing exactly once
    if ((pa >= g0) === (qa >= g0)) continue;
    const lam = qa === pa ? 0 : (g0 - pa) / (qa - pa);
    const along = horizontal ? p[0] + lam * (q[0] - p[0]) : p[1] + lam * (q[1] - p[1]);
    if (along <= lo || along >= hi) continue;
    // node along the edge, in C's increasing direction, jittered within the link's span
    const t = ((Math.atan(along) / (Math.PI / 4) + 1) / 2 - (horizontal ? r.a0 : r.b0)) / r.size;
    const span = Math.max(0, Math.min(0.25, (0.5 * ctx.linkSpacing) / edgeLen));
    const jitter = (uniform(key, L.id) - 0.5) * span;
    const node = Math.max(3, Math.min(n - 3, Math.round((t + jitter) * n)));
    const intoC = outward > 0 ? qa < g0 : qa >= g0; // flows from beyond the edge into C
    found.push({ node, area: L.area, sed: L.sed ? Array.from(L.sed) : null, inflowToCanonical: intoC, links: [L.id] });
  }
  // merge crossings closer than 4 nodes (the larger river's place wins); areas of
  // rivers flowing the same way add up
  found.sort((x, y) => y.area - x.area || x.node - y.node);
  const kept = [];
  for (const x of found) {
    const near = kept.find((k) => Math.abs(k.node - x.node) < 4);
    if (!near) { kept.push(x); continue; }
    if (near.inflowToCanonical === x.inflowToCanonical) {
      near.area += x.area; near.links.push(...x.links);
      if (near.sed && x.sed) for (let k = 0; k < near.sed.length; k++) near.sed[k] += x.sed[k];
    }
  }
  return kept.sort((x, y) => x.node - y.node);
}

/**
 * Minor crossings: streams too small for the parent's drainage tree, which would
 * otherwise turn every tile edge into a straight divide. At low points of the edge
 * profile (lowest within minorSpacing / 2 nodes either way, and minorSpacing from any
 * other crossing) a stream crosses, flowing from the side where the parent's surface
 * stands higher; it drains a strip about minorSpacing wide and as long as it is wide.
 */
function minorCrossings(ctx, C, sC, z, existing) {
  const n = ctx.n, w = Math.max(2, Math.floor(ctx.minorSpacing / 2));
  const out = [];
  const taken = existing.map((x) => x.node);
  const inward = { W: [1, 0], E: [-1, 0], S: [0, 1], N: [0, -1] }[sC];
  const step = Math.max(2, Math.round(ctx.minorSpacing / 4));
  const r = tileRect(C, ctx.layout);
  const dx = (r.size * Math.PI / 2) * ctx.radius / n;
  for (let t = w; t <= n - w; t++) {
    let low = true;
    for (let k = t - w; k <= t + w && low; k++) if (k !== t && z[k] < z[t]) low = false;
    if (!low || taken.some((u) => Math.abs(u - t) < ctx.minorSpacing)) continue;
    const [ci, cj] = sideNode(sC, t, n);
    const inside = ctx.base(nodePosition(C, n, ci + step * inward[0], cj + step * inward[1], undefined, ctx.layout));
    const outside = ctx.base(nodePosition(C, n, ci - step * inward[0], cj - step * inward[1], undefined, ctx.layout));
    if (inside === outside) continue;
    // a small stream crosses near the parent's surface, not at the bottom of a noise
    // trough (which would drop base level along the whole edge)
    const p = nodePosition(C, n, ci, cj, undefined, ctx.layout);
    z[t] = Math.fround(Math.max(z[t], ctx.base(p) - 0.5 * ctx.notchDepth(p)));
    out.push({ node: t, area: (ctx.minorSpacing * dx) ** 2, inflowToCanonical: outside > inside, links: [], minor: true, z: z[t] });
    taken.push(t);
  }
  return out;
}

/**
 * The border of a tile in its own node order: for each side, the n + 1 elevations and
 * the crossings ({node, area, inflow: true if the river flows into this tile, z, links}).
 */
export function tileBorder(ctx, addr) {
  const out = {};
  for (const side of SIDES) out[side] = sideBorder(ctx, addr, side);
  return out;
}

/** One side of a tile's border, in the tile's own node order (see tileBorder). */
export function sideBorder(ctx, addr, side) {
  const n = ctx.n;
  {
    const rec = edgeRecord(ctx, addr, side);
    // orientation: does this tile's node 0 of the side sit at the record's node 0?
    const [ci, cj] = sideNode(side, 0, n);
    const p0 = nodePosition(addr, n, ci, cj, undefined, ctx.layout);
    const d0 = dist2(p0, rec.pos[0]), dn = dist2(p0, rec.pos[n]);
    const reversed = dn < d0;
    const z = new Float32Array(n + 1);
    for (let t = 0; t <= n; t++) z[t] = rec.z[reversed ? n - t : t];
    const mine = tileKey(addr) === rec.canonical;
    const crossings = rec.crossings.map((x) => ({
      node: reversed ? n - x.node : x.node, area: x.area, sed: x.sed || null, z: x.z, links: x.links, minor: !!x.minor,
      inflow: mine ? x.inflowToCanonical : !x.inflowToCanonical,
    }));
    return { key: rec.key, z, crossings };
  }
}

function dist2(a, b) { return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2; }

/** Which tile of `level` holds a unit vector, and where in it (fractional node coordinates). */
export function locateInTile(addr, v, n, layout = WORLD_LAYOUT) {
  const f = vectorToFace(v[0], v[1], v[2]);
  if (f.face !== addr.face) return null;
  const r = tileRect(addr, layout);
  const ci = ((f.a - r.a0) / r.size) * n, cj = ((f.b - r.b0) / r.size) * n;
  if (ci < 0 || cj < 0 || ci > n || cj > n) return null;
  return [ci, cj];
}
