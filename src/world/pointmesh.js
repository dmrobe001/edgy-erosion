// A global irregular mesh, seen through windows.
//
// Every level has one point per cell of a fine grid on each cube face (the tiles' cells),
// jittered within its cell by a stateless hash of (world, level, face, cell). The point
// set is therefore defined everywhere at once, and any window onto it (a tile, a seam
// strip, a corner patch) contains exactly the same points where windows overlap, with
// the same global ids.
//
// A window is triangulated (Delaunay, Bowyer-Watson) in a stereographic projection
// centred on it. Stereographic projection maps circles on the sphere to circles in the
// plane, so the planar Delaunay triangulation is the spherical one: two windows that
// overlap agree on every triangle whose circumcircle lies inside both. Lengths, areas
// and diffusion weights are then computed in 3-D, so they do not depend on which window
// (which projection) they were computed in either.
//
// Routing on this mesh has no preferred directions: a node has about six neighbours in
// random directions, instead of the eight fixed directions of a square grid.

import { faceToVector, tileRect, toGnomonic, tilesPerSide, neighborOf, cornerTiles, tileKey, WORLD_LAYOUT } from './cubesphere.js';
import { worldKey, deriveKey, mix32 } from './hash.js';

export const JITTER = 0.8; // share of a cell the point may wander over (centred)

/** Cells per face side at a level. */
export function cellsPerFace(level, layout = WORLD_LAYOUT) {
  return tilesPerSide(level, layout) * layout.tileCells;
}

/** Global id of a cell's point. */
export function pointId(face, gi, gj, N) {
  return (face * N + gj) * N + gi;
}

/** The point of cell (gi, gj) of a face, as face coordinates [a, b]. */
export function cellPoint(key, face, gi, gj, N) {
  const u = mix32(key, face * 1048576 + gi, gj, 1) / 4294967296;
  const v = mix32(key, face * 1048576 + gi, gj, 2) / 4294967296;
  return [(gi + 0.5 + JITTER * (u - 0.5)) / N, (gj + 0.5 + JITTER * (v - 0.5)) / N];
}

export function pointKey(worldSeed, level) {
  return deriveKey(worldKey(worldSeed), level, 0x706f696e);
}

/**
 * Face coordinates of a unit vector in the frame of `face`, extended beyond the face
 * (the gnomonic plane continues past its edges), or null if it faces away.
 */
export function extendedFace(face, v) {
  const g = toGnomonic(face, v[0], v[1], v[2]);
  if (!g) return null;
  return [(Math.atan(g[0]) / (Math.PI / 4) + 1) / 2, (Math.atan(g[1]) / (Math.PI / 4) + 1) / 2];
}

/**
 * The points of a window: every point whose extended face coordinates in the frame of
 * tile `frame` fall inside the rectangle [a0, a1] x [b0, b1]. Candidates come from the
 * frame tile and its neighbours (the rectangle may not reach further than one tile out).
 * Returned in a canonical order (by global id), so every window lists shared points alike.
 */
export function windowPoints(worldSeed, level, frame, rect, layout = WORLD_LAYOUT) {
  const key = pointKey(worldSeed, level);
  const N = cellsPerFace(level, layout), n = layout.tileCells;
  const tiles = new Map([[tileKey(frame), frame]]);
  for (const s of ['W', 'E', 'S', 'N']) { const t = neighborOf(frame, s, layout); tiles.set(tileKey(t), t); }
  for (const c of ['SW', 'SE', 'NW', 'NE']) for (const t of cornerTiles(frame, c, layout)) tiles.set(tileKey(t), t);
  const ids = [], pos = [];
  const seen = new Set();
  for (const t of tiles.values()) {
    for (let cj = 0; cj < n; cj++) for (let ci = 0; ci < n; ci++) {
      const gi = t.i * n + ci, gj = t.j * n + cj;
      const id = pointId(t.face, gi, gj, N);
      if (seen.has(id)) continue;
      const [a, b] = cellPoint(key, t.face, gi, gj, N);
      const v = faceToVector(t.face, a, b);
      let fa = a, fb = b;
      if (t.face !== frame.face) {
        const e = extendedFace(frame.face, v);
        if (!e) continue;
        [fa, fb] = e;
      }
      if (fa < rect.a0 || fa > rect.a1 || fb < rect.b0 || fb > rect.b1) continue;
      seen.add(id);
      ids.push(id);
      pos.push(v);
    }
  }
  const order = ids.map((_, k) => k).sort((p, q) => ids[p] - ids[q]);
  return { ids: Float64Array.from(order, (k) => ids[k]), pos: order.map((k) => pos[k]) };
}

/**
 * Planar Delaunay triangulation (Bowyer-Watson with walking point location).
 * @param {Float64Array} xs
 * @param {Float64Array} ys
 * @returns {Int32Array} triangles, counter-clockwise vertex triples
 */
export function delaunay(xs, ys) {
  const n = xs.length;
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (let i = 0; i < n; i++) { minX = Math.min(minX, xs[i]); maxX = Math.max(maxX, xs[i]); minY = Math.min(minY, ys[i]); maxY = Math.max(maxY, ys[i]); }
  const span = Math.max(maxX - minX, maxY - minY) || 1, cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
  // vertices n, n+1, n+2: a super triangle far outside
  const X = new Float64Array(n + 3), Y = new Float64Array(n + 3);
  X.set(xs); Y.set(ys);
  X[n] = cx - 40 * span; Y[n] = cy - 30 * span;
  X[n + 1] = cx + 40 * span; Y[n + 1] = cy - 30 * span;
  X[n + 2] = cx; Y[n + 2] = cy + 40 * span;
  let cap = 8 * n + 16;
  let V = new Int32Array(3 * cap), NB = new Int32Array(3 * cap), alive = new Uint8Array(cap);
  let T = 0;
  const grow = () => {
    cap *= 2;
    const V2 = new Int32Array(3 * cap), N2 = new Int32Array(3 * cap), A2 = new Uint8Array(cap);
    V2.set(V); N2.set(NB); A2.set(alive);
    V = V2; NB = N2; alive = A2;
  };
  const add = (a, b, c) => { if (T >= cap) grow(); V[3 * T] = a; V[3 * T + 1] = b; V[3 * T + 2] = c; NB[3 * T] = NB[3 * T + 1] = NB[3 * T + 2] = -1; alive[T] = 1; return T++; };
  const orient = (a, b, px, py) => (X[b] - X[a]) * (py - Y[a]) - (Y[b] - Y[a]) * (px - X[a]);
  const inCircle = (t, px, py) => {
    const a = V[3 * t], b = V[3 * t + 1], c = V[3 * t + 2];
    const adx = X[a] - px, ady = Y[a] - py, bdx = X[b] - px, bdy = Y[b] - py, cdx = X[c] - px, cdy = Y[c] - py;
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) - (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady) > 0;
  };
  add(n, n + 1, n + 2);
  let last = 0;
  let badMark = new Int32Array(cap).fill(-1);
  const stack = [], cavity = [];
  for (let p = 0; p < n; p++) {
    const px = X[p], py = Y[p];
    // walk to the triangle containing p
    let t = last, guard = 0;
    walk: for (;;) {
      if (++guard > 4 * n + 100) throw new Error('delaunay: walk did not converge');
      for (let k = 0; k < 3; k++) {
        const kk = (k + guard) % 3;
        const a = V[3 * t + kk], b = V[3 * t + ((kk + 1) % 3)];
        if (orient(a, b, px, py) < 0) { t = NB[3 * t + kk]; continue walk; }
      }
      break;
    }
    // the cavity: triangles whose circumcircle holds p
    if (badMark.length < cap) { const m2 = new Int32Array(cap).fill(-1); m2.set(badMark); badMark = m2; }
    cavity.length = 0;
    stack.length = 0;
    stack.push(t);
    badMark[t] = p;
    while (stack.length) {
      const u = stack.pop();
      cavity.push(u);
      for (let k = 0; k < 3; k++) {
        const w = NB[3 * u + k];
        if (w >= 0 && badMark[w] !== p && inCircle(w, px, py)) { badMark[w] = p; stack.push(w); }
      }
    }
    // its boundary, re-triangulated as a fan around p
    const byStart = new Map(), made = [];
    for (const u of cavity) {
      for (let k = 0; k < 3; k++) {
        const w = NB[3 * u + k];
        if (w >= 0 && badMark[w] === p) continue;
        const a = V[3 * u + k], b = V[3 * u + ((k + 1) % 3)];
        const nt = add(a, b, p);
        NB[3 * nt] = w;
        if (w >= 0) for (let q = 0; q < 3; q++) if (NB[3 * w + q] === u) NB[3 * w + q] = nt;
        byStart.set(a, nt);
        made.push(nt);
      }
      alive[u] = 0;
    }
    for (const nt of made) {
      const b = V[3 * nt + 1];
      const next = byStart.get(b); // the fan triangle starting where this one ends
      NB[3 * nt + 1] = next;
      NB[3 * next + 2] = nt;
    }
    last = made[0];
  }
  const out = [];
  for (let t = 0; t < T; t++) {
    if (!alive[t]) continue;
    const a = V[3 * t], b = V[3 * t + 1], c = V[3 * t + 2];
    if (a >= n || b >= n || c >= n) continue;
    out.push(a, b, c);
  }
  return Int32Array.from(out);
}

/**
 * Triangulate a window of points on the sphere and compute the graph the solvers use.
 * @param {number[][]} pos unit vectors
 * @param {number} radius planet radius (m)
 * @returns {{n, tri, off, nbr, dist, lapw, area, center}}
 *   off/nbr: neighbours in compressed rows; dist: edge length (m); lapw: diffusion
 *   weight of the edge (Voronoi face length over distance, from cotangents); area:
 *   area of each node's cell (m^2, a third of each triangle around it)
 */
export function triangulateWindow(pos, radius) {
  const n = pos.length;
  // window centre and a tangent frame there
  let cx = 0, cy = 0, cz = 0;
  for (const p of pos) { cx += p[0]; cy += p[1]; cz += p[2]; }
  const L = Math.hypot(cx, cy, cz); cx /= L; cy /= L; cz /= L;
  const ref = Math.abs(cz) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  let e1 = [cy * ref[2] - cz * ref[1], cz * ref[0] - cx * ref[2], cx * ref[1] - cy * ref[0]];
  const l1 = Math.hypot(...e1); e1 = e1.map((x) => x / l1);
  const e2 = [cy * e1[2] - cz * e1[1], cz * e1[0] - cx * e1[2], cx * e1[1] - cy * e1[0]];
  // stereographic projection from the antipode of the centre (circles stay circles)
  const xs = new Float64Array(n), ys = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const p = pos[i], d = 1 + p[0] * cx + p[1] * cy + p[2] * cz;
    xs[i] = (2 * (p[0] * e1[0] + p[1] * e1[1] + p[2] * e1[2])) / d;
    ys[i] = (2 * (p[0] * e2[0] + p[1] * e2[1] + p[2] * e2[2])) / d;
  }
  const tri = delaunay(xs, ys);
  // graph from the triangles, with 3-D geometry
  const nbSets = Array.from({ length: n }, () => new Map()); // neighbour -> lap weight
  const area = new Float64Array(n);
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const crs = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  for (let t = 0; t < tri.length; t += 3) {
    const ids = [tri[t], tri[t + 1], tri[t + 2]];
    const P = ids.map((i) => pos[i]);
    const A2 = Math.hypot(...crs(sub(P[1], P[0]), sub(P[2], P[0])));
    for (let k = 0; k < 3; k++) {
      const i = ids[k], j = ids[(k + 1) % 3], o = ids[(k + 2) % 3];
      area[i] += (A2 / 6) * radius * radius;
      // cotangent of the angle at o, opposite edge ij
      const u = sub(pos[i], pos[o]), w = sub(pos[j], pos[o]);
      const cot = dot(u, w) / (Math.hypot(...crs(u, w)) || 1e-30);
      nbSets[i].set(j, (nbSets[i].get(j) || 0) + 0.5 * cot);
      nbSets[j].set(i, (nbSets[j].get(i) || 0) + 0.5 * cot);
    }
  }
  const off = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) off[i + 1] = off[i] + nbSets[i].size;
  const nbr = new Int32Array(off[n]), dist = new Float32Array(off[n]), lapw = new Float32Array(off[n]);
  for (let i = 0; i < n; i++) {
    const list = [...nbSets[i].entries()].sort((p, q) => p[0] - q[0]);
    list.forEach(([j, w], k) => {
      nbr[off[i] + k] = j;
      dist[off[i] + k] = radius * Math.hypot(...sub(pos[i], pos[j]));
      lapw[off[i] + k] = Math.max(0, w); // obtuse triangles can give negative weights: clipped (stable)
    });
  }
  return { n, tri, off, nbr, dist, lapw, area: Float32Array.from(area), center: [cx, cy, cz], xs, ys, e1, e2 };
}

/**
 * Locate a unit vector in a triangulated window: the triangle and barycentric weights,
 * computed in 3-D (the plane through the triangle's corners) so that any window holding
 * the same triangle gives the same weights. Returns null outside the triangulation.
 * `hint` speeds up a walk; `g.triOf` (built on first use) maps nodes to a triangle.
 */
export function locateInWindow(g, pos, v, hint = -1) {
  const { tri, center: c, e1, e2 } = g;
  if (!g.triAdj) buildAdjacency(g);
  const d = 1 + v[0] * c[0] + v[1] * c[1] + v[2] * c[2];
  const px = (2 * (v[0] * e1[0] + v[1] * e1[1] + v[2] * e1[2])) / d, py = (2 * (v[0] * e2[0] + v[1] * e2[1] + v[2] * e2[2])) / d;
  const X = g.xs, Y = g.ys;
  let t = hint >= 0 ? hint : g.lastTri || 0;
  for (let guard = 0; guard < 100000; guard++) {
    let moved = false;
    for (let k = 0; k < 3; k++) {
      const kk = (k + guard) % 3;
      const a = tri[3 * t + kk], b = tri[3 * t + ((kk + 1) % 3)];
      if ((X[b] - X[a]) * (py - Y[a]) - (Y[b] - Y[a]) * (px - X[a]) < 0) {
        const nt = g.triAdj[3 * t + kk];
        if (nt < 0) return null;
        t = nt; moved = true; break;
      }
    }
    if (!moved) break;
  }
  g.lastTri = t;
  const ia = tri[3 * t], ib = tri[3 * t + 1], ic = tri[3 * t + 2];
  const A = pos[ia], B = pos[ib], C = pos[ic];
  // v = wa A + wb B + wc C (up to scale): Cramer's rule with triple products
  const det = (p, q, r) => p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0]) + p[2] * (q[0] * r[1] - q[1] * r[0]);
  const D = det(A, B, C);
  let wa = det(v, B, C) / D, wb = det(A, v, C) / D, wc = det(A, B, v) / D;
  const s = wa + wb + wc;
  wa /= s; wb /= s; wc /= s;
  return { tri: t, ids: [ia, ib, ic], w: [wa, wb, wc] };
}

function buildAdjacency(g) {
  const { tri } = g;
  const T = tri.length / 3;
  const adj = new Int32Array(3 * T).fill(-1);
  const edge = new Map();
  for (let t = 0; t < T; t++) for (let k = 0; k < 3; k++) {
    const a = tri[3 * t + k], b = tri[3 * t + ((k + 1) % 3)];
    const key = a < b ? a * 4294967296 + b : b * 4294967296 + a;
    const o = edge.get(key);
    if (o === undefined) edge.set(key, 3 * t + k);
    else { adj[3 * t + k] = Math.floor(o / 3); adj[o] = t; }
  }
  g.triAdj = adj;
}

/** The rectangle (extended face coordinates in the frame of `addr`) of a tile plus a margin of `cells`. */
export function tileWindowRect(addr, cells, layout = WORLD_LAYOUT) {
  const r = tileRect(addr, layout), h = r.size / layout.tileCells;
  return { a0: r.a0 - cells * h, a1: r.a0 + r.size + cells * h, b0: r.b0 - cells * h, b1: r.b0 + r.size + cells * h };
}
