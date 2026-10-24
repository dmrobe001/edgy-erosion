// The irregular mesh at level 1: windows of one global point set, and the graph solver on it.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { windowPoints, triangulateWindow, tileWindowRect, extendedFace, locateInWindow } from '../src/world/pointmesh.js';
import { GraphSolver, FLAG_OUTLET, CLASSES } from '../src/world/graphsolver.js';
import { tileRect, neighborOf } from '../src/world/cubesphere.js';

const R = 6.371e6, layout = { rootTiles: 32, branching: 32, tileCells: 32 };

test('point mesh: overlapping windows agree on the triangulation away from their edges, across a face edge too', () => {
  for (const [A, side] of [[{ level: 1, face: 1, i: 10, j: 19 }, 'E'], [{ level: 1, face: 0, i: 31, j: 5 }, 'E']]) {
    const B = neighborOf(A, side, layout);
    const wa = windowPoints('w', 1, A, tileWindowRect(A, 16, layout), layout);
    const wb = windowPoints('w', 1, B, tileWindowRect(B, 16, layout), layout);
    const ga = triangulateWindow(wa.pos, R), gb = triangulateWindow(wb.pos, R);
    const ib = new Map(Array.from(wb.ids, (id, i) => [id, i]));
    // nodes of A's window deep inside both windows: compare neighbour id sets
    const ra = tileRect(A, layout), h = ra.size / layout.tileCells;
    const rb = tileRect(B, layout);
    let compared = 0;
    for (let i = 0; i < ga.n; i++) {
      const j = ib.get(wa.ids[i]);
      if (j === undefined) continue;
      const ea = extendedFace(A.face, wa.pos[i]), eb = extendedFace(B.face, wb.pos[j]);
      const depth = Math.min(ea[0] - (ra.a0 - 16 * h), ra.a0 + ra.size + 16 * h - ea[0], ea[1] - (ra.b0 - 16 * h), ra.b0 + ra.size + 16 * h - ea[1],
        eb[0] - (rb.a0 - 16 * h), rb.a0 + rb.size + 16 * h - eb[0], eb[1] - (rb.b0 - 16 * h), rb.b0 + rb.size + 16 * h - eb[1]) / h;
      if (depth < 7) continue;
      const na = [], nb = [];
      for (let e = ga.off[i]; e < ga.off[i + 1]; e++) na.push(wa.ids[ga.nbr[e]]);
      for (let e = gb.off[j]; e < gb.off[j + 1]; e++) nb.push(wb.ids[gb.nbr[e]]);
      assert.deepEqual(na.sort((p, q) => p - q), nb.sort((p, q) => p - q), `node ${wa.ids[i]}`);
      compared++;
    }
    assert.ok(compared > 200, `compared ${compared}`);
  }
});

test('point mesh: point location gives barycentric weights that reproduce the point', () => {
  const A = { level: 1, face: 2, i: 3, j: 7 };
  const w = windowPoints('w', 1, A, tileWindowRect(A, 4, layout), layout);
  const g = triangulateWindow(w.pos, R);
  for (let k = 0; k < 50; k++) {
    const p = w.pos[(k * 37) % g.n], q = w.pos[(k * 91 + 5) % g.n];
    const v = [0.3 * p[0] + 0.7 * q[0], 0.3 * p[1] + 0.7 * q[1], 0.3 * p[2] + 0.7 * q[2]];
    const L = Math.hypot(...v);
    for (let t = 0; t < 3; t++) v[t] /= L;
    const loc = locateInWindow(g, w.pos, v);
    if (!loc) continue;
    const s = [0, 1, 2].map((t) => loc.w[0] * w.pos[loc.ids[0]][t] + loc.w[1] * w.pos[loc.ids[1]][t] + loc.w[2] * w.pos[loc.ids[2]][t]);
    const Ls = Math.hypot(...s);
    for (let t = 0; t < 3; t++) assert.ok(Math.abs(s[t] / Ls - v[t]) < 1e-9);
    for (const x of loc.w) assert.ok(x > -1e-9);
  }
});

function slab(G) {
  const A = { level: 1, face: 1, i: 10, j: 19 };
  const W = windowPoints('w', 1, A, tileWindowRect(A, 0, layout), layout);
  const g = triangulateWindow(W.pos, R);
  const n = g.n, r = tileRect(A, layout);
  const z = new Float32Array(n), U = new Float32Array(n).fill(1e-3), K = new Float32Array(n).fill(1e-5), P = new Float32Array(n).fill(1), flags = new Uint32Array(n);
  for (let i = 0; i < n; i++) {
    const y = (extendedFace(A.face, W.pos[i])[1] - r.b0) / r.size;
    z[i] = 1 + ((i * 7919) % 13) / 13 + 50 * y;
    if (y < 0.03) { flags[i] = FLAG_OUTLET; z[i] = 0; }
  }
  const s = new GraphSolver(g, { z, U, K, P, flags });
  s.run(2500, { dt: 5e3, G: [G, G, G] });
  let err = 0, cnt = 0, out = 0, upl = 0;
  for (let i = 0; i < n; i++) {
    if (flags[i] & 1) { for (let k = 0; k < CLASSES; k++) out += s.Qs[k * n + i]; continue; }
    upl += 1e-3 * g.area[i];
    const q = s.rec[i];
    if (q < 0 || s.A[i] < 30 * g.area[i]) continue;
    err += (s.z[i] - s.z[q]) / s.recDist[i] / (100 * s.A[i] ** -0.5) - 1; cnt++;
  }
  return { err: err / cnt, balance: out / upl };
}

test('graph solver: without deposition the steady state is K A^m S = U; with it, mass still balances', () => {
  const bare = slab(0);
  assert.ok(Math.abs(bare.err) < 0.01, `slope error ${bare.err}`);
  assert.ok(Math.abs(bare.balance - 1) < 0.01, `sediment out / uplift ${bare.balance}`);
  const dep = slab(1);
  // deposition steepens the channels (Davy & Lague: S = (U + G Qs / A) / (K A^m), and Qs = U A at steady state, so S doubles at G = 1)
  assert.ok(dep.err > 0.5 && dep.err < 1.2, `slope excess ${dep.err}`);
  assert.ok(Math.abs(dep.balance - 1) < 0.01, `sediment out / uplift ${dep.balance}`);
});
