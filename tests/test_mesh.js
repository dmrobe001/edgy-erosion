// The marker core (src/mesh.js, L0_MESH_SPEC.md Stage 1).
// Part 1, the spherical triangulation: build, delete, reinsert, points on edges, small motions; every result must
// pass the full check (orientation, adjacency, Euler's formula, Delaunay) and its triangles must cover 4 pi.
// Part 2, the world: every step, the continental volume ledger closes, the mesh is valid, no edge joins two plates'
// interior markers, ids are unique; at the end, determinism (two runs, and snapshot -> restore -> replay) and the
// area budget. Options: --n=32 (grid), --t=150 (Myr), --seed=genesis.
const crypto = require("crypto");
const { loadSim, loadMesh, arg } = require("./lib");
const SIM = loadSim(), MESH = loadMesh(SIM);
let fails = 0;
const ok = (c, msg) => { console.log((c ? "ok   " : "FAIL ") + msg); if (!c) fails++; };

// ---------- part 1 ----------
{
  const N = 20000; let s = 7; const r = () => { s = (s * 16807) % 2147483647; return s / 2147483647; };
  const T = new MESH.SphereTri(N + 12000, 2 * N + 10), vs = [];
  for (let i = 0; i < N; i++) { const z = 2 * r() - 1, p = 2 * Math.PI * r(), q = Math.sqrt(1 - z * z); T.X[i] = q * Math.cos(p); T.Y[i] = q * Math.sin(p); T.Z[i] = z; vs.push(i); }
  vs.sort((a, b) => { const k = v => Math.floor((Math.atan2(T.Y[v], T.X[v]) + 4) * 20) * 1000 + Math.floor((T.Z[v] + 1) * 20); return k(a) - k(b); });
  T.build(vs); let n = N, c = T.check(n);
  let A = 0; for (let t = 0; t < T.nT; t++) if (T.tAlive[t]) A += T.area(t);
  ok(c.ok && Math.abs(A / (4 * Math.PI) - 1) < 1e-12, `build: ${N} points, ${c.F} triangles, valid, area/4pi - 1 = ${(A / (4 * Math.PI) - 1).toExponential(1)}`);
  for (let i = 0; i < N; i += 5) { T.remove(i); n--; } c = T.check(n); ok(c.ok, `delete every fifth point: valid ${c.errs.join(";")}`);
  for (let i = 0; i < N; i += 5) if (T.insert(i, T.vt[i + 1])) n++; c = T.check(n); ok(c.ok, `reinsert them: valid ${c.errs.join(";")}`);
  const edges = []; for (let t = 0; t < T.nT && edges.length < 8000; t += 7) if (T.tAlive[t]) edges.push(T.tv[3 * t], T.tv[3 * t + 1]);
  let v = N, ins = 0; for (let q = 0; q < edges.length; q += 2) { const a = edges[q], b = edges[q + 1]; let x = T.X[a] + T.X[b], y = T.Y[a] + T.Y[b], z = T.Z[a] + T.Z[b]; const l = Math.hypot(x, y, z); T.X[v] = x / l; T.Y[v] = y / l; T.Z[v] = z / l; if (T.insert(v, T.vt[a])) { ins++; n++; } v++; }
  c = T.check(n); ok(c.ok && c.minOrient > 0, `insert ${ins} exact edge midpoints (edge splits): valid, smallest orientation ${c.minOrient.toExponential(1)}`);
  // a small motion that turns no triangle over, then flips
  for (let i = 0; i < v; i++) { if (T.vt[i] < 0) continue; const e = 2e-5; let x = T.X[i] + e * (r() - .5), y = T.Y[i] + e * (r() - .5), z = T.Z[i] + e * (r() - .5); const l = Math.hypot(x, y, z); T.X[i] = x / l; T.Y[i] = y / l; T.Z[i] = z / l; }
  const pre = T.check(n, false), q = []; for (let t = 0; t < T.nT; t++) if (T.tAlive[t]) q.push(t);
  const flips = T.legalize(q); c = T.check(n);
  ok(!pre.ok || c.ok, `small motion (${pre.ok ? "no triangle turned over" : "skipped: a triangle turned over"}), ${flips} flips: Delaunay again ${c.errs.join(";")}`);
}

// ---------- part 2 ----------
const n = +arg("n", 32), Tend = +arg("t", 150), seed = arg("seed", "genesis");
const opts = { n, seed, nPlates: 10, mode: "noise", nnr: true, contTarget: 0.4 };
const hash = w => { const h = crypto.createHash("sha1"), k = w.nSlots;
  for (const a of [w.tri.X, w.tri.Y, w.tri.Z, w.mT, w.mA, w.mPlate, w.mKind, w.mB]) h.update(Buffer.from(a.buffer, a.byteOffset, a.BYTES_PER_ELEMENT * k));
  h.update(String(w.t)); return h.digest("hex").slice(0, 12); };
const contV = w => { let v = 0; for (let m = 0; m < w.nSlots; m++) if (w.mKind[m] === 1 && w.mK[m]) v += w.mA[m] * w.mT[m]; return v; };
{
  const t0 = Date.now(), w = new MESH.MeshWorld(opts);
  let worstLed = 0, worstClosure = 0, meshBad = 0, crossMax = 0, dup = 0, snap = null, snapHashAt = null, steps = 0;
  while (w.t < Tend) {
    const v0 = contV(w), l0 = { ...w.ledger }; w.step(); steps++;
    let s = 0; for (const k in w.ledger) if (k.startsWith("cont")) s += w.ledger[k] - (l0[k] || 0);
    worstLed = Math.max(worstLed, Math.abs(contV(w) - v0 - s) / v0);
    if (steps % 5 === 0) { const a = w.audit(); if (!a.mesh.ok) meshBad++; crossMax = Math.max(crossMax, a.crossEdges); dup += a.dupIds; worstClosure = Math.max(worstClosure, Math.abs(a.closure)); }
    if (!snap && w.t >= Tend / 2) { snap = w.snapshot(); snapHashAt = w.steps; }
  }
  const a = w.audit(), end = hash(w);
  console.log(`world n=${n} seed=${seed}: ${w.t.toFixed(0)} Myr in ${steps} steps, ${((Date.now() - t0) / steps).toFixed(0)} ms/step; markers ${a.markers}; ${JSON.stringify(w.mStats)}`);
  ok(worstLed < 1e-8, `continental volume ledger closes every step (worst ${worstLed.toExponential(1)} of the total)`);
  ok(meshBad === 0, `mesh valid at every check (orientation, adjacency, Euler)`);
  ok(crossMax === 0, `no edge joins two plates' interior markers (most seen: ${crossMax})`);
  ok(dup === 0, `marker ids unique`);
  ok(worstClosure < 0.02, `area budget: interior markers' area within 2% of the sphere at every check (worst ${(100 * worstClosure).toFixed(2)}%; oceanic areas are geometric, so this is continental crust's carried area against the space it covers)`);
  const w2 = new MESH.MeshWorld(opts); while (w2.t < Tend) w2.step();
  ok(hash(w2) === end, `determinism: a second run ends in the same state (${end})`);
  const w3 = new MESH.MeshWorld(opts); w3.restore(snap); while (w3.t < Tend) w3.step();
  ok(hash(w3) === end, `restore at step ${snapHashAt} and replay ends in the same state`);
}
console.log(fails ? `\n${fails} check(s) failed` : "\nall mesh checks passed"); process.exitCode = fails ? 1 : 0;
