// The marker core (src/mesh.js, L0_MESH_SPEC.md Stage 1).
// Part 1, the spherical triangulation: build, delete, reinsert, points on edges, small motions; every result must
// pass the full check (orientation, adjacency, Euler's formula, Delaunay) and its triangles must cover 4 pi.
// Part 2, the world: every step, the continental volume ledger closes, the mesh is valid, no edge joins two plates'
// interior markers, ids are unique; at the end, determinism (two runs, and snapshot -> restore -> replay) and the
// area budget. Part 3, the worker protocol (MESH.workerMain against a fake worker scope).
// Options: --n=32 (grid), --t=150 (Myr), --seed=genesis.
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
// ---------- part 2b: the forces of Stage 2 ----------
{
  const { gpe, PH } = MESH, rp = gpe(0, 7, 0) - gpe(0, 7, 100);
  ok(rp > 3.5e12 && rp < 4.5e12, `ridge push at 100 Myr from the GPE of floated columns: ${rp.toExponential(2)} N/m (Turcotte & Schubert: about 3.9e12)`);
  let mono = true; for (let a = 0; a < 200; a += 5) if (gpe(0, 7, a + 5) > gpe(0, 7, a) + 1) mono = false;
  ok(mono, `GPE falls with ocean age (ridges push)`);
  // uniform ocean, no slabs: every torque cancels, so nothing moves
  const w = new MESH.MeshWorld({ ...opts, n: 16, seed: "still", mantle: false });
  for (let m = 0; m < w.nSlots; m++) if (w.mKind[m] === 1) { w.mK[m] = 0; w.mT[m] = 7; w.mB[m] = -50; w.sA[m] = 0; }
  w.checkSolve = true; w.solveVelocities(); const vmax = Math.max(...w.plates.filter(P => P.alive).map(P => Math.hypot(...P.W))) * 6371 / 10;
  ok(vmax < 1e-9, `uniform ocean with no slabs: no plate moves (fastest ${vmax.toExponential(1)} cm/yr)`);
  // mantle flow (§6.1): basis fields orthonormal, so the tabulated RMS equals the root sum of squared coefficients;
  // a uniform ocean with no slabs moves only by the mantle's drag, so no faster than the mantle
  const wm = new MESH.MeshWorld({ ...opts, n: 16, seed: "drift" }), expect = Math.sqrt(wm.mc.reduce((a, b) => a + b * b, 0));
  ok(Math.abs(wm.mantleRmsNow / expect - 1) < 0.02, `mantle flow: RMS on the grid ${wm.mantleRmsNow.toFixed(3)} cm/yr against ${expect.toFixed(3)} from its coefficients`);
  for (let m = 0; m < wm.nSlots; m++) if (wm.mKind[m] === 1) { wm.mK[m] = 0; wm.mT[m] = 7; wm.mB[m] = -50; wm.sA[m] = 0; }
  wm.solveVelocities(); wm.updateStats(true);
  ok(wm.rmsSpeed > 0.05 && wm.rmsSpeed <= wm.mantleRmsNow, `uniform ocean over flowing mantle: plates drift at ${wm.rmsSpeed.toFixed(2)} cm/yr, under the mantle's ${wm.mantleRmsNow.toFixed(2)}`);
  // the experimental option: subduction starting at old passive margins keeps the ledger, the mesh and the area budget
  { const wi = new MESH.MeshWorld({ ...opts, seed: "margins", initiation: true }); let worst = 0, bad = 0;
    const contV = w => { let v = 0; for (let m = 0; m < w.nSlots; m++) if (w.mKind[m] === 1 && w.mK[m]) v += w.mA[m] * w.mT[m]; return v; };
    while (wi.t < 80) { const v0 = contV(wi), l0 = { ...wi.ledger }; wi.step(); let d = 0; for (const k in wi.ledger) if (k.startsWith("cont")) d += wi.ledger[k] - (l0[k] || 0);
      worst = Math.max(worst, Math.abs(contV(wi) - v0 - d) / v0); if (wi.steps % 10 === 0) { const a = wi.audit(); if (!a.mesh.ok || Math.abs(a.closure) > 0.02) bad++; } }
    ok(worst < 1e-8 && bad === 0 && (wi.mStats.initiations || 0) > 0, `subduction starting at old margins (option): ${wi.mStats.initiations || 0} starts in 80 Myr; ledger worst ${worst.toExponential(1)}; mesh and area budget held`); }
  // stress recovery (§7): converged, tangent to the sphere, zero without loads; damage within [0, 1]
  { const ws = new MESH.MeshWorld({ ...opts, seed: "stress" }); while (ws.t < 30) ws.step(); ws._loadsDue = true; ws.solveVelocities(); ws.stressStep(); // checked as solved (markers move after)
    let tan = 0, mag = 0, dOk = true;
    for (let m = 0; m < ws.nSlots; m++) { if (ws.mKind[m] !== 1) continue; const S = ws.mSig, o = 6 * m, x = ws.tri.X[m], y = ws.tri.Y[m], z = ws.tri.Z[m];
      const v = [S[o] * x + S[o + 3] * y + S[o + 4] * z, S[o + 3] * x + S[o + 1] * y + S[o + 5] * z, S[o + 4] * x + S[o + 5] * y + S[o + 2] * z];
      tan = Math.max(tan, Math.hypot(...v)); mag = Math.max(mag, Math.abs(S[o]) + Math.abs(S[o + 1]) + Math.abs(S[o + 2])); if (!(ws.mD[m] >= 0 && ws.mD[m] <= 1)) dOk = false; }
    const si = ws.stressInfo || {};
    ok(si.residual <= 1.001 * MESH.PH.stressTol && tan < 1e-6 * mag && dOk, `stress: solved to ${(si.residual || 0).toExponential(1)} in ${si.iters} iterations; normal component ${(tan / mag).toExponential(1)} of the largest; damage within [0, 1]`);
    const w0 = new MESH.MeshWorld({ ...opts, n: 16, seed: "still", mantle: false });
    for (let m = 0; m < w0.nSlots; m++) if (w0.mKind[m] === 1) { w0.mK[m] = 0; w0.mT[m] = 7; w0.mB[m] = -50; w0.sA[m] = 0; }
    for (const P of w0.plates) P.W.fill(0); w0._loadsDue = true; w0.solveVelocities(); w0.stressStep(); let mx = 0; for (let q = 0; q < 6 * w0.nSlots; q++) mx = Math.max(mx, Math.abs(w0.mSig[q]));
    ok(mx < 1, `uniform ocean, no slabs, no mantle: no stress (largest ${mx.toExponential(1)} N/m)`); }
  // cracks keep the ledger, the mesh and the area budget (weak plates, so they crack within the run)
  { const keep = MESH.PH.strScale; MESH.PH.strScale = 1; const wc = new MESH.MeshWorld({ ...opts, seed: "cracks" }); let worst = 0, bad = 0;
    const contV = w => { let v = 0; for (let m = 0; m < w.nSlots; m++) if (w.mKind[m] === 1 && w.mK[m]) v += w.mA[m] * w.mT[m]; return v; };
    while (wc.t < 80) { const v0 = contV(wc), l0 = { ...wc.ledger }; wc.step(); let d = 0; for (const k in wc.ledger) if (k.startsWith("cont")) d += wc.ledger[k] - (l0[k] || 0);
      worst = Math.max(worst, Math.abs(contV(wc) - v0 - d) / v0); if (wc.steps % 10 === 0) { const a = wc.audit(); if (!a.mesh.ok || Math.abs(a.closure) > 0.02) bad++; } }
    MESH.PH.strScale = keep;
    ok(worst < 1e-8 && bad === 0 && (wc.mStats.cracks || 0) > 0, `cracks: ${wc.mStats.cracks || 0} in 80 Myr on weak plates; ledger worst ${worst.toExponential(1)}; mesh and area budget held`); }
  const w2 = new MESH.MeshWorld({ ...opts, seed: "balance" }); w2.checkSolve = true; let worst = 0;
  while (w2.t < 60) { w2.step(); worst = Math.max(worst, w2.solveResidual); }
  ok(worst < 1e-6, `torque balance (I7): |K omega - b| / |b| below 1e-6 at every step to 60 Myr (worst ${worst.toExponential(1)})`);
}
// ---------- part 3: the worker protocol ----------
// MESH.workerMain runs against a fake worker scope: init, run to 20 Myr, one more step, ask
// for a column and the handoff; the views it posts must match a MeshWorld stepped directly.
(async () => {
  const posted = [], fake = { postMessage: m => posted.push(m), onmessage: null };
  const main = new Function("self", "MESH", "MessageChannel", "performance", "(" + MESH.workerMain + ")(); return self;");
  // a stand-in for MessageChannel on setImmediate (Node's real one would keep the process alive after the test)
  class Chan { constructor() { this.port1 = {}; this.port2 = { postMessage: () => setImmediate(() => this.port1.onmessage && this.port1.onmessage()) }; } }
  const self = main(fake, MESH, Chan, performance), send = m => self.onmessage({ data: m });
  const until = (type, test = () => true) => new Promise(res => { const poll = () => { const m = posted.find(x => x.type === type && test(x)); if (m) res(m); else setTimeout(poll, 5); }; poll(); });
  const o = { ...opts, seed: "worker" };
  send({ type: "init", opts: o }); await until("ready");
  send({ type: "run", until: 20 }); const st = await until("stopped", m => m.t >= 20 - 1e-9);
  send({ type: "step" }); const st2 = await until("stopped", m => m.steps === st.steps + 1);
  send({ type: "column", i: 123 }); const col = await until("column");
  send({ type: "extract" }); const ex = await until("extract");
  const w = new MESH.MeshWorld(o); while (w.steps < st2.steps) w.step();
  let same = true; for (const k of ["owner", "elev", "wK", "wT"]) for (let i = 0; i < w.N; i++) if (st2.view[k][i] !== w[k][i]) { same = false; break; }
  ok(same && st2.t === w.t, `worker protocol: the view after ${st2.steps} steps (run to 20 Myr, then one step) matches a direct run`);
  ok(col.units.length > 0 && Array.isArray(col.eps), `worker answers a column request (${col.units.length} units) and a handoff request (${ex.handoff.faults.length} faults)`);
  console.log(fails ? `\n${fails} check(s) failed` : "\nall mesh checks passed"); process.exitCode = fails ? 1 : 0;
})();

