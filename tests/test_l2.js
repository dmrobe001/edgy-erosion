// L2 checks, on a landmass refined at L1 (as in test_l1.js): the relief L2 starts from is
// 1. a pure function of position: a point alone equals the same point in a grid, and two builds agree;
// 2. L1's at L1's scale: averaged over each 10 km cell, it reproduces L1's surface (the added relief owns the band
//    below 20 km);
// 3. as loud as designed where the band is L2's alone: the octaves of 4 km and less, RMS per unit of L1 roughness;
// 4. stretched along strike where L1 shortened the crust: gradients across strike over along strike, about a^2
//    (on the landmass if it has shortening, and on a made-up plain under uniform shortening);
// 5. free of dams on L1's trunk rivers: the lowest pass between two points on a river, along any path at 500 m, is no
//    higher than the river's level at the upstream point;
// 6. basin outlines burnt in by strips agree with routing the whole region at 500 m;
// 7. two neighbouring selections, refined independently, agree on their shared outline;
// 8. a worker built from the module's source makes the same tiles as the main thread.
const { loadSim, loadL1, loadL2, arg } = require("./lib");
const SIM = loadSim(), L1 = loadL1(SIM), L2 = loadL2(SIM, L1);
let ok = true; const check = (cond, msg) => { if (!cond) { ok = false; console.log("FAIL " + msg); } };

const n = +arg("n", 64), W = +arg("win", 10), T0 = +arg("t0", 400), H = +arg("h", 70), target = +arg("land", 200);
const w = new SIM.World({ n, seed: "genesis", nPlates: 10, mode: "scatter", contTarget: 0.5 }), snaps = [];
while (w.t < T0) { const last = snaps[snaps.length - 1]; if (!last || w.t - last.t >= 10 - 1e-9) { snaps.push({ t: w.t, steps: w.steps, snap: w.snapshot() }); while (snaps.length > 4) snaps.shift(); } w.step(); }
let cell = -1, bd = Infinity; { const seen = new Uint8Array(w.N);
  for (let i = 0; i < w.N; i++) { if (seen[i] || w.elev[i] <= 0) continue; const lm = L1.landmassAt(w, i); for (const c of lm.cells) seen[c] = 1; const d = Math.abs(lm.land - target); if (d < bd) { bd = d; cell = i; } } }
const pick = L1.landmassAt(w, cell), entry = [...snaps].reverse().find(s => s.t <= w.t - W - 1e-9);
const rp = new L1.Replay(w, entry, W, pick); while (!rp.advance(1e9)); const plan = rp.finish();
w.restore(entry.snap); w.trackParcels(); while (w.steps < plan.t0Steps) w.step();
const sh = new L1.Sheet(w, plan, { h: H }); while (w.steps < plan.endSteps) sh.advance(); w.untrackParcels();
const B = L1.basins(sh), { nx, ny, h } = B;
let t0 = Date.now(); const R = L2.relief(sh, B), buildMs = Date.now() - t0;
const land = c => B.land[c] && B.e[c] > 0, cx = c => B.x0 + (c % nx + .5) * h, cy = c => B.y0 + (Math.floor(c / nx) + .5) * h;
console.log(`landmass: ${B.landCells} land cells of ${h} km; relief built in ${buildMs} ms, ${R.nSegments} trunk reaches, octaves ${R.octaves.join(", ")} km`);

// 1. pure function of position
{ const R2 = L2.relief(sh, B), g = R.grid(cx(0) + 3000 % (nx * h), B.y0 + ny * h / 2, 0.5, 40, 30); let same = true, alone = true;
  for (let j = 0; j < 30; j++) for (let i = 0; i < 40; i++) { const x = cx(0) + 3000 % (nx * h) + (i + .5) * 0.5, y = B.y0 + ny * h / 2 + (j + .5) * 0.5, v = R.z(x, y, 1);
    if (Math.fround(v) !== g[j * 40 + i]) alone = false; if (R2.z(x, y, 1) !== v) same = false; }
  check(alone && same, "relief: a point alone equals the grid's value there, and a second build agrees"); }

// 2. and 3. at L1's scale, and the added relief's RMS; on land cells away from the coast (whole 10 km cells at 500 m)
{ const k = 20, dx = h / k; let n2 = 0, d2 = 0, e2 = 0, a2 = 0, s2 = 0, t = Date.now(), cells = 0;
  for (let c = 0; c < nx * ny; c += 7) { if (!land(c)) continue; const i = c % nx, j = Math.floor(c / nx); let inner = true;
    for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) { const q = (j + dj) * nx + i + di; if (i + di < 0 || j + dj < 0 || i + di >= nx || j + dj >= ny || !land(q)) inner = false; }
    if (!inner) continue; let s = 0; const x0 = cx(c) - h / 2, y0 = cy(c) - h / 2;
    for (let b = 0; b < k; b++) for (let a = 0; a < k; a++) { const x = x0 + (a + .5) * dx, y = y0 + (b + .5) * dx, z = R.z(x, y, 1), r1 = R.roughAt(x, y), add = R.added(x, y, 1), short = R.band(x, y, 1) - R.band(x, y, 8);
      s += z; a2 += add * add; s2 += short * short; cells++; }
    const m = s / (k * k), d = m - B.e[c]; d2 += d * d; n2++; }
  const ms = (Date.now() - t) / cells * 1e3, sigma = Math.sqrt(R.octaves.filter(L => L <= 4).reduce((a, L) => a + Math.pow(L / R.o.refKm, 2 * R.o.H), 0)), got = Math.sqrt(s2 / cells);
  // how much L1's own relief varies within reach: RMS of L1's surface less its neighbours' mean
  for (let c = 0; c < nx * ny; c += 7) { if (!land(c)) continue; const i = c % nx, j = Math.floor(c / nx); if (i < 1 || j < 1 || i >= nx - 1 || j >= ny - 1) continue; const m = (B.e[c - 1] + B.e[c + 1] + B.e[c - nx] + B.e[c + nx]) / 4; e2 += (B.e[c] - m) ** 2; }
  const rmsD = Math.sqrt(d2 / n2), rmsL1 = Math.sqrt(e2 / n2), rmsAdd = Math.sqrt(a2 / cells);
  console.log(`L1's scale: ${n2} cells averaged at 500 m differ from L1's surface by ${(rmsD * 1000).toFixed(1)} m RMS (L1's cell-to-cell relief ${(rmsL1 * 1000).toFixed(0)} m RMS; the added relief ${(rmsAdd * 1000).toFixed(0)} m RMS); ${ms.toFixed(2)} us per point`);
  console.log(`added relief: octaves of 4 km and less ${got.toFixed(3)} RMS per unit of L1 roughness (designed ${sigma.toFixed(3)})`);
  check(rmsD < 0.25 * rmsL1, "relief: averaged over L1's cells it reproduces L1's surface");
  check(Math.abs(got / sigma - 1) < 0.1, "relief: the octaves of 4 km and less have their designed RMS"); }

// 4. anisotropy: where the lattice stretches most, the band's gradients across strike over along strike
{ let best = null; for (let c = 0; c < nx * ny; c += 3) { if (!land(c)) continue; const s = R.strikeAt(cx(c), cy(c)); if (s && (!best || s.stretch > best.s.stretch)) best = { c, s }; }
  if (!best || best.s.stretch < 1.05) console.log("anisotropy: no shortening on this landmass to test");
  else { const { c, s } = best, ux = Math.cos(s.strike), uy = Math.sin(s.strike), x0 = cx(c), y0 = cy(c), dd = 0.05; let ga = 0, gc = 0;
    // within a few km of the node, where its own frame dominates the blend
    for (let b = -20; b <= 20; b++) for (let a = -20; a <= 20; a++) { const x = x0 + a * 0.2, y = y0 + b * 0.2, f0 = R.band(x, y, 1);
      ga += ((R.band(x + dd * ux, y + dd * uy, 1) - f0) / dd) ** 2; gc += ((R.band(x - dd * uy, y + dd * ux, 1) - f0) / dd) ** 2; }
    const ratio = gc / ga; console.log(`anisotropy: stretch ${s.stretch.toFixed(2)} along strike; gradient variance across over along strike ${ratio.toFixed(2)} (a^2 = ${(s.stretch ** 2).toFixed(2)})`);
    check(ratio > 1 + 0.5 * (s.stretch ** 2 - 1), "relief: stretched along strike where L1 shortened the crust"); } }

// 4b. the same on a made-up plain (no rivers) under a uniform 30% shortening along N30E: strike N60W
{ const nx2 = 60, ny2 = 60, th = Math.PI / 3, sk = { world: { opts: { seed: "aniso" } }, nx: 1, ny: 1, strainField: () => ({ lmin: [0.7], dir: [th] }), sample: (f) => f[0] };
  const B2 = { nx: nx2, ny: ny2, h: 10, x0: 0, y0: 0, e: new Float64Array(nx2 * ny2).map((_, c) => 0.5 + 0.05 * Math.sin(c * 0.37)), land: new Uint8Array(nx2 * ny2).fill(1),
    rec: Int32Array.from({ length: nx2 * ny2 }, (_, c) => c), A: new Float64Array(nx2 * ny2).fill(100), S: { f: new Float64Array(nx2 * ny2) } };
  const R2 = L2.relief(sk, B2), st = R2.strikeAt(300, 300), ux = Math.cos(st.strike), uy = Math.sin(st.strike), dd = 0.05; let ga = 0, gc = 0;
  for (let b = 0; b < 60; b++) for (let a = 0; a < 60; a++) { const x = 200 + a * 3.1, y = 200 + b * 3.1, f0 = R2.band(x, y, 1);
    ga += ((R2.band(x + dd * ux, y + dd * uy, 1) - f0) / dd) ** 2; gc += ((R2.band(x - dd * uy, y + dd * ux, 1) - f0) / dd) ** 2; }
  const ratio = gc / ga, strikeDeg = ((st.strike * 180 / Math.PI) % 180 + 180) % 180;
  console.log(`anisotropy, made-up: stretch ${st.stretch.toFixed(2)}, strike ${strikeDeg.toFixed(0)} deg from east (expected 150); gradient variance across over along strike ${ratio.toFixed(2)} (a^2 = ${(st.stretch ** 2).toFixed(2)})`);
  check(Math.abs(strikeDeg - 150) < 1e-6 && Math.abs(ratio / st.stretch ** 2 - 1) < 0.25, "relief: stretched along strike under uniform shortening"); }

// 5. no dams on trunk rivers: chains of five reaches; at each end the point on the river's trench floor near L1's cell
//    centre (the river wanders within 3 warp RMS of it); the minimax path between them on a 500 m grid over the
//    chain's box, against the river's level at the upstream end
{ const dx = 0.5, heads = []; for (let c = 0; c < nx * ny; c += 11) if (land(c) && B.A[c] >= R.o.burnKm2 * 3) heads.push(c);
  let tested = 0, worst = -Infinity, t = Date.now();
  for (const c0 of heads.slice(0, 40)) { const chain = [c0]; for (let k = 0; k < 5; k++) { const r = B.rec[chain[chain.length - 1]]; if (r === chain[chain.length - 1] || !land(r)) break; chain.push(r); }
    if (chain.length < 4) continue; const pad = 6, xs = chain.map(cx), ys = chain.map(cy), X0 = Math.min(...xs) - pad, Y0 = Math.min(...ys) - pad, gx = Math.ceil((Math.max(...xs) + pad - X0) / dx), gy = Math.ceil((Math.max(...ys) + pad - Y0) / dx);
    // each end: the point near the reach's cell centre lowest on the river's own trench (on its floor)
    const g = R.grid(X0, Y0, dx, gx, gy, 1), tr = new Float64Array(gx * gy); for (let j = 0; j < gy; j++) for (let i = 0; i < gx; i++) tr[j * gx + i] = R.trench(X0 + (i + .5) * dx, Y0 + (j + .5) * dx);
    const low = (c) => { let b = -1; for (let j = 0; j < gy; j++) for (let i = 0; i < gx; i++) { const x = X0 + (i + .5) * dx, y = Y0 + (j + .5) * dx; if (Math.hypot(x - cx(c), y - cy(c)) > 3 * R.o.warpKm + 1) continue; const q = j * gx + i; if (b < 0 || tr[q] < tr[b]) b = q; } return b; };
    const s = low(chain[0]), e = low(chain[chain.length - 1]);
    // minimax (bottleneck) from s: the lowest possible highest point on a path to each cell
    const best = new Float64Array(gx * gy).fill(Infinity), done = new Uint8Array(gx * gy), pq = new L1.Heap(9 * gx * gy); best[s] = g[s]; pq.push(s, g[s]);
    while (pq.n) { const q = pq.pop(); if (done[q]) continue; done[q] = 1; if (q === e) break; const bh = best[q], i = q % gx, j = Math.floor(q / gx);
      for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) { const I = i + di, J = j + dj; if ((!di && !dj) || I < 0 || J < 0 || I >= gx || J >= gy) continue; const p = J * gx + I, nb = Math.max(bh, g[p]); if (nb < best[p]) { best[p] = nb; pq.push(p, nb); } } }
    worst = Math.max(worst, best[e] - Math.max(g[s], tr[s])); tested++; }
  console.log(`trunk rivers: ${tested} chains of up to 5 reaches; highest pass above the upstream end ${(worst * 1000).toFixed(1)} m (${Date.now() - t} ms)`);
  check(tested > 5 && worst <= 1e-6, "relief: no dams on L1's trunk rivers"); }

// 6. the outline burnt in by strips against routing the whole box at 500 m (the box's edge 40 km out, labelled by L1):
//    mean displacement (land assigned differently, over the outline's length), for a cluster tapped on the coast and
//    for a mid-size basin; and L1's own 10 km outline for comparison
const sizes = new Map(); for (let c = 0; c < nx * ny; c++) if (B.land[c]) sizes.set(B.basin[c], (sizes.get(B.basin[c]) || 0) + 1);
const pj = plan.proj, at = c => [cx(c), cy(c)];
let coastC = -1; for (let c = Math.floor(nx * ny / 2); c < nx * ny; c++) if (land(c) && !B.land[c + 1]) { coastC = c; break; }
let midB = -1; for (const [b, n] of sizes) if (n >= 150 && n <= 400) { midB = b; break; }
let midC = 0; while (!(B.land[midC] && B.basin[midC] === midB)) midC++;
const selA = L1.basinSelection(B, pj, ...at(coastC), 150), selM = L1.basinSelection(B, pj, ...at(midC), 0);
for (const [name, sel] of [["coastal cluster", selA], ["mid-size basin", selM]]) {
  let i0 = 1e9, i1 = -1, j0 = 1e9, j1 = -1; for (let c = 0; c < nx * ny; c++) if (sel.mask[c]) { const i = c % nx, j = Math.floor(c / nx); i0 = Math.min(i0, i); i1 = Math.max(i1, i); j0 = Math.min(j0, j); j1 = Math.max(j1, j); }
  const m = 40, dx = 0.5, X0 = B.x0 + i0 * h - m, Y0 = B.y0 + j0 * h - m, nX = Math.round(((i1 - i0 + 1) * h + 2 * m) / dx), nY = Math.round(((j1 - j0 + 1) * h + 2 * m) / dx);
  const bd = L2.boundary(R, B), t = Date.now(), ref = bd.refine(sel), ms = Date.now() - t, F = bd.full(X0, Y0, nX, nY), inS = new Set(sel.basins);
  let bad = 0, badL1 = 0, len = 0; for (let j = 0; j < nY; j++) for (let i = 0; i < nX; i++) { const q = j * nX + i; if (F.z[q] <= 0) continue; const x = X0 + (i + .5) * dx, y = Y0 + (j + .5) * dx, tr = inS.has(F.lab[q]);
    if (ref.inside(x, y) !== tr) bad += dx * dx; if (inS.has(bd.lab1[Math.floor((y - B.y0) / h) * nx + Math.floor((x - B.x0) / h)]) !== tr) badL1 += dx * dx;
    if (i + 1 < nX && F.z[q + 1] > 0 && tr !== inS.has(F.lab[q + 1])) len += dx; if (j + 1 < nY && F.z[q + nX] > 0 && tr !== inS.has(F.lab[q + nX])) len += dx; }
  console.log(`burn-in, ${name}: ${sel.basins.length} basins, ${Math.round(ref.areaKm2).toLocaleString("en-US")} km² (L1: ${Math.round(ref.areaL1Km2).toLocaleString("en-US")}); ${ref.nTiles} tiles stitched (${ref.rounds} rounds, ${ref.guessCells} L1 cells with guesses), ${ms} ms; mean displacement from full routing ${(bad / len).toFixed(3)} km (L1's outline ${(badL1 / len).toFixed(2)} km)`);
  check(bad / len < 0.25 && bad < 0.25 * badL1, "burn-in: the strips' outline is within half an L2 cell of full routing, " + name); }

// 7. seams: two neighbouring selections, refined independently, agree wherever either puts a point in one of the two,
//    and no point is in both
{ const nbr = new Map(); for (let c = 0; c < nx * ny; c++) { if (!B.land[c] || B.basin[c] !== midB) continue; for (const d of [1, -1, nx, -nx]) { const q = c + d; if (B.land[q] && B.basin[q] !== midB) nbr.set(B.basin[q], (nbr.get(B.basin[q]) || 0) + 1); } }
  const other = [...nbr].sort((a, b) => b[1] - a[1])[0][0]; let oc = 0; while (!(B.land[oc] && B.basin[oc] === other)) oc++;
  const selO = L1.basinSelection(B, pj, ...at(oc), 0), r1 = L2.boundary(R, B).refine(selM), r2 = L2.boundary(R, B).refine(selO);
  // where either refinement puts a point in one of the two selections, the other must put it in the same one (labels of
  // third basins may differ: each stitched domain ends somewhere else)
  const AB = new Set([...selM.basins, ...selO.basins]); let shared = 0, same = true, both = 0, differ = 0;
  for (const [c, L] of r1.cells) { const M = r2.cells.get(c); if (!M) continue; shared++; for (let q = 0; q < L.length; q++) if ((AB.has(L[q]) || AB.has(M[q])) && L[q] !== M[q]) { same = false; differ++; } }
  for (const [c] of r1.cells) for (let v = 0; v < 20; v++) for (let u = 0; u < 20; u++) { const x = B.x0 + (c % nx) * h + (u + .5) * 0.5, y = B.y0 + Math.floor(c / nx) * h + (v + .5) * 0.5; if (r1.inside(x, y) && r2.inside(x, y)) both++; }
  console.log(`seams: basins ${midB} and ${other} refined apart share ${shared} L1 cells; labels of the two basins there ${same ? "identical" : "DIFFER at " + differ + " points"}; ${both} points in both`);
  // stitched domains depend on the selection, so a point balanced between two exits can resolve differently: allowed at
  // 1 in 10,000 shared points
  check(shared > 0 && differ <= 1e-4 * shared * 400 && both === 0, "burn-in: neighbouring selections agree on their shared outline"); }

// 8. a worker built as the page builds it (the module's source and the relief's plain inputs) makes the same tiles
//    as the main thread, bit for bit
(async () => {
  const fs = require("fs"), path = require("path"), { Worker } = require("worker_threads"), src2 = fs.readFileSync(path.join(__dirname, "..", "src", "l2.js"), "utf8");
  const mod = src2.slice(src2.indexOf("function L2module("), src2.indexOf("const L2=L2module(SIM,L1);"));
  const code = `const {parentPort}=require("worker_threads");const self={postMessage:(m,t)=>parentPort.postMessage(m,t)};parentPort.on("message",d=>self.onmessage({data:d}));
const SIM={hashStr:${SIM.hashStr}};const L1={Heap:${L1.Heap}};${mod}\nconst L2=L2module(SIM,L1);let M=null;
self.onmessage=e=>{const m=e.data;if(m.init){M=L2.tileMaker(L2.reliefFrom(m.init.D),m.init.B,m.init.o);return;}const g=M.make(m.TI,m.TJ);self.postMessage({key:m.key,g});};`;
  const bd = L2.boundary(R, B), wi = bd.workerInputs(), D = L2.reliefInputs(sh, B), w = new Worker(code, { eval: true });
  w.postMessage({ init: { D, B: wi.B, o: wi.o } });
  let i0 = 0; while (!land(i0)) i0++; const TI = Math.floor((i0 % nx) / 4), TJ = Math.floor(Math.floor(i0 / nx) / 4), key = TI * bd.tny + TJ;
  const g = await new Promise(res => { w.on("message", m => res(m.g)); w.postMessage({ key, TI, TJ }); }); w.terminate();
  const own = bd.ensure(key); let same = g.z.length === own.z.length && g.eu.length === own.eu.length;
  for (const f of ["z", "out", "lab", "per", "eu", "ev", "ew"]) for (let q = 0; same && q < own[f].length; q++) if (own[f][q] !== g[f][q]) same = false;
  console.log(`worker: tile (${TI}, ${TJ}) from a worker ${same ? "identical to" : "DIFFERS from"} the main thread's (${own.eu.length} graph edges)`);
  check(same, "burn-in: a worker computes the same tiles as the main thread");
  console.log(ok ? "PASS L2" : "FAIL L2"); process.exitCode = ok ? 0 : 1; })();
