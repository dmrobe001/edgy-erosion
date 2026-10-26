// Experiment (README section 6, item 3; not part of run_all, about 2 minutes). Decision gate: does erodibility from rock type (and age, fault damage) change L1's basin map more than re-drawing the
// unresolved fine relief does? Adjusted Rand index (ARI) between basin labellings of the land.
const { loadSim, loadL1, arg } = require("./lib");
const SIM = loadSim(), L1 = loadL1(SIM), RK = SIM.ROCK;
const T = +arg("T", 1000), plate = +arg("plate", 44);
const w = new SIM.World({ n: 96, seed: "genesis", nPlates: 10, mode: "scatter", contTarget: 0.5 });
const snaps = []; let last = -1e9;
while (w.t < T) { if (w.t - last >= 10) { snaps.push({ t: w.t, steps: w.steps, snap: w.snapshot() }); last = w.t; while (snaps.length > 4) snaps.shift(); } w.step(); }
let cell = -1, best = 0; const seen = new Uint8Array(w.N);
for (let i = 0; i < w.N; i++) { if (seen[i] || w.elev[i] <= 0 || (plate >= 0 && w.owner[i] !== plate)) continue; const lm = L1.landmassAt(w, i); for (const c of lm.cells) seen[c] = 1; if (lm.land > best) { best = lm.land; cell = i; } }
const pick = L1.landmassAt(w, cell), entry = [...snaps].reverse().find(s => s.t <= w.t - 20 - 1e-9);
const rp = new L1.Replay(w, entry, 20, pick); while (!rp.advance(1e9)); const plan = rp.finish();
w.restore(entry.snap); w.trackParcels(); while (w.steps < plan.t0Steps) w.step();
const sh = new L1.Sheet(w, plan, {}); while (w.steps < plan.endSteps) sh.advance();
const { nx, ny } = sh, h = sh.o.h, nc = nx * ny;
// per-marker factors, log-averaged per sheet cell (weighted by marker volume)
const TYPE = { [RK.OCEAN]: 0.7, [RK.HOTSPOT]: 0.7, [RK.ARC]: 1.5, [RK.CRATON]: 0.5, [RK.OROGEN]: 0.7, [RK.SHELF]: 3, [RK.CLASTIC]: 3 };
const tmp = [0, 0, 0], b = [0, 0, 0], wv = [0, 0, 0];
const faultAt = (x, y) => { sh.proj.inv(x, y, b); const Q = sh.Qref; wv[0] = Q[0] * b[0] + Q[1] * b[1] + Q[2] * b[2]; wv[1] = Q[3] * b[0] + Q[4] * b[1] + Q[5] * b[2]; wv[2] = Q[6] * b[0] + Q[7] * b[1] + Q[8] * b[2]; const i = w.cube.lookup(wv[0], wv[1], wv[2]); const fi = w.faultCellInfo(i); return fi ? fi.slip : 0; };
function field(fn, pow = 1) { const s = new Float64Array(nc), ws = new Float64Array(nc);
  for (let k = 0; k < sh.nm; k++) { if (!sh.mIn[k]) continue; const i = Math.floor((sh.mx[k] - sh.x0) / h), j = Math.floor((sh.my[k] - sh.y0) / h), c = j * nx + i, col = sh.mCol[k]; if (!col.length) continue;
    const vol = sh.mS[k] * sh.mA[k]; s[c] += vol * Math.log(fn(col[col.length - 1], k)); ws[c] += vol; }
  const out = new Float64Array(nc); let g = 0, n = 0; for (let c = 0; c < nc; c++) if (ws[c] > 0 && sh.act[c] && sh.e[c] > 0) { g += s[c] / ws[c]; n++; } g /= n;
  for (let c = 0; c < nc; c++) out[c] = ws[c] > 0 ? Math.exp(pow * (s[c] / ws[c] - g)) : 1; return out; }
const slipC = new Float64Array(nc); for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if (sh.act[j * nx + i]) slipC[j * nx + i] = faultAt(sh.x0 + (i + .5) * h, sh.y0 + (j + .5) * h);
const typeK = field(u => TYPE[u[0]] || 1);
const ageF = u => (u[0] === RK.SHELF || u[0] === RK.CLASTIC) ? 1 + 9 * Math.exp(-(sh.world.t - u[1]) / 30) : 1;
const grounded = field((u, k) => (TYPE[u[0]] || 1) * ageF(u));
for (let c = 0; c < nc; c++) grounded[c] *= 1 + 4 * Math.min(1, slipC[c] / 300);
{ let g = 0, n = 0; for (let c = 0; c < nc; c++) if (sh.act[c] && sh.e[c] > 0) { g += Math.log(grounded[c]); n++; } g /= n; for (let c = 0; c < nc; c++) grounded[c] /= Math.exp(g); }
const strong = Float64Array.from(grounded, v => v * v);
const spread = f => { const v = []; for (let c = 0; c < nc; c++) if (sh.act[c] && sh.e[c] > 0) v.push(f[c]); v.sort((a, b) => a - b); return `K x${v[Math.floor(.05 * v.length)].toFixed(2)}..x${v[Math.floor(.95 * v.length)].toFixed(2)} (5-95%), max/min ${(v[v.length - 1] / v[0]).toFixed(0)}`; };
let faulted = 0, nl = 0; for (let c = 0; c < nc; c++) if (sh.act[c] && sh.e[c] > 0) { nl++; if (slipC[c] > 60) faulted++; }
console.log(`sheet cells on land with 60+ km of recorded fault slip: ${(100 * faulted / nl).toFixed(0)}%`);
function ari(A, B, mask) { const nij = new Map(), ai = new Map(), bj = new Map(); let n = 0;
  for (let c = 0; c < A.length; c++) { if (!mask[c]) continue; n++; const k = A[c] * 4194304 + B[c]; nij.set(k, (nij.get(k) || 0) + 1); ai.set(A[c], (ai.get(A[c]) || 0) + 1); bj.set(B[c], (bj.get(B[c]) || 0) + 1); }
  const C2 = x => x * (x - 1) / 2; let sij = 0, sa = 0, sb = 0; for (const v of nij.values()) sij += C2(v); for (const v of ai.values()) sa += C2(v); for (const v of bj.values()) sb += C2(v);
  const ex = sa * sb / C2(n), mx = (sa + sb) / 2; return (sij - ex) / (mx - ex); }
const runs = { base: {}, reseed: { salt: "b" }, reseed2: { salt: "c" }, noclimate: { climate: false }, type: { Kmul: typeK }, grounded: { Kmul: grounded }, strong: { Kmul: strong } }, res = {};
for (const [k, o] of Object.entries(runs)) { const t0 = Date.now(); res[k] = L1.basins(sh, o); console.log(`${k}: ${o.Kmul ? spread(o.Kmul) : "uniform K"}; ${res[k].nBasins} basins; ${Date.now() - t0} ms`); }
const B0 = res.base, land = B0.land, size = new Map(); for (let c = 0; c < land.length; c++) if (land[c]) size.set(B0.basin[c], (size.get(B0.basin[c]) || 0) + 1);
const major = Uint8Array.from(land, (v, c) => v && size.get(B0.basin[c]) >= 1000 ? 1 : 0); let nm = 0; for (const v of major) nm += v;
console.log(`major basins (>= 100,000 km2 in the base run) cover ${(100 * nm / B0.landCells).toFixed(0)}% of the land`);
for (const k of ["reseed", "reseed2", "noclimate", "type", "grounded", "strong"]) {
  // compare only cells that are land in both (the relief re-seed moves the coast a little)
  const both = Uint8Array.from(land, (v, c) => v && res[k].land[c] ? 1 : 0), bm = Uint8Array.from(major, (v, c) => v && res[k].land[c] ? 1 : 0);
  console.log(`base vs ${k}: ARI all land ${ari(B0.basin, res[k].basin, both).toFixed(3)}, in major basins ${ari(B0.basin, res[k].basin, bm).toFixed(3)}`); }
