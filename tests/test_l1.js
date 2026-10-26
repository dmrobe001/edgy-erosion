// L1 checks: the stream-power solver against its analytic steady state, then a real continent coupled to L0:
// crust budget closes, every land cell drains to the sea (or a recorded inland basin), runs are deterministic, peaks stay bounded.
// Hack's law and the L0-L1 divergence are printed as diagnostics. Then the suture lookahead: L0 restored exactly,
// the welding continent found, stitched on at the weld, rivers crossing the suture.
const { loadSim, loadL1, runTo, arg } = require("./lib");
const SIM = loadSim(), L1 = loadL1(SIM);
let ok = true; const check = (cond, msg) => { if (!cond) { ok = false; console.log("FAIL " + msg); } };

// 1. steady state on a planar grid: uniform uplift U, edges at sea level. The implicit scheme's fixed point is
//    exactly S * A^m = U / K along every flow path, so this is a tight test of routing and the solver.
{ const nx = 80, ny = 60, G = L1.planarGrid(nx, ny, 1), S = L1.scratch(G.nC), e = new Float64Array(G.nC), alive = new Uint8Array(G.nC).fill(1), out = new Uint8Array(G.nC);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) { const l = y * nx + x; out[l] = x == 0 || y == 0 || x == nx - 1 || y == ny - 1 ? 1 : 0; e[l] = out[l] ? 0 : 1e-3 * Math.sin(l * 12.9898) ** 2; }
  const U = 0.5, K = 2, m = 0.5, dt = 0.5;
  for (let it = 0; it < 400; it++) { for (let l = 0; l < G.nC; l++) if (!out[l]) e[l] += U * dt; L1.route(G, e, alive, out, S); e.set(L1.streamPower(e, S, K, m, dt)); }
  let dev = 0, n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let l = 0; l < G.nC; l++) { if (out[l]) continue; const s = (e[l] - e[S.rec[l]]) / S.rdist[l], A = S.A[l]; dev = Math.max(dev, Math.abs(s * Math.sqrt(A) / (U / K) - 1));
    const X = Math.log(A), Y = Math.log(s); n++; sx += X; sy += Y; sxx += X * X; sxy += X * Y; }
  const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  console.log(`steady state: slope-area exponent ${slope.toFixed(4)} (expect -0.5), max |ksn/(U/K) - 1| ${dev.toExponential(1)}`);
  check(Math.abs(slope + 0.5) < 1e-6 && dev < 1e-9, "stream-power steady state"); }

// 2. a continent coupled to L0
const n = +arg("n", 64), r = +arg("r", 4), W = +arg("win", 20), T0 = +arg("t0", 400);
function make() {
  const w = runTo(new SIM.World({ n, seed: "genesis", nPlates: 10, mode: "scatter" }), T0);
  let best = null; const seen = new Set();
  for (let i = 0; i < w.N; i++) { if (!w.wK[i] || w.elev[i] <= 0 || seen.has(w.owner[i] + ":" + w.bodyIdx[i])) continue; const s = L1.select(w, i, 3); if (!s) continue;
    for (const k of s.slots) seen.add(s.plate + ":" + k); if (!best || s.nCont > best.s.nCont) best = { i, s }; }
  return new L1.Region(w, { cell: best.i, r, window: W });
}
const R = make(); const t0 = Date.now(); while (R.t < W) R.advance(); const ms = (Date.now() - t0) / R.steps;
const { G, S, alive, outlet, cont, e } = R, h = R.history[R.history.length - 1];
console.log(`continent: ${R.nCont} L0 cells, ${G.nC} L1 cells of ${(R.world.cellKm / r).toFixed(0)} km, ${R.steps} steps, ${ms.toFixed(0)} ms/step`);
console.log(`after ${R.t.toFixed(0)} Myr: mean land elevation ${h.meanElevKm.toFixed(2)} km (start ${R.history[0].meanElevKm.toFixed(2)}), max ${h.maxElevKm.toFixed(1)} km, erosion ${Math.round(h.erosionRateKmPerMyr * 1000)} m/Myr`);
const bal = R.balance(); console.log(`crust budget: ${JSON.stringify(Object.fromEntries(Object.entries(R.budget).map(([k, v]) => [k, Math.round(v / 1e3)])))} (thousand km^3), error ${bal.error.toExponential(1)}`);
check(Math.abs(bal.error) < 1e-9, "crust budget closes");
check(h.maxElevKm < 10.5, "peaks bounded by the thickness cap plus the strength limit");
// every live land cell's receiver chain ends at the sea or at a recorded inland basin, with no cycles
let stuck = 0, inland = 0, land = 0;
for (let l = 0; l < G.nC; l++) { if (!alive[l] || outlet[l]) continue; land++; let c = l, g = 0; while (S.rec[c] !== c && g < G.nC) { c = S.rec[c]; g++; }
  if (S.sink[c]) inland++; else if (!outlet[c]) stuck++; }
console.log(`drainage: ${land} land cells; ${inland} drain to ${S.sinks} inland basin(s), ${stuck} reach neither sea nor basin`); check(stuck === 0, "all land drains to the sea or an inland basin");
// Hack's law: longest upstream flow length against drainage area, for cells draining 20+ cells
{ const Lup = new Float64Array(G.nC); for (let q = S.ns - 1; q >= 0; q--) { const l = S.stack[q], rc = S.rec[l]; if (rc !== l) Lup[rc] = Math.max(Lup[rc], Lup[l] + S.rdist[l]); }
  let n2 = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; const amin = 20 * G.area[0];
  for (let l = 0; l < G.nC; l++) { if (!alive[l] || outlet[l] || S.A[l] < amin) continue; const X = Math.log(S.A[l]), Y = Math.log(Lup[l]); n2++; sx += X; sy += Y; sxx += X * X; sxy += X * Y; }
  console.log(`Hack's law: length ~ area^${((n2 * sxy - sx * sy) / (n2 * sxx - sx * sx)).toFixed(2)} over ${n2} cells (Earth's rivers: roughly 0.5-0.6)`); }
const dv = R.divergence(); console.log(`L0 vs L1: mean thickness difference ${dv.meanThicknessDiffKm.toFixed(2)} km, rms ${dv.rmsThicknessDiffKm.toFixed(2)} km over ${dv.parcels} parcels`);
// determinism: a second run from scratch gives the same surface
{ const R2 = make(); while (R2.t < W) R2.advance(); let same = R2.G.nC === G.nC; for (let l = 0; same && l < G.nC; l++) if (R2.e[l] !== e[l]) same = false;
  check(same, "deterministic"); console.log(`determinism: ${same ? "identical" : "different"} on a second run`); }
// 3. suture lookahead: a continent that another continent welds onto within the window
{ const hash = w => { let a = 0; for (let i = 0; i < w.N; i++) a = (a * 31 + Math.round(w.elev[i] * 1e4) + w.owner[i]) | 0; return a; };
  const w = runTo(new SIM.World({ n: 64, seed: "genesis", nPlates: 10, mode: "scatter", contTarget: 0.5 }), 1000);
  let cell = -1; for (let i = 0; i < w.N && cell < 0; i++) if (w.wK[i] && w.owner[i] === 21) { const s = L1.select(w, i, 3); if (s && s.nCont === 451) cell = i; }
  check(cell >= 0, "lookahead test continent found");
  const before = hash(w), la = new L1.Lookahead(w, cell, 100); while (!la.advance(1e9)); const partners = la.finish();
  check(hash(w) === before && w.t === 1000, "lookahead restores L0 exactly");
  console.log(`lookahead: partners ${JSON.stringify(partners)}`); check(partners.length > 0, "lookahead finds the welding continent");
  const R = new L1.Region(w, { cell, r: 4, window: 100, partners }); while (R.t < 100) R.advance();
  const { G, S, alive, outlet } = R, piece = l => R.pieceOf[(l / G.rr) | 0]; let cross = 0, crossA = 0;
  for (let l = 0; l < G.nC; l++) { if (!alive[l] || outlet[l]) continue; const r = S.rec[l]; if (r !== l && piece(r) !== piece(l)) { cross++; crossA = Math.max(crossA, S.A[l]); } }
  const b2 = R.balance();
  console.log(`pieces: ${R.pieces.length}, welds ${JSON.stringify(R.welds)}; ${cross} cells drain across the suture (largest ${Math.round(crossA)} km^2); budget error ${b2.error.toExponential(1)}`);
  check(R.welds.length > 0 && R.welds[0].links > 0, "pieces stitched at the weld"); check(cross > 0, "rivers cross the suture"); check(Math.abs(b2.error) < 1e-9, "multi-piece crust budget closes"); }
console.log(ok ? "PASS L1" : "FAIL L1"); process.exitCode = ok ? 0 : 1;
