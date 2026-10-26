// L1 checks: the stream-power solver against its analytic steady state; then a real continent coupled to L0, with
// the lookahead choosing which crust is simulated: L0 restored exactly after looking ahead, crust budget closes,
// every land cell drains (to the sea, a sea edge or a recorded inland basin), no overlapping crust is left, peaks
// stay bounded, runs are deterministic; then a continent that another collides with, checking that the other is
// brought in and that crust and rivers interact across the contact. Hack's law and L0-L1 divergence are printed.
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

// 2. a continent coupled to L0, with the lookahead deciding which crust is simulated
const n = +arg("n", 64), r = +arg("r", 4), W = +arg("win", 20), T0 = +arg("t0", 400);
function make(seed, ct, t0, plate, win = W) {
  const w = runTo(new SIM.World({ n, seed, nPlates: 10, mode: "scatter", contTarget: ct }), t0);
  let cell = -1, bn = 0; for (let i = 0; i < w.N; i++) { if (!w.wK[i] || w.elev[i] <= 0 || (plate >= 0 && w.owner[i] !== plate)) continue; const s = L1.select(w, i, 3); if (s && s.nCont > bn) { bn = s.nCont; cell = i; } }
  const hash = () => { let a = 0; for (let i = 0; i < w.N; i++) a = (a * 31 + Math.round(w.elev[i] * 1e4) + w.owner[i]) | 0; return a; }, before = hash();
  const la = new L1.Lookahead(w, cell, win); while (!la.advance(1e9)); const rel = la.finish();
  return { R: new L1.Region(w, { cell, r, window: win, relevance: rel }), rel, restored: hash() === before };
}
function checks(R, label) {
  const { G, S, live, outlet } = R, h = R.history[R.history.length - 1], rr = G.rr, part = l => R.pieceOf[(l / rr) | 0] * 65536 + R.slotPlate[(l / rr) | 0];
  const bal = R.balance(); console.log(`${label}: budget ${JSON.stringify(Object.fromEntries(Object.entries(R.budget).map(([k, v]) => [k, Math.round(v / 1e3)])))} (thousand km^3), error ${bal.error.toExponential(1)}`);
  check(Math.abs(bal.error) < 1e-9, label + ": crust budget closes");
  check(h.maxElevKm < 10.5, label + ": peaks bounded by the thickness cap plus the strength limit");
  // every live land cell's receiver chain ends at the sea, a sea edge, or a recorded inland basin, with no cycles
  let stuck = 0, land = 0, cross = 0;
  for (let l = 0; l < G.nC; l++) { if (!live[l] || outlet[l]) continue; land++; let c = l, g = 0; while (S.rec[c] !== c && g < G.nC) { c = S.rec[c]; g++; }
    if (!outlet[c] && !S.sink[c] && !S.seaEdge[c]) stuck++; const rc = S.rec[l]; if (rc !== l && part(rc) !== part(l)) cross++; }
  console.log(`${label}: ${land} land cells, ${stuck} reach neither sea nor basin; ${cross} drain across a contact between rigid parts; ${h.contacts} contacts`);
  check(stuck === 0, label + ": all land drains");
  // no two live cells of different rigid parts left overlapping
  R.place(); const idx = new L1.CellIndex(G.NR, G.nC); idx.clear(); let bad = 0;
  for (let l = 0; l < G.nC; l++) if (live[l] && R._pl[l]) idx.insert(l, R.W3[3 * l], R.W3[3 * l + 1], R.W3[3 * l + 2]);
  for (let a = 0; a < G.nC; a++) { if (!live[a] || !R._pl[a]) continue; const hh = 0.5 * Math.sqrt(G.area[a]) / SIM.R_KM, na = idx.neighbours(idx.g[a]);
    for (let q = 0; q < na; q++) for (let b = idx.first(idx.around[q]); b >= 0; b = idx.next[b]) { if (b <= a || part(a) === part(b)) continue; const W3 = R.W3, dx = W3[3 * a] - W3[3 * b], dy = W3[3 * a + 1] - W3[3 * b + 1], dz = W3[3 * a + 2] - W3[3 * b + 2]; if (dx * dx + dy * dy + dz * dz < hh * hh) bad++; } }
  check(bad === 0, label + ": no overlapping crust left"); return cross;
}
{ const { R, rel, restored } = make("genesis", 0.4, T0, -1); check(restored, "lookahead restores L0 exactly");
  const t0 = Date.now(); while (R.t < W - 1e-9) R.advance(); const ms = (Date.now() - t0) / R.steps, h = R.history[R.history.length - 1];
  console.log(`continent: ${rel.target.reduce((a, v) => a + v, 0)} target parcels, ${rel.pids.length} relevant in ${R.pieces.length} piece(s), ${G0(R)} L1 cells of ${(R.world.cellKm / r).toFixed(0)} km, ${R.steps} steps, ${ms.toFixed(0)} ms/step`);
  console.log(`after ${R.t.toFixed(0)} Myr: mean land elevation ${h.meanElevKm.toFixed(2)} km (start ${R.history[0].meanElevKm.toFixed(2)}), max ${h.maxElevKm.toFixed(1)} km, erosion ${Math.round(h.erosionRateKmPerMyr * 1000)} m/Myr`);
  checks(R, "continent");
  // Hack's law: longest upstream flow length against drainage area, for cells draining 20+ cells
  { const { G, S, live, outlet } = R, Lup = new Float64Array(G.nC); for (let q = S.ns - 1; q >= 0; q--) { const l = S.stack[q], rc = S.rec[l]; if (rc !== l) Lup[rc] = Math.max(Lup[rc], Lup[l] + S.rdist[l]); }
    let n2 = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; const amin = 20 * G.area[0];
    for (let l = 0; l < G.nC; l++) { if (!live[l] || outlet[l] || S.A[l] < amin || Lup[l] <= 0) continue; const X = Math.log(S.A[l]), Y = Math.log(Lup[l]); n2++; sx += X; sy += Y; sxx += X * X; sxy += X * Y; }
    console.log(`Hack's law: length ~ area^${((n2 * sxy - sx * sy) / (n2 * sxx - sx * sx)).toFixed(2)} over ${n2} cells (Earth's rivers: roughly 0.5-0.6)`); }
  const dv = R.divergence(); console.log(`L0 vs L1: mean thickness difference ${dv.meanThicknessDiffKm.toFixed(2)} km, rms ${dv.rmsThicknessDiffKm.toFixed(2)} km over ${dv.parcels} parcels`);
  // determinism: a second run from scratch gives the same surface
  const R2 = make("genesis", 0.4, T0, -1).R; while (R2.t < W - 1e-9) R2.advance(); let same = R2.G.nC === R.G.nC; for (let l = 0; same && l < R.G.nC; l++) if (R2.e[l] !== R.e[l]) same = false;
  check(same, "deterministic"); console.log(`determinism: ${same ? "identical" : "different"} on a second run`); }
function G0(R) { return R.G.nC; }

// 3. a continent another one collides with: the lookahead brings the second one in, and they interact
{ const CW = 60, { R, rel } = make("genesis", 0.5, 1000, 21, CW); const plates = new Set([...rel.pids].map(p => (p / rel.N) | 0));
  console.log(`collision case: relevant crust from plates ${[...plates].join(", ")}; pieces ${R.pieces.map(p => p.plate0 + ":" + p.nS).join(" ")}`);
  check(plates.size > 1, "other land comes within reach");
  let crossAny = 0, contactsMax = 0; while (R.t < CW - 1e-9) { R.advance(); contactsMax = Math.max(contactsMax, R.X.count); }
  crossAny = checks(R, "collision case");
  check(contactsMax > 0 && crossAny > 0, "crust in contact across rigid parts, with rivers crossing"); }
console.log(ok ? "PASS L1" : "FAIL L1"); process.exitCode = ok ? 0 : 1;
