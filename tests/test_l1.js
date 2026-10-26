// L1 checks: the stream-power solver against its analytic steady state; the thin viscous sheet against its exact
// solutions; then a real landmass: the replay brings L0 back exactly, the sheet leaves L0 untouched, the crust budget
// closes, thickness and peaks stay bounded, runs are deterministic, and in the drainage pass every land cell drains
// to an outlet and the basin labels partition the land. Hack's law is printed.
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

// 2. the thin viscous sheet against exact solutions: a Newtonian sheet with no buoyancy, its left half dragged at +U
//    and its right half at -U, has u(x) = U sign(-x) (1 - exp(-|x| / lambda)). In a narrow strip (free sides, so the
//    across-strip stress vanishes) lambda = sqrt(3/4) ell; in the middle of a wide sheet (no sideways motion) lambda = ell.
function synthetic(nx, ny, h, ell, U0) {
  const s = Object.create(L1.Sheet.prototype);
  s.o = { h, n: 1, picard: 1, ellKm: ell, U0: 50, S0: 35, Ar: 3, cgTol: 1e-8, cgMax: 4000 }; s.nx = nx; s.ny = ny; s.x0 = -nx * h / 2; s.y0 = -ny * h / 2;
  const nc = nx * ny; s.e0 = s.o.U0 / ell; s.eta0 = 0.5 * 35 * 35 / (2 * 3 * s.e0); s.kappa = 4 * s.eta0 / (ell * ell);
  s.S = new Float64Array(nc).fill(35); s.Sref = new Float64Array(nc).fill(35); s.str = new Float64Array(nc).fill(1); s.act = new Uint8Array(nc).fill(1);
  s.eta = new Float64Array(nc).fill(s.eta0); s.eps = new Float64Array(nc); s.Ux = new Float64Array(nc); s.Uy = new Float64Array(nc);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) s.Ux[j * nx + i] = (s.x0 + (i + 0.5) * h) < 0 ? U0 : -U0;
  s.u = new Float64Array((nx + 1) * ny); s.v = new Float64Array(nx * (ny + 1)); return s;
}
for (const [label, ny, f, row] of [["narrow strip, free sides", 4, Math.sqrt(3 / 4), 2], ["middle of a wide sheet", 160, 1, 80]]) {
  const h = 20, ell = 300, U0 = 40, t = synthetic(300, ny, h, ell, U0); t.solve(); let e = 0;
  for (let i = 0; i <= t.nx; i++) { const x = t.x0 + i * h; if (Math.abs(x) > 2000) continue; const ex = U0 * Math.sign(-x) * (1 - Math.exp(-Math.abs(x) / (f * ell))); e = Math.max(e, Math.abs(t.u[row * (t.nx + 1) + i] - ex)); }
  console.log(`sheet, ${label}: max |u - exact| ${e.toFixed(3)} km/Myr (${(100 * e / U0).toFixed(2)}% of the plate speed)`);
  check(e / U0 < 0.005, "sheet: " + label);
}

// 2b. finite strain against exact solutions for uniform velocity gradients: pure shear gives F = diag(e^-et, e^et),
//     simple shear F = [[1, gt], [0, 1]], rigid rotation F = rotation by wt (the three branches of the 2x2 exponential)
for (const [label, uf, vf, exact] of [
  ["pure shear", (x, y) => -0.05 * x, (x, y) => 0.05 * y, t => [Math.exp(-0.05 * t), 0, 0, Math.exp(0.05 * t)]],
  ["simple shear", (x, y) => 0.05 * y, () => 0, t => [1, 0.05 * t, 0, 1]],
  ["rigid rotation", (x, y) => -0.05 * y, (x, y) => 0.05 * x, t => [Math.cos(0.05 * t), -Math.sin(0.05 * t), Math.sin(0.05 * t), Math.cos(0.05 * t)]]]) {
  const s = Object.create(L1.Sheet.prototype), nx = 40, ny = 40, h = 25; s.o = { h }; s.nx = nx; s.ny = ny; s.x0 = -500; s.y0 = -500; s.act = new Uint8Array(nx * ny).fill(1);
  s.u = new Float64Array((nx + 1) * ny); s.v = new Float64Array(nx * (ny + 1)); s.budget = { left: 0 };
  for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) s.u[j * (nx + 1) + i] = uf(s.x0 + i * h, s.y0 + (j + .5) * h);
  for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) s.v[j * nx + i] = vf(s.x0 + (i + .5) * h, s.y0 + j * h);
  s.nm = 3; s.mx = Float64Array.from([0, 60, -40]); s.my = Float64Array.from([0, -30, 50]); s.mS = new Float64Array(3).fill(35); s.mA = new Float64Array(3).fill(1);
  s.mIn = new Uint8Array(3).fill(1); s.mF = Float64Array.from([1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1]);
  for (let q = 0; q < 10; q++) s.advect(1);
  const ex = exact(10); let e = 0, ea = 0; for (let k = 0; k < 3; k++) { for (let q = 0; q < 4; q++) e = Math.max(e, Math.abs(s.mF[4 * k + q] - ex[q]));
    const det = s.mF[4 * k] * s.mF[4 * k + 3] - s.mF[4 * k + 1] * s.mF[4 * k + 2]; ea = Math.max(ea, Math.abs(det - s.mA[k])); }
  console.log(`finite strain, ${label}: max |F - exact| ${e.toExponential(1)}, max |det F - area ratio| ${ea.toExponential(1)}`);
  check(e < 1e-12 && ea < 1e-12, "finite strain: " + label);
}

// 2c. inland drainage on a round landmass 1 km high with a broad bowl 600 m deep in its middle:
//     - in the wet tropics its lake overflows to the sea;
//     - in the dry subtropics, standing still, sediment fills the bowl within the pass (no lake is left);
//     - a dry plain, 0.5 km high and sloping to its coast, has rivers that run dry before they reach the sea;
//     - in the dry subtropics, subsiding 0.1 km/Myr at its centre, it stays a lake with no outlet;
//     - a sea cut off from the ocean (floor 0.3 km below sea level) overflows in the wet tropics and has no outlet
//       in the dry subtropics.
function bowl(latDeg, floor, sub, plain) {
  const la = latDeg * Math.PI / 180, c = [Math.cos(la), 0, Math.sin(la)], ea = [0, 1, 0], nr = [-Math.sin(la), 0, Math.cos(la)];
  const nx = 60, ny = 60, H = 50, sheet = { nx, ny, o: { h: H }, x0: -1500, y0: -1500, e: new Float64Array(nx * ny), uRate: new Float64Array(nx * ny), proj: new L1.Projection(c, ea, nr), Qref: [1, 0, 0, 0, 1, 0, 0, 0, 1], world: { opts: { seed: "bowl" } } };
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = sheet.x0 + (i + .5) * H, y = sheet.y0 + (j + .5) * H, r = Math.hypot(x, y), g = Math.exp(-((r / 350) ** 2));
    sheet.e[j * nx + i] = r > 1200 ? -2 : plain ? 0.5 * (1 - r / 1200) + 0.02 : 1 - (1 - floor) * g; sheet.uRate[j * nx + i] = -sub * g; }
  const bs = L1.basins(sheet), mid = ((bs.ny / 2) | 0) * bs.nx + ((bs.nx / 2) | 0); let c0 = mid; while (bs.rec[c0] !== c0) c0 = bs.rec[c0];
  return { bs, mid, end: bs.land[c0] ? bs.lake[c0] : 0, lakeAtMid: bs.lake[mid] };
}
for (const [label, lat, floor, sub, plain, ok] of [["bowl, wet tropics", 5, 0.4, 0, 0, r => r.end === 0],
  ["bowl, dry subtropics, standing still", 24, 0.4, 0, 0, r => r.lakeAtMid === 0 && r.bs.nClosed === 0],
  ["dry plain", 24, 0, 0, 1, r => r.bs.nDry > 10 && r.bs.inlandFrac > 0.3],
  ["bowl, dry subtropics, subsiding", 24, 0.4, 0.1, 0, r => r.end === 2], ["cut-off sea, wet tropics", 5, -0.3, 0, 0, r => r.end === 0],
  ["cut-off sea, dry subtropics", 24, -0.3, 0, 0, r => r.end === 2]]) {
  const res = bowl(lat, floor, sub, plain), { bs, mid, end } = res, names = ["reaches the sea", "", "ends in a lake with no outlet", "", "ends where a river runs dry"];
  console.log(`inland drainage, ${label}: rain ${bs.climate.P[mid].toFixed(2)} m/yr, evaporation ${bs.climate.E[mid].toFixed(2)} m/yr; the centre ${names[end]}; ${bs.nFilled} hollows filled, ${bs.nClosed} lakes with no outlet, ${bs.nDry} rivers run dry, ${bs.nSea} cells of cut-off sea; ${(100 * bs.inlandFrac).toFixed(0)}% drains inland (${(100 * bs.areicFrac).toFixed(0)}% has no flowing water at all)`);
  check(ok(res), "inland drainage: " + label);
}

// 2d. dynamic topography: one slab bin, 5 Myr old (sunk 100 km), recorded in a single cell, must give a depression of
//     width 300 + 0.8 x 100 = 380 km: at that distance the deflection (less the far field) is exp(-1/2) of its depth.
//     The field averages zero over the sphere (return flow lifts the surface away from slabs).
{ const n = 96, cube = SIM.getCube(n), N = cube.N, nB = SIM.C.SLAB_BINS, fake = { n, N, cube, t: 100, slab: new Float32Array(nB * N), slabT: new Float64Array(nB).fill(-1e9) };
  const c0 = cube.lookup(0.3, 0.2, Math.sqrt(1 - 0.13)); fake.slabT[0] = 90; fake.slab[c0] = 1e6;
  const D = L1.dynTopo(fake), out = [0, 0, 0], p = [cube.pos[3 * c0], cube.pos[3 * c0 + 1], cube.pos[3 * c0 + 2]], at = (ang, dir) => {
    const e1 = [-p[1], p[0], 0], l = Math.hypot(...e1); e1[0] /= l; e1[1] /= l; const e2 = [p[1] * e1[2] - p[2] * e1[1], p[2] * e1[0] - p[0] * e1[2], p[0] * e1[1] - p[1] * e1[0]];
    const u = [Math.cos(dir) * e1[0] + Math.sin(dir) * e2[0], Math.cos(dir) * e1[1] + Math.sin(dir) * e2[1], Math.cos(dir) * e1[2] + Math.sin(dir) * e2[2]];
    return L1.dynAt(D, p[0] * Math.cos(ang) + u[0] * Math.sin(ang), p[1] * Math.cos(ang) + u[1] * Math.sin(ang), p[2] * Math.cos(ang) + u[2] * Math.sin(ang), out); };
  const far = at(Math.PI / 2, 0), peak = at(0, 0) - far; let r = 0; for (let k = 0; k < 8; k++) r += (at(380 / SIM.R_KM, k * Math.PI / 4) - far) / peak / 8;
  let sm = 0, sa = 0; for (let c = 0; c < D.g.cube2.N; c++) { sm += D.F[c] * D.g.area2[c]; sa += D.g.area2[c]; }
  console.log(`dynamic topography, one slab: depth ${(peak * 1000).toFixed(0)} m; at 380 km ${r.toFixed(3)} of it (Gaussian: ${Math.exp(-0.5).toFixed(3)}); mean over the sphere ${(sm / sa * 1000).toExponential(1)} m`);
  check(Math.abs(r - Math.exp(-0.5)) < 0.06 && peak < 0 && Math.abs(sm / sa) < 1e-9, "dynamic topography kernel"); }

// 3. a real landmass: snapshots every 10 Myr, pick at the present, replay, sheet, drainage pass
const n = +arg("n", 64), W = +arg("win", 10), T0 = +arg("t0", 400), H = +arg("h", 70), target = +arg("land", 200);
const w = new SIM.World({ n, seed: "genesis", nPlates: 10, mode: "scatter", contTarget: 0.5 }), snaps = [];
while (w.t < T0) { const last = snaps[snaps.length - 1]; if (!last || w.t - last.t >= 10 - 1e-9) { snaps.push({ t: w.t, steps: w.steps, snap: w.snapshot() }); while (snaps.length > 4) snaps.shift(); } w.step(); }
const hash = () => { let a = 0; for (let i = 0; i < w.N; i++) a = (a * 31 + Math.round(w.elev[i] * 1e4) + w.owner[i]) | 0; return a; }, before = hash();
let cell = -1, bd = Infinity; { const seen = new Uint8Array(w.N);
  for (let i = 0; i < w.N; i++) { if (seen[i] || w.elev[i] <= 0) continue; const lm = L1.landmassAt(w, i); for (const c of lm.cells) seen[c] = 1; const d = Math.abs(lm.land - target); if (d < bd) { bd = d; cell = i; } } }
// the selection rules: land first and all of it land; a wider buffer only adds; every island joined lies within the
// buffer of land already in; no water piece but the open ocean is left outside
{ const N = w.N, nb = w.cube.nb, pos = w.cube.pos; let prev = null, okSub = true, okLand = true, okChain = true, okSeas = true, line = [];
  for (const b of [0, 200, 500]) { const p = L1.landmassAt(w, cell, b), inP = new Uint8Array(N); for (const c of p.cells) inP[c] = 1;
    for (let h = 0; h < p.cells.length; h++) if ((w.elev[p.cells[h]] > 0) !== (h < p.land)) okLand = false;
    if (prev) for (const c of prev.cells) if (!inP[c]) okSub = false;
    // each land cell is within buffer + a cell of some other land cell of the pick, or touches one
    const lc = Array.from(p.cells.slice(0, p.land)), lim = Math.cos((b + 1.5 * w.cellKm) / 6371);
    if (p.land > 1) for (const a of lc) { let near = false; for (const c of lc) if (c !== a && pos[3 * a] * pos[3 * c] + pos[3 * a + 1] * pos[3 * c + 1] + pos[3 * a + 2] * pos[3 * c + 2] >= lim) { near = true; break; } if (!near) okChain = false; }
    const piece = new Int32Array(N).fill(-1); let np = 0; for (let i = 0; i < N; i++) { if (inP[i] || piece[i] >= 0) continue; const q = [i]; piece[i] = np; for (let h = 0; h < q.length; h++) for (let k = 0; k < 4; k++) { const j = nb[4 * q[h] + k]; if (!inP[j] && piece[j] < 0) { piece[j] = np; q.push(j); } } np++; }
    // pieces left outside: the open ocean (the largest), and otherwise only continental crust the shelf rule left out
    const size = new Array(np).fill(0); for (let i = 0; i < N; i++) if (piece[i] >= 0) size[piece[i]]++; const open = size.indexOf(Math.max(...size));
    for (let i = 0; i < N; i++) if (piece[i] >= 0 && piece[i] !== open && !w.wK[i]) okSeas = false;
    line.push(`buffer ${b} km: ${p.land} land cells (${p.islands} islands), ${p.cells.length} in all, ${p.seaCells} enclosed sea, ${np} pieces outside`); prev = p; }
  console.log("selection: " + line.join("; "));
  check(okLand && okSub && okChain && okSeas, "selection: land first, wider buffers only add, islands chained within the buffer, no sea left enclosed"); }
const pick = L1.landmassAt(w, cell), entry = [...snaps].reverse().find(s => s.t <= w.t - W - 1e-9);
function run() {
  const rp = new L1.Replay(w, entry, W, pick); while (!rp.advance(1e9)); const plan = rp.finish(), back = hash() === before;
  w.restore(entry.snap); w.trackParcels(); while (w.steps < plan.t0Steps) w.step();
  const t0 = Date.now(), sh = new L1.Sheet(w, plan, { h: H }); let maxS = 0;
  while (w.steps < plan.endSteps) { sh.advance(); for (let k = 0; k < sh.nm; k++) if (sh.mIn[k] && sh.mS[k] > maxS) maxS = sh.mS[k]; }
  w.untrackParcels(); return { plan, sh, back, maxS, ms: (Date.now() - t0) / sh.steps, after: hash() === before };
}
const A = run(), { sh, plan } = A, hl = sh.history[sh.history.length - 1], bal = sh.balance();
console.log(`landmass: ${pick.land} L0 land cells, ${plan.nObj} parcels; sheet ${sh.nx}x${sh.ny} at ${H} km, ${sh.nm} markers, ${sh.steps} steps, ${A.ms.toFixed(0)} ms/step`);
console.log(`after ${sh.t.toFixed(1)} Myr: mean land elevation ${hl.meanElevKm.toFixed(2)} km (start ${sh.history[0].meanElevKm.toFixed(2)}), max ${hl.maxElevKm.toFixed(1)} km, most thickening ${hl.maxThickeningKm.toFixed(1)} km, thickest marker ${A.maxS.toFixed(1)} km`);
console.log(`budget ${JSON.stringify(Object.fromEntries(Object.entries(sh.budget).map(([k, v]) => [k, Math.round(v / 1e3)])))} (thousand km^3), error ${bal.error.toExponential(1)}`);
check(A.back, "the replay brings L0 back to the present exactly");
check(A.after, "the sheet's pass leaves L0 exactly as it was");
check(Math.abs(bal.error) < 1e-9, "crust budget closes");
check(A.maxS < 100 && hl.maxElevKm < 11, "thickness and peaks bounded");
{ let e = 0; for (let k = 0; k < sh.nm; k++) { if (!sh.mIn[k]) continue; const F = sh.mF, det = F[4 * k] * F[4 * k + 3] - F[4 * k + 1] * F[4 * k + 2]; e = Math.max(e, Math.abs(det / (sh.mA[k] / sh.mA0[k]) - 1)); }
  const sf = sh.strainField(); let smax = 0, rmax = 0; for (let c = 0; c < sh.nx * sh.ny; c++) if (!isNaN(sf.lmin[c])) { smax = Math.max(smax, 1 - sf.lmin[c]); rmax = Math.max(rmax, Math.abs(sf.rot[c])); }
  console.log(`finite strain: most horizontal shortening ${(100 * smax).toFixed(0)}%, largest rotation ${(rmax * 180 / Math.PI).toFixed(0)} deg; det F against area ratio within ${e.toExponential(1)}`);
  check(e < 1e-9, "det F matches each marker's area change"); }
const B = run(); let same = B.sh.nm === sh.nm; for (let c = 0; same && c < sh.nx * sh.ny; c++) if (!(B.sh.e[c] === sh.e[c] || isNaN(B.sh.e[c]) && isNaN(sh.e[c]))) same = false;
check(same, "deterministic"); console.log(`determinism: ${same ? "identical" : "different"} on a second run`);
// the drainage pass
{ const t0 = Date.now(), bs = L1.basins(sh), { S } = bs, rec = bs.rec; let stuck = 0, label = 0, sum = 0;
  for (let l = 0; l < bs.land.length; l++) { if (!bs.land[l]) continue; let c = l, g = 0; while (rec[c] !== c && g < bs.land.length) { c = rec[c]; g++; } if (bs.land[c] && !S.sink[c] && bs.lake[c] < 2) stuck++; if (bs.basin[l] !== bs.basin[rec[l]] && bs.land[rec[l]]) label++; }
  for (const v of bs.sizes) sum += v;
  console.log(`drainage pass: ${bs.nx}x${bs.ny} at ${bs.h} km, ${bs.landCells} land cells, ${bs.nBasins} basins (largest ${bs.sizes.slice(0, 3).join(", ")} cells), ${Date.now() - t0} ms; ${stuck} cells drain nowhere`);
  check(stuck === 0, "all land drains to the sea or an inland sink");
  console.log(`climate: ${bs.nClosed} lakes cannot overflow, ${bs.nDry} rivers run dry, ${bs.nFilled} hollows filled; ${(100 * bs.inlandFrac).toFixed(1)}% of the land drains inland (${(100 * bs.areicFrac).toFixed(1)}% areic)`); check(label === 0 && sum === bs.landCells, "basins partition the land");
  const Lup = new Float64Array(bs.land.length); for (let q = S.ns - 1; q >= 0; q--) { const l = S.stack[q], rc = S.rec[l]; if (rc !== l) Lup[rc] = Math.max(Lup[rc], Lup[l] + S.rdist[l]); }
  let n2 = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; const amin = 20 * bs.h * bs.h;
  for (let l = 0; l < bs.land.length; l++) { if (!bs.land[l] || S.A[l] < amin || Lup[l] <= 0) continue; const X = Math.log(S.A[l]), Y = Math.log(Lup[l]); n2++; sx += X; sy += Y; sxx += X * X; sxy += X * Y; }
  console.log(`Hack's law: length ~ area^${((n2 * sxy - sx * sy) / (n2 * sxx - sx * sx)).toFixed(2)} over ${n2} cells (Earth's rivers: roughly 0.5-0.6)`);
  // the L1 -> L2 handoff: tiles partition the land into whole basins (no receiver crosses a tile edge), each with an
  // outlet; imprints come oldest first; response time grows up the network; tile fields fill their windows
  const HO = L1.handoff(sh, bs), nT = HO.tiles.length; let unassigned = 0, cross = 0, area = 0, upOk = true;
  for (let l = 0; l < bs.land.length; l++) { if (!bs.land[l]) continue; if (HO.tileOf[l] < 0) unassigned++; const r = bs.rec[l]; if (r !== l && bs.land[r] && HO.tileOf[r] !== HO.tileOf[l]) cross++; if (r !== l && bs.land[r] && !(HO.responseMyr[l] > HO.responseMyr[r])) upOk = false; }
  for (const T of HO.tiles) area += T.areaKm2;
  const sorted = HO.imprints.every((I, k) => k === 0 || HO.imprints[k - 1].startMa >= I.startMa), outlets = HO.tiles.every(T => T.outlets.length > 0);
  const td = HO.tile(0), sizes = Object.entries(td.fine).filter(([k, v]) => Array.isArray(v)).every(([k, v]) => v.length === td.fine.nx * td.fine.ny) && td.coarse.stacks.length === td.coarse.nx * td.coarse.ny;
  const same = JSON.stringify(L1.handoff(sh, bs).summary()) === JSON.stringify(HO.summary()), small = HO.tiles.filter(T => T.areaKm2 < 5000).length;
  const kinds = {}; for (const T of HO.tiles) kinds[T.kind] = (kinds[T.kind] || 0) + 1;
  console.log(`L2 handoff: ${nT} tiles (${Object.entries(kinds).map(([k, v]) => v + " " + k).join(", ")}; ${small} under 5,000 km^2), ${HO.imprints.length} imprints; ${unassigned} land cells without a tile, ${cross} receivers crossing a tile edge; summary ${(JSON.stringify(HO.summary()).length / 1e3).toFixed(0)} kB`);
  check(unassigned === 0 && cross === 0 && area === bs.landCells * bs.h * bs.h, "tiles partition the land into whole basins");
  check(outlets && sorted && upOk && sizes && same, "handoff: outlets, imprint order, response time, field sizes, determinism"); }
console.log(ok ? "PASS L1" : "FAIL L1"); process.exitCode = ok ? 0 : 1;
