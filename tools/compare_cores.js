// Stage 1 exit check (L0_MESH_SPEC.md §11): the marker core against the lattice core on the same seeds.
// Both start from the same state (the marker core converts the lattice core's start). Reported every 100 Myr and at
// the end: land and continental fraction, peak elevation, plate count, mean ocean-floor age, rifts and sutures so far,
// displayed-rock flicker (share of continental cells whose displayed rock changes in a step), and cost per step.
// Options: --n=96 --t=800 --seeds=genesis,atlas,kepler
const { loadSim, loadMesh, arg } = require("../tests/lib");
const SIM = loadSim(), MESH = loadMesh(SIM);
const n = +arg("n", 96), T = +arg("t", 800), seeds = arg("seeds", "genesis,atlas,kepler").split(",");
const fmt = (x, d = 3) => x.toFixed(d);
function run(make, label, seed) {
  const w = make({ n, seed, nPlates: 10, mode: "scatter", nnr: true, contTarget: 0.4 }), rows = [];
  let ms = 0, steps = 0, flick = 0, flickN = 0, prevR = null, next = 100;
  while (w.t < T) {
    const t0 = performance.now(); w.step(); ms += performance.now() - t0; steps++;
    if (prevR) { let ch = 0, c = 0; for (let i = 0; i < w.N; i++) if (w.wK[i]) { c++; if (w.wR[i] !== prevR[i]) ch++; } flick += c ? ch / c : 0; flickN++; }
    prevR = Uint8Array.from(w.wR);
    if (w.t >= next - 1e-9 || w.t >= T) { next += 100; w.updateStats(true);
      let peak = -1e9; for (let i = 0; i < w.N; i++) peak = Math.max(peak, w.elev[i]);
      rows.push({ t: Math.round(w.t), land: w.landFrac, cont: w.contFrac, peak, plates: w.plates.filter(p => p.alive).length, oceanAge: w.oceanAge,
        rifts: w._rifts || 0, sutures: w._sutures || 0, flicker: flickN ? flick / flickN : 0, msPerStep: ms / steps }); flick = 0; flickN = 0; }
  }
  return rows;
}
// count rifts and sutures over the whole run (the event log keeps only the last 400); MeshWorld inherits this log
for (const C of [SIM.World]) { const log = C.prototype.log; C.prototype.log = function (type, ...a) { if (type === "rift") this._rifts = (this._rifts || 0) + 1; if (type === "suture") this._sutures = (this._sutures || 0) + 1; return log.call(this, type, ...a); }; }
for (const seed of seeds) {
  const L = run(o => new SIM.World(o), "lattice", seed), M = run(o => new MESH.MeshWorld(o), "mesh", seed);
  console.log(`\nseed ${seed}, n=${n}        lattice | mesh`);
  console.log("  t     land        continent    peak km     plates  ocean age   rifts   sutures  flicker      ms/step");
  for (let k = 0; k < L.length; k++) { const a = L[k], b = M[k] || {};
    console.log(`  ${String(a.t).padStart(4)}  ${fmt(a.land)}|${fmt(b.land)}  ${fmt(a.cont)}|${fmt(b.cont)}  ${fmt(a.peak, 1)}|${fmt(b.peak, 1)}  ${String(a.plates).padStart(2)}|${String(b.plates).padStart(2)}  ${fmt(a.oceanAge, 0).padStart(4)}|${fmt(b.oceanAge, 0).padStart(4)}  ${String(a.rifts).padStart(3)}|${String(b.rifts).padStart(3)}  ${String(a.sutures).padStart(3)}|${String(b.sutures).padStart(3)}  ${fmt(a.flicker)}|${fmt(b.flicker)}  ${fmt(a.msPerStep, 0)}|${fmt(b.msPerStep, 0)}`); }
}
