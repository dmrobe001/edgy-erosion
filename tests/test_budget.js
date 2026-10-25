// Long-run crust budget: land above sea, continental crust area, peak elevation, plate count, events.
const { loadSim, runTo, arg } = require("./lib");
const SIM = loadSim(), T = +arg("t", 1600), n = +arg("n", 64), target = +arg("cont", 0.5);
let ok = true;
for (const [mode, seed] of [["scatter", "genesis"], ["super", "genesis"], ["scatter", "basalt-7"]]) {
  const w = new SIM.World({ n, seed, nPlates: 10, mode, contTarget: target }); const row = [];
  for (const t of [0, 200, 400, 800, T]) { runTo(w, t); let mx = -9; for (const e of w.elev) if (e > mx) mx = e;
    row.push(`${t}: land ${(w.landFrac * 100).toFixed(0)}% cont ${(w.contFrac * 100).toFixed(0)}% peak ${mx.toFixed(1)} km`); }
  const ev = {}; for (const e of w.events) ev[e.type] = (ev[e.type] || 0) + 1;
  console.log(`${mode}/${seed}\n  ` + row.join("\n  ") + "\n  events " + JSON.stringify(ev));
  if (w.landFrac < 0.1 || Math.abs(w.contFrac - target) > 0.08) ok = false;
}
console.log(ok ? "PASS budget within bounds" : "FAIL budget drifted (land < 10% or continental crust > 8 points off target)");
process.exitCode = ok ? 0 : 1;
