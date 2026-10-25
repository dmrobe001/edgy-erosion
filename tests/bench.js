// Milliseconds per simulation step at each resolution, plus handoff extraction cost.
const { loadSim, arg } = require("./lib");
const SIM = loadSim(), T = +arg("t", 200);
for (const n of [64, 96, 128]) {
  const w = new SIM.World({ n, seed: "genesis", nPlates: 10, mode: "scatter" });
  let s = 0; const t0 = Date.now(); while (w.t < T) { w.step(); s++; }
  const ms = (Date.now() - t0) / s, t1 = Date.now(); const h = w.extractFaults(); h.l1 = w.buildL1(); const ex = Date.now() - t1;
  const t2 = Date.now(); const json = JSON.stringify(w.exportHandoff(h)); const exp = Date.now() - t2;
  console.log(`n=${n} cells=${w.N} dt=${w.dt.toFixed(2)} Myr  ${ms.toFixed(1)} ms/step  extract ${ex} ms  export ${exp} ms, ${(json.length / 1e6).toFixed(1)} MB`);
}
