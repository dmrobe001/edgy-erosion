// Net lithospheric rotation with and without removal, and RMS plate speed.
const { loadSim, runTo, arg } = require("./lib");
const SIM = loadSim(), T = +arg("t", 800);
let ok = true;
for (const nnr of [false, true]) {
  const w = new SIM.World({ n: 64, seed: "genesis", nPlates: 10, mode: "super", nnr }); const row = [];
  for (const t of [0, 50, 200, T]) { runTo(w, t); row.push(`${t}: net ${w.netRot.toFixed(3)} deg/Myr, rms speed ${w.rmsSpeed.toFixed(1)} cm/yr`); }
  console.log((nnr ? "removed" : "free") + "\n  " + row.join("\n  "));
  if (nnr && w.netRot > 1e-6) ok = false;
}
const w = new SIM.World({ n: 64, seed: "x", nPlates: 8, mode: "scatter", nnr: false });
const rel = () => w.plates[0].W.map((v, k) => v - w.plates[1].W[k]); const a = rel(); w.removeNetRotation(); const b = rel();
const drel = Math.hypot(...a.map((v, k) => v - b[k]));
console.log(`relative motion change from removal: ${drel.toExponential(1)}; residual net rotation: ${Math.hypot(...w.netRotation()).toExponential(1)}`);
if (drel > 1e-12) ok = false;
console.log(ok ? "PASS" : "FAIL"); process.exitCode = ok ? 0 : 1;
