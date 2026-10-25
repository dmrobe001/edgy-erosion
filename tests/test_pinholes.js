// Single-cell ocean "pinholes" (ocean cells with 3+ continental neighbours): inside freshly absorbed plates vs elsewhere.
const { loadSim } = require("./lib");
const SIM = loadSim();
for (const seed of ["genesis", "basalt-7"]) {
  const w = new SIM.World({ n: 96, seed, nPlates: 10, mode: "scatter" }), nb = w.cube.nb;
  const isPin = i => { if (w.wK[i]) return false; let k = 0; for (let q = 0; q < 4; q++) if (w.wK[nb[4 * i + q]]) k++; return k >= 3; };
  const orig = w.merge.bind(w); let region = null;
  w.merge = (a, b) => { const L = a.area >= b.area ? b : a; region = []; for (let i = 0; i < w.N; i++) if (w.owner[i] === L.id) region.push(i); orig(a, b); };
  const rin = [], rout = [];
  while (w.t < 800) { region = null; w.step();
    if (region && region.length > 50) { const set = new Set(region); for (let k = 0; k < 20; k++) w.step();
      let pi = 0, ci = 0, po = 0, co = 0; for (let i = 0; i < w.N; i++) { const c = w.wK[i], p = isPin(i) ? 1 : 0; if (set.has(i)) { pi += p; ci += c; } else { po += p; co += c; } }
      rin.push(1000 * pi / Math.max(1, ci)); rout.push(1000 * po / Math.max(1, co)); } }
  const avg = a => (a.reduce((x, y) => x + y, 0) / a.length).toFixed(1);
  console.log(`${seed}: pinholes per 1000 continental cells, 20 steps after a suture: absorbed plate ${avg(rin)}, elsewhere ${avg(rout)} (${rin.length} sutures)`);
}
