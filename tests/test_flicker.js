// How often continental cells change displayed surface rock per step, split into surface processes vs resampling.
const { loadSim, runTo } = require("./lib");
const SIM = loadSim(); const w = runTo(new SIM.World({ n: 96, seed: "genesis", nPlates: 10, mode: "scatter" }), 400);
let tot = 0, flip = 0, same = 0; const N = w.N;
for (let s = 0; s < 50; s++) {
  const R0 = Uint8Array.from(w.wR), K0 = Uint8Array.from(w.wK), P0 = Int16Array.from(w.owner), B0 = Int32Array.from(w.bodyIdx); w.step();
  for (let i = 0; i < N; i++) { if (!K0[i] || !w.wK[i]) continue; tot++; if (w.wR[i] !== R0[i]) { flip++; if (w.owner[i] === P0[i] && w.bodyIdx[i] === B0[i]) same++; } }
}
let dis = 0, pairs = 0; const nb = w.cube.nb;
for (let i = 0; i < N; i++) { if (!w.wK[i]) continue; for (let k = 0; k < 4; k++) { const j = nb[4 * i + k]; if (!w.wK[j]) continue; pairs++; if (w.wR[j] !== w.wR[i]) dis++; } }
console.log(`changing surface rock per step: ${(100 * flip / tot).toFixed(1)}% of continental cells`);
console.log(`  same parcel (surface processes): ${(100 * same / flip).toFixed(0)}%, different parcel sampled (resampling): ${(100 - 100 * same / flip).toFixed(0)}%`);
console.log(`neighbouring continental cells with different surface rock: ${(100 * dis / pairs).toFixed(0)}%`);
