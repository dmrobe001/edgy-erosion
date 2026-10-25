// Where continental crust volume goes: instruments a patched copy of the sim (patterns must match src/sim.js).
const { loadSimPatched } = require("./lib");
const V = globalThis.V = { collisionConsumed: 0, collisionKept: 0, erosion: 0, redepositedOnLand: 0 };
const SIM = loadSimPatched([
  ["const amt=Pl.T[kl],W0=C.COLLW;", "const amt=Pl.T[kl],W0=C.COLLW;V.collisionConsumed+=Pl.T[kl];V.collisionKept+=amt;"],
  ["if(E[i]>0&&drop>0){const loss=C.KE*E[i]*dt;", "if(E[i]>0&&drop>0){const loss=C.KE*E[i]*dt;V.erosion+=loss;"],
  ["if(wK[j]){dT[j]+=sh;", "if(wK[j]){V.redepositedOnLand+=sh;dT[j]+=sh;"],
]);
const w = new SIM.World({ n: 64, seed: "genesis", nPlates: 10, mode: "scatter" });
while (w.t < 200) w.step(); for (const k in V) V[k] = 0;
const vol = () => { let s = 0; for (let i = 0; i < w.N; i++) if (w.wK[i]) s += w.wT[i]; return s; };
const v0 = vol(); while (w.t < 600) w.step(); const v1 = vol();
const k = x => Math.round(x / 1000) + "k";
console.log(`continental volume change 200-600 Myr: ${k(v1 - v0)} cell-km (${((v1 / v0 - 1) * 100).toFixed(0)}%)`);
for (const key in V) console.log(`  ${key}: ${k(V[key])}`);
