// Stage 2 exit checks (L0_MESH_SPEC.md §11): plate speeds from the torque balance, against Earth.
// Reported per seed, after a transient (--skip Myr; there is no spin-up, so slabs start at zero length):
//   speed: root mean square over the crust's area (cm/yr), mean and range; Earth's is about 4-5 cm/yr
//   net rotation of the crust in the mantle frame (deg/Myr); estimates for Earth 0.1-0.4
//   Forsyth & Uyeda (1975): rank correlation, across plates larger than 1% of the sphere, between a plate's speed and
//     the share of its boundary where it is the subducting side (sampled every 10 Myr; pooled, and median per sample)
//   collisions: for each plate pair converging at 1 cm/yr or more when at least --collkm (400) km of their shared
//     convergent boundary becomes continent against continent, the convergence rate then against 10 and 20 Myr later
//     (Earth: India slowed about 3x, 15 to 4-5 cm/yr)
//   tail: share of steps when any plate (larger than 1%, and any size) moves faster than 10 and 15 cm/yr, and the
//     longest such spell
// Options: --n=96 --t=600 --skip=100 --seeds=genesis,atlas,kepler --cgeom=(default PH value) --ph.<constant>=<value>
const { loadSim, loadMesh, arg } = require("../tests/lib");
const SIM = loadSim(), MESH = loadMesh(SIM);
const collKm = +arg("collkm", 400), n = +arg("n", 96), T = +arg("t", 600), skip = +arg("skip", 100), seeds = arg("seeds", "genesis,atlas,kepler").split(",");
if (arg("cgeom", "")) MESH.PH.cGeom = +arg("cgeom", "");
// any other constant of the balance, for sensitivity probes: --ph.cColl=1e23
for (const a of process.argv.slice(2)) { const m = /^--ph\.(\w+)=(.+)$/.exec(a); if (m) { if (!(m[1] in MESH.PH)) throw new Error("no constant " + m[1]); MESH.PH[m[1]] = +m[2]; console.log(`PH.${m[1]} = ${m[2]}`); } }
const spearman = (a, b) => { const rk = x => { const o = x.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]), r = new Array(x.length); o.forEach((e, k) => r[e[1]] = k); return r; };
  const ra = rk(a), rb = rk(b), m = (a.length - 1) / 2; let s = 0, sa = 0, sb = 0; for (let i = 0; i < a.length; i++) { s += (ra[i] - m) * (rb[i] - m); sa += (ra[i] - m) ** 2; sb += (rb[i] - m) ** 2; } return s / Math.sqrt(sa * sb); };
const med = x => { const s = [...x].sort((a, b) => a - b); return s.length ? s[Math.floor(s.length / 2)] : NaN; };
const all = { rms: [], net: [], rho: [], coll: [], collHalf: [] };
for (const seed of seeds) {
  const w = new MESH.MeshWorld({ n, seed, nPlates: 10, mode: "scatter", nnr: true, contTarget: 0.4 }), t0 = Date.now(), evN = {};
  { const log = w.log; w.log = function (type, ...a) { if (this.t >= skip) evN[type] = (evN[type] || 0) + 1; return log.call(this, type, ...a); }; }
  const rms = [], net = [], xs = [], ys = [], perSample = [], series = [], onset = new Map(), onsetHalf = new Map(), coll = [], collHalf = [];
  let nBig = 0, steps = 0, fast10 = 0, fast15 = 0, any10 = 0, any15 = 0, spell = 0, longest = 0, nextS = skip;
  while (w.t < T) {
    w.step(); w.updateStats(true); steps++;
    // pair convergence, for the collision check (from the start, so onsets are seen)
    // two onsets per pair: contact (collkm of continent against continent) and dominance (also half the pair's
    // convergent boundary); each needs the pair converging at 1 cm/yr or more then
    for (const [mapk, test] of [[onset, e => e[2] >= collKm], [onsetHalf, e => e[2] >= collKm && e[2] >= 0.5 * e[0]]])
      for (const [k, e] of w.fstat.pairs) { const o = mapk.get(k);
        if (!o && test(e) && e[1] >= 1 && w.t > 5) mapk.set(k, { t: w.t, v0: e[1], v: [] });
        else if (o && !o.done) { o.v.push([w.t - o.t, e[1]]); if (w.t - o.t >= 20) o.done = true; } }
    for (const mapk of [onset, onsetHalf]) for (const [k, o] of mapk) if (!o.done && !w.fstat.pairs.has(k)) { o.v.push([w.t - o.t, 0]); if (w.t - o.t >= 20) o.done = true; } // no longer converging
    if (w.t < skip) continue;
    rms.push(w.rmsSpeed); net.push(w.netRot);
    const big = w.plates.filter(P => P.alive && P.area > 0.01 * w.N), vmaxBig = Math.max(0, ...big.map(P => P.vrms)), vmaxAny = Math.max(0, ...w.plates.filter(P => P.alive && P.area > 0).map(P => P.vrms));
    nBig += big.length; if (vmaxBig > 10) fast10++; if (vmaxBig > 15) fast15++; if (vmaxAny > 10) any10++; if (vmaxAny > 15) any15++;
    spell = vmaxBig > 10 ? spell + w.dt : 0; longest = Math.max(longest, spell);
    if (w.t >= nextS) { nextS += 10; const a = [], b = [];
      for (const P of big) { a.push(P.vrms); b.push(P.perim > 0 ? P.subLen / P.perim : 0); xs.push(P.vrms); ys.push(P.perim > 0 ? P.subLen / P.perim : 0); }
      if (a.length >= 4) perSample.push(spearman(a, b)); series.push([Math.round(w.t), w.rmsSpeed, big.length]); }
  }
  // contacts during the start's transient are not counted
  for (const [mapk, out] of [[onset, coll], [onsetHalf, collHalf]])
    for (const [k, o] of mapk) { if (o.t < skip) continue; const at = d => { let best = null; for (const [dt, v] of o.v) if (dt <= d + 1e-9) best = v; return best; };
      out.push({ pair: `${Math.floor(k / 4096)}-${k % 4096}`, t: Math.round(o.t), v0: o.v0, v10: at(10), v20: at(20) }); }
  const mean = rms.reduce((a, b) => a + b, 0) / rms.length, rho = spearman(xs, ys);
  console.log(`\nseed ${seed}, n=${n}, ${skip}-${T} Myr (${((Date.now() - t0) / steps).toFixed(0)} ms/step)`);
  console.log(`  rms speed  mean ${mean.toFixed(2)} cm/yr, range ${Math.min(...rms).toFixed(2)}-${Math.max(...rms).toFixed(2)}; every 10 Myr: ${series.filter((_, i) => i % 5 === 0).map(s => `${s[0]}:${s[1].toFixed(1)}`).join(" ")}`);
  console.log(`  plates larger than 1% of the sphere: ${(nBig / rms.length).toFixed(1)} on average (Earth: about 10)`);
  const per = x => ((x || 0) * 1000 / (T - skip)).toFixed(0);
  console.log(`  events per Gyr: subduction starts ${per(evN.initiation)}, welds ${per(evN.suture)}, timed rifts ${per(evN.rift)}, plates subducted away ${per(evN.lost)}`);
  console.log(`  net rotation  median ${med(net).toFixed(3)} deg/Myr, max ${Math.max(...net).toFixed(3)}`);
  console.log(`  Forsyth-Uyeda  pooled rho ${rho.toFixed(2)} (${xs.length} plate samples); per-sample median ${med(perSample).toFixed(2)} (${perSample.length} samples)`);
  console.log(`  tail  steps with a plate >1% faster than 10 cm/yr: ${(100 * fast10 / rms.length).toFixed(1)}%, 15: ${(100 * fast15 / rms.length).toFixed(1)}%; any plate: ${(100 * any10 / rms.length).toFixed(1)}%, ${(100 * any15 / rms.length).toFixed(1)}%; longest spell >10 (large plates) ${longest.toFixed(0)} Myr`);
  const f = x => x === null ? "  -  " : x.toFixed(1).padStart(5);
  console.log(`  collisions (convergence cm/yr at contact, +10, +20 Myr; ratio contact/+20):`);
  for (const c of coll) console.log(`    plates ${c.pair.padEnd(6)} t=${String(c.t).padStart(4)}  ${f(c.v0)} ${f(c.v10)} ${f(c.v20)}  ${c.v20 !== null && c.v20 > 0.05 ? (c.v0 / c.v20).toFixed(1) + "x" : c.v20 !== null ? "stopped" : "-"}`);
  all.rms.push(mean); all.net.push(med(net)); all.rho.push(rho); all.coll.push(...coll); all.collHalf.push(...collHalf);
}
console.log(`\nall seeds: rms ${all.rms.map(x => x.toFixed(2)).join(", ")} cm/yr; net rotation ${all.net.map(x => x.toFixed(3)).join(", ")} deg/Myr; Forsyth-Uyeda rho ${all.rho.map(x => x.toFixed(2)).join(", ")}`);
for (const [label, L] of [["at contact", all.coll], ["once dominant (half the pair's convergent boundary)", all.collHalf]]) {
  const ratios = L.filter(c => c.v20 !== null).map(c => c.v20 > 0.05 ? c.v0 / c.v20 : Infinity);
  console.log(`collisions ${label}, with 20 Myr of record: ${ratios.length}; slowed 2x or more: ${ratios.filter(r => r >= 2).length}; slowed 2-4x: ${ratios.filter(r => r >= 2 && r <= 4).length}; stopped or reversed: ${ratios.filter(r => r === Infinity).length}; median ratio ${med(ratios).toFixed(1)}`); }
