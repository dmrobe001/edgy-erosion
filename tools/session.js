// Sessions exported from the page ("Export session"): a JSON log from which the whole session reruns exactly.
//   node tools/session.js session.json            summary, and which commit built the page that saved it
//   node tools/session.js session.json --rerun    rerun L0 and L1 from the log, check against the recorded values,
//                                                 and stash the results next to the log (session.stash) for reuse
// As a module:
//   const S = require("./tools/session.js"), sess = S.load(file);
//   const r = S.reproduce(sess, SIM, L1)  -> {world, entry, plan, sheet, basins, handoff, checks}; world ends at the
//                                            session's present, the L1 objects as they were when it was refined
//   S.stash(r, file); S.unstash(file, SIM, L1, sess) -> the same, from disk, without rerunning
// L0 is deterministic given its settings and step count; L1 given the selected cell, the saved state it started from
// (L0 at entrySteps), the present it ran to, its window and its settings. Reruns need the same code: sess.build is a
// hash of the page's sources, and findBuild() finds the commit with that hash.
const fs = require("fs"), path = require("path"), v8 = require("v8"), crypto = require("crypto"), { execFileSync } = require("child_process");
const ROOT = path.join(__dirname, "..");

function load(file) { return JSON.parse(fs.readFileSync(file, "utf8")); }
function l0Hash(w) { let a = 0; for (let i = 0; i < w.N; i++) a = (a * 31 + Math.round(w.elev[i] * 1e4) + w.owner[i]) | 0; return a; }
function newWorld(sess, SIM) { const s = sess.settings; return new SIM.World({ n: s.resolution, seed: s.seed, nPlates: s.initialPlates, mode: s.start, nnr: true, contTarget: s.continentalTarget }); }

function reproduce(sess, SIM, L1, log = () => {}, L2 = null) {
  const w = newWorld(sess, SIM), R = sess.reproduce, L = R.L1, out = { world: w, L1, L2 };
  const to = steps => { while (w.steps < steps) w.step(); };
  if (L) {
    to(L.entrySteps); const entry = { t: w.t, steps: w.steps, snap: w.snapshot() }; out.entry = entry; log(`L0 at the saved state: ${w.t.toFixed(2)} Myr, ${w.steps} steps`);
    if (L.presentSteps != null && L.phase !== "replay") {
      to(L.presentSteps); const pick = L1.landmassAt(w, L.cell, L.bufferKm || 0); log(`L1 present: ${w.t.toFixed(2)} Myr; landmass of ${pick.land} land cells`);
      const rp = new L1.Replay(w, entry, L.windowMyr, pick); while (!rp.advance(1e9)); const plan = rp.finish();
      w.restore(entry.snap); w.trackParcels(); to(plan.t0Steps);
      const sh = new L1.Sheet(w, plan, L.sheet); while (w.steps < plan.endSteps) sh.advance(); log(`sheet: ${sh.steps} steps over ${sh.t.toFixed(1)} Myr`);
      Object.assign(out, { plan, sheet: sh });
      if (L.phase === "done") { out.basins = L1.basins(sh); out.handoff = L1.handoff(sh, out.basins); }
      w.untrackParcels();
    }
  }
  to(R.L0.steps); log(`L0 at the present: ${w.t.toFixed(2)} Myr, ${w.steps} steps`);
  out.checks = check(sess, out); return out;
}
// compare a rerun with the values the page recorded
function check(sess, r) {
  const c = sess.checks, res = [], eq = (name, a, b) => res.push({ name, recorded: a, rerun: b, ok: a === b || (typeof a === "number" && typeof b === "number" && Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(a))) });
  eq("L0 hash", c.L0.hash, l0Hash(r.world)); eq("plates alive", c.L0.platesAlive, r.world.plates.filter(p => p.alive).length);
  if (c.L1 && r.sheet) { const sh = r.sheet; eq("sheet steps", c.L1.sheetSteps, sh.steps); eq("markers", c.L1.markers, sh.nm); eq("max elevation", c.L1.last.maxElevKm, sh.history.at(-1).maxElevKm);
    if (c.L1.basins && r.basins) { eq("basins", c.L1.basins.nBasins, r.basins.nBasins); eq("inland share", c.L1.basins.inlandFrac, r.basins.inlandFrac); }
    if (c.L1.tiles != null && r.handoff) eq("tiles", c.L1.tiles, r.handoff.tiles.length); }
  const L2 = sess.reproduce.L2;
  if (c.L2 && L2 && r.basins) { const s = r.L1.basinSelection(r.basins, r.plan.proj, L2.xKm, L2.yKm, L2.radiusKm); eq("L2 tapped basin", c.L2.tapped, s.tapped); eq("L2 basins", c.L2.basins, s.basins.length); eq("L2 cells", c.L2.cells, s.cells);
    if (c.L2.reliefAtTap != null && r.L2) eq("L2 relief at the tap", c.L2.reliefAtTap, r.L2.relief(r.sheet, r.basins).z(L2.xKm, L2.yKm)); }
  return res;
}
// Stash: the rerun's states on disk (v8 serialization), so later work can skip the rerun
function stash(r, file) {
  const strip = o => { const d = {}; for (const k of Object.keys(o)) if (!/^(world|plan|proj|dynD)$/.test(k) && k[0] !== "_") d[k] = o[k]; return d; };
  const plan = r.plan && Object.assign({}, r.plan, { proj: { c: r.plan.proj.c, ea: r.plan.proj.ea, nr: r.plan.proj.nr } });
  fs.writeFileSync(file, v8.serialize({ world: r.world.snapshot(), entry: r.entry || null, plan, sheet: r.sheet ? strip(r.sheet) : null, basins: r.basins || null, checks: r.checks }));
}
function unstash(file, SIM, L1, sess) {
  const d = v8.deserialize(fs.readFileSync(file)), w = newWorld(sess, SIM); w.restore(d.world); const out = { world: w, entry: d.entry, checks: d.checks };
  if (d.plan) out.plan = Object.assign({}, d.plan, { proj: new L1.Projection(d.plan.proj.c, d.plan.proj.ea, d.plan.proj.nr) });
  if (d.sheet) { const sh = Object.create(L1.Sheet.prototype); Object.assign(sh, d.sheet, { world: w, plan: out.plan, proj: out.plan.proj }); out.sheet = sh; }
  if (d.basins) { out.basins = d.basins; out.handoff = L1.handoff(out.sheet, d.basins); }
  return out;
}
// the commit whose sources hash to the session's build id (as build.sh computes it)
// (src/l2.js joined the page later; a commit without it hashes the other four, as build.sh did then)
function srcHash(get) { const h = crypto.createHash("sha1"); for (const f of ["src/head.html", "src/sim.js", "src/l1.js", "src/l2.js", "src/ui.js"]) { let d; try { d = get(f); } catch (_) { continue; } h.update(d); } return h.digest("hex").slice(0, 12); }
function findBuild(id) {
  if (srcHash(f => fs.readFileSync(path.join(ROOT, f))) === id) return "the working tree";
  let commits = [];
  try { commits = execFileSync("git", ["-C", ROOT, "log", "--format=%h", "-n", "200"], { encoding: "utf8" }).trim().split("\n"); } catch (_) {}
  for (const c of commits) { try { if (srcHash(f => execFileSync("git", ["-C", ROOT, "show", `${c}:${f}`], { stdio: ["ignore", "pipe", "ignore"] })) === id) return "commit " + c; } catch (_) {} }
  return null;
}
module.exports = { load, reproduce, check, stash, unstash, findBuild, l0Hash };

if (require.main === module) {
  const file = process.argv[2]; if (!file) { console.log("usage: node tools/session.js session.json [--rerun]"); process.exit(1); }
  const s = load(file), R = s.reproduce, V = s.viewports;
  console.log(`${s.format} v${s.version}, saved ${s.savedAt}; build ${s.build} = ${findBuild(s.build) || "no commit in the last 200 matches"}`);
  console.log(`settings ${JSON.stringify(s.settings)}`);
  console.log(`L0: ${R.L0.steps} steps to ${R.L0.t.toFixed(2)} Myr; selected cell ${R.selectedCell}${V.L0.selected ? " at " + V.L0.selected.latLon.join(", ") : ""}`);
  if (R.L1) console.log(`L1 ${R.L1.phase}${R.L1.stale ? " (stale)" : ""}: cell ${R.L1.cell}, from L0 at ${R.L1.entrySteps} steps to ${R.L1.presentSteps}, window ${R.L1.windowMyr} Myr, sheet ${JSON.stringify(R.L1.sheet)}`);
  console.log(`viewports: L0 ${V.L0.colorMode} on the ${V.L0.projection}; L1 ${V.L1 ? `${V.L1.colorMode} at ${V.L1.zoom.toFixed(2)}x` : "-"}; L2 tile ${V.L2 ? V.L2.selectedTile : "-"}`);
  console.log(`${s.actions.length} actions: ${s.actions.map(a => a.a).join(" ")}`);
  if (process.argv.includes("--rerun")) {
    const { loadSim, loadL1, loadL2 } = require(path.join(ROOT, "tests", "lib.js")), SIM = loadSim(), L1 = loadL1(SIM), L2 = loadL2(SIM, L1);
    const r = reproduce(s, SIM, L1, m => console.log("  " + m), L2);
    for (const c of r.checks) console.log(`  ${c.ok ? "ok  " : "DIFF"} ${c.name}: recorded ${c.recorded}, rerun ${c.rerun}`);
    const sf = file.replace(/\.json$/, "") + ".stash"; stash(r, sf); console.log(`stashed to ${sf}`);
  }
}
