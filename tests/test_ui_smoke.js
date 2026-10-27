// Runs the page's full script against a stubbed DOM: every view, both projections, inspector, extraction, export,
// and an L1 landmass: select, preview, refine (replay, a short sheet run, the drainage pass), every L1 view, zoom and
// inspect; then a session export, loaded back with tools/session.js and checked against the live state.
const fs = require("fs"), path = require("path"), os = require("os"), SESSION = require("../tools/session.js"), { loadSim, loadL1, loadL2 } = require("./lib");
const html = fs.readFileSync(path.join(__dirname, "..", "plates.html"), "utf8");
const code = html.slice(html.indexOf("<script>") + 8, html.lastIndexOf("</script>"));
const defaults = { grid: "64", seed: "genesis", nPlates: "10", mode: "scatter", speed: "3", contT: "50", l1Win: "10", l1H: "70", l1Ar: "3", l1Ell: "500" }, els = {};
function el(id) { return els[id] || (els[id] = { id, value: defaults[id] || "", checked: true, hidden: false, textContent: "", innerHTML: "", style: {}, className: "", dataset: {},
  clientWidth: 400, clientHeight: 400, width: 0, height: 0, children: [], options: [], disabled: false,
  addEventListener() {}, setAttribute() {}, setPointerCapture() {}, insertAdjacentHTML(_, h) { this.innerHTML += h; }, remove() {}, prepend(c) { this.children.unshift(c); },
  get lastChild() { return { remove() {} }; }, getBoundingClientRect() { return { left: 0, top: 0, width: 400, height: 400 }; },
  getContext() { return { createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) }), putImageData() {}, clearRect() {}, beginPath() {}, arc() {}, stroke() {}, fill() {},
    moveTo() {}, lineTo() {}, closePath() {}, setLineDash() {} }; } }); }
global.document = { getElementById: el, querySelectorAll: () => [], createElement: () => ({}), documentElement: {}, body: { classList: { add() {}, remove() {} } } };
global.setTimeout = f => { f(); return 0; }; global.clearTimeout = () => {};
global.getComputedStyle = () => ({ getPropertyValue: () => "#d2452f" });
global.window = { devicePixelRatio: 2, addEventListener() {} };
global.requestAnimationFrame = () => {};
const run = new Function("SESSION", "fs", "os", "path", "loadSim", "loadL1", "loadL2", code + `
 return (async () => {
  for (let k = 0; k < 150; k++) { world.step(); snapTick(); }
  for (const p of ["globe", "flat"]) { proj = p; pixIdx = null; for (const v of ["elev", "rock", "l1", "plates"]) { view = v; colorsDirty = true; render(); } }
  inspect({ clientX: 200, clientY: 120 }); document.getElementById("extract").onclick(); render();
  gZoom = 3; mZoom = 4; mC = [0.5, 0.2]; clampMap(); for (const p of ["globe", "flat"]) { proj = p; pixIdx = null; render(); } gZoom = 1; mZoom = 1; mC = [0, 0];
  const json = JSON.stringify(world.exportHandoff(handoff));
  if (!document.getElementById("summary").innerHTML.includes("formations")) throw new Error("summary not rendered");
  // a modest landmass (closest to 150 L0 cells of land), to keep the run short
  let land = -1, bd = Infinity; const seen = new Uint8Array(world.N);
  for (let i = 0; i < world.N; i++) { if (seen[i] || world.elev[i] <= 0) continue; const lm = L1.landmassAt(world, i); for (const c of lm.cells) seen[c] = 1; const d = Math.abs(lm.land - 150); if (d < bd) { bd = d; land = i; } }
  selCell = land; selMask = null; showColumn(); l1Select(land); render();
  if (!l1St || l1St.phase !== "picked" || !l1Pix) throw new Error("selection did not preview the landmass");
  const tNow = world.t; document.getElementById("l1Run").onclick();
  if (!l1Busy()) throw new Error("L1 did not start: " + document.getElementById("l1Info").innerHTML);
  while (l1Busy()) l1Tick();
  if (!l1St || l1St.phase !== "done") throw new Error("L1 did not finish: " + document.getElementById("l1Info").innerHTML);
  if (world.t !== tNow) throw new Error("L0 did not come back to the present");
  for (const v of ["terrain", "thick", "finite", "strain", "rock", "dyn", "rain", "basins", "tiles"]) { l1View = v; l1Legend(); l1Pix = null; l1Render(); }
  l1Zoom = 4; l1Pix = null; l1Render(); l1Inspect({ clientX: 200, clientY: 200 });
  // the L2 handoff: tap a land pixel in the tiles view, then build both downloads' contents
  l1View = "tiles"; l1Zoom = 1; l1Pix = null; l1Render(); { let q = -1; for (let k = 0; k < l1Pix.fine.length; k++) if (l1Pix.fine[k] >= 0 && l1St.bs.land[l1Pix.fine[k]]) { q = k; break; }
    l1Inspect({ clientX: (q % l1PW + 0.5) / l1PW * 400, clientY: (Math.floor(q / l1PW) + 0.5) / l1PH * 400 }); }
  if (l1TileSel < 0 || !document.getElementById("l1tile").innerHTML.includes("Tile " + l1TileSel)) throw new Error("tile not selected");
  const tj = JSON.stringify(l1St.ho.tile(l1TileSel)), sj = JSON.stringify(l1St.ho.summary());
  console.log("PASS ui smoke L2 handoff: " + l1St.ho.tiles.length + " tiles, summary " + (sj.length / 1e3).toFixed(0) + " kB, tile " + l1TileSel + " " + (tj.length / 1e3).toFixed(0) + " kB");
  // L2 selection: a tap in the basins view with a 200 km radius picks the tapped basin, the basins wholly inside, the enclosed
  l1View = "basins"; l1Pix = null; l1Render(); { const el = document.getElementById("l2Rad"); el.value = 200; el.oninput({ target: el });
    let q = -1; for (let k = Math.floor(l1Pix.fine.length / 2); k < l1Pix.fine.length; k++) if (l1Pix.fine[k] >= 0 && l1St.bs.land[l1Pix.fine[k]]) { q = k; break; }
    l1Inspect({ clientX: (q % l1PW + 0.5) / l1PW * 400, clientY: (Math.floor(q / l1PW) + 0.5) / l1PH * 400 });
    const s = l2Sel, bs = l1St.bs; if (!s || s.tapped !== bs.basin[l1Pix.fine[q]]) throw new Error("L2 tap did not select the tapped basin");
    // every cell of a basin wholly within the radius is within it (to the cell's centre)
    const p = l1St.plan.proj.inv(s.x, s.y, [0, 0, 0]), b = [0, 0, 0];
    for (let c = 0; c < bs.land.length; c++) if (s.mask[c] === 2) { l1St.plan.proj.inv(bs.x0 + (c % bs.nx + .5) * bs.h, bs.y0 + (Math.floor(c / bs.nx) + .5) * bs.h, b);
      if (Math.acos(Math.min(1, p[0] * b[0] + p[1] * b[1] + p[2] * b[2])) * 6371 > 200 + 1e-6) throw new Error("basin outside the radius selected"); }
    if (!document.getElementById("l2Info").innerHTML.includes("Selected")) throw new Error("L2 selection not reported");
    // the relief preview, drawn in one go (frames run at once here)
    // with a clock that jumps 100 ms at every reading, every frame is over its time budget: the burn-in must still route
    // a tile per frame (it once listed the pending tiles, ran out of time and never routed one)
    const raf = global.requestAnimationFrame, now = performance.now, queue = []; let fake = now.call(performance), frames = 0;
    global.requestAnimationFrame = f => queue.push(f); performance.now = () => (fake += 100);
    try { document.getElementById("l2Prev").onclick(); while (queue.length) { if (++frames > 20000) throw new Error("L2 burn-in makes no progress"); queue.shift()(); } }
    finally { global.requestAnimationFrame = raf; performance.now = now; }
    if (!l2V || !l2V.Z || !document.getElementById("l2stats").innerHTML.includes("Drawn at")) throw new Error("L2 relief preview not drawn");
    if (!s.ref || !document.getElementById("l2Info").innerHTML.includes("Outline at 500 m")) throw new Error("L2 outline not burnt in");
    console.log("PASS ui smoke L2 outline: " + Math.round(s.ref.areaKm2) + " km² at 500 m (L1: " + Math.round(s.ref.areaL1Km2) + "), " + s.ref.nTiles + " tiles (" + s.ref.evaluated + " evaluated, " + s.ref.rounds + " rounds), " + frames + " frames under a clock running 100 ms per reading");
    let zl = Infinity, zh = -Infinity; for (const v of l2V.Z) { zl = Math.min(zl, v); zh = Math.max(zh, v); }
    console.log("PASS ui smoke L2 relief preview: " + l2V.W + "x" + l2V.H + " px at " + l2V.s.toFixed(2) + " km/px, " + zl.toFixed(2) + " to " + zh.toFixed(2) + " km");
    // the L2 pass, 1 Myr, the same way
    { const raf2 = global.requestAnimationFrame, now2 = performance.now, queue2 = []; let fake2 = now2.call(performance), fr2 = 0;
      global.requestAnimationFrame = f => queue2.push(f); performance.now = () => (fake2 += 100);
      try { const el = document.getElementById("l2Myr"); el.value = 1; document.getElementById("l2Run").onclick(); while (queue2.length) { if (++fr2 > 20000) throw new Error("L2 run makes no progress"); queue2.shift()(); } }
      finally { global.requestAnimationFrame = raf2; performance.now = now2; }
      if (!s.run || !document.getElementById("l2stats").innerHTML.includes("Hack")) throw new Error("L2 run not done");
      for (const v of ["terrain", "change", "area", "order"]) { l2View = v; l2Draw(); }
      console.log("PASS ui smoke L2 run: " + s.run.nAct + " points, " + s.run.steps + " steps, deepest incision " + Math.round(s.run.stats.maxIncisionM) + " m, Hack " + s.run.stats.hack.toFixed(2)); }
    console.log("PASS ui smoke L2 selection: " + s.basins.length + " basins (" + s.nWithin + " within 200 km, " + s.nEnclosed + " enclosed), " + (s.areaKm2 / 1e3).toFixed(0) + " thousand km², " + s.ms.toFixed(0) + " ms"); }
  if (!document.getElementById("l1budget").innerHTML.includes("Drainage basins")) throw new Error("L1 summary not rendered");
  const sh = l1St.sh, bal = sh.balance(); if (Math.abs(bal.error) > 1e-9) throw new Error("L1 crust budget does not close: " + bal.error);
  console.log("PASS ui smoke L1: " + sh.nx + "x" + sh.ny + " sheet, " + sh.nm + " markers, t=" + sh.t.toFixed(1) + " Myr, " + l1St.bs.nBasins + " basins, budget error " + bal.error.toExponential(1));
  console.log("PASS ui smoke: t=" + world.t.toFixed(0) + " Myr, export " + (json.length / 1e6).toFixed(1) + " MB");
  // session export: the log reruns in a fresh simulation to the same state, checked against the recorded values
  selCell = land; l1View = "tiles"; l1TileSel = 0; const file = path.join(os.tmpdir(), "edgy-session-test.json"); fs.writeFileSync(file, JSON.stringify(sessionData()));
  const sess = SESSION.load(file), S2 = loadSim(), L2 = loadL1(S2), r = SESSION.reproduce(sess, S2, L2, () => {}, loadL2(S2, L2)), bad = r.checks.filter(c => !c.ok);
  if (bad.length || r.checks.length < 16) throw new Error("rerun differs: " + JSON.stringify(bad.length ? bad : r.checks));
  if (!sess.viewports.L1.visible || sess.viewports.L2.selectedTile !== 0 || sess.viewports.L1.colorMode !== "tiles") throw new Error("viewports not recorded");
  const sf = file.replace(".json", ".stash"); SESSION.stash(r, sf); const u = SESSION.unstash(sf, S2, L2, sess), a = u.sheet.strainField(), b = l1St.sh.strainField();
  for (let c = 0; c < a.lmin.length; c++) if (!(a.lmin[c] === b.lmin[c] || isNaN(a.lmin[c]) && isNaN(b.lmin[c]))) throw new Error("stashed sheet differs");
  if (SESSION.l0Hash(u.world) !== SESSION.l0Hash(world) || u.handoff.tiles.length !== l1St.ho.tiles.length) throw new Error("stash differs");
  console.log("PASS ui smoke session: log " + (fs.statSync(file).size / 1e3).toFixed(0) + " kB, " + sess.actions.length + " actions; rerun matches all " + r.checks.length + " checks; stash " + (fs.statSync(sf).size / 1e6).toFixed(0) + " MB restores the sheet exactly");
 })();
`);
run(SESSION, fs, os, path, loadSim, loadL1, loadL2).catch(e => { console.error(e); process.exitCode = 1; });
