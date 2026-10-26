// Runs the page's full script against a stubbed DOM: every view, both projections, inspector, extraction, export,
// and an L1 landmass: pick, replay, a short sheet run, the drainage pass, every L1 view, zoom and inspect.
const fs = require("fs"), path = require("path");
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
new Function(code + `
  for (let k = 0; k < 150; k++) { world.step(); snapTick(); }
  for (const p of ["globe", "flat"]) { proj = p; pixIdx = null; for (const v of ["elev", "rock", "faults", "l1", "plates", "age", "thick"]) { view = v; colorsDirty = true; render(); } }
  inspect({ clientX: 200, clientY: 120 }); document.getElementById("extract").onclick(); render();
  const json = JSON.stringify(world.exportHandoff(handoff));
  if (!document.getElementById("summary").innerHTML.includes("formations")) throw new Error("summary not rendered");
  // a modest landmass (closest to 150 L0 cells of land), to keep the run short
  let land = -1, bd = Infinity; const seen = new Uint8Array(world.N);
  for (let i = 0; i < world.N; i++) { if (seen[i] || world.elev[i] <= 0) continue; const lm = L1.landmassAt(world, i); for (const c of lm.cells) seen[c] = 1; const d = Math.abs(lm.land - 150); if (d < bd) { bd = d; land = i; } }
  document.getElementById("l1Pick").onclick(); l1PickAt(land); render();
  const tNow = world.t; document.getElementById("l1Run").onclick();
  if (!l1Busy()) throw new Error("L1 did not start: " + document.getElementById("l1Info").innerHTML);
  while (l1Busy()) l1Tick();
  if (!l1St || l1St.phase !== "done") throw new Error("L1 did not finish: " + document.getElementById("l1Info").innerHTML);
  if (world.t !== tNow) throw new Error("L0 did not come back to the present");
  for (const v of ["terrain", "thick", "finite", "strain", "rock", "basins"]) { l1View = v; l1Legend(); l1Pix = null; l1Render(); }
  l1Zoom = 4; l1Pix = null; l1Render(); l1Inspect({ clientX: 200, clientY: 200 });
  if (!document.getElementById("l1budget").innerHTML.includes("Drainage basins")) throw new Error("L1 summary not rendered");
  const sh = l1St.sh, bal = sh.balance(); if (Math.abs(bal.error) > 1e-9) throw new Error("L1 crust budget does not close: " + bal.error);
  console.log("PASS ui smoke L1: " + sh.nx + "x" + sh.ny + " sheet, " + sh.nm + " markers, t=" + sh.t.toFixed(1) + " Myr, " + l1St.bs.nBasins + " basins, budget error " + bal.error.toExponential(1));
  console.log("PASS ui smoke: t=" + world.t.toFixed(0) + " Myr, export " + (json.length / 1e6).toFixed(1) + " MB");
`)();
