// Runs the page's full script against a stubbed DOM: every view, both projections, inspector, extraction, export,
// and an L1 continent: pick, a short coupled run, every L1 view, zoom and inspect.
const fs = require("fs"), path = require("path");
const html = fs.readFileSync(path.join(__dirname, "..", "plates.html"), "utf8");
const code = html.slice(html.indexOf("<script>") + 8, html.lastIndexOf("</script>"));
const defaults = { grid: "64", seed: "genesis", nPlates: "10", mode: "scatter", speed: "3", contT: "50", l1Ref: "2", l1Win: "6", l1K: "-1", l1Reach: "1200" }, els = {};
function el(id) { return els[id] || (els[id] = { id, value: defaults[id] || "", checked: true, hidden: false, textContent: "", innerHTML: "", style: {}, className: "", dataset: {},
  clientWidth: 400, clientHeight: 400, width: 0, height: 0, children: [], options: id === "l1Ref" ? [{ value: "2" }, { value: "4" }] : [], disabled: false,
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
  for (let k = 0; k < 150; k++) world.step();
  for (const p of ["globe", "flat"]) { proj = p; pixIdx = null; for (const v of ["elev", "rock", "faults", "l1", "plates", "age", "thick"]) { view = v; colorsDirty = true; render(); } }
  inspect({ clientX: 200, clientY: 120 }); document.getElementById("extract").onclick(); render();
  const json = JSON.stringify(world.exportHandoff(handoff));
  if (!document.getElementById("summary").innerHTML.includes("formations")) throw new Error("summary not rendered");
  let land = -1; for (let i = 0; i < world.N; i++) if (world.wK[i] && world.elev[i] > 0.5) { land = i; break; }
  document.getElementById("l1Pick").onclick(); l1PickAt(land); render();
  document.getElementById("l1Run").onclick(); while (l1Look) l1LookTick();
  if (!region) throw new Error("L1 region not built: " + document.getElementById("l1Info").innerHTML);
  while (l1Running) l1Tick();
  for (const v of ["terrain", "erosion", "uplift", "basins"]) { l1View = v; l1Legend(); l1Render(); }
  l1Zoom = 4; l1Pix = null; l1Render(); l1Inspect({ clientX: 200, clientY: 200 });
  if (!document.getElementById("l1budget").innerHTML.includes("L0 against L1")) throw new Error("L1 summary not rendered");
  const bal = region.balance(); if (Math.abs(bal.error) > 1e-9) throw new Error("L1 crust budget does not close: " + bal.error);
  console.log("PASS ui smoke L1: " + region.G.nC + " cells, t=" + region.t.toFixed(1) + " Myr, budget error " + bal.error.toExponential(1));
  console.log("PASS ui smoke: t=" + world.t.toFixed(0) + " Myr, export " + (json.length / 1e6).toFixed(1) + " MB");
`)();
