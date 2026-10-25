// Runs the page's full script against a stubbed DOM: every view, both projections, inspector, extraction, export.
const fs = require("fs"), path = require("path");
const html = fs.readFileSync(path.join(__dirname, "..", "plates.html"), "utf8");
const code = html.slice(html.indexOf("<script>") + 8, html.lastIndexOf("</script>"));
const defaults = { grid: "64", seed: "genesis", nPlates: "10", mode: "scatter", speed: "3", contT: "50" }, els = {};
function el(id) { return els[id] || (els[id] = { id, value: defaults[id] || "", checked: true, hidden: false, textContent: "", innerHTML: "", style: {}, className: "", dataset: {},
  clientWidth: 400, clientHeight: 400, width: 0, height: 0, children: [], addEventListener() {}, setAttribute() {}, setPointerCapture() {}, prepend(c) { this.children.unshift(c); },
  get lastChild() { return { remove() {} }; }, getBoundingClientRect() { return { left: 0, top: 0, width: 400, height: 400 }; },
  getContext() { return { createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) }), putImageData() {}, clearRect() {}, beginPath() {}, arc() {}, stroke() {}, fill() {},
    moveTo() {}, lineTo() {}, closePath() {}, setLineDash() {} }; } }); }
global.document = { getElementById: el, querySelectorAll: () => [], createElement: () => ({}), documentElement: {} };
global.getComputedStyle = () => ({ getPropertyValue: () => "#d2452f" });
global.window = { devicePixelRatio: 2, addEventListener() {} };
global.requestAnimationFrame = () => {};
new Function(code + `
  for (let k = 0; k < 150; k++) world.step();
  for (const p of ["globe", "flat"]) { proj = p; pixIdx = null; for (const v of ["elev", "rock", "faults", "l1", "plates", "age", "thick"]) { view = v; colorsDirty = true; render(); } }
  inspect({ clientX: 200, clientY: 120 }); document.getElementById("extract").onclick(); render();
  const json = JSON.stringify(world.exportHandoff(handoff));
  if (!document.getElementById("summary").innerHTML.includes("formations")) throw new Error("summary not rendered");
  console.log("PASS ui smoke: t=" + world.t.toFixed(0) + " Myr, export " + (json.length / 1e6).toFixed(1) + " MB");
`)();
