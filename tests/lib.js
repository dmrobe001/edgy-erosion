// Shared helpers for the headless tests. Node.js 18+ with no dependencies.
const fs = require("fs"), path = require("path");
const SRC = path.join(__dirname, "..", "src", "sim.js");
function loadSim(file = SRC) { return new Function(fs.readFileSync(file, "utf8") + "\nreturn SIM;")(); }
// Load a patched copy of the simulation: each [from, to] pair is applied once; throws if a pattern is missing.
function loadSimPatched(patches, file = SRC) {
  let src = fs.readFileSync(file, "utf8");
  for (const [a, b] of patches) { if (!src.includes(a)) throw new Error("patch target not found: " + a.slice(0, 60)); src = src.replace(a, b); }
  return new Function(src + "\nreturn SIM;")();
}
function runTo(w, t) { while (w.t < t) w.step(); w.updateStats(true); return w; }
function arg(name, def) { const a = process.argv.find(s => s.startsWith("--" + name + "=")); return a ? a.split("=")[1] : def; }
const pct = x => (100 * x).toFixed(0) + "%";
module.exports = { loadSim, loadSimPatched, runTo, arg, pct };
