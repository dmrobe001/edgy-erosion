// Fault extraction, formations/domains, and JSON export: statistics plus basic schema checks.
const { loadSim, runTo } = require("./lib");
const SIM = loadSim(); const w = runTo(new SIM.World({ n: 96, seed: "genesis", nPlates: 10, mode: "scatter" }), 600);
const h = w.extractFaults(); h.l1 = w.buildL1(); const L = h.l1;
const kinds = {}; for (const f of h.faults) kinds[f.kindName] = (kinds[f.kindName] || 0) + 1;
console.log(`faults ${h.faults.length} ${JSON.stringify(kinds)}, active ${h.faults.filter(f => f.active).length}, mean length ${Math.round(h.faults.reduce((a, f) => a + f.lengthKm, 0) / h.faults.length)} km`);
console.log(`formations ${L.forms.length} (spanning 20+ columns: ${L.forms.filter(f => f.cells >= 20).length}), domains ${L.nDomains}`);
const out = w.exportHandoff(h), c = out.cells; let ok = true;
const check = (cond, msg) => { if (!cond) { ok = false; console.log("FAIL " + msg); } };
check(c.lat.length === w.N && c.columnOffsets.length === w.N + 1, "per-cell array lengths");
check(c.columnUnits.length === 4 * c.columnOffsets[w.N], "column units length matches offsets");
check(out.faults.every(f => f.trace.length >= 2), "every fault trace has 2+ points");
check(out.plates.length === w.plates.filter(p => p.alive).length, "plate list");
let maxF = -1; for (let k = 3; k < c.columnUnits.length; k += 4) maxF = Math.max(maxF, c.columnUnits[k]); check(maxF < out.formations.length, "formation ids in range");
const json = JSON.stringify(out); console.log(`export ${(json.length / 1e6).toFixed(1)} MB`);
console.log(ok ? "PASS handoff schema" : "FAIL handoff schema"); process.exitCode = ok ? 0 : 1;
