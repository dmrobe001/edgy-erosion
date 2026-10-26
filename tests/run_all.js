// Quick pass over every test with short run lengths. Full-length runs: call each file directly (see README).
const { execFileSync } = require("child_process"), path = require("path");
const runs = [["test_ui_smoke.js"], ["test_netrotation.js", "--t=300"], ["test_budget.js", "--t=800"], ["test_handoff.js"], ["test_l1.js"], ["test_flicker.js"], ["test_volume_budget.js"], ["bench.js", "--t=60"]];
let failed = 0;
for (const [f, ...a] of runs) {
  console.log(`\n=== ${f} ${a.join(" ")}`);
  try { process.stdout.write(execFileSync(process.execPath, [path.join(__dirname, f), ...a], { encoding: "utf8" })); }
  catch (e) { failed++; process.stdout.write((e.stdout || "") + (e.stderr || "")); }
}
console.log(failed ? `\n${failed} test(s) failed` : "\nall tests ran"); process.exitCode = failed ? 1 : 0;
