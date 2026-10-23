import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Simulation } from '../src/simulation.js';
import { zoomIn, compareWithParent, interpolator } from '../src/nest.js';

// One parent shared by the tests below: a granite dome, uplifted and dissected.
const parent = new Simulation({ preset: 'alpine' });
parent.run(300);
const g = parent.grid;
let trunk = -1; // the biggest river away from the domain edges
for (let i = 0; i < g.size; i++) {
  const c = g.col(i), r = g.row(i);
  if (c < 40 || c > 88 || r < 40 || r > 88) continue;
  if (trunk < 0 || parent.flow.area[i] > parent.flow.area[trunk]) trunk = i;
}
const zoom = (extra = {}) => zoomIn(parent, { center: [g.x[trunk], g.y[trunk]], factor: 4, cols: 96, rows: 96, ...extra });

test('interpolation reproduces a linear field exactly', () => {
  const f = Float64Array.from(g.x, (x, i) => 3 * x - 2 * g.y[i] + 7);
  const at = interpolator(g, f);
  for (const [x, y] of [[1000.3, 900.7], [2222.2, 1500.1], [3000, 2000.5]]) {
    assert.ok(Math.abs(at(x, y) - (3 * x - 2 * y + 7)) < 1e-6);
  }
});

test('the zoomed patch preserves the parent elevation averaged over each parent cell', () => {
  const child = zoom();
  const cg = child.grid;
  const s = new Float64Array(g.size), c = new Float64Array(g.size);
  for (let i = 0; i < cg.size; i++) { const p = g.cellAt(cg.x[i], cg.y[i]); s[p] += child.z[i]; c[p]++; }
  let worst = 0, covered = 0;
  for (let p = 0; p < g.size; p++) {
    if (c[p] < 0.9 * 16) continue;
    covered++;
    worst = Math.max(worst, Math.abs(s[p] / c[p] - parent.z[p]));
  }
  assert.ok(covered > 200);
  assert.ok(worst < 0.05, `worst block-mean misfit ${worst} m`);
  // and the added detail is real: the child is rougher than plain interpolation
  const smooth = zoom({ roughness: 0 });
  let diff = 0;
  for (let i = 0; i < cg.size; i++) diff = Math.max(diff, Math.abs(child.z[i] - smooth.z[i]));
  assert.ok(diff > 1, 'roughness added detail');
});

test('geology lines up across scales: same strata planes, plugs and weak zones', () => {
  const p = new Simulation({ preset: 'volcanic-necks', strata: { undulation: { amplitude: 0 }, dip: 3, dipDirection: 40 } });
  const child = zoomIn(p, { center: [0.38 * p.grid.width, 0.42 * p.grid.height], factor: 4, cols: 64, rows: 64 });
  const top = interpolator(p.grid, p.lith.columnTop); // the column top is planar here: linear interpolation is exact
  const cg = child.grid;
  let worstTop = 0, plugMismatch = 0, weakWorst = 0;
  for (let i = 0; i < cg.size; i++) {
    worstTop = Math.max(worstTop, Math.abs(child.lith.columnTop[i] - top(cg.x[i], cg.y[i])));
    const q = p.grid.cellAt(cg.x[i], cg.y[i]);
    // deep inside the plug in the parent (cell and all neighbours basalt) => basalt in the child
    let deep = p.lith.bodyRock[q] >= 0;
    for (let k = 0; k < 6; k++) { const r = p.grid.neighbors[6 * q + k]; if (r < 0 || p.lith.bodyRock[r] < 0) deep = false; }
    if (deep && child.lith.bodyRock[i] < 0) plugMismatch++;
    if (cg.x[i] === p.grid.x[q] && cg.y[i] === p.grid.y[q]) weakWorst = Math.max(weakWorst, Math.abs(child.lith.weakness[i] - p.lith.weakness[q]));
  }
  assert.ok(worstTop < 1e-6, `strata top differs by ${worstTop} m`);
  assert.equal(plugMismatch, 0);
  assert.ok(weakWorst < 1e-12);
});

test('rivers entering from the parent leave the patch with the parent discharge', () => {
  const child = zoom();
  const c = compareWithParent(child, parent);
  const { parent: pa, child: ca } = c.largestExitArea;
  assert.ok(Math.abs(ca - pa) / pa < 0.1, `largest exit: child ${ca}, parent ${pa}`);
  assert.ok(child.zoom.enteringArea > 0.5 * pa, 'most of that water entered from outside the patch');
});

test('a same-resolution patch reproduces the parent run on at the same time step', () => {
  // The control that separates handoff errors from resolution and time-step effects:
  // zoom with factor 1, run both for the same time with the same dt, compare patch means.
  const dt = 500, T = 50000;
  const child = zoomIn(parent, { center: [g.x[trunk], g.y[trunk]], factor: 1, cols: 40, rows: 40, dt, roughness: 0 });
  const twin = Simulation.fromSnapshot(JSON.parse(JSON.stringify(parent.snapshot())));
  twin.setParams({ time: { dt } });
  const b = child.grid.bounds;
  const m = 2 * g.dx; // stay away from the patch edges
  const patchMean = (sim) => {
    let s = 0, c = 0;
    for (let i = 0; i < sim.grid.size; i++) {
      const x = sim.grid.x[i], y = sim.grid.y[i];
      if (x > b.xmin + m && x < b.xmax - m && y > b.ymin + m && y < b.ymax - m) { s += sim.z[i]; c++; }
    }
    return s / c;
  };
  const z0 = patchMean(twin);
  const v0 = child.volume();
  child.run(T / dt);
  twin.run(T / dt);
  const relief = parent.stats().maxElevation - parent.stats().minElevation;
  const dChild = patchMean(child) - z0, dTwin = patchMean(twin) - z0;
  assert.ok(Math.abs(dChild - dTwin) < 0.01 * relief, `patch mean change: child ${dChild.toFixed(1)} m, parent run on ${dTwin.toFixed(1)} m`);
  // and the child's own volume budget closes (uplift and edge motion + inflow - export)
  const lhs = child.volume() - v0;
  const rhs = child.budget.uplift + child.budget.inflow - child.budget.exported;
  assert.ok(Math.abs(lhs - rhs) < 1e-8 * Math.max(Math.abs(child.budget.uplift), 1), `budget ${lhs} vs ${rhs}`);
});

test('patch edges follow the parent\'s last rate of change exactly', () => {
  const child = zoom();
  const edge = [];
  for (let i = 0; i < child.grid.size; i++) if (child.boundaryOutlet[i]) edge.push([i, child.rock[i], child.boundaryRate[i]]);
  child.run(20);
  const t = 20 * child.params.time.dt;
  for (const [i, r0, rate] of edge) assert.ok(Math.abs(child.rock[i] - (r0 + rate * t)) < 1e-9 * (1 + Math.abs(r0)));
  assert.ok(edge.some(([, , rate]) => Math.abs(rate) > 0), 'some edges move');
});

test('zooming twice works and records the lineage', () => {
  const child = zoom();
  child.run(20);
  const cg = child.grid;
  const grand = zoomIn(child, { center: [cg.x[cg.size >> 1], cg.y[cg.size >> 1]], factor: 3, cols: 48, rows: 48 });
  assert.equal(grand.zoom.level, 2);
  assert.ok(Math.abs(grand.grid.dx - parent.grid.dx / 12) < 1e-12);
  grand.run(5);
  assert.ok(grand.z.every(Number.isFinite));
});
