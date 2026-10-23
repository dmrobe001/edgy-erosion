import { test } from 'node:test';
import assert from 'node:assert/strict';
import { HexGrid } from '../src/hexgrid.js';
import { GraphMesh, checkMesh } from '../src/mesh.js';
import { Simulation } from '../src/simulation.js';
import { creep } from '../src/processes/hillslope.js';

const sum = (a, w) => a.reduce((s, v, i) => s + v * (w ? w[i] : 1), 0);

test('the hex grid satisfies the mesh contract, also when offset in world coordinates', () => {
  const g = new HexGrid(12, 9, 7, 1000, -500);
  assert.ok(checkMesh(g));
  for (let i = 0; i < g.size; i++) assert.equal(g.cellAt(g.x[i], g.y[i]), i);
  assert.ok(g.bounds.xmin < g.x[0] && g.x[0] - g.bounds.xmin === 3.5);
});

test('an irregular mesh (jittered, 2 to 6 neighbours, gaps between slots) satisfies the contract', () => {
  const m = GraphMesh.irregular(new HexGrid(30, 20, 10), { seed: 4 });
  assert.ok(checkMesh(m));
  const degrees = new Set();
  for (let i = 0; i < m.size; i++) {
    let d = 0;
    for (let k = 0; k < m.degree; k++) if (m.neighbors[m.degree * i + k] >= 0) d++;
    degrees.add(d);
  }
  assert.ok(degrees.has(2) || degrees.has(3), 'some cells lost neighbours');
  assert.ok(degrees.has(6), 'some cells kept all six');
});

const bare = {
  grid: { cols: 40, rows: 30, dx: 20 },
  initial: { type: 'dome', height: 80, noise: { amplitude: 5, wavelength: 200 }, soil: 0.5 },
  strata: { layers: [{ rock: 'siltstone', thickness: 1e9 }] },
  coast: { enabled: false },
  rockfall: { timescale: 1 },
  uplift: { rate: 2e-4 },
};

test('every process runs on an irregular mesh, water drains everywhere and the volume budget closes', () => {
  const hex = new HexGrid(bare.grid.cols, bare.grid.rows, bare.grid.dx);
  const mesh = GraphMesh.irregular(hex, { seed: 9 });
  const sim = new Simulation(bare, { mesh });
  const v0 = sim.volume();
  sim.run(40);
  for (const a of [sim.rock, sim.soil, sim.z]) assert.ok(a.every(Number.isFinite));
  assert.ok(sim.soil.every((h) => h >= 0));
  const { receiver, filled } = sim.flow;
  for (let i = 0; i < mesh.size; i++) {
    if (sim.outlet[i]) continue;
    assert.ok(receiver[i] >= 0 && filled[receiver[i]] < filled[i]);
  }
  const lhs = sim.volume() - v0;
  const rhs = sim.budget.uplift + sim.budget.inflow - sim.budget.exported;
  assert.ok(Math.abs(lhs - rhs) < 1e-8 * Math.max(sim.budget.uplift, 1), `change ${lhs} vs budget ${rhs}`);
});

test('creep conserves volume on an irregular mesh with unequal cell areas', () => {
  const mesh = GraphMesh.irregular(new HexGrid(25, 25, 10), { seed: 2, perturb: 0.6 });
  const sim = new Simulation({ ...bare, grid: { cols: 25, rows: 25, dx: 10 }, hillslope: { diffusivity: 0.5 } }, { mesh });
  sim.boundaryOutlet.fill(0);
  const v0 = sum(sim.z, mesh.area);
  creep(sim, 5000);
  assert.ok(Math.abs(sum(sim.z, mesh.area) - v0) < 1e-9 * Math.abs(v0));
  assert.ok(sim.soil.every((h) => h >= 0));
});

test('the coast process refuses a mesh without the hex lattice, with a clear message', () => {
  const mesh = GraphMesh.irregular(new HexGrid(20, 20, 10), { seed: 1 });
  const sim = new Simulation({ ...bare, grid: { cols: 20, rows: 20, dx: 10 }, boundary: { seaLevel: 0 }, coast: { enabled: true } }, { mesh });
  assert.throws(() => sim.step(), /hexagonal-lattice/);
});
