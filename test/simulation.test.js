import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Simulation, PRESETS, resolveParams } from '../src/index.js';

const small = (preset, extra = {}) => ({ preset, grid: { cols: 32, rows: 30 }, ...extra });

test('every preset runs and stays finite', () => {
  for (const name of Object.keys(PRESETS)) {
    const sim = new Simulation(small(name));
    sim.run(5);
    for (const a of [sim.rock, sim.soil, sim.z]) assert.ok(a.every(Number.isFinite), `${name}: non-finite value`);
    assert.ok(sim.soil.every((h) => h >= 0), `${name}: negative soil`);
  }
});

for (const name of ['mesa-country', 'coastal-plateau', 'alpine', 'rocky-coast', 'marine-terraces']) {
  test(`volume budget closes for ${name}: change = uplift - exported`, () => {
    const sim = new Simulation(small(name));
    const v0 = sim.volume();
    sim.run(30);
    const lhs = sim.volume() - v0;
    const rhs = sim.budget.uplift - sim.budget.exported;
    const scale = Math.max(Math.abs(sim.budget.uplift), Math.abs(sim.budget.exported), 1);
    assert.ok(Math.abs(lhs - rhs) < 1e-8 * scale, `change ${lhs} vs uplift-export ${rhs}`);
  });
}

test('same parameters give bit-identical results; a different seed does not', () => {
  const a = new Simulation(small('mesa-country')).run(15);
  const b = new Simulation(small('mesa-country')).run(15);
  const c = new Simulation(small('mesa-country', { seed: 2 })).run(15);
  assert.deepEqual(a.z, b.z);
  assert.notDeepEqual(a.z, c.z);
});

test('a snapshot restores the run exactly, including the random number stream', () => {
  const a = new Simulation(small('badlands')).run(10);
  const snap = JSON.parse(JSON.stringify(a.snapshot()));
  const b = Simulation.fromSnapshot(snap);
  a.run(10);
  b.run(10);
  assert.deepEqual(a.z, b.z);
  assert.equal(a.time, b.time);
});

test('parameters can change mid-run, but the grid cannot', () => {
  const sim = new Simulation(small('humid-hills'));
  sim.step();
  sim.setParams({ climate: { precipitation: 0.2 }, uplift: { rate: 0 } });
  sim.step();
  assert.equal(sim.params.climate.precipitation, 0.2);
  assert.throws(() => sim.setParams({ grid: { dx: 5 } }));
});

test('resolveParams layers defaults <- preset <- overrides and rejects nonsense', () => {
  const p = resolveParams({ preset: 'alpine', climate: { frost: 0.1 } });
  assert.equal(p.climate.frost, 0.1);
  assert.equal(p.climate.precipitation, PRESETS.alpine.climate.precipitation);
  assert.throws(() => resolveParams({ preset: 'nope' }), /unknown preset/);
  assert.throws(() => resolveParams({ strata: { layers: [{ rock: 'cheese', thickness: 1 }] } }), /unknown rock/);
  assert.throws(() => resolveParams({ boundary: { open: [], seaLevel: null } }), /nowhere to go/);
});

test('lattice spacing and extent are independent parameters', () => {
  const coarse = new Simulation({ grid: { cols: 20, rows: 20, dx: 50 } });
  const fine = new Simulation({ grid: { cols: 40, rows: 40, dx: 25 } });
  assert.ok(Math.abs(coarse.grid.width - fine.grid.width) / coarse.grid.width < 0.02);
  assert.equal(fine.grid.size, 4 * coarse.grid.size);
});

test('results converge as the time step shrinks', () => {
  // Mean elevation change over the same span with dt, dt/2, dt/4, dt/8: successive
  // differences must shrink. (Before "strip, then cut" incision they did not: a
  // sediment layer left by each step shielded the bed by an amount set by dt.)
  const base = new Simulation({ preset: 'alpine', grid: { cols: 48, rows: 44 } }).run(100);
  const snap = JSON.parse(JSON.stringify(base.snapshot()));
  const mean = (s) => s.z.reduce((a, v) => a + v, 0) / s.z.length;
  const T = 64000;
  const result = [2000, 1000, 500, 250].map((dt) => {
    const s = Simulation.fromSnapshot(snap);
    s.setParams({ time: { dt } });
    s.run(T / dt);
    return mean(s);
  });
  const d = [1, 2, 3].map((k) => Math.abs(result[k] - result[k - 1]));
  assert.ok(d[1] < 0.85 * d[0] && d[2] < 0.85 * d[1], `successive differences ${d.map((v) => v.toFixed(2)).join(', ')} m`);
});
