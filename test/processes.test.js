import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Simulation } from '../src/simulation.js';
import { creep, avalanche } from '../src/processes/hillslope.js';
import { rockfall, culmannHeight } from '../src/processes/rockfall.js';

const sum = (a) => a.reduce((s, v) => s + v, 0);

// Only the processes under test are switched on.
const bare = (extra = {}) => ({
  grid: { cols: 30, rows: 30, dx: 100 },
  initial: { type: 'flat', base: 0, height: 0, noise: { amplitude: 5, wavelength: 500 }, soil: 0 },
  strata: { layers: [{ rock: 'granite', thickness: 1e9 }], strengthVariation: { amplitude: 0 } },
  fluvial: { routing: 'sfd', randomReceivers: false, deposition: 0 },
  hillslope: { enabled: false },
  weathering: { enabled: false },
  rockfall: { enabled: false },
  ...extra,
});

for (const n of [1, 2]) {
  test(`stream power reaches the analytical steady state K A^m S^n = U (n = ${n})`, () => {
    const U = 1e-4;
    const sim = new Simulation(bare({ uplift: { rate: U }, time: { dt: 2e4 }, fluvial: { routing: 'sfd', randomReceivers: false, deposition: 0, n } }));
    let converged = false;
    for (let s = 0; s < 5000 && !converged; s++) {
      sim.step();
      converged = sim.dzdt.every((v) => Math.abs(v) < 1e-12);
    }
    assert.ok(converged, 'landscape reached steady state');
    const K = sim.lith.erodibility[sim.lith.names.indexOf('granite')];
    const { receiver, area } = sim.flow;
    let worst = 0;
    for (let i = 0; i < sim.grid.size; i++) {
      if (sim.outlet[i]) continue;
      const S = (sim.z[i] - sim.z[receiver[i]]) / sim.grid.dx;
      const E = K * Math.sqrt(area[i]) * Math.pow(S, n);
      worst = Math.max(worst, Math.abs(E - U) / U);
    }
    assert.ok(worst < 1e-6, `worst relative error ${worst}`);
  });
}

test('linear creep spreads a bump with variance growing at 2 D t (hex Laplacian is isotropic)', () => {
  const D = 0.05;
  const sim = new Simulation(bare({ grid: { cols: 61, rows: 61, dx: 10 }, hillslope: { enabled: true, diffusivity: D, criticalSlope: 1e9, transportDepth: 1e-9 } }));
  const g = sim.grid;
  sim.boundaryOutlet.fill(0);
  sim.veg.fill(0);
  const xc = g.x[g.index(30, 30)], yc = g.y[g.index(30, 30)];
  const s0 = 40;
  for (let i = 0; i < g.size; i++) {
    sim.rock[i] = 0;
    sim.soil[i] = 10 + Math.exp(-((g.x[i] - xc) ** 2 + (g.y[i] - yc) ** 2) / (2 * s0 * s0));
    sim.z[i] = sim.soil[i];
  }
  const moments = () => {
    let m0 = 0, mx = 0, my = 0;
    for (let i = 0; i < g.size; i++) {
      const h = sim.soil[i] - 10;
      m0 += h; mx += h * (g.x[i] - xc) ** 2; my += h * (g.y[i] - yc) ** 2;
    }
    return { m0, vx: mx / m0, vy: my / m0 };
  };
  const before = moments();
  const T = 20000;
  creep(sim, T);
  const after = moments();
  assert.ok(Math.abs(after.m0 - before.m0) < 1e-9 * before.m0, 'volume conserved');
  for (const [b, a] of [[before.vx, after.vx], [before.vy, after.vy]]) {
    const growth = a - b;
    assert.ok(Math.abs(growth - 2 * D * T) / (2 * D * T) < 0.01, `variance grew by ${growth}, expected ${2 * D * T}`);
  }
});

test('avalanching leaves no soil slope steeper than the angle of repose, and conserves volume', () => {
  const sim = new Simulation(bare({ grid: { cols: 25, rows: 25, dx: 10 }, hillslope: { enabled: true, criticalSlope: 0.7, avalancheIterations: 5000 }, vegetation: { cover: 0 } }));
  const g = sim.grid;
  sim.boundaryOutlet.fill(0);
  sim.veg.fill(0);
  sim.rock.fill(0);
  sim.soil.fill(0);
  sim.soil[g.index(12, 12)] = 200; // a tall column of debris
  for (let i = 0; i < g.size; i++) sim.z[i] = sim.soil[i];
  const v0 = sum(sim.soil);
  avalanche(sim);
  assert.ok(Math.abs(sum(sim.soil) - v0) < 1e-9 * v0);
  let steepest = 0;
  for (let i = 0; i < g.size; i++) {
    if (sim.soil[i] <= 0) continue;
    for (let k = 0; k < 6; k++) {
      const j = g.neighbors[6 * i + k];
      if (j >= 0) steepest = Math.max(steepest, (sim.z[i] - sim.z[j]) / g.dx);
    }
  }
  assert.ok(steepest <= 0.7 + 1e-4, `steepest soil slope ${steepest}`);
});

test('Culmann critical height: infinite below the friction angle, finite and decreasing above', () => {
  const phi = (35 * Math.PI) / 180;
  assert.equal(culmannHeight(phi * 0.99, 1e5, phi, 2600), Infinity);
  const h60 = culmannHeight((60 * Math.PI) / 180, 1e5, phi, 2600);
  const h85 = culmannHeight((85 * Math.PI) / 180, 1e5, phi, 2600);
  assert.ok(Number.isFinite(h60) && h85 < h60);
  // closed form check at beta = 90 deg: Hc = 4C/(rho g) * cos(phi) / (1 - sin(phi))
  const h90 = culmannHeight(Math.PI / 2, 1e5, phi, 2600);
  assert.ok(Math.abs(h90 - ((4e5 / (2600 * 9.81)) * Math.cos(phi)) / (1 - Math.sin(phi))) < 1e-9);
});

for (const [rock, shouldFail] of [['shale', true], ['granite', false]]) {
  test(`a 150 m vertical cliff of ${rock} ${shouldFail ? 'collapses' : 'stands'}`, () => {
    const sim = new Simulation(bare({
      grid: { cols: 20, rows: 10, dx: 10 },
      strata: { layers: [{ rock, thickness: 1e9 }], top: 1000, strengthVariation: { amplitude: 0 } },
      rockfall: { enabled: true, timescale: 1e-3 },
    }));
    const g = sim.grid;
    for (let i = 0; i < g.size; i++) {
      sim.rock[i] = g.col(i) >= 10 ? 150 : 0;
      sim.soil[i] = 0;
      sim.z[i] = sim.rock[i];
    }
    sim.boundaryOutlet.fill(0);
    const v0 = sum(sim.rock) + sum(sim.soil);
    const { events } = rockfall(sim, 100);
    assert.ok(Math.abs(sum(sim.rock) + sum(sim.soil) - v0) < 1e-9 * v0, 'volume conserved');
    if (shouldFail) assert.ok(events > 0 && sum(sim.soil) > 0);
    else assert.equal(events, 0);
  });
}
