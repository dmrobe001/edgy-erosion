import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Simulation } from '../src/simulation.js';
import { seaLevelAt, latticeWaveWeights } from '../src/processes/coast.js';
import { resolveParams } from '../src/params.js';

// Only the coast and uplift act; land is uniform rock with no noise.
const coastOnly = (extra = {}) => {
  const base = {
    grid: { cols: 60, rows: 6, dx: 10 },
    time: { dt: 10 },
    initial: { type: 'flat', base: 20, height: 0, noise: { amplitude: 0 }, soil: 0 },
    uplift: { rate: 0 },
    boundary: { open: ['east'], seaLevel: 0 },
    strata: { layers: [{ rock: 'sandstone', thickness: 1e9 }], strengthVariation: { amplitude: 0 } },
    fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false }, rockfall: { enabled: false },
    coast: { waves: [{ from: 0, weight: 1 }], refraction: 0, debrisRetained: 0, marineDiffusivity: 0, dissipationLength: 1e12 },
  };
  return resolveParams({ ...base, ...extra, coast: { ...base.coast, ...(extra.coast || {}) } });
};

/** Simulation of a vertical sea cliff: land at +20 m for col < 30, sea floor at -20 m beyond. */
function cliff(params) {
  const sim = new Simulation(params);
  const g = sim.grid;
  for (let i = 0; i < g.size; i++) {
    sim.rock[i] = g.col(i) < 30 ? 20 : -20;
    sim.soil[i] = 0;
    sim.z[i] = sim.rock[i];
  }
  return sim;
}

const shoreline = (sim, row = 3) => {
  const g = sim.grid;
  let c = 0;
  while (c < g.cols && sim.z[row * g.cols + c] >= sim.seaLevel) c++;
  return c * g.dx; // distance from the west edge to the first sea cell (m)
};

test('sea level follows the sinusoid, or a piecewise-linear curve when given', () => {
  const p = resolveParams({ boundary: { seaLevel: 5 }, coast: { seaLevelAmplitude: 30, seaLevelPeriod: 40000, seaLevelPhase: 90 } });
  assert.ok(Math.abs(seaLevelAt(p, 0) - 35) < 1e-9);
  assert.ok(Math.abs(seaLevelAt(p, 20000) + 25) < 1e-9);
  const q = resolveParams({ boundary: { seaLevel: 0 }, coast: { seaLevelCurve: [[0, -100], [10000, 0], [20000, 0]] } });
  assert.equal(seaLevelAt(q, 5000), -50);
  assert.equal(seaLevelAt(q, 50000), 0);
  assert.equal(seaLevelAt(resolveParams({}), 0), null);
});

test('wave climate is split onto the six lattice directions and normalised', () => {
  const w = latticeWaveWeights([{ from: 30, weight: 2 }, { from: 180, weight: 2 }]);
  assert.deepEqual(Array.from(w).map((v) => +v.toFixed(6)), [0.25, 0.25, 0, 0.5, 0, 0]);
});

test('with no wave damping, a cliff retreats at exactly its wave erodibility', () => {
  const sim = cliff(coastOnly());
  const R = sim.lith.waveErodibility[sim.lith.names.indexOf('sandstone')];
  const x0 = shoreline(sim);
  const T = 2000;
  sim.run(T / sim.params.time.dt);
  const retreat = x0 - shoreline(sim);
  // retreat is quantised to whole cells: allow one cell either way
  assert.ok(Math.abs(retreat - R * T) <= sim.grid.dx + 1e-9, `retreat ${retreat} m, expected ${R * T} m`);
});

test('a widening shallow platform slows cliff retreat (the wave-dissipation feedback)', () => {
  const free = cliff(coastOnly());
  const damped = cliff(coastOnly({ coast: { dissipationLength: 50 } }));
  free.run(300);
  damped.run(300);
  assert.ok(shoreline(damped) > shoreline(free) + 2 * free.grid.dx, 'damped cliff retreats less');
});

test('a beach in front of a cliff shields it', () => {
  const bare = cliff(coastOnly());
  const beach = cliff(coastOnly({ coast: { attrition: 0 } }));
  const g = beach.grid;
  for (let i = 0; i < g.size; i++) if (g.col(i) === 29) { beach.soil[i] = 3; beach.z[i] += 3; }
  bare.run(100);
  beach.run(100);
  assert.ok(shoreline(beach) > shoreline(bare), 'shielded cliff retreats less');
});

test('a sheltered lagoon gets less wave energy than the open coast', () => {
  const sim = new Simulation(coastOnly({ grid: { cols: 40, rows: 30, dx: 10 }, coast: { waves: [{ from: 0, weight: 1 }, { from: 60, weight: 1 }, { from: 300, weight: 1 }] } }));
  const g = sim.grid;
  // mainland at col < 10; a barrier at col 25 joined to it by land along the top and
  // bottom rows encloses a lagoon; open sea east of the barrier
  for (let i = 0; i < g.size; i++) {
    const c = g.col(i), r = g.row(i);
    const land = c < 10 || (c === 25) || ((r === 0 || r === g.rows - 1) && c <= 25);
    sim.rock[i] = land ? 5 : -10;
    sim.soil[i] = 0;
    sim.z[i] = sim.rock[i];
  }
  sim.step();
  const lagoonShore = sim.waveExposureField[g.index(9, 15)];
  const openShore = sim.waveExposureField[g.index(26, 15)];
  assert.ok(openShore > 0.9, `open shore exposure ${openShore}`);
  assert.ok(lagoonShore < 0.5 * openShore, `lagoon shore exposure ${lagoonShore}`);
});

test('marine terraces: shoreline angles sit at highstand level plus uplift since then', () => {
  const P = 41000, U = 1e-3, A = 30;
  const sim = new Simulation({
    preset: 'marine-terraces', grid: { cols: 128, rows: 6 },
    initial: { noise: { amplitude: 0 } }, strata: { strengthVariation: { amplitude: 0 } },
    uplift: { rate: U }, coast: { seaLevelAmplitude: A, seaLevelPeriod: P, seaLevelPhase: 90 },
    fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false }, rockfall: { enabled: false },
  });
  const T = 4 * P;
  sim.run(Math.round(T / sim.params.time.dt));
  const g = sim.grid;
  const z = Array.from({ length: g.cols }, (_, c) => sim.z[3 * g.cols + c]);
  // shoreline angle = foot of a riser (a drop of more than 5 m between neighbours), above sea
  const angles = [];
  for (let c = 1; c < g.cols; c++) if (z[c - 1] - z[c] > 5 && z[c] > sim.seaLevel + 5) angles.push(z[c]);
  const halfTide = sim.params.coast.tidalRange / 2;
  // the three youngest uplifted highstands (at t = T - P, T - 2P, T - 3P)
  for (let k = 1; k <= 3; k++) {
    const predicted = A - halfTide + U * k * P;
    const nearest = Math.min(...angles.map((a) => Math.abs(a - predicted)));
    // Terraces are U*P = 41 m apart. Each leaves one or two risers, and platforms are cut
    // just below the peak, so allow 6 m; a 10% error in period or uplift misses by more.
    assert.ok(nearest < 6, `no shoreline angle near ${predicted.toFixed(1)} m; found ${angles.map((a) => a.toFixed(0)).join(', ')}`);
  }
});
