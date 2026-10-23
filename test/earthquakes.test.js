import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Simulation } from '../src/simulation.js';
import { sampleMagnitude, averageDisplacement, ruptureLength, faultRate } from '../src/processes/earthquakes.js';
import { Rng } from '../src/rng.js';
import { rockfall } from '../src/processes/rockfall.js';

// Everything but the faults switched off, so displacements can be measured directly.
const quiet = {
  grid: { cols: 48, rows: 40, dx: 20 },
  initial: { type: 'flat', height: 0, base: 100, noise: { amplitude: 0 }, soil: 1 },
  uplift: { rate: 0 },
  boundary: { open: ['west', 'east'] },
  fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false },
  rockfall: { enabled: false }, coast: { enabled: false },
  strata: { strengthVariation: { amplitude: 0 } },
};

test('magnitudes follow the Gutenberg-Richter law: maximum-likelihood b recovers the parameter', () => {
  const rng = new Rng(3);
  const b = 1.1, Mmin = 4, Mmax = 9; // wide range: truncation barely biases Aki's estimator
  let sum = 0;
  const N = 40000;
  for (let k = 0; k < N; k++) sum += sampleMagnitude(rng.next(), b, Mmin, Mmax);
  // Aki 1965: b = log10(e) / (mean M - Mmin)
  const bHat = Math.LOG10E / (sum / N - Mmin);
  assert.ok(Math.abs(bHat / b - 1) < 0.03, `b estimate ${bHat}`);
});

test('rupture scaling follows Wells & Coppersmith (1994)', () => {
  // M 7: surface rupture ~ 41 km, average displacement ~ 1.1 m
  assert.ok(Math.abs(ruptureLength(7) / 1000 - 10 ** (-3.22 + 4.83)) < 1e-9);
  assert.ok(Math.abs(averageDisplacement(7) - 10 ** (-4.8 + 4.83)) < 1e-12);
});

test('over many events the fault slips at its long-term slip rate', () => {
  const slipRate = 1e-3, dip = 50;
  const sim = new Simulation({
    ...quiet,
    time: { dt: 1e4 },
    earthquakes: {
      faults: [{ x0: 0.1, y0: 0.5, x1: 0.9, y1: 0.5, dip, slipRate, upliftShare: 1, decayLength: 1e9 }],
      minMagnitude: 5.5, maxMagnitude: 6.5,
    },
  });
  const z0 = Float64Array.from(sim.z);
  const T = 1e6;
  sim.run(T / 1e4);
  // average the raised-side uplift along the whole fault (rows just north of the trace)
  const g = sim.grid;
  const L = faultRate(sim, 0).length;
  let sum = 0, cnt = 0;
  for (let i = 0; i < g.size; i++) {
    const fx = (g.x[i] - sim.frame.x) / sim.frame.width, fy = (g.y[i] - sim.frame.y) / sim.frame.height;
    if (fx > 0.1 && fx < 0.9 && fy < 0.5 && fy > 0.3) { sum += sim.z[i] - z0[i]; cnt++; }
  }
  const expected = slipRate * T * Math.sin((dip * Math.PI) / 180);
  assert.ok(sim.quakes.count > 1000, `${sim.quakes.count} events`);
  assert.ok(Math.abs(sum / cnt / expected - 1) < 0.05, `mean throw ${(sum / cnt).toFixed(1)} m, expected ${expected.toFixed(1)} m (fault ${L.toFixed(0)} m)`);
  // and nothing moved on the other side, as upliftShare = 1 asked
  for (let i = 0; i < g.size; i++) {
    const fy = (g.y[i] - sim.frame.y) / sim.frame.height;
    if (fy > 0.55) assert.equal(sim.z[i], z0[i]);
  }
});

test('the fault offsets the strata: rock carries its layers up with it', () => {
  const sim = new Simulation({
    ...quiet,
    strata: { layers: [{ rock: 'sandstone', thickness: 2 }, { rock: 'shale', thickness: 2 }], top: 101, repeat: true, strengthVariation: { amplitude: 0 } },
    earthquakes: { faults: [{ x0: -1, y0: 0.5, x1: 2, y1: 0.5, dip: 90, slipRate: 1e-3, upliftShare: 1, decayLength: 1e9 }] },
  });
  const g = sim.grid;
  const i = g.index(24, 5); // raised (north) side
  const before = sim.lith.rockAt(i, sim.rock[i] - sim.upliftCum[i]);
  const ev = sim.earthquake(0, 7.0, 1.5 * sim.frame.width); // rupture centred on the domain
  assert.ok(ev.throw > 0.5);
  const lifted = sim.upliftCum[i];
  assert.ok(lifted > 0.5, `raised ${lifted} m`);
  // the same rock sits at the surface (it rose with its layer), unlike a cell that was merely eroded
  assert.equal(sim.lith.rockAt(i, sim.rock[i] - sim.upliftCum[i]), before);
});

test('a fresh fault scarp degrades by linear creep into an error-function profile', () => {
  // Hanks, Bucknam, Lajoie & Wallace 1984 (J. Geophys. Res. 89): z = H/2 (1 + erf(d / 2 sqrt(D t)))
  const D = 1, dx = 10, T = 900;
  const sim = new Simulation({
    ...quiet,
    grid: { cols: 20, rows: 60, dx },
    boundary: { open: ['south'] }, // far from the scarp, so no edge feels it
    initial: { type: 'flat', height: 0, base: 100, noise: { amplitude: 0 }, soil: 50 },
    hillslope: { enabled: true, diffusivity: D, criticalSlope: 100, maxNonlinearity: 1, transportDepth: 1e-3 },
    vegetation: { cover: 0 },
    time: { dt: 50 },
    earthquakes: { enabled: false, faults: [{ x0: -5, y0: 0.5, x1: 5, y1: 0.5, dip: 90, slipRate: 1e-3, upliftShare: 1, decayLength: 1e9 }] },
  });
  sim.earthquake(0, 7.0, 5 * sim.frame.width);
  const g = sim.grid;
  const col = 10;
  // the step lies on the face between the last raised row and the first one left behind
  let rUp = 0;
  while (sim.upliftCum[g.index(col, rUp + 1)] > 0) rUp++;
  const yStep = 0.5 * (g.y[g.index(col, rUp)] + g.y[g.index(col, rUp + 1)]);
  const top = sim.z[g.index(col, 5)], base = sim.z[g.index(col, g.rows - 5)];
  const H = top - base;
  assert.ok(H > 0.5, `throw ${H} m`);
  sim.run(T / 50);
  const erf = (x) => { // Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
    const s = Math.sign(x); x = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * x);
    return s * (1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x));
  };
  let err = 0, m = 0;
  for (let r = 5; r < g.rows - 5; r++) {
    const i = g.index(col, r);
    const d = yStep - g.y[i]; // > 0 on the raised (north) side
    const model = base + (H / 2) * (1 + erf(d / (2 * Math.sqrt(D * T))));
    err += (sim.z[i] - model) ** 2; m++;
  }
  const rms = Math.sqrt(err / m);
  assert.ok(rms < 0.01 * H, `rms misfit ${rms.toFixed(4)} m of a ${H.toFixed(2)} m scarp`);
});

test('earthquake runs are deterministic and survive a snapshot', () => {
  const p = { preset: 'earthquake-scarp', grid: { cols: 40, rows: 40 } };
  const a = new Simulation(p);
  a.run(12);
  const snap = JSON.parse(JSON.stringify(a.snapshot()));
  a.run(12);
  const b = Simulation.fromSnapshot(snap);
  b.run(12);
  const c = new Simulation(p);
  c.run(24);
  assert.ok(a.quakes.count > 0);
  assert.deepEqual(Array.from(b.z), Array.from(a.z));
  assert.deepEqual(Array.from(c.z), Array.from(a.z));
  assert.equal(b.quakes.count, a.quakes.count);
});

test('shaking brings down cliffs that were close to failing', () => {
  // an 18 m step in weak shale over a 10 m cell: Culmann critical height ~23 m, so it
  // stands; strong shaking (weakening 0.6) cuts that to ~9 m and it fails
  const base = {
    ...quiet,
    grid: { cols: 30, rows: 30, dx: 10 },
    rockfall: { enabled: true, timescale: 1e12 }, // spontaneous failure practically never
    strata: { layers: [{ rock: 'shale', thickness: 1e9 }], strengthVariation: { amplitude: 0 } },
    earthquakes: { faults: [{ x0: -2, y0: 0.95, x1: 3, y1: 0.95, dip: 60, slipRate: 0, decayLength: 1e9 }], coseismicWeakening: 0.6 },
  };
  const make = () => {
    const s = new Simulation(base);
    for (let i = 0; i < s.grid.size; i++) if (s.grid.col(i) < 15) s.rock[i] += 18;
    s.setState({});
    return s;
  };
  const calm = make();
  assert.equal(rockfall(calm, 1).events, 0, 'stable without shaking');
  const shaken = make();
  shaken.earthquake(0, 7, 2.5 * shaken.frame.width); // shaking is 1 near the rupture, ~0.9 across the domain
  assert.ok(Math.min(...shaken.shaking) > 0.8);
  const r = rockfall(shaken, 1);
  assert.ok(r.events >= 20, `${r.events} cliff failures under shaking`);
});
