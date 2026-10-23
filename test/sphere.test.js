import { test } from 'node:test';
import assert from 'node:assert/strict';
import { SphereMesh, azimuthalProjection, lonLatToVector } from '../src/sphere.js';
import { checkMesh } from '../src/mesh.js';
import { Simulation } from '../src/simulation.js';
import { Rng } from '../src/rng.js';

test('the geodesic mesh has 10 f^2 + 2 cells, 12 pentagons, and tiles the sphere exactly', () => {
  for (const f of [3, 8]) {
    const m = new SphereMesh(f, 2);
    assert.equal(m.size, 10 * f * f + 2);
    assert.ok(checkMesh(m));
    let pent = 0;
    for (let i = 0; i < m.size; i++) {
      assert.ok(m.valence[i] === 5 || m.valence[i] === 6);
      if (m.valence[i] === 5) pent++;
    }
    assert.equal(pent, 12);
    const area = m.area.reduce((s, a) => s + a, 0);
    assert.ok(Math.abs(area / (4 * Math.PI * 4) - 1) < 1e-12, `area ${area}`);
  }
});

/** Relative RMS error of the finite-volume Laplacian against the exact -l(l+1)/R^2 f. */
function laplacianError(m, fn, l) {
  let num = 0, den = 0;
  const R = m.radius;
  for (let i = 0; i < m.size; i++) {
    let s = 0;
    for (let k = 0; k < 6; k++) {
      const e = 6 * i + k, j = m.neighbors[e];
      if (j >= 0) s += (m.face[e] * (fn(m, j) - fn(m, i))) / m.dist[e];
    }
    s /= m.area[i];
    const exact = (-l * (l + 1) * fn(m, i)) / (R * R);
    num += (s - exact) ** 2 * m.area[i];
    den += exact ** 2 * m.area[i];
  }
  return Math.sqrt(num / den);
}

test('the discrete Laplacian converges to the spherical-harmonic eigenvalues at second order', () => {
  const Y10 = (m, i) => m.pos[3 * i + 2];
  const Y22 = (m, i) => m.pos[3 * i] * m.pos[3 * i + 1];
  for (const [fn, l] of [[Y10, 1], [Y22, 2]]) {
    const e8 = laplacianError(new SphereMesh(8, 1000), fn, l);
    const e16 = laplacianError(new SphereMesh(16, 1000), fn, l);
    assert.ok(e16 < 0.005, `l=${l}: error ${e16}`);
    assert.ok(e8 / e16 > 3.3, `l=${l}: halving the spacing cut the error only ${(e8 / e16).toFixed(2)}x`);
  }
});

test('nearest-cell search by greedy walk agrees with brute force, from any start', () => {
  const m = new SphereMesh(12);
  const rng = new Rng(5);
  for (let t = 0; t < 400; t++) {
    const z = 2 * rng.next() - 1, a = 2 * Math.PI * rng.next(), r = Math.sqrt(1 - z * z);
    const v = [r * Math.cos(a), r * Math.sin(a), z];
    let best = 0, bd = -2;
    for (let i = 0; i < m.size; i++) {
      const d = m.pos[3 * i] * v[0] + m.pos[3 * i + 1] * v[1] + m.pos[3 * i + 2] * v[2];
      if (d > bd) { bd = d; best = i; }
    }
    assert.equal(m.nearest(v[0], v[1], v[2]), best);
    assert.equal(m.nearest(v[0], v[1], v[2], Math.floor(rng.next() * m.size)), best);
  }
  // map coordinates round-trip
  for (let i = 0; i < m.size; i += 37) assert.equal(m.cellAt(m.x[i], m.y[i]), i);
});

test('linear interpolation weights sum to one and reproduce linear fields to O(h^2)', () => {
  const m = new SphereMesh(16);
  const rng = new Rng(9);
  const a = [0.3, -0.5, 0.8];
  const f = (i) => a[0] * m.pos[3 * i] + a[1] * m.pos[3 * i + 1] + a[2] * m.pos[3 * i + 2];
  let worst = 0;
  for (let t = 0; t < 300; t++) {
    const z = 2 * rng.next() - 1, ang = 2 * Math.PI * rng.next(), r = Math.sqrt(1 - z * z);
    const v = [r * Math.cos(ang), r * Math.sin(ang), z];
    const { ids, weights } = m.locate(v[0], v[1], v[2]);
    assert.ok(Math.abs(weights[0] + weights[1] + weights[2] - 1) < 1e-12);
    assert.ok(Math.min(...weights) > -1e-9, 'the point lies in its triangle');
    const interp = weights.reduce((s, w, k) => s + w * f(ids[k]), 0);
    worst = Math.max(worst, Math.abs(interp - (a[0] * v[0] + a[1] * v[1] + a[2] * v[2])));
  }
  assert.ok(worst < 2e-3, `worst error ${worst}`);
});

test('the azimuthal projection keeps distances from its centre and round-trips', () => {
  const R = 6.371e6;
  const c = lonLatToVector(40, -25);
  const proj = azimuthalProjection(c, R);
  for (const [x, y] of [[0, 0], [3e5, -1e5], [-8e5, 6e5], [1.5e6, 2e6]]) {
    const v = proj.inverse(x, y);
    assert.ok(Math.abs(Math.hypot(...v) - 1) < 1e-12);
    const ang = Math.acos(Math.min(1, v[0] * c[0] + v[1] * c[1] + v[2] * c[2]));
    assert.ok(Math.abs(R * ang - Math.hypot(x, y)) < 1e-3, 'true distance from the centre');
    const back = proj.forward(v);
    assert.ok(Math.abs(back[0] - x) < 1e-3 && Math.abs(back[1] - y) < 1e-3);
  }
  // north is -y
  const north = proj.inverse(0, -1e5);
  assert.ok(north[2] > c[2]);
});

test('surface processes run on the sphere with a sea, and the volume budget closes', () => {
  const mesh = new SphereMesh(10, 1e6);
  const sim = new Simulation({
    boundary: { open: [], seaLevel: 0 }, initial: { type: 'flat', height: 0, noise: { amplitude: 0 } },
    uplift: { rate: 1e-3 }, coast: { enabled: false }, rockfall: { enabled: false },
    strata: { strengthVariation: { amplitude: 0 } }, hillslope: { diffusivity: 50 },
    time: { dt: 1e4 },
  }, { mesh });
  // a continent: everything north of the equator starts 500 m up
  for (let i = 0; i < mesh.size; i++) sim.rock[i] = mesh.pos[3 * i + 2] > 0 ? 500 + 200 * mesh.pos[3 * i] : -2000;
  sim.setState({});
  const v0 = sim.volume();
  sim.run(20);
  for (const a of [sim.rock, sim.soil, sim.z]) assert.ok(a.every(Number.isFinite));
  const lhs = sim.volume() - v0;
  const rhs = sim.budget.uplift + sim.budget.inflow - sim.budget.exported;
  assert.ok(Math.abs(lhs - rhs) < 1e-9 * Math.max(Math.abs(sim.budget.uplift), 1), `change ${lhs} vs budget ${rhs}`);
  // rivers drain to the sea: every land cell has a downhill path
  for (let i = 0; i < mesh.size; i++) if (!sim.outlet[i]) assert.ok(sim.flow.receiver[i] >= 0);
});
