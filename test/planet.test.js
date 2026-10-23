import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Planet } from '../src/planet.js';
import { zoomFromPlanet, compareWithParent } from '../src/nest.js';
import { faultMeanUplift } from '../src/processes/earthquakes.js';

const small = { subdivisions: 16, plates: { count: 6 } };

test('isostasy: reference continent, ridge depth, sqrt-age subsidence, and erosional rebound', () => {
  const pl = new Planet(small);
  const cr = pl.params.crust;
  const near = (a, b, tol, what) => assert.ok(Math.abs(a - b) < tol, `${what}: ${a} vs ${b}`);
  near(pl.elevation(cr.continentalThickness, 1, 0), cr.referenceElevation, 1e-9, 'reference continent');
  near(pl.elevation(cr.oceanicThickness, 0, 0), -cr.ridgeDepth, 1e-9, 'ridge');
  near(pl.elevation(cr.oceanicThickness, 0, 25e6), -cr.ridgeDepth - cr.subsidence * 5, 1e-9, '25 Myr sea floor');
  // old floor stops sinking once it reaches the depth of cold crust
  near(pl.elevation(cr.oceanicThickness, 0, 300e6), pl.elevation(cr.oceanicThickness, 0, 200e6), 1e-9, 'old floor');
  // removing 1 km from a continent lowers it by only (1 - rho_c / rho_m) km
  const drop = pl.elevation(40000, 1, 0) - pl.elevation(39000, 1, 0);
  near(drop, 1000 * (1 - cr.continentalDensity / cr.mantleDensity), 1e-9, 'rebound');
});

test('crust volume is conserved: change = ridges + arcs - subduction + surface, to round-off', () => {
  const pl = new Planet(small);
  pl.run(25);
  const b = pl.budget;
  const lhs = pl.crustVolume() - pl.volume0;
  const rhs = b.created + b.arc - b.subducted + b.surface;
  assert.ok(b.created > 0 && b.subducted > 0 && b.arc > 0, 'ridges spread and slabs sink');
  assert.ok(Math.abs(lhs - rhs) < 1e-10 * pl.volume0, `change ${lhs} vs budget ${rhs}`);
  assert.ok(Math.abs(b.surface) < 1e-9 * pl.volume0, 'rivers keep their sediment on the planet (deltas)');
  // and the ocean's water volume is held fixed by moving sea level
  assert.ok(Math.abs(pl._waterBelow(pl.seaLevel) / pl.waterVolume - 1) < 1e-9);
});

test('a single plate rotates rigidly: no crust made or destroyed, the continent keeps its shape', () => {
  const pl = new Planet({ ...small, plates: { count: 1 }, surface: { fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false } }, crust: { flowDiffusivity: 0 } });
  const w = 0.05 / pl.params.radius; // 5 cm/yr at the equator, about the z axis
  pl.setPlates(() => 0, [[0, 0, w]]);
  const n = pl.mesh.size;
  const contBefore = Array.from(pl.continental);
  const P3 = pl.mesh.pos;
  const centroid = (c) => {
    let x = 0, y = 0, z = 0;
    for (let i = 0; i < n; i++) if (c[i] > 0.5) { x += P3[3 * i] * pl.mesh.area[i]; y += P3[3 * i + 1] * pl.mesh.area[i]; z += P3[3 * i + 2] * pl.mesh.area[i]; }
    return Math.atan2(y, x);
  };
  const lon0 = centroid(contBefore);
  const area = (c) => c.reduce((s, v, i) => s + (v > 0.5 ? pl.mesh.area[i] : 0), 0);
  const steps = 40;
  pl.run(steps);
  assert.equal(pl.budget.created, 0);
  assert.equal(pl.budget.subducted, 0);
  // rotated by w t radians about z
  let dlon = centroid(pl.continental) - lon0;
  dlon = Math.atan2(Math.sin(dlon), Math.cos(dlon));
  const expected = w * steps * pl.params.time.dt;
  assert.ok(Math.abs(dlon - expected) < 0.25 * pl.mesh.spacing / pl.params.radius, `rotated ${dlon} rad, expected ${expected}`);
  assert.ok(Math.abs(area(pl.continental) / area(contBefore) - 1) < 0.03, 'continental area kept');
});

test('two plates spreading apart: ocean floor ages linearly away from the ridge', () => {
  const pl = new Planet({ ...small, subdivisions: 40, plates: { count: 2 }, continents: { fraction: 0 }, surface: { fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false } } });
  const R = pl.params.radius;
  const w = 0.04 / R; // each plate moves 4 cm/yr at its fastest
  // plate 0 is the x > 0 hemisphere. Rotating the plates in opposite senses about z
  // pulls them apart along the boundary where y > 0 and pushes them together where y < 0.
  pl.setPlates((i) => (pl.mesh.pos[3 * i] > 0 ? 0 : 1), [[0, 0, -w], [0, 0, w]]);
  const T = 40e6;
  pl.run(T / pl.params.time.dt);
  // along the equator near longitude 90 degrees the ridge is the meridian x = 0 and each
  // side moves away from it at w R: age = distance / (w R). Only the floor made since the
  // plates started moving counts (within w R T of the ridge, less a cell of slack).
  const pts = [];
  const reach = w * R * T - pl.mesh.spacing;
  for (let i = 0; i < pl.mesh.size; i++) {
    const lat = pl.mesh.lat[i], lon = pl.mesh.lon[i];
    const d = Math.abs(lon - Math.PI / 2) * R;
    // (cells still showing pre-existing floor are not new crust: skip them)
    if (Math.abs(lat) > 0.1 || d > reach || pl.age[i] > T) continue;
    pts.push([d, pl.age[i]]);
  }
  // least-squares slope of age against distance from the ridge
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (const [x, y] of pts) { sx += x; sy += y; sxx += x * x; sxy += x * y; }
  const N = pts.length;
  const slope = (N * sxy - sx * sy) / (N * sxx - sx * sx);
  const halfRate = 1 / slope; // m/yr
  assert.ok(N > 20);
  assert.ok(Math.abs(halfRate / (w * R) - 1) < 0.1, `half-spreading rate ${(halfRate * 100).toFixed(2)} cm/yr, expected ${(w * R * 100).toFixed(2)}`);
  // and the closing side has a subduction zone with a trench
  let trench = 0;
  for (let i = 0; i < pl.mesh.size; i++) if (pl.mesh.pos[3 * i + 1] < -0.5 && pl.boundary[i] === 1 && pl.sim.z[i] < -7000) trench++;
  assert.ok(trench > 0, 'a trench where the plates converge');
});

test('planets are deterministic', () => {
  const a = new Planet(small).run(5);
  const b = new Planet(small).run(5);
  assert.deepEqual(Array.from(a.sim.z), Array.from(b.sim.z));
  const c = new Planet({ ...small, seed: 2 }).run(5);
  assert.notDeepEqual(Array.from(a.sim.z), Array.from(c.sim.z));
});

test('a patch zoomed out of a planet keeps its block means, sea level, uplift and rainfall', () => {
  const pl = new Planet({ subdivisions: 20 });
  pl.run(10);
  // a continental cell with some uplift
  let c = -1;
  for (let i = 0; i < pl.mesh.size; i++) if (pl.sim.z[i] > pl.seaLevel + 200 && (c < 0 || pl.upliftRate[i] > pl.upliftRate[c])) c = i;
  const lon = (pl.mesh.lon[c] * 180) / Math.PI, lat = (pl.mesh.lat[c] * 180) / Math.PI;
  const child = zoomFromPlanet(pl, { lon, lat, factor: 24, cols: 96, rows: 96 });
  assert.ok(child.zoom.meanResidual < 0.5, `block-mean residual ${child.zoom.meanResidual} m`);
  const cmp = compareWithParent(child, pl);
  assert.ok(cmp.parentCells >= 2);
  assert.ok(Math.abs(cmp.meanElevation.child - cmp.meanElevation.parent) < 1, 'means agree at the hand-off');
  assert.equal(child.seaLevel, pl.seaLevel);
  // the centre cell of the patch gets the planet's uplift and rainfall (smoothed)
  const mid = child.grid.index(48, 48);
  assert.ok(child.fields.uplift && child.fields.precipitation);
  assert.ok(Math.abs(child.upliftField[mid] - pl.upliftRate[c]) < 0.5 * Math.abs(pl.upliftRate[c]) + 2e-5);
  assert.ok(Math.abs(child.precipWeight[mid] / pl.sim.precipWeight[c] - 1) < 0.2);
  child.run(5);
  for (const a of [child.rock, child.soil, child.z]) assert.ok(a.every(Number.isFinite));
});

test('a plate boundary crossing a patch becomes a fault with the plates\' convergence, and long-term uplift still matches the planet', () => {
  // two plates converging along the meridian x = 0 where y < 0 (see the spreading test)
  const pl = new Planet({ subdivisions: 32, plates: { count: 2 }, continents: { fraction: 0 } });
  const R = pl.params.radius;
  const w = 0.02 / R;
  pl.setPlates((i) => (pl.mesh.pos[3 * i] > 0 ? 0 : 1), [[0, 0, -w], [0, 0, w]]);
  pl.run(30); // long enough for subduction to start and pick the overriding plate
  // erosion off so uplift can be measured; fewer, larger events keep the test quick
  const quiet = { fluvial: { enabled: false }, hillslope: { enabled: false }, weathering: { enabled: false }, rockfall: { enabled: false }, earthquakes: { minMagnitude: 6.5 } };
  const child = zoomFromPlanet(pl, { lon: -90, lat: 0, factor: 24, cols: 80, rows: 80, params: quiet });
  const thrusts = child.zoom.boundaries.filter((b) => b.kind === 'thrust');
  assert.ok(thrusts.length >= 1, JSON.stringify(child.zoom.boundaries.map((b) => b.kind)));
  // at longitude -90 the plates close head-on at 2 w R = 4 cm/yr
  for (const t of thrusts) {
    assert.ok(Math.abs(t.closingRate / (2 * w * R) - 1) < 0.1, `closing ${t.closingRate}`);
    assert.ok(Math.abs(t.fault.slipRate - (0.1 * t.closingRate) / Math.cos(Math.PI / 6)) < 1e-12);
    assert.equal(t.raised, pl.polarity.get(0 * 2 + 1), 'the overriding plate is raised');
    // the trace runs north-south, like the boundary (the meridian through the patch centre)
    const F = child.frame;
    const dxw = (t.fault.x1 - t.fault.x0) * F.width, dyw = (t.fault.y1 - t.fault.y0) * F.height;
    assert.ok(Math.abs(dyw) > 5 * Math.abs(dxw), 'north-south trace');
  }
  // the smooth field plus what the faults deliver on average is the planet's uplift, exactly
  const total = Float64Array.from(child.upliftField);
  child.params.earthquakes.faults.forEach((f, k) => { const m = faultMeanUplift(child, k); for (let i = 0; i < total.length; i++) total[i] += m[i]; });
  for (let i = 0; i < total.length; i++) {
    if (child.boundaryOutlet[i]) continue;
    assert.ok(Math.abs(total[i] - child.zoom.planetUplift[i]) < 1e-12, `cell ${i}`);
  }
  // and a long run delivers it: mean change near the fault over 2 Myr against the planet's rate
  const z0 = Float64Array.from(child.z);
  const steps = 100;
  child.run(steps);
  const T = steps * child.params.time.dt;
  let got = 0, want = 0, cnt = 0;
  const g = child.grid;
  for (let i = 0; i < g.size; i++) {
    if (child.boundaryOutlet[i] || g.col(i) < 10 || g.col(i) > 69 || g.row(i) < 10 || g.row(i) > 69) continue;
    got += child.z[i] - z0[i]; want += child.zoom.planetUplift[i] * T; cnt++;
  }
  assert.ok(child.quakes.count > 100, `${child.quakes.count} earthquakes`);
  assert.ok(Math.abs(got - want) < 0.1 * Math.abs(want) + cnt * 1, `mean change ${(got / cnt).toFixed(1)} m vs planet ${(want / cnt).toFixed(1)} m`);
});
