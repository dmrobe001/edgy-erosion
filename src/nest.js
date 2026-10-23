// Zooming in: refine a patch of a running simulation at finer spacing.
//
// This is coarse-graining in reverse ("backmapping", as in polymer simulation):
// the parent is frozen, a patch of it is reconstructed on a finer grid, detail the
// parent could not resolve is added, and the patch then evolves under its own
// physics with the parent supplying everything that crosses the patch edges.
//
// What the child inherits from the frozen parent:
//   - topography, soil and cumulative uplift, interpolated linearly and corrected
//     so the child's average over each parent cell equals the parent value
//     ("mean-preserving"), plus zero-mean roughness scaled to the local parent relief;
//   - the parameters, the seed and the reference frame, so strata, plugs, fault lines,
//     weak zones and rainfall line up exactly with the parent's;
//   - per-cell uplift and rainfall fields, if the parent had them (see setFields);
//   - the time and the sea level;
//   - boundary conditions: every edge cell of the patch is an outlet whose elevation
//     follows the parent's last rate of change (dz/dt), and rivers that flow into the
//     patch in the parent enter the child with their drainage area and sediment flux.
//
// The parent may be a flat hex-grid Simulation (zoomIn) or a Planet (zoomFromPlanet):
// both reduce to a MAPPING from each child cell to a triangle of parent cells with
// linear-interpolation weights, and a projection of parent cell centres into child
// coordinates. For a planet the flat patch is laid on the sphere by an azimuthal
// equidistant projection about its centre (true distances and directions from the
// centre; a 500-km patch on Earth is distorted by well under 1%).
//
// What it does NOT do: anything the child does never flows back to the parent
// (one-way coupling), and the parent's rates are frozen at its last step.

import { Simulation, referenceFrame } from './simulation.js';
import { deepMerge } from './params.js';
import { fbm } from './rng.js';
import { azimuthalProjection, lonLatToVector } from './sphere.js';

/**
 * Linear interpolation of a per-cell field of a hex-lattice simulation at (x, y):
 * barycentric on the triangle of cell centres that contains the point.
 */
export function interpolator(mesh, field) {
  const w = [0, 0, 0], ids = [0, 0, 0];
  return (x, y) => {
    if (!hexTriangle(mesh, x, y, ids, w)) return NaN;
    return w[0] * field[ids[0]] + w[1] * field[ids[1]] + w[2] * field[ids[2]];
  };
}

/** Triangle of hex-lattice cell centres containing (x, y) and barycentric weights. */
function hexTriangle(mesh, x, y, ids, w) {
  if (!mesh.lattice) throw new Error('interpolation needs a hexagonal-lattice parent mesh');
  const D = mesh.degree;
  const nb = mesh.neighbors;
  const p = mesh.cellAt(x, y);
  if (p < 0) return false;
  ids[0] = ids[1] = ids[2] = p;
  w[0] = 1; w[1] = w[2] = 0;
  const dx = x - mesh.x[p];
  const dy = y - mesh.y[p];
  if (dx === 0 && dy === 0) return true;
  // slot k points along azimuth 60k degrees (y grows southwards, so north is -y)
  let az = (Math.atan2(-dy, dx) * 180) / Math.PI;
  if (az < 0) az += 360;
  const k = Math.floor(az / 60) % 6;
  const a = nb[D * p + k];
  const b = nb[D * p + ((k + 1) % 6)];
  if (a < 0 || b < 0) return true; // at the parent's edge: nearest value
  const x1 = mesh.x[a] - mesh.x[p], y1 = mesh.y[a] - mesh.y[p];
  const x2 = mesh.x[b] - mesh.x[p], y2 = mesh.y[b] - mesh.y[p];
  const det = x1 * y2 - x2 * y1;
  const wa = (dx * y2 - x2 * dy) / det;
  const wb = (x1 * dy - dx * y1) / det;
  ids[1] = a; ids[2] = b;
  w[0] = 1 - wa - wb; w[1] = wa; w[2] = wb;
  return true;
}

/**
 * Interpolate a parent field onto every child cell through a mapping: map.k parent cells
 * per child cell (ids, weights summing to 1), and parentOf, the parent cell containing it.
 */
function sample(map, field) {
  const n = map.parentOf.length;
  const K = map.k;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let v = 0;
    for (let k = K * i; k < K * i + K; k++) v += map.w[k] * field[map.ids[k]];
    out[i] = v;
  }
  return out;
}

/**
 * Create a finer simulation of a patch of `parent`.
 * @param {Simulation} parent
 * @param {object} o
 * @param {[number, number]} o.center world position of the patch centre (m)
 * @param {number} [o.factor] refinement: child dx = parent dx / factor
 * @param {number} [o.cols] [o.rows] child grid size in cells
 * @param {number} [o.dt] child time step (default parent dt / factor)
 * @param {number} [o.roughness] amplitude of added detail, as a fraction of the local
 *   parent relief (max - min over a parent cell and its neighbours)
 * @param {number} [o.correctionIterations] passes of the mean-preserving correction
 * @returns {Simulation} the child, with `child.zoom` describing where it came from
 */
export function zoomIn(parent, { center, factor = 4, cols = 128, rows = 128, dt, roughness = 0.3, correctionIterations = 16 } = {}) {
  const pg = parent.grid;
  if (!pg.lattice) throw new Error('zooming needs a hexagonal-lattice parent');
  const dx = pg.dx / factor;
  const rowSpacing = (Math.sqrt(3) / 2) * dx;
  // centre the child grid on `center`, then keep it inside the parent's domain
  const halfW = (dx * (cols - 1 + 0.5)) / 2;
  const halfH = (rowSpacing * (rows - 1)) / 2;
  const b = pg.bounds;
  const cx = Math.min(Math.max(center[0], b.xmin + halfW + pg.dx), b.xmax - halfW - pg.dx);
  const cy = Math.min(Math.max(center[1], b.ymin + halfH + pg.dx), b.ymax - halfH - pg.dx);
  if (!(cx - halfW > b.xmin && cx + halfW < b.xmax && cy - halfH > b.ymin && cy + halfH < b.ymax)) {
    throw new Error('the requested patch does not fit inside the parent domain; use fewer cells or a larger factor');
  }
  const origin = [cx - halfW, cy - halfH];

  const pp = parent.params;
  const params = deepMerge(pp, {
    grid: { cols, rows, dx, origin },
    frame: referenceFrame(pp, pg),
    time: { dt: dt ?? pp.time.dt / factor },
    boundary: { open: ['north', 'south', 'east', 'west'], seaLevel: pp.boundary.seaLevel },
    initial: { type: 'flat', base: 0, height: 0, noise: { amplitude: 0 }, soil: 0 },
    zoom: { level: (pp.zoom?.level ?? 0) + 1, factor, center: [cx, cy], parentDx: pg.dx },
  });
  delete params.preset; // the parameters are now fully explicit
  const child = new Simulation(params);
  const cg = child.grid;
  const n = cg.size;
  const map = { k: 3, ids: new Int32Array(3 * n), w: new Float64Array(3 * n), parentOf: new Int32Array(n) };
  const ids = [0, 0, 0], w = [0, 0, 0];
  for (let i = 0; i < n; i++) {
    hexTriangle(pg, cg.x[i], cg.y[i], ids, w);
    for (let k = 0; k < 3; k++) { map.ids[3 * i + k] = ids[k]; map.w[3 * i + k] = w[k]; }
    map.parentOf[i] = pg.cellAt(cg.x[i], cg.y[i]);
  }
  refine(parentView(parent), child, map, (p) => [pg.x[p], pg.y[p]], { roughness, correctionIterations, wavelength: pg.dx, seed: (pp.seed | 0) + 4049 });
  Object.assign(child.zoom, { level: params.zoom.level, factor, center: [cx, cy], parentDx: pg.dx });
  return child;
}

/**
 * Create a flat hex-grid patch of a Planet, laid on the sphere around (lon, lat).
 * The patch inherits the planet's topography, soil, sea level, rainfall and its
 * tectonic uplift rate (isostasy included) as a per-cell uplift field; its edges follow
 * the planet's surface change rate; rivers entering from the planet bring their
 * drainage area and sediment. Erosion parameters are rescaled for the patch's spacing
 * (regionalErosionParams). From the patch, zoomIn refines further.
 * @param {import('./planet.js').Planet} planet
 * @param {object} o
 * @param {number} o.lon longitude of the patch centre (degrees)
 * @param {number} o.lat latitude of the patch centre (degrees)
 * @param {number} [o.dx] child spacing (m); default planet spacing / factor
 * @param {number} [o.factor] used when dx is not given (default 32)
 * @param {number} [o.cols] [o.rows] child grid size
 * @param {number} [o.dt] child time step (default planet dt x dx / planet spacing)
 * @param {number} [o.smooth] width of the interpolation kernel in planet cell spacings;
 *   0 = linear on the triangles of planet cells (shows their facets at the patch scale)
 * @param {object} [o.params] parameter overrides for the child (e.g. erodibilities)
 * @returns {Simulation}
 */
export function zoomFromPlanet(planet, { lon, lat, dx, factor = 32, cols = 128, rows = 128, dt, roughness = 0.3, smooth = 0.4, correctionIterations = 48, params: over = {} } = {}) {
  const mesh = planet.mesh;
  dx = dx ?? mesh.spacing / factor;
  const rowSpacing = (Math.sqrt(3) / 2) * dx;
  const halfW = (dx * (cols - 1 + 0.5)) / 2;
  const halfH = (rowSpacing * (rows - 1)) / 2;
  if (Math.max(halfW, halfH) > 0.5 * planet.params.radius) throw new Error('the requested patch is too large to lay flat on the planet; use fewer cells or a smaller dx');
  const proj = azimuthalProjection(lonLatToVector(lon, lat), planet.params.radius);
  const fac = mesh.spacing / dx;
  let params = deepMerge(planet.sim.params, {
    grid: { cols, rows, dx, origin: [-halfW, -halfH] },
    frame: null,
    time: { dt: dt ?? planet.params.time.dt / fac },
    boundary: { open: ['north', 'south', 'east', 'west'], seaLevel: planet.seaLevel },
    initial: { type: 'flat', base: 0, height: 0, noise: { amplitude: 0 }, soil: 0 },
    uplift: { rate: 0, pattern: 'uniform' },
    // waves and cliff collapse are calibrated for cells of tens of metres
    coast: { enabled: false },
    rockfall: { enabled: true },
    earthquakes: { enabled: true, faults: [] },
    zoom: { level: 1, factor: fac, center: [0, 0], parentDx: mesh.spacing, planet: { lon, lat } },
  });
  params = deepMerge(deepMerge(params, regionalErosionParams(dx, planet.params)), over);
  delete params.preset;
  const child = new Simulation(params);
  const cg = child.grid;
  const n = cg.size;
  const map = smooth > 0 ? kernelMapping(mesh, cg, proj, smooth * mesh.spacing) : linearMapping(mesh, cg, proj);
  const P3 = mesh.pos;
  const project = (p) => proj.forward([P3[3 * p], P3[3 * p + 1], P3[3 * p + 2]]);
  refine(parentView(planet), child, map, project, { roughness, correctionIterations, wavelength: mesh.spacing, seed: (planet.params.seed | 0) + 4049 });
  Object.assign(child.zoom, { level: 1, factor: fac, center: [0, 0], parentDx: mesh.spacing, planet: { lon, lat } });
  return child;
}

/** Linear interpolation on the triangles of sphere-mesh cell centres. */
function linearMapping(mesh, cg, proj) {
  const n = cg.size;
  const map = { k: 3, ids: new Int32Array(3 * n), w: new Float64Array(3 * n), parentOf: new Int32Array(n) };
  let hint = -1;
  for (let i = 0; i < n; i++) {
    const v = proj.inverse(cg.x[i], cg.y[i]);
    const loc = mesh.locate(v[0], v[1], v[2], hint);
    hint = loc.cell;
    map.parentOf[i] = loc.cell;
    for (let k = 0; k < 3; k++) { map.ids[3 * i + k] = loc.ids[k]; map.w[3 * i + k] = loc.weights[k]; }
  }
  return map;
}

/**
 * Smooth interpolation: a Gaussian-weighted average (width sigma, m) over the planet cells
 * within two rings of the nearest one. It does not pass through the parent values, but
 * the mean-preserving correction in refine restores every parent cell's average.
 */
function kernelMapping(mesh, cg, proj, sigma) {
  const n = cg.size;
  const K = 19; // a cell, its 6 neighbours and the 12 beyond (fewer next to pentagons)
  const map = { k: K, ids: new Int32Array(K * n), w: new Float64Array(K * n), parentOf: new Int32Array(n) };
  const P3 = mesh.pos;
  const R = mesh.radius;
  let hint = -1;
  const ring = [];
  for (let i = 0; i < n; i++) {
    const v = proj.inverse(cg.x[i], cg.y[i]);
    const p = mesh.nearest(v[0], v[1], v[2], hint);
    hint = p;
    map.parentOf[i] = p;
    ring.length = 0;
    ring.push(p);
    for (let r = 0, start = 0; r < 2; r++) {
      const end = ring.length;
      for (let a = start; a < end; a++) {
        for (let k = 0; k < 6; k++) {
          const q = mesh.neighbors[6 * ring[a] + k];
          if (q >= 0 && !ring.includes(q)) ring.push(q);
        }
      }
      start = end;
    }
    let sum = 0;
    for (let k = 0; k < K; k++) {
      if (k >= ring.length) { map.ids[K * i + k] = p; map.w[K * i + k] = 0; continue; }
      const q = ring[k];
      const d = R * Math.acos(Math.min(1, v[0] * P3[3 * q] + v[1] * P3[3 * q + 1] + v[2] * P3[3 * q + 2]));
      const wk = Math.exp(-(d * d) / (2 * sigma * sigma));
      map.ids[K * i + k] = q;
      map.w[K * i + k] = wk;
      sum += wk;
    }
    for (let k = 0; k < K; k++) map.w[K * i + k] /= sum;
  }
  return map;
}

/**
 * Erosion parameters for a patch of spacing dx zoomed out of a planet. Stream power and
 * creep coefficients calibrated on a grid are effective values that depend on what the
 * grid resolves: the landscape defaults suit ~30-m cells and the planet's ~150-km cells.
 * [adapted] Between them, interpolate log-linearly in log(dx). A rule of thumb, not a
 * derivation; override with `params` in zoomFromPlanet when you have better values.
 */
export function regionalErosionParams(dx, planetParams, { landscapeDx = 30, planetDx = 1.5e5 } = {}) {
  const t = Math.min(1, Math.max(0, Math.log(dx / landscapeDx) / Math.log(planetDx / landscapeDx)));
  const mix = (a, b) => Math.exp((1 - t) * Math.log(a) + t * Math.log(b));
  const s = planetParams.surface;
  const crust = s.rocks.crust;
  return {
    rocks: { crust: { ...crust, erodibility: mix(1e-5, crust.erodibility), cohesion: mix(1e6, crust.cohesion), weathering: mix(1e-4, crust.weathering) } },
    fluvial: { soilErodibility: mix(4e-4, s.fluvial.soilErodibility) },
    hillslope: { diffusivity: mix(0.01, s.hillslope.diffusivity), transportDepth: mix(0.3, s.hillslope.transportDepth) },
    weathering: { productionDepth: mix(0.5, s.weathering.productionDepth) },
  };
}

/** What refine and compareWithParent need from a parent, Simulation or Planet. */
function parentView(parent) {
  if (parent.sim) {
    const sim = parent.sim;
    return {
      mesh: parent.mesh, z: sim.z, soil: sim.soil, upliftCum: new Float64Array(parent.mesh.size),
      dzdt: parent.dzdt, flow: sim.flow, qsOut: sim.qsOut, time: parent.time,
      uplift: parent.upliftRate, precipitation: sim.precipWeight,
    };
  }
  return {
    mesh: parent.grid, z: parent.z, soil: parent.soil, upliftCum: parent.upliftCum,
    dzdt: parent.dzdt, flow: parent.flow, qsOut: parent.qsOut, time: parent.time,
    uplift: parent.fields.uplift ? parent.upliftField : null,
    precipitation: parent.fields.precipitation ? parent.precipWeight : null,
  };
}

/** Fill a freshly made child from its parent through a mapping (see the header). */
function refine(pv, child, map, project, { roughness, correctionIterations, wavelength, seed }) {
  const pg = pv.mesh;
  const cg = child.grid;
  const n = cg.size;
  const pz = pv.z;

  // --- topography: linear interpolation, roughness, mean-preserving correction ---
  const relief = new Float64Array(pg.size);
  for (let p = 0; p < pg.size; p++) {
    let lo = pz[p], hi = pz[p];
    for (let k = 0; k < pg.degree; k++) {
      const q = pg.neighbors[pg.degree * p + k];
      if (q >= 0) { lo = Math.min(lo, pz[q]); hi = Math.max(hi, pz[q]); }
    }
    relief[p] = hi - lo;
  }
  const z = sample(map, pz);
  const reliefC = sample(map, relief);
  const soil = sample(map, pv.soil);
  const upliftCum = sample(map, pv.upliftCum);
  const rate = sample(map, pv.dzdt);
  const boundaryRate = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const detail = fbm(cg.x[i], cg.y[i], { wavelength, octaves: 4, seed });
    z[i] += roughness * 0.5 * reliefC[i] * detail;
    if (soil[i] < 0) soil[i] = 0;
    if (child.boundaryOutlet[i]) boundaryRate[i] = rate[i];
  }
  // Correct so each parent cell's child-average equals the parent value: interpolate the
  // per-parent-cell misfit and add it, a few times (a fixed-point iteration).
  const parentOf = map.parentOf;
  const sum = new Float64Array(pg.size);
  const cnt = new Float64Array(pg.size);
  const misfit = new Float64Array(pg.size);
  const childArea = cg.area[0];
  let residual = 0;
  for (let it = 0; it <= correctionIterations; it++) {
    sum.fill(0); cnt.fill(0);
    for (let i = 0; i < n; i++) { sum[parentOf[i]] += z[i]; cnt[parentOf[i]]++; }
    residual = 0;
    for (let p = 0; p < pg.size; p++) {
      misfit[p] = cnt[p] ? pz[p] - sum[p] / cnt[p] : 0;
      // judge convergence on parent cells the patch covers completely
      if (cnt[p] * childArea >= 0.9 * pg.area[p]) residual = Math.max(residual, Math.abs(misfit[p]));
    }
    if (it === correctionIterations) break;
    const corr = sample(map, misfit);
    for (let i = 0; i < n; i++) z[i] += corr[i];
  }

  // --- rivers entering the patch from the parent ---
  const inflowArea = new Float64Array(n);
  const inflowSediment = new Float64Array(n);
  const where = new Map();
  const childCell = (p) => {
    if (where.has(p)) return where.get(p);
    const xy = project(p);
    const c = xy ? cg.cellAt(xy[0], xy[1]) : -1;
    const r = c >= 0 && !child.boundaryOutlet[c] ? c : -1;
    where.set(p, r);
    return r;
  };
  const pf = pv.flow;
  const pD = pg.degree;
  const mfd = pf.flowExponent > 0;
  let enteringArea = 0;
  for (let p = 0; p < pg.size; p++) {
    if (pf.receiver[p] < 0 || childCell(p) >= 0) continue;
    // sediment follows the single receiver
    const r = pf.receiver[p];
    const tr = childCell(r);
    if (tr >= 0) inflowSediment[tr] += pv.qsOut[p];
    // water: all to the receiver, or shared by slope^p as in the parent's MFD routing
    if (!mfd) {
      if (tr >= 0) { inflowArea[tr] += pf.area[p]; enteringArea += pf.area[p]; }
      continue;
    }
    let wsum = 0;
    const w = [];
    for (let k = 0; k < pD; k++) {
      const e = pD * p + k;
      const q = pg.neighbors[e];
      const d = q >= 0 ? pf.filled[p] - pf.filled[q] : 0;
      const wk = d > 0 ? Math.pow(d / pg.dist[e], pf.flowExponent) : 0;
      w.push([q, wk]);
      wsum += wk;
    }
    for (const [q, wk] of w) {
      if (wk <= 0) continue;
      const t = childCell(q);
      if (t >= 0) { const a = (pf.area[p] * wk) / wsum; inflowArea[t] += a; enteringArea += a; }
    }
  }

  // per-cell fields the parent had (a planet's uplift and rainfall) carry over
  const fields = {};
  if (pv.uplift) fields.uplift = sample(map, pv.uplift);
  if (pv.precipitation) fields.precipitation = sample(map, pv.precipitation);
  if (fields.uplift || fields.precipitation) child.setFields(fields);

  // open edges carry no soil; their elevation is the surface itself
  const rock = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (child.boundaryOutlet[i]) soil[i] = 0;
    rock[i] = z[i] - soil[i];
  }
  child.setState({ rock, soil, upliftCum, time: pv.time, inflowArea, inflowSediment, boundaryRate });
  child.stepCount = 0;
  child.zoom = { bounds: { ...cg.bounds }, meanResidual: residual, enteringArea, parentTime: pv.time, parentOf, project };
}

/**
 * Compare a zoomed child with its frozen parent (a Simulation or a Planet) over the
 * parent cells the patch covers completely: mean elevation, mean rate of change, and the
 * discharge of the largest river leaving the patch. These are the cross-scale consistency
 * checks. The parent is frozen, so its "expected" elevation is extrapolated at its last
 * rate of change over the time the child has run.
 */
export function compareWithParent(child, parent) {
  const pv = parentView(parent);
  const pg = pv.mesh;
  const cg = child.grid;
  const parentOf = child.zoom.parentOf;
  const sumZ = new Float64Array(pg.size);
  const sumR = new Float64Array(pg.size);
  const cnt = new Float64Array(pg.size);
  for (let i = 0; i < cg.size; i++) {
    if (child.boundaryOutlet[i]) continue;
    const p = parentOf[i];
    sumZ[p] += child.z[i]; sumR[p] += child.dzdt[i]; cnt[p]++;
  }
  const childArea = cg.area[0];
  let nP = 0, parentZ = 0, childZ = 0, parentRate = 0, childRate = 0;
  for (let p = 0; p < pg.size; p++) {
    if (cnt[p] * childArea < 0.9 * pg.area[p]) continue;
    nP++;
    parentZ += pv.z[p]; childZ += sumZ[p] / cnt[p];
    parentRate += pv.dzdt[p]; childRate += sumR[p] / cnt[p];
  }
  // biggest river leaving: child = largest drainage area arriving at an edge cell;
  // parent = largest area at a covered cell whose receiver lies outside the patch
  let childExit = 0;
  for (let i = 0; i < cg.size; i++) if (child.boundaryOutlet[i]) childExit = Math.max(childExit, child.flow.area[i]);
  let parentExit = 0;
  const b = cg.bounds;
  const within = (p) => {
    const xy = child.zoom.project(p);
    return !!xy && xy[0] > b.xmin && xy[0] < b.xmax && xy[1] > b.ymin && xy[1] < b.ymax;
  };
  for (let p = 0; p < pg.size; p++) {
    const r = pv.flow.receiver[p];
    if (r < 0 || !within(p) || within(r)) continue;
    parentExit = Math.max(parentExit, pv.flow.area[p]);
  }
  const elapsed = child.time - child.zoom.parentTime;
  const m = Math.max(nP, 1);
  return {
    parentCells: nP,
    elapsed,
    meanElevation: { parent: parentZ / m, expected: (parentZ + parentRate * elapsed) / m, child: childZ / m },
    meanRate: { parent: parentRate / m, child: childRate / m },
    largestExitArea: { parent: parentExit, child: childExit },
  };
}
