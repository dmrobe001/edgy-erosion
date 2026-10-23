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
//   - the time and the sea level;
//   - boundary conditions: every edge cell of the patch is an outlet whose elevation
//     follows the parent's last rate of change (dz/dt), and rivers that flow into the
//     patch in the parent enter the child with their drainage area and sediment flux.
//
// What it does NOT do: anything the child does never flows back to the parent
// (one-way coupling), and the parent's rates are frozen at its last step.

import { Simulation, referenceFrame } from './simulation.js';
import { deepMerge } from './params.js';
import { fbm } from './rng.js';

/**
 * Linear interpolation of a per-cell field of a hex-lattice simulation at (x, y):
 * barycentric on the triangle of cell centres that contains the point.
 */
export function interpolator(mesh, field) {
  if (!mesh.lattice) throw new Error('interpolation needs a hexagonal-lattice parent mesh');
  const D = mesh.degree;
  const nb = mesh.neighbors;
  return (x, y) => {
    const p = mesh.cellAt(x, y);
    if (p < 0) return NaN;
    const dx = x - mesh.x[p];
    const dy = y - mesh.y[p];
    if (dx === 0 && dy === 0) return field[p];
    // slot k points along azimuth 60k degrees (y grows southwards, so north is -y)
    let az = (Math.atan2(-dy, dx) * 180) / Math.PI;
    if (az < 0) az += 360;
    const k = Math.floor(az / 60) % 6;
    const a = nb[D * p + k];
    const b = nb[D * p + ((k + 1) % 6)];
    if (a < 0 || b < 0) return field[p]; // at the parent's edge: nearest value
    // barycentric coordinates of (x, y) in triangle (p, a, b)
    const x1 = mesh.x[a] - mesh.x[p], y1 = mesh.y[a] - mesh.y[p];
    const x2 = mesh.x[b] - mesh.x[p], y2 = mesh.y[b] - mesh.y[p];
    const det = x1 * y2 - x2 * y1;
    const wa = (dx * y2 - x2 * dy) / det;
    const wb = (x1 * dy - dx * y1) / det;
    return (1 - wa - wb) * field[p] + wa * field[a] + wb * field[b];
  };
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

  // --- topography: linear interpolation, roughness, mean-preserving correction ---
  const pz = parent.z;
  const zAt = interpolator(pg, pz);
  const soilAt = interpolator(pg, parent.soil);
  const upAt = interpolator(pg, parent.upliftCum);
  const rateAt = interpolator(pg, parent.dzdt);
  const relief = new Float64Array(pg.size);
  for (let p = 0; p < pg.size; p++) {
    let lo = pz[p], hi = pz[p];
    for (let k = 0; k < pg.degree; k++) {
      const q = pg.neighbors[pg.degree * p + k];
      if (q >= 0) { lo = Math.min(lo, pz[q]); hi = Math.max(hi, pz[q]); }
    }
    relief[p] = hi - lo;
  }
  const reliefAt = interpolator(pg, relief);
  const parentOf = new Int32Array(n);
  const z = new Float64Array(n);
  const soil = new Float64Array(n);
  const upliftCum = new Float64Array(n);
  const boundaryRate = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const x = cg.x[i], y = cg.y[i];
    parentOf[i] = pg.cellAt(x, y);
    const detail = fbm(x, y, { wavelength: pg.dx, octaves: 4, seed: (pp.seed | 0) + 4049 });
    z[i] = zAt(x, y) + roughness * 0.5 * reliefAt(x, y) * detail;
    soil[i] = Math.max(0, soilAt(x, y));
    upliftCum[i] = upAt(x, y);
    if (child.boundaryOutlet[i]) boundaryRate[i] = rateAt(x, y);
  }
  // Correct so each parent cell's child-average equals the parent value: interpolate the
  // per-parent-cell misfit and add it, a few times (a fixed-point iteration).
  const sum = new Float64Array(pg.size);
  const cnt = new Float64Array(pg.size);
  const misfit = new Float64Array(pg.size);
  const misfitAt = interpolator(pg, misfit);
  let residual = 0;
  for (let it = 0; it <= correctionIterations; it++) {
    sum.fill(0); cnt.fill(0);
    for (let i = 0; i < n; i++) { sum[parentOf[i]] += z[i]; cnt[parentOf[i]]++; }
    residual = 0;
    for (let p = 0; p < pg.size; p++) {
      misfit[p] = cnt[p] ? pz[p] - sum[p] / cnt[p] : 0;
      // judge convergence on parent cells the patch covers completely
      if (cnt[p] >= 0.9 * factor * factor) residual = Math.max(residual, Math.abs(misfit[p]));
    }
    if (it === correctionIterations) break;
    for (let i = 0; i < n; i++) z[i] += misfitAt(cg.x[i], cg.y[i]);
  }

  // --- rivers entering the patch from the parent ---
  const inflowArea = new Float64Array(n);
  const inflowSediment = new Float64Array(n);
  const inside = (x, y) => { const c = cg.cellAt(x, y); return c >= 0 && !child.boundaryOutlet[c] ? c : -1; };
  const pf = parent.flow;
  const pD = pg.degree;
  const mfd = pf.flowExponent > 0;
  let enteringArea = 0;
  for (let p = 0; p < pg.size; p++) {
    if (pf.receiver[p] < 0 || inside(pg.x[p], pg.y[p]) >= 0) continue;
    // sediment follows the single receiver
    const r = pf.receiver[p];
    const tr = inside(pg.x[r], pg.y[r]);
    if (tr >= 0) inflowSediment[tr] += parent.qsOut[p];
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
      const t = inside(pg.x[q], pg.y[q]);
      if (t >= 0) { const a = (pf.area[p] * wk) / wsum; inflowArea[t] += a; enteringArea += a; }
    }
  }

  // open edges carry no soil; their elevation is the surface itself
  const rock = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (child.boundaryOutlet[i]) soil[i] = 0;
    rock[i] = z[i] - soil[i];
  }
  child.setState({ rock, soil, upliftCum, time: parent.time, inflowArea, inflowSediment, boundaryRate });
  child.stepCount = 0;
  child.zoom = {
    level: params.zoom.level, factor, center: [cx, cy], parentDx: pg.dx,
    bounds: { ...cg.bounds }, meanResidual: residual, enteringArea, parentTime: parent.time,
  };
  return child;
}

/**
 * Compare a zoomed child with its frozen parent over the parent cells the patch covers
 * completely: mean elevation, mean rate of change, and the discharge of the largest
 * river leaving the patch. These are the cross-scale consistency checks. The parent is
 * frozen, so its "expected" elevation is extrapolated at its last rate of change over
 * the time the child has run.
 */
export function compareWithParent(child, parent) {
  const pg = parent.grid;
  const cg = child.grid;
  const factor = child.zoom.factor;
  const sumZ = new Float64Array(pg.size);
  const sumR = new Float64Array(pg.size);
  const cnt = new Float64Array(pg.size);
  for (let i = 0; i < cg.size; i++) {
    if (child.boundaryOutlet[i]) continue;
    const p = pg.cellAt(cg.x[i], cg.y[i]);
    sumZ[p] += child.z[i]; sumR[p] += child.dzdt[i]; cnt[p]++;
  }
  let nP = 0, parentZ = 0, childZ = 0, parentRate = 0, childRate = 0;
  for (let p = 0; p < pg.size; p++) {
    if (cnt[p] < 0.9 * factor * factor) continue;
    nP++;
    parentZ += parent.z[p]; childZ += sumZ[p] / cnt[p];
    parentRate += parent.dzdt[p]; childRate += sumR[p] / cnt[p];
  }
  // biggest river leaving: child = largest drainage area arriving at an edge cell;
  // parent = largest area at a covered cell whose receiver lies outside the patch
  let childExit = 0;
  for (let i = 0; i < cg.size; i++) if (child.boundaryOutlet[i]) childExit = Math.max(childExit, child.flow.area[i]);
  let parentExit = 0;
  const b = cg.bounds;
  const within = (x, y) => x > b.xmin && x < b.xmax && y > b.ymin && y < b.ymax;
  for (let p = 0; p < pg.size; p++) {
    const r = parent.flow.receiver[p];
    if (r < 0 || !within(pg.x[p], pg.y[p]) || within(pg.x[r], pg.y[r])) continue;
    parentExit = Math.max(parentExit, parent.flow.area[p]);
  }
  const elapsed = child.time - child.zoom.parentTime;
  return {
    parentCells: nP,
    elapsed,
    meanElevation: { parent: parentZ / nP, expected: (parentZ + parentRate * elapsed) / nP, child: childZ / nP },
    meanRate: { parent: parentRate / nP, child: childRate / nP },
    largestExitArea: { parent: parentExit, child: childExit },
  };
}
