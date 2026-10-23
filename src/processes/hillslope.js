// Hillslopes: soil creep and soil/talus avalanching.
//
// CREEP. Soil moves downslope with a flux per unit width
//
//   q = D * S / (1 - (S/Sc)^2) * (1 - exp(-h / h0))
//
// the nonlinear creep law of Roering, Kirchner & Dietrich 1999 (Water Resour. Res. 35)
// with a depth dependence (thin soil creeps less; bare rock not at all). Only soil
// moves, so outflow from a cell is capped at the soil it has. The enhancement factor
// 1/(1-(S/Sc)^2) is capped at `maxNonlinearity`; anything steeper than Sc is handled
// by avalanching instead.
//
// On the hex lattice, the flux across the edge shared by cells i and j changes the
// height of each by
//     dh = q * edgeLength / cellArea * dt = q * 2 / (3 dx) * dt,   with S = dz / dx,
// which works out to D * f * dz * 2 / (3 dx^2) * dt -- the standard hex-lattice
// Laplacian. Forward Euler is stable for dt < dx^2 / (3 D f); we substep.
//
// AVALANCHING. Soil steeper than the angle of repose (slope Sc) slides to lower
// neighbours until it is not, conserving volume (a "thermal erosion" relaxation as
// in Musgrave, Kolb & Mace 1989, SIGGRAPH). This is what builds talus slopes under
// cliffs from rockfall debris.

const LAPLACE = 2 / 3; // hex-lattice Laplacian coefficient (divided by dx^2)

/** Effective critical slope for cell i, raised by vegetation (root cohesion). */
function criticalSlope(sim, i) {
  return sim.params.hillslope.criticalSlope * (1 + sim.params.vegetation.slopeBoost * sim.veg[i]);
}

/**
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} dt years
 * @returns {{substeps:number}}
 */
export function creep(sim, dt) {
  const hp = sim.params.hillslope;
  const vp = sim.params.vegetation;
  const { grid, soil, z, veg } = sim;
  const n = grid.size;
  const nb = grid.neighbors;
  const dx = grid.dx;
  const D0 = hp.diffusivity;
  const fmax = Math.max(1, hp.maxNonlinearity);
  const h0 = hp.transportDepth;
  if (!(D0 > 0)) return { substeps: 0 };

  let vmax = 0;
  for (let i = 0; i < n; i++) if (veg[i] > vmax) vmax = veg[i];
  const Dmax = D0 * (1 + vp.diffusivityBoost * vmax);
  const dtStable = (0.5 * dx * dx) / (3 * Dmax * fmax);
  const substeps = Math.max(1, Math.ceil(dt / dtStable));
  const h = dt / substeps;
  const coef = (LAPLACE * h) / (dx * dx);

  const flux = sim.scratch.edgeFlux; // 3 per cell: directions 0 (E), 1 (NE), 2 (NW)
  const outTotal = sim.scratch.outTotal;
  for (let sub = 0; sub < substeps; sub++) {
    outTotal.fill(0);
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < 3; k++) {
        const j = nb[6 * i + k];
        const e = 3 * i + k;
        if (j < 0) { flux[e] = 0; continue; }
        const dz = z[i] - z[j];
        const up = dz > 0 ? i : j;
        const S = Math.abs(dz) / dx;
        const Sc = criticalSlope(sim, up);
        const r = S / Sc;
        const f = r < 1 ? Math.min(fmax, 1 / (1 - r * r)) : fmax;
        const D = D0 * (1 + vp.diffusivityBoost * veg[up]);
        const depth = 1 - Math.exp(-soil[up] / h0);
        const q = D * f * depth * dz * coef; // height moved i -> j (negative: j -> i)
        flux[e] = q;
        outTotal[up] += Math.abs(q);
      }
    }
    // scale factor so no cell exports more soil than it has
    for (let i = 0; i < n; i++) {
      outTotal[i] = outTotal[i] > soil[i] ? soil[i] / outTotal[i] : 1;
    }
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < 3; k++) {
        const q0 = flux[3 * i + k];
        if (q0 === 0) continue;
        const j = nb[6 * i + k];
        const q = q0 > 0 ? q0 * outTotal[i] : q0 * outTotal[j];
        soil[i] -= q;
        soil[j] += q;
        z[i] -= q;
        z[j] += q;
      }
    }
    for (let i = 0; i < n; i++) if (soil[i] < 0) { z[i] -= soil[i]; soil[i] = 0; } // round-off only
  }
  return { substeps };
}

/**
 * Relax soil slopes steeper than the critical slope.
 * The first sweep visits every cell; later sweeps only revisit cells next to
 * where something moved (a worklist), alternating direction to avoid bias.
 * @returns {number} total volume moved (m^3)
 */
export function avalanche(sim) {
  const { grid, soil, z, boundaryOutlet } = sim;
  const n = grid.size;
  const nb = grid.neighbors;
  const dx = grid.dx;
  const maxIter = sim.params.hillslope.avalancheIterations;
  const excess = new Float64Array(6);
  const ws = sim.scratch.avalanche || (sim.scratch.avalanche = {
    list: new Int32Array(n), next: new Int32Array(n), mark: new Uint32Array(n), stamp: 0,
  });
  let list = ws.list;
  let next = ws.next;
  let count = n;
  for (let i = 0; i < n; i++) list[i] = i;
  let moved = 0;
  for (let it = 0; it < maxIter && count > 0; it++) {
    const stamp = ++ws.stamp;
    let nextCount = 0;
    const forward = (it & 1) === 0;
    for (let s = 0; s < count; s++) {
      const i = list[forward ? s : count - 1 - s];
      if (soil[i] <= 0 || boundaryOutlet[i]) continue;
      const lim = criticalSlope(sim, i) * dx;
      let sum = 0;
      let emax = 0;
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        let e = 0;
        if (j >= 0) {
          e = z[i] - z[j] - lim;
          if (e < 0) e = 0;
        }
        excess[k] = e;
        sum += e;
        if (e > emax) emax = e;
      }
      if (emax <= 1e-6) continue;
      // Move half the largest excess (the amount that would level the steepest pair),
      // shared among all over-steep neighbours in proportion to their excess.
      const M = Math.min(soil[i], 0.5 * emax);
      for (let k = 0; k < 6; k++) {
        if (excess[k] === 0) continue;
        const j = nb[6 * i + k];
        const dm = (M * excess[k]) / sum;
        soil[j] += dm;
        z[j] += dm;
      }
      soil[i] -= M;
      z[i] -= M;
      moved += M;
      // i and its neighbours may now be unstable
      if (ws.mark[i] !== stamp) { ws.mark[i] = stamp; next[nextCount++] = i; }
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j >= 0 && ws.mark[j] !== stamp) { ws.mark[j] = stamp; next[nextCount++] = j; }
      }
    }
    [list, next] = [next, list];
    count = nextCount;
  }
  return moved * grid.cellArea;
}
