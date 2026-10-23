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
// Finite volumes (see mesh.js): the flux across the face shared by cells i and j moves
// volume q * face * dt, with S = dz / dist. On the regular hex grid this is
// D * f * dz * 2 / (3 dx^2) * dt per cell, the standard isotropic hex Laplacian.
// Explicit substeps keep forward Euler stable.
//
// AVALANCHING. Soil steeper than the angle of repose (slope Sc) slides to lower
// neighbours until it is not, conserving volume (a "thermal erosion" relaxation as
// in Musgrave, Kolb & Mace 1989, SIGGRAPH). This is what builds talus slopes under
// cliffs from rockfall debris.

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
  const { grid: mesh, soil, z, veg } = sim;
  const n = mesh.size;
  const Dg = mesh.degree;
  const nb = mesh.neighbors;
  const { dist, face, area } = mesh;
  const D0 = hp.diffusivity;
  const fmax = Math.max(1, hp.maxNonlinearity);
  const h0 = hp.transportDepth;
  if (!(D0 > 0)) return { substeps: 0 };

  let vmax = 0;
  for (let i = 0; i < n; i++) if (veg[i] > vmax) vmax = veg[i];
  const Dmax = D0 * (1 + vp.diffusivityBoost * vmax) * fmax;
  // Forward Euler is stable when dt * sum_j(D face / (dist area_i)) <= 1 for every cell
  // (Gershgorin); keep a safety factor.
  let worst = 0;
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let k = 0; k < Dg; k++) { const e = Dg * i + k; if (nb[e] >= 0) sum += face[e] / dist[e]; }
    worst = Math.max(worst, sum / area[i]);
  }
  const dtStable = 0.75 / (Dmax * worst);
  const substeps = Math.max(1, Math.ceil(dt / dtStable));
  const h = dt / substeps;

  const flux = sim.scratch.edgeFlux; // volume moved per edge per substep, indexed like neighbors
  const outTotal = sim.scratch.outTotal;
  for (let sub = 0; sub < substeps; sub++) {
    outTotal.fill(0);
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < Dg; k++) {
        const e = Dg * i + k;
        const j = nb[e];
        flux[e] = 0;
        if (j <= i) continue; // each edge once, from its lower-numbered cell
        const dz = z[i] - z[j];
        const up = dz > 0 ? i : j;
        const S = Math.abs(dz) / dist[e];
        const r = S / criticalSlope(sim, up);
        const f = r < 1 ? Math.min(fmax, 1 / (1 - r * r)) : fmax;
        const D = D0 * (1 + vp.diffusivityBoost * veg[up]);
        const depth = 1 - Math.exp(-soil[up] / h0);
        const q = D * f * depth * (dz / dist[e]) * face[e] * h; // volume moved i -> j (negative: j -> i)
        flux[e] = q;
        outTotal[up] += Math.abs(q);
      }
    }
    // scale factor so no cell exports more soil than it has
    for (let i = 0; i < n; i++) {
      const have = soil[i] * area[i];
      outTotal[i] = outTotal[i] > have ? have / outTotal[i] : 1;
    }
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < Dg; k++) {
        const e = Dg * i + k;
        const q0 = flux[e];
        if (q0 === 0) continue;
        const j = nb[e];
        const q = q0 > 0 ? q0 * outTotal[i] : q0 * outTotal[j];
        soil[i] -= q / area[i];
        z[i] -= q / area[i];
        soil[j] += q / area[j];
        z[j] += q / area[j];
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
  const { grid: mesh, soil, z, boundaryOutlet } = sim;
  const n = mesh.size;
  const Dg = mesh.degree;
  const nb = mesh.neighbors;
  const { dist, area } = mesh;
  const maxIter = sim.params.hillslope.avalancheIterations;
  const ws = sim.scratch.avalanche || (sim.scratch.avalanche = {
    list: new Int32Array(n), next: new Int32Array(n), mark: new Uint32Array(n), stamp: 0,
    excess: new Float64Array(Dg),
  });
  const excess = ws.excess;
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
      const Sc = criticalSlope(sim, i);
      let sum = 0;
      let emax = 0;
      let kmax = -1;
      for (let k = 0; k < Dg; k++) {
        const e = Dg * i + k;
        const j = nb[e];
        let ex = 0;
        if (j >= 0) {
          ex = z[i] - z[j] - Sc * dist[e];
          if (ex < 0) ex = 0;
        }
        excess[k] = ex;
        sum += ex;
        if (ex > emax) { emax = ex; kmax = k; }
      }
      if (emax <= 1e-6) continue;
      // Move the volume that would level the steepest pair's excess halfway (half the
      // excess for equal areas), shared among all over-steep neighbours by excess.
      const jm = nb[Dg * i + kmax];
      const want = (emax * area[i] * area[jm]) / (area[i] + area[jm]);
      const M = Math.min(soil[i] * area[i], want); // volume
      for (let k = 0; k < Dg; k++) {
        if (excess[k] === 0) continue;
        const j = nb[Dg * i + k];
        const dv = (M * excess[k]) / sum;
        soil[j] += dv / area[j];
        z[j] += dv / area[j];
      }
      if (M === soil[i] * area[i]) { z[i] -= soil[i]; soil[i] = 0; } // exact, no round-off sliver
      else { soil[i] -= M / area[i]; z[i] -= M / area[i]; }
      moved += M;
      // i and its neighbours may now be unstable
      if (ws.mark[i] !== stamp) { ws.mark[i] = stamp; next[nextCount++] = i; }
      for (let k = 0; k < Dg; k++) {
        const j = nb[Dg * i + k];
        if (j >= 0 && ws.mark[j] !== stamp) { ws.mark[j] = stamp; next[nextCount++] = j; }
      }
    }
    [list, next] = [next, list];
    count = nextCount;
  }
  return moved;
}
