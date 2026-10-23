// Rockfall and cliff collapse, using the Culmann slope-stability criterion.
//
// A planar rock face of slope angle beta, made of material with cohesion C,
// friction angle phi and unit weight rho*g, can stand to a maximum height
//
//   Hc = (4 C / (rho g)) * sin(beta) cos(phi) / (1 - cos(beta - phi))
//
// and fails along a plane at angle (beta + phi) / 2. (Culmann 1866; used in
// landscape evolution models by Densmore, Ellis & Anderson 1998, JGR 103.) Hc grows
// without bound as beta -> phi, so faces no steeper than phi never fail, and
// strong rock (high C) supports tall vertical cliffs while weak rock does not.
//
// On the grid, each cell is compared with its steepest-down neighbour. The drop
// between them is the face height and atan(drop / distance) its angle. The face strength is the
// thickness-weighted mean over the layers it cuts, so a soft layer at the foot of
// a cliff weakens it -- undercutting, the engine of escarpment retreat.
//
// An unstable face fails stochastically with rate (1/timescale)(1 + frostFactor*frost),
// i.e. probability 1 - exp(-rate dt) per step. The failed wedge is removed down to
// the failure plane (but never more than half the drop, so the two cells cannot
// swap order) and dumped as loose debris (soil) on the lower cell, where the
// avalanche step spreads it into a talus slope.
//
// Grid-scale caveat: cliffs lower than about dx * tan(phi) cannot be represented as
// unstable faces, so rockfall only "switches on" when dx is small compared with the
// cliff heights the rock can sustain.

const G_ACCEL = 9.81;

/** Culmann critical height (m) for a face at angle beta (rad). Infinity if beta <= phi. */
export function culmannHeight(beta, cohesion, phi, density) {
  if (beta <= phi) return Infinity;
  return ((4 * cohesion) / (density * G_ACCEL)) * (Math.sin(beta) * Math.cos(phi)) / (1 - Math.cos(beta - phi));
}

/**
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} dt years
 * @returns {{events:number, volume:number}}
 */
export function rockfall(sim, dt) {
  const rp = sim.params.rockfall;
  const { grid, lith, rock, soil, z, boundaryOutlet, upliftCum, rng, rockfallActivity } = sim;
  const n = grid.size;
  const Dg = grid.degree;
  const nb = grid.neighbors;
  const { dist, area } = grid;
  const rate = (1 / rp.timescale) * (1 + rp.frostFactor * sim.params.climate.frost);
  const pFail = 1 - Math.exp(-rate * dt);
  let events = 0;
  let volume = 0;

  const forward = (sim.stepCount & 1) === 0;
  for (let s = 0; s < n; s++) {
    const i = forward ? s : n - 1 - s;
    if (boundaryOutlet[i]) continue;
    // the steepest face: largest drop / distance to a neighbour
    let j = -1;
    let L = 0;
    let best = 0;
    for (let k = 0; k < Dg; k++) {
      const e = Dg * i + k;
      const c = nb[e];
      if (c < 0) continue;
      const S = (z[i] - z[c]) / dist[e];
      if (S > best) { best = S; j = c; L = dist[e]; }
    }
    if (j < 0) continue;
    const zj = z[j];
    const drop = z[i] - zj;
    // only the rock part of the face carries load; the soil cap is just weight
    const rockFace = rock[i] - zj;
    if (rockFace <= 0) continue;
    const beta = Math.atan(drop / L);
    const { cohesion, phi } = lith.meanStrength(i, zj - upliftCum[i], rock[i] - upliftCum[i]);
    if (drop <= culmannHeight(beta, cohesion, phi, rp.rockDensity)) continue;
    if (rng.next() >= pFail) continue;

    const stableDrop = L * Math.tan(0.5 * (beta + phi));
    // never so much that the lower cell ends up higher (half the drop for equal areas)
    const remove = Math.min(drop - stableDrop, (drop * area[j]) / (area[i] + area[j]));
    if (remove <= 0) continue;
    const fromSoil = Math.min(soil[i], remove);
    soil[i] -= fromSoil;
    rock[i] -= remove - fromSoil;
    z[i] -= remove;
    const v = remove * area[i];
    soil[j] += v / area[j];
    z[j] += v / area[j];
    events++;
    volume += v;
    rockfallActivity[i] += remove;
  }
  return { events, volume };
}
