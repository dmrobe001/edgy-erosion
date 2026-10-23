// Rivers: bedrock/soil incision by the stream power law, then sediment transport
// and deposition downstream.
//
// INCISION. Channel lowering rate E = K * A^m * S^n, where A is the effective
// drainage area (upstream area weighted by precipitation / 1 m/yr) and S the slope
// to the receiver. Rock has its own K; loose sediment (soil, alluvium, talus) has
// K_soil, reduced by vegetation. Within a step the river first strips any sediment
// on its bed at the soil rate, and only then cuts rock, for the remaining time.
//
// Why "strip, then cut" rather than a cover factor exp(-soil/H*) taken from the
// start of the step (as in SPACE, Shobe, Tucker & Barnhart 2017, Geosci. Model Dev. 10):
// each step leaves a sediment layer whose thickness scales with dt (from weathering
// and deposition), and a cover factor computed from it made incision depend strongly
// on dt -- quartering dt changed 150 kyr of evolution by tens of metres, without
// converging. Stripping costs a fraction of the step that does not depend on dt.
//
// Each phase is solved implicitly, downstream first (Braun & Willett 2013,
// Geomorphology 180-181): unconditionally stable for n = 1; Newton's method for n != 1.
//
// DEPOSITION. Eroded material travels downstream as a volumetric flux Qs (m^3/yr).
// Each cell deposits at rate D = G * Qs / A ("xi-q" erosion-deposition model of
// Davy & Lague 2009, in the dimensionless form of Yuan et al. 2019, JGR Earth Surface
// 124). G = 0 recovers pure detachment-limited incision; G ~ 1 gives alluvial
// valleys and fans. Lakes and the sea trap sediment up to their water level.
//
/**
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} dt years
 */
export function fluvial(sim, dt) {
  const p = sim.params.fluvial;
  const { grid, lith, rock, soil, z, outlet, boundaryOutlet, veg } = sim;
  const { filled, receiver, receiverDist, order, area } = sim.flow;
  const n = grid.size;
  const cellArea = grid.area;
  const m = p.m;
  const nExp = p.n;
  const Ks = p.soilErodibility;
  const G = p.deposition;
  const vegRed = sim.params.vegetation.erodibilityReduction;
  const sea = sim.seaLevel;
  const lakeTol = 1e-3;
  const eroded = sim.scratch.eroded; // lowering depth this step (m)
  const zNew = sim.scratch.zNew;

  // --- incision, downstream first ---
  for (let s = 0; s < n; s++) {
    const i = order[s];
    eroded[i] = 0;
    zNew[i] = z[i];
    if (outlet[i]) continue;
    if (filled[i] - z[i] > lakeTol) continue; // under a lake: no incision
    const r = receiver[i];
    const zr = zNew[r];
    if (z[i] <= zr) continue;
    const L = receiverDist[i];
    const Am = Math.pow(area[i], m);
    const kSoil = Ks * Math.exp(-vegRed * veg[i]);
    const kRock = lith.erodibility[sim.surfaceRock[i]] * lith.weakness[i];
    // Strip, then cut: loose sediment on the bed is entrained at the soil rate first;
    // only once it is gone does the river cut rock, for the rest of the step.
    let fromSoil = 0;
    let fromRock = 0;
    const zSoil = implicitStep(z[i], zr, kSoil * Am, dt, L, nExp);
    const soilWanted = z[i] - zSoil;
    if (soilWanted <= soil[i]) {
      fromSoil = soilWanted;
    } else {
      fromSoil = soil[i];
      // time spent stripping, in proportion to the depth stripped
      const tRest = dt * (1 - soil[i] / soilWanted);
      const z1 = z[i] - soil[i];
      if (tRest > 0 && z1 > zr) fromRock = z1 - implicitStep(z1, zr, kRock * Am, tRest, L, nExp);
    }
    soil[i] -= fromSoil;
    rock[i] -= fromRock;
    eroded[i] = fromSoil + fromRock;
    zNew[i] = z[i] - eroded[i];
  }

  // --- sediment routing and deposition, upstream first ---
  const qs = sim.scratch.qs; // incoming sediment flux (m^3/yr)
  if (sim.inflowSediment) qs.set(sim.inflowSediment); // rivers entering from outside the domain
  else qs.fill(0);
  const qsOut = sim.qsOut;
  qsOut.fill(0);
  let exported = 0;
  for (let s = n - 1; s >= 0; s--) {
    const i = order[s];
    const ac = cellArea[i];
    const supply = qs[i] + (eroded[i] * ac) / dt; // m^3/yr
    if (supply <= 0) continue;
    const maxDepth = (supply * dt) / ac; // all of it, as a thickness on this cell
    let d;
    if (outlet[i]) {
      d = !boundaryOutlet[i] && sea != null ? clamp(sea - zNew[i], 0, maxDepth) : 0;
    } else if (filled[i] - z[i] > lakeTol) {
      // a lake (judged on the pre-incision surface, which z still holds): fill to the water level
      d = Math.min(filled[i] - zNew[i], maxDepth);
    } else if (G > 0) {
      d = ((G * supply) / (area[i] + G * ac)) * dt;
    } else {
      d = 0;
    }
    soil[i] += d;
    zNew[i] += d;
    const out = supply - (d * ac) / dt;
    qsOut[i] = out;
    const r = receiver[i];
    if (r >= 0) qs[r] += out;
    else exported += out * dt;
  }
  sim.qs = qs;
  for (let i = 0; i < n; i++) z[i] = rock[i] + soil[i];
  return exported;
}

function clamp(v, lo, hi) {
  return v < lo ? lo : v > hi ? hi : v;
}

/**
 * One implicit (backward Euler) step of dz/dt = -KA (z - zr)^n / L^n for a cell whose
 * receiver has already been updated to zr (Braun & Willett 2013). Exact algebra for
 * n = 1; Newton's method otherwise. Returns the new elevation, never below zr.
 */
function implicitStep(z0, zr, KA, dt, L, n) {
  if (n === 1) {
    const F = (KA * dt) / L;
    return (z0 + F * zr) / (1 + F);
  }
  // f(h) = h - h0 + F h^n = 0 with h = z - zr, h0 = z0 - zr
  const F = (KA * dt) / Math.pow(L, n);
  const h0 = z0 - zr;
  let h = h0;
  for (let it = 0; it < 30; it++) {
    const hn1 = Math.pow(h, n - 1);
    const f = h - h0 + F * hn1 * h;
    const df = 1 + n * F * hn1;
    const hNext = Math.max(h - f / df, h * 1e-3);
    if (Math.abs(hNext - h) < 1e-9 * (1 + h0)) { h = hNext; break; }
    h = hNext;
  }
  return zr + h;
}
