// Rivers: bedrock/soil incision by the stream power law, then sediment transport
// and deposition downstream.
//
// INCISION. Channel lowering rate E = K * A^m * S^n, where A is the effective
// drainage area (upstream area weighted by precipitation / 1 m/yr), S the slope to
// the receiver, and K an erodibility that blends between soil and bedrock:
//
//   K = K_soil * (1 - c) + K_rock * c,      c = exp(-soil / soilCoverScale)
//
// (thick soil -> soil erodibility; bare rock -> rock erodibility; after the SPACE
// model of Shobe, Tucker & Barnhart 2017, Geosci. Model Dev. 10). Vegetation
// reduces K_soil. We solve implicitly in time, downstream first, following
// Braun & Willett 2013 (Geomorphology 180-181), which is unconditionally stable for
// n = 1; for n != 1 each cell's implicit equation is solved by Newton's method.
//
// DEPOSITION. Eroded material travels downstream as a volumetric flux Qs (m^3/yr).
// Each cell deposits at rate D = G * Qs / A ("xi-q" erosion-deposition model of
// Davy & Lague 2009, in the dimensionless form of Yuan et al. 2019, JGR Earth Surface
// 124). G = 0 recovers pure detachment-limited incision; G ~ 1 gives alluvial
// valleys and fans. Lakes and the sea trap sediment up to their water level.
//
// As in SPACE, the soil term only removes soil and the rock term only removes rock.
// Simplification: the implicit solve uses the blended K for the total lowering and
// then splits it in proportion; if the soil runs out mid-step the soil share is
// simply truncated (the cell lowers a little less than the implicit solution).

/**
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} dt years
 */
export function fluvial(sim, dt) {
  const p = sim.params.fluvial;
  const { grid, lith, rock, soil, z, outlet, boundaryOutlet, veg } = sim;
  const { filled, receiver, order, area } = sim.flow;
  const n = grid.size;
  const L = grid.dx;
  const ac = grid.cellArea;
  const m = p.m;
  const nExp = p.n;
  const Ks = p.soilErodibility;
  const Hs = p.soilCoverScale;
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
    const cover = Math.exp(-soil[i] / Hs);
    const kSoil = Ks * Math.exp(-vegRed * veg[i]) * (1 - cover);
    const kRock = lith.erodibility[sim.surfaceRock[i]] * lith.weakness[i] * cover;
    const K = kSoil + kRock;
    const Am = Math.pow(area[i], m);
    let zi;
    if (nExp === 1) {
      const F = (K * Am * dt) / L;
      zi = (z[i] + F * zr) / (1 + F);
    } else {
      // f(h) = h - h0 + F h^n = 0 with h = zi - zr, h0 = z[i] - zr
      const F = (K * Am * dt) / Math.pow(L, nExp);
      const h0 = z[i] - zr;
      let h = h0;
      for (let it = 0; it < 30; it++) {
        const hn1 = Math.pow(h, nExp - 1);
        const f = h - h0 + F * hn1 * h;
        const df = 1 + nExp * F * hn1;
        const hNext = Math.max(h - f / df, h * 1e-3);
        if (Math.abs(hNext - h) < 1e-9 * (1 + h0)) { h = hNext; break; }
        h = hNext;
      }
      zi = zr + h;
    }
    // Split the lowering between soil entrainment and bedrock incision in proportion to
    // their terms in K. Soil removal cannot exceed the soil present; it never eats into rock.
    const total = z[i] - zi;
    const fromSoil = Math.min(soil[i], (total * kSoil) / K);
    const fromRock = (total * kRock) / K;
    soil[i] -= fromSoil;
    rock[i] -= fromRock;
    eroded[i] = fromSoil + fromRock;
    zNew[i] = z[i] - eroded[i];
  }

  // --- sediment routing and deposition, upstream first ---
  const qs = sim.scratch.qs; // incoming sediment flux (m^3/yr)
  qs.fill(0);
  let exported = 0;
  for (let s = n - 1; s >= 0; s--) {
    const i = order[s];
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
