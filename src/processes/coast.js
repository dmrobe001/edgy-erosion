// Coasts: wave exposure, sea-cliff retreat, shore-platform lowering, and the
// reworking of sediment near the shore.
//
// SEA LEVEL may vary in time (seaLevelAt). Every cell below it is sea.
//
// WAVE EXPOSURE. For each cell near the shore we follow a ray towards each direction
// waves come from (the wave climate) across the sea until it hits land or leaves the
// grid. Rays run along the six lattice directions (see waveExposure). Along each ray:
//   generation   = 1 - exp(-fetch / fetchScale)       waves grow with open water upwind
//   dissipation  = exp(-shallowLength / dissipationLength)
// where shallowLength is how far the ray runs over water shallower than the breaking
// depth (waveHeight / 0.78). Exposure is the weighted mean of generation x dissipation
// over the wave climate, between 0 (sheltered) and 1 (fully exposed). Dissipation is
// the feedback that limits platform width: a wide shallow platform breaks the waves
// before they reach the cliff (Trenhaile 2000, Marine Geology 166; Anderson, Densmore &
// Ellis 1999, Basin Research 11). [adapted: a cheap stand-in for wave refraction and
// shoaling, which this does not model]
//
// CLIFF RETREAT. Any land cell with sea beside it is attacked from the side: a sea
// cliff, or on a gentle coast a low shore that will become one. Waves cut a notch at
// its foot; the notch deepens horizontally at
//   R = waveErodibility(rock at sea level) x weakness x exposure x exp(-soil / beachShielding)   (m/yr)
// (the beach-shielding exponential follows the SCAPE model, Walkden & Hall 2005,
// Coastal Engineering 52). When the notch reaches dx the cell has been eaten
// through: the column collapses to low-tide level. A cliff can retreat at most one
// cell per step, so keep dt < dx / (fastest retreat rate). A fraction `debrisRetained`
// of the collapsed volume stays as loose debris on the cell; the rest (fines carried
// off in suspension) leaves the domain. The debris shields the next cliff until the
// sea clears it. [adapted: notch-and-collapse is our grid version of lateral retreat]
//
// PLATFORM LOWERING. Submerged cells lower vertically at
//   platformRatio x R x f(depth),  f = 1 within the tidal range, decaying as
//   exp(-(depth - tide/2) / (depthScale x waveHeight)) below it
// (depth-decay of wave erosion: Sunamura 1992, "Geomorphology of Rocky Coasts";
// Anderson et al. 1999). Eroded rock becomes soil in the same proportion as above.
//
// ATTRITION. Shore sediment is ground down to fines and carried away at
//   attrition x exposure x f(depth)   (m/yr; f = 1 in the swash and tidal zone)
// so debris from cliff falls is cleared in time and beaches survive mainly where
// waves are weak or supply is high. [heuristic]
//
// SHORE SEDIMENT. Soil near the shore diffuses with diffusivity
//   D = marineDiffusivity x exposure x exp(-depth / closureDepth)
// down to well below the closure depth and up to the run-up limit above sea level,
// which builds beaches, shorefaces and delta fronts and parks sand in sheltered bays.
// (Depth-decaying marine diffusion as in basin models, e.g. Yuan et al. 2019.) Explicit
// substeps are capped; if the cap binds, D is reduced to the stable value and the
// step reports it.

const BREAK_RATIO = 0.78; // breaking wave height / water depth

/** Sea level (m) at time t (yr), or null when there is no sea. */
export function seaLevelAt(params, t) {
  const base = params.boundary.seaLevel;
  if (base == null) return null;
  const c = params.coast;
  if (Array.isArray(c.seaLevelCurve) && c.seaLevelCurve.length) {
    const pts = c.seaLevelCurve;
    if (t <= pts[0][0]) return base + pts[0][1];
    for (let k = 1; k < pts.length; k++) {
      if (t <= pts[k][0]) {
        const [t0, z0] = pts[k - 1];
        const [t1, z1] = pts[k];
        return base + z0 + ((z1 - z0) * (t - t0)) / (t1 - t0);
      }
    }
    return base + pts[pts.length - 1][1];
  }
  if (!c.seaLevelAmplitude) return base;
  return base + c.seaLevelAmplitude * Math.sin((2 * Math.PI * t) / c.seaLevelPeriod + ((c.seaLevelPhase || 0) * Math.PI) / 180);
}

/** Lattice direction k points at azimuth 60k degrees (0 = east, 60 = north-east, ...). */
const LATTICE_AZIMUTH = [0, 60, 120, 180, 240, 300];

/**
 * Split a wave climate onto the six lattice directions: each wave direction's weight
 * goes to the two nearest lattice directions, linearly in angle.
 * @returns {Float64Array} six weights summing to 1
 */
export function latticeWaveWeights(waves) {
  const w = new Float64Array(6);
  let total = 0;
  for (const { from, weight } of waves) {
    const a = (((from % 360) + 360) % 360) / 60;
    const k0 = Math.floor(a) % 6;
    const f = a - Math.floor(a);
    w[k0] += weight * (1 - f);
    w[(k0 + 1) % 6] += weight * f;
    total += weight;
  }
  if (total > 0) for (let k = 0; k < 6; k++) w[k] /= total;
  return w;
}

/** For each lattice direction, cell indices sorted so a cell comes after its neighbour in that direction. */
function directionOrders(grid) {
  if (grid._dirOrders) return grid._dirOrders;
  const orders = LATTICE_AZIMUTH.map((az) => {
    const ux = Math.cos((az * Math.PI) / 180);
    const uy = -Math.sin((az * Math.PI) / 180);
    const idx = Array.from({ length: grid.size }, (_, i) => i);
    const proj = (i) => grid.x[i] * ux + grid.y[i] * uy;
    idx.sort((a, b) => proj(b) - proj(a)); // furthest along the direction first
    return Int32Array.from(idx);
  });
  grid._dirOrders = orders;
  return orders;
}

/**
 * Wave exposure in [0, 1] for every cell in `cells`, written into `out`.
 *
 * Rays run along the six lattice directions, which land exactly on cell centres, so
 * each ray's fetch and shallow-water length follow from its neighbour's by dynamic
 * programming: one sweep of the grid per direction, O(n), instead of marching rays.
 * The price is 60-degree angular resolution: wave directions in between are shared
 * between the two nearest lattice directions.
 * @param {import('../simulation.js').Simulation} sim
 */
export function waveExposure(sim, sea, cells, count, out) {
  const cp = sim.params.coast;
  const { grid, z } = sim;
  const n = grid.size;
  const nb = grid.neighbors;
  const dx = grid.dx;
  const breakDepth = cp.waveHeight / BREAK_RATIO;
  // Refraction and diffraction bend waves round headlands and towards the shore, so a
  // share of the energy attacks from every direction. [heuristic]
  const W = latticeWaveWeights(cp.waves);
  const r = Math.min(1, Math.max(0, cp.refraction ?? 0));
  for (let k = 0; k < 6; k++) W[k] = (1 - r) * W[k] + r / 6;
  const orders = directionOrders(grid);
  const ws = sim.scratch.waves || (sim.scratch.waves = { fetch: new Float64Array(n), shallow: new Float64Array(n) });
  const { fetch, shallow } = ws;
  for (let c = 0; c < count; c++) out[cells[c]] = 0;
  for (let k = 0; k < 6; k++) {
    if (W[k] === 0) continue;
    const order = orders[k];
    for (let s = 0; s < n; s++) {
      const i = order[s];
      const j = nb[6 * i + k];
      if (j < 0) { fetch[i] = Infinity; shallow[i] = 0; continue; } // open sea beyond the grid
      const depth = sea - z[j];
      if (depth <= 0) { fetch[i] = dx; shallow[i] = 0; continue; } // land upwind
      fetch[i] = dx + fetch[j];
      shallow[i] = shallow[j] + (depth < breakDepth ? dx : 0);
    }
    for (let c = 0; c < count; c++) {
      const i = cells[c];
      const own = sea - z[i];
      const sh = shallow[i] + (own > 0 && own < breakDepth ? dx / 2 : 0);
      const gen = fetch[i] === Infinity ? 1 : 1 - Math.exp(-fetch[i] / cp.fetchScale);
      out[i] += W[k] * gen * Math.exp(-sh / cp.dissipationLength);
    }
  }
}

/**
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} dt years
 * @returns {{exported:number, collapses:number, cliffRetreat:number, substeps:number, diffusivityCapped:boolean}}
 */
export function coast(sim, dt) {
  const cp = sim.params.coast;
  const sea = sim.seaLevel;
  const { grid, lith, rock, soil, z, boundaryOutlet, upliftCum } = sim;
  const n = grid.size;
  const nb = grid.neighbors;
  const dx = grid.dx;
  const ac = grid.cellArea;
  const exposure = sim.waveExposureField;
  exposure.fill(0);
  const result = { exported: 0, collapses: 0, cliffRetreat: 0, substeps: 0, diffusivityCapped: false };
  if (sea == null) return result;

  const H = cp.waveHeight;
  const halfTide = cp.tidalRange / 2;
  const runup = halfTide + H; // top of the swash zone above mean sea level
  const activeDepth = 3 * cp.closureDepth;
  const keep = cp.debrisRetained;

  // --- the coastal band: shallow sea, and land touching the sea ---
  const band = sim.scratch.coastBand || (sim.scratch.coastBand = new Int32Array(n));
  let count = 0;
  for (let i = 0; i < n; i++) {
    const depth = sea - z[i];
    if (depth > 0) {
      if (depth < activeDepth) band[count++] = i;
      continue;
    }
    for (let k = 0; k < 6; k++) {
      const j = nb[6 * i + k];
      if (j >= 0 && z[j] < sea) { band[count++] = i; break; }
    }
  }
  waveExposure(sim, sea, band, count, exposure);

  // --- cliff notching and platform lowering ---
  const notch = sim.notch;
  for (let c = 0; c < count; c++) {
    const i = band[c];
    if (boundaryOutlet[i]) continue;
    const ex = exposure[i];
    if (ex <= 0) continue;
    if (z[i] >= sea) {
      // land touching the sea (a sea cliff, or a low shore): waves cut a notch at its foot
      const k = lith.waveErodibility[lith.rockAt(i, sea - upliftCum[i])] * lith.weakness[i];
      const R = k * ex * Math.exp(-soil[i] / cp.beachShielding);
      notch[i] += R * dt;
      result.cliffRetreat += R * dt;
      if (notch[i] >= dx) {
        const newRock = Math.min(rock[i], sea - halfTide);
        const removed = rock[i] - newRock;
        rock[i] = newRock;
        soil[i] += keep * removed;
        result.exported += (1 - keep) * removed * ac;
        // carry the remainder (at most one more cell's worth): retreat is limited to one
        // cell per step, so keep dt below dx / (fastest retreat rate)
        notch[i] = Math.min(notch[i] - dx, dx);
        result.collapses++;
      }
    } else {
      // submerged platform or shoreface: vertical lowering of the rock beneath any soil
      const depth = sea - rock[i];
      const f = depth <= halfTide ? 1 : Math.exp(-(depth - halfTide) / (cp.depthScale * H));
      const k = lith.waveErodibility[sim.surfaceRock[i]] * lith.weakness[i];
      const lower = cp.platformRatio * k * ex * f * Math.exp(-soil[i] / cp.beachShielding) * dt;
      rock[i] -= lower;
      soil[i] += keep * lower;
      result.exported += (1 - keep) * lower * ac;
      notch[i] = 0;
    }
    // attrition: waves grind shore sediment to fines, fastest in the swash zone
    if (soil[i] > 0 && cp.attrition > 0) {
      const d = sea - z[i];
      const g = d < -runup ? 0 : d <= halfTide ? 1 : Math.exp(-(d - halfTide) / (cp.depthScale * H));
      const loss = Math.min(soil[i], cp.attrition * ex * g * dt);
      soil[i] -= loss;
      result.exported += loss * ac;
    }
    z[i] = rock[i] + soil[i];
  }

  // --- shore sediment diffusion (soil only), explicit substeps on the band ---
  const D0 = cp.marineDiffusivity;
  if (D0 > 0 && count > 0) {
    const maxSub = cp.maxSubsteps;
    const dtStable = (0.5 * dx * dx) / (3 * D0); // exposure and depth factors are <= 1
    let substeps = Math.ceil(dt / dtStable);
    let Dscale = 1;
    if (substeps > maxSub) {
      Dscale = (maxSub * dtStable) / dt;
      substeps = maxSub;
      result.diffusivityCapped = true;
    }
    result.substeps = substeps;
    const h = dt / substeps;
    const coef = ((2 / 3) * h * D0 * Dscale) / (dx * dx);
    const inBand = sim.scratch.inBand || (sim.scratch.inBand = new Uint8Array(n));
    for (let c = 0; c < count; c++) inBand[band[c]] = 1;
    // Edge list with conductances, fixed for the step (they change slowly with z).
    const es = sim.scratch.coastEdges || (sim.scratch.coastEdges = {
      a: new Int32Array(6 * n), b: new Int32Array(6 * n), g: new Float64Array(6 * n),
    });
    let ne = 0;
    for (let c = 0; c < count; c++) {
      const i = band[c];
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j < 0 || (inBand[j] && j < i)) continue; // each edge once
        const d = sea - Math.max(z[i], z[j]); // depth over the shallower of the pair
        if (d < -runup) continue;
        const g = coef * 0.5 * (exposure[i] + exposure[j]) * Math.exp(-Math.max(d, 0) / cp.closureDepth);
        if (g <= 0) continue;
        es.a[ne] = i; es.b[ne] = j; es.g[ne] = g; ne++;
      }
    }
    const flux = sim.scratch.coastFlux || (sim.scratch.coastFlux = new Float64Array(6 * n));
    const outTotal = sim.scratch.outTotal;
    const { a: ea, b: eb, g: eg } = es;
    for (let sub = 0; sub < substeps; sub++) {
      for (let e = 0; e < ne; e++) { outTotal[ea[e]] = 0; outTotal[eb[e]] = 0; }
      for (let e = 0; e < ne; e++) {
        const i = ea[e], j = eb[e];
        const q = eg[e] * (z[i] - z[j]); // positive: i -> j
        const up = q > 0 ? i : j;
        // only band cells give sediment, and only if they have some
        if (!inBand[up] || soil[up] <= 0) { flux[e] = 0; continue; }
        flux[e] = q;
        outTotal[up] += Math.abs(q);
      }
      // apply with the soil limit: no cell gives more than it has
      for (let e = 0; e < ne; e++) {
        const q0 = flux[e];
        if (q0 === 0) continue;
        const i = ea[e], j = eb[e];
        const up = q0 > 0 ? i : j;
        const q = outTotal[up] > soil[up] ? (q0 * soil[up]) / outTotal[up] : q0;
        soil[i] -= q; z[i] -= q;
        soil[j] += q; z[j] += q;
      }
      for (let c = 0; c < count; c++) {
        const i = band[c];
        if (soil[i] < 0) { z[i] -= soil[i]; soil[i] = 0; } // round-off only
      }
    }
    for (let c = 0; c < count; c++) inBand[band[c]] = 0;
  }
  return result;
}
