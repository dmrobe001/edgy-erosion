// Earthquakes: faults that move in sudden slips instead of smooth uplift.
//
// A fault is a straight surface trace (in fractions of the reference frame, so it lines
// up across zoom levels and may run off the domain: the domain is then a window onto
// a longer fault) with a dip angle and a long-term slip rate. The side to the LEFT of
// the trace's direction (seen on the map, north up) is the raised side, as for the
// `fault` uplift pattern.
//
// MAGNITUDES follow a truncated Gutenberg-Richter distribution: the number of events of
// magnitude >= M falls as 10^(-b M) between minMagnitude and maxMagnitude
// (Gutenberg & Richter 1944, Bull. Seismol. Soc. Am. 34). Below about magnitude 5.5
// ruptures rarely reach the surface, so smaller events are left out.
//
// RUPTURE SIZE from magnitude, by the empirical regressions of Wells & Coppersmith 1994
// (Bull. Seismol. Soc. Am. 84), all slip types:
//   surface rupture length   log10 L(km) = -3.22 + 0.69 M
//   average displacement     log10 D(m)  = -4.80 + 0.69 M
//   down-dip rupture width   log10 W(km) = -1.01 + 0.32 M
// Ruptures longer than the fault are cut to its length. Slip along a rupture has an
// elliptical profile with the regression's mean, zero at the tips.
//
// HOW OFTEN. Events are a Poisson process whose rate makes the long-term slip match the
// slip rate (a potency balance, the slip-times-length version of a moment balance):
//   rate = slipRate * faultLength / E[D(M) L(M)]
// with the expectation over the magnitude distribution. So over many events the fault
// slips at its slip rate on average, while any one place sees a staircase of offsets.
//
// SURFACE DISPLACEMENT. Vertical throw = slip x sin(dip). The raised side rises by
// upliftShare of it and the other side drops by the rest (on normal faults the
// hanging wall typically drops several times more than the footwall rises). Both decay
// away from the trace as exp(-distance / W), a stand-in for the elastic dislocation
// solution (Okada 1985, Bull. Seismol. Soc. Am. 75) [adapted: exponential decay instead
// of the full half-space solution; horizontal motion is ignored]. Rock carries its
// strata with it, so the fault offsets the layers.
//
// SHAKING. Each event shakes the ground within about 10^(0.5 M - 2.5) km of the rupture
// [adapted: a one-parameter attenuation, not a ground-motion model]. For that step
// cliffs fail more easily (their critical height is cut by coseismicWeakening x shaking)
// and unstable ones fail at once (see rockfall.js): earthquake-triggered rockfall.

const sinD = (deg) => Math.sin((deg * Math.PI) / 180);

/** Surface rupture length (m) for magnitude M (Wells & Coppersmith 1994). */
export function ruptureLength(M) { return 1000 * Math.pow(10, -3.22 + 0.69 * M); }
/** Average surface displacement (m) for magnitude M (Wells & Coppersmith 1994). */
export function averageDisplacement(M) { return Math.pow(10, -4.80 + 0.69 * M); }
/** Down-dip rupture width (m) for magnitude M (Wells & Coppersmith 1994). */
export function ruptureWidth(M) { return 1000 * Math.pow(10, -1.01 + 0.32 * M); }
/** Distance (m) over which shaking from a magnitude-M rupture falls by e. */
export function shakingLength(M) { return 1000 * Math.pow(10, 0.5 * M - 2.5); }

/** Magnitude from a uniform random number u in [0, 1): inverse CDF of the truncated G-R law. */
export function sampleMagnitude(u, b, Mmin, Mmax) {
  const c = 1 - Math.pow(10, -b * (Mmax - Mmin));
  return Mmin - Math.log10(1 - u * c) / b;
}

/** Per-fault geometry and event rate, cached on the simulation until parameters change. */
function faultModel(sim, k) {
  const ep = sim.params.earthquakes;
  const fault = ep.faults[k];
  const F = sim.frame;
  const key = JSON.stringify([fault, ep.bValue, ep.minMagnitude, ep.maxMagnitude, F]);
  const cache = sim.scratch.faults || (sim.scratch.faults = []);
  if (cache[k] && cache[k].key === key) return cache[k];
  const g = sim.grid;
  const x0 = F.x + fault.x0 * F.width, y0 = F.y + fault.y0 * F.height;
  const x1 = F.x + fault.x1 * F.width, y1 = F.y + fault.y1 * F.height;
  const Lf = Math.hypot(x1 - x0, y1 - y0);
  if (!(Lf > 0)) throw new Error(`earthquakes.faults[${k}] has zero length`);
  const sx = (x1 - x0) / Lf, sy = (y1 - y0) / Lf;
  const along = new Float64Array(g.size);
  const across = new Float64Array(g.size); // > 0 on the raised (left) side
  for (let i = 0; i < g.size; i++) {
    const rx = g.x[i] - x0, ry = g.y[i] - y0;
    along[i] = rx * sx + ry * sy;
    across[i] = rx * sy - ry * sx; // left of the strike, with y pointing south
  }
  // expected potency (slip x rupture length) per event, by quadrature over the
  // magnitude distribution's quantiles
  const Q = 400;
  let potency = 0;
  for (let q = 0; q < Q; q++) {
    const M = sampleMagnitude((q + 0.5) / Q, ep.bValue, ep.minMagnitude, ep.maxMagnitude);
    potency += averageDisplacement(M) * Math.min(ruptureLength(M), Lf);
  }
  potency /= Q;
  // cells sorted by position along the fault, so an event visits only the cells near it
  const order = Int32Array.from({ length: g.size }, (_, i) => i).sort((a, b) => along[a] - along[b]);
  const alongSorted = Float64Array.from(order, (i) => along[i]);
  const model = { key, x0, y0, sx, sy, length: Lf, along, across, order, alongSorted, rate: (fault.slipRate * Lf) / potency, meanPotency: potency };
  cache[k] = model;
  return model;
}

/**
 * Slip one fault now.
 * @param {import('../simulation.js').Simulation} sim
 * @param {number} k index into params.earthquakes.faults
 * @param {number} M magnitude
 * @param {number} [center] position of the rupture centre along the fault (m from its
 *   first point); default random, keeping the rupture on the fault
 * @returns {object} the event, also appended to sim.quakes.events
 */
export function triggerEarthquake(sim, k, M, center, fm = faultModel(sim, k)) {
  const ep = sim.params.earthquakes;
  const fault = ep.faults[k];
  const { rock, upliftCum, boundaryOutlet, shaking } = sim;
  const area = sim.grid.area;
  const L = Math.min(ruptureLength(M), fm.length);
  const slip = averageDisplacement(M);
  const c = center ?? L / 2 + sim.quakeRng.next() * (fm.length - L);
  const half = L / 2;
  const peak = (4 / Math.PI) * slip; // elliptical profile with mean `slip`
  const throwPeak = peak * sinD(fault.dip ?? 60);
  const share = fault.upliftShare ?? ep.upliftShare;
  const decay = fault.decayLength ?? ruptureWidth(M);
  const r0 = shakingLength(M);
  let volume = 0;
  // cells within the rupture's along-strike span, plus a margin where shaking matters
  // (only to rockfall: without it, shaking is not worth computing)
  const shake = sim.params.rockfall.enabled;
  const reach = half + (shake ? 5 * r0 : 0);
  const first = lowerBound(fm.alongSorted, c - reach);
  const last = lowerBound(fm.alongSorted, c + reach);
  for (let k = first; k < last; k++) {
    const i = fm.order[k];
    const s = fm.along[i] - c;
    const d = fm.across[i];
    // shaking: by distance to the rupture segment
    if (shake) {
      const ds = Math.max(0, Math.abs(s) - half);
      const sh = Math.exp(-Math.hypot(ds, d) / r0);
      if (sh > shaking[i]) shaking[i] = sh;
    }
    if (boundaryOutlet[i] || Math.abs(s) >= half) continue;
    const t = throwPeak * Math.sqrt(1 - (s / half) * (s / half)) * Math.exp(-Math.abs(d) / decay);
    const u = d > 0 ? share * t : -(1 - share) * t;
    rock[i] += u;
    upliftCum[i] += u;
    volume += u * area[i];
  }
  sim.budget.uplift += volume;
  const ev = {
    time: sim.time, fault: k, magnitude: M, length: L, slip, throw: slip * sinD(fault.dip ?? 60),
    center: [fm.x0 + c * fm.sx, fm.y0 + c * fm.sy],
  };
  const q = sim.quakes;
  q.events.push(ev);
  if (q.events.length > 1000) q.events.splice(0, q.events.length - 1000);
  q.count++;
  q.lastStep++;
  return ev;
}

/** First index of a sorted array whose value is >= v. */
function lowerBound(a, v) {
  let lo = 0, hi = a.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m; }
  return lo;
}

/**
 * Long-term mean rate of surface uplift (m/yr) that fault k's earthquakes produce at each
 * cell: the event rate times the expected vertical displacement per event, averaged over
 * magnitudes (quantile quadrature) and over rupture positions (uniform along the fault,
 * so slip tapers towards the fault's ends). The integral of the elliptical slip profile
 * over rupture centres is done in closed form. Used to keep a zoomed patch's long-term
 * uplift equal to its planet's while earthquakes deliver part of it.
 */
export function faultMeanUplift(sim, k) {
  const ep = sim.params.earthquakes;
  const fault = ep.faults[k];
  const fm = faultModel(sim, k);
  const Q = 200;
  const Lf = fm.length;
  const share = fault.upliftShare ?? ep.upliftShare;
  const sinDip = sinD(fault.dip ?? 60);
  // integral of sqrt(1 - t^2) from 0 to t
  const G = (t) => { t = Math.max(-1, Math.min(1, t)); return 0.5 * (t * Math.sqrt(1 - t * t) + Math.asin(t)); };
  const out = new Float64Array(sim.grid.size);
  for (let q = 0; q < Q; q++) {
    const M = sampleMagnitude((q + 0.5) / Q, ep.bValue, ep.minMagnitude, ep.maxMagnitude);
    const L = Math.min(ruptureLength(M), Lf);
    const half = L / 2;
    const peak = (4 / Math.PI) * averageDisplacement(M);
    const W = fault.decayLength ?? ruptureWidth(M);
    const cmin = half, cmax = Lf - half;
    for (let i = 0; i < out.length; i++) {
      const s = fm.along[i];
      let slip; // expected slip at s per event of this magnitude
      if (cmax - cmin < 1e-9 * Lf) {
        const t = (s - Lf / 2) / half;
        slip = Math.abs(t) < 1 ? peak * Math.sqrt(1 - t * t) : 0;
      } else {
        const lo = Math.max(cmin, s - half), hi = Math.min(cmax, s + half);
        slip = hi > lo ? (peak * half * (G((hi - s) / half) - G((lo - s) / half))) / (cmax - cmin) : 0;
      }
      if (slip === 0) continue;
      const d = fm.across[i];
      out[i] += slip * Math.exp(-Math.abs(d) / W) * (d > 0 ? share : -(1 - share));
    }
  }
  const scale = (fm.rate * sinDip) / Q;
  for (let i = 0; i < out.length; i++) out[i] *= scale;
  return out;
}

/**
 * Draw this step's earthquakes on every fault and apply them.
 * @returns {number} events this step
 */
export function earthquakes(sim, dt) {
  const ep = sim.params.earthquakes;
  sim.quakes.lastStep = 0;
  for (let k = 0; k < ep.faults.length; k++) {
    const fm = faultModel(sim, k);
    if (!(fm.rate > 0)) continue;
    // Poisson process: exponential waiting times
    let t = -Math.log(1 - sim.quakeRng.next()) / fm.rate;
    while (t < dt) {
      const M = sampleMagnitude(sim.quakeRng.next(), ep.bValue, ep.minMagnitude, ep.maxMagnitude);
      triggerEarthquake(sim, k, M, undefined, fm);
      t += -Math.log(1 - sim.quakeRng.next()) / fm.rate;
    }
  }
  return sim.quakes.lastStep;
}

/** Long-term event rate (events/yr) and mean potency of fault k, for checks and display. */
export function faultRate(sim, k) {
  const fm = faultModel(sim, k);
  return { rate: fm.rate, length: fm.length, meanPotency: fm.meanPotency };
}
