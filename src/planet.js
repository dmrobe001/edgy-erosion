// A planet: kinematic plate tectonics on the icosahedral sphere, isostasy, and the
// same surface processes (rivers, creep, weathering) as the landscape scale.
//
// PLATES are rigid caps that rotate about Euler poles at constant angular velocity
// ("kinematic": the motions are prescribed, not computed from forces). Each plate keeps
// its crust in its OWN frame of reference, on its own copy of the mesh: crust thickness,
// continental fraction, formation time and soil. The world map is a view that, for each
// world cell, looks back along every plate's rotation to the plate-frame cell that has
// moved there (nearest neighbour). Rigid motion is therefore exact: nothing is
// interpolated, so coastlines do not blur however far a plate travels.
//
// Where plates meet, after every small rotation:
//   no plate covers a cell   a GAP opens: new oceanic crust (age zero) is added to the
//                            plate that was there before. This is sea-floor spreading;
//                            ridges form and ocean floor ages away from them.
//   two or more plates       they CONVERGE. The most buoyant column (the one that would
//                            float highest) overrides. A continent beneath it is too
//                            buoyant to sink and is scraped off onto the overriding plate
//                            (collision: crust thickens into mountains); ocean floor
//                            subducts (is destroyed) except for part of its sediment, while
//                            a fraction `arcFraction` of the subducted volume returns as
//                            new continental crust on the overriding plate (volcanic arcs;
//                            this is how island arcs grow into continents).
// Colliding continents lock together: their relative rotation decays with timescale
// `collisionLocking` while continental crust is being scraped off along their contact
// (an adapted stand-in for the forces that stall real collisions). Optionally, every
// `reorganizeEvery` years all plates get new random poles (plate reorganisation).
// Thick crust spreads sideways where it exceeds `flowThreshold` (crustal flow, the
// gravitational spreading that keeps orogenic plateaus flat).
//
// ISOSTASY. Each column floats on the mantle (Airy isostasy, no flexure). A column of
// thickness h and density rho has a "dry" freeboard  e = h (1 - rho / rho_m) - K,  and
// below sea level the water's weight deepens it to  e rho_m / (rho_m - rho_w).  K is set so
// that continental crust of `continentalThickness` stands at `referenceElevation`. The
// density mixes oceanic and continental crust by the continental fraction. Oceanic
// lithosphere cools and sinks with age t:  depth = ridgeDepth + subsidence sqrt(t / 1 Myr)
// until it reaches the depth of old cold crust (Parsons & Sclater 1977, J. Geophys. Res. 82).
// Erosion therefore triggers rebound: removing 1 m of rock lowers a continent by only
// 1 - rho_c / rho_m ~ 0.15 m, which is why old mountain belts outlive their erosion.
//
// SEA LEVEL is set every step so the ocean's water volume stays constant: young, buoyant
// ocean floor or shrinking continents push it up; old deep basins pull it down.
//
// SURFACE processes are a Simulation on the same mesh with planet-scale parameters
// (rivers, creep, weathering; no rockfall or waves at 100-km cells). Its erosion and
// deposition change crust thickness, written back to the plate frames.
//
// Conservation: in the plate frames, total crust volume changes by exactly
//   created (ridges) + arcs - subducted + surface change (sediment lost to the sea).
// The world view samples the plate frames by nearest neighbour; where neighbouring
// world cells look back to the same plate cell (a few per cent of cells, along seams
// where the rotated mesh does not line up), the map shows that column twice and skips
// another for that step.

import { SphereMesh, lonLatToVector } from './sphere.js';
import { Simulation } from './simulation.js';
import { deepMerge } from './params.js';
import { Rng, fbm3 } from './rng.js';
import { MinHeap } from './flow.js';

const YR_PER_MYR = 1e6;

/**
 * Crust material classes, which seed lithology at finer levels. A column's class is
 * decided from what it is made of (its arc and accreted shares) and its setting.
 */
export const MATERIAL = {
  oceanic: 0,         // basaltic ocean floor
  craton: 1,          // the original continental crust, never reworked
  terrane: 2,         // crust scraped off or collided onto a plate (accreted terrane)
  arc: 3,             // volcanic arc crust grown above a subduction zone
  basin: 4,           // continental crust that is mostly sediment, or drowned: sedimentary basin
};
export const MATERIAL_NAMES = Object.keys(MATERIAL);

export const PLANET_DEFAULTS = {
  seed: 1,
  subdivisions: 40,
  surfaceProcesses: true,        // false: no river routing; land is worn down by `denudation` instead
  denudation: {
    timescale: 1.5e8,            // yr: with surfaceProcesses false, land higher than `floor` above sea level
    floor: 400,                  // m      loses height at (z - sea level - floor) / timescale
  },              // icosahedron edge cuts: 10 f^2 + 2 cells (40 -> 16 002, ~180 km apart on Earth)
  radius: 6.371e6,               // m
  time: { dt: 5e5 },             // years per step
  plates: {
    count: 10,
    minSpeed: 0.01,              // m/yr, fastest point of the slowest plate
    maxSpeed: 0.08,              // m/yr, fastest point of the fastest plate
    sizeVariation: 1.0,          // spread of plate sizes (0 = all similar)
    collisionLocking: 2e7,       // yr: relative motion of colliding continents decays over this
    lockingLength: 1.5e6,        // m of collision front over which locking reaches full strength
    reorganizeEvery: 0,          // yr between plate reorganisations (new poles); 0 = never
  },
  continents: {
    mode: 'scattered',           // 'scattered' (noise blobs) | 'supercontinent' (one cap)
    fraction: 0.35,              // share of the surface that starts continental
    wavelength: 1.2,             // noise wavelength in planet radii (smaller = more, smaller continents)
    thicknessNoise: 4000,        // m
  },
  crust: {
    continentalThickness: 35000, // m
    oceanicThickness: 7000,      // m
    continentalDensity: 2800,    // kg/m^3
    oceanicDensity: 2900,
    mantleDensity: 3300,
    waterDensity: 1030,
    referenceElevation: 600,     // m: where continental crust of continentalThickness floats
    ridgeDepth: 2600,            // m below sea level at a spreading ridge
    subsidence: 350,             // m per sqrt(Myr) of oceanic subsidence with age
    trenchDepth: 2500,           // m: extra depth of a subducting plate's edge (ocean trenches)
    arcFraction: 0.1,            // volume of new arc crust per volume of oceanic crust subducted
    sedimentAccretion: 0.5,      // share of a sinking ocean floor's sediment scraped off at the trench
    arcWidth: 250000,            // m: width over which arc and scraped-off crust spreads behind the margin
    polarityFlip: 2000,          // m: buoyancy excess of the sinking column that flips a subduction zone
    flowThreshold: 45000,        // m: crust thicker than this spreads sideways
    flowDiffusivity: 5e4,        // m^2/yr: how fast it spreads at the threshold...
    flowScale: 10000,            // m: ...rising as (1 + (excess / flowScale)^2) for thicker crust
    maxFlowSubsteps: 400,        // cap on explicit flow substeps; when it binds, the fastest flow is capped instead
    maxThickness: 0,             // m: crust thicker than this sheds the excess outward, like arc crust; 0 = no limit
  },
  events: {
    log: true,                   // record collisions and continental rifting as dated events
    minCells: 2,                 // boundary cells of a kind, per plate pair and step, for an event to register
    gap: 1e7,                    // yr: an event quiet for longer than this has ended; shorter pauses join up
  },
  climate: {
    zonal: true,                 // wet tropics, dry subtropics, wet mid-latitudes, dry poles
    precipitation: 1.0,          // m/yr, scales the whole pattern
  },
  // Surface processes: parameters for the Simulation on the sphere (see params.js).
  // Erodibilities are effective values for ~100-km cells, far below landscape-scale ones:
  // the cells hide all the relief that does the eroding at small scales.
  surface: {
    rocks: { crust: { color: '#9a8a74', erodibility: 1.5e-7, cohesion: 1e8, frictionAngle: 45, weathering: 1e-4, waveErodibility: 0 } },
    strata: { layers: [{ rock: 'crust', thickness: 1e9 }], basement: 'crust', strengthVariation: { amplitude: 0 }, undulation: { amplitude: 0 } },
    fluvial: { soilErodibility: 1e-6, deposition: 0.5 },
    hillslope: { diffusivity: 1000, transportDepth: 20 },
    weathering: { productionDepth: 30, moistureExponent: 0.5 },
    vegetation: { cover: 'auto' },
    climate: { frost: 0.3 },
  },
};

export const PLANET_PRESETS = {
  'earth-like': {
    description: 'About a third continental crust in scattered blocks, ten plates moving 1-8 cm/yr: ridges, trenches, arcs and collisions.',
  },
  'supercontinent': {
    description: 'One supercontinent over several plates, which tear it apart: rifts open into new oceans while its far side crumples.',
    continents: { mode: 'supercontinent', fraction: 0.3 },
    plates: { count: 8 },
  },
  'archipelago': {
    description: 'Mostly ocean with small continental fragments and many plates: island arcs grow, collide and weld into new land.',
    continents: { fraction: 0.12, wavelength: 0.5 },
    plates: { count: 14 },
  },
};

/** DEFAULTS <- preset <- overrides, like resolveParams for landscapes. */
export function resolvePlanetParams(input = {}) {
  const over = typeof input === 'string' ? { preset: input } : input;
  let p = deepMerge(PLANET_DEFAULTS, {});
  if (over.preset) {
    const pre = PLANET_PRESETS[over.preset];
    if (!pre) throw new Error(`unknown planet preset "${over.preset}"; try one of: ${Object.keys(PLANET_PRESETS).join(', ')}`);
    const { description, ...rest } = pre;
    p = deepMerge(p, rest);
  }
  p = deepMerge(p, over);
  const fail = (m) => { throw new Error(`invalid planet parameters: ${m}`); };
  if (!(p.subdivisions >= 2 && p.subdivisions <= 160)) fail('subdivisions must be 2..160');
  if (!(p.plates.count >= 1)) fail('plates.count must be at least 1');
  if (!(p.continents.fraction >= 0 && p.continents.fraction <= 1)) fail('continents.fraction must be 0..1');
  if (!(p.time.dt > 0)) fail('time.dt must be positive');
  return p;
}

/** The Simulation parameters for the surface processes on the sphere. */
function surfaceParams(p, seaLevel) {
  return deepMerge(p.surface, {
    seed: p.seed,
    time: { dt: p.time.dt },
    boundary: { open: [], seaLevel },
    initial: { type: 'flat', base: 0, height: 0, noise: { amplitude: 0 }, soil: 0 },
    uplift: { rate: 0, pattern: 'uniform' },
    coast: { enabled: false },
    rockfall: { enabled: false },
    earthquakes: { enabled: false, faults: [] },
  });
}

/** Rotation matrix (row-major 3x3) for angle theta about unit axis (x, y, z): Rodrigues' formula. */
function axisAngle(x, y, z, theta, out) {
  const c = Math.cos(theta), s = Math.sin(theta), C = 1 - c;
  out[0] = c + x * x * C; out[1] = x * y * C - z * s; out[2] = x * z * C + y * s;
  out[3] = y * x * C + z * s; out[4] = c + y * y * C; out[5] = y * z * C - x * s;
  out[6] = z * x * C - y * s; out[7] = z * y * C + x * s; out[8] = c + z * z * C;
  return out;
}

export class Planet {
  /** @param {object|string} [params] planet parameter overrides (optionally `preset`), or a preset name */
  constructor(params = {}) {
    this.params = resolvePlanetParams(params);
    const p = this.params;
    this.mesh = new SphereMesh(p.subdivisions, p.radius);
    const mesh = this.mesh;
    const n = mesh.size;
    this.rng = new Rng(p.seed);
    this.time = 0;
    this.stepCount = 0;
    this.seaLevel = 0;
    this.sim = new Simulation(surfaceParams(p, 0), { mesh });

    const P = p.plates.count;
    this.plateCount = P;
    this.rot = new Float64Array(9 * P);
    for (let q = 0; q < P; q++) { this.rot[9 * q] = this.rot[9 * q + 4] = this.rot[9 * q + 8] = 1; }
    this.omega = new Float64Array(3 * P); // angular velocity vectors (rad/yr)
    // plate-frame state, P copies of the mesh, index q * n + s
    this.bOwned = new Uint8Array(P * n);
    this.bConsumed = new Uint8Array(P * n);
    this.bStamp = new Int32Array(P * n);
    this.bThick = new Float64Array(P * n);
    this.bCont = new Float64Array(P * n);
    this.bBirth = new Float64Array(P * n);
    this.bSoil = new Float64Array(P * n);
    // share of each column's crust that is volcanic arc, and that was scraped off or
    // collided onto it (accreted terrane); the rest of its continental part is original craton
    this.bArc = new Float64Array(P * n);
    this.bAcc = new Float64Array(P * n);
    this.bSed = new Float64Array(P * n); // share that is sediment worn off land (denudation mode)
    this.bodyOf = new Int32Array(P * n);
    for (let q = 0; q < P; q++) for (let i = 0; i < n; i++) this.bodyOf[q * n + i] = i;
    this.stamp = 0;
    // world view
    this.owner = new Int32Array(n);
    this.prevOwner = new Int32Array(n);
    this.thickness = new Float64Array(n);
    this.continental = new Float64Array(n);
    this.age = new Float64Array(n);
    this.arcShare = new Float64Array(n);
    this.accretedShare = new Float64Array(n);
    this.sedimentShare = new Float64Array(n);
    /** 0 interior, 1 convergent, 2 divergent, 3 transform (largest-rate neighbour) */
    this.boundary = new Uint8Array(n);
    /** isostatic + tectonic elevation before surface processes, this step (m) */
    this.zTectonic = new Float64Array(n);
    /** rock uplift rate from tectonics and isostasy over the last step, following the plate (m/yr) */
    this.upliftRate = new Float64Array(n);
    /** total surface change rate over the last step, following the plate (m/yr) */
    this.dzdt = new Float64Array(n);
    /** crust (m^3, continental) scraped off or grown as arcs this step, not yet placed */
    this.pending = new Float64Array(n);
    /** the arc part of `pending` (m^3) */
    this.pendingArc = new Float64Array(n);
    /** dated tectonic events: {kind, start, end, center [x, y, z], radius (rad), plates [a, b], cells} */
    this.events = [];
    this._openEvents = new Map();
    /** plate pair -> the plate that overrides at their convergent boundary */
    this.polarity = new Map();
    this.budget = { created: 0, arc: 0, subducted: 0, accreted: 0, surface: 0 };
    this.last = { created: 0, subducted: 0, arc: 0, accreted: 0, substeps: 0, collisions: 0 };

    this._setupPlates();
    this._setupCrust();
    this._buildView();
    this._initialAges();
    this._buildView();
    this._precipitation();
    // water volume: whatever lies below sea level 0 at the start
    this._isostasy();
    this.waterVolume = this._waterBelow(0);
    this.seaLevel = 0;
    this._applyToSurface();
    /** surface elevation of every plate column when it was last seen (m) */
    this.bZ = new Float64Array(P * n);
    for (let i = 0; i < n; i++) this.bZ[this.owner[i] * n + i] = this.sim.z[i];
    this.volume0 = this.crustVolume();
  }

  /** Total crust volume held by all plates, in their own frames (m^3). */
  crustVolume() {
    const n = this.mesh.size;
    const A = this.mesh.area;
    let v = 0;
    for (let q = 0; q < this.plateCount; q++) {
      for (let s = 0; s < n; s++) if (this.bOwned[q * n + s]) v += this.bThick[q * n + s] * A[s];
    }
    return v;
  }

  /** Isostatic surface elevation (m, relative to the datum sea level 0). */
  elevation(h, c, ageYr) {
    const cr = this.params.crust;
    const rm = cr.mantleDensity;
    const rho = cr.oceanicDensity + c * (cr.continentalDensity - cr.oceanicDensity);
    const K = cr.continentalThickness * (1 - cr.continentalDensity / rm) - cr.referenceElevation;
    const dry = h * (1 - rho / rm) - K;
    let e = dry >= 0 ? dry : (dry * rm) / (rm - cr.waterDensity);
    if (c < 1) e += (1 - c) * this._thermal(ageYr);
    return e;
  }

  /** Thermal buoyancy of young oceanic lithosphere (m, water-loaded), zero once old. */
  _thermal(ageYr) {
    const cr = this.params.crust;
    if (this._theta0 === undefined || this._theta0key !== cr) {
      // so that new crust sits at ridgeDepth and old crust at its cold isostatic depth
      this._theta0key = cr;
      const rm = cr.mantleDensity;
      const K = cr.continentalThickness * (1 - cr.continentalDensity / rm) - cr.referenceElevation;
      const dry = cr.oceanicThickness * (1 - cr.oceanicDensity / rm) - K;
      const old = dry >= 0 ? dry : (dry * rm) / (rm - cr.waterDensity);
      this._theta0 = Math.max(0, -cr.ridgeDepth - old);
    }
    return Math.max(0, this._theta0 - this.params.crust.subsidence * Math.sqrt(Math.max(0, ageYr) / YR_PER_MYR));
  }

  /** Surface velocity (m/yr, 3-vector) of plate q at unit position (x, y, z). */
  velocity(q, x, y, z, out = [0, 0, 0]) {
    const w = this.omega, R = this.params.radius;
    const wx = w[3 * q], wy = w[3 * q + 1], wz = w[3 * q + 2];
    out[0] = R * (wy * z - wz * y); out[1] = R * (wz * x - wx * z); out[2] = R * (wx * y - wy * x);
    return out;
  }

  /** Advance one step (default params.time.dt years). */
  step(dt = this.params.time.dt) {
    const p = this.params;
    const mesh = this.mesh;
    const n = mesh.size;
    const R = p.radius;
    // tectonics, in substeps that move no point more than half a cell
    let vmax = 0;
    for (let q = 0; q < this.plateCount; q++) vmax = Math.max(vmax, R * Math.hypot(this.omega[3 * q], this.omega[3 * q + 1], this.omega[3 * q + 2]));
    const sub = Math.max(1, Math.ceil((vmax * dt) / (0.5 * mesh.spacing)));
    const h = dt / sub;
    const last = this.last;
    last.created = last.subducted = last.arc = last.accreted = last.collisions = 0;
    for (let k = 0; k < sub; k++) this._tectonicSubstep(h, this.time + (k + 1) * h);
    last.substeps = sub;
    this._maybeReorganize(dt);

    this._spreadAdditions();
    this._buildView();
    this._crustalFlow(dt);
    this._buildView();
    this._classifyBoundaries();
    this._noteRifts();
    this._closeStepEvents(this.time, this.time + dt);
    this._isostasy();
    this.seaLevel = this._seaLevelForVolume(this.waterVolume);
    this._applyToSurface();

    // surface processes
    const sim = this.sim;
    const zT = this.zTectonic;
    if (p.surfaceProcesses) sim.step(dt);
    else this._denude(dt);
    // write erosion and deposition back to the plates
    const A = mesh.area;
    let surf = 0;
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      const s = this.bodyOf[q * n + i];
      const b = q * n + s;
      // sediment a river carries past its first sea cell (once that is filled to sea level)
      // stays in that column and builds a delta, instead of leaving the planet; isostasy
      // places it next step
      const delta = sim.exportedAt[i] / A[i];
      const d = sim.z[i] - zT[i] + delta;
      const dv = d * A[i];
      surf += dv;
      const hOld = this.bThick[b];
      const hNew = Math.max(0, hOld + dv / A[s]);
      if (d > 0) this.bCont[b] = (this.bCont[b] * hOld + (hNew - hOld)) / hNew; // sediment is continental material
      this.bThick[b] = hNew;
      this.bSoil[b] = Math.min(sim.soil[i] + delta, hNew);
    }
    this.budget.surface += surf;
    // Rates follow the rock, not the map: compare each column with itself a step ago, in
    // its plate's frame. (Measured at fixed map points, a mountain range sliding past would
    // read as uplift on one flank and collapse on the other.)
    const tStart = this.time;
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      const b = q * n + this.bodyOf[q * n + i];
      const prev = this.bBirth[b] > tStart ? zT[i] : this.bZ[b]; // new crust: no history
      this.dzdt[i] = (sim.z[i] - prev) / dt;
      this.upliftRate[i] = (zT[i] - prev) / dt;
    }
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      this.bZ[q * n + this.bodyOf[q * n + i]] = sim.z[i];
    }
    this.time += dt;
    this.stepCount++;
    this._buildView();
    return this;
  }

  run(steps, onStep) {
    for (let k = 0; k < steps; k++) { this.step(); if (onStep) onStep(this); }
    return this;
  }

  /** Change parameters mid-run (crust, climate, surface, plate locking/reorganisation). */
  setParams(over) {
    const next = resolvePlanetParams(deepMerge(this.params, over));
    for (const k of ['subdivisions', 'radius', 'seed']) if (next[k] !== this.params[k]) throw new Error(`planet ${k} cannot change mid-run; make a new Planet`);
    if (next.plates.count !== this.params.plates.count) throw new Error('plate count cannot change mid-run; make a new Planet');
    this.params = next;
    this.sim.setParams(deepMerge(surfaceParams(next, this.seaLevel), {}));
    this._precipitation();
    return this;
  }

  /**
   * Replace the plate layout, keeping the crust where it is: `owner(i)` gives the plate
   * of world cell i (0..plates.count-1) and `omega[q]` plate q's angular velocity vector
   * (rad/yr; its direction is the Euler pole). For custom set-ups and tests.
   */
  setPlates(owner, omega) {
    const n = this.mesh.size;
    const P = this.plateCount;
    this._buildView();
    const soil = new Float64Array(n), birth = new Float64Array(n), zNow = Float64Array.from(this.sim.z);
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      soil[i] = this.bSoil[q * n + this.bodyOf[q * n + i]];
      birth[i] = this.bBirth[q * n + this.bodyOf[q * n + i]];
    }
    const arc = Float64Array.from(this.arcShare), acc = Float64Array.from(this.accretedShare), sed = Float64Array.from(this.sedimentShare);
    for (const a of [this.bOwned, this.bConsumed, this.bThick, this.bCont, this.bBirth, this.bSoil, this.bZ, this.bArc, this.bAcc, this.bSed]) a.fill(0);
    this.rot.fill(0);
    for (let q = 0; q < P; q++) { this.rot[9 * q] = this.rot[9 * q + 4] = this.rot[9 * q + 8] = 1; }
    for (let q = 0; q < P; q++) for (let i = 0; i < n; i++) this.bodyOf[q * n + i] = i;
    for (let i = 0; i < n; i++) {
      const q = owner(i);
      if (!(q >= 0 && q < P)) throw new Error(`setPlates: cell ${i} given plate ${q}`);
      const b = q * n + i;
      this.owner[i] = this.prevOwner[i] = q;
      this.bOwned[b] = 1;
      this.bThick[b] = this.thickness[i];
      this.bCont[b] = this.continental[i];
      this.bBirth[b] = birth[i];
      this.bSoil[b] = soil[i];
      this.bZ[b] = zNow[i];
      this.bArc[b] = arc[i];
      this.bAcc[b] = acc[i];
      this.bSed[b] = sed[i];
    }
    omega.forEach((w, q) => { this.omega[3 * q] = w[0]; this.omega[3 * q + 1] = w[1]; this.omega[3 * q + 2] = w[2]; });
    this.polarity.clear();
    this._buildView();
    this.volume0 = this.crustVolume() - (this.budget.created + this.budget.arc - this.budget.subducted + this.budget.surface);
    return this;
  }

  /** Summary numbers. */
  stats() {
    const n = this.mesh.size;
    const A = this.mesh.area;
    const z = this.sim.z;
    let land = 0, total = 0, zmin = Infinity, zmax = -Infinity, cont = 0;
    const alive = new Set();
    for (let i = 0; i < n; i++) {
      total += A[i];
      if (z[i] > this.seaLevel) land += A[i];
      if (this.continental[i] > 0.5) cont += A[i];
      zmin = Math.min(zmin, z[i]);
      zmax = Math.max(zmax, z[i]);
      alive.add(this.owner[i]);
    }
    return {
      time: this.time, step: this.stepCount, plates: alive.size, seaLevel: this.seaLevel,
      landFraction: land / total, continentalFraction: cont / total, minElevation: zmin, maxElevation: zmax,
      created: this.last.created, subducted: this.last.subducted, arc: this.last.arc, accreted: this.last.accreted,
    };
  }

  // ------------------------------------------------------------------ set-up

  _setupPlates() {
    const p = this.params;
    const mesh = this.mesh;
    const n = mesh.size;
    const P = this.plateCount;
    const rng = this.rng;
    // seeds spread out: rejection of candidates too close to earlier seeds
    const minSep = Math.acos(1 - 2 / P) * 0.6;
    const seeds = [];
    for (let tries = 0; seeds.length < P && tries < 10000; tries++) {
      const c = Math.floor(rng.next() * n);
      const ok = seeds.every((s) => {
        const d = mesh.pos[3 * s] * mesh.pos[3 * c] + mesh.pos[3 * s + 1] * mesh.pos[3 * c + 1] + mesh.pos[3 * s + 2] * mesh.pos[3 * c + 2];
        return Math.acos(Math.min(1, d)) > minSep;
      });
      if (ok || tries > 5000) seeds.push(c);
    }
    // grow plates outward at different speeds over a noisy cost: irregular, varied plates
    const speed = seeds.map(() => Math.exp(p.plates.sizeVariation * (rng.next() - 0.5) * 1.5));
    const cost = new Float64Array(n);
    const P3 = mesh.pos;
    for (let i = 0; i < n; i++) {
      cost[i] = 1 + 0.6 * fbm3(P3[3 * i], P3[3 * i + 1], P3[3 * i + 2], { wavelength: 0.35, octaves: 4, seed: p.seed + 11 });
    }
    const label = new Int32Array(n).fill(-1);
    const best = new Float64Array(n).fill(Infinity);
    const heap = new MinHeap(6 * n + P);
    seeds.forEach((s, q) => { best[s] = 0; heap.push(s, 0); label[s] = q; });
    const done = new Uint8Array(n);
    while (heap.n) {
      const i = heap.pop();
      if (done[i]) continue;
      done[i] = 1;
      const q = label[i];
      for (let k = 0; k < 6; k++) {
        const e = 6 * i + k;
        const j = mesh.neighbors[e];
        if (j < 0 || done[j]) continue;
        const d = best[i] + (mesh.dist[e] * 0.5 * (cost[i] + cost[j])) / speed[q];
        if (d < best[j]) { best[j] = d; label[j] = q; heap.push(j, d); }
      }
    }
    for (let i = 0; i < n; i++) {
      const q = label[i];
      this.owner[i] = q;
      this.prevOwner[i] = q;
      this.bOwned[q * n + i] = 1;
    }
    for (let q = 0; q < P; q++) this._randomMotion(q);
  }

  _randomMotion(q) {
    const p = this.params.plates;
    const rng = this.rng;
    // uniform random pole on the sphere
    const z = 2 * rng.next() - 1;
    const a = 2 * Math.PI * rng.next();
    const r = Math.sqrt(1 - z * z);
    const v = p.minSpeed + (p.maxSpeed - p.minSpeed) * rng.next();
    const w = v / this.params.radius;
    this.omega[3 * q] = w * r * Math.cos(a);
    this.omega[3 * q + 1] = w * r * Math.sin(a);
    this.omega[3 * q + 2] = w * z;
  }

  _setupCrust() {
    const p = this.params;
    const mesh = this.mesh;
    const n = mesh.size;
    const P3 = mesh.pos;
    const cp = p.continents;
    const cr = p.crust;
    const field = new Float64Array(n);
    if (cp.mode === 'supercontinent') {
      const [cx, cy, cz] = lonLatToVector(360 * this.rng.next() - 180, 60 * this.rng.next() - 30);
      for (let i = 0; i < n; i++) {
        const d = P3[3 * i] * cx + P3[3 * i + 1] * cy + P3[3 * i + 2] * cz;
        field[i] = d + 0.25 * fbm3(P3[3 * i], P3[3 * i + 1], P3[3 * i + 2], { wavelength: 0.6, octaves: 5, seed: p.seed + 23 });
      }
    } else if (cp.mode === 'scattered') {
      for (let i = 0; i < n; i++) field[i] = fbm3(P3[3 * i], P3[3 * i + 1], P3[3 * i + 2], { wavelength: cp.wavelength, octaves: 5, seed: p.seed + 23 });
    } else {
      throw new Error(`unknown continents.mode "${cp.mode}"`);
    }
    // threshold at the area quantile that gives the requested continental fraction
    const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => field[b] - field[a]);
    const total = mesh.area.reduce((s, a) => s + a, 0);
    const cont = new Float64Array(n);
    let acc = 0;
    for (const i of order) { if (acc >= cp.fraction * total) break; cont[i] = 1; acc += mesh.area[i]; }
    const thick = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const noise = cp.thicknessNoise * fbm3(P3[3 * i], P3[3 * i + 1], P3[3 * i + 2], { wavelength: 0.3, octaves: 4, seed: p.seed + 29 });
      thick[i] = cont[i] ? cr.continentalThickness + noise : cr.oceanicThickness;
    }
    // soften margins into shelves and slopes: two passes of neighbour averaging
    for (let pass = 0; pass < 2; pass++) {
      const t0 = Float64Array.from(thick), c0 = Float64Array.from(cont);
      for (let i = 0; i < n; i++) {
        let st = t0[i], sc = c0[i], m = 1;
        for (let k = 0; k < 6; k++) { const j = mesh.neighbors[6 * i + k]; if (j >= 0) { st += t0[j]; sc += c0[j]; m++; } }
        // keep interiors as they are; only mixed neighbourhoods change
        thick[i] = 0.5 * t0[i] + 0.5 * (st / m);
        cont[i] = 0.5 * c0[i] + 0.5 * (sc / m);
      }
    }
    for (let i = 0; i < n; i++) {
      const b = this.owner[i] * n + i;
      this.bThick[b] = thick[i];
      this.bCont[b] = cont[i];
      this.bBirth[b] = -100 * YR_PER_MYR; // replaced by _initialAges
    }
  }

  /** Initial ocean-floor ages: travel time from the nearest spreading boundary of the plate. */
  _initialAges() {
    const mesh = this.mesh;
    const n = mesh.size;
    this._classifyBoundaries();
    const best = new Float64Array(n).fill(Infinity);
    const heap = new MinHeap(6 * n + n);
    const v = [0, 0, 0];
    for (let i = 0; i < n; i++) {
      if (this.boundary[i] !== 2) continue;
      best[i] = 0;
      heap.push(i, 0);
    }
    const done = new Uint8Array(n);
    while (heap.n) {
      const i = heap.pop();
      if (done[i]) continue;
      done[i] = 1;
      const q = this.owner[i];
      const P3 = mesh.pos;
      this.velocity(q, P3[3 * i], P3[3 * i + 1], P3[3 * i + 2], v);
      const speed = Math.max(0.005, Math.hypot(v[0], v[1], v[2])); // half-spreading rate proxy
      for (let k = 0; k < 6; k++) {
        const e = 6 * i + k;
        const j = mesh.neighbors[e];
        if (j < 0 || done[j] || this.owner[j] !== q) continue;
        const d = best[i] + mesh.dist[e] / speed;
        if (d < best[j]) { best[j] = d; heap.push(j, d); }
      }
    }
    const maxAge = 180 * YR_PER_MYR;
    for (let i = 0; i < n; i++) {
      const b = this.owner[i] * n + i;
      this.bBirth[b] = -Math.min(maxAge, Number.isFinite(best[i]) ? best[i] : maxAge);
    }
  }

  _precipitation() {
    const c = this.params.climate;
    const mesh = this.mesh;
    const P = new Float64Array(mesh.size);
    for (let i = 0; i < mesh.size; i++) {
      const lat = Math.abs((mesh.lat[i] * 180) / Math.PI);
      // [adapted] a smooth caricature of Earth's zonal-mean rainfall
      const f = c.zonal ? 0.25 + 1.5 * Math.exp(-((lat / 12) ** 2)) + 0.65 * Math.exp(-(((lat - 50) / 13) ** 2)) : 1;
      P[i] = c.precipitation * f;
    }
    this.sim.setFields({ precipitation: P });
  }

  // ------------------------------------------------------------------ tectonics

  _tectonicSubstep(h, tNow) {
    const mesh = this.mesh;
    const n = mesh.size;
    const P = this.plateCount;
    const P3 = mesh.pos;
    const { rot, omega, bodyOf, bOwned, bStamp } = this;
    // 1. rotate every plate a little further
    const Q = new Float64Array(9), M = new Float64Array(9);
    const alive = new Uint8Array(P);
    for (let i = 0; i < n; i++) alive[this.owner[i]] = 1;
    for (let q = 0; q < P; q++) {
      const wx = omega[3 * q], wy = omega[3 * q + 1], wz = omega[3 * q + 2];
      const w = Math.hypot(wx, wy, wz);
      if (!(w > 0)) continue;
      axisAngle(wx / w, wy / w, wz / w, w * h, Q);
      const Rq = rot.subarray(9 * q, 9 * q + 9);
      for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) M[3 * r + c] = Q[3 * r] * Rq[c] + Q[3 * r + 1] * Rq[3 + c] + Q[3 * r + 2] * Rq[6 + c];
      // re-orthonormalise (Gram-Schmidt on the rows) so round-off cannot shear the plate
      for (let r = 0; r < 3; r++) {
        for (let r2 = 0; r2 < r; r2++) {
          const d = M[3 * r] * M[3 * r2] + M[3 * r + 1] * M[3 * r2 + 1] + M[3 * r + 2] * M[3 * r2 + 2];
          for (let c = 0; c < 3; c++) M[3 * r + c] -= d * M[3 * r2 + c];
        }
        const L = Math.hypot(M[3 * r], M[3 * r + 1], M[3 * r + 2]);
        for (let c = 0; c < 3; c++) M[3 * r + c] /= L;
      }
      Rq.set(M);
    }
    // 2. look back along each plate's rotation: which plate cell is now at world cell i?
    //    A substep moves no point more than half a cell, so only the plates that owned i
    //    or one of its neighbours can have arrived there: those are the plates looked up
    //    (`near`, up to 7 per cell); bodyOf is kept current for them only.
    const stamp = ++this.stamp;
    this.prevOwner.set(this.owner);
    const nbr = mesh.neighbors;
    const near = this._near || (this._near = new Int32Array(7 * n));
    const nearCount = this._nearCount || (this._nearCount = new Uint8Array(n));
    for (let i = 0; i < n; i++) {
      let m = 0;
      const q0 = this.owner[i];
      near[7 * i] = q0; m = 1;
      for (let k = 0; k < 6; k++) {
        const j = nbr[6 * i + k];
        if (j < 0) continue;
        const q = this.owner[j];
        let seen = false;
        for (let c = 0; c < m; c++) if (near[7 * i + c] === q) { seen = true; break; }
        if (seen) continue;
        near[7 * i + m++] = q;
        // warm start from the neighbour, whose look-up for q is current
        bodyOf[q * n + i] = bodyOf[q * n + j];
      }
      nearCount[i] = m;
      for (let c = 0; c < m; c++) {
        const q = near[7 * i + c];
        if (!alive[q]) continue;
        const Rq = rot.subarray(9 * q, 9 * q + 9);
        const x = P3[3 * i], y = P3[3 * i + 1], z = P3[3 * i + 2];
        // body = R^T world
        const bx = Rq[0] * x + Rq[3] * y + Rq[6] * z;
        const by = Rq[1] * x + Rq[4] * y + Rq[7] * z;
        const bz = Rq[2] * x + Rq[5] * y + Rq[8] * z;
        const s = mesh.nearest(bx, by, bz, bodyOf[q * n + i]);
        bodyOf[q * n + i] = s;
        bStamp[q * n + s] = stamp;
      }
    }
    // 3. resolve each world cell: gap, single owner, or convergence
    const cand = new Int32Array(P);
    const A = mesh.area;
    const cr = this.params.crust;
    const last = this.last;
    const collide = new Map();
    for (let i = 0; i < n; i++) {
      let m = 0;
      for (let c = 0; c < nearCount[i]; c++) { const q = near[7 * i + c]; if (alive[q] && bOwned[q * n + bodyOf[q * n + i]]) cand[m++] = q; }
      if (m === 1) { this.owner[i] = cand[0]; continue; }
      if (m === 0) {
        // spreading: new oceanic crust joins the plate that was here
        let q = this.prevOwner[i];
        if (!alive[q]) for (let k = 0; k < 6; k++) { const j = mesh.neighbors[6 * i + k]; if (j >= 0 && alive[this.prevOwner[j]]) { q = this.prevOwner[j]; break; } }
        const s = bodyOf[q * n + i];
        const b = q * n + s;
        bOwned[b] = 1; this.bConsumed[b] = 0;
        this.bThick[b] = cr.oceanicThickness; this.bCont[b] = 0; this.bBirth[b] = this._birthTime(q, s, i, tNow); this.bSoil[b] = 0; this.bArc[b] = this.bAcc[b] = this.bSed[b] = 0;
        last.created += cr.oceanicThickness * A[s];
        this.budget.created += cr.oceanicThickness * A[s];
        this.owner[i] = q;
        continue;
      }
      // convergence: decide which plate overrides (see _overrider)
      let w = cand[0];
      for (let c = 1; c < m; c++) w = this._overrider(w, cand[c], i, tNow);
      for (let c = 0; c < m; c++) {
        const q = cand[c];
        if (q === w) continue;
        const sq = bodyOf[q * n + i];
        // a continent (not a sliver of arc crust) arriving at the margin: a collision
        const continent = this.bCont[q * n + sq] > 0.5;
        this._consume(q, sq, i, stamp);
        // also take plate-frame cells next to it that no world cell looked at this substep,
        // so subduction leaves no stray patches behind in the sinking plate
        for (let k = 0; k < 6; k++) {
          const t = mesh.neighbors[6 * sq + k];
          if (t >= 0 && bOwned[q * n + t] && bStamp[q * n + t] !== stamp) this._consume(q, t, i, stamp);
        }
        if (continent) {
          const key = Math.min(w, q) * P + Math.max(w, q);
          collide.set(key, (collide.get(key) || 0) + 1);
          this._noteEvent('collision', w, q, i);
        }
      }
      this.owner[i] = w;
    }
    // A world cell resolved early may look back to a plate cell that a later convergence
    // consumed (two world cells can sample the same plate cell): hand it to a plate that
    // still has crust there, or treat it as a gap.
    for (let i = 0; i < n; i++) {
      const q0 = this.owner[i];
      if (bOwned[q0 * n + bodyOf[q0 * n + i]]) continue;
      let q1 = -1, e1 = -Infinity;
      for (let c = 0; c < nearCount[i]; c++) {
        const q = near[7 * i + c];
        if (!alive[q]) continue;
        const b = q * n + bodyOf[q * n + i];
        if (!bOwned[b]) continue;
        const e = this.elevation(this.bThick[b], this.bCont[b], tNow - this.bBirth[b]);
        if (e > e1) { e1 = e; q1 = q; }
      }
      if (q1 >= 0) { this.owner[i] = q1; continue; }
      const b = q0 * n + bodyOf[q0 * n + i];
      const s = bodyOf[q0 * n + i];
      bOwned[b] = 1; this.bConsumed[b] = 0;
      this.bThick[b] = cr.oceanicThickness; this.bCont[b] = 0; this.bBirth[b] = this._birthTime(q0, s, i, tNow); this.bSoil[b] = 0; this.bArc[b] = this.bAcc[b] = this.bSed[b] = 0;
      last.created += cr.oceanicThickness * A[s];
      this.budget.created += cr.oceanicThickness * A[s];
    }
    // 4. close pinholes in the plate frames (cells skipped by the nearest-neighbour
    //    sampling while crust was being added around them)
    for (let q = 0; q < P; q++) {
      if (!alive[q]) continue;
      for (let s = 0; s < n; s++) {
        const b = q * n + s;
        if (bOwned[b] || this.bConsumed[b]) continue;
        let own = 0, th = 0, co = 0, bi = 0, so = 0, ar = 0, ac = 0, se = 0;
        const val = mesh.valence[s];
        for (let k = 0; k < val; k++) {
          const t = mesh.neighbors[6 * s + k];
          const bt = q * n + t;
          if (bOwned[bt]) { own++; th += this.bThick[bt]; co += this.bCont[bt]; bi += this.bBirth[bt]; so += this.bSoil[bt]; ar += this.bArc[bt]; ac += this.bAcc[bt]; se += this.bSed[bt]; }
        }
        if (own < val - 1) continue;
        bOwned[b] = 1;
        this.bThick[b] = th / own; this.bCont[b] = co / own; this.bBirth[b] = bi / own; this.bSoil[b] = so / own;
        this.bArc[b] = ar / own; this.bAcc[b] = ac / own; this.bSed[b] = se / own;
        last.created += this.bThick[b] * A[s];
        this.budget.created += this.bThick[b] * A[s];
      }
    }
    // 5. colliding continents lock together
    const pl = this.params.plates;
    if (pl.collisionLocking > 0 && collide.size) {
      const area = new Float64Array(P);
      for (let i = 0; i < n; i++) area[this.owner[i]] += A[i];
      const ref = Math.max(1, pl.lockingLength / mesh.spacing);
      for (const [key, cells] of collide) {
        const a = Math.floor(key / P), b = key % P;
        const f = Math.min(1, cells / ref);
        const keep = Math.exp((-h / pl.collisionLocking) * f);
        const wa = area[b] / (area[a] + area[b] || 1), wb = 1 - wa;
        for (let c = 0; c < 3; c++) {
          const d = (omega[3 * a + c] - omega[3 * b + c]) * (1 - keep);
          omega[3 * a + c] -= wa * d;
          omega[3 * b + c] += wb * d;
        }
        last.collisions += cells;
      }
    }
  }

  /**
   * Formation time of new crust in plate q's cell s, which fills the gap at world cell i.
   * The nearest-neighbour view opens gaps a whole cell at a time, so the new cell is
   * usually some way past the gap, or short of it; date it from when it passes:
   * signed distance beyond i along the plate's motion, over the plate's speed.
   */
  _birthTime(q, s, i, tNow) {
    const P3 = this.mesh.pos;
    const R = this.rot.subarray(9 * q, 9 * q + 9);
    const bx = P3[3 * s], by = P3[3 * s + 1], bz = P3[3 * s + 2];
    const wx = R[0] * bx + R[1] * by + R[2] * bz, wy = R[3] * bx + R[4] * by + R[5] * bz, wz = R[6] * bx + R[7] * by + R[8] * bz;
    const v = this.velocity(q, P3[3 * i], P3[3 * i + 1], P3[3 * i + 2]);
    const speed = Math.hypot(v[0], v[1], v[2]);
    if (!(speed > 0)) return tNow;
    const ahead = (this.params.radius * ((wx - P3[3 * i]) * v[0] + (wy - P3[3 * i + 1]) * v[1] + (wz - P3[3 * i + 2]) * v[2])) / speed;
    // behind the gap (not yet at it) counts as forming a little later
    return tNow - Math.max(-this.mesh.spacing, Math.min(ahead, this.mesh.spacing)) / speed;
  }

  /**
   * Which of plates a and b overrides the other at world cell i. Each pair keeps the
   * polarity it took at first contact (the more buoyant column then), because a
   * subduction zone does not flip from one cell to the next; but if the column meant to
   * sink stands higher than the overriding one by more than polarityFlip (a continent
   * arriving at the trench), the polarity flips, as in arc-continent collisions.
   */
  _overrider(a, b, i, tNow) {
    const n = this.mesh.size;
    const ba = a * n + this.bodyOf[a * n + i], bb = b * n + this.bodyOf[b * n + i];
    const ea = this.elevation(this.bThick[ba], this.bCont[ba], tNow - this.bBirth[ba]);
    const eb = this.elevation(this.bThick[bb], this.bCont[bb], tNow - this.bBirth[bb]);
    const key = Math.min(a, b) * this.plateCount + Math.max(a, b);
    let top = this.polarity.get(key);
    if (top === undefined) top = ea >= eb ? a : b;
    else if ((top === a ? eb - ea : ea - eb) > this.params.crust.polarityFlip) top = top === a ? b : a;
    this.polarity.set(key, top);
    return top;
  }

  /**
   * Plate q's frame cell sq goes under the overriding plate at world cell i: its
   * continental part is scraped off and its oceanic part subducts, feeding arc crust.
   * Both additions wait in `pending` until _spreadAdditions places them.
   */
  _consume(q, sq, i, stamp) {
    const n = this.mesh.size;
    const A = this.mesh.area;
    const b = q * n + sq;
    const vol = this.bThick[b] * A[sq];
    // a continent is too buoyant to sink and is scraped off whole; on sinking ocean floor
    // only part of the sediment it carries is scraped off (accretionary wedge), the rest
    // goes down with the slab
    const c = this.bCont[b];
    const vc = c * vol * (c > 0.5 ? 1 : this.params.crust.sedimentAccretion);
    const vo = vol - vc;
    const arc = this.params.crust.arcFraction * vo;
    this.pending[i] += vc + arc;
    this.pendingArc[i] += arc;
    this.bOwned[b] = 0;
    this.bConsumed[b] = 1;
    this.bStamp[b] = stamp;
    this.last.subducted += vo;
    this.last.arc += arc;
    this.last.accreted += vc;
    this.budget.subducted += vo;
    this.budget.arc += arc;
    this.budget.accreted += vc;
  }

  _maybeReorganize(dt) {
    const every = this.params.plates.reorganizeEvery;
    if (!(every > 0)) return;
    if (Math.floor((this.time + dt) / every) > Math.floor(this.time / every)) {
      for (let q = 0; q < this.plateCount; q++) this._randomMotion(q);
    }
  }

  /** Copy the plate-frame state into the world view. */
  _buildView() {
    const n = this.mesh.size;
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      const b = q * n + this.bodyOf[q * n + i];
      this.thickness[i] = this.bThick[b];
      this.continental[i] = this.bCont[b];
      this.age[i] = Math.max(0, this.time - this.bBirth[b]);
      this.arcShare[i] = this.bArc[b];
      this.accretedShare[i] = this.bAcc[b];
      this.sedimentShare[i] = this.bSed[b];
    }
  }

  /**
   * Place this step's scraped-off and arc crust: spread it over the overriding plate
   * with a diffusion of width ~arcWidth (arcs stand a few hundred km behind trenches;
   * collision belts are as wide), then add it to the plate frames as continental crust.
   */
  _spreadAdditions() {
    const mesh = this.mesh;
    const n = mesh.size;
    const { neighbors: nb, dist, face, area } = mesh;
    const vol = this.pending, volArc = this.pendingArc;
    let any = false;
    for (let i = 0; i < n; i++) if (vol[i] > 0) { any = true; break; }
    if (!any) return;
    let worst = 0;
    for (let i = 0; i < n; i++) {
      let t = 0;
      for (let k = 0; k < 6; k++) { const e = 6 * i + k; if (nb[e] >= 0) t += face[e] / dist[e]; }
      worst = Math.max(worst, t / area[i]);
    }
    // each pass adds about 0.45 / worst x 2 to the variance (m^2); worst ~ 4 / spacing^2
    const w = this.params.crust.arcWidth;
    const passes = Math.min(200, Math.ceil((w * w * worst) / 0.9));
    const kdt = 0.45 / worst;
    // only cells within 3 widths of a source take part (beyond, the spread is negligible);
    // edges leading out of that set carry nothing, so volume is still conserved exactly
    const rings = Math.ceil((3 * w) / mesh.spacing) + 1;
    const inSet = new Uint8Array(n);
    let set = [];
    for (let i = 0; i < n; i++) if (vol[i] > 0) { inSet[i] = 1; set.push(i); }
    for (let r = 0, front = set.slice(); r < rings; r++) {
      const next = [];
      for (const i of front) for (let k = 0; k < 6; k++) { const j = nb[6 * i + k]; if (j >= 0 && !inSet[j] && this.owner[j] === this.owner[i]) { inSet[j] = 1; next.push(j); } }
      set = set.concat(next);
      front = next;
    }
    set.sort((a, b) => a - b);
    const u = new Float64Array(n), ua = new Float64Array(n);
    for (let it = 0; it < passes; it++) {
      for (const i of set) { u[i] = vol[i] / area[i]; ua[i] = volArc[i] / area[i]; }
      for (const i of set) {
        for (let k = 0; k < 6; k++) {
          const e = 6 * i + k;
          const j = nb[e];
          if (j <= i || !inSet[j] || this.owner[j] !== this.owner[i]) continue;
          const g = (kdt * face[e]) / dist[e];
          const f = g * (u[i] - u[j]), fa = g * (ua[i] - ua[j]);
          vol[i] -= f; vol[j] += f;
          volArc[i] -= fa; volArc[j] += fa;
        }
      }
    }
    for (const i of set) {
      if (!(vol[i] > 0)) continue;
      const q = this.owner[i];
      const s = this.bodyOf[q * n + i];
      const b = q * n + s;
      const add = vol[i] / area[s];
      const addArc = Math.min(add, Math.max(0, volArc[i] / area[s]));
      const h0 = this.bThick[b];
      this.bThick[b] = h0 + add;
      this.bCont[b] = (this.bCont[b] * h0 + add) / (h0 + add);
      this.bArc[b] = (this.bArc[b] * h0 + addArc) / (h0 + add);
      this.bAcc[b] = (this.bAcc[b] * h0 + add - addArc) / (h0 + add);
      this.bSed[b] = (this.bSed[b] * h0) / (h0 + add);
    }
    vol.fill(0);
    volArc.fill(0);
  }

  /** Crust thicker than flowThreshold spreads to neighbours on the same plate. */
  _crustalFlow(dt) {
    const cr = this.params.crust;
    const D = cr.flowDiffusivity;
    if (!(D > 0)) return;
    const mesh = this.mesh;
    const n = mesh.size;
    const { neighbors: nb, dist, face, area } = mesh;
    const h0 = cr.flowThreshold;
    const th = Float64Array.from(this.thickness);
    let worst = 0;
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let k = 0; k < 6; k++) { const e = 6 * i + k; if (nb[e] >= 0) s += face[e] / dist[e]; }
      worst = Math.max(worst, s / area[i]);
    }
    // thicker crust is hotter and weaker at depth: it flows faster, as (1 + (excess / e0)^2)
    const e0 = cr.flowScale;
    let emax = 0;
    for (let i = 0; i < n; i++) emax = Math.max(emax, th[i] - h0);
    let fmax = 1 + (Math.max(0, emax) / e0) ** 2;
    const cap = Math.max(1, cr.maxFlowSubsteps);
    // when the substep cap binds, cap the flow law's speed-up so each substep stays stable
    if ((dt * D * fmax * worst) / 0.5 > cap) fmax = Math.max(1, (0.5 * cap) / (dt * D * worst));
    const subs = Math.min(cap, Math.max(1, Math.ceil((dt * D * fmax * worst) / 0.5)));
    const hstep = dt / subs;
    const out = new Float64Array(n);
    const flux = new Float64Array(6 * n);
    for (let it = 0; it < subs; it++) {
      out.fill(0);
      for (let i = 0; i < n; i++) {
        const ei = th[i] - h0;
        for (let k = 0; k < 6; k++) {
          const e = 6 * i + k;
          const j = nb[e];
          flux[e] = 0;
          if (j <= i || this.owner[j] !== this.owner[i]) continue;
          const ej = th[j] - h0;
          if (ei <= 0 && ej <= 0) continue;
          const em = Math.max(ei, ej);
          const q = ((D * Math.min(fmax, 1 + (em / e0) ** 2) * (Math.max(ei, 0) - Math.max(ej, 0))) / dist[e]) * face[e] * hstep; // volume i -> j
          flux[e] = q;
          out[q > 0 ? i : j] += Math.abs(q);
        }
      }
      for (let i = 0; i < n; i++) {
        const avail = Math.max(0, th[i] - h0) * area[i];
        out[i] = out[i] > avail ? avail / out[i] : 1;
      }
      for (let i = 0; i < n; i++) {
        for (let k = 0; k < 6; k++) {
          const e = 6 * i + k;
          const q0 = flux[e];
          if (q0 === 0) continue;
          const j = nb[e];
          const q = q0 > 0 ? q0 * out[i] : q0 * out[j];
          th[i] -= q / area[i];
          th[j] += q / area[j];
        }
      }
    }
    // crust beyond maxThickness is shed: it joins the crust waiting to be spread outward
    const hmax = cr.maxThickness;
    let shed = false;
    if (hmax > 0) {
      for (let i = 0; i < n; i++) {
        if (th[i] <= hmax) continue;
        this.pending[i] += (th[i] - hmax) * area[i];
        th[i] = hmax;
        shed = true;
      }
    }
    // write the change back to the plate frames
    for (let i = 0; i < n; i++) {
      const d = th[i] - this.thickness[i];
      if (d === 0) continue;
      const q = this.owner[i];
      const s = this.bodyOf[q * n + i];
      this.bThick[q * n + s] += (d * area[i]) / area[s];
    }
    if (shed) this._spreadAdditions();
  }

  /**
   * Cheap stand-in for erosion at the planet scale (used when surfaceProcesses is off):
   * mean denudation rate proportional to height above a floor a few hundred metres above
   * sea level (Ahnert 1970, Am. J. Sci. 268, found it proportional to relief, and the
   * relief of old low continents is small; the floor keeps them from all being planed
   * down to the sea). Isostatic rebound returns part of every metre
   * removed, so the crust loses dz / (1 - rho_c / rho_m) to lower the surface by dz. The
   * sediment is laid down in the nearest sea cell (by steps across the mesh) as
   * continental crust: margins build out into shelves and sedimentary basins, and the
   * crust budget still closes.
   */
  _denude(dt) {
    const { timescale: tau, floor } = this.params.denudation;
    if (!(tau > 0)) return;
    const cr = this.params.crust;
    const sim = this.sim;
    const n = this.mesh.size;
    const f = 1 - Math.exp(-dt / tau);
    const A = this.mesh.area, nb = this.mesh.neighbors;
    // nearest sea cell of every cell: breadth-first from the sea
    const target = this._seaTarget || (this._seaTarget = new Int32Array(n));
    target.fill(-1);
    const queue = new Int32Array(n);
    let head = 0, tail = 0;
    for (let i = 0; i < n; i++) if (this.zTectonic[i] < this.seaLevel) { target[i] = i; queue[tail++] = i; }
    while (head < tail) {
      const i = queue[head++];
      for (let k = 0; k < 6; k++) { const j = nb[6 * i + k]; if (j >= 0 && target[j] < 0) { target[j] = target[i]; queue[tail++] = j; } }
    }
    const sed = this._sed || (this._sed = new Float64Array(n));
    sed.fill(0);
    for (let i = 0; i < n; i++) {
      const above = this.zTectonic[i] - this.seaLevel - floor;
      if (!(above > 0) || target[i] < 0) continue;
      const rho = this.continental[i] * cr.continentalDensity + (1 - this.continental[i]) * cr.oceanicDensity;
      const dh = Math.min(this.thickness[i] - 1000, (f * above) / (1 - rho / cr.mantleDensity));
      if (!(dh > 0)) continue;
      sim.rock[i] -= dh; sim.z[i] -= dh; // the write-back takes the thickness change from sim.z
      sed[target[i]] += dh * A[i];
    }
    let laid = 0;
    for (let i = 0; i < n; i++) {
      if (!(sed[i] > 0)) continue;
      const q = this.owner[i];
      const s = this.bodyOf[q * n + i];
      const b = q * n + s;
      const add = sed[i] / A[s];
      const h0 = this.bThick[b];
      this.bThick[b] = h0 + add;
      this.bCont[b] = (this.bCont[b] * h0 + add) / (h0 + add);
      this.bSed[b] = (this.bSed[b] * h0 + add) / (h0 + add);
      this.bArc[b] = (this.bArc[b] * h0) / (h0 + add);
      this.bAcc[b] = (this.bAcc[b] * h0) / (h0 + add);
      laid += sed[i];
    }
    this.budget.surface += laid;
  }

  /** Record that world cell i took part in an event of `kind` between plates a and b this step. */
  _noteEvent(kind, a, b, i) {
    if (!this.params.events.log) return;
    const key = `${kind}:${Math.min(a, b)}:${Math.max(a, b)}`;
    if (!this._stepCells) this._stepCells = new Map();
    let list = this._stepCells.get(key);
    if (!list) this._stepCells.set(key, (list = []));
    list.push(i);
  }

  /** Continental rifting: a divergent boundary with continental crust on both sides. */
  _noteRifts() {
    if (!this.params.events.log) return;
    const nb = this.mesh.neighbors;
    for (let i = 0; i < this.mesh.size; i++) {
      if (this.boundary[i] !== 2 || this.continental[i] < 0.5) continue;
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j >= 0 && this.owner[j] !== this.owner[i] && this.continental[j] >= 0.5) { this._noteEvent('rifting', this.owner[i], this.owner[j], i); break; }
      }
    }
  }

  /**
   * Turn this step's event cells into dated events. An event that continues from the
   * previous step is extended (its region grows to cover both); one that did not occur
   * this step is closed.
   */
  _closeStepEvents(t0, t1) {
    if (!this.params.events.log) return;
    const P3 = this.mesh.pos;
    const step = this._stepCells || new Map();
    const min = this.params.events.minCells;
    const gap = this.params.events.gap;
    for (const [key, ev] of this._openEvents) if (t1 - ev.end > gap && !(step.has(key) && step.get(key).length >= min)) { this.events.push(ev); this._openEvents.delete(key); }
    for (const [key, cells] of step) {
      if (cells.length < min) continue;
      let ev = this._openEvents.get(key);
      const [kind, a, b] = key.split(':');
      if (!ev) {
        ev = { kind, start: t0, end: t1, plates: [+a, +b], cells: 0, sum: [0, 0, 0], members: new Set() };
        this._openEvents.set(key, ev);
      }
      ev.end = t1;
      for (const i of cells) {
        if (ev.members.has(i)) continue;
        ev.members.add(i);
        ev.sum[0] += P3[3 * i]; ev.sum[1] += P3[3 * i + 1]; ev.sum[2] += P3[3 * i + 2];
      }
      ev.cells = ev.members.size;
    }
    this._stepCells = new Map();
  }

  /**
   * All events so far, open ones included, with their regions as spherical caps:
   * {kind, start, end (model years), center (unit vector), radius (rad), plates, cells}.
   */
  eventLog() {
    const P3 = this.mesh.pos;
    const out = [];
    for (const ev of [...this.events, ...this._openEvents.values()]) {
      const L = Math.hypot(...ev.sum) || 1;
      const c = ev.sum.map((v) => v / L);
      let r = 0;
      for (const i of ev.members) r = Math.max(r, Math.acos(Math.min(1, c[0] * P3[3 * i] + c[1] * P3[3 * i + 1] + c[2] * P3[3 * i + 2])));
      out.push({ kind: ev.kind, start: ev.start, end: ev.end, center: c, radius: r + (0.5 * this.mesh.spacing) / this.params.radius, plates: ev.plates, cells: ev.cells });
    }
    return out.sort((x, y) => x.start - y.start || (x.kind < y.kind ? -1 : x.kind > y.kind ? 1 : 0) || x.plates[0] - y.plates[0] || x.plates[1] - y.plates[1]);
  }

  /** Material class of every world cell (see MATERIAL). */
  material(out = new Uint8Array(this.mesh.size)) {
    const cr = this.params.crust;
    for (let i = 0; i < this.mesh.size; i++) {
      const c = this.continental[i];
      if (c < 0.5) { out[i] = MATERIAL.oceanic; continue; }
      const arc = this.arcShare[i], acc = this.accretedShare[i];
      if (arc > 0.3 && arc >= acc) out[i] = MATERIAL.arc;
      else if (acc > 0.3) out[i] = MATERIAL.terrane;
      else if (this.sedimentShare[i] > 0.3 || this.zTectonic[i] < this.seaLevel) out[i] = MATERIAL.basin;
      else out[i] = MATERIAL.craton;
    }
    return out;
  }

  _classifyBoundaries() {
    const mesh = this.mesh;
    const n = mesh.size;
    const P3 = mesh.pos;
    const vi = [0, 0, 0], vj = [0, 0, 0];
    this.boundary.fill(0);
    for (let i = 0; i < n; i++) {
      const a = this.owner[i];
      let bestRate = 0, type = 0;
      for (let k = 0; k < 6; k++) {
        const j = mesh.neighbors[6 * i + k];
        if (j < 0 || this.owner[j] === a) continue;
        const x = P3[3 * i], y = P3[3 * i + 1], z = P3[3 * i + 2];
        this.velocity(a, x, y, z, vi);
        this.velocity(this.owner[j], x, y, z, vj);
        // unit direction from i towards j, in the tangent plane
        let dx = P3[3 * j] - x, dy = P3[3 * j + 1] - y, dz = P3[3 * j + 2] - z;
        const L = Math.hypot(dx, dy, dz);
        dx /= L; dy /= L; dz /= L;
        const rx = vi[0] - vj[0], ry = vi[1] - vj[1], rz = vi[2] - vj[2];
        const closing = rx * dx + ry * dy + rz * dz; // > 0: plates approach
        const rate = Math.hypot(rx, ry, rz);
        if (rate <= bestRate) continue;
        bestRate = rate;
        type = Math.abs(closing) < 0.5 * rate ? 3 : closing > 0 ? 1 : 2;
      }
      this.boundary[i] = type;
    }
  }

  /** Elevations from isostasy (plus trenches) into zTectonic. */
  _isostasy() {
    const n = this.mesh.size;
    const zT = this.zTectonic;
    for (let i = 0; i < n; i++) zT[i] = this.elevation(this.thickness[i], this.continental[i], this.age[i]);
    const td = this.params.crust.trenchDepth;
    if (td > 0) {
      // trenches: sea floor at the edge of a plate sinking under its neighbour
      const nb = this.mesh.neighbors;
      const P3 = this.mesh.pos;
      const P = this.plateCount;
      const vi = [0, 0, 0], vj = [0, 0, 0];
      for (let i = 0; i < n; i++) {
        if (zT[i] > 0) continue;
        const a = this.owner[i];
        for (let k = 0; k < 6; k++) {
          const j = nb[6 * i + k];
          if (j < 0) continue;
          const b = this.owner[j];
          if (b === a || this.polarity.get(Math.min(a, b) * P + Math.max(a, b)) !== b) continue;
          const x = P3[3 * i], y = P3[3 * i + 1], z = P3[3 * i + 2];
          this.velocity(a, x, y, z, vi);
          this.velocity(b, x, y, z, vj);
          const closing = (vi[0] - vj[0]) * (P3[3 * j] - x) + (vi[1] - vj[1]) * (P3[3 * j + 1] - y) + (vi[2] - vj[2]) * (P3[3 * j + 2] - z);
          if (closing > 0) { zT[i] -= td; break; }
        }
      }
    }
  }

  _waterBelow(level) {
    const z = this.zTectonic;
    const A = this.mesh.area;
    let v = 0;
    for (let i = 0; i < this.mesh.size; i++) if (z[i] < level) v += (level - z[i]) * A[i];
    return v;
  }

  /** Sea level that holds `volume` of water over the current topography (bisection). */
  _seaLevelForVolume(volume) {
    let lo = -12000, hi = 12000;
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      if (this._waterBelow(mid) < volume) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /** Put the tectonic surface into the surface simulation. */
  _applyToSurface() {
    const n = this.mesh.size;
    const sim = this.sim;
    for (let i = 0; i < n; i++) {
      const q = this.owner[i];
      const soil = Math.min(this.bSoil[q * n + this.bodyOf[q * n + i]], this.thickness[i]);
      sim.soil[i] = soil;
      sim.rock[i] = this.zTectonic[i] - soil;
    }
    sim.params.boundary.seaLevel = this.seaLevel; // the planet owns sea level
    sim.seaLevel = this.seaLevel;
    if (this.params.surfaceProcesses) sim._updateDerivedState();
    else for (let i = 0; i < n; i++) sim.z[i] = sim.rock[i] + sim.soil[i]; // nothing else is read without surface processes
  }
}
