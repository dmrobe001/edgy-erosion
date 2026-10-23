// Parameters: defaults, presets, and merging.
//
// Units throughout: metres (m), years (yr), pascals (Pa), kilograms (kg).
// Parameters are plain JSON so that a run is reproducible from its parameter
// object alone (see Simulation.snapshot()).

/**
 * Rock library. Values are order-of-magnitude, chosen to be plausible at the
 * landscape scale rather than for any particular formation.
 *
 *  erodibility    K in the stream power law  E = K * A^m * S^n  (for m = 0.5, n = 1:
 *                 units yr^-1). Larger = channels cut it faster.
 *  cohesion       rock-MASS cohesion (Pa) in the Culmann cliff-stability criterion.
 *                 Much lower than intact-sample strength because of joints.
 *  frictionAngle  internal friction angle (degrees). Cliffs no steeper than this
 *                 never fail, whatever their height.
 *  weathering     maximum soil (regolith) production rate on bare rock (m/yr).
 *  color          for display only.
 */
export const ROCKS = {
  granite:   { color: '#b9a99a', erodibility: 1e-5, cohesion: 1.0e6, frictionAngle: 45, weathering: 5e-5 },
  basalt:    { color: '#5d5a5c', erodibility: 1e-5, cohesion: 8.0e5, frictionAngle: 42, weathering: 6e-5 },
  sandstone: { color: '#c98a5a', erodibility: 2e-5, cohesion: 4.0e5, frictionAngle: 40, weathering: 1e-4 },
  limestone: { color: '#d8d2bd', erodibility: 1.5e-5, cohesion: 5.0e5, frictionAngle: 40, weathering: 8e-5 },
  schist:    { color: '#7f8a86', erodibility: 3e-5, cohesion: 2.0e5, frictionAngle: 35, weathering: 1.5e-4 },
  siltstone: { color: '#a8846a', erodibility: 8e-5, cohesion: 6.0e4, frictionAngle: 31, weathering: 3e-4 },
  shale:     { color: '#8f6f6a', erodibility: 1.5e-4, cohesion: 3.0e4, frictionAngle: 28, weathering: 4e-4 },
};

export const DEFAULTS = {
  seed: 1,
  grid: {
    cols: 128,      // cells per row
    rows: 128,      // number of rows
    dx: 30,         // centre-to-centre lattice spacing (m)
  },
  time: {
    dt: 2000,       // years per step
  },
  boundary: {
    // Which grid edges are fixed-elevation outlets (base level). Others are closed (no flux).
    open: ['north', 'south', 'east', 'west'],
    // Cells below this elevation are sea: water and sediment sinks. null = no sea.
    seaLevel: null,
  },
  initial: {
    // 'flat' | 'plateau' | 'dome' | 'ramp' | 'ridge'
    type: 'plateau',
    base: 0,         // elevation of the surroundings (m)
    height: 300,     // relief of the initial landform (m)
    margin: 0.15,    // plateau: fraction of the domain left low around the edge
    edgeWidth: 0.05, // plateau: width of the initial slope, as a fraction of domain width
    edgeNoise: 0,    // plateau: irregularity of the outline, as a fraction of domain width
    noise: { amplitude: 15, wavelength: 1000, octaves: 5, gain: 0.5 },
    soil: 0.5,       // initial soil (regolith) thickness (m)
  },
  uplift: {
    rate: 1e-4,      // rock uplift rate relative to base level (m/yr)
    // 'uniform' | 'dome' | 'tilt' | 'fault'
    pattern: 'uniform',
    // dome: gaussian of this width (fraction of domain width)
    radius: 0.3,
    // fault: line through the domain centre at this angle (degrees from +x); the
    // side to the left of the direction of strike is uplifted. offset shifts the line
    // (fraction of domain width).
    angle: 0,
    offset: 0,
    // tilt: rate goes from 0 at one side to `rate` at the other, along `angle`.
  },
  rocks: ROCKS,
  strata: {
    // Layers from the top of the stratigraphic column downwards.
    layers: [{ rock: 'sandstone', thickness: 1e9 }],
    top: 1000,          // elevation of the top of the column at the domain centre (m)
    repeat: false,      // repeat the layer sequence indefinitely downwards
    basement: 'granite',// rock below the last layer (when not repeating)
    dip: 0,             // dip angle of the layers (degrees)
    dipDirection: 0,    // azimuth the layers descend towards (degrees, 0 = east, 90 = north)
    undulation: { amplitude: 0, wavelength: 2000 }, // gentle folding/irregularity of layer surfaces (m)
    // Random spatial variation of rock strength: erodibility x f, cohesion / f, with
    // ln f a smooth noise field of this amplitude (natural-log units) and wavelength (m).
    strengthVariation: { amplitude: 0.5, wavelength: 300 },
  },
  // Vertical bodies of rock that cut across the strata: volcanic plugs, intrusions, dikes.
  // Positions and sizes are fractions of the domain width/height so presets survive
  // changes of grid size. `top` (optional) is the elevation the body reaches up to.
  // e.g. { shape: 'disc', x: 0.5, y: 0.5, radius: 0.08, rock: 'basalt', top: 250 }
  //      { shape: 'dike', x0: 0.1, y0: 0.2, x1: 0.9, y1: 0.7, width: 0.01, rock: 'basalt' }
  bodies: [],
  climate: {
    precipitation: 1.0,          // mean effective runoff-producing precipitation (m/yr)
    precipitationGradient: 0,    // fractional change from west edge to east edge (e.g. -0.5)
    frost: 0.3,                  // 0..1, freeze-thaw intensity: speeds rockfall and weathering
  },
  vegetation: {
    // 'auto' derives the maximum cover from precipitation; a number 0..1 fixes it.
    cover: 'auto',
    halfSaturation: 0.5,         // precipitation (m/yr) giving 50% cover when 'auto'
    soilDepthScale: 0.3,         // cover falls off on soils thinner than this (m): bare rock is bare
    erodibilityReduction: 2.0,   // soil erodibility multiplied by exp(-this * cover)
    diffusivityBoost: 1.0,       // creep diffusivity multiplied by (1 + this * cover) (bioturbation)
    slopeBoost: 0.3,             // critical soil slope multiplied by (1 + this * cover) (root cohesion)
    productionBoost: 0.5,        // soil production multiplied by (1 + this * cover)
  },
  fluvial: {
    enabled: true,
    m: 0.5,                      // drainage-area exponent
    n: 1.0,                      // slope exponent
    soilErodibility: 4e-4,       // K for loose soil/sediment (same units as rock K)
    soilCoverScale: 0.5,         // soil thickness (m) that shields bedrock from incision (e-folding)
    // 'mfd': drainage area spread over all lower neighbours (smoother hillslopes, no
    // grid-aligned streaks); 'sfd': all to the steepest neighbour (classic D6).
    routing: 'mfd',
    flowExponent: 4,             // MFD: share proportional to drop^this (higher = more channelised)
    // Pick each cell's receiver at random among lower neighbours (weights drop^randomExponent),
    // re-drawn every step. Breaks up channels that would otherwise lock onto the lattice axes.
    randomReceivers: true,
    randomExponent: 2,
    deposition: 0.5,             // dimensionless deposition coefficient G; 0 = detachment-limited
  },
  hillslope: {
    enabled: true,
    diffusivity: 0.01,           // soil creep diffusivity (m^2/yr)
    criticalSlope: 0.8,          // gradient where creep diverges; also the soil angle of repose (~39 deg)
    maxNonlinearity: 10,         // cap on the nonlinear enhancement factor 1/(1-(S/Sc)^2)
    transportDepth: 0.3,         // creep falls off on soils thinner than this (m)
    avalancheIterations: 40,     // max relaxation sweeps per step for over-steep soil/talus
  },
  weathering: {
    enabled: true,
    productionDepth: 0.5,        // soil production e-folding depth h* (m)
    frostBoost: 1.0,             // production multiplied by (1 + this * frost)
    moistureExponent: 1.0,       // production multiplied by (precipitation / 1 m/yr)^this
  },
  rockfall: {
    enabled: true,
    rockDensity: 2600,           // kg/m^3
    timescale: 300,              // mean waiting time (yr) for an unstable cliff to fail
    frostFactor: 5,              // failure rate multiplied by (1 + this * frost)
  },
};

/**
 * Presets: partial parameter objects layered over DEFAULTS. Each aims at a
 * recognisable landscape. Spatial positions are fractions of the domain, so they
 * work for any grid size, but process lengths (layer thickness, relief) are in
 * metres, so the look depends on dx.
 */
export const PRESETS = {
  'mesa-country': {
    description: 'Flat-lying hard sandstone over soft shale in a dry climate: stepped escarpments, mesas, buttes and talus.',
    time: { dt: 5000 },
    initial: { type: 'plateau', height: 350, margin: 0.06, edgeWidth: 0.04, edgeNoise: 0.08, noise: { amplitude: 25, wavelength: 1500 }, soil: 0.2 },
    uplift: { rate: 3e-4 },
    strata: {
      layers: [
        { rock: 'sandstone', thickness: 40 },
        { rock: 'shale', thickness: 90 },
        { rock: 'sandstone', thickness: 25 },
        { rock: 'shale', thickness: 70 },
      ],
      top: 345, repeat: true, undulation: { amplitude: 5, wavelength: 3000 },
    },
    climate: { precipitation: 0.3, frost: 0.4 },
    vegetation: { halfSaturation: 1.0 },
    fluvial: { deposition: 0 },
    hillslope: { diffusivity: 0.003 },
  },
  'alpine': {
    description: 'Fast uplift of strong granite under a cold, wet climate: sharp ridges, steep valley walls, frequent rockfall.',
    time: { dt: 2000 },
    initial: { type: 'dome', height: 300, noise: { amplitude: 120, wavelength: 2500 }, soil: 0.1 },
    uplift: { rate: 5e-4, pattern: 'uniform' },
    strata: { layers: [{ rock: 'granite', thickness: 1e9 }], undulation: { amplitude: 0 } },
    bodies: [],
    climate: { precipitation: 2.0, frost: 0.9 },
    vegetation: { cover: 0.3 },
    hillslope: { diffusivity: 0.005 },
  },
  'humid-hills': {
    description: 'Slow uplift of deeply weathered schist under a warm, wet, forested climate: rounded soil-mantled hills.',
    time: { dt: 5000 },
    initial: { type: 'flat', height: 0, base: 50, noise: { amplitude: 60, wavelength: 2000 }, soil: 2 },
    uplift: { rate: 2e-4 },
    strata: { layers: [{ rock: 'schist', thickness: 1e9 }] },
    climate: { precipitation: 1.8, frost: 0 },
    vegetation: { cover: 'auto' },
    hillslope: { diffusivity: 0.02 },
    weathering: { productionDepth: 0.8 },
  },
  'cuestas': {
    description: 'Gently dipping alternating hard and soft layers: asymmetric ridges (cuestas) with a steep scarp and a gentle dip slope.',
    time: { dt: 4000 },
    initial: { type: 'ramp', height: 200, angle: 0, noise: { amplitude: 10 }, soil: 0.3 },
    uplift: { rate: 2e-4, pattern: 'uniform' },
    boundary: { open: ['east', 'west'] },
    strata: {
      layers: [
        { rock: 'limestone', thickness: 35 },
        { rock: 'shale', thickness: 80 },
        { rock: 'sandstone', thickness: 30 },
        { rock: 'siltstone', thickness: 60 },
      ],
      top: 260, repeat: true, dip: 4, dipDirection: 0,
    },
    climate: { precipitation: 0.7, frost: 0.3 },
  },
  'fault-scarp': {
    description: 'A normal-fault block rising out of the plain: a fresh fault scarp degrading into triangular facets and fans.',
    time: { dt: 2000 },
    initial: { type: 'flat', height: 0, base: 0, noise: { amplitude: 5, wavelength: 800 }, soil: 0.5 },
    uplift: { rate: 8e-4, pattern: 'fault', angle: 10, offset: 0.05 },
    boundary: { open: ['south'] },
    strata: { layers: [{ rock: 'limestone', thickness: 1e9 }] },
    climate: { precipitation: 0.4, frost: 0.5 },
    vegetation: { halfSaturation: 0.8 },
  },
  'badlands': {
    description: 'Weak shale and siltstone, no vegetation, dry climate with intense storms: dense gullies and sharp knife-edge ridges.',
    time: { dt: 1000 },
    initial: { type: 'plateau', height: 80, margin: 0.1, noise: { amplitude: 4, wavelength: 300 }, soil: 0.05 },
    uplift: { rate: 5e-4 },
    strata: {
      layers: [{ rock: 'siltstone', thickness: 6 }, { rock: 'shale', thickness: 10 }],
      top: 80, repeat: true,
    },
    climate: { precipitation: 0.25, frost: 0.2 },
    vegetation: { cover: 0 },
    hillslope: { diffusivity: 0.002 },
    rockfall: { timescale: 100 },
  },
  'volcanic-necks': {
    description: 'Soft sediments cut by basalt dikes and plugs: as the land wears down the resistant bodies stand out as walls and towers.',
    time: { dt: 4000 },
    initial: { type: 'flat', base: 150, height: 0, noise: { amplitude: 20, wavelength: 1500 }, soil: 0.5 },
    uplift: { rate: 3e-4 },
    strata: { layers: [{ rock: 'siltstone', thickness: 30 }, { rock: 'shale', thickness: 50 }], top: 200, repeat: true },
    bodies: [
      { shape: 'disc', x: 0.38, y: 0.42, radius: 0.045, rock: 'basalt' },
      { shape: 'disc', x: 0.66, y: 0.62, radius: 0.03, rock: 'basalt', top: 120 },
      { shape: 'dike', x0: 0.38, y0: 0.42, x1: 0.85, y1: 0.12, width: 0.012, rock: 'basalt' },
      { shape: 'dike', x0: 0.38, y0: 0.42, x1: 0.12, y1: 0.9, width: 0.01, rock: 'basalt' },
    ],
    climate: { precipitation: 0.4, frost: 0.4 },
  },
  'coastal-plateau': {
    description: 'A sandstone-capped plateau rising from the sea: sea cliffs, deltas and drowned valleys.',
    time: { dt: 2000 },
    initial: { type: 'ramp', base: -60, height: 260, angle: 180, noise: { amplitude: 15 }, soil: 0.3 },
    uplift: { rate: 3e-5 },
    boundary: { open: [], seaLevel: 0 },
    strata: {
      layers: [{ rock: 'sandstone', thickness: 30 }, { rock: 'siltstone', thickness: 50 }],
      top: 220, repeat: true, dip: 1, dipDirection: 180,
    },
    climate: { precipitation: 1.0, frost: 0.2 },
  },
};

const isPlainObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

/** Recursive merge: objects merge key by key, everything else (including arrays) replaces. */
export function deepMerge(base, over) {
  if (!isPlainObject(over)) return over === undefined ? clone(base) : clone(over);
  const out = isPlainObject(base) ? clone(base) : {};
  for (const [k, v] of Object.entries(over)) {
    out[k] = isPlainObject(v) && isPlainObject(out[k]) ? deepMerge(out[k], v) : clone(v);
  }
  return out;
}

function clone(v) {
  return v === undefined ? undefined : JSON.parse(JSON.stringify(v));
}

/**
 * Build a complete parameter object: DEFAULTS <- PRESETS[preset] <- overrides.
 * The returned object records which preset it came from but is self-contained.
 *
 * @param {object|string} [input] overrides (may include `preset: name`), or a preset name
 */
export function resolveParams(input = {}) {
  const overrides = typeof input === 'string' ? { preset: input } : input;
  let p = clone(DEFAULTS);
  if (overrides.preset) {
    const preset = PRESETS[overrides.preset];
    if (!preset) throw new Error(`unknown preset "${overrides.preset}"; try one of: ${Object.keys(PRESETS).join(', ')}`);
    const { description, ...rest } = preset;
    p = deepMerge(p, rest);
  }
  p = deepMerge(p, overrides);
  validate(p);
  return p;
}

/** Throw a descriptive error for parameter combinations that cannot work. */
export function validate(p) {
  const fail = (msg) => { throw new Error(`invalid parameters: ${msg}`); };
  if (!(p.grid.cols >= 3 && p.grid.rows >= 3)) fail('grid must be at least 3 x 3');
  if (!(p.grid.dx > 0)) fail('grid.dx must be positive');
  if (!(p.time.dt > 0)) fail('time.dt must be positive');
  if (!(p.fluvial.n > 0)) fail('fluvial.n must be positive');
  const rockNames = new Set(Object.keys(p.rocks));
  const needRock = (name, where) => { if (!rockNames.has(name)) fail(`${where} refers to unknown rock "${name}"`); };
  if (!p.strata.layers.length) fail('strata.layers must not be empty');
  p.strata.layers.forEach((l, k) => {
    needRock(l.rock, `strata.layers[${k}]`);
    if (!(l.thickness > 0)) fail(`strata.layers[${k}].thickness must be positive`);
  });
  needRock(p.strata.basement, 'strata.basement');
  p.bodies.forEach((b, k) => needRock(b.rock, `bodies[${k}]`));
  const sides = ['north', 'south', 'east', 'west'];
  for (const s of p.boundary.open) if (!sides.includes(s)) fail(`boundary.open: unknown side "${s}"`);
  if (!p.boundary.open.length && p.boundary.seaLevel == null) {
    fail('need at least one open boundary or a sea level, or water has nowhere to go');
  }
}
