// Public API. Works unchanged in browsers (as an ES module) and in Node >= 18.
export { Simulation } from './simulation.js';
export { HexGrid } from './hexgrid.js';
export { Lithology } from './lithology.js';
export { DEFAULTS, PRESETS, ROCKS, resolveParams, deepMerge } from './params.js';
export { routeFlow } from './flow.js';
export { culmannHeight } from './processes/rockfall.js';
export { VIEWS, cellColors, hillshade, pixelMap, paint, crossSection } from './render.js';
export { Rng, fbm } from './rng.js';
