// Public API. Works unchanged in browsers (as an ES module) and in Node >= 18.
export { Simulation, referenceFrame } from './simulation.js';
export { HexGrid } from './hexgrid.js';
export { Lithology } from './lithology.js';
export { DEFAULTS, PRESETS, ROCKS, resolveParams, deepMerge } from './params.js';
export { routeFlow } from './flow.js';
export { culmannHeight } from './processes/rockfall.js';
export { ruptureLength, averageDisplacement, ruptureWidth, sampleMagnitude, faultRate } from './processes/earthquakes.js';
export { VIEWS, cellColors, hillshade, pixelMap, paint, crossSection, PLANET_VIEWS, planetColors, spherePixelMap, planetSection, plateColor } from './render.js';
export { Rng, fbm, fbm3 } from './rng.js';
export { checkMesh, GraphMesh } from './mesh.js';
export { SphereMesh, azimuthalProjection, lonLatToVector, vectorToLonLat } from './sphere.js';
export { Planet, PLANET_DEFAULTS, PLANET_PRESETS, resolvePlanetParams } from './planet.js';
export { zoomIn, zoomFromPlanet, compareWithParent, interpolator, regionalErosionParams } from './nest.js';
