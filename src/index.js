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
export { Planet, PLANET_DEFAULTS, PLANET_PRESETS, resolvePlanetParams, MATERIAL, MATERIAL_NAMES } from './planet.js';
export { zoomIn, zoomFromPlanet, compareWithParent, interpolator, regionalErosionParams } from './nest.js';
export { Detail, DETAIL_DEFAULTS, COVER, detailColors } from './detail.js';
export { migrate } from './meander.js';
// Hierarchical generator (docs/DESIGN.md, docs/PLAN.md)
export { philox4x32, worldKey, deriveKey, keyHex, uniform, normal, mix32, fieldHash } from './world/hash.js';
export {
  WORLD_LAYOUT, faceToVector, vectorToFace, tilesPerSide, tileKey, parseTileKey, tileAt, tileRect, cellCenter,
  parentOf, childrenOf, neighborOf, cornerTiles, edgeKey, cornerKey,
} from './world/cubesphere.js';
export { generatePlanet, L0_DEFAULTS, GENERATOR_VERSION, seaLevelTemperature } from './world/l0.js';
export { RasterSolver } from './world/fastscape.js';
