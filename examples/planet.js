// Run a planet for 100 Myr, save maps, then zoom from the planet to a region and on
// to a landscape:   node examples/planet.js [preset] [out-dir]
import { mkdirSync, writeFileSync } from 'node:fs';
import {
  Planet, zoomFromPlanet, zoomIn, compareWithParent,
  planetColors, spherePixelMap, cellColors, pixelMap, paint,
} from '../src/index.js';
import { encodePng } from '../bin/png.js';

const [preset = 'earth-like', out = 'out/planet'] = process.argv.slice(2);
mkdirSync(out, { recursive: true });

const planet = new Planet({ preset, subdivisions: 40 });
const W = 1200, H = 600;
const pm = spherePixelMap(planet.mesh, W, H);
const saveMap = (name, view) => {
  const rgba = new Uint8ClampedArray(4 * W * H);
  paint(pm, planetColors(planet, view), rgba, [0, 0, 0]);
  writeFileSync(`${out}/${name}.png`, encodePng(W, H, rgba));
};
const saveSim = (name, sim) => {
  const w = 700, h = Math.round((w * sim.grid.height) / sim.grid.width);
  const rgba = new Uint8ClampedArray(4 * w * h);
  paint(pixelMap(sim.grid, w, h), cellColors(sim, 'terrain'), rgba);
  writeFileSync(`${out}/${name}.png`, encodePng(w, h, rgba));
};

for (let k = 1; k <= 4; k++) {
  planet.run(50); // 25 Myr at dt = 0.5 Myr
  const s = planet.stats();
  console.log(`${(s.time / 1e6).toFixed(0)} Myr: ${s.plates} plates, land ${(100 * s.landFraction).toFixed(0)}%, ` +
    `highest ${(s.maxElevation / 1000).toFixed(1)} km, sea level ${s.seaLevel.toFixed(0)} m`);
}
for (const v of ['terrain', 'plates', 'age', 'crust']) saveMap(`planet-${v}`, v);

// zoom onto the fastest-rising land
let c = -1;
for (let i = 0; i < planet.mesh.size; i++) {
  if (planet.sim.z[i] > planet.seaLevel + 300 && (c < 0 || planet.upliftRate[i] > planet.upliftRate[c])) c = i;
}
const lon = (planet.mesh.lon[c] * 180) / Math.PI, lat = (planet.mesh.lat[c] * 180) / Math.PI;
const region = zoomFromPlanet(planet, { lon, lat, factor: 32, cols: 128, rows: 128 });
console.log(`region at ${lat.toFixed(1)}°, ${lon.toFixed(1)}°: dx ${(region.grid.dx / 1000).toFixed(1)} km, dt ${region.params.time.dt.toFixed(0)} yr`);
region.run(80);
const cr = compareWithParent(region, planet);
console.log(`  after ${(cr.elapsed / 1e6).toFixed(2)} Myr: mean ${cr.meanElevation.child.toFixed(0)} m (planet extrapolated ${cr.meanElevation.expected.toFixed(0)} m)`);
saveSim('region', region);

const g = region.grid;
let trunk = -1;
for (let i = 0; i < g.size; i++) {
  const col = g.col(i), row = g.row(i);
  if (col < 40 || col > 88 || row < 40 || row > 88 || region.z[i] < region.seaLevel) continue;
  if (trunk < 0 || region.flow.area[i] > region.flow.area[trunk]) trunk = i;
}
if (trunk >= 0) {
  const land = zoomIn(region, { center: [g.x[trunk], g.y[trunk]], factor: 8, cols: 96, rows: 96 });
  console.log(`landscape: dx ${land.grid.dx.toFixed(0)} m`);
  land.run(60);
  saveSim('landscape', land);
}
console.log(`wrote ${out}/`);
