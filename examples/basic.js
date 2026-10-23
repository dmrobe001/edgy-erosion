// Programmatic use, no browser needed:   node examples/basic.js
//
// Builds a custom landscape (a layered plateau with a volcanic plug), runs it while
// logging a few numbers, then changes the climate halfway through.

import { writeFileSync } from 'node:fs';
import { Simulation } from '../src/index.js';

const sim = new Simulation({
  preset: 'mesa-country', // start from a preset, then override what you like
  seed: 42,
  grid: { cols: 96, rows: 80, dx: 25 },
  bodies: [{ shape: 'disc', x: 0.5, y: 0.5, radius: 0.05, rock: 'basalt' }],
  climate: { precipitation: 0.25 },
});

console.log('time_kyr,max_elevation_m,mean_soil_m,rockfalls');
for (let s = 1; s <= 200; s++) {
  if (s === 100) sim.setParams({ climate: { precipitation: 0.8, frost: 0.8 } }); // wetter, colder
  sim.step();
  if (s % 20 === 0) {
    const st = sim.stats();
    console.log([st.time / 1000, st.maxElevation.toFixed(1), st.meanSoil.toFixed(3), st.rockfallEvents].join(','));
  }
}

// The raw arrays are ordinary typed arrays: sim.z, sim.rock, sim.soil, sim.flow.area, ...
// and sim.grid.x / sim.grid.y give the cell centres in metres.
writeFileSync('example-snapshot.json', JSON.stringify(sim.snapshot()));
console.log('wrote example-snapshot.json (reload with Simulation.fromSnapshot, or bin/erode.js --params)');
