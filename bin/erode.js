#!/usr/bin/env node
// Command-line runner: simulate without a browser and write images and data.
//
//   node bin/erode.js --preset mesa-country --steps 400 --out out/mesa
//   node bin/erode.js --params my-run.json --steps 1000 --views terrain,geology --section 64
//   node bin/erode.js --list-presets
//
// Run with --help for every option.

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { Simulation, PRESETS, VIEWS, cellColors, pixelMap, paint, crossSection } from '../src/index.js';
import { encodePng } from './png.js';

const HELP = `Usage: node bin/erode.js [options]

  --preset NAME       start from a preset (see --list-presets)
  --params FILE       JSON file of parameter overrides (applied after --preset);
                      may also be a snapshot written by --snapshot, to continue a run
  --steps N           number of time steps (default 200)
  --dt YEARS          time step length
  --cols N --rows N   grid size in cells
  --dx METRES         lattice spacing
  --seed N            random seed
  --out DIR           output directory (default out)
  --views LIST        comma-separated views to render (default terrain,geology)
                      available: ${Object.keys(VIEWS).join(', ')}
  --width PX          image width in pixels (default 800)
  --section ROW       also render a cross-section along this grid row ("mid" = middle)
  --every N           print stats every N steps (default: 10 lines total)
  --frames N          also save a terrain image every N steps
  --snapshot          write snapshot.json (full state; reload with --params)
  --heightmap         write heightmap.pgm (16-bit, resampled to square pixels) + heightmap.json
  --list-presets      list presets and exit
`;

function parseArgs(argv) {
  const a = {};
  for (let k = 0; k < argv.length; k++) {
    const key = argv[k];
    if (!key.startsWith('--')) throw new Error(`unexpected argument "${key}"`);
    const name = key.slice(2);
    const flag = ['snapshot', 'heightmap', 'list-presets', 'help'].includes(name);
    a[name] = flag ? true : argv[++k];
    if (!flag && a[name] === undefined) throw new Error(`${key} needs a value`);
  }
  return a;
}

const args = parseArgs(process.argv.slice(2));
if (args.help) { process.stdout.write(HELP); process.exit(0); }
if (args['list-presets']) {
  for (const [k, v] of Object.entries(PRESETS)) console.log(`${k.padEnd(18)} ${v.description}`);
  process.exit(0);
}

let sim;
let overrides = {};
if (args.params) {
  const loaded = JSON.parse(readFileSync(args.params, 'utf8'));
  if (loaded.format && loaded.format.startsWith('edgy-erosion-snapshot')) {
    sim = Simulation.fromSnapshot(loaded);
  } else {
    overrides = loaded;
  }
}
if (!sim) {
  if (args.preset) overrides.preset = args.preset;
  const num = (k) => (args[k] !== undefined ? Number(args[k]) : undefined);
  const grid = { cols: num('cols'), rows: num('rows'), dx: num('dx') };
  overrides.grid = { ...overrides.grid, ...Object.fromEntries(Object.entries(grid).filter(([, v]) => v !== undefined)) };
  if (args.dt) overrides.time = { ...overrides.time, dt: Number(args.dt) };
  if (args.seed) overrides.seed = Number(args.seed);
  sim = new Simulation(overrides);
} else if (args.dt) {
  sim.setParams({ time: { dt: Number(args.dt) } });
}

const outDir = args.out || 'out';
mkdirSync(outDir, { recursive: true });
const steps = Number(args.steps ?? 200);
const every = Number(args.every ?? Math.max(1, Math.round(steps / 10)));
const frames = args.frames ? Number(args.frames) : 0;
const views = (args.views || 'terrain,geology').split(',').filter(Boolean);
const g = sim.grid;
const width = Number(args.width ?? 800);
const height = Math.round((width * g.height) / g.width);
const pm = pixelMap(g, width, height);
const rgba = new Uint8ClampedArray(width * height * 4);

function writeView(view, file) {
  paint(pm, cellColors(sim, view), rgba);
  writeFileSync(file, encodePng(width, height, rgba));
}

const p = sim.params;
console.log(`grid ${g.cols} x ${g.rows} hexes, dx = ${g.dx} m, extent ${(g.width / 1000).toFixed(2)} x ${(g.height / 1000).toFixed(2)} km, dt = ${p.time.dt} yr`);
const fmt = (s) =>
  `t=${(s.time / 1000).toFixed(1).padStart(8)} kyr  z=[${s.minElevation.toFixed(1)}, ${s.maxElevation.toFixed(1)}] m  ` +
  `mean ${s.meanElevation.toFixed(1)} m  soil ${s.meanSoil.toFixed(2)} m  bare ${(100 * s.bareRockFraction).toFixed(0)}%  ` +
  `rockfalls ${s.rockfallEvents}  creep substeps ${s.creepSubsteps}`;

const t0 = Date.now();
for (let s = 1; s <= steps; s++) {
  sim.step();
  if (s % every === 0 || s === steps) console.log(fmt(sim.stats()));
  if (frames && s % frames === 0) writeView('terrain', join(outDir, `frame-${String(s).padStart(5, '0')}.png`));
}
const secs = (Date.now() - t0) / 1000;
console.log(`${steps} steps in ${secs.toFixed(2)} s (${((1000 * secs) / Math.max(1, steps)).toFixed(1)} ms/step)`);

for (const v of views) {
  if (!VIEWS[v]) throw new Error(`unknown view "${v}"`);
  writeView(v, join(outDir, `${v}.png`));
}
if (args.section !== undefined) {
  const row = args.section === 'mid' ? g.rows >> 1 : Number(args.section);
  const sw = width, sh = Math.round(width / 3);
  const buf = new Uint8ClampedArray(sw * sh * 4);
  const info = crossSection(sim, row, sw, sh, buf);
  writeFileSync(join(outDir, `section-row${row}.png`), encodePng(sw, sh, buf));
  console.log(`section row ${row}: z ${info.zmin.toFixed(0)}..${info.zmax.toFixed(0)} m, vertical exaggeration ${info.exaggeration.toFixed(1)}x`);
}
if (args.snapshot) writeFileSync(join(outDir, 'snapshot.json'), JSON.stringify(sim.snapshot()));
if (args.heightmap) {
  let lo = Infinity, hi = -Infinity;
  for (const v of sim.z) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  const header = Buffer.from(`P5\n${width} ${height}\n65535\n`, 'ascii');
  const body = Buffer.alloc(width * height * 2);
  pm.map.forEach((i, k) => body.writeUInt16BE(i < 0 ? 0 : Math.round((65535 * (sim.z[i] - lo)) / Math.max(hi - lo, 1e-9)), 2 * k));
  writeFileSync(join(outDir, 'heightmap.pgm'), Buffer.concat([header, body]));
  writeFileSync(join(outDir, 'heightmap.json'), JSON.stringify({ minElevation: lo, maxElevation: hi, metresPerPixel: 1 / pm.scale, width, height }, null, 2));
}
writeFileSync(join(outDir, 'params.json'), JSON.stringify(sim.params, null, 2));
console.log(`wrote ${outDir}/`);
