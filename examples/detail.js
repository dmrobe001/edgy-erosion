// From a landscape simulation to a walking-scale map:  node examples/detail.js [out-dir]
// Runs the river-coast preset (25 m cells), then details it: meandering rivers with
// widths, point bars, oxbows, floodplains, beaches and dunes, tree crowns. Writes the
// whole patch at 4 m per pixel and some close-ups at 1 m and 0.5 m.
import { mkdirSync, writeFileSync } from 'node:fs';
import { Simulation, Detail, detailColors } from '../src/index.js';
import { encodePng } from '../bin/png.js';

const out = process.argv[2] || 'out/detail';
mkdirSync(out, { recursive: true });

const sim = new Simulation('river-coast');
sim.run(200); // 20 kyr at dt = 100 yr
const detail = new Detail(sim);
const big = detail.rivers[detail.rivers.length - 1];
console.log(`${detail.rivers.length} rivers, the largest ${big.wmax.toFixed(0)} m wide` +
  `${big.meandering ? ` with sinuosity ${big.sinuosity.toFixed(2)}` : ''}; ${detail.oxbowCount} cut-off bends`);

const save = (name, x0, y0, width, height, pixel) => {
  const t0 = Date.now();
  const tile = detail.render({ x0, y0, width, height, pixel });
  writeFileSync(`${out}/${name}.png`, encodePng(width, height, detailColors(detail, tile)));
  console.log(`${name}: ${width} x ${height} px at ${pixel} m (${Date.now() - t0} ms)`);
};
const b = detail.bounds;
save('patch-4m', b.xmin, b.ymin, 1000, Math.round((1000 * (b.ymax - b.ymin)) / (b.xmax - b.xmin)), (b.xmax - b.xmin) / 1000);
// close-ups: the tightest bend of the largest river, and where it meets the sea
let bend = big.nodes[0];
for (const nd of big.nodes) if (Math.abs(nd.c) > Math.abs(bend.c) && nd.x > b.xmin + 500 && nd.x < b.xmax - 500) bend = nd;
save('river-1m', bend.x - 450, bend.y - 350, 900, 700, 1);
save('river-50cm', bend.x - 225, bend.y - 175, 900, 700, 0.5);
const mouth = big.nodes.find((nd) => detail.z.sample(nd.x, nd.y) < detail.sea) || big.nodes[big.nodes.length - 1];
save('coast-1m', mouth.x - 450, mouth.y - 500, 900, 700, 1);
console.log(`wrote ${out}/`);
