// Zoom into a patch of a running landscape at four times finer spacing, and compare the
// refined patch with its frozen parent:   node examples/zoom.js
import { Simulation, zoomIn, compareWithParent } from '../src/index.js';

const parent = new Simulation({ preset: 'alpine' });
parent.run(300); // 600 kyr at dt = 2000 yr

// centre the patch on the largest river away from the domain edges
const g = parent.grid;
let trunk = -1;
for (let i = 0; i < g.size; i++) {
  const c = g.col(i), r = g.row(i);
  if (c < 40 || c > 88 || r < 40 || r > 88) continue;
  if (trunk < 0 || parent.flow.area[i] > parent.flow.area[trunk]) trunk = i;
}

const child = zoomIn(parent, { center: [g.x[trunk], g.y[trunk]], factor: 4, cols: 96, rows: 96 });
console.log(`parent dx ${g.dx} m -> child dx ${child.grid.dx} m, dt ${child.params.time.dt} yr`);
console.log(`drainage area entering the patch from outside: ${(child.zoom.enteringArea / 1e6).toFixed(3)} km²`);

const report = () => {
  const c = compareWithParent(child, parent);
  console.log(
    `after ${(c.elapsed / 1000).toFixed(0).padStart(3)} kyr: mean z ${c.meanElevation.child.toFixed(1)} m ` +
    `(parent at its last rate ${c.meanElevation.expected.toFixed(1)} m), ` +
    `largest river out ${(c.largestExitArea.child / 1e6).toFixed(3)} km² (parent ${(c.largestExitArea.parent / 1e6).toFixed(3)} km²)`,
  );
};
report();
for (let k = 0; k < 3; k++) { child.run(100); report(); }
