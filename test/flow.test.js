import { test } from 'node:test';
import assert from 'node:assert/strict';
import { HexGrid } from '../src/hexgrid.js';
import { routeFlow } from '../src/flow.js';
import { Rng } from '../src/rng.js';

function setup(cols, rows, zfun, opts = {}) {
  const g = new HexGrid(cols, rows, 10);
  const n = g.size;
  const z = new Float64Array(n);
  const outlet = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    z[i] = zfun(i, g);
    if (g.edge[i]) outlet[i] = 1;
  }
  const out = {
    filled: new Float64Array(n), receiver: new Int32Array(n), order: new Int32Array(n),
    area: new Float64Array(n), visited: new Uint8Array(n), ...opts,
  };
  routeFlow(g, z, outlet, new Float64Array(n).fill(1), out);
  return { g, z, outlet, out };
}

const rough = (seed) => { const r = new Rng(seed); return () => 100 * r.next(); };

for (const [name, opts] of [['single flow', {}], ['multiple flow', { flowExponent: 1.1 }], ['random receivers', { randomExponent: 2, rng: new Rng(5) }]]) {
  test(`${name}: no pits remain, receivers are lower, order is topological, area is conserved`, () => {
    const f = rough(3);
    const { g, z, outlet, out } = setup(25, 20, f, opts);
    const pos = new Int32Array(g.size);
    out.order.forEach((i, s) => { pos[i] = s; });
    let outletArea = 0;
    for (let i = 0; i < g.size; i++) {
      assert.ok(out.filled[i] >= z[i]);
      if (outlet[i]) { outletArea += out.area[i]; assert.equal(out.receiver[i], -1); continue; }
      const r = out.receiver[i];
      assert.ok(r >= 0, 'every interior cell drains somewhere');
      assert.ok(out.filled[r] < out.filled[i], 'receiver is lower on the filled surface');
      assert.ok(pos[r] < pos[i], 'receiver comes earlier in the stack');
    }
    // all rain ends up at an outlet
    assert.ok(Math.abs(outletArea - g.size * g.cellArea) < 1e-6 * g.size * g.cellArea);
  });
}

test('a closed basin fills to its spill point', () => {
  // a bowl: high rim at distance, deep centre; outlets at the edge with z = 50
  const { g, z, out } = setup(15, 15, (i, g) => {
    if (g.edge[i]) return 50;
    const dx = g.x[i] - g.x[g.index(7, 7)], dy = g.y[i] - g.y[g.index(7, 7)];
    return Math.hypot(dx, dy) < 30 ? 0 : 80;
  });
  const centre = g.index(7, 7);
  // the lake surface sits at the lowest rim level (80), not at the bottom (0)
  assert.ok(out.filled[centre] >= 80 && out.filled[centre] < 80.01);
  assert.equal(z[centre], 0);
});
