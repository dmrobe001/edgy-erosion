import { test } from 'node:test';
import assert from 'node:assert/strict';
import { migrate, curvature } from '../src/meander.js';
import { Simulation } from '../src/simulation.js';
import { Detail, COVER } from '../src/detail.js';

test('meander migration grows bends from a straight channel, cuts off necks and keeps its ends fixed', () => {
  const nodes = [];
  for (let k = 0; k <= 120; k++) nodes.push({ x: k * 40, y: 0, w: 40, d: 2, zw: 0 });
  const belt = 600;
  const r = migrate(nodes, { iterations: 800, maxSinuosity: 9, seed: 2, erodibility: (x, y) => (Math.abs(y) < belt ? 1 : 0) });
  assert.ok(r.sinuosity > 1.5, `sinuosity ${r.sinuosity}`);
  assert.ok(r.oxbows.length > 0, 'loops were cut off into oxbows');
  const first = r.nodes[0], last = r.nodes[r.nodes.length - 1];
  assert.deepEqual([first.x, first.y], [0, 0]);
  assert.deepEqual([last.x, last.y], [4800, 0]);
  // no step moved the line out of its belt by more than a width or so
  assert.ok(r.nodes.every((nd) => Math.abs(nd.y) < belt + 2 * nd.w));
  // and no two stretches of the live channel overlap (necks were cut before they could cross)
  const n = r.nodes;
  for (let i = 0; i < n.length; i += 3) {
    for (let j = i + 20; j < n.length; j += 3) assert.ok(Math.hypot(n[i].x - n[j].x, n[i].y - n[j].y) > 0.5 * n[i].w);
  }
  assert.ok(Math.max(...curvature(n).map(Math.abs)) > 0);
});

// one coarse run shared by the tests below
const sim = new Simulation({ preset: 'river-coast', grid: { cols: 96, rows: 96, dx: 25 } });
sim.run(120);
const detail = new Detail(sim);

test('rivers become curves with widths from discharge, wider downstream, reaching the sea', () => {
  assert.ok(detail.rivers.length >= 3);
  const big = detail.rivers[detail.rivers.length - 1]; // sorted by width: the largest last
  assert.ok(big.wmax > 30, `widest river ${big.wmax} m`);
  // hydraulic geometry: W = aW Q^0.5
  const h = detail.hydraulics(3e9);
  assert.ok(Math.abs(h.w - detail.o.widthCoefficient * Math.sqrt(h.Q)) < 1e-9);
  for (const r of detail.rivers) {
    for (let k = 1; k < r.nodes.length; k++) {
      assert.ok(r.nodes[k].w >= r.nodes[k - 1].w - 1e-9, 'width never shrinks downstream');
      assert.ok(r.nodes[k].zw <= r.nodes[k - 1].zw + 1e-9, 'water surface never rises downstream');
    }
  }
  assert.ok(detail.rivers.some((r) => r.meandering), 'a big lowland river meanders');
});

test('a rendered map has wide rivers, sand and water, and tiles rendered apart join exactly', () => {
  const b = detail.bounds;
  const x0 = b.xmin + 300, y0 = b.ymin + 300;
  const whole = detail.render({ x0, y0, width: 240, height: 160, pixel: 4 });
  const counts = new Array(14).fill(0);
  for (const c of whole.cover) counts[c]++;
  assert.ok(counts[COVER.river] > 50 || counts[COVER.lake] > 50, 'river water is several pixels wide');
  // the same window as two halves
  const left = detail.render({ x0, y0, width: 120, height: 160, pixel: 4 });
  const right = detail.render({ x0: x0 + 120 * 4, y0, width: 120, height: 160, pixel: 4 });
  let worst = 0, coverDiff = 0;
  for (let j = 0; j < 160; j++) {
    for (let i = 0; i < 240; i++) {
      const t = i < 120 ? left : right;
      const k = j * 120 + (i % 120);
      worst = Math.max(worst, Math.abs(t.z[k] - whole.z[j * 240 + i]));
      if (t.cover[k] !== whole.cover[j * 240 + i]) coverDiff++;
    }
  }
  assert.ok(worst < 1e-3, `elevations differ by up to ${worst} m across the seam`);
  assert.ok(coverDiff === 0, `${coverDiff} pixels change cover`);
});
