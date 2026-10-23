import { test } from 'node:test';
import assert from 'node:assert/strict';
import { HexGrid, OPPOSITE } from '../src/hexgrid.js';

test('every cell centre is dx from each of its neighbours, and neighbourship is symmetric', () => {
  const g = new HexGrid(9, 7, 2.5);
  for (let i = 0; i < g.size; i++) {
    for (let k = 0; k < 6; k++) {
      const j = g.neighbors[6 * i + k];
      if (j < 0) continue;
      const d = Math.hypot(g.x[j] - g.x[i], g.y[j] - g.y[i]);
      assert.ok(Math.abs(d - g.dx) < 1e-9, `distance ${d} between ${i} and ${j}`);
      assert.equal(g.neighbors[6 * j + OPPOSITE[k]], i);
    }
  }
});

test('interior cells have six neighbours; corners fewer', () => {
  const g = new HexGrid(6, 6, 1);
  const count = (i) => [0, 1, 2, 3, 4, 5].filter((k) => g.neighbors[6 * i + k] >= 0).length;
  assert.equal(count(g.index(2, 2)), 6);
  assert.equal(count(g.index(3, 3)), 6);
  assert.ok(count(g.index(0, 0)) < 6);
});

test('cell area and edge length match hexagon geometry', () => {
  const g = new HexGrid(4, 4, 10);
  // hexagon of side a has area (3 sqrt3 / 2) a^2
  const a = g.edgeLength;
  assert.ok(Math.abs(g.cellArea - 1.5 * Math.sqrt(3) * a * a) < 1e-9);
  assert.ok(Math.abs(a - 10 / Math.sqrt(3)) < 1e-12);
});

test('cellAt finds the cell under its own centre and corners stay inside', () => {
  const g = new HexGrid(11, 8, 3);
  for (let i = 0; i < g.size; i++) {
    assert.equal(g.cellAt(g.x[i], g.y[i]), i);
    // a point 90% of the way to each corner is still in cell i
    for (const [cx, cy] of g.corners(i)) {
      assert.equal(g.cellAt(g.x[i] + 0.9 * (cx - g.x[i]), g.y[i] + 0.9 * (cy - g.y[i])), i);
    }
  }
  assert.equal(g.cellAt(-100, -100), -1);
});
