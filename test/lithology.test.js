import { test } from 'node:test';
import assert from 'node:assert/strict';
import { HexGrid } from '../src/hexgrid.js';
import { Lithology } from '../src/lithology.js';
import { resolveParams } from '../src/params.js';

const grid = new HexGrid(21, 21, 10);
const centre = grid.index(10, 10);
const lith = (strata, bodies = []) =>
  new Lithology(resolveParams({ strata: { strengthVariation: { amplitude: 0 }, ...strata }, bodies }), grid);

test('horizontal layers are found at the right depths', () => {
  const L = lith({ top: 100, layers: [{ rock: 'sandstone', thickness: 10 }, { rock: 'shale', thickness: 20 }], basement: 'granite' });
  const name = (z) => L.names[L.rockAt(centre, z)];
  assert.equal(name(105), 'sandstone'); // above the column top: top layer extends up
  assert.equal(name(95), 'sandstone');
  assert.equal(name(85), 'shale');
  assert.equal(name(71), 'shale');
  assert.equal(name(69), 'granite');
});

test('repeating layers cycle downwards', () => {
  const L = lith({ top: 0, repeat: true, layers: [{ rock: 'limestone', thickness: 5 }, { rock: 'shale', thickness: 5 }] });
  const name = (z) => L.names[L.rockAt(centre, z)];
  assert.equal(name(-2), 'limestone');
  assert.equal(name(-7), 'shale');
  assert.equal(name(-12), 'limestone');
  assert.equal(name(-1007), 'shale');
});

test('dipping layers descend towards the dip direction', () => {
  const L = lith({ top: 0, dip: 45, dipDirection: 0, layers: [{ rock: 'limestone', thickness: 1e9 }] });
  const east = grid.index(15, 10), west = grid.index(5, 10);
  assert.ok(L.columnTop[east] < L.columnTop[centre] && L.columnTop[centre] < L.columnTop[west]);
  // 45 degrees: one metre down per metre east
  assert.ok(Math.abs((L.columnTop[west] - L.columnTop[east]) - (grid.x[east] - grid.x[west])) < 1e-9);
});

test('a plug overrides strata below its top only', () => {
  const L = lith({ top: 0, layers: [{ rock: 'shale', thickness: 1e9 }] },
    [{ shape: 'disc', x: 0.5, y: 0.5, radius: 0.1, rock: 'basalt', top: -10 }]);
  const name = (i, z) => L.names[L.rockAt(i, z)];
  assert.equal(name(centre, -20), 'basalt');
  assert.equal(name(centre, -5), 'shale');
  assert.equal(name(grid.index(1, 1), -20), 'shale');
});

test('face strength averages the layers a cliff cuts through', () => {
  const L = lith({ top: 0, layers: [{ rock: 'sandstone', thickness: 50 }, { rock: 'shale', thickness: 1e9 }] });
  const s = L.meanStrength(centre, -100, 0, 100);
  const expected = 0.5 * (L.cohesion[L.names.indexOf('sandstone')] + L.cohesion[L.names.indexOf('shale')]);
  assert.ok(Math.abs(s.cohesion - expected) / expected < 1e-9);
});
