// The mesh contract: what the processes are allowed to know about the grid.
//
// Every process is written in finite-volume form against this interface, so a
// mesh does not have to be hexagonal, regular or flat. A mesh provides:
//
//   size            number of cells
//   degree          stride of the neighbour arrays (the maximum number of neighbours)
//   neighbors       Int32Array(degree * size); slot k of cell i is neighbors[degree*i + k],
//                   or -1 if unused (slots may be empty anywhere, not only at the end)
//   dist            Float64Array(degree * size): distance between the two cell centres (m)
//   face            Float64Array(degree * size): length of the face the two cells share (m)
//   area            Float64Array(size): cell area (m^2)
//   x, y            Float64Array(size): cell centres in world coordinates (m)
//   edge            Uint8Array(size): which domain sides a cell lies on
//                   (1 west, 2 east, 4 north, 8 south; 0 for interior cells)
//   bounds          { xmin, ymin, xmax, ymax } of the domain (m)
//
// Optional capabilities:
//   cellAt(x, y)    cell containing a point, or -1 (needed for drawing and zooming)
//   lattice         true if slot k always points along azimuth 60k degrees, with
//                   equal spacing (the hexagonal lattice); the wave-exposure sweep needs it
//
// The finite-volume rule every diffusive process uses: a flux q (m^2/yr, per unit face
// length) across the face between i and j moves volume q * face * dt, which changes
// the heights by -q*face*dt/area[i] and +q*face*dt/area[j]. For a gradient flux
// q = D (z_i - z_j)/dist this is the standard two-point flux approximation; on the
// regular hex grid it reduces to the isotropic hex Laplacian.

import { Rng } from './rng.js';
import { OPPOSITE } from './hexgrid.js';

/** Throw if a mesh breaks the contract (sizes, symmetry, positivity). */
export function checkMesh(m) {
  const fail = (msg) => { throw new Error(`invalid mesh: ${msg}`); };
  const { size: n, degree: D } = m;
  if (!(n > 0 && D > 0)) fail('size and degree must be positive');
  for (const [k, len] of [['neighbors', D * n], ['dist', D * n], ['face', D * n], ['area', n], ['x', n], ['y', n], ['edge', n]]) {
    if (!m[k] || m[k].length !== len) fail(`${k} must have length ${len}`);
  }
  for (let i = 0; i < n; i++) {
    if (!(m.area[i] > 0)) fail(`cell ${i} has non-positive area`);
    for (let k = 0; k < D; k++) {
      const j = m.neighbors[D * i + k];
      if (j < 0) continue;
      if (j >= n || j === i) fail(`cell ${i} slot ${k} points at ${j}`);
      if (!(m.dist[D * i + k] > 0 && m.face[D * i + k] > 0)) fail(`edge ${i}-${j} has non-positive dist or face`);
      let back = -1;
      for (let q = 0; q < D; q++) if (m.neighbors[D * j + q] === i) back = D * j + q;
      if (back < 0) fail(`edge ${i}->${j} has no reverse edge`);
      if (Math.abs(m.dist[back] - m.dist[D * i + k]) > 1e-9 * m.dist[back] ||
          Math.abs(m.face[back] - m.face[D * i + k]) > 1e-9 * m.face[back]) fail(`edge ${i}-${j} is not symmetric`);
    }
  }
  return true;
}

/**
 * A mesh given directly as arrays. `GraphMesh.irregular(hex, ...)` derives a
 * deliberately messy one from a hex grid for testing: jittered centres, missing
 * edges (so neighbour counts vary from 2 to 6, with gaps between slots), and
 * randomly perturbed face widths and areas. It satisfies the contract but is NOT
 * a geometrically consistent Voronoi tessellation, so it checks that processes
 * conserve mass and route water on any mesh, not that they converge to the PDE.
 */
export class GraphMesh {
  constructor({ size, degree, neighbors, dist, face, area, x, y, edge, bounds, cellAt }) {
    Object.assign(this, { size, degree, neighbors, dist, face, area, x, y, edge, bounds });
    if (cellAt) this.cellAt = cellAt;
    this.lattice = false;
    checkMesh(this);
  }

  static irregular(hex, { jitter = 0.25, drop = 0.15, perturb = 0.3, seed = 1 } = {}) {
    const rng = new Rng(seed);
    const n = hex.size;
    const D = hex.degree;
    const x = Float64Array.from(hex.x, (v) => v + jitter * hex.dx * (rng.next() - 0.5));
    const y = Float64Array.from(hex.y, (v) => v + jitter * hex.dx * (rng.next() - 0.5));
    const neighbors = Int32Array.from(hex.neighbors);
    const degreeOf = (i) => { let c = 0; for (let k = 0; k < D; k++) if (neighbors[D * i + k] >= 0) c++; return c; };
    const slot = (i, j) => { for (let k = 0; k < D; k++) if (neighbors[D * i + k] === j) return D * i + k; return -1; };
    // drop interior edges, keeping every cell at two or more neighbours
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < D; k++) {
        const j = neighbors[D * i + k];
        if (j <= i || rng.next() >= drop) continue;
        if (degreeOf(i) <= 2 || degreeOf(j) <= 2) continue;
        neighbors[D * i + k] = -1;
        neighbors[slot(j, i)] = -1;
      }
    }
    // restore edges until every cell is connected to cell 0
    for (;;) {
      const seen = new Uint8Array(n);
      const stack = [0];
      seen[0] = 1;
      while (stack.length) {
        const i = stack.pop();
        for (let k = 0; k < D; k++) { const j = neighbors[D * i + k]; if (j >= 0 && !seen[j]) { seen[j] = 1; stack.push(j); } }
      }
      const lost = seen.indexOf(0);
      if (lost < 0) break;
      // give the unreached cell back its original hex edges (slot k pairs with OPPOSITE[k])
      for (let k = 0; k < D; k++) {
        const j = hex.neighbors[D * lost + k];
        if (j >= 0) { neighbors[D * lost + k] = j; neighbors[D * j + OPPOSITE[k]] = lost; }
      }
    }
    const dist = new Float64Array(D * n);
    const face = new Float64Array(D * n);
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < D; k++) {
        const j = neighbors[D * i + k];
        if (j < i) continue;
        const d = Math.hypot(x[j] - x[i], y[j] - y[i]);
        const f = hex.edgeLength * (1 + perturb * (rng.next() - 0.5));
        const back = slot(j, i);
        dist[D * i + k] = dist[back] = d;
        face[D * i + k] = face[back] = f;
      }
    }
    const area = Float64Array.from(hex.area, (a) => a * (1 + perturb * (rng.next() - 0.5)));
    return new GraphMesh({
      size: n, degree: D, neighbors, dist, face, area, x, y,
      edge: Uint8Array.from(hex.edge), bounds: { ...hex.bounds }, cellAt: (px, py) => hex.cellAt(px, py),
    });
  }
}
