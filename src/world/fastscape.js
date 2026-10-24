// Implicit stream power and hillslope diffusion on a square raster (level 1).
//
//   dz/dt = U - K A^m S^n + D lap(z)
//
// Fluvial term: the O(n) implicit scheme of Braun & Willett (2013, Geomorphology 180):
// each cell drains to one receiver; visiting cells from base level upstream, each new
// elevation depends only on its receiver's new elevation, already known:
//   n = 1:  z_i' = (z_i + U dt + F z_r') / (1 + F),   F = K dt A^m / L
// (for n != 1 the same equation, with (z_i' - z_r')^n, is solved by Newton's method).
// It is unconditionally stable, so a step can be 10^4-10^5 years at kilometre cells.
//
// Routing: a priority flood (Barnes, Lehman & Mulla 2014, Computers & Geosciences 62)
// from the outlet cells gives every cell a path to an outlet, depressions included, in
// a single pass: cells come off the queue in an order in which every cell's receiver
// came off before it, which is exactly the visiting order the implicit scheme needs.
// A cell's receiver is its steepest-descent neighbour among cells already visited, or,
// in a depression or flat, the cell the flood reached it from. A cell in a depression
// (receiver higher than itself) does not erode; it is a lake.
//
// Hillslope term: implicit diffusion split by direction (an implicit solve along every
// row, then along every column, each a tridiagonal system): first order in time and
// unconditionally stable.
//
// Fixed cells (Dirichlet) keep their elevation: tile edges take their profile from the
// parent and are not drains. Outlet cells are the only sinks: designated river exits
// and the sea.

import { MinHeap } from '../flow.js';

const D8 = [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, 1], [1, -1], [-1, -1]];

export class RasterSolver {
  /**
   * @param {number} nx columns
   * @param {number} ny rows
   * @param {number} dx cell size (m)
   */
  constructor(nx, ny, dx) {
    this.nx = nx; this.ny = ny; this.dx = dx;
    const n = nx * ny;
    this.size = n;
    this.receiver = new Int32Array(n);
    this.recDist = new Float64Array(n);
    this.order = new Int32Array(n);
    this.orderCount = 0;
    this.area = new Float64Array(n);
    this.level = new Float64Array(n);
    this.visited = new Uint8Array(n);
    this.heap = new MinHeap(n + 8);
    // neighbour offsets and distances
    this.nbOff = D8.map(([a, b]) => a + b * nx);
    this.nbDist = D8.map(([a, b]) => dx * Math.hypot(a, b));
    const m = Math.max(nx, ny);
    this._a = new Float64Array(m); this._b = new Float64Array(m); this._c = new Float64Array(m); this._r = new Float64Array(m);
  }

  /**
   * Receivers and visiting order. Fast path: steepest-descent receivers and a
   * breadth-first walk up the drainage tree from the outlets, O(n). If some cells are not
   * reached (they drain into a depression, or a cycle through uphill-routed fixed cells),
   * the priority flood is used instead. `this.flooded` says which ran.
   *
   * Fixed cells that are not outlets are walls, not drains: no water flows into them.
   * Their own (small) runoff goes to their lowest open neighbour, uphill if need be.
   * @param {Float64Array} z
   * @param {Uint8Array} outlet 1 where water leaves
   * @param {Uint8Array} [fixed]
   */
  route(z, outlet, fixed) {
    if (this.routeDescent(z, outlet, fixed)) { this.flooded = false; return this; }
    this.flooded = true;
    return this.routeFlood(z, outlet, fixed);
  }

  /** Steepest-descent receivers; returns false if some cell has no path to an outlet. */
  routeDescent(z, outlet, fixed) {
    const { nx, ny, receiver, recDist, order, nbOff, nbDist } = this;
    const n = this.size;
    const count = this._count || (this._count = new Int32Array(n + 1));
    const donors = this._donors || (this._donors = new Int32Array(n));
    count.fill(0);
    for (let cj = 0; cj < ny; cj++) {
      const edgeRow = cj === 0 || cj === ny - 1;
      for (let ci = 0; ci < nx; ci++) {
        const i = cj * nx + ci;
        receiver[i] = i; recDist[i] = 0;
        if (outlet[i]) continue;
        const wall = fixed && fixed[i];
        let best = wall ? -Infinity : 0;
        const interior = !edgeRow && ci > 0 && ci < nx - 1;
        for (let k = 0; k < 8; k++) {
          if (!interior) {
            const a = ci + D8[k][0], b = cj + D8[k][1];
            if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
          }
          const j = i + nbOff[k];
          if (fixed && fixed[j] && !outlet[j]) continue;
          const sl = (z[i] - z[j]) / nbDist[k];
          if (sl > best) { best = sl; receiver[i] = j; recDist[i] = nbDist[k]; }
        }
        if (receiver[i] === i) return false; // a pit
        count[receiver[i] + 1]++;
      }
    }
    // donors in compressed rows: donors of r are donors[count[r] .. count[r + 1])
    for (let i = 0; i < n; i++) count[i + 1] += count[i];
    const fill = this._fill || (this._fill = new Int32Array(n));
    fill.set(count.subarray(0, n));
    for (let i = 0; i < n; i++) { const r = receiver[i]; if (r !== i) donors[fill[r]++] = i; }
    let head = 0, tail = 0;
    for (let i = 0; i < n; i++) if (outlet[i]) order[tail++] = i;
    while (head < tail) {
      const r = order[head++];
      for (let k = count[r]; k < count[r + 1]; k++) order[tail++] = donors[k];
    }
    this.orderCount = tail;
    return tail === n;
  }

  /**
   * Receivers and visiting order by priority flood from the outlets.
   * @param {Float64Array} z
   * @param {Uint8Array} outlet 1 where water leaves
   */
  routeFlood(z, outlet, fixed) {
    const { nx, ny, receiver, recDist, order, level, visited, heap, nbOff, nbDist } = this;
    const n = this.size;
    visited.fill(0);
    heap.n = 0;
    let count = 0;
    for (let i = 0; i < n; i++) {
      if (!outlet[i]) continue;
      visited[i] = 1;
      level[i] = z[i];
      receiver[i] = i;
      recDist[i] = 0;
      heap.push(i, z[i]);
    }
    const done = this._done || (this._done = new Uint8Array(n));
    done.fill(0);
    while (heap.n) {
      const i = heap.pop();
      if (done[i]) continue;
      done[i] = 1;
      order[count++] = i;
      const ci = i % nx, cj = (i - ci) / nx;
      // steepest descent among cells already settled (if any lies below)
      if (!outlet[i]) {
        let best = 0, br = receiver[i], bd = recDist[i];
        for (let k = 0; k < 8; k++) {
          const a = ci + D8[k][0], b = cj + D8[k][1];
          if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
          const j = i + nbOff[k];
          if (!done[j] || j === i || (fixed && fixed[j] && !outlet[j])) continue;
          const s = (z[i] - z[j]) / nbDist[k];
          if (s > best) { best = s; br = j; bd = nbDist[k]; }
        }
        receiver[i] = br; recDist[i] = bd;
      }
      if (fixed && fixed[i] && !outlet[i]) continue; // walls pass no water on
      for (let k = 0; k < 8; k++) {
        const a = ci + D8[k][0], b = cj + D8[k][1];
        if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
        const j = i + nbOff[k];
        if (visited[j]) continue;
        visited[j] = 1;
        receiver[j] = i; // provisional: the flood parent
        recDist[j] = nbDist[k];
        level[j] = Math.max(z[j], level[i]);
        heap.push(j, level[j]);
      }
    }
    this.orderCount = count;
    return this;
  }

  /**
   * Drainage area (m^2, or discharge-weighted area with `weight`), plus any inflow.
   * @param {Float64Array} [weight] per-cell multiplier (relative runoff)
   * @param {Float64Array} [inflow] extra area entering at each cell from outside the tile
   */
  accumulate(weight, inflow) {
    const { area, receiver, order } = this;
    const cellArea = this.dx * this.dx;
    for (let i = 0; i < this.size; i++) area[i] = cellArea * (weight ? weight[i] : 1) + (inflow ? inflow[i] : 0);
    for (let k = this.orderCount - 1; k >= 0; k--) {
      const i = order[k];
      const r = receiver[i];
      if (r !== i) area[r] += area[i];
    }
    return area;
  }

  /**
   * One implicit stream-power step (after route and accumulate).
   * @param {Float64Array} z modified in place
   * @param {object} p {dt, K (number|array), m, n, U (number|array), fixed (Uint8Array), fillLakes}
   */
  erode(z, { dt, K, m = 0.5, n = 1, U = 0, fixed, fillLakes = true }) {
    const { receiver, recDist, order, area } = this;
    const Kv = typeof K === 'number' ? null : K, Uv = typeof U === 'number' ? null : U;
    for (let k = 0; k < this.orderCount; k++) {
      const i = order[k];
      if (fixed && fixed[i]) continue;
      const r = receiver[i];
      if (r === i) continue; // an outlet: base level stays put
      const z0 = z[i] + dt * (Uv ? Uv[i] : U);
      const zr = z[r];
      if (z0 <= zr) {
        // a lake: the receiver stands higher. With fillLakes, sediment fills it to just
        // above its outlet (lakes fill within a step at these time scales); this keeps the
        // surface free of depressions, so routing stays on the fast path.
        z[i] = fillLakes ? zr + 1e-6 * recDist[i] : z0;
        continue;
      }
      const F = ((Kv ? Kv[i] : K) * dt * area[i] ** m) / recDist[i] ** n;
      if (n === 1) { z[i] = (z0 + F * zr) / (1 + F); continue; }
      // Newton on g(x) = x - (z0 - zr) + F x^n, x = z' - zr > 0
      let x = z0 - zr;
      for (let it = 0; it < 20; it++) {
        const g = x - (z0 - zr) + F * x ** n;
        const dg = 1 + n * F * x ** (n - 1);
        const nx = Math.max(1e-12, x - g / dg);
        if (Math.abs(nx - x) < 1e-9) { x = nx; break; }
        x = nx;
      }
      z[i] = zr + x;
    }
  }

  /**
   * Implicit linear diffusion, rows then columns. Fixed cells are Dirichlet values.
   * @param {Float64Array} z
   * @param {{dt:number, D:number|Float64Array, fixed?:Uint8Array}} p
   */
  diffuse(z, { dt, D, fixed }) {
    const { nx, ny } = this;
    const Dv = typeof D === 'number' ? null : D;
    const sweep = (count, len, idx) => {
      const { _a: a, _b: b, _c: c, _r: r } = this;
      for (let line = 0; line < count; line++) {
        for (let t = 0; t < len; t++) {
          const i = idx(line, t);
          const d = (Dv ? Dv[i] : D) * dt / (this.dx * this.dx);
          if ((fixed && fixed[i]) || t === 0 || t === len - 1) { a[t] = 0; b[t] = 1; c[t] = 0; r[t] = z[i]; continue; }
          a[t] = -d; b[t] = 1 + 2 * d; c[t] = -d; r[t] = z[i];
        }
        // Thomas algorithm
        for (let t = 1; t < len; t++) {
          const w = a[t] / b[t - 1];
          b[t] -= w * c[t - 1];
          r[t] -= w * r[t - 1];
        }
        let prev = r[len - 1] / b[len - 1];
        z[idx(line, len - 1)] = prev;
        for (let t = len - 2; t >= 0; t--) {
          prev = (r[t] - c[t] * prev) / b[t];
          z[idx(line, t)] = prev;
        }
      }
    };
    sweep(ny, nx, (row, t) => row * nx + t);
    sweep(nx, ny, (col, t) => t * nx + col);
  }
}
