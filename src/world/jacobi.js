// Stream power for the GPU: every cell updates at once from the previous pass.
//
// The implicit solver in fastscape.js visits cells in drainage order, which is
// sequential by nature. Here an iteration is a fixed sequence of data-parallel passes
// over the raster ("Jacobi" style: every cell reads only the previous pass's values),
// so each pass maps one-to-one onto a GPU compute kernel. gpu.js runs the same passes
// in WGSL; this file is the CPU reference, in 32-bit floats like the GPU.
//
//   0. uplift    z += U dt (not outlets or fixed cells); the fill surface w moves along.
//   1. fill      w' = max(z, min(w, min_j (w_j + eps d_j)))       x fillSweeps
//                The depression-filling iteration of Planchon & Darboux (2001, Catena
//                46), lowered from above: where drainage is unresolved w starts at
//                +infinity, and after as many sweeps as the longest flow path w is z with
//                every depression filled to its spill point plus a gradient eps.
//   2. receiver  steepest descent on w.
//   3. area      A'_i = a_i + sum of A_j over neighbours j whose receiver is i
//                                                                   x areaSweeps
//                Area advances areaSweeps cells per iteration and is carried over
//                between iterations (amortised, as in Schott, Paris, Fournier, Guerin &
//                Galin 2023, ACM Trans. Graph. 42: terrain changes slowly compared with
//                how fast area settles).
//   4. erode     z'_i = (w_i + F z_r) / (1 + F),   F = K dt A_i^m / L_ir   x sweeps
//                The implicit update of fastscape.js, with the receiver's value z_r from
//                the previous sweep (the first guess is its value before this
//                iteration's uplift). Block Jacobi iteration towards the implicit
//                solution: each sweep is a weighted mean of the cell and its receiver,
//                so it is unconditionally stable, and the steady state is exactly the
//                implicit scheme's, K A^m S = U. A knickpoint moves at most `sweeps`
//                cells per iteration. Lakes (w > z) fill with sediment to their spill
//                surface.
//   5. creep     z' = z + D dt lap(z)  (explicit, 4 neighbours).
//   6. restart   the next fill starts from the eroded z where drainage resolved, and
//                from +infinity where it did not. (Starting from the old fill instead
//                refills fresh incision to the previous surface: a spurious slope floor
//                of U dt / L.)
//
// Outlet cells (the sea, designated river exits) hold their elevation and absorb water.
// Fixed cells (tile edges) hold their elevation, which comes from the parent; unless they
// are also outlets they are walls: no water flows into them.

export const OFF8 = [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, 1], [1, -1], [-1, -1]];
// q drains to c when q's receiver direction is the reverse of the direction c -> q
export const OPPOSITE = OFF8.map(([a, b]) => OFF8.findIndex(([c, d]) => c === -a && d === -b));
export const FLAG_OUTLET = 1, FLAG_FIXED = 2;
export const BIG = 3.0e38;
export const JACOBI_DEFAULTS = { m: 0.5, D: 0, eps: 1e-4, sweeps: 8, fillSweeps: 4, areaSweeps: 4 };

const f = Math.fround;

export class JacobiSolver {
  /**
   * @param {number} nx
   * @param {number} ny
   * @param {number} dx cell size (m)
   * @param {object} fields {z, U, K, P (runoff weight), flags} typed arrays of nx*ny
   */
  constructor(nx, ny, dx, { z, U, K, P, flags }) {
    const n = nx * ny;
    this.nx = nx; this.ny = ny; this.dx = dx; this.size = n;
    this.z = Float32Array.from(z);
    this.U = Float32Array.from(U);
    this.K = Float32Array.from(K);
    this.P = Float32Array.from(P);
    this.flags = Uint32Array.from(flags);
    this.w = new Float32Array(n);
    for (let i = 0; i < n; i++) this.w[i] = this.flags[i] & (FLAG_OUTLET | FLAG_FIXED) ? this.z[i] : BIG;
    this.w2 = new Float32Array(n);
    this.rec = new Int32Array(n).fill(-1);
    this.A = new Float32Array(n);
    this.A2 = new Float32Array(n);
    this.zr = new Float32Array(n);
    this.zn = new Float32Array(n);
    this.iterations = 0;
  }

  _blocked(q) { return (this.flags[q] & FLAG_FIXED) && !(this.flags[q] & FLAG_OUTLET); }

  /** One fill sweep from src into dst. */
  _fill(src, dst, eps) {
    const { nx, ny, dx, z, flags } = this;
    const step = [f(eps * dx), f(eps * f(dx * Math.SQRT2))];
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const c = j * nx + i;
        if (flags[c] & (FLAG_OUTLET | FLAG_FIXED)) { dst[c] = z[c]; continue; }
        let lo = src[c];
        for (let k = 0; k < 8; k++) {
          const a = i + OFF8[k][0], b = j + OFF8[k][1];
          if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
          const q = b * nx + a;
          if (this._blocked(q)) continue;
          const v = f(src[q] + step[k < 4 ? 0 : 1]);
          if (v < lo) lo = v;
        }
        dst[c] = Math.max(z[c], lo);
      }
    }
  }

  _receivers(w) {
    const { nx, ny, dx, flags, rec } = this;
    const diag = f(dx * Math.SQRT2);
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const c = j * nx + i;
        rec[c] = -1;
        if (flags[c] & (FLAG_OUTLET | FLAG_FIXED) || w[c] >= BIG) continue;
        let best = 0;
        for (let k = 0; k < 8; k++) {
          const a = i + OFF8[k][0], b = j + OFF8[k][1];
          if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
          const q = b * nx + a;
          if (this._blocked(q)) continue;
          const s = f(f(w[c] - w[q]) / (k < 4 ? dx : diag));
          if (s > best) { best = s; rec[c] = k; }
        }
      }
    }
  }

  _area(src, dst) {
    const { nx, ny, dx, P, rec } = this;
    const cellArea = f(dx * dx);
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const c = j * nx + i;
        let s = f(cellArea * P[c]);
        for (let k = 0; k < 8; k++) {
          const a = i + OFF8[k][0], b = j + OFF8[k][1];
          if (a < 0 || a >= nx || b < 0 || b >= ny) continue;
          const q = b * nx + a;
          if (rec[q] === OPPOSITE[k]) s = f(s + src[q]);
        }
        dst[c] = s;
      }
    }
  }

  /**
   * Run iterations.
   * @param {number} count
   * @param {{dt:number, m?:number, D?:number, eps?:number, sweeps?:number, fillSweeps?:number, areaSweeps?:number}} p
   */
  run(count, p) {
    const { dt, m, D, eps, sweeps, fillSweeps, areaSweeps } = { ...JACOBI_DEFAULTS, ...p };
    const { nx, ny, dx, z, U, K, flags, rec, zr, zn } = this;
    const n = nx * ny;
    const diag = f(dx * Math.SQRT2);
    for (let it = 0; it < count; it++) {
      // 0. uplift, carrying the fill surface along
      for (let c = 0; c < n; c++) {
        if (flags[c] & (FLAG_OUTLET | FLAG_FIXED)) continue;
        const up = f(U[c] * dt);
        z[c] = f(z[c] + up);
        if (this.w[c] < BIG) this.w[c] = f(this.w[c] + up);
      }
      // 1. fill (ping-pong; the result ends up in this.w)
      for (let s = 0; s < fillSweeps; s++) { this._fill(this.w, this.w2, eps); [this.w, this.w2] = [this.w2, this.w]; }
      const w = this.w;
      // 2. receivers
      this._receivers(w);
      // 3. area (ping-pong; the result ends up in this.A)
      for (let s = 0; s < areaSweeps; s++) { this._area(this.A, this.A2); [this.A, this.A2] = [this.A2, this.A]; }
      const A = this.A;
      // 4. erode in sweeps; zr holds each receiver's value from the sweep before
      for (let c = 0; c < n; c++) zr[c] = flags[c] & (FLAG_OUTLET | FLAG_FIXED) || rec[c] < 0 ? w[c] : f(w[c] - f(U[c] * dt));
      for (let s = 0; s < sweeps; s++) {
        for (let c = 0; c < n; c++) {
          if (flags[c] & (FLAG_OUTLET | FLAG_FIXED) || rec[c] < 0) { zn[c] = z[c]; continue; }
          const k = rec[c];
          const r = c + OFF8[k][0] + OFF8[k][1] * nx;
          const F = f(f(f(K[c] * dt) * f(Math.pow(A[c], m))) / (k < 4 ? dx : diag));
          zn[c] = f(f(w[c] + f(F * zr[r])) / f(1 + F));
        }
        if (s < sweeps - 1) zr.set(zn);
      }
      // 5. creep
      if (D > 0) {
        const kd = f((D * dt) / (dx * dx));
        for (let j = 0; j < ny; j++) {
          for (let i = 0; i < nx; i++) {
            const c = j * nx + i;
            if (flags[c] & (FLAG_OUTLET | FLAG_FIXED) || i === 0 || j === 0 || i === nx - 1 || j === ny - 1) { z[c] = zn[c]; continue; }
            z[c] = f(zn[c] + f(kd * f(f(f(zn[c - 1] + zn[c + 1]) + f(zn[c - nx] + zn[c + nx])) - f(4 * zn[c]))));
          }
        }
      } else {
        z.set(zn);
      }
      // 6. restart the fill
      for (let c = 0; c < n; c++) w[c] = flags[c] & (FLAG_OUTLET | FLAG_FIXED) || rec[c] >= 0 ? z[c] : BIG;
      this.iterations++;
    }
    return this;
  }

  /** Receiver as a cell index (-1: none), from the direction code. */
  receiverIndex(c) {
    const k = this.rec[c];
    return k < 0 ? -1 : c + OFF8[k][0] + OFF8[k][1] * this.nx;
  }

  receiverDistance(c) {
    return this.rec[c] < 4 ? this.dx : this.dx * Math.SQRT2;
  }
}
