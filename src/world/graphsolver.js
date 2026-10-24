// The parallel stream-power solver of jacobi.js on an irregular graph (pointmesh.js),
// with sediment transport and deposition.
//
// The passes per iteration are those of jacobi.js (uplift, fill, receivers, area,
// erosion sweeps, creep, fill restart), with neighbours from compressed rows instead of
// eight fixed offsets, plus sediment:
//
//   dz/dt = U - K A^m S + sum_k G_k Qs_k / A        (Davy & Lague 2009, J. Geophys. Res.
//                                                   114; implicit form: Yuan, Braun,
//                                                   Guerit et al. 2019, J. Geophys. Res.
//                                                   Earth Surf. 124)
//
// Qs_k is the flux of grain class k (m^3/yr) arriving at a node from upstream and from
// outside the window. A node lays down the share 1 - exp(-G_k a / A) of what arrives (a:
// the node's own area, A: its drainage area), so a class travels a characteristic
// distance set by G_k: small G (fine sediment) travels far, large G (gravel) settles
// soon. What leaves is the rest plus what the node erodes, split among the classes by
// the rock's fractions. Like drainage area, Qs is gathered over a few sweeps per
// iteration and carried between iterations. Deposition enters the erosion update
// explicitly: z' = (w + dt D + F z_r) / (1 + F), D the deposition rate (m/yr).
//
// All arithmetic in 32-bit floats, as on the GPU.

export const FLAG_OUTLET = 1, FLAG_FIXED = 2;
export const BIG = 3.0e38;
export const CLASSES = 3; // fine (mostly suspended), sand, gravel
export const GRAPH_DEFAULTS = { m: 0.5, D: 0, eps: 1e-5, sweeps: 8, fillSweeps: 4, areaSweeps: 4, sedSweeps: 4, G: [0.05, 0.6, 2.5] };

const f = Math.fround;

export class GraphSolver {
  /**
   * @param {object} g graph: {n, off, nbr, dist, lapw, area}
   * @param {object} fields {z, U, K, P (runoff weight), flags, inflow? (m^2 of area entering),
   *   sedIn? (CLASSES x n, m^3/yr entering), frac? (CLASSES x n, share of eroded rock per class)}
   */
  constructor(g, fields) {
    const n = g.n;
    this.g = g; this.n = n;
    this.z = Float32Array.from(fields.z);
    this.U = Float32Array.from(fields.U);
    this.K = Float32Array.from(fields.K);
    this.P = Float32Array.from(fields.P);
    this.flags = Uint32Array.from(fields.flags);
    this.inflow = fields.inflow ? Float32Array.from(fields.inflow) : new Float32Array(n);
    this.sedIn = fields.sedIn ? Float32Array.from(fields.sedIn) : new Float32Array(CLASSES * n);
    this.frac = fields.frac ? Float32Array.from(fields.frac) : Float32Array.from({ length: CLASSES * n }, (_, k) => (k < n ? 0.4 : k < 2 * n ? 0.4 : 0.2));
    this.w = new Float32Array(n);
    for (let i = 0; i < n; i++) this.w[i] = this.flags[i] & (FLAG_OUTLET | FLAG_FIXED) ? this.z[i] : BIG;
    this.w2 = new Float32Array(n);
    this.rec = new Int32Array(n).fill(-1);
    this.recDist = new Float32Array(n);
    this.A = new Float32Array(n); this.A2 = new Float32Array(n);
    this.Qs = new Float32Array(CLASSES * n); this.Qs2 = new Float32Array(CLASSES * n);
    this.prod = new Float32Array(CLASSES * n); // m^3/yr eroded per class, last iteration
    this.dep = new Float32Array(n);            // deposition rate, m/yr, last sediment pass
    this.zr = new Float32Array(n); this.zn = new Float32Array(n);
    this.iterations = 0;
  }

  _blocked(q) { return (this.flags[q] & FLAG_FIXED) && !(this.flags[q] & FLAG_OUTLET); }

  run(count, p) {
    const { dt, m, D, eps, sweeps, fillSweeps, areaSweeps, sedSweeps, G } = { ...GRAPH_DEFAULTS, ...p };
    const { n, z, U, K, P, flags, rec, recDist, zr, zn, inflow, sedIn, frac, prod, dep } = this;
    const { off, nbr, dist, lapw, area } = this.g;
    const held = (i) => flags[i] & (FLAG_OUTLET | FLAG_FIXED);
    for (let it = 0; it < count; it++) {
      // 0. uplift, carrying the fill surface along
      for (let i = 0; i < n; i++) {
        if (held(i)) continue;
        const up = f(U[i] * dt);
        z[i] = f(z[i] + up);
        if (this.w[i] < BIG) this.w[i] = f(this.w[i] + up);
      }
      // 1. fill
      for (let s = 0; s < fillSweeps; s++) {
        const src = this.w, dst = this.w2;
        for (let i = 0; i < n; i++) {
          if (held(i)) { dst[i] = z[i]; continue; }
          let lo = src[i];
          for (let e = off[i]; e < off[i + 1]; e++) {
            const j = nbr[e];
            if (this._blocked(j)) continue;
            const v = f(src[j] + f(eps * dist[e]));
            if (v < lo) lo = v;
          }
          dst[i] = Math.max(z[i], lo);
        }
        this.w = dst; this.w2 = src;
      }
      const w = this.w;
      // 2. receivers: steepest descent
      for (let i = 0; i < n; i++) {
        rec[i] = -1;
        if (held(i) || w[i] >= BIG) continue;
        let best = 0;
        for (let e = off[i]; e < off[i + 1]; e++) {
          const j = nbr[e];
          if (this._blocked(j)) continue;
          const s = f(f(w[i] - w[j]) / dist[e]);
          if (s > best) { best = s; rec[i] = j; recDist[i] = dist[e]; }
        }
      }
      // 3. drainage area, gathered from donors
      for (let s = 0; s < areaSweeps; s++) {
        const src = this.A, dst = this.A2;
        for (let i = 0; i < n; i++) {
          let a = f(f(area[i] * P[i]) + inflow[i]);
          for (let e = off[i]; e < off[i + 1]; e++) { const j = nbr[e]; if (rec[j] === i) a = f(a + src[j]); }
          dst[i] = a;
        }
        this.A = dst; this.A2 = src;
      }
      const A = this.A;
      // 4. sediment: what arrives, what is laid down, what leaves
      for (let s = 0; s < sedSweeps; s++) {
        const src = this.Qs, dst = this.Qs2;
        for (let i = 0; i < n; i++) {
          let d = 0;
          for (let k = 0; k < CLASSES; k++) {
            let q = sedIn[k * n + i];
            for (let e = off[i]; e < off[i + 1]; e++) { const j = nbr[e]; if (rec[j] === i) q = f(q + src[k * n + j]); }
            const laid = held(i) ? 0 : f(q * f(1 - Math.exp(-G[k] * area[i] / Math.max(A[i], area[i]))));
            d = f(d + laid);
            dst[k * n + i] = f(f(q - laid) + prod[k * n + i]);
          }
          dep[i] = f(d / area[i]);
        }
        this.Qs = dst; this.Qs2 = src;
      }
      // 5. erode in sweeps (deposition enters explicitly)
      for (let i = 0; i < n; i++) zr[i] = held(i) || rec[i] < 0 ? w[i] : f(w[i] - f(U[i] * dt));
      for (let s = 0; s < sweeps; s++) {
        for (let i = 0; i < n; i++) {
          if (held(i) || rec[i] < 0) { zn[i] = z[i]; continue; }
          const F = f(f(f(K[i] * dt) * f(Math.pow(A[i], m))) / recDist[i]);
          zn[i] = f(f(f(w[i] + f(dt * dep[i])) + f(F * zr[rec[i]])) / f(1 + F));
        }
        if (s < sweeps - 1) zr.set(zn);
      }
      // what each node eroded this iteration (detachment, K A^m S), split by rock fractions
      for (let i = 0; i < n; i++) {
        const r = rec[i];
        const E = held(i) || r < 0 ? 0 : Math.max(0, f(f(f(K[i] * f(Math.pow(A[i], m))) * f(zn[i] - zn[r])) / recDist[i]));
        for (let k = 0; k < CLASSES; k++) prod[k * n + i] = f(f(E * area[i]) * frac[k * n + i]);
      }
      // 6. creep (explicit, cotangent weights)
      if (D > 0) {
        for (let i = 0; i < n; i++) {
          if (held(i)) { z[i] = zn[i]; continue; }
          let s = 0;
          for (let e = off[i]; e < off[i + 1]; e++) s = f(s + f(lapw[e] * f(zn[nbr[e]] - zn[i])));
          z[i] = f(zn[i] + f(f(D * dt / area[i]) * s));
        }
      } else z.set(zn);
      // 7. restart the fill from the eroded surface (or from above where unresolved)
      for (let i = 0; i < n; i++) w[i] = held(i) || rec[i] >= 0 ? z[i] : BIG;
      this.iterations++;
    }
    return this;
  }

  /** Sediment flux of class k leaving node i (m^3/yr). */
  sediment(k, i) { return this.Qs[k * this.n + i]; }
}
