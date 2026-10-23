// Water routing: depression filling, flow directions, and discharge.
//
// 1. Priority-flood (Barnes, Lehman & Mulla 2014, "Priority-flood: An optimal
//    depression-filling and watershed-labeling algorithm", Computers & Geosciences 62)
//    grows inward from the outlets, always expanding the lowest known cell. Any
//    cell reached from a higher one is a pit and is raised to just above it; the
//    raised surface `filled` is the water surface of lakes. A tiny epsilon makes
//    lake surfaces drain in a definite direction.
// 2. Each cell drains to its steepest-descent neighbour on the filled surface (single
//    flow direction).
//    Drainage AREA can optionally be spread over all lower neighbours (multiple flow
//    direction), which avoids grid-aligned streaks on planar hillslopes; incision and
//    sediment still follow the single steepest receiver.
// 3. Cells pop off the priority queue in increasing filled elevation, so every
//    cell appears after its receiver. That pop order is the "stack": iterate it
//    forwards to go downstream-first, backwards to go upstream-first.

export class MinHeap {
  constructor(capacity) {
    this.idx = new Int32Array(capacity);
    this.key = new Float64Array(capacity);
    this.n = 0;
  }
  push(i, k) {
    let p = this.n++;
    const { idx, key } = this;
    while (p > 0) {
      const parent = (p - 1) >> 1;
      if (key[parent] <= k) break;
      idx[p] = idx[parent];
      key[p] = key[parent];
      p = parent;
    }
    idx[p] = i;
    key[p] = k;
  }
  pop() {
    const { idx, key } = this;
    const top = idx[0];
    const last = --this.n;
    const li = idx[last];
    const lk = key[last];
    let p = 0;
    for (;;) {
      let c = 2 * p + 1;
      if (c >= last) break;
      if (c + 1 < last && key[c + 1] < key[c]) c++;
      if (key[c] >= lk) break;
      idx[p] = idx[c];
      key[p] = key[c];
      p = c;
    }
    idx[p] = li;
    key[p] = lk;
    return top;
  }
}

/** Lake surfaces rise by this much (m) per cell away from their outlet. */
export const FLAT_EPSILON = 1e-6;

/**
 * @param {object} mesh any mesh satisfying the contract in mesh.js
 * @param {Float64Array} z surface elevation
 * @param {Uint8Array} outlet 1 where water leaves the domain (base level / sea)
 * @param {Float64Array} weight per-cell runoff multiplier (precipitation / 1 m/yr)
 * @param {object} out preallocated arrays: filled, receiver, receiverDist, order, area,
 *   visited (scratch). Optional: out.flowExponent > 0 switches drainage area to multiple
 *   flow direction; out.randomExponent > 0 with out.rng draws random receivers;
 *   out.inflow (Float64Array, m^2) adds drainage area entering from outside the domain.
 */
export function routeFlow(mesh, z, outlet, weight, out) {
  const n = mesh.size;
  const D = mesh.degree;
  const nb = mesh.neighbors;
  const dist = mesh.dist;
  const { filled, receiver, order, area } = out;
  const rdist = out.receiverDist;
  const visited = out.visited;
  visited.fill(0);
  const heap = out.heap || (out.heap = new MinHeap(n));
  heap.n = 0;

  let anyOutlet = false;
  for (let i = 0; i < n; i++) {
    if (outlet[i]) {
      filled[i] = z[i];
      receiver[i] = -1;
      visited[i] = 1;
      heap.push(i, z[i]);
      anyOutlet = true;
    }
  }
  if (!anyOutlet) throw new Error('routeFlow: no outlet cells');

  let count = 0;
  while (heap.n > 0) {
    const c = heap.pop();
    order[count++] = c;
    const fc = filled[c];
    for (let k = 0; k < D; k++) {
      const j = nb[D * c + k];
      if (j < 0 || visited[j]) continue;
      visited[j] = 1;
      const raised = fc + FLAT_EPSILON;
      filled[j] = z[j] > raised ? z[j] : raised;
      heap.push(j, filled[j]);
    }
  }
  if (count < n) throw new Error('routeFlow: some cells are not connected to any outlet');

  // Receivers: the steepest descent on the filled surface (drop / distance), or with
  // randomExponent a random lower neighbour weighted by slope^p (the hex analogue of
  // the "Rho8" method of Fairfield & Leymarie 1991, Water Resour. Res. 27). Re-drawn
  // every step, random receivers remove the bias towards the lattice directions while
  // keeping a single receiver for the implicit solver.
  const rng = out.rng;
  const rp = out.randomExponent ?? 0;
  for (let i = 0; i < n; i++) {
    if (outlet[i]) { if (rdist) rdist[i] = 0; continue; }
    const fi = filled[i];
    let best = -1;
    let bestS = 0;
    let bestD = 0;
    let sum = 0;
    for (let k = 0; k < D; k++) {
      const e = D * i + k;
      const j = nb[e];
      if (j < 0 || filled[j] >= fi) continue;
      const S = (fi - filled[j]) / dist[e];
      if (S > bestS) { bestS = S; best = j; bestD = dist[e]; }
      if (rp > 0) sum += Math.pow(S, rp);
    }
    if (rng && rp > 0) {
      let u = rng.next() * sum;
      for (let k = 0; k < D; k++) {
        const e = D * i + k;
        const j = nb[e];
        if (j < 0 || filled[j] >= fi) continue;
        u -= Math.pow((fi - filled[j]) / dist[e], rp);
        if (u <= 0) { best = j; bestD = dist[e]; break; }
      }
    }
    receiver[i] = best;
    if (rdist) rdist[i] = bestD;
  }

  const inflow = out.inflow;
  for (let i = 0; i < n; i++) area[i] = weight[i] * mesh.area[i] + (inflow ? inflow[i] : 0);
  const p = out.flowExponent ?? 0;
  if (!(p > 0)) {
    // single flow direction: everything goes to the receiver
    for (let s = n - 1; s >= 0; s--) {
      const i = order[s];
      const r = receiver[i];
      if (r >= 0) area[r] += area[i];
    }
    return;
  }
  // Multiple flow direction (Freeman 1991, Computers & Geosciences 17): area is shared
  // among all lower neighbours in proportion to slope^p. Large p tends to single flow.
  // Every lower neighbour has lower filled elevation, so it comes earlier in `order`
  // and the upstream-first sweep stays valid.
  const w = out.mfdWeights || (out.mfdWeights = new Float64Array(D));
  for (let s = n - 1; s >= 0; s--) {
    const i = order[s];
    if (receiver[i] < 0) continue;
    const fi = filled[i];
    let sum = 0;
    for (let k = 0; k < D; k++) {
      const e = D * i + k;
      const j = nb[e];
      const d = j >= 0 ? fi - filled[j] : 0;
      w[k] = d > 0 ? Math.pow(d / dist[e], p) : 0;
      sum += w[k];
    }
    const ai = area[i] / sum;
    for (let k = 0; k < D; k++) if (w[k] > 0) area[nb[D * i + k]] += ai * w[k];
  }
}
