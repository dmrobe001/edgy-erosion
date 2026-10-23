// Water routing: depression filling, flow directions, and discharge.
//
// 1. Priority-flood (Barnes, Lehman & Mulla 2014, "Priority-flood: An optimal
//    depression-filling and watershed-labeling algorithm", Computers & Geosciences 62)
//    grows inward from the outlets, always expanding the lowest known cell. Any
//    cell reached from a higher one is a pit and is raised to just above it; the
//    raised surface `filled` is the water surface of lakes. A tiny epsilon makes
//    lake surfaces drain in a definite direction.
// 2. Each cell drains to its lowest neighbour on the filled surface (single flow
//    direction; on a hex grid all six neighbours are equidistant, so lowest = steepest).
//    Drainage AREA can optionally be spread over all lower neighbours (multiple flow
//    direction), which avoids grid-aligned streaks on planar hillslopes; incision and
//    sediment still follow the single steepest receiver.
// 3. Cells pop off the priority queue in increasing filled elevation, so every
//    cell appears after its receiver. That pop order is the "stack": iterate it
//    forwards to go downstream-first, backwards to go upstream-first.

class MinHeap {
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
 * @param {import('./hexgrid.js').HexGrid} grid
 * @param {Float64Array} z surface elevation
 * @param {Uint8Array} outlet 1 where water leaves the domain (base level / sea)
 * @param {Float64Array} weight per-cell runoff multiplier (precipitation / 1 m/yr)
 * @param {object} out preallocated arrays: filled, receiver, order, area, visited (scratch);
 *   optional out.flowExponent > 0 switches drainage area to multiple flow direction
 */
export function routeFlow(grid, z, outlet, weight, out) {
  const n = grid.size;
  const nb = grid.neighbors;
  const { filled, receiver, order, area } = out;
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
    for (let k = 0; k < 6; k++) {
      const j = nb[6 * c + k];
      if (j < 0 || visited[j]) continue;
      visited[j] = 1;
      const raised = fc + FLAT_EPSILON;
      filled[j] = z[j] > raised ? z[j] : raised;
      heap.push(j, filled[j]);
    }
  }

  const rng = out.rng;
  const rp = out.randomExponent ?? 0;
  for (let i = 0; i < n; i++) {
    if (outlet[i]) continue;
    let best = -1;
    let bestZ = filled[i];
    for (let k = 0; k < 6; k++) {
      const j = nb[6 * i + k];
      if (j >= 0 && filled[j] < bestZ) { bestZ = filled[j]; best = j; }
    }
    if (rng && rp > 0) {
      // Random receiver among the lower neighbours, chosen with probability proportional
      // to drop^rp (the hex analogue of the "Rho8" method of Fairfield & Leymarie 1991,
      // Water Resour. Res. 27). Re-drawn every step, it removes the bias towards the six
      // lattice directions while keeping a single receiver for the implicit solver.
      const fi = filled[i];
      let sum = 0;
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j >= 0 && filled[j] < fi) sum += Math.pow(fi - filled[j], rp);
      }
      let u = rng.next() * sum;
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j >= 0 && filled[j] < fi) {
          u -= Math.pow(fi - filled[j], rp);
          if (u <= 0) { best = j; break; }
        }
      }
    }
    receiver[i] = best;
  }

  const a = grid.cellArea;
  for (let i = 0; i < n; i++) area[i] = weight[i] * a;
  const p = out.flowExponent ?? 0;
  if (!(p > 0)) {
    // single flow direction: everything goes to the steepest receiver
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
  const w = new Float64Array(6);
  for (let s = n - 1; s >= 0; s--) {
    const i = order[s];
    if (receiver[i] < 0) continue;
    const fi = filled[i];
    let sum = 0;
    for (let k = 0; k < 6; k++) {
      const j = nb[6 * i + k];
      const d = j >= 0 ? fi - filled[j] : 0;
      w[k] = d > 0 ? Math.pow(d, p) : 0;
      sum += w[k];
    }
    const ai = area[i] / sum;
    for (let k = 0; k < 6; k++) if (w[k] > 0) area[nb[6 * i + k]] += ai * w[k];
  }
}
