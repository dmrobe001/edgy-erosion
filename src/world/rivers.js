// Drainage graph extraction and the plausibility statistics of a drainage network.
//
// Input is a drainage raster: each cell's receiver (-1 for outlets and undrained cells),
// the distance to it, drainage area, and elevation. Output:
//
//   river graph   nodes at channel heads, confluences and mouths; edges are the channel
//                 reaches between them, with their cell paths, discharge and Strahler
//                 order (the design's RiverGraph)
//   Hack's law    main-stream length L (longest flow path to the divide) against
//                 drainage area A: L = c A^h, h ~ 0.5-0.6 in real networks (Hack 1957)
//   slope-area    channel slope S against area: S = k A^-theta, theta ~ 0.4-0.6
//                 (Flint 1974; the concavity index)
//
// Both fits are least squares on log-binned medians (equal weight per decade of area,
// not per cell, so the many small channels do not swamp the fit).

/**
 * Cells in upstream-to-downstream order (every cell before its receiver): a
 * topological order of the drainage forest (Kahn's algorithm on donor counts).
 * @param {Int32Array} receiver
 */
export function upstreamOrder(receiver) {
  const n = receiver.length;
  const donors = new Int32Array(n);
  for (let i = 0; i < n; i++) if (receiver[i] >= 0) donors[receiver[i]]++;
  const order = new Int32Array(n);
  let head = 0, tail = 0;
  for (let i = 0; i < n; i++) if (donors[i] === 0) order[tail++] = i;
  while (head < tail) {
    const i = order[head++];
    const r = receiver[i];
    if (r >= 0 && --donors[r] === 0) order[tail++] = r;
  }
  if (tail !== n) throw new Error('drainage has a cycle');
  return order;
}

/**
 * Longest flow path from each cell up to the divide (m): Hack's main-stream length.
 * @param {Int32Array} receiver
 * @param {ArrayLike<number>} dist distance to the receiver
 */
export function mainStreamLength(receiver, dist, order = upstreamOrder(receiver)) {
  const L = new Float64Array(receiver.length);
  for (const i of order) {
    const r = receiver[i];
    if (r >= 0) L[r] = Math.max(L[r], L[i] + dist[i]);
  }
  return L;
}

/** Least-squares slope and intercept of y on x. */
function fit(xs, ys) {
  const n = xs.length;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let k = 0; k < n; k++) { sx += xs[k]; sy += ys[k]; sxx += xs[k] * xs[k]; sxy += xs[k] * ys[k]; }
  const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const intercept = (sy - slope * sx) / n;
  let ss = 0, st = 0;
  const my = sy / n;
  for (let k = 0; k < n; k++) { ss += (ys[k] - intercept - slope * xs[k]) ** 2; st += (ys[k] - my) ** 2; }
  return { slope, intercept, r2: 1 - ss / (st || 1), points: n };
}

/** Medians of y in bins of log10(x) `perDecade` to a decade, keeping bins with >= minCount. */
function binnedMedians(x, y, perDecade = 4, minCount = 5) {
  const bins = new Map();
  for (let k = 0; k < x.length; k++) {
    if (!(x[k] > 0 && y[k] > 0)) continue;
    const b = Math.floor(Math.log10(x[k]) * perDecade);
    if (!bins.has(b)) bins.set(b, []);
    bins.get(b).push(y[k]);
  }
  const xs = [], ys = [];
  for (const [b, list] of [...bins].sort((p, q) => p[0] - q[0])) {
    if (list.length < minCount) continue;
    list.sort((p, q) => p - q);
    xs.push((b + 0.5) / perDecade);
    ys.push(Math.log10(list[Math.floor(list.length / 2)]));
  }
  return { xs, ys };
}

/**
 * Hack's law and the slope-area relation on the channel cells of a drainage raster.
 * @param {object} d {receiver, dist, area, z, mask?} mask: 1 for cells to include
 * @param {{channelArea:number, flatSlope?:number}} o channels drain at least channelArea
 *   (m^2); cells with slopes at or below flatSlope are lakes and flats (filled
 *   depressions), left out of the slope-area fit as is usual
 */
export function drainageStatistics(d, { channelArea, flatSlope = 0 }) {
  const { receiver, dist, area, z } = d;
  const order = upstreamOrder(receiver);
  const L = mainStreamLength(receiver, dist, order);
  const aH = [], lH = [], aS = [], sS = [];
  for (let i = 0; i < receiver.length; i++) {
    const r = receiver[i];
    if (r < 0 || area[i] < channelArea || (d.mask && !d.mask[i])) continue;
    aH.push(area[i]); lH.push(L[i]);
    const s = (z[i] - z[r]) / dist[i];
    if (s > flatSlope) { aS.push(area[i]); sS.push(s); }
  }
  const h = binnedMedians(aH, lH);
  const s = binnedMedians(aS, sS);
  const hack = fit(h.xs, h.ys);
  const sa = fit(s.xs, s.ys);
  return {
    hack: { exponent: hack.slope, coefficient: 10 ** hack.intercept, r2: hack.r2, bins: hack.points },
    slopeArea: { concavity: -sa.slope, steepness: 10 ** sa.intercept, r2: sa.r2, bins: sa.points },
    channelCells: aH.length,
    areaRange: aH.reduce((r, v) => [Math.min(r[0], v), Math.max(r[1], v)], aH.length ? [Infinity, -Infinity] : [0, 0]),
  };
}

/**
 * The river graph: channel cells (area >= channelArea) as a forest of reaches.
 * @param {object} d {receiver, dist, area, z, nx, dx, runoff?} runoff in m/yr per cell (for discharge)
 * @returns {{nodes: {cell:number, i:number, j:number, z:number, kind:string}[],
 *            edges: {from:number, to:number, path:number[], discharge:number, area:number, order:number, length:number}[],
 *            crossings: object[]}}
 */
export function extractRiverGraph(d, { channelArea, runoff = 0.5 }) {
  const { receiver, dist, area, z, nx } = d;
  const n = receiver.length;
  const isChan = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (area[i] >= channelArea) isChan[i] = 1;
  const upCount = new Int32Array(n);
  for (let i = 0; i < n; i++) if (isChan[i] && receiver[i] >= 0 && isChan[receiver[i]]) upCount[receiver[i]]++;
  const order = upstreamOrder(receiver);
  // Strahler order
  const strahler = new Int32Array(n), top = new Int32Array(n), topCount = new Int32Array(n);
  for (const i of order) {
    if (!isChan[i]) continue;
    strahler[i] = top[i] === 0 ? 1 : topCount[i] > 1 ? top[i] + 1 : top[i];
    const r = receiver[i];
    if (r < 0 || !isChan[r]) continue;
    if (strahler[i] > top[r]) { top[r] = strahler[i]; topCount[r] = 1; } else if (strahler[i] === top[r]) topCount[r]++;
  }
  // nodes: heads (no channel donors), confluences (2+), mouths (receiver not a channel)
  const nodeOf = new Int32Array(n).fill(-1);
  const nodes = [];
  const addNode = (c, kind) => {
    if (nodeOf[c] >= 0) return nodeOf[c];
    nodeOf[c] = nodes.length;
    nodes.push({ cell: c, i: c % nx, j: Math.floor(c / nx), z: z[c], kind });
    return nodeOf[c];
  };
  for (const c of order) {
    if (!isChan[c]) continue;
    const r = receiver[c];
    if (upCount[c] === 0 && (r < 0 || !isChan[r])) continue; // an isolated cell is not a river
    if (upCount[c] === 0) addNode(c, 'head');
    else if (upCount[c] > 1) addNode(c, 'confluence');
    if (r < 0 || !isChan[r]) addNode(c, 'mouth');
  }
  // edges: walk down from every node to the next node
  const edges = [];
  const secPerYr = 3.15576e7;
  for (const nd of nodes) {
    const c0 = nd.cell;
    if (nd.kind === 'mouth' && (receiver[c0] < 0 || !isChan[receiver[c0]])) continue;
    const path = [c0];
    let c = receiver[c0], length = dist[c0];
    while (c >= 0 && isChan[c] && nodeOf[c] < 0) { path.push(c); length += dist[c]; c = receiver[c]; }
    if (c < 0 || !isChan[c]) continue;
    path.push(c);
    const last = path[path.length - 2];
    edges.push({ from: nodeOf[c0], to: nodeOf[c], path, area: area[last], discharge: (area[last] * runoff) / secPerYr, order: strahler[last], length });
  }
  return { nodes, edges, crossings: [] };
}
