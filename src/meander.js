// Meandering rivers as moving centre lines.
//
// A river wider than a few grid cells is better represented as a curve than as cells:
// its centre line is a polyline whose nodes carry width W, depth D and water-surface
// elevation. Bends migrate by the kinematic model of Howard & Knutson (1984, Water
// Resour. Res. 20), itself a simplification of Ikeda, Parker & Sawai (1981, J. Fluid
// Mech. 112): the bank at a node erodes at a rate set by the local curvature and, with
// more weight, by the curvature upstream, because the flow's momentum carries the
// fastest water to the outer bank some distance downstream of each bend's apex:
//
//   R1(s) = Omega C(s) W + Gamma * integral_0^inf C(s - xi) W G(xi) dxi / integral G
//   G(xi) = exp(-2 Cf xi / D),   Omega = -1, Gamma = 2.5, Cf ~ 0.01
//
// Nodes move sideways towards the outer bank by (migration rate) x R1 x (bank
// erodibility). This makes bends grow, skew downstream and migrate down-valley. When
// a loop's neck closes (two parts of the line come within about one width of each
// other) the loop is cut off and left behind as an oxbow lake.
//
// Bank erodibility comes from a callback, so the river can be confined by valley walls,
// stopped at the sea, or slowed on steep (bedrock) reaches.

import { Rng, fbm } from './rng.js';

/** Resample a polyline at roughly even spacing ds(node), interpolating every attribute. */
export function resample(nodes, spacing) {
  const out = [nodes[0]];
  let carry = 0;
  for (let k = 1; k < nodes.length; k++) {
    const a = nodes[k - 1], b = nodes[k];
    const L = Math.hypot(b.x - a.x, b.y - a.y);
    if (L === 0) continue;
    let s = spacing(a) - carry;
    while (s < L) {
      const t = s / L;
      out.push(lerpNode(a, b, t));
      s += spacing(out[out.length - 1]);
    }
    carry = L - (s - spacing(out[out.length - 1]));
  }
  const last = nodes[nodes.length - 1];
  const prev = out[out.length - 1];
  if (Math.hypot(last.x - prev.x, last.y - prev.y) < 0.3 * spacing(last) && out.length > 1) out[out.length - 1] = last;
  else out.push(last);
  return out;
}

function lerpNode(a, b, t) {
  const n = {};
  for (const key of Object.keys(a)) n[key] = typeof a[key] === 'number' ? a[key] + (b[key] - a[key]) * t : a[key];
  return n;
}

/** Chaikin corner cutting, keeping the end points: removes the zigzag of grid paths. */
export function chaikin(nodes, iterations = 3) {
  let p = nodes;
  for (let it = 0; it < iterations; it++) {
    if (p.length < 3) return p;
    const q = [p[0]];
    for (let k = 0; k < p.length - 1; k++) {
      q.push(lerpNode(p[k], p[k + 1], 0.25), lerpNode(p[k], p[k + 1], 0.75));
    }
    q.push(p[p.length - 1]);
    p = q;
  }
  return p;
}

/** Signed curvature (1/m) at each node; > 0 turning left (counter-clockwise in x-right, y-up). */
export function curvature(nodes) {
  const n = nodes.length;
  const c = new Float64Array(n);
  for (let k = 1; k < n - 1; k++) {
    const ax = nodes[k].x - nodes[k - 1].x, ay = nodes[k].y - nodes[k - 1].y;
    const bx = nodes[k + 1].x - nodes[k].x, by = nodes[k + 1].y - nodes[k].y;
    const ang = Math.atan2(ax * by - ay * bx, ax * bx + ay * by);
    c[k] = ang / (0.5 * (Math.hypot(ax, ay) + Math.hypot(bx, by)) || 1);
  }
  return c;
}

/**
 * Migrate a river centre line.
 * @param {object[]} nodes [{x, y, w, d, zw}] from upstream to downstream; ends stay fixed
 * @param {object} o
 * @param {number} [o.iterations]
 * @param {number} [o.step] node displacement per iteration, in widths, per unit of R1
 * @param {(x:number, y:number, node:object) => number} [o.erodibility] bank erodibility 0..1
 * @param {number} [o.seed]
 * @param {number} [o.heterogeneity] amplitude (natural-log units) of random bank strength
 * @param {number} [o.kick] size of the random bank collapse each iteration, in widths
 * @returns {{nodes: object[], oxbows: {nodes: object[], iteration: number}[], sinuosity: number}}
 */
export function migrate(nodes, { iterations = 600, step = 0.5, erodibility = () => 1, seed = 1, cf = 0.025, maxSinuosity = 2.2, heterogeneity = 1.0, kick = 0.3 } = {}) {
  const rng = new Rng(seed);
  const spacing = (nd) => 0.5 * nd.w;
  let line = resample(nodes, spacing);
  // seed bends: a few random waves along the line, 8-30 widths long, 0.3 widths high
  // (the unstable band of the model; noise at node scale would only be damped)
  const waves = Array.from({ length: 4 }, () => ({ k: (2 * Math.PI) / (8 + 22 * rng.next()), ph: 2 * Math.PI * rng.next() }));
  let arc = 0;
  for (let k = 1; k < line.length - 1; k++) {
    const a = line[k - 1], b = line[k + 1];
    arc += Math.hypot(line[k].x - a.x, line[k].y - a.y) / line[k].w;
    const tx = b.x - a.x, ty = b.y - a.y, L = Math.hypot(tx, ty) || 1;
    let r = 0;
    for (const wv of waves) r += Math.sin(wv.k * arc + wv.ph);
    r *= 0.15 * line[k].w * Math.min(1, k / 8, (line.length - 1 - k) / 8);
    line[k] = { ...line[k], x: line[k].x - (ty / L) * r, y: line[k].y + (tx / L) * r };
  }
  const startLen = pathLength(line);
  const straight = Math.hypot(line[line.length - 1].x - line[0].x, line[line.length - 1].y - line[0].y) || 1;
  const oxbows = [];
  for (let it = 0; it < iterations; it++) {
    line = resample(line, spacing);
    const n = line.length;
    if (n < 5) break;
    const c0 = curvature(line);
    // smooth the curvature (a few [1 2 1] passes): damps node-scale wiggles, which the
    // explicit update would otherwise amplify, and lets the step be large
    let c = c0;
    for (let pass = 0; pass < 4; pass++) {
      const s2 = new Float64Array(n);
      for (let k = 1; k < n - 1; k++) s2[k] = 0.25 * c[k - 1] + 0.5 * c[k] + 0.25 * c[k + 1];
      c = s2;
    }
    // upstream-weighted curvature (the convolution, truncated where G < 1%)
    const move = new Float64Array(n);
    for (let k = 1; k < n - 1; k++) {
      const nd = line[k];
      const ds = 0.5 * nd.w;
      const alpha = (2 * cf) / Math.max(nd.d, 0.1);
      let sum = 0, norm = 0;
      for (let j = 1; j <= k; j++) {
        const g = Math.exp(-alpha * j * ds);
        if (g < 0.01) break;
        sum += c[k - j] * line[k - j].w * g;
        norm += g;
      }
      const r1 = -1 * c[k] * nd.w + (norm > 0 ? (2.5 * sum) / norm : 0);
      const taper = Math.min(1, k / 4, (n - 1 - k) / 4); // fixed ends
      // banks are not uniform (old fills, clay plugs, tree roots): a fixed random field of
      // bank strength, a few widths across, keeps seeding new bends along the whole river
      const het = Math.exp(heterogeneity * fbm(nd.x, nd.y, { wavelength: 4 * nd.w, octaves: 3, seed: seed + 71 }));
      // capped at 0.15 widths per iteration, so that no two stretches of river can jump
      // across each other between neck-cutoff checks
      const m = step * nd.w * r1 * taper * het * erodibility(nd.x, nd.y, nd);
      move[k] = Math.max(-0.15 * nd.w, Math.min(0.15 * nd.w, m));
    }
    // move towards the outer bank: for a left turn (c > 0) that is the right-hand side
    const next = line.map((nd) => ({ ...nd }));
    for (let k = 1; k < n - 1; k++) {
      if (move[k] === 0) continue;
      const tx = line[k + 1].x - line[k - 1].x, ty = line[k + 1].y - line[k - 1].y;
      const L = Math.hypot(tx, ty) || 1;
      // right-hand normal of (tx, ty) in x-right/y-up axes is (ty, -tx)
      next[k].x += (ty / L) * move[k];
      next[k].y += (-tx / L) * move[k];
    }
    // bank collapses: a random small sideways kick somewhere along the river each
    // iteration keeps bends forming upstream (the model's instability is convective:
    // bends grow as they travel downstream, so without new disturbances the upper
    // reach would straighten)
    {
      const k0 = 1 + Math.floor(rng.next() * (n - 2));
      const amp = kick * line[k0].w * (rng.next() < 0.5 ? -1 : 1);
      const span = 3;
      for (let k = Math.max(1, k0 - 3 * span); k < Math.min(n - 1, k0 + 3 * span); k++) {
        const g = amp * Math.exp(-(((k - k0) / span) ** 2));
        const tx = line[k + 1].x - line[k - 1].x, ty = line[k + 1].y - line[k - 1].y;
        const L = Math.hypot(tx, ty) || 1;
        next[k].x += (ty / L) * g * erodibility(next[k].x, next[k].y, next[k]);
        next[k].y += (-tx / L) * g * erodibility(next[k].x, next[k].y, next[k]);
      }
    }
    line = next;
    // neck cutoffs
    line = cutoffs(line, oxbows, it);
    if (pathLength(line) / straight > maxSinuosity) break;
  }
  // round off the corners the last cutoffs left (two passes of [1 2 1] on the positions)
  for (let pass = 0; pass < 2; pass++) {
    const prev = line.map((nd) => ({ ...nd }));
    for (let k = 1; k < line.length - 1; k++) {
      line[k].x = 0.25 * prev[k - 1].x + 0.5 * prev[k].x + 0.25 * prev[k + 1].x;
      line[k].y = 0.25 * prev[k - 1].y + 0.5 * prev[k].y + 0.25 * prev[k + 1].y;
    }
  }
  return { nodes: line, oxbows, sinuosity: pathLength(line) / straight, lengthening: pathLength(line) / startLen };
}

function pathLength(line) {
  let L = 0;
  for (let k = 1; k < line.length; k++) L += Math.hypot(line[k].x - line[k - 1].x, line[k].y - line[k - 1].y);
  return L;
}

/** Cut loops whose necks have closed to within a width; the loops become oxbows. */
function cutoffs(line, oxbows, iteration) {
  let n = line.length;
  let wmax = 0;
  for (const nd of line) wmax = Math.max(wmax, nd.w);
  const cell = wmax;
  const grid = new Map();
  const keyOf = (x, y) => `${Math.floor(x / cell)},${Math.floor(y / cell)}`;
  for (let k = 0; k < n; k++) {
    const kk = keyOf(line[k].x, line[k].y);
    if (!grid.has(kk)) grid.set(kk, []);
    grid.get(kk).push(k);
  }
  for (let k = 0; k < n; k++) {
    const nd = line[k];
    const gx = Math.floor(nd.x / cell), gy = Math.floor(nd.y / cell);
    let cut = -1;
    for (let dx = -1; dx <= 1 && cut < 0; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        const list = grid.get(`${gx + dx},${gy + dy}`);
        if (!list) continue;
        for (const j of list) {
          // far apart along the river (more than ~3 widths of channel between them) but close in space
          if (j <= k + 12) continue;
          if (Math.hypot(line[j].x - nd.x, line[j].y - nd.y) < 0.9 * nd.w && j > cut) cut = j;
        }
      }
    }
    if (cut > 0) {
      oxbows.push({ nodes: line.slice(k, cut + 1), iteration });
      return line.slice(0, k + 1).concat(line.slice(cut));
    }
  }
  return line;
}
