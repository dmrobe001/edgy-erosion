// Detail: from a landscape simulation to a map at walking scale (metres per pixel).
//
// Below a few tens of metres the grid simulation is the wrong tool: every river is one
// cell wide and runs along the lattice, and beaches, bars and levees are narrower than
// a cell. This module keeps the simulation as the large-scale truth (relief, drainage,
// soils, vegetation, rock, sea level, wave exposure) and adds what it cannot resolve:
//
//   rivers     traced from the simulation's drainage network as curves, given a width
//              and depth from their discharge (downstream hydraulic geometry, Leopold &
//              Maddock 1953: W = aW Q^0.5, D = aD Q^0.35 at bankfull), smoothed out of
//              the lattice, and, where wide and gentle, left to meander (meander.js),
//              which leaves oxbow lakes behind;
//   channels   carved into the surface with an asymmetric cross-section (deep at the
//              outer bank of a bend, a point bar on the inside), low-flow water that
//              exposes the bars, natural levees and a flattened floodplain around
//              meandering rivers;
//   beaches    on sandy, gentle, wave-exposed coasts: an equilibrium shoreface offshore
//              (h = A y^(2/3), Dean 1977), a beach face rising to a berm, and dunes
//              behind, their width growing with wave exposure;
//   texture    micro-relief scaled to the local relief, and land cover (forest with
//              crowns, grass, scrub, bare rock, cliffs, marsh) from the simulation's
//              vegetation, soil and slope.
// [heuristic] All of this is procedural detail conditioned on the simulation, chosen to
// look plausible, not computed physics (the meander model aside).
//
// Everything a pixel gets depends only on its world position (noise is seeded by world
// coordinates; the river network is built once for the whole patch), so tiles rendered
// separately join up exactly: render any window, at any pixel size, and pan.

import { interpolator } from './nest.js';
import { fbm } from './rng.js';
import { chaikin, resample, migrate, curvature } from './meander.js';

const SECONDS_PER_YEAR = 3.156e7;
const FLOODPLAIN = 254; // provisional cover while rendering

/** Land-cover classes of a rendered tile. */
export const COVER = {
  deepWater: 0, shallowWater: 1, river: 2, lake: 3, wetSand: 4, sand: 5, gravel: 6, marsh: 7,
  duneGrass: 8, grass: 9, scrub: 10, forest: 11, rock: 12, cliff: 13,
};

export const DETAIL_DEFAULTS = {
  seed: 1,
  minWidth: 0.6,          // m: narrowest stream drawn (bankfull width)
  runoff: 0.5,            // share of precipitation that becomes river flow
  bankfullFactor: 4,      // bankfull discharge / mean discharge
  widthCoefficient: 4,    // aW in W = aW Q^0.5 (m, Q in m^3/s)
  depthCoefficient: 0.4,  // aD in D = aD Q^0.35
  meanderMinWidth: 6,     // m: rivers at least this wide may meander...
  meanderMaxSlope: 0.012, // ...where gentler than this (generous: coarse valleys are steep)
  belt: 7,                // meander belt half-width, in channel widths: the floodplain the detail opens up
  warp: 1.5,              // bend the coarse grid's straight valleys, by this many coarse spacings
  meanderIterations: 2500,
  sinuosity: 2.7,         // meandering stops at this sinuosity, if neck cutoffs have not held it lower
  oxbowsShown: 5,         // newest cut-off bends drawn per river: half as lakes, half as marsh
  wiggle: 0.3,            // sideways noise on smaller streams, in coarse cell spacings
  microRelief: 0.4,       // m: smallest micro-relief amplitude
  beachWidth: [15, 60],   // m: beach width on sheltered .. fully exposed sandy coast
  duneWidth: [20, 70],    // m
  duneHeight: [1.5, 6],   // m
  shoreface: 0.12,        // Dean's A (m^(1/3)), fine sand
};

/** A field of a hex simulation resampled on a regular grid and read back bicubically. */
class Raster {
  constructor(sim, field) {
    const g = sim.grid;
    const b = g.bounds;
    this.h = g.dx;
    this.x0 = b.xmin;
    this.y0 = b.ymin;
    this.nx = Math.ceil((b.xmax - b.xmin) / this.h) + 1;
    this.ny = Math.ceil((b.ymax - b.ymin) / this.h) + 1;
    this.data = new Float64Array(this.nx * this.ny);
    const at = interpolator(g, field);
    for (let j = 0; j < this.ny; j++) {
      for (let i = 0; i < this.nx; i++) {
        // clamp into the cell-centre hull so edge nodes take the nearest values
        const x = Math.min(Math.max(this.x0 + i * this.h, g.x[0]), g.x[g.size - 1]);
        const y = Math.min(Math.max(this.y0 + j * this.h, g.y[0]), g.y[g.size - 1]);
        let v = at(x, y);
        if (!Number.isFinite(v)) v = field[Math.max(0, g.cellAt(x, y))];
        this.data[j * this.nx + i] = v;
      }
    }
  }

  /** Catmull-Rom bicubic value at world (x, y): smooth, so the hexagons disappear. */
  sample(x, y) {
    if (this.warp) { const w = this.warp(x, y); x += w[0]; y += w[1]; }
    return this._sample(x, y);
  }

  /** As sample, at coordinates already warped. */
  _sample(x, y) {
    const fx = (x - this.x0) / this.h, fy = (y - this.y0) / this.h;
    const ix = Math.floor(fx), iy = Math.floor(fy);
    const tx = fx - ix, ty = fy - iy;
    const { nx, ny, data } = this;
    let v = 0;
    for (let m = -1; m <= 2; m++) {
      const jy = Math.min(ny - 1, Math.max(0, iy + m));
      const wy = cr(ty, m);
      let row = 0;
      for (let k = -1; k <= 2; k++) {
        const jx = Math.min(nx - 1, Math.max(0, ix + k));
        row += cr(tx, k) * data[jy * nx + jx];
      }
      v += wy * row;
    }
    return v;
  }

  /** Bilinear value (for fields where overshoot would matter). */
  linear(x, y) {
    if (this.warp) { const w = this.warp(x, y); x += w[0]; y += w[1]; }
    return this._linear(x, y);
  }

  /** As linear, at coordinates already warped. */
  _linear(x, y) {
    const fx = Math.min(this.nx - 1.001, Math.max(0, (x - this.x0) / this.h));
    const fy = Math.min(this.ny - 1.001, Math.max(0, (y - this.y0) / this.h));
    const ix = Math.floor(fx), iy = Math.floor(fy), tx = fx - ix, ty = fy - iy;
    const d = this.data, n = this.nx;
    return (1 - ty) * ((1 - tx) * d[iy * n + ix] + tx * d[iy * n + ix + 1]) + ty * ((1 - tx) * d[(iy + 1) * n + ix] + tx * d[(iy + 1) * n + ix + 1]);
  }

  /** Slope magnitude from the bilinear field (central differences at the grid spacing). */
  slope(x, y) {
    if (this.warp) { const w = this.warp(x, y); x += w[0]; y += w[1]; }
    return Math.hypot(...this._gradient(x, y));
  }

  /** Gradient [d/dx, d/dy] at coordinates already warped. */
  _gradient(x, y) {
    const h = this.h;
    return [(this._linear(x + h, y) - this._linear(x - h, y)) / (2 * h), (this._linear(x, y + h) - this._linear(x, y - h)) / (2 * h)];
  }
}

/** Catmull-Rom weight of tap k (-1..2) at fraction t. */
function cr(t, k) {
  const t2 = t * t, t3 = t2 * t;
  switch (k) {
    case -1: return 0.5 * (-t3 + 2 * t2 - t);
    case 0: return 0.5 * (3 * t3 - 5 * t2 + 2);
    case 1: return 0.5 * (-3 * t3 + 4 * t2 + t);
    default: return 0.5 * (t3 - t2);
  }
}

function smooth01(a, b, x) {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
}

export class Detail {
  /**
   * @param {import('./simulation.js').Simulation} sim a hex-lattice landscape simulation
   * @param {object} [options] see DETAIL_DEFAULTS
   */
  constructor(sim, options = {}) {
    if (!sim.grid.lattice) throw new Error('detailing needs a hex-lattice simulation');
    this.sim = sim;
    this.o = { ...DETAIL_DEFAULTS, ...options };
    this.seed = (this.o.seed ?? sim.params.seed) | 0;
    this.bounds = { ...sim.grid.bounds };
    const g0 = sim.grid;
    this.hull = { x0: g0.x[0], x1: g0.x[g0.size - 1], y0: g0.y[0], y1: g0.y[g0.size - 1] };
    this.sea = sim.seaLevel;
    const t0 = Date.now();
    this.z = new Raster(sim, sim.z);
    this.water = new Raster(sim, sim.flow.filled);
    this.soil = new Raster(sim, sim.soil);
    this.veg = new Raster(sim, sim.veg);
    this.waves = new Raster(sim, sim.waveExposureField);
    const lakeDepth = Float64Array.from(sim.flow.filled, (f, i) => f - sim.z[i]);
    this.lake = new Raster(sim, lakeDepth);
    // A gentle warp of where the coarse fields are read: straight, lattice-parallel valleys
    // bend (amplitude warp x coarse spacing, wavelength 12 spacings). Rivers get the
    // inverse warp, so they stay in their valleys.
    const A = this.o.warp * sim.grid.dx, Lw = 12 * sim.grid.dx, ws = this.seed;
    this.warp = A > 0 ? (x, y) => [A * fbm(x, y, { wavelength: Lw, octaves: 2, seed: ws + 301 }), A * fbm(x, y, { wavelength: Lw, octaves: 2, seed: ws + 302 })] : null;
    for (const r of [this.z, this.water, this.soil, this.veg, this.waves, this.lake]) r.warp = this.warp;
    this.rockColor = sim.lith.color.map((c) => { const v = parseInt(c.slice(1), 16); return [(v >> 16) & 255, (v >> 8) & 255, v & 255]; });
    this._buildRivers();
    this.buildTime = Date.now() - t0;
  }

  /** Bankfull hydraulics from a drainage area (m^2, precipitation-weighted). */
  hydraulics(area) {
    const o = this.o;
    const Q = (o.bankfullFactor * o.runoff * area) / SECONDS_PER_YEAR; // m^3/s
    return { Q, w: o.widthCoefficient * Math.sqrt(Q), d: o.depthCoefficient * Math.pow(Q, 0.35) };
  }

  _buildRivers() {
    const sim = this.sim, g = sim.grid, o = this.o;
    const { receiver, filled, order } = sim.flow;
    const n = g.size;
    // drainage area along the single-receiver tree (the simulation may spread area over
    // several lower neighbours, which would start wide "rivers" with no stream feeding them)
    const area = new Float64Array(n);
    for (let i = 0; i < n; i++) area[i] = sim.precipWeight[i] * g.area[i] + sim.inflowArea[i];
    for (let s2 = n - 1; s2 >= 0; s2--) { const i = order[s2]; if (receiver[i] >= 0) area[receiver[i]] += area[i]; }
    // channel cells: wide enough to draw, and not sea or open edge
    const minQ = (o.minWidth / o.widthCoefficient) ** 2;
    const minArea = (minQ * SECONDS_PER_YEAR) / (o.bankfullFactor * o.runoff);
    const isChannel = new Uint8Array(n);
    for (let i = 0; i < n; i++) if (area[i] >= minArea && !sim.outlet[i]) isChannel[i] = 1;
    // each channel cell's main donor: the upstream channel cell with the most area
    const main = new Int32Array(n).fill(-1);
    for (let i = 0; i < n; i++) {
      if (!isChannel[i]) continue;
      const r = receiver[i];
      if (r >= 0 && isChannel[r] && (main[r] < 0 || area[i] > area[main[r]])) main[r] = i;
    }
    // polylines start at channel heads (no channel upstream) and run downstream while
    // they are the main donor of the next cell; a tributary ends where it joins a bigger
    // river, whose own line continues through the junction
    const lines = [];
    const lineOf = new Int32Array(n).fill(-1);
    const starts = [];
    for (let i = 0; i < n; i++) if (isChannel[i] && main[i] < 0) starts.push(i);
    starts.sort((a, b) => area[b] - area[a]); // big rivers first
    for (const s of starts) {
      const cells = [s];
      let c = s, joins = -1;
      const estuary = () => {
        // at the sea: carry the channel on downhill across the shallows (an estuary, or a
        // channel across the shoreface) until the water is deep enough to hide it
        // Follow the lowest neighbour that leads away from where the river met the sea,
        // so it runs out along a flooded inlet rather than along the shore.
        const hd = this.hydraulics(area[c]).d;
        const sx = g.x[c], sy = g.y[c];
        for (let step = 0; step < 60 && this.sea != null && sim.z[c] > this.sea - 2 * hd - 0.5; step++) {
          const here = Math.hypot(g.x[c] - sx, g.y[c] - sy);
          let lo = -1;
          for (let q = 0; q < 6; q++) {
            const j = g.neighbors[6 * c + q];
            if (j < 0 || cells.includes(j) || Math.hypot(g.x[j] - sx, g.y[j] - sy) <= here) continue;
            if (lo < 0 || sim.z[j] < sim.z[lo]) lo = j;
          }
          if (lo < 0 || sim.z[lo] > sim.z[c] + 0.3) break;
          cells.push(lo);
          c = lo;
        }
      };
      for (;;) {
        const r = receiver[c];
        if (r < 0) { if (sim.outlet[c] && !sim.boundaryOutlet[c]) estuary(); break; }
        cells.push(r); // the receiver ends this line (sea, edge or junction) or continues it
        if (!isChannel[r]) {
          if (sim.outlet[r] && !sim.boundaryOutlet[r]) { c = r; estuary(); }
          break;
        }
        if (main[r] !== c) { joins = r; break; }
        c = r;
      }
      const id = lines.length;
      for (const cell of cells) if (lineOf[cell] < 0) lineOf[cell] = id;
      lines.push({ cells, joins });
    }
    // geometry: coarse cell centres, smoothed out of the lattice, with a little sideways
    // noise, and hydraulics from drainage area
    const rivers = [];
    for (const L of lines) {
      if (L.cells.length < 2) continue;
      let zw = Infinity, a = 0;
      let nodes = L.cells.map((c) => {
        zw = Math.min(zw, filled[c]); // water surface never rises downstream...
        a = Math.max(a, area[c]);     // ...and discharge never falls (sea cells carry none of their own)
        const h = this.hydraulics(a);
        return { x: g.x[c], y: g.y[c], w: h.w, d: h.d, q: h.Q, zw };
      });
      nodes = chaikin(nodes, 3);
      const coarse = g.dx;
      nodes = resample(nodes, () => coarse / 3);
      // sideways noise along the line, fading where the river is wide
      let arc = 0;
      for (let k = 1; k < nodes.length - 1; k++) {
        arc += Math.hypot(nodes[k].x - nodes[k - 1].x, nodes[k].y - nodes[k - 1].y);
        const tx = nodes[k + 1].x - nodes[k - 1].x, ty = nodes[k + 1].y - nodes[k - 1].y;
        const len = Math.hypot(tx, ty) || 1;
        const amp = o.wiggle * coarse * Math.min(1, k / 3, (nodes.length - 1 - k) / 3) / (1 + nodes[k].w / coarse);
        const off = amp * fbm(arc, L.cells[0], { wavelength: 2.5 * coarse, octaves: 3, seed: this.seed + 5 });
        nodes[k] = { ...nodes[k], x: nodes[k].x + (ty / len) * off, y: nodes[k].y - (tx / len) * off };
      }
      if (this.warp) {
        nodes = nodes.map((nd) => {
          const w1 = this.warp(nd.x, nd.y);
          const w2 = this.warp(nd.x - w1[0], nd.y - w1[1]);
          return { ...nd, x: nd.x - w2[0], y: nd.y - w2[1] };
        });
      }
      rivers.push({ nodes, joins: L.joins, first: L.cells[0], oxbows: [], meandering: false });
    }
    // meander the wide, gentle ones, inside a belt around the valley line
    const sea = this.sea;
    const b = this.bounds;
    rivers.forEach((r, idx) => {
      const nodes = r.nodes;
      const wmax = Math.max(...nodes.map((nd) => nd.w));
      if (wmax < o.meanderMinWidth) return;
      // local slope of the water surface, over about five widths
      for (let k = 0; k < nodes.length; k++) {
        const a = nodes[Math.max(0, k - 3)], c = nodes[Math.min(nodes.length - 1, k + 3)];
        let len = 0;
        for (let m = Math.max(1, k - 2); m <= Math.min(nodes.length - 1, k + 3); m++) len += Math.hypot(nodes[m].x - nodes[m - 1].x, nodes[m].y - nodes[m - 1].y);
        nodes[k].s = (a.zw - c.zw) / Math.max(len, 1);
      }
      if (!nodes.some((nd) => nd.w >= o.meanderMinWidth && nd.s < o.meanderMaxSlope)) return;
      // the valley line the belt is centred on, in a spatial hash
      const cell = 2 * wmax;
      const hash = new Map();
      for (const nd of nodes) {
        const key = `${Math.floor(nd.x / cell)},${Math.floor(nd.y / cell)}`;
        if (!hash.has(key)) hash.set(key, []);
        hash.get(key).push(nd);
      }
      const reach = Math.ceil((o.belt * wmax) / cell) + 1;
      const distToValley = (x, y) => {
        const gx = Math.floor(x / cell), gy = Math.floor(y / cell);
        let d2 = Infinity;
        for (let dx = -reach; dx <= reach; dx++) {
          for (let dy = -reach; dy <= reach; dy++) {
            const list = hash.get(`${gx + dx},${gy + dy}`);
            if (list) for (const v of list) d2 = Math.min(d2, (v.x - x) ** 2 + (v.y - y) ** 2);
          }
        }
        return Math.sqrt(d2);
      };
      const erodibility = (x, y, nd) => {
        if (x < b.xmin + nd.w || x > b.xmax - nd.w || y < b.ymin + nd.w || y > b.ymax - nd.w) return 0;
        if (nd.w < o.meanderMinWidth) return 0;
        if (sea != null && this.z.sample(x, y) < sea - 0.5) return 0; // not out into open water
        // inside a belt around the valley line, whatever the coarse terrain says: meandering
        // stands in for the lateral erosion the grid model lacks, and the belt becomes a
        // floodplain when rendered [heuristic]
        const belt = o.belt * nd.w;
        return (1 - smooth01(0.6 * belt, belt, distToValley(x, y))) * (1 - smooth01(0.5 * o.meanderMaxSlope, o.meanderMaxSlope, nd.s));
      };
      const res = migrate(nodes, { iterations: o.meanderIterations, maxSinuosity: o.sinuosity, erodibility, seed: this.seed + 31 * idx });
      r.nodes = res.nodes;
      // keep only the parts of cut-off bends that the live channel has not since crossed
      const live = new Map();
      for (const nd of r.nodes) {
        const key = `${Math.floor(nd.x / cell)},${Math.floor(nd.y / cell)}`;
        if (!live.has(key)) live.set(key, []);
        live.get(key).push(nd);
      }
      const nearLive = (p) => {
        const gx = Math.floor(p.x / cell), gy = Math.floor(p.y / cell);
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) {
          const list = live.get(`${gx + dx},${gy + dy}`);
          if (list) for (const v of list) if ((v.x - p.x) ** 2 + (v.y - p.y) ** 2 < (1.3 * Math.max(v.w, p.w)) ** 2) return true;
        }
        return false;
      };
      r.oxbows = [];
      for (const ox of res.oxbows) {
        let piece = [];
        for (const nd of ox.nodes) {
          if (nearLive(nd)) { if (piece.length >= 6) r.oxbows.push({ nodes: piece, iteration: ox.iteration }); piece = []; } else piece.push(nd);
        }
        if (piece.length >= 6) r.oxbows.push({ nodes: piece, iteration: ox.iteration });
      }
      r.meandering = true;
      r.sinuosity = res.sinuosity;
    });
    // tributaries re-attach to wherever their river now runs
    for (const r of rivers) {
      if (r.joins < 0) continue;
      // the nearest node of any bigger river to the tributary's end
      const end = r.nodes[r.nodes.length - 1];
      let best = null, bd = Infinity;
      for (const h of rivers) {
        if (h === r || h.nodes[0].w < end.w * 0.5) continue;
        for (const nd of h.nodes) {
          const d = (nd.x - end.x) ** 2 + (nd.y - end.y) ** 2;
          if (d < bd) { bd = d; best = nd; }
        }
      }
      if (!best || bd > (o.belt * best.w + 3 * this.sim.grid.dx) ** 2) continue;
      // the creek winds across the floodplain to wherever the river now runs
      const L = Math.sqrt(bd);
      const steps = Math.max(1, Math.ceil(L / (this.sim.grid.dx / 3)));
      const ux = (best.x - end.x) / (L || 1), uy = (best.y - end.y) / (L || 1);
      const extra = [];
      for (let q = 1; q <= steps; q++) {
        const t = q / steps;
        const off = q < steps ? 0.15 * L * Math.sin(Math.PI * t) * fbm(t * 3, r.first, { wavelength: 1, octaves: 2, seed: this.seed + 9 }) : 0;
        extra.push({ ...end, x: end.x + ux * L * t - uy * off, y: end.y + uy * L * t + ux * off, zw: Math.min(end.zw, best.zw) });
      }
      r.nodes = r.nodes.concat(extra);
    }
    for (const r of rivers) {
      const c = curvature(r.nodes);
      r.nodes.forEach((nd, k) => { nd.c = c[k]; });
      r.wmax = Math.max(...r.nodes.map((nd) => nd.w));
    }
    // small rivers first, so bigger ones win where they overlap
    rivers.sort((a, b2) => a.wmax - b2.wmax);
    this.rivers = rivers;
    this.oxbowCount = rivers.reduce((s, r) => s + r.oxbows.length, 0);
  }

  /**
   * Render a window of the map.
   * @param {object} o
   * @param {number} o.x0 world x of the window's left edge (m)
   * @param {number} o.y0 world y of the window's top edge (m)
   * @param {number} o.width pixels
   * @param {number} o.height pixels
   * @param {number} o.pixel metres per pixel
   * @returns {{x0, y0, width, height, pixel, z: Float32Array, level: Float32Array, cover: Uint8Array, rock: Uint8Array, shadeZ: Float32Array}}
   */
  render({ x0, y0, width, height, pixel }) {
    const o = this.o;
    // beaches need the distance to the shore, so work on a padded window
    const pad = Math.ceil((o.beachWidth[1] + o.duneWidth[1] + 40) / pixel) + 2;
    const W = width + 2 * pad, H = height + 2 * pad;
    const N = W * H;
    const X0 = x0 - pad * pixel, Y0 = y0 - pad * pixel;
    const z = new Float32Array(N), level = new Float32Array(N).fill(-Infinity);
    const cover = new Uint8Array(N).fill(255), soil = new Float32Array(N), veg = new Float32Array(N), slopeC = new Float32Array(N), aspectC = new Float32Array(N);
    const sea = this.sea;
    const coarse = this.sim.grid.dx;
    const seed = this.seed;
    // 1. base surface: smooth interpolation plus micro-relief scaled to the local relief
    // (octaves down to about two pixels: finer would only alias)
    const octaves = Math.max(2, Math.min(7, Math.round(Math.log2((1.5 * coarse) / (2 * pixel))) + 1));
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const k = j * W + i;
        const x = X0 + (i + 0.5) * pixel, y = Y0 + (j + 0.5) * pixel;
        let xw = x, yw = y;
        if (this.warp) { const w = this.warp(x, y); xw += w[0]; yw += w[1]; }
        const gr = this.z._gradient(xw, yw);
        const s = Math.hypot(gr[0], gr[1]);
        const zc = this.z._sample(xw, yw);
        // micro-relief, scaled to the local relief; low flats by the sea stay flat (tidal
        // flats, estuaries) rather than breaking up into speckled islands
        let amp = o.microRelief + 0.25 * coarse * Math.min(s, 1);
        if (sea != null) amp = Math.min(amp, 0.1 + amp * smooth01(0.5, 5, zc - sea));
        z[k] = zc + amp * fbm(x, y, { wavelength: 1.5 * coarse, octaves, seed: seed + 101 });
        soil[k] = Math.max(0, this.soil._linear(xw, yw));
        if (x < this.hull.x0 || x > this.hull.x1 || y < this.hull.y0 || y > this.hull.y1) soil[k] = Math.max(soil[k], 0.1); // ragged hex edge
        veg[k] = Math.min(1, Math.max(0, this.veg._linear(xw, yw)));
        slopeC[k] = s;
        aspectC[k] = (gr[0] - gr[1]) * Math.SQRT1_2; // > 0 facing south-west (uphill to the north-east)
        // lakes the simulation filled
        if (this.lake._linear(xw, yw) > 0.3) {
          const lv = this.water._linear(xw, yw);
          if (lv > z[k]) { level[k] = lv; cover[k] = COVER.lake; }
        }
      }
    }
    // 2. beaches and shorefaces on sandy coasts
    if (sea != null) this._beaches({ z, level, cover, soil, slopeC, W, H, X0, Y0, pixel });
    // 3. rivers and oxbows (river: 0..1 nearness to a river, for riparian vegetation)
    const river = this._rivers({ z, level, cover, W, H, X0, Y0, pixel });
    // 4. land cover and the sea
    const tide = this.sim.params.coast?.tidalRange ?? 1;
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const k = j * W + i;
        if (sea != null && z[k] < sea && level[k] < sea) level[k] = sea;
        if (level[k] > z[k]) {
          if (cover[k] === 255 || cover[k] === FLOODPLAIN || cover[k] === COVER.wetSand || cover[k] === COVER.sand) cover[k] = level[k] - z[k] > 3 ? COVER.deepWater : COVER.shallowWater;
          continue;
        }
        if (cover[k] === FLOODPLAIN) {
          const x = X0 + (i + 0.5) * pixel, y = Y0 + (j + 0.5) * pixel;
          const f = fbm(x, y, { wavelength: 90, octaves: 3, seed: seed + 41 });
          cover[k] = f > 0.25 ? COVER.forest : f < -0.35 ? COVER.marsh : COVER.grass;
          continue;
        }
        if (cover[k] !== 255) continue; // already sand, gravel, marsh...
        const x = X0 + (i + 0.5) * pixel, y = Y0 + (j + 0.5) * pixel;
        // fine-scale slope, from this tile's surface
        const zl = z[j * W + Math.max(0, i - 1)], zr = z[j * W + Math.min(W - 1, i + 1)];
        const zu = z[Math.max(0, j - 1) * W + i], zd = z[Math.min(H - 1, j + 1) * W + i];
        const s = Math.hypot(zr - zl, zd - zu) / (2 * pixel);
        if (sea != null && z[k] < sea + tide / 2 && soil[k] > 0.2) { cover[k] = COVER.wetSand; continue; }
        if (s > 1.2) { cover[k] = COVER.cliff; continue; }
        // bare rock on steep ground, or thin soil on slopes; gentle bare valley floors are
        // alluvium the coarse grid could not hold, so they get cover like their neighbours
        if (s > 0.8 || (soil[k] < 0.08 && slopeC[k] > 0.08)) { cover[k] = COVER.rock; continue; }
        // vegetation: the simulation's cover, broken into stands and clearings, wetter
        // (more trees) along rivers and drier on steep, sunny (south-west facing) slopes
        const patch = 0.35 * fbm(x, y, { wavelength: 260, octaves: 4, seed: seed + 7 });
        const riparian = river[k] > 0 ? 0.25 * river[k] : 0;
        const aspect = 0.6 * Math.max(0, aspectC[k]);
        const v = veg[k] + patch + riparian - aspect;
        cover[k] = v > 0.6 ? COVER.forest : v > 0.3 ? COVER.grass : COVER.scrub;
      }
    }
    // rock type at the surface, for colour
    const rock = new Uint8Array(N);
    const g = this.sim.grid;
    for (let j = 0; j < H; j += 1) {
      for (let i = 0; i < W; i++) {
        // clamped into the cell-centre hull, so the ragged hex edge takes its neighbours' rock
        let px = X0 + (i + 0.5) * pixel, py = Y0 + (j + 0.5) * pixel;
        if (this.warp) { const w = this.warp(px, py); px += w[0]; py += w[1]; }
        const c = g.cellAt(Math.min(Math.max(px, g.x[0]), g.x[g.size - 1]), Math.min(Math.max(py, g.y[0]), g.y[g.size - 1]));
        // the strata at this pixel's own elevation, so outcrops follow the contours
        rock[j * W + i] = c >= 0 ? this.sim.lith.rockAt(c, z[j * W + i] - this.sim.upliftCum[c] - 0.5) : 0;
      }
    }
    // crop the padding
    const crop = (a, T) => { const out = new T(width * height); for (let j = 0; j < height; j++) out.set(a.subarray((j + pad) * W + pad, (j + pad) * W + pad + width), j * width); return out; };
    return {
      x0, y0, width, height, pixel,
      z: crop(z, Float32Array), level: crop(level, Float32Array), cover: crop(cover, Uint8Array), rock: crop(rock, Uint8Array),
    };
  }

  _beaches({ z, level, cover, soil, slopeC, W, H, X0, Y0, pixel }) {
    const o = this.o;
    const sea = this.sea;
    const N = W * H;
    const wet = new Uint8Array(N);
    for (let k = 0; k < N; k++) wet[k] = z[k] < sea && level[k] === -Infinity ? 1 : 0;
    const dLand = edt(wet, W, H, 1); // land pixels: distance to the sea (pixels)
    const dSea = edt(wet, W, H, 0);  // sea pixels: distance to land
    const cp = this.sim.params.coast || {};
    const berm = (cp.tidalRange ?? 1) / 2 + 0.6 * (cp.waveHeight ?? 1);
    const maxReach = (o.beachWidth[1] + o.duneWidth[1]) / pixel;
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const k = j * W + i;
        const isSea = wet[k];
        const dpx = Math.sqrt(isSea ? dSea[k] : dLand[k]);
        if (!isSea && dpx > maxReach) continue;
        if (isSea && dpx * pixel > 500) continue;
        const x = X0 + (i + 0.5) * pixel, y = Y0 + (j + 0.5) * pixel;
        // sandy where the simulation has sediment on gentle ground
        const sand = smooth01(0.1, 1.0, this.soil.linear(x, y)) * (1 - smooth01(0.1, 0.3, slopeC[k]));
        if (sand <= 0.02) continue;
        const expo = Math.min(1, Math.max(0.2, this.waves.linear(x, y) * 1.5));
        const yd = dpx * pixel;
        if (isSea) {
          // equilibrium shoreface, blending back into the simulated sea floor offshore
          const dean = sea - o.shoreface * Math.pow(yd, 2 / 3);
          const w = sand * (1 - smooth01(150, 500, yd));
          z[k] = z[k] + (dean - z[k]) * w;
          continue;
        }
        const bw = o.beachWidth[0] + (o.beachWidth[1] - o.beachWidth[0]) * expo;
        const dw = (o.duneWidth[0] + (o.duneWidth[1] - o.duneWidth[0]) * expo) * sand;
        const zBeach = sea + Math.min(berm, 0.07 * yd);
        if (yd < bw) {
          z[k] = z[k] + (zBeach - z[k]) * sand;
          if (sand > 0.4) cover[k] = COVER.sand;
        } else if (yd < bw + dw) {
          const t = (yd - bw) / dw;
          const ridges = Math.max(0, fbm(x, y, { wavelength: 45, octaves: 3, seed: this.seed + 211 }) + 0.25);
          const dh = (o.duneHeight[0] + (o.duneHeight[1] - o.duneHeight[0]) * expo) * ridges * Math.sin(Math.PI * Math.min(1, t * 1.3));
          const target = sea + berm + dh;
          const wgt = sand * (1 - smooth01(0.6, 1, t));
          z[k] = z[k] + (target - z[k]) * wgt;
          if (wgt > 0.3) cover[k] = dh > 1 && fbm(x, y, { wavelength: 25, octaves: 2, seed: this.seed + 17 }) > 0 ? COVER.sand : COVER.duneGrass;
        }
      }
    }
  }

  _rivers({ z, level, cover, W, H, X0, Y0, pixel }) {
    const N = W * H;
    // nearest river influence per pixel, in units of that river's half-width
    const best = new Float32Array(N).fill(Infinity);
    const half = new Float32Array(N), depth = new Float32Array(N), zw = new Float32Array(N);
    const side = new Float32Array(N), curv = new Float32Array(N), reach = new Float32Array(N), kind = new Uint8Array(N);
    const xmax = X0 + W * pixel, ymax = Y0 + H * pixel;
    const stamp = (nodes, F, k0, scaleD, age) => {
      for (let s = 1; s < nodes.length; s++) {
        const a = nodes[s - 1], b = nodes[s];
        const hw = 0.5 * Math.max(a.w, b.w);
        const R = hw * (1 + F) + pixel;
        const bx0 = Math.min(a.x, b.x) - R, bx1 = Math.max(a.x, b.x) + R;
        const by0 = Math.min(a.y, b.y) - R, by1 = Math.max(a.y, b.y) + R;
        if (bx1 < X0 || bx0 > xmax || by1 < Y0 || by0 > ymax) continue;
        const i0 = Math.max(0, Math.floor((bx0 - X0) / pixel)), i1 = Math.min(W - 1, Math.ceil((bx1 - X0) / pixel));
        const j0 = Math.max(0, Math.floor((by0 - Y0) / pixel)), j1 = Math.min(H - 1, Math.ceil((by1 - Y0) / pixel));
        const ex = b.x - a.x, ey = b.y - a.y, L2 = ex * ex + ey * ey || 1;
        for (let j = j0; j <= j1; j++) {
          const py = Y0 + (j + 0.5) * pixel;
          for (let i = i0; i <= i1; i++) {
            const px = X0 + (i + 0.5) * pixel;
            let t = ((px - a.x) * ex + (py - a.y) * ey) / L2;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
            const cx = a.x + t * ex, cy = a.y + t * ey;
            const d = Math.hypot(px - cx, py - cy);
            const hwt = 0.5 * (a.w + t * (b.w - a.w));
            const u = d / hwt;
            const k = j * W + i;
            if (u >= 1 + F || u >= best[k]) continue;
            best[k] = u;
            half[k] = hwt;
            depth[k] = scaleD * (a.d + t * (b.d - a.d));
            zw[k] = a.zw + t * (b.zw - a.zw);
            // which side of the centre line (left = +1), for the asymmetric section
            side[k] = (ex * (py - a.y) - ey * (px - a.x)) > 0 ? 1 : -1;
            curv[k] = a.c + t * (b.c - a.c);
            reach[k] = F;
            kind[k] = k0 + (age > 0 ? Math.min(age, 3) : 0);
          }
        }
      }
    };
    // floodplains of meandering rivers first: distance to the nearest centre line by a
    // feature transform on a grid extended by the belt width (so tiles join up)
    const near = new Float32Array(N);
    this._floodplains({ z, cover, near, W, H, X0, Y0, pixel });
    for (const r of this.rivers) {
      const F = 1.2;
      // oxbows: the newest are lakes, older ones marsh; older still have filled in and
      // are not drawn (on the ground they would be faint scroll bars in the floodplain)
      r.oxbows.forEach((ox, n) => {
        const ageRank = r.oxbows.length - n; // 1 = newest
        if (ageRank > this.o.oxbowsShown) return;
        const age = ageRank <= this.o.oxbowsShown / 2 ? 1 : 2;
        stamp(ox.nodes, 0.4, 1, age === 1 ? 0.7 : 0.35, age);
      });
      stamp(r.nodes, F, 0, 1, 0);
    }
    const sea = this.sea;
    for (let k = 0; k < N; k++) {
      const u = best[k];
      if (u === Infinity) continue;
      if (sea != null && z[k] < sea - 0.3) continue; // out at sea the channel is under water anyway
      const D = depth[k], zk = zw[k];
      if (kind[k] === 0) {
        if (u < 1) {
          // asymmetric bed: deeper towards the outer bank of a bend, a point bar inside.
          // For a left turn (curvature > 0) the outer bank is on the right (side = -1).
          const bend = Math.max(-1, Math.min(1, curv[k] * half[k] * 10));
          const outer = -side[k] * bend * u; // > 0 towards the outer bank
          const shape = Math.max(0, (1 - u * u) * (1 + 0.8 * outer));
          const bed = zk - D * (0.15 + shape);
          z[k] = Math.min(z[k], bed);
          const low = zk - 0.55 * D; // low-flow stage: bars stand out of the water
          if (low > z[k]) { if (low > level[k]) level[k] = low; cover[k] = COVER.river; } else cover[k] = half[k] > 8 ? COVER.sand : COVER.gravel;
        } else {
          // banks, natural levees, and the floodplain the river has levelled
          const F = reach[k];
          const levee = 0.3 * D * Math.exp(-(u - 1) * 2.5) * (half[k] > 0.5 * this.o.meanderMinWidth ? 1 : 0.4);
          const target = zk + levee;
          const wgt = 1 - smooth01(0, F, u - 1);
          z[k] = z[k] > target ? z[k] + (target - z[k]) * wgt : z[k] + (target - z[k]) * wgt * 0.5;
          // fresh sand on the inside of bends (the top of the point bar); the outer bank is
          // a cut bank, vegetated to its edge
          const bendB = Math.max(-1, Math.min(1, curv[k] * half[k] * 10));
          const inner = side[k] * bendB > 0.15;
          if (half[k] > 0.5 * this.o.meanderMinWidth && u < 1 + 0.8 * Math.abs(bendB) && inner && (cover[k] === 255 || cover[k] === FLOODPLAIN)) cover[k] = COVER.sand;
        }
      } else if (u < 1) {
        // oxbow: an abandoned bend, deep (lake), filling (marsh) or nearly gone (scar)
        const age = kind[k] - 1;
        const bed = zk - D * (1 - u * u) - 0.2;
        z[k] = Math.min(z[k], bed);
        if (age <= 1) {
          const lv = zk - 0.1;
          if (lv > z[k]) { level[k] = Math.max(level[k], lv); cover[k] = COVER.lake; } else cover[k] = COVER.marsh;
        } else cover[k] = COVER.marsh;
      }
    }
    for (let k = 0; k < N; k++) if (best[k] < Infinity) near[k] = Math.max(near[k], 1 - Math.min(1, best[k] / (1 + reach[k])));
    return near;
  }

  /** Level the meander belts into floodplains (meadow, wet woodland, marsh). */
  _floodplains({ z, cover, near, W, H, X0, Y0, pixel }) {
    const rivers = this.rivers.filter((r) => r.meandering);
    if (!rivers.length) return;
    const F = 2 * this.o.belt; // in half-widths
    let wmax = 0;
    for (const r of rivers) wmax = Math.max(wmax, r.wmax);
    const M = Math.min(1500, Math.ceil((0.5 * wmax * (1 + F)) / pixel) + 1);
    const We = W + 2 * M, He = H + 2 * M;
    const site = new Int32Array(We * He).fill(-1);
    const sZw = [], sHalf = [], sDepth = [];
    const ex0 = X0 - M * pixel, ey0 = Y0 - M * pixel;
    for (const r of rivers) {
      for (let s2 = 1; s2 < r.nodes.length; s2++) {
        const a = r.nodes[s2 - 1], b = r.nodes[s2];
        const L = Math.hypot(b.x - a.x, b.y - a.y);
        const steps = Math.max(1, Math.ceil(L / (0.5 * pixel)));
        for (let q = 0; q <= steps; q++) {
          const t = q / steps;
          const i = Math.floor((a.x + t * (b.x - a.x) - ex0) / pixel), j = Math.floor((a.y + t * (b.y - a.y) - ey0) / pixel);
          if (i < 0 || j < 0 || i >= We || j >= He) continue;
          const k = j * We + i;
          if (site[k] >= 0) continue;
          site[k] = sZw.length;
          sZw.push(a.zw + t * (b.zw - a.zw));
          sHalf.push(0.5 * (a.w + t * (b.w - a.w)));
          sDepth.push(a.d + t * (b.d - a.d));
        }
      }
    }
    if (!sZw.length) return;
    const { d2, nearest } = featureTransform(site, We, He);
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const ke = (j + M) * We + (i + M);
        const id = nearest[ke];
        if (id < 0) continue;
        const u = (Math.sqrt(d2[ke]) * pixel) / sHalf[id];
        if (u >= 1 + F) continue;
        const k = j * W + i;
        // levelled, not planed flat: some of the old relief survives, and the plain rises
        // gently away from the river
        // and bounded by the valley walls: ground well above the river is left alone
        const reachUp = sDepth[id] + 0.6 * sHalf[id] + 1;
        // (but right beside the channel the ground is always levelled: where a bend swings
        // into high ground it cuts a bluff, not a slot)
        const far = smooth01(2, 4, u);
        const walls = 1 - far * smooth01(0.5 * reachUp, 1.5 * reachUp, z[k] - sZw[id]);
        const wgt = (1 - 0.15 * far) * (1 - smooth01(0.3 * F, F, u - 1)) * walls;
        const target = sZw[id] + 0.1 + 0.002 * u * sHalf[id];
        if (z[k] > target) z[k] += (target - z[k]) * wgt;
        if (wgt > 0.45 && cover[k] === 255) cover[k] = FLOODPLAIN;
        near[k] = Math.max(near[k], 0.6 * wgt);
      }
    }
  }
}

/**
 * Squared Euclidean distance (pixels^2) from every pixel to the nearest pixel whose mask
 * value equals `target`: the exact separable transform of Felzenszwalb & Huttenlocher
 * (2012, Theory of Computing 8).
 */
function edt(mask, W, H, target) {
  const INF = 1e20;
  const f = new Float64Array(Math.max(W, H));
  const d = new Float64Array(W * H);
  for (let k = 0; k < W * H; k++) d[k] = mask[k] === target ? 0 : INF;
  const v = new Int32Array(Math.max(W, H));
  const zz = new Float64Array(Math.max(W, H) + 1);
  const out = new Float64Array(Math.max(W, H));
  const pass = (n, get, set) => {
    for (let q = 0; q < n; q++) f[q] = get(q);
    let k = 0;
    v[0] = 0; zz[0] = -INF; zz[1] = INF;
    for (let q = 1; q < n; q++) {
      let s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
      while (s <= zz[k]) { k--; s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]); }
      k++; v[k] = q; zz[k] = s; zz[k + 1] = INF;
    }
    k = 0;
    for (let q = 0; q < n; q++) { while (zz[k + 1] < q) k++; out[q] = (q - v[k]) ** 2 + f[v[k]]; }
    for (let q = 0; q < n; q++) set(q, out[q]);
  };
  for (let x = 0; x < W; x++) pass(H, (q) => d[q * W + x], (q, val) => { d[q * W + x] = val; });
  for (let y = 0; y < H; y++) pass(W, (q) => d[y * W + q], (q, val) => { d[y * W + q] = val; });
  return d;
}

/**
 * Exact Euclidean feature transform: for every pixel, the squared distance (pixels^2) to
 * the nearest site (site >= 0) and that site's id. Felzenszwalb & Huttenlocher's
 * separable algorithm, carrying the site along.
 */
function featureTransform(site, W, H) {
  const INF = 1e20;
  const n = W * H;
  const d2 = new Float64Array(n);
  const nearest = new Int32Array(n);
  for (let k = 0; k < n; k++) { d2[k] = site[k] >= 0 ? 0 : INF; nearest[k] = site[k]; }
  const len = Math.max(W, H);
  const f = new Float64Array(len), id = new Int32Array(len), v = new Int32Array(len), zz = new Float64Array(len + 1);
  const outD = new Float64Array(len), outI = new Int32Array(len);
  const pass = (count, base, stride) => {
    for (let q = 0; q < count; q++) { f[q] = d2[base + q * stride]; id[q] = nearest[base + q * stride]; }
    let k = 0;
    v[0] = 0; zz[0] = -INF; zz[1] = INF;
    for (let q = 1; q < count; q++) {
      let s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
      while (s <= zz[k]) { k--; s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]); }
      k++; v[k] = q; zz[k] = s; zz[k + 1] = INF;
    }
    k = 0;
    for (let q = 0; q < count; q++) { while (zz[k + 1] < q) k++; outD[q] = (q - v[k]) ** 2 + f[v[k]]; outI[q] = id[v[k]]; }
    for (let q = 0; q < count; q++) { d2[base + q * stride] = outD[q]; nearest[base + q * stride] = outI[q]; }
  };
  for (let x = 0; x < W; x++) pass(H, x, W);
  for (let y = 0; y < H; y++) pass(W, y * W, 1);
  for (let k = 0; k < n; k++) if (d2[k] >= INF) nearest[k] = -1;
  return { d2, nearest };
}

const COVER_RGB = [
  [22, 52, 98], [52, 102, 150], [58, 108, 140], [52, 96, 128], [176, 160, 120], [222, 206, 158], [168, 160, 146], [104, 116, 82],
  [168, 172, 112], [128, 152, 84], [150, 140, 96], [58, 92, 50], [150, 142, 132], [104, 96, 90],
];

/**
 * Colour a rendered tile: land cover, lit by a hillshade of the detailed surface; water
 * by depth; forest canopies broken into crowns.
 * @returns {Uint8ClampedArray} RGBA
 */
export function detailColors(detail, tile, rgba) {
  const { width: W, height: H, pixel, z, level, cover, rock, x0, y0 } = tile;
  rgba = rgba || new Uint8ClampedArray(4 * W * H);
  const L = [-0.5, -0.5, Math.SQRT1_2];
  const seed = detail.seed;
  for (let j = 0; j < H; j++) {
    for (let i = 0; i < W; i++) {
      const k = j * W + i;
      const zl = z[j * W + Math.max(0, i - 1)], zr = z[j * W + Math.min(W - 1, i + 1)];
      const zu = z[Math.max(0, j - 1) * W + i], zd = z[Math.min(H - 1, j + 1) * W + i];
      const gx = (zr - zl) / (2 * pixel), gy = (zd - zu) / (2 * pixel);
      const shade = Math.max(0, (-gx * L[0] - gy * L[1] + L[2]) / Math.sqrt(gx * gx + gy * gy + 1));
      const c = cover[k];
      let col;
      if (level[k] > z[k]) {
        const dep = level[k] - z[k];
        const base = COVER_RGB[c] || COVER_RGB[COVER.shallowWater];
        const t = Math.min(1, dep / (c === COVER.river ? 4 : c === COVER.lake ? 2 : 12));
        // shallow water shows a sandy bed
        col = [0, 1, 2].map((m) => base[m] * t + COVER_RGB[COVER.sand][m] * (1 - t) * 0.6 + base[m] * (1 - t) * 0.4);
        col = col.map((v) => v * (0.9 + 0.15 * shade));
      } else {
        col = c === COVER.rock || c === COVER.cliff ? detail.rockColor[rock[k]].map((v) => v * (c === COVER.cliff ? 0.75 : 1)) : (COVER_RGB[c] || COVER_RGB[COVER.grass]).slice();
        if (c === COVER.forest) {
          // tree crowns: cellular (Worley-like) noise, crowns ~6 m across
          const wx = (x0 + (i + 0.5) * pixel) / 6, wy = (y0 + (j + 0.5) * pixel) / 6;
          const f = crown(wx, wy, seed);
          col = col.map((v) => v * (0.7 + 0.55 * f));
        } else if (c === COVER.grass || c === COVER.scrub || c === COVER.duneGrass) {
          const n = fbm(x0 + i * pixel, y0 + j * pixel, { wavelength: 18, octaves: 2, seed: seed + 3 });
          col = col.map((v) => v * (0.93 + 0.12 * n));
        }
        col = col.map((v) => v * (0.45 + 0.7 * shade));
      }
      rgba[4 * k] = col[0]; rgba[4 * k + 1] = col[1]; rgba[4 * k + 2] = col[2]; rgba[4 * k + 3] = 255;
    }
  }
  return rgba;
}

/** Brightness of a tree crown at (x, y) in crown units: 1 at a crown's centre, 0 between. */
function crown(x, y, seed) {
  const ix = Math.floor(x), iy = Math.floor(y);
  let d1 = 9;
  for (let a = -1; a <= 1; a++) {
    for (let b = -1; b <= 1; b++) {
      const cx = ix + a, cy = iy + b;
      let h = Math.imul(cx, 0x27d4eb2d) ^ Math.imul(cy, 0x165667b1) ^ Math.imul(seed, 0x9e3779b1);
      h = Math.imul(h ^ (h >>> 15), 0x85ebca6b);
      const px = cx + ((h >>> 8) & 255) / 255, py = cy + ((h >>> 16) & 255) / 255;
      const d = (px - x) ** 2 + (py - y) ** 2;
      if (d < d1) d1 = d;
    }
  }
  return Math.max(0, 1 - Math.sqrt(d1) * 1.4);
}
