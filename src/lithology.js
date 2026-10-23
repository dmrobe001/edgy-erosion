// Lithology: which rock lies at a given place and depth.
//
// Rock type is a function of position (x, y) and of elevation in the rock's own
// reference frame. Uplift carries the strata up with it, so we look rock up at
//   z_frame = z - cumulativeUplift[i]
// i.e. the elevation the rock had at t = 0. The stratigraphic column is a stack of
// layers whose top surface is a plane (optionally dipping, optionally undulating):
//
//   columnTop(x, y) = top - tan(dip) * s(x, y) + undulation(x, y)
//
// where s is the distance along the dip direction from the domain centre. Vertical
// bodies (plugs, dikes) override the strata below their `top` elevation.

import { fbm } from './rng.js';

export class Lithology {
  /**
   * @param {object} params resolved parameters
   * @param {import('./hexgrid.js').HexGrid} grid
   */
  constructor(params, grid) {
    this.grid = grid;
    const names = Object.keys(params.rocks);
    this.names = names;
    const R = names.length;
    this.erodibility = new Float64Array(R);
    this.cohesion = new Float64Array(R);
    this.frictionAngle = new Float64Array(R); // radians
    this.weathering = new Float64Array(R);
    this.color = [];
    names.forEach((name, k) => {
      const r = params.rocks[name];
      this.erodibility[k] = r.erodibility;
      this.cohesion[k] = r.cohesion;
      this.frictionAngle[k] = (r.frictionAngle * Math.PI) / 180;
      this.weathering[k] = r.weathering;
      this.color.push(r.color || '#888888');
    });
    const idx = (name) => names.indexOf(name);

    const st = params.strata;
    this.layerRock = Int32Array.from(st.layers.map((l) => idx(l.rock)));
    // layerBottom[k] = depth below column top of the base of layer k
    this.layerBottom = new Float64Array(st.layers.length);
    let d = 0;
    st.layers.forEach((l, k) => { d += l.thickness; this.layerBottom[k] = d; });
    this.columnThickness = d;
    this.repeat = !!st.repeat;
    this.basement = idx(st.basement);

    const n = grid.size;
    const W = grid.width;
    const H = grid.height;
    const xc = W / 2;
    const yc = H / 2;
    const tanDip = Math.tan(((st.dip || 0) * Math.PI) / 180);
    const az = ((st.dipDirection || 0) * Math.PI) / 180;
    // screen y points south, so "north" is -y
    const ux = Math.cos(az);
    const uy = -Math.sin(az);
    const und = st.undulation || { amplitude: 0 };
    this.columnTop = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const s = (grid.x[i] - xc) * ux + (grid.y[i] - yc) * uy;
      let top = st.top - tanDip * s;
      if (und.amplitude) top += und.amplitude * fbm(grid.x[i], grid.y[i], { wavelength: und.wavelength, octaves: 3, seed: (params.seed | 0) + 77 });
      this.columnTop[i] = top;
    }

    // Spatial variation in rock strength (fracture density, cementation...): a smooth
    // random field f(x, y) with E[ln f] = 0. Erodibility is multiplied by f, cohesion
    // divided by it. Breaks up the perfectly regular patterns uniform rock would give.
    const sv = params.strata.strengthVariation || { amplitude: 0 };
    this.weakness = new Float64Array(n).fill(1);
    if (sv.amplitude) {
      for (let i = 0; i < n; i++) {
        this.weakness[i] = Math.exp(sv.amplitude * fbm(grid.x[i], grid.y[i], { wavelength: sv.wavelength, octaves: 4, seed: (params.seed | 0) + 191 }));
      }
    }

    // Vertical bodies
    this.bodyRock = new Int32Array(n).fill(-1);
    this.bodyTop = new Float64Array(n).fill(Infinity);
    for (const b of params.bodies || []) {
      const rock = idx(b.rock);
      const top = b.top ?? Infinity;
      for (let i = 0; i < n; i++) {
        const px = grid.x[i] / W;
        const py = grid.y[i] / H;
        let inside = false;
        if (b.shape === 'disc') {
          // radius is a fraction of width; measure in width units on both axes
          const ddx = px - b.x;
          const ddy = (py - b.y) * (H / W);
          inside = ddx * ddx + ddy * ddy <= b.radius * b.radius;
        } else if (b.shape === 'dike') {
          inside = segmentDistance(px, py * (H / W), b.x0, b.y0 * (H / W), b.x1, b.y1 * (H / W)) <= b.width / 2;
        } else {
          throw new Error(`unknown body shape "${b.shape}"`);
        }
        if (inside) {
          this.bodyRock[i] = rock;
          this.bodyTop[i] = top;
        }
      }
    }
  }

  /**
   * Rock index at cell i and frame elevation z (elevation corrected for uplift).
   */
  rockAt(i, z) {
    if (this.bodyRock[i] >= 0 && z <= this.bodyTop[i]) return this.bodyRock[i];
    let depth = this.columnTop[i] - z;
    if (depth < 0) return this.layerRock[0]; // above the column: top layer extends upward
    if (depth >= this.columnThickness) {
      if (!this.repeat) return this.basement;
      depth %= this.columnThickness;
    }
    const lb = this.layerBottom;
    for (let k = 0; k < lb.length; k++) if (depth < lb[k]) return this.layerRock[k];
    return this.layerRock[lb.length - 1];
  }

  /**
   * Thickness-weighted mean cohesion and tan(friction angle) of the rock between
   * frame elevations zLow and zHigh in cell i. This is the strength of a cliff
   * face cutting through several layers: a failure plane must shear through all of them.
   * Samples at `samples` evenly spaced elevations (cheap, and exact enough for display-scale layers).
   */
  meanStrength(i, zLow, zHigh, samples = 8) {
    if (zHigh <= zLow) {
      const k = this.rockAt(i, zHigh);
      return { cohesion: this.cohesion[k] / this.weakness[i], phi: this.frictionAngle[k] };
    }
    let c = 0;
    let t = 0;
    const h = (zHigh - zLow) / samples;
    for (let s = 0; s < samples; s++) {
      const k = this.rockAt(i, zLow + (s + 0.5) * h);
      c += this.cohesion[k];
      t += Math.tan(this.frictionAngle[k]);
    }
    return { cohesion: c / samples / this.weakness[i], phi: Math.atan(t / samples) };
  }
}

function segmentDistance(px, py, x0, y0, x1, y1) {
  const vx = x1 - x0;
  const vy = y1 - y0;
  const len2 = vx * vx + vy * vy;
  let t = len2 > 0 ? ((px - x0) * vx + (py - y0) * vy) / len2 : 0;
  t = Math.max(0, Math.min(1, t));
  const dx = px - (x0 + t * vx);
  const dy = py - (y0 + t * vy);
  return Math.sqrt(dx * dx + dy * dy);
}
