// The simulation: state, the order of processes in a time step, and bookkeeping.
//
// STATE (per cell, all in metres):
//   rock       elevation of the bedrock surface
//   soil       thickness of loose material on top (regolith, alluvium, talus)
//   z          surface elevation = rock + soil (kept in sync)
//   upliftCum  total rock uplift so far (to look up strata in the rock's own frame)
//
// ONE STEP of length dt:
//   1. uplift       rock rises by U(x, y) dt (except at fixed base-level outlets)
//   2. weathering   bedrock converts to soil at P0 exp(-soil / h*)   (Heimsath et al. 1997, Nature 388)
//   3. vegetation   cover from climate and soil depth
//   4. flow         lakes, flow directions, drainage area
//   5. rivers       stream power incision + sediment deposition
//   6. creep        nonlinear soil creep
//   7. rockfall     Culmann cliff failure -> debris
//   8. avalanche    over-steep soil and talus slide to the angle of repose
//   9. boundaries   material reaching an open edge leaves the domain
//
// Operator splitting like this is first-order accurate in dt; keep dt small
// compared with the timescales you care about (e.g. the rockfall timescale).

import { HexGrid } from './hexgrid.js';
import { Rng, fbm } from './rng.js';
import { resolveParams, deepMerge, validate } from './params.js';
import { Lithology } from './lithology.js';
import { routeFlow } from './flow.js';
import { fluvial } from './processes/fluvial.js';
import { creep, avalanche } from './processes/hillslope.js';
import { rockfall } from './processes/rockfall.js';

export class Simulation {
  /**
   * @param {object|string} [params] parameter overrides (see params.js), optionally
   *   with `preset: name`; or just a preset name.
   */
  constructor(params = {}) {
    this.params = resolveParams(params);
    const g = this.params.grid;
    this.grid = new HexGrid(g.cols, g.rows, g.dx);
    const n = this.grid.size;
    this.rock = new Float64Array(n);
    this.soil = new Float64Array(n);
    this.z = new Float64Array(n);
    this.upliftCum = new Float64Array(n);
    this.veg = new Float64Array(n);
    this.surfaceRock = new Int32Array(n);
    this.outlet = new Uint8Array(n);
    this.boundaryOutlet = new Uint8Array(n);
    this.upliftField = new Float64Array(n);
    this.precipWeight = new Float64Array(n);
    this.vegMax = new Float64Array(n);
    /** Net surface change rate over the last step (m/yr), uplift included. */
    this.dzdt = new Float64Array(n);
    /** Recent rockfall thickness removed per cell (m), decaying over a few steps; for display. */
    this.rockfallActivity = new Float64Array(n);
    this.flow = {
      filled: new Float64Array(n),
      receiver: new Int32Array(n),
      order: new Int32Array(n),
      area: new Float64Array(n),
      visited: new Uint8Array(n),
      flowExponent: 0,
    };
    this.qs = new Float64Array(n);
    this.scratch = {
      eroded: new Float64Array(n),
      zNew: new Float64Array(n),
      qs: new Float64Array(n),
      edgeFlux: new Float64Array(3 * n),
      outTotal: new Float64Array(n),
      zStart: new Float64Array(n),
    };
    this.time = 0;
    this.stepCount = 0;
    this.rng = new Rng(this.params.seed);
    this.budget = { uplift: 0, exported: 0 };
    this.last = { substeps: 0, rockfallEvents: 0, rockfallVolume: 0, avalancheVolume: 0, exported: 0 };

    this._buildBoundaries();
    this._buildDerived();
    this._initialTerrain();
    this._updateDerivedState();
    routeFlow(this.grid, this.z, this.outlet, this.precipWeight, this.flow);
  }

  /** Rebuild a simulation exactly from snapshot() output. */
  static fromSnapshot(snap) {
    const sim = new Simulation(snap.params);
    sim.rock.set(snap.rock);
    sim.soil.set(snap.soil);
    sim.upliftCum.set(snap.upliftCum);
    sim.time = snap.time;
    sim.stepCount = snap.stepCount;
    sim.budget = { ...snap.budget };
    sim._updateDerivedState();
    routeFlow(sim.grid, sim.z, sim.outlet, sim.precipWeight, sim.flow);
    sim.rng.state = snap.rngState >>> 0; // after routing, which may draw random numbers
    return sim;
  }

  /** Everything needed to reproduce and continue this run, as plain JSON-able data. */
  snapshot() {
    return {
      format: 'edgy-erosion-snapshot/1',
      params: this.params,
      time: this.time,
      stepCount: this.stepCount,
      rngState: this.rng.state,
      budget: { ...this.budget },
      rock: Array.from(this.rock),
      soil: Array.from(this.soil),
      upliftCum: Array.from(this.upliftCum),
    };
  }

  /**
   * Change parameters mid-run (climate, rates, rocks, strata...). Grid geometry and
   * the initial condition cannot change; make a new Simulation for that.
   */
  setParams(overrides) {
    const next = deepMerge(this.params, overrides);
    for (const k of ['cols', 'rows', 'dx']) {
      if (next.grid[k] !== this.params.grid[k]) throw new Error('grid cannot change mid-run; create a new Simulation');
    }
    validate(next);
    this.params = next;
    this._buildBoundaries();
    this._buildDerived();
    this._updateDerivedState();
  }

  /** Advance one step. */
  step(dt = this.params.time.dt) {
    const p = this.params;
    const { grid, rock, soil, z, boundaryOutlet } = this;
    const n = grid.size;
    const ac = grid.cellArea;
    this.scratch.zStart.set(z);

    // 1. uplift
    let upliftVol = 0;
    for (let i = 0; i < n; i++) {
      if (boundaryOutlet[i]) continue;
      const u = this.upliftField[i] * dt;
      rock[i] += u;
      this.upliftCum[i] += u;
      upliftVol += u * ac;
    }
    this.budget.uplift += upliftVol;
    this._updateSurfaceRock();

    // 2. weathering
    if (p.weathering.enabled) {
      const wp = p.weathering;
      // climate: frost cracking speeds production, and so does moisture (chemical weathering)
      const frost = 1 + wp.frostBoost * p.climate.frost;
      const wet = Math.pow(Math.max(p.climate.precipitation, 1e-6), wp.moistureExponent ?? 0);
      const vb = p.vegetation.productionBoost;
      const hs = wp.productionDepth;
      for (let i = 0; i < n; i++) {
        if (boundaryOutlet[i] || soil[i] > 40 * hs) continue; // production negligible under deep soil
        // dh/dt = P exp(-h/h*) integrates exactly to h(t+dt) = h* ln(exp(h/h*) + P dt / h*)
        const P = this.lith.weathering[this.surfaceRock[i]] * frost * wet * (1 + vb * this.veg[i]);
        const h1 = hs * Math.log(Math.exp(soil[i] / hs) + (P * dt) / hs);
        const d = h1 - soil[i];
        if (!(d > 0)) continue;
        rock[i] -= d;
        soil[i] = h1;
      }
      this._updateSurfaceRock();
    }
    for (let i = 0; i < n; i++) z[i] = rock[i] + soil[i];

    // 3. vegetation, 4. flow
    this._updateVegetation();
    this._updateOutlets();
    routeFlow(grid, z, this.outlet, this.precipWeight, this.flow);

    // 5. rivers
    let exported = 0;
    if (p.fluvial.enabled) exported += fluvial(this, dt);

    // 6. creep
    this.last.substeps = p.hillslope.enabled ? creep(this, dt).substeps : 0;

    // 7. rockfall
    const act = this.rockfallActivity;
    for (let i = 0; i < n; i++) act[i] *= 0.7;
    if (p.rockfall.enabled) {
      const rf = rockfall(this, dt);
      this.last.rockfallEvents = rf.events;
      this.last.rockfallVolume = rf.volume;
    } else {
      this.last.rockfallEvents = 0;
      this.last.rockfallVolume = 0;
    }

    // 8. avalanche
    this.last.avalancheVolume = p.hillslope.enabled ? avalanche(this) : 0;

    // 9. open boundaries: loose material leaves; base level stays fixed
    for (let i = 0; i < n; i++) {
      if (boundaryOutlet[i] && soil[i] > 0) {
        exported += soil[i] * ac;
        soil[i] = 0;
      }
      z[i] = rock[i] + soil[i];
    }
    this.budget.exported += exported;
    this.last.exported = exported;

    this.time += dt;
    this.stepCount++;
    const zs = this.scratch.zStart;
    for (let i = 0; i < n; i++) this.dzdt[i] = (z[i] - zs[i]) / dt;
    this._updateDerivedState();
    // Re-route on the final surface so lakes, drainage area etc. match `z` for display
    // and analysis. (The next step routes again after uplift.)
    routeFlow(grid, z, this.outlet, this.precipWeight, this.flow);
    return this;
  }

  /** Advance `steps` steps; optional callback(sim) after each. */
  run(steps, onStep) {
    for (let s = 0; s < steps; s++) {
      this.step();
      if (onStep) onStep(this);
    }
    return this;
  }

  /** Total volume of rock + soil above z = 0 datum (m^3), for budget checks. */
  volume() {
    let v = 0;
    for (let i = 0; i < this.grid.size; i++) v += this.rock[i] + this.soil[i];
    return v * this.grid.cellArea;
  }

  /** Summary numbers for logging. */
  stats() {
    const { z, soil, grid } = this;
    let zmin = Infinity, zmax = -Infinity, zsum = 0, ssum = 0, bare = 0, lake = 0;
    const f = this.flow.filled;
    for (let i = 0; i < grid.size; i++) {
      if (z[i] < zmin) zmin = z[i];
      if (z[i] > zmax) zmax = z[i];
      zsum += z[i];
      ssum += soil[i];
      if (soil[i] < 0.05) bare++;
      if (f[i] - z[i] > 0.1) lake++;
    }
    return {
      time: this.time,
      step: this.stepCount,
      minElevation: zmin,
      maxElevation: zmax,
      meanElevation: zsum / grid.size,
      meanSoil: ssum / grid.size,
      bareRockFraction: bare / grid.size,
      lakeFraction: lake / grid.size,
      rockfallEvents: this.last.rockfallEvents,
      rockfallVolume: this.last.rockfallVolume,
      exportedVolume: this.last.exported,
      creepSubsteps: this.last.substeps,
    };
  }

  // ---------------------------------------------------------------- internals

  _buildBoundaries() {
    const open = new Set(this.params.boundary.open);
    const g = this.grid;
    for (let i = 0; i < g.size; i++) {
      const e = g.edge[i];
      this.boundaryOutlet[i] =
        ((e & 1) && open.has('west')) || ((e & 2) && open.has('east')) ||
        ((e & 4) && open.has('north')) || ((e & 8) && open.has('south')) ? 1 : 0;
    }
  }

  _buildDerived() {
    const p = this.params;
    const g = this.grid;
    this.flow.flowExponent = p.fluvial.routing === 'mfd' ? p.fluvial.flowExponent : 0;
    this.flow.randomExponent = p.fluvial.randomReceivers ? p.fluvial.randomExponent : 0;
    this.flow.rng = this.rng;
    this.lith = new Lithology(p, g);
    const W = g.width, H = g.height;
    const u = p.uplift;
    const ang = ((u.angle || 0) * Math.PI) / 180;
    const ux = Math.cos(ang), uy = -Math.sin(ang);
    for (let i = 0; i < g.size; i++) {
      const fx = g.x[i] / W - 0.5;
      const fy = (g.y[i] / W) - 0.5 * (H / W); // width units, centred
      let f = 1;
      switch (u.pattern) {
        case 'uniform': f = 1; break;
        case 'dome': f = Math.exp(-(fx * fx + fy * fy) / (2 * u.radius * u.radius)); break;
        case 'tilt': f = clamp01(0.5 + fx * ux + fy * uy); break;
        case 'fault': f = -fx * uy + fy * ux + (u.offset || 0) < 0 ? 1 : 0; break;
        default: throw new Error(`unknown uplift pattern "${u.pattern}"`);
      }
      this.upliftField[i] = this.boundaryOutlet[i] ? 0 : u.rate * f;

      const P = Math.max(0, p.climate.precipitation * (1 + p.climate.precipitationGradient * (fx)));
      this.precipWeight[i] = P; // relative to a 1 m/yr reference
      const vc = p.vegetation.cover;
      this.vegMax[i] = vc === 'auto' ? P / (P + p.vegetation.halfSaturation) : Math.max(0, Math.min(1, +vc));
    }
  }

  _initialTerrain() {
    const p = this.params.initial;
    const g = this.grid;
    const W = g.width, H = g.height;
    const ang = ((p.angle || 0) * Math.PI) / 180;
    const ux = Math.cos(ang), uy = -Math.sin(ang);
    const noise = p.noise || { amplitude: 0 };
    for (let i = 0; i < g.size; i++) {
      const fx = g.x[i] / W;
      const fy = g.y[i] / H;
      let shape = 0;
      switch (p.type) {
        case 'flat': shape = 0; break;
        case 'plateau': {
          // distance inside the plateau rectangle, in width units
          let d = Math.min(fx - p.margin, 1 - p.margin - fx, (fy - p.margin) * (H / W), (1 - p.margin - fy) * (H / W));
          if (p.edgeNoise) d += p.edgeNoise * fbm(g.x[i], g.y[i], { wavelength: 0.3 * W, octaves: 4, seed: (this.params.seed | 0) + 31 });
          shape = smoothstep(clamp01(d / p.edgeWidth + 0.5));
          break;
        }
        case 'dome': {
          const rx = fx - 0.5, ry = (fy - 0.5) * (H / W);
          shape = Math.exp(-(rx * rx + ry * ry) / (2 * 0.22 * 0.22));
          break;
        }
        case 'ramp': shape = clamp01(0.5 + (fx - 0.5) * ux + (fy - 0.5) * (H / W) * uy); break;
        case 'ridge': shape = clamp01(1 - 2 * Math.abs((fx - 0.5) * -uy + (fy - 0.5) * (H / W) * ux)); break;
        default: throw new Error(`unknown initial terrain type "${p.type}"`);
      }
      let zi = p.base + p.height * shape;
      if (noise.amplitude) {
        zi += noise.amplitude * fbm(g.x[i], g.y[i], { wavelength: noise.wavelength, octaves: noise.octaves ?? 5, gain: noise.gain ?? 0.5, seed: this.params.seed | 0 });
      }
      this.rock[i] = zi - (this.boundaryOutlet[i] ? 0 : p.soil);
      this.soil[i] = this.boundaryOutlet[i] ? 0 : p.soil;
    }
    // Open edges are base level: no higher than `base`, so rivers have somewhere to go.
    for (let i = 0; i < g.size; i++) if (this.boundaryOutlet[i]) this.rock[i] = Math.min(this.rock[i], p.base);
  }

  _updateSurfaceRock() {
    const { rock, upliftCum, lith, surfaceRock } = this;
    // sample just below the surface so a cell sitting exactly on a boundary uses the rock beneath
    for (let i = 0; i < this.grid.size; i++) surfaceRock[i] = lith.rockAt(i, rock[i] - upliftCum[i] - 1e-6);
  }

  _updateVegetation() {
    const hv = this.params.vegetation.soilDepthScale;
    for (let i = 0; i < this.grid.size; i++) this.veg[i] = this.vegMax[i] * (1 - Math.exp(-this.soil[i] / hv));
  }

  _updateOutlets() {
    const sea = this.params.boundary.seaLevel;
    for (let i = 0; i < this.grid.size; i++) {
      this.outlet[i] = this.boundaryOutlet[i] || (sea != null && this.z[i] < sea) ? 1 : 0;
    }
  }

  _updateDerivedState() {
    for (let i = 0; i < this.grid.size; i++) this.z[i] = this.rock[i] + this.soil[i];
    this._updateSurfaceRock();
    this._updateVegetation();
    this._updateOutlets();
  }
}

function clamp01(v) {
  return v < 0 ? 0 : v > 1 ? 1 : v;
}

function smoothstep(t) {
  return t * t * (3 - 2 * t);
}
