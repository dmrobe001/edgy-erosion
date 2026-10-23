// Turning simulation state into pixels. Pure functions on typed arrays, so they
// work the same in a browser (into ImageData) and in Node (into a PNG file).

export const VIEWS = {
  terrain: 'Terrain (rock, soil, vegetation, water)',
  elevation: 'Elevation (hypsometric tint)',
  geology: 'Geology (rock at the surface)',
  soil: 'Soil thickness',
  vegetation: 'Vegetation cover',
  discharge: 'Drainage area',
  erosion: 'Erosion / deposition rate',
  rockfall: 'Recent rockfall',
  waves: 'Wave exposure',
};

export function hexToRgb(hex) {
  const v = parseInt(hex.replace('#', ''), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

/** Piecewise-linear colormap; stops = [[t, [r,g,b]], ...] with t ascending in [0, 1]. */
function ramp(stops, t) {
  if (!(t > 0)) return stops[0][1];
  if (t >= 1) return stops[stops.length - 1][1];
  for (let k = 1; k < stops.length; k++) {
    if (t <= stops[k][0]) {
      const [t0, c0] = stops[k - 1];
      const [t1, c1] = stops[k];
      const u = (t - t0) / (t1 - t0);
      return [c0[0] + u * (c1[0] - c0[0]), c0[1] + u * (c1[1] - c0[1]), c0[2] + u * (c1[2] - c0[2])];
    }
  }
  return stops[stops.length - 1][1];
}

const HYPSO = [[0, [62, 110, 70]], [0.25, [120, 150, 85]], [0.5, [190, 175, 115]], [0.75, [150, 110, 80]], [0.92, [200, 190, 185]], [1, [250, 250, 250]]];
const SEQ = [[0, [250, 245, 225]], [0.33, [220, 180, 100]], [0.66, [160, 95, 50]], [1, [70, 40, 25]]];
const GREEN = [[0, [235, 225, 200]], [0.5, [150, 185, 95]], [1, [35, 95, 40]]];
const BLUE = [[0, [245, 245, 240]], [0.4, [150, 190, 220]], [1, [20, 60, 140]]];
const HEAT = [[0, [0, 0, 0]], [0.3, [200, 40, 20]], [0.7, [255, 170, 30]], [1, [255, 255, 200]]];
const SOIL_RGB = [176, 152, 112];
const SAND_RGB = [222, 204, 150];
const VEG_RGB = [72, 112, 52];
const LAKE_RGB = [70, 120, 175];
const SEA_RGB = [40, 80, 140];

/**
 * Lambertian hillshade per cell, light from the north-west at 45 degrees.
 * Gradient by least squares over the available neighbours.
 * @returns {Float32Array} values in [0, 1]
 */
export function hillshade(sim, exaggeration = 1) {
  const { grid, z } = sim;
  const n = grid.size;
  const nb = grid.neighbors;
  const D = grid.degree;
  const out = new Float32Array(n);
  const L = [-0.5, -0.5, Math.SQRT1_2]; // unit vector towards the light (x east, y south, z up)
  const off = [0, 0];
  for (let i = 0; i < n; i++) {
    let sxx = 0, sxy = 0, syy = 0, sxz = 0, syz = 0;
    for (let k = 0; k < D; k++) {
      const j = nb[D * i + k];
      if (j < 0) continue;
      // meshes with curved geometry (the sphere) give displacements in the local tangent plane
      if (grid.offset) grid.offset(i, j, off);
      const ddx = grid.offset ? off[0] : grid.x[j] - grid.x[i];
      const ddy = grid.offset ? off[1] : grid.y[j] - grid.y[i];
      const dz = (z[j] - z[i]) * exaggeration;
      sxx += ddx * ddx; sxy += ddx * ddy; syy += ddy * ddy;
      sxz += ddx * dz; syz += ddy * dz;
    }
    const det = sxx * syy - sxy * sxy;
    let gx = 0, gy = 0;
    if (det > 1e-12) {
      gx = (syy * sxz - sxy * syz) / det;
      gy = (sxx * syz - sxy * sxz) / det;
    }
    const norm = Math.sqrt(gx * gx + gy * gy + 1);
    out[i] = Math.max(0, (-gx * L[0] - gy * L[1] + L[2]) / norm);
  }
  return out;
}

/**
 * Colour every cell for the given view.
 * @param {import('./simulation.js').Simulation} sim
 * @param {string} view one of VIEWS
 * @param {Uint8ClampedArray} [out] 3 bytes per cell
 */
export function cellColors(sim, view = 'terrain', out) {
  const { grid, z, soil, lith, surfaceRock, veg } = sim;
  const n = grid.size;
  out = out || new Uint8ClampedArray(3 * n);
  const shade = hillshade(sim, view === 'terrain' || view === 'geology' ? 2 : 1.5);
  const rockRgb = lith.color.map(hexToRgb);
  const { filled, area } = sim.flow;
  const sea = sim.seaLevel;
  const cellA = grid.cellArea;

  let zmin = Infinity, zmax = -Infinity, amax = 0, emax = 1e-9;
  for (let i = 0; i < n; i++) {
    if (z[i] < zmin) zmin = z[i];
    if (z[i] > zmax) zmax = z[i];
    if (area[i] > amax) amax = area[i];
    const e = Math.abs(sim.dzdt[i]);
    if (e > emax) emax = e;
  }
  const zr = Math.max(1e-6, zmax - zmin);
  const logAmax = Math.log(Math.max(amax / cellA, 2));
  const cp = sim.params.coast;
  const beachTop = sea != null && cp ? sea + cp.tidalRange / 2 + cp.waveHeight : null;
  const riverStart = Math.min(Math.log(150), 0.6 * logAmax); // contributing cells at which rivers show

  for (let i = 0; i < n; i++) {
    let c;
    const s = shade[i];
    const lit = 0.35 + 0.75 * s;
    const water = sea != null && z[i] < sea ? 2 : filled[i] - z[i] > 0.05 ? 1 : 0;
    switch (view) {
      case 'elevation':
        c = ramp(HYPSO, (z[i] - zmin) / zr).map((v) => v * lit);
        break;
      case 'geology':
        c = rockRgb[surfaceRock[i]].map((v) => v * lit);
        break;
      case 'soil': {
        const t = (Math.log10(Math.max(soil[i], 0.01)) + 2) / 3; // 1 cm .. 10 m
        c = ramp(SEQ, t).map((v) => v * (0.6 + 0.4 * s));
        break;
      }
      case 'vegetation':
        c = ramp(GREEN, veg[i]).map((v) => v * (0.6 + 0.4 * s));
        break;
      case 'discharge': {
        const t = Math.log(Math.max(area[i] / cellA, 1)) / logAmax;
        c = ramp(BLUE, t).map((v) => v * (0.75 + 0.25 * s));
        break;
      }
      case 'erosion': {
        // red = lowering, blue = rising (uplift included); square-root stretch
        const r = sim.dzdt[i] / emax;
        const t = Math.sign(r) * Math.sqrt(Math.abs(r));
        c = t < 0 ? [255, 255 * (1 + t), 255 * (1 + t)] : [255 * (1 - t), 255 * (1 - t), 255];
        c = c.map((v) => v * (0.7 + 0.3 * s));
        break;
      }
      case 'rockfall': {
        const g = 90 + 120 * s;
        const t = Math.min(1, Math.sqrt(sim.rockfallActivity[i] / 10));
        c = t > 0.01 ? ramp(HEAT, 0.3 + 0.7 * t) : [g, g, g];
        break;
      }
      case 'waves': {
        const g = 80 + 110 * s;
        const e = sim.waveExposureField[i];
        c = water === 2 && e === 0 ? [g * 0.5, g * 0.55, g * 0.7] : e > 0 ? ramp(HEAT, 0.25 + 0.75 * e) : [g, g, g];
        if (sim.notch[i] > 0) c = [255, 255, 255];
        break;
      }
      case 'terrain':
      default: {
        const rc = rockRgb[surfaceRock[i]];
        const cover = 1 - Math.exp(-soil[i] / 0.3);
        c = [0, 1, 2].map((k) => rc[k] + cover * (SOIL_RGB[k] - rc[k]));
        c = c.map((v, k) => v + veg[i] * (VEG_RGB[k] - v));
        if (beachTop !== null && soil[i] > 0.2 && z[i] < beachTop) c = SAND_RGB.slice(); // beach sand
        c = c.map((v) => v * lit);
        const t = (Math.log(Math.max(area[i] / cellA, 1)) - riverStart) / (logAmax - riverStart);
        if (t > 0 && !water) c = c.map((v, k) => v + Math.min(1, 0.4 + t) * (LAKE_RGB[k] - v));
        if (water === 1) c = LAKE_RGB.map((v) => v * (0.85 + 0.15 * s));
        if (water === 2) {
          // shallow water shows the bed (platforms, shoals); deep water darkens
          const depth = Math.min(1, (sea - z[i]) / 40);
          const bed = soil[i] > 0.2 ? SAND_RGB : rockRgb[surfaceRock[i]];
          const shallow = Math.exp(-(sea - z[i]) / 4);
          c = SEA_RGB.map((v, k) => (v * (1 - 0.5 * depth) + 40 * (1 - depth)) * (1 - 0.5 * shallow) + 0.5 * shallow * bed[k] * (0.7 + 0.3 * s));
        }
      }
    }
    out[3 * i] = c[0];
    out[3 * i + 1] = c[1];
    out[3 * i + 2] = c[2];
  }
  return out;
}

/**
 * For each pixel of a width x height image, the cell it shows (-1 = none).
 * The grid is fitted inside the image preserving aspect ratio.
 */
export function pixelMap(grid, width, height) {
  if (!grid.cellAt) throw new Error('drawing a map needs a mesh with cellAt()');
  const b = grid.bounds;
  const x0 = b.xmin;
  const y0 = b.ymin;
  const scale = Math.min(width / (b.xmax - b.xmin), height / (b.ymax - b.ymin)); // pixels per metre
  const map = new Int32Array(width * height);
  for (let py = 0; py < height; py++) {
    for (let px = 0; px < width; px++) {
      map[py * width + px] = grid.cellAt(x0 + (px + 0.5) / scale, y0 + (py + 0.5) / scale);
    }
  }
  return { map, width, height, scale, x0, y0 };
}

/** Fill an RGBA buffer from per-cell colours through a pixel map. */
export function paint(pm, colors, rgba, background = [24, 24, 28]) {
  const { map } = pm;
  for (let p = 0; p < map.length; p++) {
    const i = map[p];
    const o = 4 * p;
    if (i < 0) {
      rgba[o] = background[0]; rgba[o + 1] = background[1]; rgba[o + 2] = background[2];
    } else {
      rgba[o] = colors[3 * i]; rgba[o + 1] = colors[3 * i + 1]; rgba[o + 2] = colors[3 * i + 2];
    }
    rgba[o + 3] = 255;
  }
}

/**
 * Vertical cross-section along one grid row: sky, water, soil, and the rock
 * layers below the surface, coloured by lithology.
 * @returns {{zmin:number, zmax:number, exaggeration:number}}
 */
export function crossSection(sim, row, width, height, rgba) {
  const { grid, z, rock, lith, upliftCum } = sim;
  const filled = sim.flow.filled;
  const sea = sim.seaLevel;
  const cols = grid.cols;
  row = Math.max(0, Math.min(grid.rows - 1, row | 0));
  let zmin = Infinity, zmax = -Infinity;
  for (let c = 0; c < cols; c++) {
    const i = row * cols + c;
    zmin = Math.min(zmin, rock[i]);
    zmax = Math.max(zmax, z[i], sea ?? -Infinity);
  }
  const relief = Math.max(zmax - zmin, 10);
  zmin -= 0.35 * relief;
  zmax += 0.08 * relief;
  const rockRgb = lith.color.map(hexToRgb);
  const sky = [34, 36, 44];
  const water = [70, 120, 175];
  const soilC = [196, 172, 128];
  for (let px = 0; px < width; px++) {
    const c = Math.min(cols - 1, Math.floor((px / width) * cols));
    const i = row * cols + c;
    const wl = Math.max(filled[i], sea ?? -Infinity);
    for (let py = 0; py < height; py++) {
      const e = zmax - ((py + 0.5) / height) * (zmax - zmin);
      let col;
      if (e > z[i]) col = e <= wl ? water : sky;
      else if (e > rock[i]) col = soilC;
      else col = rockRgb[lith.rockAt(i, e - upliftCum[i])];
      const o = 4 * (py * width + px);
      rgba[o] = col[0]; rgba[o + 1] = col[1]; rgba[o + 2] = col[2]; rgba[o + 3] = 255;
    }
  }
  const horiz = grid.width / width; // metres per pixel horizontally
  const vert = (zmax - zmin) / height;
  return { zmin, zmax, exaggeration: horiz / vert };
}

// ------------------------------------------------------------------ planets

export const PLANET_VIEWS = {
  terrain: 'Terrain (land and sea floor)',
  plates: 'Plates and boundaries',
  age: 'Ocean-floor age',
  crust: 'Crust thickness',
  discharge: 'Drainage area',
  erosion: 'Erosion / deposition rate',
};

const LAND = [[0, [72, 120, 72]], [0.06, [110, 150, 85]], [0.18, [185, 175, 110]], [0.38, [160, 120, 80]], [0.62, [130, 105, 90]], [0.8, [225, 225, 225]], [1, [255, 255, 255]]];
const OCEAN = [[0, [150, 205, 225]], [0.06, [90, 160, 205]], [0.35, [45, 95, 165]], [0.7, [25, 55, 115]], [1, [10, 25, 60]]];
const AGE = [[0, [200, 40, 40]], [0.15, [235, 140, 40]], [0.3, [240, 220, 70]], [0.5, [90, 190, 110]], [0.75, [60, 120, 200]], [1, [70, 40, 140]]];
const CRUST = [[0, [20, 30, 70]], [0.12, [50, 90, 160]], [0.35, [120, 170, 120]], [0.5, [220, 200, 110]], [0.75, [200, 100, 50]], [1, [110, 20, 30]]];
const BOUNDARY_RGB = [null, [255, 70, 60], [255, 235, 90], [235, 235, 235]];

/** A distinct colour for plate q (golden-angle hues). */
export function plateColor(q) {
  const h = (q * 137.508) % 360;
  const s = 0.45, l = 0.55;
  const f = (k) => {
    const a = s * Math.min(l, 1 - l);
    const kk = (k + h / 30) % 12;
    return 255 * (l - a * Math.max(-1, Math.min(kk - 3, 9 - kk, 1)));
  };
  return [f(0), f(8), f(4)];
}

/**
 * Colour every cell of a planet for the given view (see PLANET_VIEWS).
 * @param {import('./planet.js').Planet} planet
 * @param {Uint8ClampedArray} [out] 3 bytes per cell
 */
export function planetColors(planet, view = 'terrain', out) {
  const sim = planet.sim;
  const mesh = planet.mesh;
  const n = mesh.size;
  out = out || new Uint8ClampedArray(3 * n);
  const z = sim.z;
  const sea = planet.seaLevel;
  const shade = hillshade(sim, 25);
  const area = sim.flow.area;
  let amax = 0, emax = 1e-12;
  for (let i = 0; i < n; i++) {
    if (area[i] > amax) amax = area[i];
    if (z[i] > sea) emax = Math.max(emax, Math.abs(sim.dzdt[i]));
  }
  const logAmax = Math.log(Math.max(amax / mesh.cellArea, 2));
  const riverStart = Math.min(Math.log(40), 0.55 * logAmax);
  for (let i = 0; i < n; i++) {
    const s = shade[i];
    const lit = 0.45 + 0.65 * s;
    const wet = z[i] < sea;
    let c;
    switch (view) {
      case 'plates': {
        c = plateColor(planet.owner[i]);
        const k = 0.55 + 0.45 * planet.continental[i];
        c = c.map((v) => v * k * (0.8 + 0.2 * s));
        const b = planet.boundary[i];
        if (b) c = BOUNDARY_RGB[b];
        break;
      }
      case 'age': {
        if (planet.continental[i] > 0.5) { const g = 70 + 90 * s; c = [g, g, g]; break; }
        c = ramp(AGE, planet.age[i] / 2e8).map((v) => v * (0.8 + 0.2 * s));
        break;
      }
      case 'crust':
        c = ramp(CRUST, planet.thickness[i] / 80000).map((v) => v * (0.75 + 0.25 * s));
        break;
      case 'discharge': {
        const t = Math.log(Math.max(area[i] / mesh.cellArea, 1)) / logAmax;
        c = wet ? [30, 40, 60] : ramp(BLUE, t).map((v) => v * (0.75 + 0.25 * s));
        break;
      }
      case 'erosion': {
        if (wet) { c = [40, 50, 70]; break; }
        const r = sim.dzdt[i] / emax;
        const t = Math.sign(r) * Math.sqrt(Math.abs(r));
        c = t < 0 ? [255, 255 * (1 + t), 255 * (1 + t)] : [255 * (1 - t), 255 * (1 - t), 255];
        c = c.map((v) => v * (0.7 + 0.3 * s));
        break;
      }
      case 'terrain':
      default: {
        if (wet) {
          c = ramp(OCEAN, (sea - z[i]) / 8000).map((v) => v * (0.85 + 0.25 * s));
        } else {
          c = ramp(LAND, (z[i] - sea) / 6000).map((v) => v * lit);
          const t = (Math.log(Math.max(area[i] / mesh.cellArea, 1)) - riverStart) / (logAmax - riverStart);
          if (t > 0) c = c.map((v, k) => v + Math.min(1, 0.3 + t) * (LAKE_RGB[k] - v));
        }
      }
    }
    out[3 * i] = c[0];
    out[3 * i + 1] = c[1];
    out[3 * i + 2] = c[2];
  }
  return out;
}

/**
 * Pixel map for a sphere mesh: 'map' (equirectangular, the whole planet, 2:1) or
 * 'globe' (orthographic view centred on lon0, lat0 in degrees). Pixels off the globe are -1.
 */
export function spherePixelMap(mesh, width, height, { projection = 'map', lon0 = 0, lat0 = 0 } = {}) {
  const map = new Int32Array(width * height);
  let hint = 0;
  if (projection === 'globe') {
    const l0 = (lon0 * Math.PI) / 180, p0 = (lat0 * Math.PI) / 180;
    // camera frame: c = centre, e = east, u = north
    const c = [Math.cos(p0) * Math.cos(l0), Math.cos(p0) * Math.sin(l0), Math.sin(p0)];
    const e = [-Math.sin(l0), Math.cos(l0), 0];
    const u = [c[1] * e[2] - c[2] * e[1], c[2] * e[0] - c[0] * e[2], c[0] * e[1] - c[1] * e[0]];
    const r = 0.5 * Math.min(width, height) - 1;
    for (let py = 0; py < height; py++) {
      for (let px = 0; px < width; px++) {
        const a = (px + 0.5 - width / 2) / r, b = -(py + 0.5 - height / 2) / r;
        const d = 1 - a * a - b * b;
        if (d < 0) { map[py * width + px] = -1; continue; }
        const w = Math.sqrt(d);
        hint = mesh.nearest(w * c[0] + a * e[0] + b * u[0], w * c[1] + a * e[1] + b * u[1], w * c[2] + a * e[2] + b * u[2], hint);
        map[py * width + px] = hint;
      }
    }
    return { map, width, height, projection, lon0, lat0, radiusPx: r };
  }
  for (let py = 0; py < height; py++) {
    const lat = Math.PI / 2 - ((py + 0.5) / height) * Math.PI;
    for (let px = 0; px < width; px++) {
      const lon = ((px + 0.5) / width) * 2 * Math.PI - Math.PI;
      hint = mesh.nearest(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat), hint);
      map[py * width + px] = hint;
    }
  }
  return { map, width, height, projection, scale: width / mesh.width, x0: 0, y0: 0 };
}

/**
 * Cross-section of a planet along a parallel (constant latitude, degrees): the surface,
 * sea, and the crust down to its base, coloured continental or oceanic, floating on
 * the mantle. Shows isostasy at work: thick crust has deep roots under high ground.
 */
export function planetSection(planet, latDeg, width, height, rgba) {
  const mesh = planet.mesh;
  const z = planet.sim.z;
  const sea = planet.seaLevel;
  const lat = (latDeg * Math.PI) / 180;
  const cells = new Int32Array(width);
  let hint = 0;
  for (let px = 0; px < width; px++) {
    const lon = ((px + 0.5) / width) * 2 * Math.PI - Math.PI;
    hint = mesh.nearest(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat), hint);
    cells[px] = hint;
  }
  const top = 10000, bottom = -90000;
  const sky = [34, 36, 44], water = [55, 100, 165], mantle = [120, 60, 50];
  const contC = [196, 160, 120], oceC = [70, 75, 85];
  for (let px = 0; px < width; px++) {
    const i = cells[px];
    const base = z[i] - planet.thickness[i];
    const cc = planet.continental[i];
    const crustC = [0, 1, 2].map((k) => oceC[k] + cc * (contC[k] - oceC[k]));
    for (let py = 0; py < height; py++) {
      const e = top - ((py + 0.5) / height) * (top - bottom);
      const col = e > z[i] ? (e <= sea ? water : sky) : e > base ? crustC : mantle;
      const o = 4 * (py * width + px);
      rgba[o] = col[0]; rgba[o + 1] = col[1]; rgba[o + 2] = col[2]; rgba[o + 3] = 255;
    }
  }
  return { zmin: bottom, zmax: top, exaggeration: (2 * Math.PI * planet.params.radius * Math.cos(lat) / width) / ((top - bottom) / height) };
}
