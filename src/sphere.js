// A whole planet: the icosahedral geodesic mesh.
//
// Start from an icosahedron (12 vertices, 20 triangles), split every triangle into
// f x f smaller ones, and push the new vertices out onto the sphere. The vertices are
// the cell centres; the cells are their spherical Voronoi polygons: 12 pentagons
// (at the icosahedron's corners) and 10 f^2 - 10 hexagons, 10 f^2 + 2 cells in all.
// This is the grid of several global atmosphere models (e.g. Heikes & Randall 1995,
// Monthly Weather Review 123) and of geodesic domes.
//
//   f        cells     mean spacing on Earth
//   16       2 562     ~ 450 km
//   32      10 242     ~ 220 km
//   48      23 042     ~ 150 km
//   64      40 962     ~ 110 km
//
// It satisfies the mesh contract (mesh.js) with real geometry:
//   dist   great-circle distance between neighbouring centres (R x angle)
//   face   length of the shared Voronoi edge: the arc between the circumcentres of
//          the two triangles on either side of the edge
//   area   spherical area of the Voronoi polygon (sums to 4 pi R^2 to round-off)
// so the finite-volume processes are consistent discretisations on the sphere: the
// tests check the discrete Laplacian against spherical harmonics.
//
// Map coordinates. The contract wants x, y in metres for drawing; here they are
// equirectangular: x = R (longitude + pi) eastwards, y = R (pi/2 - latitude) southwards,
// so the map is 2 pi R by pi R. Geometry never uses them (they jump at the date line
// and squash towards the poles); `pos` holds the true unit vectors.
//
// There are no edges: `edge` is all zero, and water reaches the sea or a closed basin.

const PHI = (1 + Math.sqrt(5)) / 2;

const ICO_VERTICES = [
  [-1, PHI, 0], [1, PHI, 0], [-1, -PHI, 0], [1, -PHI, 0],
  [0, -1, PHI], [0, 1, PHI], [0, -1, -PHI], [0, 1, -PHI],
  [PHI, 0, -1], [PHI, 0, 1], [-PHI, 0, -1], [-PHI, 0, 1],
];
const ICO_FACES = [
  [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
  [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
  [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
  [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
];

const dot3 = (a, i, b, j) => a[i] * b[j] + a[i + 1] * b[j + 1] + a[i + 2] * b[j + 2];

/** Angle (rad) between unit vectors a[i..] and b[j..], accurate at small angles. */
function angleBetween(a, i, b, j) {
  const cx = a[i + 1] * b[j + 2] - a[i + 2] * b[j + 1];
  const cy = a[i + 2] * b[j] - a[i] * b[j + 2];
  const cz = a[i] * b[j + 1] - a[i + 1] * b[j];
  return Math.atan2(Math.sqrt(cx * cx + cy * cy + cz * cz), dot3(a, i, b, j));
}

/** Solid angle of the spherical triangle with unit-vector corners (Van Oosterom & Strackee 1983). */
function solidAngle(a, i, b, j, c, k) {
  const det = a[i] * (b[j + 1] * c[k + 2] - b[j + 2] * c[k + 1]) -
    a[i + 1] * (b[j] * c[k + 2] - b[j + 2] * c[k]) +
    a[i + 2] * (b[j] * c[k + 1] - b[j + 1] * c[k]);
  const den = 1 + dot3(a, i, b, j) + dot3(b, j, c, k) + dot3(c, k, a, i);
  return Math.abs(2 * Math.atan2(det, den));
}

export class SphereMesh {
  /**
   * @param {number} subdivisions f: each icosahedron edge is cut into f pieces
   * @param {number} [radius] planet radius (m), default Earth's
   */
  constructor(subdivisions, radius = 6.371e6) {
    const f = subdivisions | 0;
    if (!(f >= 1)) throw new Error('subdivisions must be at least 1');
    if (!(radius > 0)) throw new Error('radius must be positive');
    this.subdivisions = f;
    this.radius = radius;
    const R = radius;

    // --- vertices and triangles: subdivide each face, sharing vertices along edges ---
    const pos = [];
    const key = new Map();
    const vertex = (x, y, z) => {
      const L = Math.hypot(x, y, z);
      x /= L; y /= L; z /= L;
      const k = `${Math.round(x * 1e9)},${Math.round(y * 1e9)},${Math.round(z * 1e9)}`;
      let id = key.get(k);
      if (id === undefined) { id = pos.length / 3; pos.push(x, y, z); key.set(k, id); }
      return id;
    };
    const tris = [];
    for (const [ia, ib, ic] of ICO_FACES) {
      const A = ICO_VERTICES[ia], B = ICO_VERTICES[ib], C = ICO_VERTICES[ic];
      const idx = [];
      for (let i = 0; i <= f; i++) {
        idx.push([]);
        for (let j = 0; j <= f - i; j++) {
          const u = i / f, v = j / f, w = 1 - u - v;
          idx[i].push(vertex(w * A[0] + u * B[0] + v * C[0], w * A[1] + u * B[1] + v * C[1], w * A[2] + u * B[2] + v * C[2]));
        }
      }
      for (let i = 0; i < f; i++) {
        for (let j = 0; j < f - i; j++) {
          tris.push(idx[i][j], idx[i + 1][j], idx[i][j + 1]);
          if (j < f - i - 1) tris.push(idx[i + 1][j], idx[i + 1][j + 1], idx[i][j + 1]);
        }
      }
    }
    const n = pos.length / 3;
    this.size = n;
    this.pos = Float64Array.from(pos);
    const P = this.pos;
    // orient every triangle counter-clockwise seen from outside
    const T = tris.length / 3;
    for (let t = 0; t < T; t++) {
      const a = 3 * tris[3 * t], b = 3 * tris[3 * t + 1], c = 3 * tris[3 * t + 2];
      const ux = P[b] - P[a], uy = P[b + 1] - P[a + 1], uz = P[b + 2] - P[a + 2];
      const vx = P[c] - P[a], vy = P[c + 1] - P[a + 1], vz = P[c + 2] - P[a + 2];
      const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      if (nx * P[a] + ny * P[a + 1] + nz * P[a + 2] < 0) { const s = tris[3 * t + 1]; tris[3 * t + 1] = tris[3 * t + 2]; tris[3 * t + 2] = s; }
    }
    this.triangles = Int32Array.from(tris);

    // circumcentres of the triangles = Voronoi vertices
    const cc = new Float64Array(3 * T);
    for (let t = 0; t < T; t++) {
      const a = 3 * tris[3 * t], b = 3 * tris[3 * t + 1], c = 3 * tris[3 * t + 2];
      const ux = P[b] - P[a], uy = P[b + 1] - P[a + 1], uz = P[b + 2] - P[a + 2];
      const vx = P[c] - P[a], vy = P[c + 1] - P[a + 1], vz = P[c + 2] - P[a + 2];
      let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      const L = Math.hypot(nx, ny, nz);
      cc[3 * t] = nx / L; cc[3 * t + 1] = ny / L; cc[3 * t + 2] = nz / L;
    }

    // --- neighbours in counter-clockwise order, from "after b around a comes c" ---
    const D = 6;
    this.degree = D;
    const nextAround = new Map(); // key a*n+b -> [c, triangle]
    for (let t = 0; t < T; t++) {
      const v = [tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]];
      for (let s = 0; s < 3; s++) nextAround.set(v[s] * n + v[(s + 1) % 3], [v[(s + 2) % 3], t]);
    }
    const first = new Int32Array(n).fill(-1);
    for (let t = 0; t < T; t++) for (let s = 0; s < 3; s++) { const a = tris[3 * t + s]; if (first[a] < 0) first[a] = tris[3 * t + (s + 1) % 3]; }
    this.neighbors = new Int32Array(D * n).fill(-1);
    /** triangle between neighbour slot k and slot k+1 (cyclically), -1 for unused slots */
    this.cornerTriangle = new Int32Array(D * n).fill(-1);
    this.valence = new Uint8Array(n);
    for (let a = 0; a < n; a++) {
      let b = first[a];
      let k = 0;
      do {
        if (k >= D) throw new Error('geodesic mesh: vertex with more than six neighbours');
        const [c, t] = nextAround.get(a * n + b);
        this.neighbors[D * a + k] = b;
        this.cornerTriangle[D * a + k] = t;
        k++;
        b = c;
      } while (b !== first[a]);
      this.valence[a] = k;
    }

    // --- geometry: distances, Voronoi face lengths, cell areas ---
    this.dist = new Float64Array(D * n);
    this.face = new Float64Array(D * n);
    this.area = new Float64Array(n);
    for (let a = 0; a < n; a++) {
      const d = this.valence[a];
      let omega = 0;
      for (let k = 0; k < d; k++) {
        const e = D * a + k;
        const b = this.neighbors[e];
        this.dist[e] = R * angleBetween(P, 3 * a, P, 3 * b);
        // the edge a-b is flanked by the corner triangles before and after slot k
        const tPrev = this.cornerTriangle[D * a + ((k + d - 1) % d)];
        const tNext = this.cornerTriangle[e];
        this.face[e] = R * angleBetween(cc, 3 * tPrev, cc, 3 * tNext);
        omega += solidAngle(P, 3 * a, cc, 3 * tPrev, cc, 3 * tNext);
      }
      this.area[a] = R * R * omega;
    }
    this.circumcenters = cc;
    this.edge = new Uint8Array(n);
    this.lattice = false;
    this.sphere = true;
    let s = 0, m = 0;
    for (let e = 0; e < D * n; e++) if (this.neighbors[e] >= 0) { s += this.dist[e]; m++; }
    /** mean centre-to-centre spacing (m) */
    this.spacing = s / m;
    this.cellArea = (4 * Math.PI * R * R) / n; // mean, for display scales

    // --- map coordinates ---
    this.lon = new Float64Array(n);
    this.lat = new Float64Array(n);
    this.x = new Float64Array(n);
    this.y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      this.lon[i] = Math.atan2(P[3 * i + 1], P[3 * i]);
      this.lat[i] = Math.asin(Math.max(-1, Math.min(1, P[3 * i + 2])));
      this.x[i] = R * (this.lon[i] + Math.PI);
      this.y[i] = R * (Math.PI / 2 - this.lat[i]);
    }
    this.width = 2 * Math.PI * R;
    this.height = Math.PI * R;
    this.bounds = { xmin: 0, ymin: 0, xmax: this.width, ymax: this.height };

    // coarse longitude-latitude buckets give every search a nearby starting cell
    this._bw = 4 * f;
    this._bh = 2 * f;
    this._bucket = new Int32Array(this._bw * this._bh);
    let hint = 0;
    for (let r = 0; r < this._bh; r++) {
      for (let c = 0; c < this._bw; c++) {
        const lon = ((c + 0.5) / this._bw) * 2 * Math.PI - Math.PI;
        const lat = Math.PI / 2 - ((r + 0.5) / this._bh) * Math.PI;
        hint = this.nearest(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat), hint);
        this._bucket[r * this._bw + c] = hint;
      }
    }
  }

  /**
   * Cell whose centre is nearest the unit vector (vx, vy, vz): greedy walk from `start`
   * towards the point. On a Delaunay triangulation the walk cannot get stuck short of
   * the nearest vertex.
   */
  nearest(vx, vy, vz, start = -1) {
    const P = this.pos;
    const nb = this.neighbors;
    let i = start >= 0 ? start : this._startFor(vx, vy, vz);
    let best = P[3 * i] * vx + P[3 * i + 1] * vy + P[3 * i + 2] * vz;
    for (;;) {
      let next = -1;
      for (let k = 0; k < 6; k++) {
        const j = nb[6 * i + k];
        if (j < 0) continue;
        const d = P[3 * j] * vx + P[3 * j + 1] * vy + P[3 * j + 2] * vz;
        if (d > best) { best = d; next = j; }
      }
      if (next < 0) return i;
      i = next;
    }
  }

  _startFor(vx, vy, vz) {
    if (!this._bucket) return 0;
    const lon = Math.atan2(vy, vx);
    const lat = Math.asin(Math.max(-1, Math.min(1, vz)));
    const c = Math.min(this._bw - 1, Math.floor(((lon + Math.PI) / (2 * Math.PI)) * this._bw));
    const r = Math.min(this._bh - 1, Math.floor(((Math.PI / 2 - lat) / Math.PI) * this._bh));
    return this._bucket[r * this._bw + c];
  }

  /** Cell at map coordinates (x, y) (see the header), or -1 off the map. */
  cellAt(x, y) {
    if (!(x >= 0 && x <= this.width && y >= 0 && y <= this.height)) return -1;
    const lon = x / this.radius - Math.PI;
    const lat = Math.PI / 2 - y / this.radius;
    return this.nearest(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat));
  }

  /**
   * Linear interpolation weights at the unit vector v: the triangle of cell centres that
   * contains it and barycentric weights (summing to 1) in the plane of that triangle.
   * @returns {{cell:number, ids:number[], weights:number[]}} cell = nearest centre
   */
  locate(vx, vy, vz, start = -1) {
    const P = this.pos;
    const tris = this.triangles;
    const p = this.nearest(vx, vy, vz, start);
    let bestT = -1, bestMin = -Infinity, bestW = null;
    for (let k = 0; k < this.valence[p]; k++) {
      const t = this.cornerTriangle[6 * p + k];
      const a = 3 * tris[3 * t], b = 3 * tris[3 * t + 1], c = 3 * tris[3 * t + 2];
      // v = wa A + wb B + wc C solved by triple products (Cramer's rule)
      const det = triple(P, a, b, c, P[a], P[a + 1], P[a + 2], 0);
      const wa = triple(P, -1, b, c, vx, vy, vz, 1) / det;
      const wb = triple(P, a, -1, c, vx, vy, vz, 2) / det;
      const wc = triple(P, a, b, -1, vx, vy, vz, 3) / det;
      const lo = Math.min(wa, wb, wc);
      if (lo > bestMin) { bestMin = lo; bestT = t; bestW = [wa, wb, wc]; }
      if (lo >= 0) break;
    }
    const s = bestW[0] + bestW[1] + bestW[2];
    return {
      cell: p,
      ids: [tris[3 * bestT], tris[3 * bestT + 1], tris[3 * bestT + 2]],
      weights: bestW.map((w) => w / s),
    };
  }

  /**
   * Displacement from cell i to cell j in i's local tangent plane (m): x east, y south,
   * matching the map convention. Used for gradients (hillshading).
   */
  offset(i, j, out = [0, 0]) {
    const P = this.pos;
    const px = P[3 * i], py = P[3 * i + 1], pz = P[3 * i + 2];
    // east = z x p (normalised); north = p x east
    let ex = -py, ey = px;
    const el = Math.hypot(ex, ey) || 1;
    ex /= el; ey /= el;
    const nx = -pz * ey, ny = pz * ex, nz = px * ey - py * ex;
    const dx = P[3 * j] - px, dy = P[3 * j + 1] - py, dz = P[3 * j + 2] - pz;
    out[0] = this.radius * (dx * ex + dy * ey);
    out[1] = -this.radius * (dx * nx + dy * ny + dz * nz);
    return out;
  }
}

/** Triple product with one of the three columns optionally replaced by (vx, vy, vz). */
function triple(P, a, b, c, vx, vy, vz, which) {
  const A = which === 1 ? [vx, vy, vz] : [P[a], P[a + 1], P[a + 2]];
  const B = which === 2 ? [vx, vy, vz] : [P[b], P[b + 1], P[b + 2]];
  const C = which === 3 ? [vx, vy, vz] : [P[c], P[c + 1], P[c + 2]];
  return A[0] * (B[1] * C[2] - B[2] * C[1]) - A[1] * (B[0] * C[2] - B[2] * C[0]) + A[2] * (B[0] * C[1] - B[1] * C[0]);
}

/** Unit vector of a longitude and latitude given in degrees. */
export function lonLatToVector(lonDeg, latDeg) {
  const lon = (lonDeg * Math.PI) / 180, lat = (latDeg * Math.PI) / 180;
  return [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
}

/** Longitude and latitude (degrees) of a unit vector. */
export function vectorToLonLat(v) {
  return [(Math.atan2(v[1], v[0]) * 180) / Math.PI, (Math.asin(Math.max(-1, Math.min(1, v[2]))) * 180) / Math.PI];
}

/**
 * Azimuthal equidistant projection centred on the unit vector `center`: distances and
 * directions from the centre are true. x points east and y south (the grid convention).
 * Used to lay a flat hex patch onto the sphere when zooming in from a planet.
 */
export function azimuthalProjection(center, radius) {
  const [cx, cy, cz] = center;
  let ex = -cy, ey = cx, ez = 0; // east = z x c
  let el = Math.hypot(ex, ey, ez);
  if (el < 1e-9) { ex = 0; ey = 1; ez = 0; el = 1; } // at a pole: pick y as "east"
  ex /= el; ey /= el; ez /= el;
  const nx = cy * ez - cz * ey, ny = cz * ex - cx * ez, nz = cx * ey - cy * ex; // north = c x east
  return {
    /** unit vector -> [x, y] (m), or null on the far side of the planet */
    forward(v) {
      const d = v[0] * cx + v[1] * cy + v[2] * cz;
      const e = v[0] * ex + v[1] * ey + v[2] * ez;
      const nn = v[0] * nx + v[1] * ny + v[2] * nz;
      const s = Math.hypot(e, nn);
      const c = Math.atan2(s, d);
      if (c > Math.PI * 0.75) return null;
      if (s < 1e-15) return [0, 0];
      return [(radius * c * e) / s, (-radius * c * nn) / s];
    },
    /** [x, y] (m) -> unit vector */
    inverse(x, y) {
      const c = Math.hypot(x, y) / radius;
      if (c < 1e-15) return [cx, cy, cz];
      const se = x / (c * radius), sn = -y / (c * radius); // direction cosines east and north
      const sc = Math.sin(c), co = Math.cos(c);
      return [co * cx + sc * (se * ex + sn * nx), co * cy + sc * (se * ey + sn * ny), co * cz + sc * (se * ez + sn * nz)];
    },
  };
}
