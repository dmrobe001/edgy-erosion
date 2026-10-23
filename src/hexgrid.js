// Hexagonal lattice geometry.
//
// Layout: "pointy-top" hexagons stored in rows, with odd rows shifted right by
// half a cell ("odd-r" offset coordinates). Cell index i = row * cols + col.
//
//   row 0:   0   1   2   3
//   row 1:     4   5   6   7
//   row 2:   8   9  10  11
//
// Every cell centre is exactly `dx` from each of its six neighbours, so the
// cell centres form a triangular lattice and the cells are its Voronoi
// (hexagonal) tiles:
//
//   centre of (col, row):  x = dx * (col + (row & 1) / 2),  y = dx * (sqrt(3)/2) * row
//   hexagon side length (= circumradius):  a = dx / sqrt(3)
//   cell area:                             (sqrt(3) / 2) * dx^2
//
// Internally we also use axial coordinates (q, r) where the six neighbour
// directions are simple constant offsets; see https://www.redblobgames.com/grids/hexagons/
//
// HexGrid implements the general mesh contract in mesh.js (degree, area, dist,
// face, edge flags, bounds), so the processes never assume hexagons. It adds
// `lattice = true`: neighbour slot k always points along azimuth 60k degrees, which
// the wave-exposure sweep relies on. `originX/originY` place cell (0, 0) in world
// coordinates, so nested patches share their parent's coordinate system.

const SQRT3 = Math.sqrt(3);

/** Axial direction vectors, counter-clockwise starting from east. */
export const AXIAL_DIRECTIONS = [
  [1, 0], [1, -1], [0, -1], [-1, 0], [-1, 1], [0, 1],
];

/** Index of the opposite direction: OPPOSITE[k] points back along direction k. */
export const OPPOSITE = [3, 4, 5, 0, 1, 2];

export class HexGrid {
  /**
   * @param {number} cols number of cells per row
   * @param {number} rows number of rows
   * @param {number} dx centre-to-centre lattice spacing (m)
   * @param {number} [originX] world x of the centre of cell (0, 0) (m)
   * @param {number} [originY] world y of the centre of cell (0, 0) (m)
   */
  constructor(cols, rows, dx, originX = 0, originY = 0) {
    if (!(cols >= 3 && rows >= 3)) throw new Error('grid must be at least 3 x 3');
    if (!(dx > 0)) throw new Error('lattice spacing dx must be positive');
    this.cols = cols | 0;
    this.rows = rows | 0;
    this.dx = dx;
    this.size = this.cols * this.rows;
    this.rowSpacing = (SQRT3 / 2) * dx;
    this.cellArea = (SQRT3 / 2) * dx * dx;
    this.edgeLength = dx / SQRT3;
    /** Extent of the cell centres, plus half a cell of padding on each side. */
    this.width = dx * (this.cols + 0.5);
    this.height = this.rowSpacing * (this.rows - 1) + 2 * this.edgeLength;
    this.originX = originX;
    this.originY = originY;
    this.bounds = {
      xmin: originX - dx / 2, ymin: originY - this.edgeLength,
      xmax: originX - dx / 2 + this.width, ymax: originY - this.edgeLength + this.height,
    };
    this.degree = 6;
    this.lattice = true;

    const n = this.size;
    this.x = new Float64Array(n);
    this.y = new Float64Array(n);
    /** neighbours[6 * i + k] = index of the neighbour in direction k, or -1 outside the grid. */
    this.neighbors = new Int32Array(6 * n);
    /** Bit flags: 1 = west edge, 2 = east edge, 4 = north (row 0), 8 = south (last row). */
    this.edge = new Uint8Array(n);

    for (let r = 0; r < this.rows; r++) {
      for (let c = 0; c < this.cols; c++) {
        const i = r * this.cols + c;
        this.x[i] = originX + dx * (c + 0.5 * (r & 1));
        this.y[i] = originY + this.rowSpacing * r;
        let flags = 0;
        if (c === 0) flags |= 1;
        if (c === this.cols - 1) flags |= 2;
        if (r === 0) flags |= 4;
        if (r === this.rows - 1) flags |= 8;
        this.edge[i] = flags;
        const [q, rr] = offsetToAxial(c, r);
        for (let k = 0; k < 6; k++) {
          const [nc, nr] = axialToOffset(q + AXIAL_DIRECTIONS[k][0], rr + AXIAL_DIRECTIONS[k][1]);
          this.neighbors[6 * i + k] =
            nc >= 0 && nc < this.cols && nr >= 0 && nr < this.rows ? nr * this.cols + nc : -1;
        }
      }
    }

    // mesh contract: per-cell area, per-edge centre distance and shared face width
    this.area = new Float64Array(n).fill(this.cellArea);
    this.dist = new Float64Array(6 * n);
    this.face = new Float64Array(6 * n);
    for (let e = 0; e < 6 * n; e++) {
      if (this.neighbors[e] >= 0) { this.dist[e] = dx; this.face[e] = this.edgeLength; }
    }
  }

  index(col, row) {
    return row * this.cols + col;
  }

  col(i) {
    return i % this.cols;
  }

  row(i) {
    return (i / this.cols) | 0;
  }

  /**
   * Cell containing the world point (x, y), or -1 if it lies outside the grid.
   * Uses fractional axial coordinates and cube rounding.
   */
  cellAt(x, y) {
    x -= this.originX;
    y -= this.originY;
    const a = this.edgeLength;
    const qf = ((SQRT3 / 3) * x - y / 3) / a;
    const rf = ((2 / 3) * y) / a;
    // cube rounding: (q, r, s) with q + r + s = 0
    const sf = -qf - rf;
    let q = Math.round(qf);
    let r = Math.round(rf);
    const s = Math.round(sf);
    const dq = Math.abs(q - qf);
    const dr = Math.abs(r - rf);
    const ds = Math.abs(s - sf);
    if (dq > dr && dq > ds) q = -r - s;
    else if (dr > ds) r = -q - s;
    const [c, row] = axialToOffset(q, r);
    if (c < 0 || c >= this.cols || row < 0 || row >= this.rows) return -1;
    return row * this.cols + c;
  }

  /** The six corner points of cell i, for drawing. */
  corners(i) {
    const a = this.edgeLength;
    const out = [];
    for (let k = 0; k < 6; k++) {
      const ang = (Math.PI / 180) * (60 * k - 30);
      out.push([this.x[i] + a * Math.cos(ang), this.y[i] + a * Math.sin(ang)]);
    }
    return out;
  }
}

export function offsetToAxial(col, row) {
  return [col - ((row - (row & 1)) >> 1), row];
}

export function axialToOffset(q, r) {
  return [q + ((r - (r & 1)) >> 1), r];
}
