// Tile addresses on a cube-sphere: the sphere is the six faces of a cube blown up onto
// it, and each face is cut into a square grid of tiles at every level.
//
// Level 0 is the whole planet (one tile, on the icosahedral mesh: face -1). From level 1
// on, a face is cut into `rootTiles` x `rootTiles` tiles, and every tile at level L is
// cut into `branching` x `branching` children at level L + 1. Each tile is a raster of
// `tileCells` x `tileCells` cells, so the cell size shrinks by `branching` per level and
// a child tile spans tileCells / branching of its parent's cells.
//
// Face coordinates (a, b) run over [0, 1] x [0, 1]. The map from face to sphere uses
// the equal-angle ("gnomonic with tangent warp") projection: a is linear in the angle
// about the face's axis, not in the cube's plane, which keeps cells within 1.41x of
// each other in width over a face (the plain cube projection would reach 1.73x...5.2x
// in area). Tile i runs along a, tile j along b.
//
// Neighbours, edges and corners are found by geometry, not by tables: step a little
// beyond the tile's edge on the sphere and ask which tile contains that point. That
// handles the twelve cube edges, where the neighbour's axes are rotated, and the eight
// cube corners, where only three tiles meet, without special cases.

/**
 * Default layout, on Earth: level-1 tiles ~310 km across with 1.2-km cells, level 2 ~9.8 km
 * with 38-m cells, level 3 ~305 m with 1.2-m cells. The design's cell sizes with 256-cell
 * tiles instead of 1024 (four times as many tiles per level, each 16x cheaper; see
 * docs/PLAN.md for why).
 */
export const WORLD_LAYOUT = { rootTiles: 32, branching: 32, tileCells: 256 };

// face: normal n, and the axes along which a and b increase
const FACES = [
  { n: [1, 0, 0], s: [0, 1, 0], t: [0, 0, 1] },
  { n: [-1, 0, 0], s: [0, -1, 0], t: [0, 0, 1] },
  { n: [0, 1, 0], s: [-1, 0, 0], t: [0, 0, 1] },
  { n: [0, -1, 0], s: [1, 0, 0], t: [0, 0, 1] },
  { n: [0, 0, 1], s: [0, 1, 0], t: [-1, 0, 0] },
  { n: [0, 0, -1], s: [0, 1, 0], t: [1, 0, 0] },
];

/** Unit vector at face coordinates (a, b) (values a little outside [0, 1] extend the face). */
export function faceToVector(face, a, b, out = [0, 0, 0]) {
  const F = FACES[face];
  const u = Math.tan((Math.PI / 4) * (2 * a - 1));
  const v = Math.tan((Math.PI / 4) * (2 * b - 1));
  const x = F.n[0] + u * F.s[0] + v * F.t[0];
  const y = F.n[1] + u * F.s[1] + v * F.t[1];
  const z = F.n[2] + u * F.s[2] + v * F.t[2];
  const L = Math.hypot(x, y, z);
  out[0] = x / L; out[1] = y / L; out[2] = z / L;
  return out;
}

/** Face and face coordinates of a direction (need not be unit length). */
export function vectorToFace(x, y, z) {
  let face = 0, best = -Infinity;
  for (let f = 0; f < 6; f++) {
    const n = FACES[f].n;
    const d = n[0] * x + n[1] * y + n[2] * z;
    if (d > best) { best = d; face = f; }
  }
  const F = FACES[face];
  const u = (F.s[0] * x + F.s[1] * y + F.s[2] * z) / best;
  const v = (F.t[0] * x + F.t[1] * y + F.t[2] * z) / best;
  const a = (Math.atan(u) / (Math.PI / 4) + 1) / 2;
  const b = (Math.atan(v) / (Math.PI / 4) + 1) / 2;
  return { face, a: Math.min(1, Math.max(0, a)), b: Math.min(1, Math.max(0, b)) };
}

/** Tiles per face side at a level. */
export function tilesPerSide(level, layout = WORLD_LAYOUT) {
  if (level < 1) throw new Error('level 0 is the whole planet; cube tiles start at level 1');
  return layout.rootTiles * layout.branching ** (level - 1);
}

export function tileKey(addr) {
  return `${addr.level}/${addr.face}/${addr.i}/${addr.j}`;
}

export function parseTileKey(key) {
  const [level, face, i, j] = key.split('/').map(Number);
  return { level, face, i, j };
}

/** Tile containing a direction at a level. */
export function tileAt(level, x, y, z, layout = WORLD_LAYOUT) {
  const { face, a, b } = vectorToFace(x, y, z);
  const N = tilesPerSide(level, layout);
  return { level, face, i: Math.min(N - 1, Math.floor(a * N)), j: Math.min(N - 1, Math.floor(b * N)) };
}

/** Face-coordinate rectangle of a tile: {a0, b0, size}. */
export function tileRect(addr, layout = WORLD_LAYOUT) {
  const N = tilesPerSide(addr.level, layout);
  return { a0: addr.i / N, b0: addr.j / N, size: 1 / N };
}

/**
 * Unit vector at the centre of cell (ci, cj) of a tile with n x n cells. Cells outside
 * 0..n-1 (a halo) extend the face plane; they are well defined up to the face's edge and
 * a little beyond.
 */
export function cellCenter(addr, n, ci, cj, out, layout = WORLD_LAYOUT) {
  const r = tileRect(addr, layout);
  return faceToVector(addr.face, r.a0 + ((ci + 0.5) / n) * r.size, r.b0 + ((cj + 0.5) / n) * r.size, out);
}

/** The parent tile (level 1 tiles have the planet, level 0, as parent). */
export function parentOf(addr, layout = WORLD_LAYOUT) {
  if (addr.level <= 1) return { level: 0, face: -1, i: 0, j: 0 };
  const b = layout.branching;
  return { level: addr.level - 1, face: addr.face, i: Math.floor(addr.i / b), j: Math.floor(addr.j / b) };
}

/** The branching^2 children. */
export function childrenOf(addr, layout = WORLD_LAYOUT) {
  const out = [];
  if (addr.level === 0) {
    const N = layout.rootTiles;
    for (let f = 0; f < 6; f++) for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) out.push({ level: 1, face: f, i, j });
    return out;
  }
  const b = layout.branching;
  for (let j = 0; j < b; j++) for (let i = 0; i < b; i++) out.push({ level: addr.level + 1, face: addr.face, i: addr.i * b + i, j: addr.j * b + j });
  return out;
}

const SIDES = { W: [-1, 0], E: [1, 0], S: [0, -1], N: [0, 1] };

/** The tile across one side ('W' = i - 1, 'E' = i + 1, 'S' = j - 1, 'N' = j + 1), across cube edges too. */
export function neighborOf(addr, side, layout = WORLD_LAYOUT) {
  const [di, dj] = SIDES[side];
  const N = tilesPerSide(addr.level, layout);
  const i = addr.i + di, j = addr.j + dj;
  if (i >= 0 && i < N && j >= 0 && j < N) return { level: addr.level, face: addr.face, i, j };
  // off the face: the point a quarter-tile beyond the edge midpoint, on the sphere
  const r = tileRect(addr, layout);
  const a = r.a0 + r.size * (0.5 + 0.75 * di), b = r.b0 + r.size * (0.5 + 0.75 * dj);
  const v = faceToVector(addr.face, a, b);
  return tileAt(addr.level, v[0], v[1], v[2], layout);
}

/** Tiles sharing a corner of a tile ('SW', 'SE', 'NW', 'NE'): four, or three at a cube corner. Sorted by key. */
export function cornerTiles(addr, corner, layout = WORLD_LAYOUT) {
  const r = tileRect(addr, layout);
  const a = r.a0 + (corner[1] === 'E' ? r.size : 0), b = r.b0 + (corner[0] === 'N' ? r.size : 0);
  const p = faceToVector(addr.face, a, b);
  // probe around the corner point in its tangent plane
  const e = 0.25 * r.size * (Math.PI / 2);
  const ref = Math.abs(p[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  const s = normalize(cross(p, ref)), t = cross(p, s);
  const seen = new Map();
  for (let k = 0; k < 8; k++) {
    const ang = (2 * Math.PI * (k + 0.5)) / 8;
    const c = Math.cos(ang) * e, d = Math.sin(ang) * e;
    const tile = tileAt(addr.level, p[0] + c * s[0] + d * t[0], p[1] + c * s[1] + d * t[1], p[2] + c * s[2] + d * t[2], layout);
    seen.set(tileKey(tile), tile);
  }
  return [...seen.keys()].sort().map((k) => seen.get(k));
}

/** Canonical identity of the edge between a tile and its neighbour across `side`: the same from both sides. */
export function edgeKey(addr, side, layout = WORLD_LAYOUT) {
  const keys = [tileKey(addr), tileKey(neighborOf(addr, side, layout))].sort();
  return `E:${keys[0]}|${keys[1]}`;
}

/** Canonical identity of a corner: the sorted keys of the tiles that share it. */
export function cornerKey(addr, corner, layout = WORLD_LAYOUT) {
  return `C:${cornerTiles(addr, corner, layout).map(tileKey).join('|')}`;
}

function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
function normalize(a) {
  const L = Math.hypot(a[0], a[1], a[2]);
  return [a[0] / L, a[1] / L, a[2] / L];
}
