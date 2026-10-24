// Stateless lattice noise in 3-D, keyed: a pure function of (key, position).
//
// Noise for a tile is evaluated at each cell's position on the sphere (metres along the
// unit vector times the radius), never in tile coordinates, so two tiles that share an
// edge, or a tile and its halo, see exactly the same noise there.

import { mix32 } from './hash.js';

const fade = (t) => t * t * t * (t * (t * 6 - 15) + 10);
const lattice = (key, x, y, z) => mix32(key, x, y, z) / 2147483647.5 - 1;

/** Value noise on the unit cubic lattice, C2 interpolation, range about [-1, 1]. */
export function valueNoise3k(key, x, y, z) {
  const x0 = Math.floor(x), y0 = Math.floor(y), z0 = Math.floor(z);
  const u = fade(x - x0), v = fade(y - y0), w = fade(z - z0);
  const c000 = lattice(key, x0, y0, z0), c100 = lattice(key, x0 + 1, y0, z0);
  const c010 = lattice(key, x0, y0 + 1, z0), c110 = lattice(key, x0 + 1, y0 + 1, z0);
  const c001 = lattice(key, x0, y0, z0 + 1), c101 = lattice(key, x0 + 1, y0, z0 + 1);
  const c011 = lattice(key, x0, y0 + 1, z0 + 1), c111 = lattice(key, x0 + 1, y0 + 1, z0 + 1);
  const a = c000 + (c100 - c000) * u, b = c010 + (c110 - c010) * u;
  const c = c001 + (c101 - c001) * u, d = c011 + (c111 - c011) * u;
  const e = a + (b - a) * v, f = c + (d - c) * v;
  return e + (f - e) * w;
}

/**
 * Fractal noise at a point given in metres: octaves from `wavelength` down by halves,
 * amplitude falling by `gain` per octave, normalised to about [-1, 1].
 */
export function fbm3k(key, x, y, z, wavelength, octaves = 4, gain = 0.5) {
  let sum = 0, amp = 1, norm = 0, f = 1 / wavelength;
  for (let k = 0; k < octaves; k++) {
    // a different key per octave, and an offset so lattice points of octaves never align
    const kk = [key[0] ^ Math.imul(k + 1, 0x9e3779b9), key[1]];
    sum += amp * valueNoise3k(kk, x * f + 17.3 * k, y * f - 5.1 * k, z * f + 9.7 * k);
    norm += amp;
    amp *= gain;
    f *= 2;
  }
  return sum / norm;
}
