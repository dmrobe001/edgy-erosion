// Seeded pseudo-random numbers and smooth noise.
//
// Everything stochastic in the simulation draws from these, so a run is fully
// determined by its parameters (including `seed`). Math.random is never used.

/**
 * Mulberry32: a tiny 32-bit generator. Not cryptographic; plenty for terrain.
 * The whole state is one uint32, so it can be saved in a snapshot and resumed.
 */
export class Rng {
  constructor(seed = 1) {
    this.state = (seed >>> 0) || 0x9e3779b9;
  }

  /** Uniform float in [0, 1). */
  next() {
    let t = (this.state = (this.state + 0x6d2b79f5) >>> 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /** Uniform float in [lo, hi). */
  uniform(lo, hi) {
    return lo + (hi - lo) * this.next();
  }
}

/** Integer hash of a lattice point and a seed, returning a float in [-1, 1]. */
function hash2(ix, iy, seed) {
  let h = Math.imul(ix, 0x27d4eb2d) ^ Math.imul(iy, 0x165667b1) ^ Math.imul(seed, 0x9e3779b1);
  h = Math.imul(h ^ (h >>> 15), 0x85ebca6b);
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35);
  h ^= h >>> 16;
  return (h >>> 0) / 2147483647.5 - 1;
}

const fade = (t) => t * t * t * (t * (t * 6 - 15) + 10);

/** Value noise on the unit square lattice, smooth (C2) interpolation, range ~[-1, 1]. */
export function valueNoise(x, y, seed) {
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const u = fade(x - x0);
  const v = fade(y - y0);
  const a = hash2(x0, y0, seed);
  const b = hash2(x0 + 1, y0, seed);
  const c = hash2(x0, y0 + 1, seed);
  const d = hash2(x0 + 1, y0 + 1, seed);
  return a + (b - a) * u + (c - a) * v + (a - b - c + d) * u * v;
}

/**
 * Fractal (fractional Brownian motion) noise: a sum of octaves of value noise,
 * each with half the wavelength and `gain` times the amplitude of the last.
 * Normalised so the result stays roughly in [-1, 1].
 *
 * @param {number} x world coordinate (m)
 * @param {number} y world coordinate (m)
 * @param {{wavelength:number, octaves?:number, gain?:number, seed?:number}} o
 */
export function fbm(x, y, { wavelength, octaves = 4, gain = 0.5, seed = 0 }) {
  let sum = 0;
  let amp = 1;
  let norm = 0;
  let f = 1 / wavelength;
  for (let k = 0; k < octaves; k++) {
    sum += amp * valueNoise(x * f, y * f, seed + 1013 * k);
    norm += amp;
    amp *= gain;
    f *= 2;
  }
  return sum / norm;
}
