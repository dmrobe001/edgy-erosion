// Stateless random numbers: every draw is a pure function of (key, counter).
//
// The world must be a pure function of its seed and a tile's address, whatever order
// tiles are generated in. A stateful generator (one that advances with every draw)
// breaks that the moment two tiles are made in a different order, so here nothing has
// state: a key is derived by hashing the world seed with an entity's address, and a
// draw is the hash of that key with a counter (which cell, which draw).
//
// philox4x32 is the counter-based generator Philox4x32-10 of Salmon, Moraes, Dror & Shaw
// (2011, "Parallel random numbers: as easy as 1, 2, 3", Proc. SC11). It passes the
// BigCrush statistical tests and is checked here against the Random123 known-answer
// vectors. It costs ten 32x32-bit multiplies per four words, so lattice noise, which
// needs millions of draws per tile, uses the cheaper integer mix `mix32` instead.
//
// Everything is 32-bit integer arithmetic (Math.imul, shifts, xor), which JavaScript
// evaluates exactly on every engine. Floating-point results downstream are reproducible
// across machines only to the extent the engine's Math.exp, Math.sin etc. are; V8 and
// SpiderMonkey both use fdlibm-derived implementations, but the standard does not
// promise it. Golden tests (test/world.test.js) catch any drift.

const M0 = 0xd2511f53, M1 = 0xcd9e8d57, W0 = 0x9e3779b9, W1 = 0xbb67ae85;

/** High 32 bits of the 64-bit product of two uint32s, exactly. */
function mulhi32(a, b) {
  const al = a & 0xffff, ah = a >>> 16, bl = b & 0xffff, bh = b >>> 16;
  const ll = al * bl, lh = al * bh, hl = ah * bl, hh = ah * bh;
  const mid = (ll >>> 16) + (lh & 0xffff) + (hl & 0xffff);
  return (hh + (lh >>> 16) + (hl >>> 16) + Math.floor(mid / 65536)) >>> 0;
}

/**
 * Philox4x32-10: four uint32 counter words and two key words to four uint32 outputs.
 * @param {Uint32Array} [out]
 */
export function philox4x32(c0, c1, c2, c3, k0, k1, out = new Uint32Array(4)) {
  c0 >>>= 0; c1 >>>= 0; c2 >>>= 0; c3 >>>= 0; k0 >>>= 0; k1 >>>= 0;
  for (let r = 0; r < 10; r++) {
    if (r > 0) { k0 = (k0 + W0) >>> 0; k1 = (k1 + W1) >>> 0; }
    const hi0 = mulhi32(M0, c0), lo0 = Math.imul(M0, c0) >>> 0;
    const hi1 = mulhi32(M1, c2), lo1 = Math.imul(M1, c2) >>> 0;
    const n0 = (hi1 ^ c1 ^ k0) >>> 0, n2 = (hi0 ^ c3 ^ k1) >>> 0;
    c0 = n0; c1 = lo1; c2 = n2; c3 = lo0;
  }
  out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
  return out;
}

/** Split an integer (|v| < 2^53) into two uint32 words, low first; negatives are distinct from positives. */
function words(v) {
  if (!Number.isInteger(v)) throw new Error(`hash parts must be integers, got ${v}`);
  const neg = v < 0;
  const a = Math.abs(v);
  const lo = a >>> 0, hi = Math.floor(a / 4294967296) >>> 0;
  return [lo, (hi ^ (neg ? 0x80000000 : 0)) >>> 0];
}

/** A world seed (integer, string or bigint) as a two-word key. */
export function worldKey(seed) {
  if (typeof seed === 'bigint') return [Number(seed & 0xffffffffn), Number((seed >> 32n) & 0xffffffffn)];
  if (typeof seed === 'string') {
    // FNV-1a over the UTF-16 code units, two lanes
    let a = 0x811c9dc5, b = 0x01000193 ^ seed.length;
    for (let k = 0; k < seed.length; k++) {
      a = Math.imul(a ^ seed.charCodeAt(k), 0x01000193);
      b = Math.imul(b ^ seed.charCodeAt(k), 0x5bd1e995) ^ (b >>> 15);
    }
    return [a >>> 0, b >>> 0];
  }
  return words(seed);
}

/**
 * Key of an entity: hash(parent key, integer parts...). The parts are absorbed one at a
 * time through Philox, so [a, b] and [b, a] (and [a] and [a, 0]) give different keys.
 * @param {number[]} key two uint32 words
 * @param {...number} parts integers, e.g. level, face, i, j
 */
export function deriveKey(key, ...parts) {
  let k0 = key[0] >>> 0, k1 = key[1] >>> 0;
  const out = new Uint32Array(4);
  for (let p = 0; p < parts.length; p++) {
    const [lo, hi] = words(parts[p]);
    philox4x32(lo, hi, p, 0x6b657921 /* "key!" */, k0, k1, out);
    k0 = out[0]; k1 = out[1];
  }
  philox4x32(parts.length, 0, 0, 0x656e6421 /* "end!" */, k0, k1, out);
  return [out[0], out[1]];
}

/** Key as a 16-digit hex string (for logs and cache keys). */
export function keyHex(key) {
  return key.map((w) => (w >>> 0).toString(16).padStart(8, '0')).join('');
}

const scratch = new Uint32Array(4);

/**
 * Uniform draw in [0, 1) with 53 random bits, a pure function of (key, a, b, c).
 * Typical use: a = cell index, b = draw number, c = purpose tag.
 */
export function uniform(key, a = 0, b = 0, c = 0) {
  philox4x32(a, b, c, 0, key[0], key[1], scratch);
  return (scratch[0] * 2097152 + (scratch[1] >>> 11)) / 9007199254740992;
}

/** Standard normal draw (Box-Muller on two uniforms), a pure function of (key, a, b, c). */
export function normal(key, a = 0, b = 0, c = 0) {
  philox4x32(a, b, c, 1, key[0], key[1], scratch);
  const u1 = (scratch[0] + 1) / 4294967297, u2 = scratch[1] / 4294967296;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

/**
 * Fast integer mix of a key and three integers to a uint32 (the murmur3 finalizer applied
 * to a combination). Not as strong as Philox; for lattice noise, where it is evaluated
 * millions of times per tile.
 */
export function mix32(key, a, b, c = 0) {
  let h = (key[0] ^ Math.imul(a | 0, 0x27d4eb2d)) >>> 0;
  h = Math.imul(h ^ (h >>> 15), 0x85ebca6b) ^ Math.imul(b | 0, 0x165667b1) ^ key[1];
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35) ^ Math.imul(c | 0, 0x1b873593);
  h = Math.imul(h ^ (h >>> 16), 0x85ebca6b);
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35);
  return (h ^ (h >>> 16)) >>> 0;
}

/**
 * Hash of a typed array's bytes (FNV-1a, 32 bits, as 8 hex digits): golden tests store
 * these to detect any change in generated fields.
 * @param {ArrayBufferView} arr
 */
export function fieldHash(arr) {
  const bytes = new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength);
  let h = 0x811c9dc5;
  for (let k = 0; k < bytes.length; k++) h = Math.imul(h ^ bytes[k], 0x01000193);
  return (h >>> 0).toString(16).padStart(8, '0');
}
