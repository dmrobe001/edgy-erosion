// The generation pipeline: any tile on demand, parents first, with caching and prefetch.
//
// A tile is a pure function of (world seed, address), so caching is only a speed-up
// (design: caching): results are kept in a small in-memory LRU and, if a store is given,
// in a persistent store (IndexedDB in the browser) under (generator version, world seed,
// tile address). Stored tiles are structured clones of the originals, so a cache hit
// returns exactly what generating would have.
//
// Prefetch: tiles queued with prefetch() are generated one at a time whenever nothing
// else is asked for, most recently queued first.

import { generatePlanet, GENERATOR_VERSION } from './l0.js';
import { regionInput, runRegion, cpuSolver, regionSeams } from './region.js';
import { landscapeSeams, landscapeInput, runLandscape } from './landscape.js';
import { siteSeams, buildSite } from './site.js';
import { parentOf, tileKey, WORLD_LAYOUT } from './cubesphere.js';

/** What a level keeps of a result (enough to render it and to build its children). */
function keepOf(level, r) {
  if (level === 1) {
    return {
      addr: r.addr, n: r.n, nodes: r.nodes, dx: r.dx, fields: r.fields, stats: r.stats,
      drainage: { receiver: r.drainage.receiver, area: r.drainage.area, nx: r.drainage.nx, dx: r.drainage.dx },
      crossings: r.rivers.crossings, riverReaches: r.rivers.edges.length, kFactor: r.kFactor,
    };
  }
  if (level === 2) {
    return {
      addr: r.addr, n: r.n, nodes: r.nodes, dx: r.dx, fields: r.fields, stats: r.stats, record: r.record,
      drainage: { receiver: r.drainage.receiver, area: r.drainage.area, nx: r.drainage.nx, dx: r.drainage.dx },
    };
  }
  return { addr: r.addr, n: r.n, nodes: r.nodes, dx: r.dx, fields: r.fields, stats: r.stats, objects: r.objects };
}

export class World {
  /**
   * @param {number|string} seed
   * @param {object} [o]
   * @param {object} [o.layout]
   * @param {(input:object) => object} [o.regionSolver] solver factory for level 1 (default: CPU reference)
   * @param {{get(key:string):Promise<any>, put(key:string, v:any):Promise<void>}} [o.store]
   * @param {object} [o.planetOptions]
   * @param {number} [o.memory] tiles kept in memory
   */
  constructor(seed, o = {}) {
    this.seed = String(seed);
    this.layout = o.layout || WORLD_LAYOUT;
    this.regionSolver = o.regionSolver || cpuSolver;
    this.store = o.store || null;
    this.planetOptions = o.planetOptions || {};
    this.memory = o.memory || 24;
    this.lru = new Map();
    this.pending = new Map();
    this.queue = [];
    this.busy = 0;
    this.stats = { generated: 0, memoryHits: 0, storeHits: 0, timings: [] };
    this._seams = new Map();
  }

  key(addr) { return `${GENERATOR_VERSION}|${this.seed}|${tileKey(addr)}`; }

  /** The planet (level 0), built once. */
  planet() {
    if (!this._planet) {
      const t0 = now();
      this._planet = generatePlanet(this.seed, this.planetOptions);
      this.stats.timings.push({ tile: '0', seconds: (now() - t0) / 1000, source: 'generated' });
    }
    return this._planet;
  }

  _remember(k, v) {
    this.lru.delete(k);
    this.lru.set(k, v);
    while (this.lru.size > this.memory) this.lru.delete(this.lru.keys().next().value);
  }

  /** A tile of level 1, 2 or 3; its parents are built (or fetched) first. */
  async tile(addr) {
    const k = this.key(addr);
    if (this.lru.has(k)) { this.stats.memoryHits++; const v = this.lru.get(k); this._remember(k, v); return v; }
    if (this.pending.has(k)) return this.pending.get(k);
    const job = (async () => {
      this.busy++;
      try {
        if (this.store) {
          const v = await this.store.get(k);
          if (v) { this.stats.storeHits++; this._remember(k, v); this.stats.timings.push({ tile: tileKey(addr), seconds: 0, source: 'store' }); return v; }
        }
        const t0 = now();
        const v = await this._generate(addr);
        this.stats.generated++;
        this.stats.timings.push({ tile: tileKey(addr), seconds: (now() - t0) / 1000, source: 'generated' });
        this._remember(k, v);
        if (this.store) await this.store.put(k, v);
        return v;
      } finally {
        this.busy--;
        this.pending.delete(k);
        this._pump();
      }
    })();
    this.pending.set(k, job);
    return job;
  }

  async _generate(addr) {
    const planet = this.planet();
    const L = this.layout;
    if (addr.level === 1) {
      const input = regionInput(planet, addr, {}, L, this._regionSeams());
      const solver = await this.regionSolver(input);
      const r = await runRegion(input, solver);
      solver.destroy?.();
      return keepOf(1, r);
    }
    const parentAddr = parentOf(addr, L);
    const parent = await this.tile(parentAddr);
    if (addr.level === 2) {
      const ctx = this._landscapeSeams(parent);
      return keepOf(2, runLandscape(landscapeInput(parent, ctx, addr, planet.events, {}, L)));
    }
    if (addr.level === 3) {
      const grand = await this.tile(parentOf(parentAddr, L));
      const ctx2 = this._landscapeSeams(grand);
      return keepOf(3, buildSite(parent, this._siteSeams(parent, ctx2), addr, {}, L));
    }
    throw new Error(`no level ${addr.level}`);
  }

  _regionSeams() {
    if (!this._s1) this._s1 = regionSeams(this.planet(), {}, this.layout);
    return this._s1;
  }

  _landscapeSeams(l1) {
    const k = `2:${tileKey(l1.addr)}`;
    if (!this._seams.has(k)) this._seams.set(k, landscapeSeams(l1, this._regionSeams(), this.planet().worldSeed, {}, this.layout));
    return this._seams.get(k);
  }

  _siteSeams(l2, ctx2) {
    const k = `3:${tileKey(l2.addr)}`;
    if (!this._seams.has(k)) this._seams.set(k, siteSeams(l2, ctx2, this.planet().worldSeed, {}, this.layout));
    return this._seams.get(k);
  }

  /** Queue tiles to build in the background when nothing else is running. */
  prefetch(addrs) {
    for (const a of addrs) {
      const k = this.key(a);
      if (this.lru.has(k) || this.pending.has(k)) continue;
      // newest first: what was asked for last is nearest to where the viewer is now
      this.queue = this.queue.filter((q) => this.key(q) !== k);
      this.queue.unshift(a);
    }
    this._pump();
  }

  _pump() {
    if (this.busy || !this.queue.length) return;
    const next = this.queue.shift();
    // run on a later turn, so a foreground request issued now goes first
    setTimeout(() => { if (!this.busy) this.tile(next).catch(() => {}); else this.queue.unshift(next); }, 0);
  }
}

const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

/** An IndexedDB-backed store for World (browser only). */
export function indexedDBStore(name = 'hierarchical-terrain') {
  const open = new Promise((resolve, reject) => {
    const req = indexedDB.open(name, 1);
    req.onupgradeneeded = () => req.result.createObjectStore('tiles');
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
  const tx = async (mode, fn) => {
    const db = await open;
    return new Promise((resolve, reject) => {
      const t = db.transaction('tiles', mode);
      const r = fn(t.objectStore('tiles'));
      t.oncomplete = () => resolve(r.result);
      t.onerror = () => reject(t.error);
    });
  };
  return {
    get: (k) => tx('readonly', (s) => s.get(k)),
    put: (k, v) => tx('readwrite', (s) => s.put(v, k)),
    clear: () => tx('readwrite', (s) => s.clear()),
  };
}

/** A store in memory (for tests and Node). */
export function memoryStore() {
  const m = new Map();
  return { get: async (k) => (m.has(k) ? structuredClone(m.get(k)) : undefined), put: async (k, v) => { m.set(k, structuredClone(v)); }, map: m };
}
