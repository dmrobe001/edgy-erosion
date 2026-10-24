// The world, generated in a worker: tiles on request, rendered to RGBA images.
import { World, indexedDBStore } from '../../src/world/pipeline.js';
import { getDevice, gpuSolver } from '../../src/world/gpu.js';
import { cpuSolver } from '../../src/world/region.js';
import { tileAt, parseTileKey, neighborOf, tileKey, childrenOf, WORLD_LAYOUT } from '../../src/world/cubesphere.js';
import { spherePixelMap } from '../../src/render.js';
import { SITE_COVER } from '../../src/world/site.js';
import { GENERATOR_VERSION } from '../../src/world/l0.js';

let world = null, solverName = 'CPU', pixelMap = null;
const W0 = 1024, H0 = 512;

async function init(seed, useStore) {
  const device = await getDevice().catch(() => null);
  solverName = device ? `WebGPU (${device.adapterInfo?.architecture || device.adapterInfo?.vendor || 'gpu'})` : 'CPU reference (no WebGPU: slow)';
  world = new World(seed, {
    regionSolver: device ? (input) => gpuSolver(device, input) : cpuSolver,
    store: useStore ? indexedDBStore() : null,
  });
}

const lerp = (a, b, t) => a.map((v, k) => v + (b[k] - v) * Math.max(0, Math.min(1, t)));
function hypso(z) {
  if (z <= 0) return lerp([70, 120, 185], [15, 35, 90], -z / 4000);
  if (z < 600) return lerp([95, 145, 80], [165, 170, 105], z / 600);
  if (z < 2200) return lerp([165, 170, 105], [150, 120, 90], (z - 600) / 1600);
  return lerp([150, 120, 90], [245, 245, 245], (z - 2200) / 2500);
}

function renderPlanet() {
  const p = world.planet();
  pixelMap ||= spherePixelMap(p.mesh, W0, H0);
  const rgba = new Uint8ClampedArray(4 * W0 * H0);
  for (let px = 0; px < W0 * H0; px++) {
    const i = pixelMap.map[px];
    const c = i < 0 ? [0, 0, 0] : hypso(p.fields.elevation[i]);
    rgba.set([c[0], c[1], c[2], 255], 4 * px);
  }
  return { w: W0, h: H0, rgba };
}

const COVER_COLORS = {
  sea: [40, 80, 150], water: [60, 115, 195], gravel: [195, 185, 165], meadow: [150, 180, 95], forest: [55, 100, 50],
  shrub: [130, 145, 90], rock: [140, 135, 128], talus: [172, 162, 150], snow: [242, 246, 250], bare: [205, 185, 145], beach: [228, 212, 165],
};

function renderTile(t, level) {
  const m = t.nodes, e = t.fields.elevation, rgba = new Uint8ClampedArray(4 * m * m);
  const A = t.fields.drainageArea;
  const shade = level === 1 ? 6 : level === 2 ? 1.5 : 0.6;
  for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) {
    const q = (m - 1 - j) * m + i, o = 4 * (j * m + i);
    let col;
    let h = e[q], hl = e[(m - 1 - j) * m + Math.max(0, i - 1)], hu = e[Math.max(0, m - 2 - j) * m + i];
    if (level === 3) {
      col = COVER_COLORS[SITE_COVER[t.fields.cover[q]]];
      h += t.fields.canopy[q]; hl += t.fields.canopy[(m - 1 - j) * m + Math.max(0, i - 1)]; hu += t.fields.canopy[Math.max(0, m - 2 - j) * m + i];
    } else {
      col = hypso(e[q]);
      if (level === 2 && t.record?.ice[q]) col = lerp(col, [225, 235, 245], 0.55);
      if (level === 2 && t.record?.deposit[q] > 0.5) col = lerp(col, [200, 182, 130], 0.5);
      const thr = level === 1 ? 2e8 : 3e6;
      if (e[q] > 0 && A[q] > thr) col = lerp([90, 150, 215], [30, 80, 190], Math.log10(A[q] / thr) / 3);
    }
    const sh = Math.max(0.5, Math.min(1.4, 1 + ((h - hl) + (h - hu)) / (2 * t.dx) * shade));
    rgba[o] = col[0] * sh; rgba[o + 1] = col[1] * sh; rgba[o + 2] = col[2] * sh; rgba[o + 3] = 255;
  }
  return { w: m, h: m, rgba };
}

function describe(t, level, source, seconds) {
  const s = t.stats || {};
  const parts = [`level ${level}`, `${(t.dx >= 1000 ? (t.dx / 1000).toFixed(2) + ' km' : t.dx.toFixed(1) + ' m')} per node`, `${source}${seconds != null ? ` in ${seconds.toFixed(2)} s` : ''}`];
  if (s.consistency != null) parts.push(`block means within ${s.consistency.toFixed(2)} m of the parent`);
  if (level === 1 && s.hack) parts.push(`Hack ${s.hack.exponent.toFixed(2)}, concavity ${s.slopeArea.concavity.toFixed(2)}`);
  if (level === 3) parts.push(`${s.trees} trees, ${s.boulders} boulders, ${s.channels} channels`);
  return parts.join(' · ');
}

self.onmessage = async (ev) => {
  const msg = ev.data;
  try {
    if (msg.type === 'init') {
      await init(msg.seed, msg.store !== false);
      const t0 = performance.now();
      const img = renderPlanet();
      self.postMessage({ type: 'planet', ...img, info: `planet · ${world.planet().stats.cells} cells · built in ${((performance.now() - t0) / 1000).toFixed(1)} s · ${GENERATOR_VERSION} · level 1 on ${solverName}` }, [img.rgba.buffer]);
      return;
    }
    if (msg.type === 'pick') {
      // a click on the planet map: the level-1 tile there
      const p = world.planet();
      const i = pixelMap.map[msg.y * W0 + msg.x];
      if (i < 0) return;
      const P = p.mesh.pos;
      self.postMessage({ type: 'picked', key: tileKey(tileAt(1, P[3 * i], P[3 * i + 1], P[3 * i + 2])) });
      return;
    }
    if (msg.type === 'tile') {
      const addr = parseTileKey(msg.key);
      const hits0 = world.stats.memoryHits + world.stats.storeHits, gen0 = world.stats.generated;
      const t0 = performance.now();
      const t = await world.tile(addr);
      const seconds = (performance.now() - t0) / 1000;
      const source = world.stats.generated > gen0 ? `generated (${world.stats.generated - gen0} tiles)` : world.stats.memoryHits + world.stats.storeHits > hits0 ? 'from cache' : 'from cache';
      const img = renderTile(t, addr.level);
      self.postMessage({ type: 'tile', key: msg.key, seconds, source, ...img, info: describe(t, addr.level, source, seconds), n: t.n }, [img.rgba.buffer]);
      // prefetch the neighbours (and, below level 3, the child under the centre)
      const pre = ['W', 'E', 'S', 'N'].map((s) => neighborOf(addr, s));
      if (addr.level < 3) pre.push(childrenOf(addr)[Math.floor(WORLD_LAYOUT.branching / 2) * (WORLD_LAYOUT.branching + 1)]);
      if (msg.prefetch !== false) world.prefetch(pre);
      return;
    }
    if (msg.type === 'suggest') {
      // a level-1 tile over the highest temperate, wet land (for the gate's zoom path)
      const p = world.planet(), e = p.fields.elevation, T = p.fields.temperature, P = p.mesh.pos;
      let c = -1;
      for (let i = 0; i < e.length; i++) if (e[i] > 300 && e[i] < 2000 && T[i] > 4 && T[i] < 16 && p.fields.precip[i] > 700 && (c < 0 || e[i] > e[c])) c = i;
      if (c < 0) for (let i = 0; i < e.length; i++) if (c < 0 || e[i] > e[c]) c = i;
      self.postMessage({ type: 'suggested', key: tileKey(tileAt(1, P[3 * c], P[3 * c + 1], P[3 * c + 2])) });
      return;
    }
    if (msg.type === 'stats') { self.postMessage({ type: 'stats', stats: world.stats, queue: world.queue.length, busy: world.busy }); return; }
  } catch (err) {
    self.postMessage({ type: 'error', message: String(err && err.stack || err) });
  }
};
