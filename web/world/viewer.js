// The zooming viewer: planet -> region -> landscape -> site, generated in a worker.
import { parseTileKey, tileKey, parentOf, neighborOf, WORLD_LAYOUT } from '../../src/world/cubesphere.js';

const params = new URLSearchParams(location.search);
const seed = params.get('seed') || 'example';
const worker = new Worker('./worker.js', { type: 'module' });
const canvas = document.getElementById('view'), ctx = canvas.getContext('2d');
const info = document.getElementById('info'), status = document.getElementById('status'), crumbs = document.getElementById('crumbs');
let current = null; // null = planet, else a tile key
let planetImage = null;
const waiting = new Map();

function show(img) {
  canvas.width = img.w; canvas.height = img.h;
  ctx.putImageData(new ImageData(new Uint8ClampedArray(img.rgba.buffer || img.rgba), img.w, img.h), 0, 0);
}

function drawCrumbs() {
  crumbs.textContent = '';
  const chain = [];
  for (let k = current; k; ) { chain.unshift(k); const a = parseTileKey(k); k = a.level > 1 ? tileKey(parentOf(a)) : null; }
  const add = (label, key) => {
    const b = document.createElement('button');
    b.className = 'crumb' + (key === current ? ' current' : '');
    b.textContent = label;
    b.onclick = () => go(key);
    crumbs.append(b);
  };
  add('Planet', null);
  for (const k of chain) add(`L${parseTileKey(k).level} ${k}`, k);
}

function request(msg) {
  return new Promise((resolve) => { waiting.set(msg.key ?? msg.type, resolve); worker.postMessage(msg); });
}

async function go(key) {
  current = key;
  drawCrumbs();
  if (key === null) { show(planetImage.img); info.textContent = planetImage.info; status.textContent = ''; return planetImage; }
  status.textContent = `building ${key}…`;
  const r = await request({ type: 'tile', key });
  if (current !== key) return r;
  show(r);
  info.textContent = r.info;
  status.textContent = '';
  return r;
}

worker.onmessage = (ev) => {
  const m = ev.data;
  if (m.type === 'planet') { planetImage = { img: m, info: m.info }; if (current === null) { show(m); info.textContent = m.info; status.textContent = ''; } waiting.get('init')?.(m); return; }
  if (m.type === 'picked') { waiting.get('pick')?.(m); return; }
  if (m.type === 'tile') { waiting.get(m.key)?.(m); waiting.delete(m.key); return; }
  if (m.type === 'stats') { waiting.get('stats')?.(m); return; }
  if (m.type === 'suggested') { waiting.get('suggest')?.(m); return; }
  if (m.type === 'error') { status.textContent = m.message; window.gateResult = { error: m.message }; }
};

canvas.addEventListener('click', async (e) => {
  const r = canvas.getBoundingClientRect();
  const x = Math.floor(((e.clientX - r.left) / r.width) * canvas.width), y = Math.floor(((e.clientY - r.top) / r.height) * canvas.height);
  if (current === null) {
    const p = await request({ type: 'pick', x, y });
    go(p.key);
    return;
  }
  const a = parseTileKey(current);
  if (a.level >= 3) return;
  const n = canvas.width - 1, b = WORLD_LAYOUT.branching;
  const ci = Math.min(b - 1, Math.floor((x / n) * b)), cj = Math.min(b - 1, Math.floor(((n - y) / n) * b));
  go(tileKey({ level: a.level + 1, face: a.face, i: a.i * b + ci, j: a.j * b + cj }));
});
document.getElementById('up').onclick = () => { if (!current) return; const a = parseTileKey(current); go(a.level > 1 ? tileKey(parentOf(a)) : null); };
for (const btn of document.querySelectorAll('[data-side]')) btn.onclick = () => { if (current) go(tileKey(neighborOf(parseTileKey(current), btn.dataset.side))); };
document.addEventListener('keydown', (e) => {
  const side = { ArrowLeft: 'W', ArrowRight: 'E', ArrowUp: 'N', ArrowDown: 'S' }[e.key];
  if (side && current) { e.preventDefault(); go(tileKey(neighborOf(parseTileKey(current), side))); }
});

await request({ type: 'init', seed, store: params.get('store') !== '0', region: params.get('region') === 'mesh' ? 'mesh' : 'raster' });
drawCrumbs();

// ?auto=1: the browser gate's zoom path (planet -> a region -> landscape -> site), then
// revisits (cache) and a neighbour (prefetch)
if (params.get('auto')) {
  const out = { steps: [] };
  const step = async (label, key) => { const r = await go(key); out.steps.push({ label, key, seconds: r.seconds, source: r.source, info: r.info }); return r; };
  const start = params.get('start');
  const p = start ? { key: start } : await request({ type: 'suggest' });
  const k1 = p.key, a1 = parseTileKey(k1), b = WORLD_LAYOUT.branching;
  await step('region', k1);
  const k2 = tileKey({ level: 2, face: a1.face, i: a1.i * b + b / 2, j: a1.j * b + b / 2 });
  await step('landscape', k2);
  const a2 = parseTileKey(k2);
  const k3 = tileKey({ level: 3, face: a2.face, i: a2.i * b + b / 2, j: a2.j * b + b / 2 });
  await step('site', k3);
  await step('site again (cache)', k3);
  await step('landscape again (cache)', k2);
  // wait for the site's neighbours to be prefetched, then step east
  const deadline = performance.now() + 60000;
  for (;;) { const s = await request({ type: 'stats' }); if ((!s.queue && !s.busy) || performance.now() > deadline) break; await new Promise((r) => setTimeout(r, 200)); }
  await new Promise((r) => setTimeout(r, 1500));
  await step('site neighbour (prefetched)', tileKey(neighborOf(parseTileKey(k3), 'E')));
  window.gateResult = out;
}
