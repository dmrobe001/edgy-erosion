// Browser front end. All physics lives in ../src; this file only wires controls,
// the animation loop and drawing.

import { Simulation, PRESETS, VIEWS, resolveParams, cellColors, pixelMap, paint, crossSection } from '../src/index.js';

const $ = (id) => document.getElementById(id);

// Live-adjustable parameters. `log` sliders move in powers of ten.
const SLIDERS = [
  { group: 'Tectonics' },
  { path: 'uplift.rate', label: 'Uplift rate', unit: 'm/yr', min: 1e-6, max: 3e-3, log: true },
  { group: 'Climate' },
  { path: 'climate.precipitation', label: 'Precipitation', unit: 'm/yr', min: 0.05, max: 4, log: true },
  { path: 'climate.frost', label: 'Frost (freeze-thaw)', unit: '', min: 0, max: 1 },
  { group: 'Vegetation' },
  { path: 'vegetation.cover', label: 'Max cover (auto = from rain)', unit: '', min: 0, max: 1, auto: true },
  { group: 'Rivers' },
  { path: 'fluvial.deposition', label: 'Deposition G', unit: '', min: 0, max: 3 },
  { path: 'fluvial.soilErodibility', label: 'Soil erodibility', unit: '', min: 1e-6, max: 1e-2, log: true },
  { group: 'Hillslopes' },
  { path: 'hillslope.diffusivity', label: 'Creep diffusivity', unit: 'm²/yr', min: 1e-4, max: 0.5, log: true },
  { path: 'hillslope.criticalSlope', label: 'Soil repose slope', unit: '', min: 0.2, max: 1.5 },
  { path: 'weathering.enabled', label: 'Weathering (soil production)', toggle: true },
  { group: 'Rockfall' },
  { path: 'rockfall.enabled', label: 'Rockfall enabled', toggle: true },
  { path: 'rockfall.timescale', label: 'Failure timescale', unit: 'yr', min: 10, max: 1e5, log: true },
];

let sim;
let running = false;
let pm;
let mapImage;
let sectionRow = null;
const mapCanvas = $('map');
const mapCtx = mapCanvas.getContext('2d');
const secCanvas = $('section');
const secCtx = secCanvas.getContext('2d');

const get = (obj, path) => path.split('.').reduce((o, k) => o?.[k], obj);
const nest = (path, value) => path.split('.').reduceRight((acc, k) => ({ [k]: acc }), value);
const fmt = (v) => {
  if (typeof v !== 'number') return String(v);
  if (v === 0) return '0';
  const a = Math.abs(v);
  return a >= 1000 || a < 0.01 ? v.toExponential(1) : +v.toPrecision(3) + '';
};

// ------------------------------------------------------------------ setup

function buildPresetMenu() {
  const sel = $('preset');
  for (const [name, p] of Object.entries(PRESETS)) {
    const o = document.createElement('option');
    o.value = name;
    o.textContent = name;
    o.title = p.description;
    sel.appendChild(o);
  }
  const custom = document.createElement('option');
  custom.value = '';
  custom.textContent = '(custom)';
  sel.appendChild(custom);
  sel.addEventListener('change', () => {
    if (!sel.value) return;
    const params = resolveParams({ preset: sel.value });
    fillGridInputs(params);
    restart(params);
  });
}

function buildViewMenu() {
  const sel = $('view');
  for (const [k, label] of Object.entries(VIEWS)) {
    const o = document.createElement('option');
    o.value = k;
    o.textContent = label;
    sel.appendChild(o);
  }
  sel.addEventListener('change', draw);
}

function buildSliders() {
  const host = $('sliders');
  host.innerHTML = '';
  for (const s of SLIDERS) {
    if (s.group) {
      const h = document.createElement('h2');
      h.textContent = s.group;
      h.style.marginTop = '0.7rem';
      host.appendChild(h);
      continue;
    }
    if (s.toggle) {
      const lab = document.createElement('label');
      lab.className = 'check';
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.dataset.path = s.path;
      cb.addEventListener('change', () => sim.setParams(nest(s.path, cb.checked)));
      lab.append(cb, document.createTextNode(s.label));
      host.appendChild(lab);
      continue;
    }
    const wrap = document.createElement('div');
    wrap.className = 'slider';
    const head = document.createElement('div');
    head.className = 'head';
    const name = document.createElement('span');
    name.textContent = s.label;
    const val = document.createElement('span');
    val.className = 'val';
    head.append(name, val);
    const input = document.createElement('input');
    input.type = 'range';
    input.min = 0;
    input.max = 1000;
    input.dataset.path = s.path;
    const toValue = (t) => (s.log ? s.min * Math.pow(s.max / s.min, t / 1000) : s.min + ((s.max - s.min) * t) / 1000);
    const toSlider = (v) => (s.log ? (1000 * Math.log(v / s.min)) / Math.log(s.max / s.min) : (1000 * (v - s.min)) / (s.max - s.min));
    input.addEventListener('input', () => {
      const v = toValue(+input.value);
      val.textContent = `${fmt(v)} ${s.unit}`;
      sim.setParams(nest(s.path, v));
      if (autoBox) autoBox.checked = false;
    });
    let autoBox = null;
    wrap._sync = () => {
      const v = get(sim.params, s.path);
      if (s.auto && v === 'auto') {
        autoBox.checked = true;
        val.textContent = 'auto';
      } else {
        input.value = toSlider(v);
        val.textContent = `${fmt(v)} ${s.unit}`;
        if (autoBox) autoBox.checked = false;
      }
    };
    wrap.append(head, input);
    if (s.auto) {
      const lab = document.createElement('label');
      lab.className = 'check';
      autoBox = document.createElement('input');
      autoBox.type = 'checkbox';
      autoBox.addEventListener('change', () => {
        sim.setParams(nest(s.path, autoBox.checked ? 'auto' : toValue(+input.value)));
        wrap._sync();
      });
      lab.append(autoBox, document.createTextNode('auto'));
      wrap.appendChild(lab);
    }
    host.appendChild(wrap);
  }
}

function syncControls() {
  for (const el of $('sliders').children) if (el._sync) el._sync();
  for (const cb of $('sliders').querySelectorAll('input[type=checkbox][data-path]')) cb.checked = !!get(sim.params, cb.dataset.path);
  $('dt').value = sim.params.time.dt;
  $('preset').value = sim.params.preset || '';
  $('preset-desc').textContent = sim.params.preset ? PRESETS[sim.params.preset]?.description ?? '' : 'Custom parameters.';
  refreshJson();
  updateExtent();
}

function fillGridInputs(params) {
  $('cols').value = params.grid.cols;
  $('rows').value = params.grid.rows;
  $('dx').value = params.grid.dx;
  $('seed').value = params.seed;
  $('dt').value = params.time.dt;
}

function updateExtent() {
  const cols = +$('cols').value, rows = +$('rows').value, dx = +$('dx').value;
  const w = dx * (cols + 0.5) / 1000, h = (dx * Math.sqrt(3) / 2 * (rows - 1) + 2 * dx / Math.sqrt(3)) / 1000;
  $('extent').textContent = `${w.toFixed(2)} × ${h.toFixed(2)} km, ${cols * rows} cells`;
}

function gridOverrides() {
  return {
    seed: +$('seed').value,
    grid: { cols: +$('cols').value, rows: +$('rows').value, dx: +$('dx').value },
    time: { dt: +$('dt').value },
  };
}

function restart(params) {
  try {
    sim = new Simulation(params);
  } catch (e) {
    $('json-error').textContent = e.message;
    return;
  }
  $('json-error').textContent = '';
  sectionRow = sim.grid.rows >> 1;
  resize();
  syncControls();
  draw();
}

// ------------------------------------------------------------------ drawing

function resize() {
  const wrap = $('map-wrap');
  const cssW = Math.max(200, wrap.clientWidth);
  const g = sim.grid;
  const aspect = g.height / g.width;
  const maxH = window.innerHeight * 0.62;
  let w = cssW;
  let h = w * aspect;
  if (h > maxH) { h = maxH; w = h / aspect; }
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  // never draw at more than ~3 device pixels per cell row: costs time, adds nothing
  const W = Math.round(Math.min(w * dpr, g.cols * 12));
  const H = Math.round(W * aspect);
  mapCanvas.width = W;
  mapCanvas.height = H;
  mapCanvas.style.width = `${w}px`;
  mapCanvas.style.height = `${w * aspect}px`;
  pm = pixelMap(g, W, H);
  mapImage = mapCtx.createImageData(W, H);
  secCanvas.width = W;
  secCanvas.height = Math.round(W / 3.2);
  secCanvas.style.width = `${w}px`;
  secCanvas.style.height = `${w / 3.2}px`;
}

function draw() {
  const view = $('view').value;
  paint(pm, cellColors(sim, view), mapImage.data);
  mapCtx.putImageData(mapImage, 0, 0);
  // cross-section marker
  const y = (sim.grid.y[sectionRow * sim.grid.cols] - pm.y0) * pm.scale;
  mapCtx.strokeStyle = 'rgba(255,255,255,0.6)';
  mapCtx.setLineDash([6, 4]);
  mapCtx.lineWidth = 1;
  mapCtx.beginPath();
  mapCtx.moveTo(0, y + 0.5);
  mapCtx.lineTo(mapCanvas.width, y + 0.5);
  mapCtx.stroke();
  mapCtx.setLineDash([]);

  const img = secCtx.createImageData(secCanvas.width, secCanvas.height);
  const info = crossSection(sim, sectionRow, secCanvas.width, secCanvas.height, img.data);
  secCtx.putImageData(img, 0, 0);
  $('section-info').textContent = `row ${sectionRow}  z ${info.zmin.toFixed(0)}…${info.zmax.toFixed(0)} m  vertical exaggeration ${info.exaggeration.toFixed(1)}×`;

  const s = sim.stats();
  $('stats').textContent =
    `t = ${(s.time / 1000).toFixed(1)} kyr  step ${s.step}  z ${s.minElevation.toFixed(0)}…${s.maxElevation.toFixed(0)} m  ` +
    `soil ${s.meanSoil.toFixed(2)} m  rockfalls ${s.rockfallEvents}`;
  drawLegend(view);
}

function drawLegend(view) {
  const host = $('legend');
  if (view !== 'geology' && view !== 'terrain') { host.innerHTML = ''; return; }
  const used = new Set();
  for (let i = 0; i < sim.grid.size; i++) used.add(sim.surfaceRock[i]);
  for (const k of sim.lith.layerRock) used.add(k);
  host.innerHTML = [...used].sort().map((k) => `<span style="--c:${sim.lith.color[k]}">${sim.lith.names[k]}</span>`).join('') +
    (view === 'terrain' ? '<span style="--c:#b09870">soil / talus</span><span style="--c:#487034">vegetation</span><span style="--c:#4678af">water</span>' : '');
}

function loop() {
  if (running) {
    const t0 = performance.now();
    do { sim.step(); } while (performance.now() - t0 < 30);
    draw();
  }
  requestAnimationFrame(loop);
}

// ------------------------------------------------------------------ events

function cellAtEvent(ev) {
  const r = mapCanvas.getBoundingClientRect();
  const px = Math.floor(((ev.clientX - r.left) / r.width) * mapCanvas.width);
  const py = Math.floor(((ev.clientY - r.top) / r.height) * mapCanvas.height);
  if (px < 0 || py < 0 || px >= pm.width || py >= pm.height) return -1;
  return pm.map[py * pm.width + px];
}

mapCanvas.addEventListener('mousemove', (ev) => {
  const i = cellAtEvent(ev);
  if (i < 0) { $('hover').innerHTML = '&nbsp;'; return; }
  const g = sim.grid;
  const lake = sim.flow.filled[i] - sim.z[i];
  $('hover').textContent =
    `cell (${g.col(i)}, ${g.row(i)})  z ${sim.z[i].toFixed(1)} m  soil ${sim.soil[i].toFixed(2)} m  ` +
    `rock ${sim.lith.names[sim.surfaceRock[i]]}  drainage ${(sim.flow.area[i] / 1e6).toFixed(3)} km²  ` +
    `vegetation ${(100 * sim.veg[i]).toFixed(0)}%  dz/dt ${(sim.dzdt[i] * 1000).toFixed(3)} mm/yr` +
    (lake > 0.05 ? `  lake ${lake.toFixed(1)} m deep` : '');
});
mapCanvas.addEventListener('click', (ev) => {
  const i = cellAtEvent(ev);
  if (i < 0) return;
  sectionRow = sim.grid.row(i);
  draw();
});

$('run').addEventListener('click', () => {
  running = !running;
  $('run').textContent = running ? 'Pause' : 'Run';
});
$('step').addEventListener('click', () => { sim.step(); draw(); });
$('reset').addEventListener('click', () => {
  const base = { ...sim.params };
  restart(resolveParams({ ...base, ...gridOverrides(), grid: gridOverrides().grid, time: { ...base.time, dt: +$('dt').value } }));
});
$('dt').addEventListener('change', () => sim.setParams({ time: { dt: +$('dt').value } }));
for (const id of ['cols', 'rows', 'dx']) $(id).addEventListener('input', updateExtent);

function refreshJson() {
  $('json').value = JSON.stringify(sim.params, null, 2);
}
$('json-refresh').addEventListener('click', refreshJson);
$('json-apply').addEventListener('click', () => {
  try {
    const p = resolveParams(JSON.parse($('json').value));
    fillGridInputs(p);
    restart(p);
  } catch (e) {
    $('json-error').textContent = e.message;
  }
});

function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
$('save-png').addEventListener('click', () => {
  mapCanvas.toBlob((b) => download(`edgy-${$('view').value}-t${Math.round(sim.time)}.png`, b));
});
$('save-snapshot').addEventListener('click', () => {
  download(`edgy-snapshot-t${Math.round(sim.time)}.json`, new Blob([JSON.stringify(sim.snapshot())], { type: 'application/json' }));
});
$('load').addEventListener('change', async (ev) => {
  const file = ev.target.files[0];
  if (!file) return;
  try {
    const data = JSON.parse(await file.text());
    if (data.format?.startsWith('edgy-erosion-snapshot')) {
      sim = Simulation.fromSnapshot(data);
      fillGridInputs(sim.params);
      sectionRow = sim.grid.rows >> 1;
      resize();
      syncControls();
      draw();
    } else {
      const p = resolveParams(data);
      fillGridInputs(p);
      restart(p);
    }
  } catch (e) {
    $('json-error').textContent = e.message;
  }
});

let resizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { resize(); draw(); }, 150);
});

// ------------------------------------------------------------------ start

buildPresetMenu();
buildViewMenu();
buildSliders();
const initial = resolveParams({ preset: 'mesa-country' });
fillGridInputs(initial);
restart(initial);
requestAnimationFrame(loop);

// handy for poking at the simulation from the browser console
window.edgy = { get sim() { return sim; }, Simulation, resolveParams };
