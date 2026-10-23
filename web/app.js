// Browser front end. All physics lives in ../src; this file only wires controls,
// the animation loop and drawing.
//
// `sim` is the level being shown: a landscape Simulation, or a Planet. Zooming in
// stacks the frozen levels above it.

import {
  Simulation, PRESETS, VIEWS, resolveParams, cellColors, pixelMap, paint, crossSection, zoomIn, compareWithParent,
  Planet, PLANET_PRESETS, PLANET_VIEWS, resolvePlanetParams, planetColors, spherePixelMap, planetSection, plateColor,
  zoomFromPlanet, azimuthalProjection, lonLatToVector,
} from '../src/index.js';

const $ = (id) => document.getElementById(id);

// Live-adjustable parameters. `log` sliders move in powers of ten.
const SLIDERS = [
  { group: 'Tectonics' },
  { path: 'uplift.rate', label: 'Uplift rate', unit: 'm/yr', min: 1e-6, max: 3e-3, log: true },
  { group: 'Earthquakes (on the faults in the parameters)' },
  { path: 'earthquakes.enabled', label: 'Faults slip in earthquakes', toggle: true },
  { path: 'earthquakes.maxMagnitude', label: 'Largest magnitude', unit: '', min: 5.6, max: 8.5 },
  { path: 'earthquakes.bValue', label: 'Gutenberg-Richter b', unit: '', min: 0.5, max: 1.5 },
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
  { group: 'Coast (when there is a sea)' },
  { path: 'coast.waveHeight', label: 'Wave height', unit: 'm', min: 0.2, max: 8, log: true },
  { path: 'coast.seaLevelAmplitude', label: 'Sea-level swing', unit: 'm', min: 0, max: 120 },
  { path: 'coast.seaLevelPeriod', label: 'Sea-level period', unit: 'yr', min: 5000, max: 200000, log: true },
  { path: 'coast.marineDiffusivity', label: 'Shore sediment mobility', unit: 'm²/yr', min: 0.1, max: 100, log: true },
  { group: 'Rockfall' },
  { path: 'rockfall.enabled', label: 'Rockfall enabled', toggle: true },
  { path: 'rockfall.timescale', label: 'Failure timescale', unit: 'yr', min: 10, max: 1e5, log: true },
];

const PLANET_SLIDERS = [
  { group: 'Plates' },
  { path: 'plates.collisionLocking', label: 'Collision locking time', unit: 'yr', min: 1e6, max: 2e8, log: true },
  { path: 'plates.reorganizeEvery', label: 'Reorganise plates every (0 = never)', unit: 'yr', min: 0, max: 2e8 },
  { group: 'Crust' },
  { path: 'crust.arcFraction', label: 'Arc growth per subducted crust', unit: '', min: 0, max: 0.5 },
  { path: 'crust.flowDiffusivity', label: 'Thick-crust flow', unit: 'm²/yr', min: 1e3, max: 1e6, log: true },
  { path: 'crust.trenchDepth', label: 'Trench depth', unit: 'm', min: 0, max: 5000 },
  { group: 'Climate' },
  { path: 'climate.precipitation', label: 'Precipitation (scale)', unit: 'm/yr', min: 0.1, max: 4, log: true },
  { path: 'climate.zonal', label: 'Zonal rain belts (wet tropics, dry subtropics)', toggle: true },
  { group: 'Erosion (effective values for planet-scale cells)' },
  { path: 'surface.rocks.crust.erodibility', label: 'Rock erodibility', unit: '', min: 1e-8, max: 3e-6, log: true },
  { path: 'surface.hillslope.diffusivity', label: 'Creep diffusivity', unit: 'm²/yr', min: 10, max: 3e4, log: true },
];

let sim;
let running = false;
let pm;
let mapImage;
let sectionRow = null;
let sectionLat = 20;
// Zoom: the frozen ancestors of the current simulation, innermost last, each with the
// options that made its child (so Reset can re-make the child from scratch).
let stack = [];
let zoomPending = false;
let lastPatch = null; // { parent, bounds | outline } of the most recently left child, drawn on its parent
let globe = { lon: 0, lat: 15 };
const mapCanvas = $('map');
const mapCtx = mapCanvas.getContext('2d');
const secCanvas = $('section');
const secCtx = secCanvas.getContext('2d');

const isPlanet = (s) => s instanceof Planet;
const get = (obj, path) => path.split('.').reduce((o, k) => o?.[k], obj);
const nest = (path, value) => path.split('.').reduceRight((acc, k) => ({ [k]: acc }), value);
const fmt = (v) => {
  if (typeof v !== 'number') return String(v);
  if (v === 0) return '0';
  const a = Math.abs(v);
  return a >= 1000 || a < 0.01 ? v.toExponential(1) : +v.toPrecision(3) + '';
};
const deg = (r) => (r * 180) / Math.PI;

// ------------------------------------------------------------------ setup

function buildPresetMenu() {
  const sel = $('preset');
  const group = (label, entries, prefix) => {
    const g = document.createElement('optgroup');
    g.label = label;
    for (const [name, p] of Object.entries(entries)) {
      const o = document.createElement('option');
      o.value = prefix + name;
      o.textContent = name;
      o.title = p.description;
      g.appendChild(o);
    }
    sel.appendChild(g);
  };
  group('Landscapes', PRESETS, '');
  group('Planets (plates, ~150-km cells)', PLANET_PRESETS, 'planet:');
  const custom = document.createElement('option');
  custom.value = '';
  custom.textContent = '(custom)';
  sel.appendChild(custom);
  sel.addEventListener('change', () => {
    if (!sel.value) return;
    stack = []; lastPatch = null;
    if (sel.value.startsWith('planet:')) {
      const params = resolvePlanetParams({ preset: sel.value.slice(7) });
      fillPlanetInputs(params);
      restart(params);
      return;
    }
    const params = resolveParams({ preset: sel.value });
    fillGridInputs(params);
    restart(params);
  });
}

function buildViewMenu() {
  const sel = $('view');
  const views = isPlanet(sim) ? PLANET_VIEWS : VIEWS;
  const keep = sel.value;
  sel.innerHTML = '';
  for (const [k, label] of Object.entries(views)) {
    const o = document.createElement('option');
    o.value = k;
    o.textContent = label;
    sel.appendChild(o);
  }
  sel.value = views[keep] ? keep : 'terrain';
}
$('view').addEventListener('change', draw);

function buildSliders() {
  const host = $('sliders');
  host.innerHTML = '';
  for (const s of isPlanet(sim) ? PLANET_SLIDERS : SLIDERS) {
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

function buildZoomFactors() {
  const sel = $('zoom-factor');
  const opts = isPlanet(sim) ? [8, 16, 32, 64] : [2, 4, 8];
  const def = isPlanet(sim) ? 32 : 4;
  sel.innerHTML = opts.map((f) => `<option${f === def ? ' selected' : ''}>${f}</option>`).join('');
}

function syncControls() {
  const planet = isPlanet(sim);
  $('planet-grid').hidden = !planet;
  $('landscape-grid').hidden = planet;
  $('projection-wrap').hidden = !planet;
  $('save-snapshot').disabled = planet;
  $('section-hint').textContent = planet ? '(along a parallel; click the map to choose the latitude)' : '(click the map to choose a row)';
  $('quake').hidden = planet || !sim.params.earthquakes?.faults?.length;
  buildViewMenu();
  buildSliders();
  buildZoomFactors();
  for (const el of $('sliders').children) if (el._sync) el._sync();
  for (const cb of $('sliders').querySelectorAll('input[type=checkbox][data-path]')) cb.checked = !!get(sim.params, cb.dataset.path);
  if (planet) {
    $('planet-dt').value = sim.params.time.dt;
    $('preset').value = sim.params.preset ? `planet:${sim.params.preset}` : '';
    $('preset-desc').textContent = sim.params.preset ? PLANET_PRESETS[sim.params.preset]?.description ?? '' : 'Custom planet.';
  } else {
    $('dt').value = sim.params.time.dt;
    $('preset').value = sim.params.preset || '';
    $('preset-desc').textContent = sim.params.preset ? PRESETS[sim.params.preset]?.description ?? '' : 'Custom parameters.';
  }
  refreshJson();
  updateExtent();
  syncZoom();
}

function levelName(s) {
  if (isPlanet(s)) return `planet: ~${(s.mesh.spacing / 1000).toFixed(0)} km cells`;
  const dx = s.grid.dx;
  return dx >= 1000 ? `dx ${+(dx / 1000).toPrecision(3)} km` : `dx ${+dx.toPrecision(3)} m`;
}

function syncZoom() {
  const levels = [...stack.map((f) => f.sim), sim];
  $('zoom-path').textContent = levels.map((s, k) => `level ${k}: ${levelName(s)}`).join(' › ') +
    (stack.length ? '' : '. Zoom in to refine a patch; the level above freezes and feeds its edges and rivers.');
  $('zoom-out').disabled = !stack.length;
  for (const id of ['cols', 'rows', 'dx', 'seed', 'preset', 'subdivisions', 'plates', 'planet-seed']) $(id).disabled = stack.length > 0;
  $('zoom-in').textContent = zoomPending ? 'Click the map…' : 'Zoom in…';
}

function setSim(next) {
  sim = next;
  $('hover').innerHTML = '&nbsp;';
  if (isPlanet(sim)) fillPlanetInputs(sim.params); else fillGridInputs(sim.params);
  sectionRow = isPlanet(sim) ? null : sim.grid.rows >> 1;
  syncControls();
  resize();
  draw();
}

function fillGridInputs(params) {
  $('cols').value = params.grid.cols;
  $('rows').value = params.grid.rows;
  $('dx').value = params.grid.dx;
  $('seed').value = params.seed;
  $('dt').value = params.time.dt;
}

function fillPlanetInputs(params) {
  $('subdivisions').value = params.subdivisions;
  $('plates').value = params.plates.count;
  $('planet-seed').value = params.seed;
  $('planet-dt').value = params.time.dt;
}

function updateExtent() {
  const cols = +$('cols').value, rows = +$('rows').value, dx = +$('dx').value;
  const w = dx * (cols + 0.5) / 1000, h = (dx * Math.sqrt(3) / 2 * (rows - 1) + 2 * dx / Math.sqrt(3)) / 1000;
  $('extent').textContent = `${w.toFixed(2)} × ${h.toFixed(2)} km, ${cols * rows} cells`;
  const f = +$('subdivisions').value;
  const cells = 10 * f * f + 2;
  $('planet-extent').textContent = `${cells} cells, ~${((6371 * Math.sqrt((4 * Math.PI) / cells)) * 1.07).toFixed(0)} km apart`;
}

function gridOverrides() {
  return {
    seed: +$('seed').value,
    grid: { cols: +$('cols').value, rows: +$('rows').value, dx: +$('dx').value },
    time: { dt: +$('dt').value },
  };
}

function planetOverrides() {
  return {
    seed: +$('planet-seed').value,
    subdivisions: +$('subdivisions').value,
    plates: { count: +$('plates').value },
    time: { dt: +$('planet-dt').value },
  };
}

function restart(params) {
  try {
    sim = params.subdivisions != null ? new Planet(params) : new Simulation(params);
  } catch (e) {
    $('json-error').textContent = e.message;
    return;
  }
  $('json-error').textContent = '';
  sectionRow = isPlanet(sim) ? null : sim.grid.rows >> 1;
  syncControls();
  resize();
  draw();
}

// ------------------------------------------------------------------ drawing

function resize() {
  const wrap = $('map-wrap');
  const cssW = Math.max(200, wrap.clientWidth);
  const planet = isPlanet(sim);
  const globeView = planet && $('projection').value === 'globe';
  const aspect = planet ? (globeView ? 1 : 0.5) : sim.grid.height / sim.grid.width;
  const maxH = window.innerHeight * 0.62;
  let w = cssW;
  let h = w * aspect;
  if (h > maxH) { h = maxH; w = h / aspect; }
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  // never draw at more than ~3 device pixels per cell row: costs time, adds nothing
  const W = Math.round(Math.min(w * dpr, planet ? 1400 : sim.grid.cols * 12));
  const H = Math.round(W * aspect);
  mapCanvas.width = W;
  mapCanvas.height = H;
  mapCanvas.style.width = `${w}px`;
  mapCanvas.style.height = `${w * aspect}px`;
  pm = planet ? spherePixelMap(sim.mesh, W, H, { projection: globeView ? 'globe' : 'map', lon0: globe.lon, lat0: globe.lat }) : pixelMap(sim.grid, W, H);
  mapImage = mapCtx.createImageData(W, H);
  secCanvas.width = W;
  secCanvas.height = Math.round(W / 3.2);
  secCanvas.style.width = `${w}px`;
  secCanvas.style.height = `${w / 3.2}px`;
}

/** Pixel position of a longitude/latitude (degrees) on the current planet map, or null. */
function lonLatToPixel(lon, lat) {
  if (pm.projection === 'globe') {
    const v = lonLatToVector(lon, lat);
    const c = lonLatToVector(pm.lon0, pm.lat0);
    if (v[0] * c[0] + v[1] * c[1] + v[2] * c[2] < 0) return null;
    const l0 = (pm.lon0 * Math.PI) / 180;
    const e = [-Math.sin(l0), Math.cos(l0), 0];
    const u = [c[1] * e[2] - c[2] * e[1], c[2] * e[0] - c[0] * e[2], c[0] * e[1] - c[1] * e[0]];
    const a = v[0] * e[0] + v[1] * e[1] + v[2] * e[2];
    const b = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
    return [pm.width / 2 + a * pm.radiusPx, pm.height / 2 - b * pm.radiusPx];
  }
  return [((lon + 180) / 360) * pm.width, ((90 - lat) / 180) * pm.height];
}

function draw() {
  const view = $('view').value;
  const planet = isPlanet(sim);
  paint(pm, planet ? planetColors(sim, view) : cellColors(sim, view), mapImage.data, planet ? [0, 0, 0] : undefined);
  mapCtx.putImageData(mapImage, 0, 0);
  mapCtx.setLineDash([6, 4]);
  mapCtx.lineWidth = 1;
  mapCtx.strokeStyle = 'rgba(255,255,255,0.6)';
  if (planet) {
    if (pm.projection !== 'globe') { // the section's parallel
      const y = ((90 - sectionLat) / 180) * pm.height;
      mapCtx.beginPath(); mapCtx.moveTo(0, y + 0.5); mapCtx.lineTo(mapCanvas.width, y + 0.5); mapCtx.stroke();
    }
  } else {
    const y = (sim.grid.y[sectionRow * sim.grid.cols] - pm.y0) * pm.scale;
    mapCtx.beginPath(); mapCtx.moveTo(0, y + 0.5); mapCtx.lineTo(mapCanvas.width, y + 0.5); mapCtx.stroke();
  }
  mapCtx.setLineDash([]);
  if (!planet) drawFaults();
  if (lastPatch && lastPatch.parent === sim) {
    mapCtx.strokeStyle = 'rgba(255,255,255,0.9)';
    mapCtx.lineWidth = 2;
    if (lastPatch.outline) {
      mapCtx.beginPath();
      let pen = false;
      for (const [lon, lat] of lastPatch.outline) {
        const p = lonLatToPixel(lon, lat);
        if (!p) { pen = false; continue; }
        if (pen && pm.projection !== 'globe' && Math.abs(p[0] - pen[0]) > pm.width / 2) pen = false; // date line
        if (pen) mapCtx.lineTo(p[0], p[1]); else mapCtx.moveTo(p[0], p[1]);
        pen = p;
      }
      mapCtx.stroke();
    } else {
      const b = lastPatch.bounds;
      mapCtx.strokeRect((b.xmin - pm.x0) * pm.scale, (b.ymin - pm.y0) * pm.scale, (b.xmax - b.xmin) * pm.scale, (b.ymax - b.ymin) * pm.scale);
    }
  }

  const img = secCtx.createImageData(secCanvas.width, secCanvas.height);
  if (planet) {
    const info = planetSection(sim, sectionLat, secCanvas.width, secCanvas.height, img.data);
    $('section-info').textContent = `latitude ${sectionLat.toFixed(1)}°  ${(info.zmin / 1000).toFixed(0)}…${(info.zmax / 1000).toFixed(0)} km  vertical exaggeration ${info.exaggeration.toFixed(0)}×  (tan: continental crust, grey: oceanic, red: mantle)`;
  } else {
    const info = crossSection(sim, sectionRow, secCanvas.width, secCanvas.height, img.data);
    $('section-info').textContent = `row ${sectionRow}  z ${info.zmin.toFixed(0)}…${info.zmax.toFixed(0)} m  vertical exaggeration ${info.exaggeration.toFixed(1)}×`;
  }
  secCtx.putImageData(img, 0, 0);

  if (planet) {
    const s = sim.stats();
    $('stats').textContent =
      `t = ${(s.time / 1e6).toFixed(1)} Myr  step ${s.step}  plates ${s.plates}  land ${(100 * s.landFraction).toFixed(0)}%  ` +
      `z ${(s.minElevation / 1000).toFixed(1)}…${(s.maxElevation / 1000).toFixed(1)} km  sea level ${s.seaLevel.toFixed(0)} m`;
  } else {
    const s = sim.stats();
    const q = sim.quakes;
    const lastQ = q.events.length ? q.events[q.events.length - 1] : null;
    $('stats').textContent =
      `t = ${(s.time / 1000).toFixed(1)} kyr  step ${s.step}  z ${s.minElevation.toFixed(0)}…${s.maxElevation.toFixed(0)} m  ` +
      `soil ${s.meanSoil.toFixed(2)} m  rockfalls ${s.rockfallEvents}` +
      (s.seaLevel != null ? `  sea level ${s.seaLevel.toFixed(1)} m  cliff collapses ${s.cliffCollapses}` : '') +
      (q.count ? `  earthquakes ${q.count} (last M${lastQ.magnitude.toFixed(1)}, ${((sim.time - lastQ.time) / 1000).toFixed(1)} kyr ago)` : '');
  }
  drawLegend(view);
  if (stack.length) {
    const parent = stack[stack.length - 1].sim;
    const c = compareWithParent(sim, parent);
    const km2 = (a) => (a / 1e6 >= 1e4 ? (a / 1e6).toExponential(2) : (a / 1e6).toFixed(a / 1e6 > 100 ? 0 : 3));
    $('zoom-compare').textContent =
      `Against the frozen level above, over ${c.parentCells} of its cells: mean elevation ${c.meanElevation.child.toFixed(1)} m ` +
      `(level above, extrapolated at its last rate: ${c.meanElevation.expected.toFixed(1)} m); largest river leaving ` +
      `${km2(c.largestExitArea.child)} km² (level above: ${km2(c.largestExitArea.parent)} km²).` + boundarySummary();
  } else {
    $('zoom-compare').textContent = '';
  }
}

/** The plate boundaries a planet patch inherited, and the faults made from them. */
function boundarySummary() {
  const bs = sim.zoom?.boundaries;
  if (!bs) return '';
  if (!bs.length) return ' No plate boundary crosses this patch.';
  return ' Plate boundaries here: ' + bs.map((b) => {
    const rate = `${(Math.abs(b.closingRate) * 100).toFixed(1)} cm/yr ${b.closingRate > 0 ? 'closing' : 'opening'}`;
    return b.fault
      ? `${b.kind} fault (plates ${b.plates.join('/')}, ${rate}, slip ${(b.fault.slipRate * 1000).toFixed(2)} mm/yr, plate ${b.raised} side raised)`
      : `transform (plates ${b.plates.join('/')}, sliding ${(b.relativeSpeed * 100).toFixed(1)} cm/yr; no fault: only vertical slip is modelled)`;
  }).join('; ') + '.';
}

/** Fault traces (dashed), and this step's ruptures (bright), on a landscape map. */
function drawFaults() {
  const faults = sim.params.earthquakes?.faults || [];
  if (!faults.length) return;
  const F = sim.frame;
  const px = (fx, fy) => [(F.x + fx * F.width - pm.x0) * pm.scale, (F.y + fy * F.height - pm.y0) * pm.scale];
  mapCtx.save();
  mapCtx.beginPath();
  mapCtx.rect(0, 0, mapCanvas.width, mapCanvas.height);
  mapCtx.clip();
  mapCtx.setLineDash([3, 4]);
  mapCtx.strokeStyle = 'rgba(255,90,70,0.8)';
  mapCtx.lineWidth = 1.5;
  for (const f of faults) {
    const a = px(f.x0, f.y0), b = px(f.x1, f.y1);
    mapCtx.beginPath(); mapCtx.moveTo(a[0], a[1]); mapCtx.lineTo(b[0], b[1]); mapCtx.stroke();
  }
  mapCtx.setLineDash([]);
  const recent = sim.quakes.events.slice(-sim.quakes.lastStep || sim.quakes.events.length);
  if (sim.quakes.lastStep) {
    mapCtx.strokeStyle = 'rgba(255,240,120,0.95)';
    mapCtx.lineWidth = 4;
    for (const ev of recent) {
      const f = faults[ev.fault];
      const dxw = (f.x1 - f.x0) * F.width, dyw = (f.y1 - f.y0) * F.height;
      const L = Math.hypot(dxw, dyw);
      const ux = dxw / L, uy = dyw / L;
      const a = [(ev.center[0] - (ux * ev.length) / 2 - pm.x0) * pm.scale, (ev.center[1] - (uy * ev.length) / 2 - pm.y0) * pm.scale];
      const b = [(ev.center[0] + (ux * ev.length) / 2 - pm.x0) * pm.scale, (ev.center[1] + (uy * ev.length) / 2 - pm.y0) * pm.scale];
      mapCtx.beginPath(); mapCtx.moveTo(a[0], a[1]); mapCtx.lineTo(b[0], b[1]); mapCtx.stroke();
    }
  }
  mapCtx.restore();
}

function drawLegend(view) {
  const host = $('legend');
  if (isPlanet(sim)) {
    const items = {
      plates: '<span style="--c:#ff463c">convergent</span><span style="--c:#ffeb5a">divergent (ridge)</span><span style="--c:#ebebeb">transform</span><span>plate colours: lighter = continental</span>',
      age: '<span style="--c:#c82828">0 Myr</span><span style="--c:#f0dc46">60</span><span style="--c:#5abe6e">100</span><span style="--c:#3c78c8">150</span><span style="--c:#46288c">200+</span><span style="--c:#999">continent</span>',
      crust: '<span style="--c:#141e46">0 km</span><span style="--c:#78aa78">28</span><span style="--c:#dcc86e">40</span><span style="--c:#c86432">60</span><span style="--c:#6e141e">80+</span>',
    };
    host.innerHTML = items[view] || '';
    if (view === 'plates') {
      const alive = [...new Set(sim.owner)].sort((a, b) => a - b);
      host.innerHTML += alive.map((q) => `<span style="--c:rgb(${plateColor(q).map(Math.round).join(',')})">plate ${q}</span>`).join('');
    }
    return;
  }
  if (view !== 'geology' && view !== 'terrain') { host.innerHTML = ''; return; }
  const used = new Set();
  for (let i = 0; i < sim.grid.size; i++) used.add(sim.surfaceRock[i]);
  for (const k of sim.lith.layerRock) used.add(k);
  host.innerHTML = [...used].sort().map((k) => `<span style="--c:${sim.lith.color[k]}">${sim.lith.names[k]}</span>`).join('') +
    (view === 'terrain' ? '<span style="--c:#b09870">soil / talus</span><span style="--c:#decc96">beach sand</span><span style="--c:#487034">vegetation</span><span style="--c:#4678af">water</span>' : '');
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

const BOUNDARY = ['', 'convergent', 'divergent', 'transform'];

mapCanvas.addEventListener('mousemove', (ev) => {
  if (drag) return;
  const i = cellAtEvent(ev);
  if (i < 0) { $('hover').innerHTML = '&nbsp;'; return; }
  if (isPlanet(sim)) {
    const m = sim.mesh, s = sim.sim;
    const wet = s.z[i] < sim.seaLevel;
    $('hover').textContent =
      `${deg(m.lat[i]).toFixed(1)}°, ${deg(m.lon[i]).toFixed(1)}°  z ${s.z[i].toFixed(0)} m${wet ? ` (${(sim.seaLevel - s.z[i]).toFixed(0)} m deep)` : ''}  ` +
      `crust ${(sim.thickness[i] / 1000).toFixed(1)} km, ${(100 * sim.continental[i]).toFixed(0)}% continental` +
      (sim.continental[i] < 0.5 ? `, ${(sim.age[i] / 1e6).toFixed(0)} Myr old` : '') +
      `  plate ${sim.owner[i]}${sim.boundary[i] ? ` (${BOUNDARY[sim.boundary[i]]} boundary)` : ''}  ` +
      `uplift ${(sim.upliftRate[i] * 1000).toFixed(3)} mm/yr  rain ${s.precipWeight[i].toFixed(2)} m/yr`;
    return;
  }
  const g = sim.grid;
  const lake = sim.flow.filled[i] - sim.z[i];
  $('hover').textContent =
    `cell (${g.col(i)}, ${g.row(i)})  z ${sim.z[i].toFixed(1)} m  soil ${sim.soil[i].toFixed(2)} m  ` +
    `rock ${sim.lith.names[sim.surfaceRock[i]]}  drainage ${(sim.flow.area[i] / 1e6).toFixed(3)} km²  ` +
    `vegetation ${(100 * sim.veg[i]).toFixed(0)}%  dz/dt ${(sim.dzdt[i] * 1000).toFixed(3)} mm/yr` +
    (lake > 0.05 ? `  lake ${lake.toFixed(1)} m deep` : '') +
    (sim.seaLevel != null && sim.z[i] < sim.seaLevel ? `  sea ${(sim.seaLevel - sim.z[i]).toFixed(1)} m deep` : '') +
    (sim.waveExposureField[i] > 0 ? `  waves ${(100 * sim.waveExposureField[i]).toFixed(0)}%` : '') +
    (sim.notch[i] > 0 ? `  notch ${sim.notch[i].toFixed(1)} m` : '');
});

// dragging turns the globe
let drag = null;
mapCanvas.addEventListener('pointerdown', (ev) => {
  if (!isPlanet(sim) || pm.projection !== 'globe' || zoomPending) return;
  drag = { x: ev.clientX, y: ev.clientY, lon: globe.lon, lat: globe.lat, moved: false };
  mapCanvas.setPointerCapture(ev.pointerId);
});
mapCanvas.addEventListener('pointermove', (ev) => {
  if (!drag) return;
  const r = mapCanvas.getBoundingClientRect();
  const k = 180 / r.width;
  const dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
  if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
  globe.lon = ((drag.lon - dx * k + 540) % 360) - 180;
  globe.lat = Math.max(-89, Math.min(89, drag.lat + dy * k));
  if (!drag.pending) {
    drag.pending = true;
    requestAnimationFrame(() => { if (drag) drag.pending = false; resize(); draw(); });
  }
});
mapCanvas.addEventListener('pointerup', () => { if (drag) { drag.wasDrag = drag.moved; setTimeout(() => { drag = null; }, 0); } });

mapCanvas.addEventListener('click', (ev) => {
  if (drag && drag.wasDrag) return;
  const i = cellAtEvent(ev);
  if (i < 0) return;
  if (zoomPending) {
    zoomPending = false;
    try {
      let child, opts;
      if (isPlanet(sim)) {
        const n = +$('zoom-cells').value;
        opts = { lon: deg(sim.mesh.lon[i]), lat: deg(sim.mesh.lat[i]), factor: +$('zoom-factor').value, cols: n, rows: n };
        child = zoomFromPlanet(sim, opts);
      } else {
        opts = { center: [sim.grid.x[i], sim.grid.y[i]], factor: +$('zoom-factor').value, cols: +$('zoom-cells').value, rows: +$('zoom-cells').value };
        child = zoomIn(sim, opts);
      }
      stack.push({ sim, opts });
      setSim(child);
    } catch (e) {
      $('zoom-compare').textContent = e.message;
      syncZoom();
    }
    return;
  }
  if (isPlanet(sim)) sectionLat = deg(sim.mesh.lat[i]);
  else sectionRow = sim.grid.row(i);
  draw();
});

/** Outline of a planet patch as longitude/latitude points, for drawing on the planet. */
function patchOutline(child) {
  const { lon, lat } = child.zoom.planet;
  const proj = azimuthalProjection(lonLatToVector(lon, lat), child.zoom.planetRadius ?? 6.371e6);
  const b = child.grid.bounds;
  const pts = [];
  const N = 24;
  const edge = (x0, y0, x1, y1) => { for (let k = 0; k < N; k++) pts.push([x0 + ((x1 - x0) * k) / N, y0 + ((y1 - y0) * k) / N]); };
  edge(b.xmin, b.ymin, b.xmax, b.ymin); edge(b.xmax, b.ymin, b.xmax, b.ymax);
  edge(b.xmax, b.ymax, b.xmin, b.ymax); edge(b.xmin, b.ymax, b.xmin, b.ymin);
  pts.push(pts[0]);
  return pts.map(([x, y]) => {
    const v = proj.inverse(x, y);
    return [deg(Math.atan2(v[1], v[0])), deg(Math.asin(Math.max(-1, Math.min(1, v[2]))))];
  });
}

$('zoom-in').addEventListener('click', () => {
  zoomPending = !zoomPending;
  syncZoom();
});
$('zoom-out').addEventListener('click', () => {
  if (!stack.length) return;
  const { sim: parent } = stack.pop();
  lastPatch = isPlanet(parent)
    ? { parent, outline: patchOutline(Object.assign(sim, { zoom: { ...sim.zoom, planetRadius: parent.params.radius } })) }
    : { parent, bounds: { ...sim.grid.bounds } };
  setSim(parent);
});

$('run').addEventListener('click', () => {
  running = !running;
  $('run').textContent = running ? 'Pause' : 'Run';
});
$('step').addEventListener('click', () => { sim.step(); draw(); });
$('quake').addEventListener('click', () => {
  if (isPlanet(sim) || !sim.params.earthquakes.faults.length) return;
  sim.earthquake(0, 7);
  sim.quakes.lastStep = 1; // highlight it
  draw();
});
$('reset').addEventListener('click', () => {
  if (stack.length) { // re-make this patch from its frozen parent
    const top = stack[stack.length - 1];
    setSim(isPlanet(top.sim) ? zoomFromPlanet(top.sim, top.opts) : zoomIn(top.sim, top.opts));
    return;
  }
  if (isPlanet(sim)) {
    restart(resolvePlanetParams({ ...sim.params, ...planetOverrides(), plates: { ...sim.params.plates, count: +$('plates').value } }));
    return;
  }
  const base = { ...sim.params };
  restart(resolveParams({ ...base, ...gridOverrides(), grid: gridOverrides().grid, time: { ...base.time, dt: +$('dt').value } }));
});
$('dt').addEventListener('change', () => { if (!isPlanet(sim)) sim.setParams({ time: { dt: +$('dt').value } }); });
$('planet-dt').addEventListener('change', () => { if (isPlanet(sim)) sim.setParams({ time: { dt: +$('planet-dt').value } }); });
for (const id of ['cols', 'rows', 'dx', 'subdivisions']) $(id).addEventListener('input', updateExtent);
$('projection').addEventListener('change', () => { resize(); draw(); });

function refreshJson() {
  $('json').value = JSON.stringify(sim.params, null, 2);
}
$('json-refresh').addEventListener('click', refreshJson);
$('json-apply').addEventListener('click', () => {
  try {
    const raw = JSON.parse($('json').value);
    const p = raw.subdivisions != null ? resolvePlanetParams(raw) : resolveParams(raw);
    stack = []; lastPatch = null;
    if (raw.subdivisions != null) fillPlanetInputs(p); else fillGridInputs(p);
    restart(p);
  } catch (e) {
    $('json-error').textContent = e.message;
  }
});

$('copy-params').addEventListener('click', async () => {
  const text = JSON.stringify(sim.params, null, 2);
  try {
    await navigator.clipboard.writeText(text);
    $('copy-status').textContent = 'Parameters copied. Paste them into the JSON box (or a file) to reproduce this run.';
  } catch {
    // clipboard refused: open the JSON panel with the text selected instead
    refreshJson();
    $('json').closest('details').open = true;
    $('json').select();
    $('copy-status').textContent = 'Copying was blocked; the parameters are selected in the JSON box above.';
  }
});
// Single-file builds published where downloads are blocked set this flag.
if (window.EDGY_NO_DOWNLOADS) for (const el of document.querySelectorAll('[data-download]')) el.hidden = true;

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
  if (isPlanet(sim)) return;
  download(`edgy-snapshot-t${Math.round(sim.time)}.json`, new Blob([JSON.stringify(sim.snapshot())], { type: 'application/json' }));
});
$('load').addEventListener('change', async (ev) => {
  const file = ev.target.files[0];
  if (!file) return;
  try {
    const data = JSON.parse(await file.text());
    stack = []; lastPatch = null; // a loaded run stands alone (a saved patch keeps its edges and inflows)
    if (data.format?.startsWith('edgy-erosion-snapshot')) {
      setSim(Simulation.fromSnapshot(data));
    } else if (data.subdivisions != null) {
      const p = resolvePlanetParams(data);
      fillPlanetInputs(p);
      restart(p);
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
const initial = resolveParams({ preset: 'mesa-country' });
fillGridInputs(initial);
fillPlanetInputs(resolvePlanetParams({}));
restart(initial);
requestAnimationFrame(loop);

// handy for poking at the simulation from the browser console
window.edgy = { get sim() { return sim; }, Simulation, Planet, resolveParams, resolvePlanetParams };
