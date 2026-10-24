# Edgy Erosion

A landscape-evolution (erosion) simulator on a **hexagonal grid**, for generating terrain
with cliffs, escarpments, mesas, outcrops and talus slopes, plus a **planet** scale above
it: plate tectonics on an icosahedral sphere that you can zoom into, region by region,
down to cliffs. It runs:

- **in a browser**, with live controls, map views and a geological cross-section, and
- **headless in Node.js**, for scripting, batch runs and tests. No browser is needed.

It is plain JavaScript (ES modules) with **no dependencies** and no build step. The same
`src/` files run in both places.

The model combines tectonic uplift, faults that slip in earthquakes, soil production,
river incision and deposition, soil creep, cliff collapse (rockfall, also shaken loose
by earthquakes), talus avalanching, and wave attack on coasts (sea cliffs, shore
platforms, beaches, marine terraces under changing sea level). It works over layered,
dipping and intruded rock, and climate and vegetation modulate the rates. At planet
scale, rigid plates spread, subduct and collide, crust floats isostatically, and the
same rivers and creep wear the continents down. The physics is
documented, with citations and an honest list of simplifications, in
**[docs/MODEL.md](docs/MODEL.md)**.

## Direction: a hierarchical terrain generator

The project is being rebuilt around a new design, **[docs/DESIGN.md](docs/DESIGN.md)**:
a planet built once from plate tectonics, then refined on demand, level by level, down
to metre-scale sites. Every level must agree with the coarser one above it, and the
result is a pure function of the world seed and the tile address. Progress, decisions
and measured budgets are in **[docs/PLAN.md](docs/PLAN.md)**; the new code lives in
`src/world/`:

```js
import { generatePlanet } from './src/world/l0.js';
const planet = generatePlanet('my world'); // ~10 s: fields, river graph, dated events
```

Everything below describes the first version, which stays in the tree and is tagged
`v1-hex-erosion`.

## Quick start

You need [Node.js](https://nodejs.org) 18 or newer. Nothing to install beyond that.

```bash
git clone <this repo> && cd edgy-erosion

npm test                 # run the test suite (about 20 s)
npm run serve            # then open http://localhost:8000/
node bin/erode.js --list-presets
node bin/erode.js --preset mesa-country --steps 300 --section mid --out out/mesa
```

The browser app needs a local web server (`npm run serve`), not a double-clicked file,
because browsers refuse to load ES modules from `file://`. The repository root can also
be published as-is with GitHub Pages.

For a version you can just double-click or email, `npm run build` bundles everything
into one self-contained file, `dist/edgy-erosion.html` (about 250 KB, works offline).
Rebuild it after changing `src/` or `web/`.

## In the browser

- **Preset** picks a landscape type. **Run/Pause**, **Step** and **Reset** control time.
- **Grid**: `cols × rows` cells with spacing `dx` metres. Extent and resolution are
  independent. These, and the seed, take effect on Reset.
- **Sliders** change uplift, climate, vegetation, river, hillslope and rockfall parameters
  *while it runs*.
- **View** switches between terrain, elevation, surface geology, soil depth, vegetation,
  drainage area, erosion rate, recent rockfall and wave exposure. Hover for cell values. Click to move
  the cross-section, which shows the strata beneath the surface.
- **All parameters (JSON)** shows the complete parameter set. That JSON *is* the
  experiment: paste it back, or load it from a file, to reproduce a run exactly.
- **Zoom** refines a patch: choose a refinement factor and size, press *Zoom in…* and
  click the map. The level above freezes and supplies the patch's edges and inflowing
  rivers. The panel compares the patch with the level above it. *Zoom out* returns to
  the level above and outlines the patch. You can zoom several levels deep.
- **Planets** are at the bottom of the preset list. Views show terrain, plates and their
  boundaries (convergent, divergent, transform), ocean-floor age and crust thickness, on
  a map or a globe you can drag. The cross-section runs along a parallel and shows crustal
  roots under mountains. *Zoom in…* on a planet lays a flat regional patch (factor 8–64
  finer) on it, inheriting its topography, uplift and rainfall. Plate boundaries crossing
  the patch become faults (thrusts where plates converge, normal faults where they
  part) that slip in earthquakes. Zoom again to reach landscape scale.
- **Walking-scale detail** (below the cross-section): *Build detail* turns the current
  landscape into a map at 0.5–4 m per pixel. Rivers become curves with real widths that
  meander, with point bars, oxbows and floodplains; sandy coasts get beaches and dunes;
  forests get tree crowns. Drag to pan; click the map above to jump. Try the
  `river-coast` preset.
- **Earthquake!** (for presets with faults, e.g. `earthquake-scarp`) slips the first
  fault in a magnitude-7 event; ruptures flash yellow on the map.
- **Export** saves the current view as PNG, or a snapshot (full state plus parameters)
  that can be resumed in the browser or the CLI.

## From code

```js
import { Simulation } from './src/index.js';

const sim = new Simulation({
  preset: 'cuestas',              // optional starting point
  seed: 7,
  grid: { cols: 200, rows: 160, dx: 20 },
  climate: { precipitation: 0.5, frost: 0.6 },
});

sim.run(500);                     // or sim.step() in your own loop
sim.setParams({ uplift: { rate: 0 } });   // change conditions mid-run
sim.run(500);

sim.z;            // Float64Array of surface elevations (m), one per cell
sim.soil;         // soil/talus thickness (m)
sim.flow.area;    // drainage area (m², precipitation-weighted)
sim.grid.x[i], sim.grid.y[i];     // cell centre (m)
sim.stats();      // summary numbers
const snap = sim.snapshot();      // JSON-able; Simulation.fromSnapshot(snap) resumes bit-for-bit
```

Zoom into a patch at finer resolution:

```js
import { zoomIn, compareWithParent } from './src/index.js';
const child = zoomIn(sim, { center: [2000, 1500], factor: 4, cols: 128, rows: 128 });
child.run(200);                   // the parent stays frozen
compareWithParent(child, sim);    // mean elevation, rates, river discharge vs the parent
```

Processes run on any mesh that satisfies the contract in `src/mesh.js`:
`new Simulation(params, { mesh })`.

A planet, and the way down from it:

```js
import { Planet, zoomFromPlanet, zoomIn } from './src/index.js';
const planet = new Planet({ preset: 'earth-like', subdivisions: 40, seed: 3 });
planet.run(200);                                  // 100 Myr at dt = 0.5 Myr
planet.sim.z; planet.owner; planet.age; planet.thickness;   // per cell of planet.mesh
const region = zoomFromPlanet(planet, { lon: 20, lat: -10, factor: 32 });  // ~6 km cells
const land = zoomIn(region, { center: [0, 0], factor: 8 });                  // ~750 m cells
```

Earthquakes: give faults in the parameters (positions are fractions of the domain and
may lie outside it):

```js
const sim = new Simulation({ preset: 'humid-hills',
  earthquakes: { faults: [{ x0: -1, y0: 0.5, x1: 2, y1: 0.4, dip: 60, slipRate: 2e-3 }] } });
sim.run(200);
sim.earthquake(0, 7.2);        // or trigger one now
sim.quakes.events;             // the catalogue
```

Walking-scale maps from a landscape run:

```js
import { Simulation, Detail, detailColors } from './src/index.js';
const sim = new Simulation('river-coast');
sim.run(200);
const detail = new Detail(sim);                // rivers as curves; big ones meander
const tile = detail.render({ x0: 600, y0: 900, width: 900, height: 700, pixel: 1 });
const rgba = detailColors(detail, tile);       // tile.z, tile.cover, tile.level as arrays too
```

See `examples/basic.js`, `examples/zoom.js`, `examples/planet.js` and `examples/detail.js` for runnable scripts. `src/render.js` turns state into RGBA
pixels without a canvas, and `bin/erode.js` uses it to write PNGs from Node.

### Command line

```
node bin/erode.js [--preset NAME] [--params FILE] [--steps N] [--dt YEARS]
                  [--cols N --rows N --dx METRES] [--seed N] [--out DIR]
                  [--views terrain,geology,...] [--section ROW|mid] [--frames N]
                  [--snapshot] [--heightmap]
```

It writes PNG images, `params.json` (the exact resolved parameters), and optionally a
snapshot and a 16-bit PGM heightmap. `--params` accepts either a parameter file or a
snapshot, so a run can be continued. Run `node bin/erode.js --help` for details.

## Presets

| name | what it shows |
|---|---|
| `mesa-country` | hard sandstone over soft shale, dry: stepped escarpments, mesas, buttes, talus |
| `alpine` | fast uplift of granite, cold and wet: steep ridges and frost-driven rockfall |
| `humid-hills` | weathered schist under forest: rounded, soil-mantled hills |
| `cuestas` | gently dipping hard/soft layers: asymmetric ridges with scarp and dip slopes |
| `fault-scarp` | a block rising along a fault: a degrading scarp, facets, alluvial fans |
| `earthquake-scarp` | a normal fault slipping in M 5.5–7.3 earthquakes: a stepped range front, offset strata, shaken-down cliffs |
| `badlands` | weak, bare, thinly layered rock: dense gullies |
| `volcanic-necks` | soft sediment cut by basalt dikes and plugs: walls and towers emerge as the land wears down |
| `coastal-plateau` | layered plateau meeting the sea: deltas and drowned valleys |
| `rocky-coast` | bands of hard and soft rock meeting the sea end-on: headlands, bays, sea cliffs, platforms |
| `river-coast` | a humid coastal plain with two large rivers from inland and a sandy coast; made for walking-scale detail |
| `marine-terraces` | a coast rising 1 mm/yr under 41-kyr sea-level cycles: a staircase of uplifted wave-cut platforms |

Planet presets (`src/planet.js`): `earth-like` (a third continental crust in scattered
blocks, ten plates), `supercontinent` (one continent across several plates, which tear
it apart), `archipelago` (small fragments, many plates, growing island arcs).

Presets are only starting points. Every number in them lives in `src/params.js`, and the
units and meanings are commented there.

## Parameters at a glance

| group | key parameters | what they control |
|---|---|---|
| `grid` | `cols`, `rows`, `dx` | domain size and lattice spacing (m) |
| `time` | `dt` | years per step |
| `boundary` | `open`, `seaLevel` | which edges are base level; sea |
| `initial` | `type`, `height`, `noise`, `soil` | starting landform |
| `uplift` | `rate`, `pattern` | rock uplift (m/yr): uniform, dome, tilt or fault |
| `rocks` | per rock: `erodibility`, `cohesion`, `frictionAngle`, `weathering` | how each rock erodes, weathers and holds cliffs |
| `strata` | `layers`, `top`, `dip`, `dipDirection`, `repeat`, `undulation`, `strengthVariation` | layered geology |
| `bodies` | `disc` / `dike` shapes | plugs and dikes cutting the strata |
| `climate` | `precipitation`, `precipitationGradient`, `frost` | runoff; freeze–thaw |
| `vegetation` | `cover` (`'auto'` or 0–1) and its effects | protects soil, thickens it, steepens stable slopes |
| `fluvial` | `m`, `n`, `soilErodibility`, `deposition`, `routing` | rivers |
| `hillslope` | `diffusivity`, `criticalSlope` | soil creep and angle of repose |
| `weathering` | `productionDepth`, `frostBoost`, `moistureExponent` | soil production |
| `rockfall` | `enabled`, `timescale`, `frostFactor` | cliff collapse |
| `earthquakes` | `faults` (trace, `dip`, `slipRate`), `bValue`, `min/maxMagnitude`, `upliftShare`, `coseismicWeakening` | faults that slip in events |
| `coast` | `seaLevelAmplitude`, `seaLevelPeriod`, `seaLevelCurve`, `waveHeight`, `waves`, `tidalRange`, `marineDiffusivity`, `attrition` | sea level through time, wave climate, cliff retreat, shore sediment |

**Scale matters.** Cliffs and rockfall need `dx` small compared with the cliff heights
the rock can support: tens of metres or less. On a coarse grid, cliffs smear into slopes.
Changing `dx` changes the landscape, not just its resolution. See "Known limitations" in
[docs/MODEL.md](docs/MODEL.md).

## Reproducibility

- All randomness comes from a seeded generator (`src/rng.js`); `Math.random` is never used.
- Same parameters give bit-identical results within one JavaScript engine. Different
  engines (Chrome vs. Firefox vs. Node) may differ in the last bits of `Math.exp` and
  similar functions, so they can drift apart over long runs.
- `params.json` / the JSON panel / snapshots fully specify a run.
- `npm test` checks the numerics against analytical results: the stream-power steady
  state, isotropic diffusion, Culmann's closed form and the volume budget. It does not
  check geological realism; nobody has calibrated these landscapes against real terrain.

## Layout

```
src/                 the simulation (browser + Node)
  simulation.js      state, step order, bookkeeping, snapshots
  mesh.js            the mesh contract every process uses; an irregular test mesh
  hexgrid.js         hex lattice geometry and neighbours
  nest.js            zooming in: refine a patch of a frozen parent (or of a planet)
  sphere.js          the icosahedral geodesic mesh, projections
  detail.js          walking-scale maps: rivers as curves, channels, floodplains, beaches, cover
  meander.js         meander migration of river centre lines (Howard & Knutson)
  planet.js          plates, isostasy, sea level, surface processes on the sphere
  flow.js            lakes, flow directions, drainage area
  lithology.js       rock type at any place and depth
  processes/         fluvial.js, hillslope.js, rockfall.js, coast.js, earthquakes.js
  params.js          defaults, rock library, presets
  render.js          colouring, rasterising, cross-sections
web/                 browser front end (app.js, style.css); index.html at the root
bin/erode.js         command-line runner (+ png.js encoder)
scripts/             serve.js (dev server), build-single.js (one-file bundle)
dist/                the built single-file app
examples/            small scripts
test/                node:test suites
docs/MODEL.md        the physics, its sources and its limits
```
