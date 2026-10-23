# Edgy Erosion

A landscape-evolution (erosion) simulator on a **hexagonal grid**, for generating terrain
with cliffs, escarpments, mesas, outcrops and talus slopes. It runs:

- **in a browser**, with live controls, map views and a geological cross-section, and
- **headless in Node.js**, for scripting, batch runs and tests. No browser is needed.

It is plain JavaScript (ES modules) with **no dependencies** and no build step. The same
`src/` files run in both places.

The model combines tectonic uplift, soil production, river incision and deposition,
soil creep, cliff collapse (rockfall) and talus avalanching. It works over layered,
dipping and intruded rock, and climate and vegetation modulate the rates. The physics is
documented, with citations and an honest list of simplifications, in
**[docs/MODEL.md](docs/MODEL.md)**.

## Quick start

You need [Node.js](https://nodejs.org) 18 or newer. Nothing to install beyond that.

```bash
git clone <this repo> && cd edgy-erosion

npm test                 # run the test suite (about 10 s)
npm run serve            # then open http://localhost:8000/
node bin/erode.js --list-presets
node bin/erode.js --preset mesa-country --steps 300 --section mid --out out/mesa
```

The browser app needs a local web server (`npm run serve`), not a double-clicked file,
because browsers refuse to load ES modules from `file://`. The repository root can also
be published as-is with GitHub Pages.

## In the browser

- **Preset** picks a landscape type. **Run/Pause**, **Step** and **Reset** control time.
- **Grid**: `cols × rows` cells with spacing `dx` metres. Extent and resolution are
  independent. These, and the seed, take effect on Reset.
- **Sliders** change uplift, climate, vegetation, river, hillslope and rockfall parameters
  *while it runs*.
- **View** switches between terrain, elevation, surface geology, soil depth, vegetation,
  drainage area, erosion rate and recent rockfall. Hover for cell values. Click to move
  the cross-section, which shows the strata beneath the surface.
- **All parameters (JSON)** shows the complete parameter set. That JSON *is* the
  experiment: paste it back, or load it from a file, to reproduce a run exactly.
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

See `examples/basic.js` for a runnable script. `src/render.js` turns state into RGBA
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
| `badlands` | weak, bare, thinly layered rock: dense gullies |
| `volcanic-necks` | soft sediment cut by basalt dikes and plugs: walls and towers emerge as the land wears down |
| `coastal-plateau` | layered plateau meeting the sea: deltas and drowned valleys |

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
  hexgrid.js         hex lattice geometry and neighbours
  flow.js            lakes, flow directions, drainage area
  lithology.js       rock type at any place and depth
  processes/         fluvial.js, hillslope.js, rockfall.js
  params.js          defaults, rock library, presets
  render.js          colouring, rasterising, cross-sections
web/                 browser front end (app.js, style.css); index.html at the root
bin/erode.js         command-line runner (+ png.js encoder)
examples/            small scripts
test/                node:test suites
docs/MODEL.md        the physics, its sources and its limits
```
