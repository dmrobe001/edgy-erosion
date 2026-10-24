# Hierarchical generator: plan, decisions and status

[DESIGN.md](DESIGN.md) is the design. This file records how the code base gets there:
what is reused from the first version, the decisions taken on the design's open
questions, the milestones, and measured numbers against the design's budgets.

The first version (a hex-grid erosion simulator with a planet, earthquakes, nested zoom
and a walking-scale detail layer) is tagged **`v1-hex-erosion`** (commit `566acc5`). It
stays in the tree and keeps working; the hierarchical generator grows beside it in
`src/world/`.

## What carries over

| Existing code | Role in the hierarchy |
| --- | --- |
| `src/sphere.js` (icosahedral Voronoi mesh) | Level-0 mesh, unchanged. |
| `src/planet.js` (kinematic plates, isostasy) | **The level-0 engine.** Extended with material tracking (arc and accreted shares of each crust column, giving the design's material classes), a dated event log (collisions, continental rifting), and a lean mode without river routing. Made about 3x faster (see below). |
| `src/flow.js` (`MinHeap`) | Priority queue for floods at every level. |
| `src/nest.js` (block-mean-preserving prolongation) | To become the level-to-level prolongation, with restriction R = block average checked after every child run, as the design asks. Milestone 2. |
| `src/processes/rockfall.js`, `hillslope.js`, `fluvial.js` (deposition), `coast.js` | Physics to port to level 2 on square rasters: Culmann-threshold slope failure is the design's threshold landsliding; fans are the deposition term. |
| `src/processes/earthquakes.js`, the plate-boundary faults in `nest.js` | Fault traces from the event log and plate boundaries at levels 1-2. |
| `src/meander.js`, `src/detail.js` | Level 3: channel width and meanders from discharge, beaches, tree crowns. They need to honour pinned river crossings. |
| `src/lithology.js` | Turning material classes into strata at levels 1-2. |
| `src/simulation.js` on the hex grid, `web/` | Version 1, kept as is. The new levels use square rasters, because seams (edge records, crossings, halos) are far simpler on rows of cells. A new viewer comes in milestone 5. |
| `src/rng.js` (stateful generator) | Not used by new code. Everything in `src/world/` draws from the stateless `src/world/hash.js`. The planet still uses its own seeded generator internally; that stays deterministic because level 0 is one sequential computation seeded from the world key. |

## Decisions on the open questions (provisional)

- **Sphere grid.** Icosahedral at level 0, where the existing plate model lives and
  cells are near-uniform. Cube-sphere quadtree-style addresses from level 1 down, where
  tiles must be square rasters with simple edges (`src/world/cubesphere.js`,
  equal-angle projection: cell widths within 1.41x over a face).
- **Refinement and tile size.** 32x per level in cell size, as designed, but tiles of
  **256 x 256** cells instead of 1024 x 1024. The first level has 32 tiles per face
  side, so cell sizes match the design's: 1.2 km, 38 m, 1.2 m. Tile size in cells is a
  performance knob, independent of the cell sizes (`WORLD_LAYOUT`).
- **L3.** One level for now.
- **Time.** Present only.

## Milestones

Each milestone closes with an automated gate. Node gates run with `node --test`; browser
gates with `node scripts/browser-gate.js` (headless Chromium through Playwright).

1. **Planet (L0)**: plates, crustal thickening, isostasy, event log. *Gate: under 10 s
   in the browser; golden-seed hashes stable.* **Passed.** Level 0 runs in a Web Worker
   in headless Chromium in 8.5-9.4 s, and its field hashes match the golden ones Node
   checks (`test/golden.json`, now `hg-0.2`).
2. **Region (L1) benchmark**: GPU stream power erosion and drainage graph extraction.
   *Gate: each tile under 1 s; Hack's law and slope-area checks pass.* **Partly
   passed.** The WGSL solver runs and agrees with its CPU reference; Hack's law passes
   on every tile tried (exponents 0.50-0.67); slope-area concavity falls within
   0.3-0.7 on most tiles (0.33-0.59 in the Node runs) but not all (below). **The 1 s
   budget cannot be checked here**: this container's only WebGPU adapter is
   SwiftShader, a CPU emulation, where a tile takes about 30 s, the same as the
   JavaScript reference. A real GPU should do it in a fraction of a second (200
   iterations of about 20 kernel dispatches over 66 000 nodes, plus about 0.3 s of CPU
   work per tile), but that is an estimate, not a measurement. To measure it, run
   `npm run serve` and open `/web/gate/region.html` in a browser with a GPU.
3. **Seams**: corner and edge records, river crossings, halos. *Gate: the
   order-independence test produces bitwise identical seams.* **Passed**
   (`test/seams.test.js`, and across a level-1 seam at level 2 in
   `test/landscape.test.js`). Halos: with edges pinned to shared records, a halo cannot
   affect a tile's interior at level 1, so seam mode uses none (the halo option remains
   for tiles without seams).
4. **Landscape (L2)**: landslides, alluvial fans, glaciation, lithology contrasts.
   *Gate: local maps vary beyond noise; consistency holds within tolerance at every
   level.* **Passed** (`test/landscape.test.js`). Across rock and climate histories,
   7 of 9 summary statistics spread more than twice as much as across noise seeds,
   most by 10-100x. Block means hold within 2 m at level 1 and within 1 m at level 2.
5. **Site (L3) and full zoom**: feature placement, rendering, caching, prefetch. *Gate:
   zooming from planet to 1 m produces no contradictions across levels.* **Passed in
   Node** (`test/zoom.test.js`): block means within tolerance at every level, coasts
   agreeing where the parent is clearly land or sea, every big parent river carried by
   the child within a parent cell, snow only where cold, site seams identical. The
   viewer (`/web/world/`) zooms planet to 1.2 m with IndexedDB caching and neighbour
   prefetch (`test/pipeline.test.js` checks that a cached tile equals a regenerated one).

## Measured against the budgets

From `node scripts/bench-world.js`, on the development container (Node 22, one core).

**Level 0, 10 s budget: met, at coarser cells than designed.** It takes about 9 s in
Node and 8.5-9.4 s in a Chromium worker for 23 042 cells about 160 km apart, 500 Myr
in 250 steps of 2 Myr. The design asks for
30-50 km cells. That is about 16x the cells at 40 km, and about 4x the steps (plates may
move only half a cell per substep), so roughly 60x the work: far out of reach on one
CPU core. Each level-1 tile still sees about 2 x 2 level-0 cells. The coarse level-0
map then leaves more of the regional structure to level 1, which is where it is
cheapest to make.

The speed-ups that got level 0 to about 9 s (from about 50 s), none of which changes
the model:
- each world cell looks up only the plates that owned it or a neighbour last substep,
  not all of them;
- arc crust is spread only over cells within three arc widths of where it was made,
  and crustal flow runs only within reach of crust above the flow threshold;
- the explicit crustal-flow substeps are capped, and the flow law's speed-up is capped
  with them;
- river routing is replaced by elevation-dependent denudation (Ahnert 1970) above a
  400 m floor. The eroded sediment is laid down in the nearest sea cell as continental
  crust, so margins build shelves and sedimentary basins, and the crust budget still
  closes.

Level-0 settings that shape the result (`L0_DEFAULTS` in `src/world/l0.js`):
- plates get new poles every 100 Myr, giving supercontinent cycles; without it, every
  continent ends up in one polar mass after 500 Myr;
- crust spreads sideways above 55 km and is capped at 70 km, so plateaus stand at
  about 3-6 km;
- 40% of the surface starts continental.

With these, three seeds gave 33-36% land, land elevations with median 0.6-1.0 km,
90th percentile about 3 km and 99th percentile 3.3-4.6 km, and 150-200 dated events.

**Pinholes.** The plate model leaves about 1.3% of cells as pinholes of ocean floor
inside continents, where its nearest-neighbour look-up consumes a plate cell inside a
collision zone. Since `hg-0.2` they are closed when the planet is packed: a sea cell
with at most two sea neighbours takes its land neighbours' mean, over three passes.
Two of about 8 000 land cells remain. The plate model itself still makes them.

**Level 1.** The implicit CPU solver (`fastscape.js`, priority-flood routing):

| Tile | ms per step (flood routing) | of which routing | fast routing alone |
| --- | --- | --- | --- |
| 128 x 128 | 10 | 7 | 4 |
| 256 x 256 | 28 | 20 | 6 |
| 512 x 512 | 115 | 82 | 23 |
| 1024 x 1024 | 455 | 310 | 100 |

The GPU-shaped solver (`jacobi.js`, `gpu.js`) runs every pass data-parallel. Per
iteration it does an uplift pass, 4 fill sweeps, one receiver pass, 4 area sweeps,
8 erosion sweeps, creep and a restart. A 257 x 257 tile takes 200 iterations: about
12-30 s in the JavaScript reference and about 30 s on SwiftShader. Its steady state is
the implicit scheme's exactly (`test/region.test.js`). Two numerical flaws turned up on
the way and were fixed:
- refilling fresh incision to the old surface made a slope floor of U dt / L;
- a one-cell-per-iteration limit held back knickpoints.

**Level 2:** 1.5-2 s per tile on the CPU (budget 1 s), most of it in routing filled
lakes. **Level 3:** about 0.5 s per tile.

## Other deviations and caveats

- **GPU results are not bitwise the CPU's.** WGSL's `pow` and division need not round
  exactly like JavaScript's. Near-ties in receiver choice then flip, and elevations
  differ by about 1-8 m rms. The design's cross-hardware determinism therefore needs
  one canonical path per world. Fixed-point arithmetic would make the GPU path exact.
  Shared seam records never touch the solver, so seams stay bitwise everywhere.
- **Tile edges show.** Edges are walls except where a parent river crosses, so they
  tend to become straight divides. Minor stream crossings at low points of each edge
  profile are implemented (`minorSpacing`) but off by default: they raised local base
  levels and flattened the tiles' drainage (concavity fell from 0.56 to 0.1 on the
  test tile). The design accepts some of this ("a child-scale channel running parallel
  to an edge has no clean representation").
- **Slope-area concavity varies between tiles.** It is 0.47-0.59 on uplifted tiles.
  Tiles whose boundary hands them a river high above the sea, or a large closed basin,
  are still adjusting after 5 Myr and give 0.14-0.3. Natural landscapes out of
  equilibrium vary as much.
- **Smoothing the planet.** Level 1 samples the planet with a Gaussian kernel (linear
  interpolation left 160-km creases), and a first-order correction restores most of
  what smoothing takes off peaks. The rest is a per-block target bias; the erodibility
  controller and 40 final correction passes steer the child onto it.
- **Large rivers at the site level** follow the parent's 38-m steepest-descent path
  smoothed, so a 170-m-wide river bends at angles. A river parallel to a tile edge
  leaves a strip of uncut land along the edge.

- **Determinism across machines.** Integer hashing is exact everywhere. Floating-point
  results also depend on the engine's `Math.exp`, `Math.sin` etc.; V8 and SpiderMonkey
  both derive theirs from fdlibm, but the language standard does not promise it. The
  golden hashes were stored under Node 22 (V8); running them in Firefox is part of
  milestone 5.
- **Event regions** are spherical caps (centre and angular radius), not polygons.
- **Regolith at level 0** is a placeholder rule (warm, wet and low means deeper), not a
  model. Level 1 should produce its own.
- **Climate at level 0** is zonal: a temperature caricature with a 6.5 K/km lapse rate,
  and the planet's zonal rainfall bands. Orographic precipitation belongs to level 1, as
  designed.
