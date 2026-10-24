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

Each milestone closes with an automated gate.

1. **Planet (L0)**: plates, crustal thickening, isostasy, event log. *Gate: under 10 s
   in the browser; golden-seed hashes stable.* **Passed**: `node scripts/browser-gate.js`
   runs level 0 in a Web Worker in headless Chromium in 8.5-9.4 s, and its golden hashes
   match the ones `node --test` checks (`test/golden.json`).
2. **Region (L1) benchmark**: GPU stream power erosion and drainage graph extraction.
   *Gate: each tile under 1 s; Hack's law and slope-area checks pass.*
3. **Seams**: corner and edge records, river crossings, halos. *Gate: the
   order-independence test produces bitwise identical seams.*
4. **Landscape (L2)**: landslides, alluvial fans, glaciation, lithology contrasts.
   *Gate: local maps vary beyond noise; consistency holds within tolerance at every
   level.*
5. **Site (L3) and full zoom**: feature placement, rendering, caching, prefetch. *Gate:
   zooming from planet to 1 m produces no contradictions across levels.*

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

**Known artifact.** About 1.3% of cells are pinholes of ocean floor inside continents.
They appear where the nearest-neighbour look-up leaves a plate cell consumed inside a
collision zone. Two quick fixes made things worse and were reverted: filling interior
gaps from neighbours misfires at ridges, and refilling consumed cells creates crust
without bound. It needs a proper fix in the plate model, or a clean-up when level 1
reads level 0.

**Level 1, under 1 s per tile: met at 256 x 256, not at 1024 x 1024.**

| Tile | ms per step (flood routing) | of which routing | fast routing alone |
| --- | --- | --- | --- |
| 128 x 128 | 10 | 7 | 4 |
| 256 x 256 | 28 | 20 | 6 |
| 512 x 512 | 115 | 82 | 23 |
| 1024 x 1024 | 455 | 310 | 100 |

A level-1 run needs tens of implicit steps: a few drainage response times at
0.1-0.5 Myr per step. At 256 cells that is well under a second; at 1024 it is 10-20 s.
Routing dominates. The fast path (steepest descent and an O(n) walk up the drainage
tree) is 3-4x cheaper than the priority flood, but it applies only to surfaces without
depressions. Implicit diffusion leaves millimetre pits on filled-lake flats, so in
practice the flood runs every step.

**Next step for routing:** resolve depressions on the basin graph rather than flooding
the whole tile (Cordonnier, Bovy & Braun 2019, *Earth Surface Dynamics* 7), which keeps
routing O(n).

## Other deviations and caveats

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
