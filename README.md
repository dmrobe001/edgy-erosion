# Procedural Earth-like worlds from geological history: L0 plate tectonics, L1 landmasses

This repository holds the first two levels of a multi-level procedural generator for Earth-like game worlds, plus sketches of the levels below them. L0 is a planet-scale plate tectonics simulation. L1 takes one landmass at the present and replays the last 10–30 Myr of its history at finer resolution, with its crust deforming as a thin viscous sheet driven by L0's plates, so collisions crumple it into mountain belts and plateaus; it then settles the landmass's drainage basins on a 10 km grid, the tiles for L2. Both run in a browser (`plates.html`, a single self-contained page) and in Node.js for headless testing.

The goal is diversity that comes from **material and history** rather than from noise layered on hand-made rules. Each level simulates or procedurally generates the structure that is plausible at its own scale. It hands down only what it can actually resolve, plus a record of what happened, so finer levels inherit causes as well as shapes.

Contents: the big picture; quick start; L0 in detail; the L1 handoff format; L1 landmass; sketches of the levels below; known limitations and next steps. `DEVLOG.md` explains how the current version came about, including what didn't work. `ONTOLOGY.md` labels every structure in the code as field, material, record, process, source, event or pattern, and lists where the labels show a problem.

---

## 1. The big picture

### Levels as conditional generators

Each level takes the coarser level's output as a condition and produces finer structure consistent with it. A level does not need to reach steady state. It only needs to replace the interpolation it inherits with structure that is plausible at its own resolution. Discrete, procedurally defined events are welcome wherever they produce richness that would be expensive to simulate physically.

| Level | Cell size | What it owns | Representative methods |
|---|---|---|---|
| L0 planet | ~80–160 km | Plates, crust type and age, rock columns, fault history, events | Rigid plates on a sphere (this repository) |
| L1 landmass | ~10–50 km | Crumple zones, drainage basins, climate, biomes, fault strands, glaciation | Thin viscous sheet, stream-power erosion, orographic rainfall |
| L2 region | ~500 m | River networks, folded and tilted strata, water table | Flow routing, layer geometry |
| L3 landscape | ~50 m | Hillslopes, landslides, terraces, cave networks, springs | Diffusion, event catalogues, graph search |
| L4 site | ~5 m | Cliff bands, talus, meanders, vegetation instances | Differential erosion, point processes |
| L5 detail | ~1 m | Small cliffs, stream beds, local texture | Lazy, per-tile, seeded |

A ratio of roughly 5–10× between levels gives five or six levels from planet to metre scale.

### Three principles that shaped L0

1. **Each level owns a band of wavelengths.** No grid model resolves features smaller than two or three of its cells; detail below that is numerical noise. A handoff passes down only the band the coarser level resolves, and the finer level regenerates everything below it. This is the same logic as multigrid methods, and the reason the L1 handoff contains *smoothed* fields.
2. **Pass down history, not just state.** A snapshot of elevation cannot tell a finer level whether a mountain belt came from a head-on or an oblique collision, or which valleys follow dead faults. L0 records rock columns, fault episodes, a dated event catalogue, and current plate motions, because those are what finer levels need to make cliffs, caves, springs and fault-controlled valleys.
3. **Categorical fields travel as geometry or proportions, not labels.** A rock type assigned per 100 km cell would become a giant checkerboard one level down. The handoff correlates rock units into formations whose *thickness* can be interpolated, so rock boundaries at L1 emerge where L1's own terrain cuts through smooth contacts. This is coarse-graining as in soft matter: hand down an order-parameter field, not particle identities.

### Reproducibility

Everything is deterministic given the seed and settings (seed text, resolution, initial plate count, start mode, continental crust target, net-rotation option). The same inputs always produce the same world.

**What that is for, and what it is not.** Determinism is a tool: it lets a session be rerun, lets L1 replay L0's history exactly, and lets tests catch accidental changes. It is not a goal in itself. Getting the physics right comes first, and a change that improves the physics is welcome even when it changes every number downstream. Tests should check physical expectations and invariants (budgets that close, ranges against Earth) more than exact matches to an earlier build; a bit-for-bit match to a wrong answer is worth nothing. Sessions saved by an earlier build rerun with that build's sources (`tools/session.js` finds the commit).

---

## 2. Quick start

**In a browser:** open `plates.html`. It needs no server and no build step; it loads only Google Fonts. Press **Run**: L0 advances by the **Gyr** slider (default 1) and stops; **Run** again continues for that long. Explanations of every view and control are under **Info** at the bottom of the page.

- **Views:** Elevation, Rock type, L1 preview, Plates. **Marks** shows hotspots and plate velocity arrows.
- **Projections:** a globe (drag to orbit) or an interrupted Goode homolosine map (drag to pan), whose cuts are placed through ocean for the current world: the map edge on the meridian crossing the least land, one more cut in the north and two in the south, each on the least-land meridian of its hemisphere at least 50° from the others. Scroll or pinch to zoom on both, and on the L1 map.
- **Tap or click without dragging:** selects that point. Its landmass is highlighted, the panel under the map shows its rock column, fault episodes and plate motion, and the L1 panel shows the landmass at once, as L0 has it now.
- **Refine selected landmass** (top of the L1 panel) replays the last 10–30 Myr of the selected landmass from a saved state of L0 (kept every 10 Myr while L0 runs, so run L0 for a while first). L0 replays the window twice (tracing the crust, then driving the sheet) and ends exactly where it started. The L1 map then shows terrain, thickening, finite strain, strain rate, surface rock, mantle pull, rainfall, drainage basins and the L2 tiles.
- **Extract handoff** traces the fault network, builds formations and domains, and overlays fault traces; **Download JSON** saves that L0 handoff.
- **Export session** saves a small JSON log from which the session reruns exactly (section 2.1). Downloads use the Claude artifact runtime's downloads capability; opened as a plain local file, they are unavailable.
- **Core** (under the world settings, applied by **Rebuild world**): the lattice core, or the marker core (section 3.14, Stage 1 of `L0_MESH_SPEC.md`), which carries the crust on markers on a Delaunay mesh. On the marker core, **Markers** draws every marker, with interface markers coloured by the kind of boundary; L1 refinement still needs the lattice core.
- Net rotation is always removed (section 3.8). The statistics the page no longer shows (plate count, continental fraction, mean seafloor age, units per column, net rotation, simulated rate) are in the session export.

**Rebuild the page** from its parts after editing `src/`:

```sh
./build.sh        # concatenates src/head.html + src/sim.js + src/mesh.js + src/l1.js + src/l2.js + src/ui.js into plates.html
```

**Headless tests** (Node.js 18+, no dependencies):

```sh
node tests/run_all.js              # quick pass over everything (about two minutes)
node tests/test_l1.js --n=96 --h=50 --win=20 --t0=1000 --land=3000   # L1 checks at the browser's defaults
node tests/test_budget.js --t=2400 # full-length crust budget run
node tests/bench.js --t=300        # ms per step at each resolution
node tests/test_mesh.js --n=64 --t=300       # the marker core's invariants (section 3.14)
node tools/compare_cores.js --n=96 --t=800   # the marker core against the lattice core, three seeds (about 15 minutes)
```

### 2.1 Session exports

**Export session** saves a JSON log, not the state: L0 is deterministic given its settings and step count, and L1 given the selected cell, the saved state it started from, the present it ran to, its window and its settings. A full state dump of a 1 Gyr session with a refined 42 M km² landmass was 240 MB raw and 37 MB compressed, while the log is a few hundred kB. The log holds:

- `settings`: seed, resolution, initial plates, start, continental target, steps per frame, Gyr, L1 settings.
- `reproduce`: L0's step count; the selected cell; for L1, the cell, `entrySteps` (the saved state), `presentSteps`, window and sheet settings, and whether L0 has moved on since (`stale`).
- `checks`: values to confirm a rerun: a hash of L0's elevation and plate ownership, plate count, the sheet's steps, markers and highest point, the basin count and inland share, the tile count.
- `viewports`: what is on screen and how. L0: colour mode, globe or map, camera (yaw, pitch, zoom), map zoom and centre, the Goode cuts, the selected cell and the last inspected point. L1: colour mode, zoom, pan, the visible extent in map km, the map frame, the last inspected point. L2: whether the tiles view is shown, the selected tile and its summary.
- `stats`, `actions` (build, run, stop, step, select, refine, extract, with step counts), `events`, the landmass's cells, and the L2 handoff summary.
- `build`: a hash of the page's sources (`build.sh` writes it); reruns need the same code.

`tools/session.js` reads it:

```sh
node tools/session.js session.json          # summary, and the commit whose sources match the build hash
node tools/session.js session.json --rerun  # rerun, compare every check, stash the results as session.stash
```

As a module, `reproduce(sess, SIM, L1)` returns the world at the session's present and L1's plan, sheet, drainage pass and handoff as they were. `stash` and `unstash` save and restore those (Node's v8 serialization, tens of MB) so later work can skip the rerun. The UI smoke test exports a session, reruns it in a fresh simulation, requires every check to match, and round-trips a stash.

### Repository layout

```
plates.html            the built single-file page (what the artifact serves)
build.sh               rebuilds plates.html from src/
src/head.html          markup, styles, explanatory notes; ends with an opening <script> tag
src/sim.js             the simulation core: a self-contained SIM module, usable in browser and Node
src/mesh.js            the marker core (L0_MESH_SPEC.md Stage 1): spherical Delaunay triangulation, MeshWorld (MESH module)
src/l1.js              L1: landmass pick, replay, thin viscous sheet with markers, flexure, crust budget, drainage pass, L2 handoff, basin selection (L1 module)
src/l2.js              L2: the relief L2 starts from, as a function of position; basin outlines burnt in at 500 m; the L2 pass (L2 module)
src/ui.js              rendering, controls, inspector, handoff panel, L1 and L2 panels
tests/                 headless tests and benchmarks (see section 3.13)
tools/session.js       reruns, checks and stashes a session exported from the page (section 2.1)
tools/compare_cores.js the marker core against the lattice core on the same seeds (section 3.14)
history/               simulation cores from earlier milestones, referenced by DEVLOG.md
ONTOLOGY.md            what each structure in the code is (field, material, pattern, ...), and what that shows
L0_MESH_SPEC.md        proposal: L0 rebuilt on Lagrangian markers with a torque balance, mantle flow, plumes and cracking
docs/                  diagnostic images produced during development
```

---

## 3. L0 in detail

### 3.1 Grid

The planet is an **equal-angle cube-sphere**: six faces of n × n cells projected onto a sphere of radius 6,371 km. Resolutions: n = 64 (24,576 cells, ~160 km), 96 (55,296, ~105 km), 128 (98,304, ~78 km).

- Lookup from any 3D direction to its cell is constant-time: pick the dominant axis, take the arctangent of the two minor coordinates, floor. `lookupFrac` also returns the fractional position within the cell, for bilinear interpolation.
- Each cell stores four neighbours, found by stepping one cell in each face direction (stepping past a face edge lands on the adjacent face).
- Cell areas vary by about ±25% and are computed exactly. They are used where integrals matter, notably net rotation.

### 3.2 Plates: body frames and coverage

Each plate stores its crust on **its own copy of the grid, in a frame that moves with the plate**. A plate has:

- an angular velocity vector Ω (rad/Myr), which defines its Euler pole and rotation rate;
- a cumulative rotation matrix Q, so a body-frame point b sits at world position Q·b;
- per-cell arrays: thickness, birth time, ocean sediment, crust type (continental or oceanic), a pointer into the rock-column pool, and two bookkeeping stamps.

Each step (Δt = min(1, 0.7 × cell angle / max angular speed) Myr, so plates move under one cell per step):

1. **Rotate** every plate: Q ← R(Ω, Δt) · Q, using the Rodrigues formula, re-orthonormalised every 64 steps.
2. **Coverage.** Each fixed cell asks the plates owning it and its four neighbours on the previous step whether they have crust at this location now, by mapping its position into each plate's frame with Qᵀ and looking up the body cell.
3. **Overlaps** (two or more plates cover the cell) are resolved by buoyancy. Continental beats oceanic, thicker continental beats thinner, and younger oceanic beats older. The loser's body cell is marked "subducted this step", with negative thickness so it can be restored (see 3.3). Boundary effects are applied (see 3.4).
4. **Gaps** (no plate covers the cell) are classified, in order:
   - a duplicate sample of crust subducted this step elsewhere: restore it, since it is still at the surface here;
   - a hole inside the departing plate (it has crust on three or more of the body cell's grid neighbours): fill from a neighbour;
   - a true trailing edge: new oceanic crust of age zero, assigned to the plate that just left. This is how ridges form;
   - otherwise, an interior hole surrounded by another plate: extend that plate, copying a neighbour whose crust type matches what was here.
5. **Gather** the owning plates' crust into fixed-grid arrays, run surface processes (3.6), write back (**scatter**), then **sync unsampled cells** (3.3).

### 3.3 The resampling problem and its patches

Nearest-cell lookup between two grids of equal resolution, rotated relative to each other, is not one-to-one. About 18% of fixed cells sample a different body cell each step, some body cells are sampled twice, and some not at all. Most of L0's numerical artifacts trace back to this. The current patches:

- **Restorable subduction:** a cleared body cell can be restored the same step when another fixed cell still sees it. Without this, crust "subducted" by one fixed cell came back as fresh crust, and continents grew to 90% of the planet.
- **Sync unsampled cells:** a body cell no fixed cell sampled this step, surrounded by sampled neighbours, takes their crust type (copying a neighbour's rock column if the type differs) and their mean thickness. Without this, old seafloor survived inside continents that grew around it and blinked into view as single-cell "pinholes".
- **Hole-filling** prefers sources matching the cell's previous crust type, which removed a bias toward continental growth.
- **Welding pulls rather than pushes:** see 3.5.

The proper fix is structural: a body grid finer than the fixed grid, or crust carried on moving points (as in Cortial et al. 2019). See section 7.

### 3.4 Plate boundaries

- **Ocean–continent and ocean–ocean subduction:** the overriding plate thickens 1–3 cells inland along the direction of relative motion (a volcanic arc; weights 0.5, 0.3, 0.2 of 0.6 km per event plus 40% of the subducting crust's sediment). A trench is marked on the down-going side. Thick enough oceanic arcs accrete as continental crust.
- **Continent–continent collision:** all of the consumed crust is kept. Its distribution depends on the across-boundary share f of the relative motion, measured against the boundary normal recorded the previous step. Head-on (f ≈ 1), it spreads over the winner's first four cells inland, imitating underthrusting and crustal stacking. Sideways (f ≈ 0), it piles up at the boundary cell. Continental slices pushed onto seafloor accrete as continent. Collision cells push a collisional-metamorphic rock unit.
- **Welding (sutures):** each pair of plates keeps a collision counter that decays by 3% per Myr. Crossing 30 logs "continents begin to collide"; crossing 180 welds the plates (3.5).
- **Rifting:** roughly every 40 Myr (faster when below the target plate count), a large plate, chosen with probability weighted by area^1.5, splits along a wiggly great circle through its centroid. The halves receive opposite extra rotations about an axis lying in the cutting plane, so they separate perpendicular to the cut at speeds of 1.5–4 cm/yr. Splits with too few cells on either side are reverted.

### 3.5 Welding plates

The smaller plate's crust is transferred into the larger plate's frame. Every cell of the surviving plate's grid asks which absorbed-plate cell lies beneath it. This is a **pull**, and it overwrites stale leftovers on the absorbed plate's territory. The angular velocities are averaged by area. Fault-record strike vectors are rotated into the surviving plate's frame by Q_Wᵀ·Q_L.

### 3.6 Elevation and surface processes

Elevation is derived from crust, never stored:

- Continental: e = 0.16 × (thickness − 33 km), Airy isostasy.
- Oceanic: e = −2.6 − 0.33·√age (km, Myr), from half-space cooling, plus 0.5 × excess volcanic thickness and 0.3 × sediment; floor −7 km.
- Sea level is fixed at 0.

Each step, on the fixed grid:

- **Hotspots** (3–6, fixed in the mantle frame) thicken the cell above them. Moving plates drag curved island chains off them.
- **Diffusion** of continental thickness between continental neighbours (0.03 per Myr per neighbour).
- **Erosion** of emergent continental cells: thickness loss 0.03 × elevation per Myr, which gives exponential decay of excess thickness with a time constant of about 210 Myr. Sediment is spread over **all lower neighbours in proportion to their drop** (multiple-flow routing). On land it pushes a clastic-sediment unit; in the ocean, half becomes shelf sediment.
- **Volcanic decay:** excess oceanic thickness erodes at 2% per Myr.
- **Shelf growth:** ocean cells accumulating enough sediment become continental (a marine-shelf unit on top).
- **Orogenic collapse:** continental crust thicker than 48 km occasionally spreads onto an adjacent seafloor cell of the same plate, creating a 28 km continental cell (basement plus clastic cover).
- **Caps:** continental thickness is held within 10–75 km, oceanic within 5–30 km.

### 3.7 Crust-budget thermostat (a deliberate fudge)

A feedback factor fb = (continental fraction / target)⁶ scales:

- the sediment threshold for shelf growth (by fb);
- the arc-accretion threshold (by fb);
- the probability of orogenic collapse (by 1/fb).

Above target, drowned coastal continental cells (elevation under 0.3 km) founder back to seafloor with probability 0.1 × Δt × (fb − 1). The target comes from the **Continental crust** slider (default 50%), which also sets the starting fraction. This is labelled in the page as a fudge, not physics.

### 3.8 Plate kinematics: speed relaxation and net rotation

- **Random walk with speed relaxation:** each Ω receives Gaussian kicks, and its magnitude relaxes toward the plate's own preferred speed over 40 Myr, within a hard cap of 8 cm/yr and a floor of 0.5 cm/yr. An unconstrained 3D random walk drifts outward in magnitude until it hits the cap.
- **Net-rotation removal (default on):** independent random motions sum to a shared spin of the whole lithosphere that grows over time, like the "flying ice cube" artifact in molecular dynamics. Each step computes ω_net = M⁻¹ Σ_cells (I − x xᵀ) Ω_owner dA, where M = Σ_cells (I − x xᵀ) dA is precomputed from the exact discrete cell areas. Subtracting ω_net from every plate zeroes net rotation to machine precision. Relative motions, and therefore all boundary behaviour, are untouched; only hotspot tracks and continent latitudes change. Earth's lithosphere does drift slowly relative to hotspots (roughly 0.1–0.4° per Myr by different estimates), so turning removal off is not wrong, only unconstrained.

### 3.9 Rock columns

Every crust parcel points into a shared **column pool**. A column is a stack of up to 8 units (bottom first), each with rock type, formation time, thickness, and **origin plate**. Rock types:

1. ocean-ridge basalt
2. hotspot or flood basalt
3. arc volcanics
4. cratonic basement
5. collisional metamorphics
6. marine shelf sediments
7. clastic sediments

Operations:

- **push:** deposition. Merges into the top unit if it has the same rock and origin plate and formed within 150 Myr. Thin veneers under 50 m of a different rock are absorbed first. When the stack is full, the two deepest units fold together into the basement.
- **erode:** shaves the top and pops units, exposing older rock.
- **reconcile:** after each step, any mismatch between stack thickness and crust thickness goes into (or comes out of) the basement unit. This represents thickening or thinning at depth. Mean mismatch is about 0.005 km.
- **clone:** columns are copied whenever crust is copied, and garbage-collected every 32 steps (with their fault records) once no live plate cell points at them.

The **Rock type** map shows the rock with the most thickness in the top kilometre (map generalisation), not the topmost veneer.

### 3.10 Fault records

Each step, every plate-boundary cell splits the relative motion of each differing neighbour into across-boundary slip (convergence positive) and along-boundary slip. A column's optional **fault record** keeps two episodes (previous and current). Each episode stores start, end, across slip, along slip, and a slip-weighted strike vector **in the plate's own frame**, so dead faults ride along with the crust.

- A new episode starts after 10 Myr of inactivity, or after 10 Myr of a different style (thrust, normal or strike-slip, classified by the ratio of along to across slip).
- Sense of strike-slip: right-lateral if the far side moves to the observer's right.

**Extraction** (at handoff) takes continental cells with at least 60 km of total recorded slip and groups neighbours whose strikes agree within about 45°. It then traces each group into polylines by greedy chaining along strike and drops traces shorter than three cells. Each fault gets a kind and sense from its summed slip, an active flag, a start and end time, mean slip, and length.

### 3.11 Formations, domains and the L1 preview

**Formations.** At handoff, units in neighbouring continental columns are linked (union-find) when rock type and origin plate match and ages agree within max(40 Myr, 10% of age). Both stacks are walked bottom to top, so order is preserved. Each linked group is a formation.

**Domains.** Neighbours on the same plate whose stacks share at least 40% of their thickness through linked units join a domain. Domains under 8 columns are absorbed into the neighbouring domain they share the most edges with. Interpolation never crosses a domain boundary, which keeps sutures and terrane edges sharp.

**Band-limited fields.**
- Elevation is smoothed with three passes of 5-point averaging.
- Formation thicknesses are averaged over neighbours in the same domain for two passes. Absent counts as zero, so isolated specks fade and edges pinch out.
- Formations are ordered by thickness-weighted age, with basal formations forced to the bottom.

**L1 preview.** The view interpolates the smoothed fields bilinearly on the cube grid, per pixel. The coast is the zero contour. A ridged noise field stands in for L1 relief; valleys cut deeper into the formation stack and expose older units. The noise is a placeholder, not an L1 model.

### 3.12 Default constants

| Constant | Value | Meaning |
|---|---|---|
| `ISO`, `CONT_BASE` | 0.16, 33 km | continental elevation = ISO × (thickness − CONT_BASE) |
| `OCEAN_T` | 7 km | new oceanic crust thickness |
| `CAP` | 75 km | continental thickness ceiling (excess is deleted; see section 7) |
| `ARC`, `ARCW` | 0.6 km, [0, .5, .3, .2] | arc thickening per subduction event and inland weights |
| `COLLW` | [.35, .3, .2, .15] | head-on collision thickening weights (all consumed crust kept) |
| `KE`, `KD` | 0.03, 0.03 per Myr | erosion and continental diffusion rates |
| `SED_CONV`, `ARC_CONV`, `COLLAPSE` | 5, 17, 48 km | shelf, arc-accretion and collapse thresholds (thermostat-scaled) |
| `HOT` | 0.7 km/Myr | hotspot thickening rate |
| `CONT_TARGET` | 0.40 (UI default 0.50) | continental crust fraction target |
| `MERGE`, `COLL_LOG` | 180, 30 | collision-counter thresholds for welding and logging |
| `SPLIT_TAU` | 40 Myr | mean interval between rifts |
| `VMAX`, `VMIN`, `SPEED_TAU` | 8, 0.5 cm/yr; 40 Myr | plate speed cap, floor and relaxation time |
| `UNIT_MERGE`, `COVER_DEPTH` | 150 Myr, 1 km | rock unit merging window; depth for displayed surface rock |
| `FAULT_GAP`, `FAULT_SWITCH` | 10, 10 Myr | fault episode boundaries |

All are tuned for plausible appearance, not fitted to Earth data.

### 3.13 Tests and measured behaviour

| Test | What it checks |
|---|---|
| `test_ui_smoke.js` | runs the whole page script against a stubbed DOM: every view, both projections, inspector, extraction, export, and an L1 landmass from pick to drainage basins |
| `test_budget.js` | land, continental crust and peak elevation over time for three seeds; fails if land falls under 10% or crust drifts more than 8 points from target |
| `test_netrotation.js` | net rotation with removal on and off; removal must zero it and leave relative motion unchanged |
| `test_handoff.js` | fault, formation and domain statistics; schema consistency of the exported JSON |
| `test_flicker.js` | rock-type flicker per step, split into surface processes versus resampling |
| `test_pinholes.js` | single-cell ocean pinholes inside freshly absorbed plates versus elsewhere |
| `test_volume_budget.js` | continental volume sources and sinks, via a patched copy of the sim |
| `test_l1.js` | L1: stream-power steady state and the viscous sheet against exact solutions; on a real landmass, exact replay, L0 untouched by the sheet, crust budget closure, bounded thickness and peaks, determinism; in the drainage pass, every land cell drains and basins partition the land; prints Hack's law |
| `test_mesh.js` | the marker core: spherical triangulation (build, delete, reinsert, points on edges, flips); per step, the continental volume ledger, mesh validity, no edge between two plates' interior markers, unique ids; area budget; determinism and snapshot replay |
| `bench.js` | ms per step, extraction and export cost at each resolution |

Representative numbers from development, in Node.js on a single core:

- **Cost per step:** about 14 ms at n = 64 and 22–25 ms at n = 96 after warm-up. The first tens of steps are slower; short benchmarks overstate cost. Extraction takes 0.2–0.5 s; the export is 2–7 MB.
- **Budget:** with a 50% target, land levels off at roughly 22–34% after the first few hundred Myr, continental crust holds at 46–50%, and peaks stay at 7–8 km.
- **Plate count and seafloor:** plates stay between 9 and 13; mean seafloor age wanders between about 70 and 220 Myr (Earth: about 60), so seafloor recycles somewhat slowly.
- **Events:** roughly 20–30 rifts and 20–30 sutures per 800 Myr, plus hotspot islands.
- **Rock flicker:** about 5% of continental cells change displayed rock per step, about 90% of it from resampling.

### 3.14 The marker core (Stage 1 of `L0_MESH_SPEC.md`)

An alternative L0 core in `src/mesh.js`, chosen with **Core** before **Rebuild world**. It keeps the lattice core's kinematics and surface rules, but the crust lives on markers instead of plate lattices, so nothing is resampled.

- **Start.** The lattice core builds the starting world; one marker is put in each of its cells (jittered within the cell by a fixed hash), carrying the cell's crust, rock column and fault record. The plates' lattices are then dropped.
- **Mesh.** A Delaunay triangulation of the markers on the sphere (on a sphere it is their convex hull), kept up to date by insertion, deletion and edge flips. Markers inside a rigid plate keep their triangles; only triangles touching two plates change.
- **Interfaces.** Wherever an edge joins two plates' markers, an interface marker goes at its midpoint, with one side per plate (three at a triple junction). Each step it is classed as ridge, convergent or transform from the plates' relative velocity across the boundary. At a convergent boundary it rides the winning plate; each connected stretch of one plate pair's convergent boundary takes its markers' majority verdict on which plate wins (continent, then thicker, then younger: the lattice core's rule).
- **Motion.** Plates rotate rigidly. The step is cut into substeps of at most a quarter spacing of relative motion across any interface; a substep that turns a triangle over is undone and the culprit removed, or halved.
- **Subduction and collision.** A losing marker within 0.3 spacings of a convergent interface marker is consumed, with the lattice core's side effects: an arc on the winner 1–3 spacings inboard and a slab record for an oceanic marker; for a continental one, its whole volume moved onto the winner's markers (head-on: inland; glancing: at the boundary) and the plate pair's collision count raised.
- **Ridges.** Where a plate's nearest marker is more than 1.3 spacings from an interface marker that is not convergent, a new oceanic marker of that plate is placed between them.
- **Area and volume.** Oceanic markers take their area from the geometry each step (their share of the triangles around them and of the interface markers beside them). Continental markers carry area and volume; every change to continental volume is booked by cause, and the ledger closes to about 10⁻⁹ every step.
- **Surface processes.** The lattice core's, per marker, over each marker's neighbours: its mesh neighbours, plus the markers of other plates around each interface marker beside it, so diffusion, erosion and sediment cross plate boundaries. Diffusion uses cotangent weights, which reduce to the lattice's on a square grid.
- **The page** reads a fixed-grid view: each cell shows the nearest marker. **Markers** draws them. L1 refinement still reads the lattice core's plate lattices and asks for it.
- **In a Web Worker.** The page runs the marker core in a worker built from the modules' own sources (`SIMmodule`, `MESHmodule`, `MESH.workerMain`) and keeps a mirror: after each step the worker posts the fixed-grid view (copies, transferred), plates, statistics and new events; the inspector's column, the handoff, its export and the L1 preview fields are asked for and arrive later. If a worker cannot start, the page steps the core itself. In Chromium at n = 96 the page keeps drawing at 60 frames a second while it runs (frames 17 ms apart, against 350 ms when the core steps on the page), but the simulation is no faster (7.6 against 7.3 Myr a second; the lattice core runs 17.8): one worker still does the whole step. The worker's world matches a headless run exactly.
- **A second worker for the view.** The page has no shared memory in claude.ai (`crossOriginIsolated` is false and `SharedArrayBuffer` is undefined; the statistics row shows this), so workers only exchange copies. Of the step, only the fixed-grid view is separable: the dynamics read the markers, not the view (the hotspot island check now reads the marker its cell shows). So the step worker hands a copy of the markers, their neighbour lists and the rock columns to a helper worker and starts the next step; the view comes back and goes to the page with that step's summary. Measured in Chromium at n = 96, 25 s each, alternating: 7.5 and 8.1 Myr a second with the helper, 7.2 and 7.3 without: about 8% faster, and noisy. The rest of the step is a chain of operations on one triangulation, so more workers would not help without shared memory.

**Measured** (`tools/compare_cores.js`, three seeds, n = 96, 800 Myr, against the lattice core on the same seeds; details in `L0_MESH_SPEC.md` §11.1 and `DEVLOG.md`):

- The invariants hold: volume ledger to about 10⁻⁹, mesh valid at every check, ids unique, runs deterministic and replayable from a snapshot (`tests/test_mesh.js`).
- Resampling flicker is gone from the material: displayed rock changes on 0.55% of continental markers a step, the lattice core's surface-process share.
- Continental crust settles about 4 points lower (0.35–0.36 against 0.39–0.40) and seafloor is older (236–268 Myr against 124–206 at 800 Myr). The lattice core counts about 3.5 times as many oceanic contacts as the marker core consumes, and each deposits arc crust; seafloor turnover is the same on both, so the lattice's extra arc crust comes from resampling.
- Land fraction, plate count, rifts and sutures stay in the lattice core's ranges, which vary more between seeds than between cores.
- It cost 100–126 ms a step at n = 96 against 29–34 ms for the comparison; caching what rigid plate interiors cannot change has since brought it to about 103 ms headless. In Chromium, in its worker with the view helper, it simulates 9.6 Myr a second against the lattice core's 18.3.

---

## 4. The L1 handoff format

`Download JSON` writes one object:

```text
meta        seed, settings, snapshot time, units, grid definition (cell index = face*n*n + j*n + i),
            rock type names, fault kind names, notes
plates      [{id, eulerPoleLat, eulerPoleLon, rateDegPerMyr, maxSpeedCmPerYr,
              areaFraction, continentalFraction}]
events      [{ma, type, text, lat, lon}]   types: rift, collide, suture, lost, island
faults      [{id, kind, sense, active, startMa, endMa, slipAcrossKm, slipAlongKm,
              lengthKm, trace: [[lat, lon], ...]}]
formations  [{id, rock, formedFromMa, formedToMa, originPlate, columns, meanThicknessKm}]
domains     number of domains
cells       per-cell arrays: lat, lon, plate, elevationKm, smoothedElevationKm, continental,
            domain, fault (id or -1), columnOffsets (N+1),
            columnUnits (flat, 4 per unit: rock, formed_ma, thickness_km, formation_id)
```

Times are "Myr before the snapshot". Velocity anywhere is Ω × position, from each plate's Euler pole and rate.

**How L1 should use it:**
- Interpolate `smoothedElevationKm`, never the raw elevation, then threshold the coast at L1's resolution.
- Interpolate each formation's thickness within its domain only.
- Treat fault traces as fault *zones* to be refined into strands, not exact lines.

---

## 5. L1: landmass

L1 takes **one landmass at the present**: a single piece of land bounded by its coastline at the end of L1. It replays the last 10–30 Myr of that landmass's history (default 20) at finer resolution, with the crust deforming instead of riding rigid plates. Then it settles the landmass's drainage basins, which become the tiles for L2.

The window is short because of the time scales involved:
- Rivers respond to a change in uplift in about 1–10 Myr, so the drainage network matures almost at once and needs no long run.
- Collision zones build mountain belts and plateaus over 10–50 Myr. That is what L0 cannot resolve: its plates are rigid, so a collision becomes thickening of a few 100 km cells.

So the window is spent on crumple zones. Tectonics then stops, and a short erosion pass on a fine grid draws the rivers.

An earlier version of L1 (DEVLOG, "L1 v1") ran one continent for ~100 Myr on a refined cube-sphere with fine stream-power erosion all the way. It handled meeting continents with a world-space overlap layer. That is replaced by the design below.

### 5.1 Picking at the end, with replay

- **Snapshots.** While L0 runs in the page, it saves its whole state every 10 Myr and keeps the last four. A saved state includes the random generator's state. It costs 14 MB at n = 64, 29 MB at n = 96 and 59 MB at n = 128, and about 50 ms to take at n = 96.
- **Pick.** At the present, a tap on land selects its landmass, grown from the tapped cell (`landmassAt(world, cell, bufferKm)`):
  1. *Land:* the land cells connected to it, corners included (two cells meeting at a corner are one coast at L0's resolution).
  2. *Buffer:* the sea within `bufferKm` of that land's coast (great-circle distance between cell centres, less one cell). An island with a cell in the buffer joins, and the buffer grows around it, until nothing more joins. This is single-linkage clustering of the land at the buffer distance: an archipelago whose gaps are all under the buffer comes as one landmass. The buffer always includes the ring of continental crust around the land, so the coast can move. Default 200 km; 0 to 1,000 km on the slider.
  3. *Enclosed seas:* water the region surrounds (every piece of the rest of the sphere except the largest, the open ocean) joins, with any islands in it. A strait narrower than about twice the buffer is closed by the buffer, so a sea behind one counts as enclosed (in image-processing terms, a morphological closing).
  4. *Crust:* of the region's continental crust, only what connects to the landmass's land through the region's continental crust is kept, so a sliver of another continent's shelf that the buffer reaches stays out.

  The sheet carries:
  - the region's continental crust, including crust that was still oceanic when the window began (accreted during it);
  - the landmass's land on oceanic crust (islands built by arcs and hotspots, 1–3% of a large landmass's land);
  - the *coastal rim*: oceanic cells touching the land, corners included, on the land cell's own plate. L0's coast is where its bilinear surface crosses zero, up to half a cell out into such a cell. Without the rim, the edge of the sheet's crust cut the coast short in straight lines. Ocean on another plate (across a trench) stays out, so the sheet doesn't drive the subducting plate into the coast.

  The rest of the region's oceanic floor is not modelled at L1. Once refined, the L1 map draws it as plain ocean, without L0's surface: L1 moves the coast (the sheet deforms where L0's plates are rigid, and adds the mantle pull), and two shorelines side by side were confusing. The plate carrying most of the land is the reference plate.

  Measured before the coastal rim (seed b, 1 Gyr, 96-cell grid, 10 Myr window, largest landmass): with a 300 km buffer, 15 islands joined and 4,873 of the 4,889 land cells lie under the sheet's crust at the end. The other 16 cells do carry markers, but the sheet moved the coast a cell away from where L0 has it. Land in the L1 frame without sheet crust is now the other landmasses, beyond the buffer.
- **Pass 1, tracing.** L0 restores the newest saved state older than the window and tags every crust parcel with an identity (plate × N + body cell). L0 carries that identity wherever it moves crust, through splits, welds and hole fills. It then replays to the present, recording which parcel sits on each fixed cell at each step of the window. L0 is deterministic, so the replay lands exactly on the present; this is tested.

  The landmass's crust is two sets of parcels:
  - every parcel that ends up in the landmass;
  - every continental parcel that touches that crust at some step and is gone by the end. These were overrun in a collision, or relabelled by L0's resampling. That crust was squeezed into the landmass and is what its crumple zones are made of.

  A continental parcel with at least three traced neighbours also counts as traced, which fills pinholes from L0's resampling.
- **Map frame.** The map is an azimuthal equidistant projection about the landmass's centroid, co-moving with the reference plate. If the reference plate split off another plate during the window, the frame follows its ancestor before the split. The footprint is everywhere the crust has been during the window, plus a 400 km margin.
- **Pass 2, the sheet.** L0 restores the same state again, replays to the window's start, and from there runs alongside the sheet, one L0 step per sheet step. It ends exactly where it started.

### 5.2 The thin viscous sheet

The crust is treated as a thin viscous sheet, the vertically averaged model of continental deformation of England and McKenzie (1982):

∂ⱼσᵢⱼ − κ (uᵢ − Uᵢ) = ∂ᵢΦ,  with σᵢⱼ = 2η (ε̇ᵢⱼ + δᵢⱼ ε̇ₖₖ) and Φ = ½ (S² − S₀²)

- **Symbols.** u is the sheet's horizontal velocity and ε̇ its strain rate. The δᵢⱼ ε̇ₖₖ term is the vertical strain, from incompressibility. S is crust thickness and S₀ the thickness at the window's start.
- **Buoyancy.** The right-hand side is the gradient of the crust's excess gravitational potential energy (constants absorbed into η and κ). Thickened crust pushes outward, so plateaus spread.
- **Basal drag.** κ (u − U) drags the sheet toward the velocity U of the plate that L0 says carries each marker's parcel at that step, relative to the reference plate. This is how L0's plate motions drive the sheet: where two plates converge, the crust between them shortens and thickens.
- **Viscosity.** Power law with exponent 3: η = η₀ · strength · (ε̇/ε̇₀)^(1/3 − 1), floored at 2% of ε̇₀. Old cratons (basement older than 500 Myr) are ten times stronger. The craton share is smoothed on L0's grid (three passes, continental crust with continental crust) and sampled bilinearly at each marker. Taken parcel by parcel, its edges were L0's cell edges, and the weak crust between cratons crumpled along those staircases into square ridges.
- **Parameters.** Two dimensionless numbers set the behaviour:
  - the **Argand number** Ar (default 3) sets η₀ at the reference strain rate ε̇₀ = 50 km/Myr ÷ ℓ. Higher means weaker crust that spreads into broader plateaus.
  - the **coupling length** ℓ = √(4η₀/κ) (default 500 km) sets how far a plate's drag reaches into the sheet.
- **Numerics.**
  - Staggered (Marker-and-Cell, MAC) grid, 50 km cells by default.
  - The operator is the gradient of a discrete dissipation function, so it is symmetric positive definite and conjugate gradients apply (Jacobi preconditioner, relative tolerance 10⁻⁴, at most 200 iterations), with two Picard iterations on the viscosity per solve.
  - Edges of the crust are traction-free; cells hold crust where the markers cover at least a quarter of them.
- **Markers.** There are 6 × 6 markers per L0 parcel, one at a random point (a fixed hash) in each of 6 × 6 sub-cells, so marker rows don't alias with the grid or bunch into lines under shortening. Each carries its parcel's identity, plate, rock column (as fractions of thickness), strength, the thickness and elevation it started with (L0's elevation smoothed over three passes, land only with land and sea only with sea so the coast stays where L0 has it, sampled bilinearly), and its own thickness and area.
  - Markers move with the sheet (midpoint rule, substeps limited to half a cell per substep).
  - They thicken by the local convergence: S → S·e^(−div·dt), A → A·e^(div·dt). Each marker's volume S·A is therefore exactly conserved.
  - Grid fields are bilinear weighted averages of marker values.
- **Finite strain.** Each marker carries a 2×2 deformation gradient F, updated every substep as F ← exp(L·dt)·F. L is the sheet's velocity gradient at the marker's midpoint: du/dx and dv/dy at cell centres, du/dy and dv/dx at cell corners. Its trace is the same divergence that thickens the marker, so det F equals the marker's area change exactly. The map's finite-strain view shows, per cell, the principal stretches of the marker-averaged B = F·Fᵀ: total horizontal shortening, its direction (ticks) and vertical-axis rotation. This is what L2 needs to build thrust sheets and folds. Under shortening with free sides, crust also stretches at right angles (squeezed out sideways), so the view colours shortening first and shows extension only where there is next to no shortening.
- **Dynamic topography.** Slabs sinking into the mantle pull the surface above them down.
  - **Recording (L0, recording only).** L0 records the oceanic area consumed at each trench. Each area is weighted by the slab's age (older is denser: √(age / 80 Myr), at most 1) and deposited about 300 km beneath the overriding plate. The record goes on its fixed grid, which is the mantle frame, in 10 Myr bins, keeping the last 250 Myr.
  - **Sinking.** A bin of age a has sunk to depth d(a): 20 km/Myr to 660 km, then 10 km/Myr.
  - **Deflection.** Each bin's pull spreads as a Gaussian of width 300 km + 0.8·d and weakens as e^(−d/1200 km). The deflection is −0.6 km times the summed, spread slab density (consumed area per area), less its mean over the sphere: return flow lifts the surface away from slabs.
  - **Computation.** On a cube grid of half L0's resolution (at least 48 cells a face, about 210 km), as lazy diffusion steps applied oldest bin first, so each bin gets exactly its own width (Gaussian widths add in quadrature). It is recomputed every sheet step (about 90 ms) and sampled under each cell as the map moves over the mantle. It is added to the sheet's elevation, so it also enters the present rates the drainage pass uses.
  - **Status.** A calibrated stand-in, not mantle flow. Over the continents of the 1 Gyr test world it ranges from −1.0 km (1st percentile) to +0.37 km (90th), the range estimated for Earth.
- **Volcanism.** L0 records the crust its arcs and hotspots add. Its collisions and collapse are not passed on, because the sheet makes those itself. Each marker thickens by what L0 gained there (interpolated bilinearly between L0's cells), as a new top unit of arc rock.
- **Elevation.** As in the first version: the inherited elevation plus the thickness change, minus a flexural root, e = e₀ + ΔS − 0.84·F[ΔS]. F is a Gaussian of width 1.18 × the 40 km flexural length, and the plate yields beyond 3 km from local balance. The inherited state counts as equilibrium.
- **Erosion during the sheet** is a smooth stand-in for rivers, which a 50 km grid cannot draw. Denudation is 0.5 km/Myr × (e / 5 km)^1.4: 0.5 km/Myr in ranges at 5 km and about 20 m/Myr on lowlands at 0.5 km. These are the magnitudes of Ahnert's (1970) relation between denudation and relief. Each marker loses the depth interpolated bilinearly to where it is. Eroded rock is stripped from the top of the rock columns and leaves the system. A first version ran stream power on the sheet grid; its 50 km gullies were aligned with the grid, and the drainage pass inherited them as parallel stripe-shaped basins.
- **Budget.** A volume budget tracks initial crust, arcs and hotspots, erosion and crust carried off the grid. It closes to about 10⁻¹⁴ of the initial volume; this is tested.

### 5.3 The drainage pass

With tectonics frozen at the end, the sheet's surface is refined to 10 km cells: bilinear, plus fixed fractal relief of ±(0.15 + 0.15·e) km in the band the sheet cannot resolve (wavelengths from 150 km down), seeded from the world seed. Then 5 Myr of stream-power erosion runs (E = K A^m S, m = 0.5, n = 1, K = 0.1, implicit solver of Braun and Willett 2013, Priority-Flood+ε depression filling of Barnes et al. 2014 in the Zhou et al. 2016 variant, D8 receivers). That is long enough for the network to organise, but not to change the landscape.

As in L1's first version, each cell's node (used for distances, hence slopes) is displaced within the cell by up to 0.4 cell widths by a fixed hash. Neighbour directions then vary from cell to cell, so rivers are not locked to the grid's 0° and 45° directions.

**Climate.** A deliberately crude present-day climate, enough to tell wet from dry:
- **Rain by latitude.** Away from relief, rain follows latitude as on Earth: wet equator, dry subtropics, wet westerly belts, dry poles.
- **Winds.** They are zonal: easterly trades below 30°, westerlies to 60°, polar easterlies beyond. The belts blend over 15° at their edges.
- **Moisture.** It blows in from the sea and rains out with distance (e-folding over 3,000 km; long because vegetation recycles rain) and with uplift (e-folding per 3 km climbed). Rain is proportional to the moisture lost, at most 5× the latitude value. Windward slopes are therefore wet and lee sides lie in a rain shadow. Uplift is measured on the surface smoothed over 30 km, since rain responds to relief tens of km across, and the rain field is smoothed over the same width.
- **Evaporation.** Potential evaporation peaks in the clear-sky subtropics (2 m/yr at sea level, 1.4 at the equator) and falls with elevation.
- **Runoff.** Land loses water by Fu's form of Budyko's curve, a standard relation between rain, potential evaporation and actual evaporation.
- **Effect on erosion.** Runoff weights the discharge that drives erosion (normalised by its mean over the land), so dry regions erode less.

**Present tectonic rates.** The sheet keeps a running average, over the last 3 Myr, of its rate of elevation change from tectonics alone: deformation, volcanism and their isostasy, not erosion. The drainage pass keeps applying it, so the pass approximates the next few Myr at present rates rather than a frozen landscape.

**Lakes, sediment and closed basins.**
- **Sediment.** Lakes trap sediment. A depression that its catchment's sediment can fill within the pass becomes a flat that water crosses. Sediment is river erosion plus background denudation, at 20 m/Myr × the Langbein–Schumm curve: sediment yield peaks in semi-arid climates at about 0.3 m/yr of rain and falls in deserts and under forest. Hollows in the fine relief fill this way. Big depressions, and those subsiding faster than they fill, survive.
- **Water balance.** A surviving lake overflows only if its inflow beats evaporation minus rain over its full surface. Otherwise its basin ends there, and the lake shrinks to the area its inflow can supply, leaving a dry lake floor.
- **Cut-off seas.** Sea cut off from the open ocean is inland water under the same rules.

**Rivers that run dry.** In dry country, flood pulses soak into channel beds and floodplains. Rivers lose about the same discharge per km whatever their size: 0.3 m/yr·km² per km (about 0.01 m³/s per km) at the driest, falling to nothing at an aridity index (rain ÷ potential evaporation) of 0.65. So small desert rivers die within a few hundred km, while a Nile-sized one crosses 2,000 km of desert losing about 1%. Where a river of at least 1 m³/s mean flow runs dry, its basin ends at a terminal fan.

**Inland drainage.** Land drains inland if its basin ends in a lake with no outlet or where a river runs dry (endorheic). It is also inland if its water never reaches the coast because the wash it feeds dies on the way (areic: desert with no flowing water).

**Basins.** Every land cell is then labelled by its outlet: a point on the coast, a lake with no outlet, or the place a river runs dry. The map colours basins of at least 10,000 km² (inland-draining ones in muted earth tones, areic desert sand-coloured); smaller coastal basins are grey. Lakes are blue if they overflow, pale blue if they have no outlet, and white for dry lake floors.

### 5.4 Measured behaviour

Measured on seed "genesis" in Node.js, unless stated:

- **Sheet solver against exact solutions.** A Newtonian sheet with no buoyancy has its left half dragged at +U and its right half at −U. The exact profile is U(1 − e^(−|x|/λ)). The solver matches it to 0.06% in a narrow strip with free sides (λ = √(3/4)·ℓ) and 0.03% in the middle of a wide sheet (λ = ℓ).
- **Finite strain** is exact to rounding (about 10⁻¹⁵) for uniform pure shear, simple shear and rigid rotation, the three branches of the 2×2 matrix exponential. On a real landmass det F matches each marker's area change to 3 × 10⁻¹⁵. A 110 M km² supercontinent's collision belt reached 84% horizontal shortening in 20 Myr.
- **Stream-power solver.** It reaches the analytic steady state S·A^m = U/K at every cell to 5 × 10⁻¹⁵.
- **How much the simulation determines drainage.** Re-drawing only the synthetic fine relief (two other seeds) moves each land cell's outlet a median of 117 km on the 1 Gyr landmass and 125–130 km on the 300 Myr supercontinent. About 70% of land keeps its outlet within 300 km. So L1 determines which stretch of coast water reaches; only the exact divides between neighbouring rivers are contingent. Label agreement (adjusted Rand index of 0.36 and 0.52 for major basins) penalises every shifted divide and understates this.
- **Dynamic topography** (−0.45 to +0.39 km on the 1 Gyr landmass, −0.63 to +0.19 km on the supercontinent, deepest along the supercontinent's inner sea above the subducted ocean) reorganises basins toward its depressions. It raised the 1 Gyr landmass's mean elevation from 1.37 to 1.49 km and its land from 36.3 to 37.8 M km². I expected it to make major basins much less contingent. It does so only slightly: the median outlet shift falls from 117 to 111 km and from about 127 to 109 km, and label agreement is unchanged. The prediction largely failed.
- **Inland drainage, tested** on a round landmass 1 km high with a broad bowl 600 m deep:
  - At 5° latitude its lake overflows to the sea.
  - At 24°, standing still, sediment fills the bowl within the pass.
  - At 24°, subsiding 0.1 km/Myr at its centre, it stays a lake with no outlet.
  - A dry plain at 24° has rivers running dry before the coast, and 58% of it drains inland.
  - A sea cut off from the ocean overflows at 5° and has no outlet at 24°. The share draining inland by aridity index is the same on two real landmasses: 88–89% of arid land (below 0.2), 19–25% of semi-arid land (0.2–0.5), 1% of dry subhumid land, and none of the humid. The totals then follow from how much land sits in dry belts: 14% for the 1 Gyr landmass (19% in the browser's world at 1 Gyr), 5% for the mostly humid 300 Myr supercontinent; Earth drains about 20% inland. Present tectonic rates decide where lakes persist rather than how much land drains inland. With them switched off, the supercontinent keeps 1 lake with no outlet instead of 22, with a similar inland share (5.0% against 4.7%).
- **Hack's law** in the drainage pass comes out at 0.52–0.55, with no tuning (Earth's rivers: roughly 0.5–0.6).
- **Lattice alignment in the drainage pass** (landmass below, channels draining 5,000 km² or more): node jitter cuts the share of channel steps in straight runs of four or more identical steps from 52% to 33%. The share of 6-step reaches within 5° of a grid direction falls from 25% to 18% (uniform would be 22%).
- **A landmass at 1 Gyr** (medium grid; 41.5 M km² of land on plate 44, three lobes, two plates at the window's start): pass 1 replays 27 Myr in 2.1 s and traces 4,224 parcels. The sheet is 223 × 185 cells of 50 km with 152,928 markers, at 1.65 s per L0 step (22 steps). The replay and the sheet both leave L0 bit-identical. After 20 Myr: the most thickening is 13 km, in parallel bands at two collision fronts; the highest point is 6.9 km; the strain rate falls from about 0.3/Myr during the collision to 0.02/Myr once L0 has welded the plates. Erosion removed 71 M km³ (1.7 km averaged over the land), arcs and hotspots added 1.5 M km³, and the budget closes to 9 × 10⁻¹⁶. The drainage pass (1115 × 925 cells) finds 5,203 basins, of which about 380 are 10,000 km² or larger; the largest is 1.4 M km². It takes 2.2 s.
- **Solver tolerance.** On the same landmass, capping conjugate gradients at 200 iterations leaves a residual of about 5 × 10⁻³ and costs 1.66 s per step. A fully converged solve (10⁻⁶, four Picard iterations) costs 14.8 s per step. The final surfaces differ by 8 m on average (at most 0.35 km), against a mean thickness change of 1.5 km. A 30 km grid (5.4 s per step, never converging within the cap) gives the same peaks to within 0.15 km and the same most-thickening (12.5 km). The 50 km default is a cost choice, and the answer is not sensitive to it at that level.
- **The thickening bands are contacts, not grid.** At the strongest collision front, two bands about 300 km apart appear at sheet grids of 70, 50 and 30 km, and with or without marker jitter. What follows the grid is their width, about 1.5 cells each time: the shortening localises more narrowly than the grid resolves, as power-law crust tends to. Tracing the markers back shows the bands are two separate contacts. One is the suture between the two plates' crust; the other is a promontory of one margin that hit first. Their spacing comes from the shape of the coastlines L0 brought together.
- **A thickness spike** in an early run (one marker at 431 km for one step) came from passing on L0's orogenic collapse as volcanism: L0 turns an ocean cell into 28 km of continent in one step, and that volume landed on the few markers in the cell. Only arcs and hotspots are passed on now, as thickness per marker.
- **Cost in Chromium:** a 110 M km² supercontinent (365 × 366 sheet cells, 381k markers) took 5.9 s per step, about 2.5 minutes in all, with the drainage pass on 1825 × 1830 cells taking 8.8 s. Small landmasses take seconds.

### 5.5 The L2 handoff format

`L1.handoff(sheet, basins)` cuts the landmass into tiles and packages what L2 needs to refine each tile on its own. The page's **L2 tiles** view shows the tiles and their divides; tapping one shows its summary, and two buttons download the landmass summary or that tile's fields as JSON.

**Tiles.**
- Each drainage basin of at least 10,000 km² is a tile.
- Smaller basins (short coastal streams, small inland sinks) are grouped with their neighbours, breadth first, into tiles of about that size. A group boxed in below half that size joins the neighbouring tile it shares the longest boundary with. Isolated small islands stay as they are.
- A union of whole basins is closed to water, so a river leaves a tile only at one of its outlets. This is tested.
- Kinds: *drains to the sea*, *drains inland* (its largest outlet is a lake with no outlet or a river that runs dry), *coastal strip* and *inland group* (grouped small basins).

**Divides: firm or contingent.** A tile edge is *firm* where its local relief exceeds the peak-to-peak amplitude of the synthetic fine relief, 2 × 0.15 × (1 + e) km. Local relief here is the height above the lowest land within one sheet cell. Elsewhere the edge is *contingent*, and L2 may move it. 18% of the length of divides between tiles is firm on the 1 Gyr landmass.

How firm a divide is turns out to be graded, not sharp.
- **Test.** On the 1 Gyr landmass, the synthetic relief is re-seeded (two seeds). The question is whether a divide still lies between two points 30 km either side of each original divide. The baseline is random pairs of points as far apart, which straddle some divide anyway, more often in rugged country.

| Local relief (km) | Divides still there | Random pairs |
|---|---|---|
| 0–0.2 | 45% | 25% |
| 0.2–0.4 | 54% | 38% |
| 0.4–0.8 | 67% | 45% |
| 0.8–1.6 | 76% | 50% |
| over 1.6 | 90% | 62% |

- **Result.** Firm divides persist 77% of the time, contingent ones 51–53%. The table travels with the handoff (`meta.dividePersistence`), so L2 can weigh it.
- **Rejected predictors.**
  - Crest height: high plateaus' divides were the least stable.
  - The convexity of the sheet's smooth surface: weaker.
  - Divides of the drainage routed over the smooth surface alone: worse than chance, because its near-flat lowlands are full of arbitrary divides.
- **Escarpments.** Steep coastal escarpments count as relief, so the short divides between their parallel streams come out firm.

**Fields at the resolution they were computed.**

| Grid | Per cell |
|---|---|
| Fine (the drainage pass, 10 km) | `elevationM`; `receiver` (the steepest-descent direction; −1 marks a basin end); `areaKm2`; `dischargeM3s`; `lake` and `drainage` codes; `divide` (0 none, 1 contingent, 2 firm); `rainMmYr`; `evaporationMmYr`; `upliftMmYr` (present tectonic rate); `responseMyr` (below); `inTile` |
| Coarse (the sheet, 50 km) | `thicknessKm`; `thickeningKm` in the window; finite strain (`shortening`, `stretch`, `shorteningAzimuthDeg` clockwise from north, `rotationDeg`); `strainRatePerMyr` now; `mantlePullKm`; `denudationKmPerMyr` (present); `arcAddedKm` in the window; `stacks` |

`stacks` lists the L0 parcels present in each coarse cell, up to three, largest first. Each entry gives the parcel's share of the cell's crust, its thickness, and its units bottom to top as [rock, formed Ma, thickness km]. Parcels lie side by side, because the sheet does not stack crust vertically.

**Imprints, oldest first.** Each imprint has a class, a start and end age (Ma, before the present), whether it is ongoing, what it leaves (structure, rock, relief, water) and the tiles it touches. L2 applies them in order, so younger structures cut older ones: Steno's and Hutton's cross-cutting principle, which is how geologists read the order of events from outcrops.
- **Classes:**
  - *fault* (L0's fault records: kind, slip, a trace in map km);
  - *crumpling* (the sheet's finite strain over the window);
  - *volcanism* (arc and hotspot rock added in the window);
  - *mantle pull* (dynamic topography from slabs subducted over the last 250 Myr);
  - *drainage* (the network the pass settled).
- **Structure versus relief:** structural imprints (rock juxtaposed at a fault, strain fabric) last until erosion strips the rock carrying them. Topographic ones decay.

**Lifetimes of topographic imprints, by erasure mechanism:**
- **Hillslope blurring.** Width σ = √(2Dt), with D = 0.01 km²/Myr, which equals 0.01 m²/yr. This is a literature value, not simulated: soil-mantled slopes range about 0.001–0.05 m²/yr. A 100 m wide scarp lasts about 0.5 Myr, a 1 km wide feature about 50 Myr.
- **Knickpoint sweep.** A change at an outlet (a fall in base level, a capture, a scarp across a river) travels up the network at K·A_w^m, where A_w is runoff-weighted drainage area. Celerities are 1 km/Myr at 100 km², 10 at 10,000 km² and 100 at 1 M km². A cell's response time, `responseMyr`, is the travel time from its outlet; a base-level change older than that has passed the cell. K = 0.1 is the drainage pass's own value, not calibrated to absolute rates, so these times are relative.
- **Relief decay.** Half-life under L1's own denudation law:
  - Wide relief, isostatically compensated: 95 Myr from 1 km, 61 Myr from 3 km.
  - Narrow relief, flexurally supported: 15 and 10 Myr.
- **Lake filling.** Each lake's volume over the sediment reaching it gives `fillMyr`.
- **Stripping.** A surface unit lasts its thickness over `denudationKmPerMyr`.

**Frame.** Map km, x east and y north at the projection centre, which is given in L0's present latitude and longitude. Outlets and tile centres also carry latitude and longitude.

**Measured on the 1 Gyr landmass:**
- **Tiles:** 646 in all, with a median of 15,800 km² and the largest 1.26 M km². 353 drain to the sea, 14 drain inland, 266 are coastal strips and 13 inland groups.
- **Imprints:** 57, of which 53 are L0 faults. The oldest started at the start of L0's run.
- **Response time:** a median of 54 Myr across the land; 90% of cells are under 206 Myr. Headwaters in desert, with almost no runoff, reach over 1,000 Myr.
- **Sizes and cost:** the summary is 0.83 MB of JSON and the largest tile 1.7 MB. Building it all takes 0.8 s.

### 5.6 Selecting basins for L2

L2 refines drainage basins, not a landmass. A tap on the refined L1 map (in any view but L2 tiles) selects them about the tapped point, with the **Selection radius** slider in the L2 panel (default 150 km, 0 to 1,000 km) (`basinSelection(B, proj, x, y, radiusKm)`):
1. *The tapped basin*, whole, whatever its size. A tap at sea selects none.
2. *Basins wholly within the radius*: every cell centre within the radius of the tap (great-circle distance). This picks up a cluster of small coastal basins without the large basin behind them, which a buffer around the tapped basin would take.
3. *Enclosed basins*: land the rest of the map cannot reach without crossing the selection (8-connected, from the grid's edge and the open sea), as enclosed seas join a landmass at L1.

The page draws the circle and fades the land outside the selection, and the panel reports the basins, their area and their cell count at L2's planned 500 m. The tap and radius go into the session log, and a rerun checks the selection (tapped basin, basins, cells). Measured on the genesis landmass at 1 Gyr: a tap on the east coast with a 300 km radius selects 54 basins, 54,200 km² (the tapped one 100 km², 53 wholly within the radius); a tap inland selects its basin of 605,000 km² alone. A selection takes about 20–60 ms.

### 5.7 What it does not do yet

- **Coupling to L0 is one-way.** L0 does not see the sheet's deformation.
- **Sutures stop the crumpling.** When L0 welds two plates, their markers get the same velocity, so a collision crumples the crust only until L0 welds it. That is L0's rule, and L1 inherits it.
- **Crust is seeded from what exists at the window's start.** Every parcel that ends up in the landmass is seeded, including continental crust that was still oceanic then (accreted during the window). Crust that is new during the window comes only from L0's arcs and hotspots, added to markers already there.
- **Most of the ocean is outside the sheet.** The sheet carries continental crust, land on oceanic crust and the coastal rim of oceanic cells on the land's own plate (section 5.1). The rest of the region's ocean floor is not modelled at L1, and the refined map draws it as plain ocean.
- **Erosion during the sheet** is a smooth function of elevation, and eroded rock leaves the system instead of filling basins.
- **The drainage pass** has uniform erodibility (by design: section 6, item 3).
- **Fault traces in the handoff** are L0's, at 100 km cells, where L0's rigid plates put them. In crumple zones the sheet has moved the crust by up to the shortening.
- **L0's event catalogue is not in the handoff.** Its locations are where events happened at the time, and the continents have moved since.
- **Cost.** The conjugate-gradient solve is the cost: a multigrid preconditioner, or a Web Worker, would make large landmasses quicker.

## 6. Sketches of the levels below

These are design sketches, not implementations (L1 above is partly built; its sketch lists what remains). Each level should follow the same pattern as L0: deterministic given the seed and its parent's output; owns one band of wavelengths; hands down band-limited fields plus history.

### L1: landmass: what remains

Built so far: section 5. The current plan, organised around the L1→L2 handoff (each item either feeds it, decides whether something belongs in it, or is parked):

**A0. Inland drainage (next)**

Roughly a fifth of Earth's land drains inland, much of it without any depression. L1 currently produces very little.
- ~~**a. Rivers that dry up.**~~ Done (section 5.3). Channels in dry country lose water to evaporation and infiltration, at evaporation minus rain over the channel's area. Channel width scales with the square root of discharge. Where discharge reaches zero the basin ends, at a terminal fan or inland delta.
- ~~**b. Enclosed seas.**~~ Done. Water below sea level that is cut off from the open ocean (a Caspian), or a depression below sea level (a Dead Sea), is a lake with a water balance, not an ocean outlet.
- ~~**c. Present rates in the drainage pass.**~~ Done. Instead of freezing tectonics, the pass applies the sheet's present-day rate of elevation change. A depression stays closed if it subsides faster than its catchment's sediment fills it, and if its rim rises faster than the outlet river can cut it down.
- ~~**d. Mantle-driven topography from L0's subduction history.**~~ Done (section 5.2); see 5.4 for the result. A separate piece of L0 and L1 work. Slabs sinking under continents pull the surface down by up to about 1 km over 1,000–3,000 km: sag basins and continent-wide tilts. L0 records subducted volume per mantle-frame location, as recording only, and L1 applies the smooth field. It should also make major basins depend more on simulated history and less on unresolved relief, which the relief re-seed test can measure.

Measure each step by the inland share in dry belts, the major-basin re-seed score, and the bowl and enclosed-sea tests.

**A. L1 work, in order**
1. ~~**Finite strain per marker.**~~ Done (section 5.2): total horizontal shortening, its direction and vertical-axis rotation, with a map view.
2. ~~**Closed basins.**~~ Done (section 5.3). A water balance on the depressions the drainage pass finds, with a crude aridity field (latitude belts plus a rain shadow from the sheet's topography and a prevailing wind): a hollow overflows only if inflow beats evaporation over its lake. This decides basins that drain inland versus to the sea.
3. ~~**Erodibility test (a decision gate).**~~ Done; **detailed erodibility goes to L2.** `tests/exp_erodibility.js` compares L1's basin map under uniform K with K from:
   - rock type (×0.9–2.2 between the 5th and 95th percentiles);
   - type plus consolidation with age and fault damage from L0's fault records (×0.5–3.5, 40× overall);
   - the same squared (×0.25–12, 1,600× overall).

   On the 1 Gyr test landmass, the adjusted Rand index against the uniform-K map (1 = identical) is 0.98, 0.92 and 0.89; switching off the climate weighting gives 0.94. Re-drawing only the synthetic fine relief with another seed gives 0.35 (two seeds), and the same within the major basins.

   So L1's lowland drainage divides are set mainly by relief at scales the sheet cannot resolve, not by erodibility. That relief is not too loud: its 53 m RMS on lowlands continues the sheet's own spectrum (75 m at 100–200 km, 109 m at 200–400 km). Lowland divides are genuinely contingent; coasts and range crests are what L1 determines. The `Kmul` and `salt` options of the drainage pass remain for such tests.
4. **Glacial overprint (optional).** Gated by a world setting: whether there was recent ice, and how long ago it left. Only what is regional and permanent: the glaciated-region mask with its age, river diversions, large troughs and lake basins. It shares the lake logic of item 2 and gets the same basin-map test. No glacial buzzsaw: it is indistinguishable from older mountains.

**B. The L1→L2 handoff**

5. ~~**Schema.**~~ Done (section 5.5), without the glaciation record (item 4 is not built). Passes are ordered by the age of the event they represent, oldest first. Every imprint carries its event class and age, plus a table of lifetimes by erasure mechanism (hillslope blurring, σ = √(2Dt); knickpoints sweeping up the network; lakes filling). Each basin tile carries:
   - its outlet and drainage graph;
   - whether it drains inland or to the sea;
   - rock stacks;
   - finite strain and vergence in crumple zones;
   - its glaciation record, if any.
6. ~~**Export** and a basin-tile view in the page.~~ Done: the **L2 tiles** view, a tile summary on tap, and downloads of the landmass summary and of one tile.

**C. L0 additions (recording only)**

7. **Climate history per parcel:** time spent in each climate belt, and glaciation exposure. This gives weathering maturity for later soils and erodibility.
8. **Icehouse index** from polar land, ocean gateways and weathering. It replaces item 4's setting.

**D. Performance, when needed**

9. A multigrid preconditioner or a Web Worker for large landmasses.

**E. Parked or decided against**
- Glacial buzzsaw.
- Jitter in L0 (aliasing would become noise, not go away).
- Thrust-sheet and fold generation (balanced cross-sections, critical taper, Biot buckling wavelengths). This belongs to L2; item 1 is shaped to serve it.
- Limits accepted in L1 v2:
  - crumpling stops once L0 welds two plates;
  - crust is seeded from what exists at the window's start (crust that is new during the window comes only from arcs and hotspots);
  - the coupling is one-way.

- **Inputs:** smoothed elevation, formation stacks and domains, fault zones, events, plate velocities (for present-day uplift).
- **Terrain:** uplift field from present plate motion and active faults; stream-power erosion (river incision ∝ drainage area^m × slope^n) with an O(n) implicit solver, as in Braun & Willett's FastScape. Run only long enough to relax wavelengths between L0's and L1's resolution.
- **Climate:** latitude bands and prevailing winds; orographic precipitation with rain shadows (e.g. the Smith & Barstad linear model); temperature from latitude and lapse rate; biomes by Köppen class.
- **Geology:** refine fault zones into strands. Strike sets orientation, slip sets spacing, and the across-to-along ratio sets style; oblique zones split into parallel thrust and strike-slip faults. Re-segment rifted coasts with transform offsets guided by fault records. Glaciation episodes enter as discrete events.
- **Outputs:** drainage graph with outlets pinned at tile edges; precipitation and temperature fields; formation contact surfaces; fault strands; event catalogue of glaciations, major floods and large landslides.
- **Tools:** Landlab (Python, from CSDMS, the Community Surface Dynamics Modeling System); fastscapelib; goSPL for a global version.

### L2: region (~500 m)

Built so far: the selection (section 5.6), the relief L2 starts from, the selection's outline burnt in at 500 m, and the L2 pass (all below). The plan for the boundaries, as first written (the burn-in below tests and revises it):
- **What moves divides between L1 and L2 is mostly resolution, not time.** Redrawing only L1's unresolved relief with another seed changed the basin map drastically (adjusted Rand index 0.35); a 1,600-fold range of erodibility barely changed it (0.89–0.98); divide persistence rose with local relief (section 5.4). Divide migration from unequal erosion across a divide (Willett et al. 2014) is estimated, not measured, at 0.1–1 km per Myr: a cell or two of L2 over an L2 run.
- **So L2's added relief should be a deterministic function of position** (fixed to the planet, from the seed and L1's fields), as L1's unresolved relief already is. It can then be evaluated anywhere without simulating the region around it.
- **Burn in the outer boundary in strips.** A point's basin is a property of its downhill path, so the high-resolution divide is found by routing flow in a strip around each L1 divide on the selection's outer edge: a path that leaves the strip into a basin's interior belongs to that basin. The strip's half-width is the unresolved relief's amplitude over the resolved slope across the divide: narrow on crests, wide on lowlands. Divides inside the selection are L2's own to resolve. Because the strip depends only on position, neighbouring selections compute the same outline and their seams agree. The relief generator must be fixed before any outline is computed.
- **During the L2 run** the burnt-in outline is a fixed water divide, open only at the outlets. The same idea, solving tiles independently and stitching only their perimeter cells, is behind Barnes's parallel depression filling and flow accumulation for DEMs (digital elevation models) of a trillion cells (2016, 2017).
- **Test before trusting it:** route a whole test region at 500 m (200 × 200 km, 160,000 cells) and compare its outline with the strip-only outline.

**The relief L2 starts from** (`src/l2.js`, `L2.relief(sheet, B)`) is a function z(x, y) on the L1 sheet's map (km, fixed to the reference plate), deterministic in the world seed. Any point can be evaluated without the rest, in about 1–3 µs. It is the sum of three parts:
1. **Base:** L1's surface from the drainage pass (10 km, rivers and lakes included), Catmull-Rom bicubic. L1's values are cell means, so the values the interpolant passes through are corrected first (four iterations of e′ += e − M e′, M the interpolant's cell mean, a 5-point separable filter). Averaged over each L1 cell, the base is L1's surface.
2. **Added relief** in the band L1 cannot resolve: gradient-noise octaves of lattice spacing 16, 8, 4, 2 and 1 km, each turned by its own angle (a lattice's axes show faintly in its noise, and octaves sharing them added up to a plaid).
   - *Amplitude.* Each octave's RMS is L1's local roughness r₁ continued down the spectrum, r₁ (λ / 40 km)^H with H = 0.5. r₁ is the RMS, over a 30 km Gaussian window, of L1's land surface less its blur at one cell (relief at about 20–60 km), on land only, so the step from land to L1's flat sea is not roughness. H = 0.5 because L1's own spectrum falls by about 0.7 per octave (section 5.4), as Earth's topography does at these scales (Huang and Turcotte 1989).
   - *Band ownership.* Noise of lattice spacing 16 km still has power at wavelengths L1 resolves. The added relief's own L1 cell means (8 × 8 samples per cell of the octaves of 4 km and up, computed when first needed) are interpolated like the base and subtracted. Averaged over L1's cells the added relief is then about zero; this removes most of the 16 km octave, so L2's relief is in practice about 1–10 km.
   - *Structure.* Where L1 shortened the crust, the relief is stretched along strike by a = 1 + 5 × shortening (at most 4), as ridges, valleys and so divides run along strike in fold-thrust belts (trellis drainage). Direction and amount come from a lattice of nodes 20 km apart, each with its own rotated, stretched copy of the noise, blended with smooth weights; where the nodes agree (no shortening) the copies are one field, otherwise they are independent and the blend is renormalized to keep the RMS.
3. **L1's trunk rivers** (draining at least 1,000 km²) burnt in, so the added relief cannot dam them. The surface is at most the river's level (L1's filled surface, linear along each 10 km reach) on a floor 0.6 km either side of the river, rising beyond it at the steeper of L1's valley side and 3 RMS of the added relief per km: a narrow trench, which doesn't flatten the valley. L1 knows a river's path only to a cell, and its reaches run at 0° and 45°, so the river is drawn where the reach is in coordinates warped by a smooth noise (RMS 1.5 km, 16 km lattice): it wanders within L1's cell instead of running as a canal. A river entering the sea is held at sea level. Small streams are L2's own.

Measured (`tests/test_l2.js`, on the L1 test landmass; and the genesis landmass at 1 Gyr):
- *Pure function of position:* a point evaluated alone equals the same point in a grid, and a second build agrees.
- *L1's scale:* averaged over L1's 10 km cells, the 500 m surface differs from L1's surface by 4.2 m RMS, against 31 m RMS of L1's cell-to-cell relief (genesis: 12 m against 99 m; before the base correction and the band's high-pass, 53 m).
- *Amplitude:* the octaves of 4 km and less have 0.419 RMS per unit of L1 roughness, against 0.418 designed.
- *Structure:* on a made-up plain under 30% shortening (a = 2.5), the strike comes out at the designed angle, and gradients across strike vary 6.21 times as much as along it (a² = 6.25). The test landmasses have too little shortening to test (under 1%, and 6% on genesis).
- *Rivers:* on 34 chains of five reaches, the lowest pass between two points on a river's floor, along any path on a 500 m grid, is no higher than the river at the upstream point (0.0 m). With a V-shaped trench and no floor it was 58 m, because a 500 m grid may have no cell on the river's line.
- *Not yet addressed:* on coastal lowland a few tens of metres above sea level, hollows of the added relief dip below sea level and show as ponds and lagoons, a few kilometres inland at most. L1 treats land below sea level as inland water, and so will L2 unless this is changed.

**The outline burnt in at 500 m** (`L2.boundary(relief, B)`):
- *Membership.* A 500 m point belongs to the L1 basin of the place where its water leaves the land: the sea, an inland sink (an L1 lake with no outlet, a dry lake floor or a river's dry end; or land below sea level), or the edge of the window routed, where it takes the L1 label of its cell. Depressions are filled (priority flood, as at L1). On a window, the flood starts from those outlets, lowest first, and each cell takes the label of the cell it was reached from (Barnes et al. 2014). L1's sea cells take the basin of the nearest L1 land cell, so a coast that L2 moves seaward keeps its basins. A point is in the selection if its label is a selected basin.
- *Stitched tiles.* Tiles are 4 × 4 L1 cells (40 km, 80 × 80 points) on a fixed lattice. Each is evaluated once, ever, and reduced to a small graph (pass A): flooded from every edge point and every outlet, each its own seed, keeping for each pair of seeds whose catchments touch the lowest level at which water crosses between them. Minimax ("lowest highest point") distances between seeds are the same on this graph as on the points (Barnes et al. 2014, 2016). For a selection, a domain of tiles is routed:
  - *Stitch:* the tiles' graphs are joined across their shared edges (touching edge points, crossing at the greater of the two heights) and flooded as one graph, from the outlets and from the domain's outer edge (labelled by L1). This gives the exact level, label and origin arriving at every tile's edge points.
  - *Pass B:* each tile is flooded again from its outlets and its edge points at those levels.
  - *Ties at the seams:* a point on a slope is reached at its own height from every lower neighbour; routing the whole domain gives it the label of the lowest one, where its water goes. The stitched graph knows the level but not the neighbour, so edge points then take the label of their lowest neighbour (across seams too), and tiles whose edge labels changed are flooded again until none change.
  - *How far to route:* a label is a guess if its water reached the domain's outer edge where L1's label is not sure (an L1 cell within a cell of another basin). The domain starts as the tiles of the L1 cells within 2 cells of the selection's outer edge, and grows by the tiles beyond the edge points that guesses came from, until no point in those cells is a guess. While it grows, only the tiles near the outline are flooded again, and the relabelling at seams waits for the last round. Every L1 cell of the domain with no guess gets L2's labels; the rest keep L1's.
- *Measured against routing the whole box at once* (the selection's box at 500 m, its edge 40 km out; widening that margin to 120 or 200 km changed nothing), as mean displacement of the outline (land assigned differently, over the outline's length), on one core:

  | Selection | L1's 10 km outline | Stitched | Tiles, rounds, time |
  |---|---|---|---|
  | Genesis, 1 Gyr: 54 coastal basins, 54,200 km² | 2.51 km | 0 km | 109 tiles, 4 rounds, 2.3 s |
  | Genesis: a mid-size coastal basin, 36,300 km² at 500 m | 4.10 km | 0 km | 75, 3, 1.3 s |
  | Genesis, another landmass: 7 basins, 301,500 km² at L1 (1.2 million points) | 5.12 km | 0.002 km (9 km², 36 points) | 184, 4, 4.9 s |
  | Test landmass: 28 coastal basins | 2.93 km | 0 km | 42, 2, 0.8 s |
  | Test landmass: a lowland basin (18,300 km² at L1, 14,600 at 500 m) | 5.29 km | 0 km | 72, 5, 1.7 s |

  L1's outline is off by 2.5–5 km on average, a quarter to a half of an L1 cell, and a basin's area at 500 m can differ from L1's by 20% in flat country. The 36 points left are near-ties: points balanced between two exits at the same level.
- *Seams.* A stitched domain depends on the selection (it ends wherever the selection's water needed it to), so neighbouring selections are not guaranteed to agree point for point. Tested: two neighbouring basins refined independently disagree on 2 of the 287,000 points in their 717 shared L1 cells (a point one calls the other's basin and the other a third's), and no point falls in both. The test allows 1 in 10,000.
- *The first version* routed each tile on its own with a halo of up to 24 L1 cells (the halo doubling while a label was a guess). It was exact and selection-independent, but in flat country every 40 km tile routed 100–200 km of neighbourhood itself: 20 s for the lowland basin above, 21.5 s for the 301,500 km² selection. Stitching evaluates each tile once.
- *Several cores.* Pass A (a tile's surface and graph) is independent per tile, so the page hands it to Web Workers, one fewer than the cores (at most 8). Each is built from the L2 module's own source and the relief's plain inputs (`reliefInputs`), so it computes exactly what the main thread would (tested bit for bit). Workers start when a basin is tapped (each takes about a second to build the relief) and are kept for the L1 run. The stitching stays on the main thread. Measured in Chromium on the 301,500 km² selection, with 3 workers: tiles 1.3 s (2.6 s on the main thread alone), stitching 2.0 s, 3.4 s in all against 4.4 s; a second preview reuses its tiles (2.1 s). Without workers (or if one fails) the main thread evaluates the tiles, at least one per frame.

**The L2 pass** (`L2.runPass`, **Run L2** in the L2 panel) evolves the selected basins at 500 m for a run length (default 2 Myr, steps of 0.25 Myr) at present rates: L2's counterpart of L1's drainage pass.
- *Domain:* the selected land points (inside the burnt-in outline). Base level is the sea (land below sea level included) at sea level, and L1's inland sinks at their own height, as in the burn-in. The rest of the outline is a water divide held fixed, as planned: divides would migrate a cell or two in a run.
- *Each step:* uplift at L1's present tectonic rate; routing by priority flood with an epsilon gradient, then steepest descent between nodes jittered within their cells (so rivers aren't locked to the grid's 0° and 45° directions); discharge, drainage area weighted by L1's runoff over its land mean; incision by stream power, implicit, n = 1: dh/dt = U − K k_rock A_w^m S (Braun and Willett 2013), K = 0.1 and m = 0.5 as in L1's pass. A point in a filled depression is left alone.
- *Erodibility* (the decision of section 6, item 3, that detailed erodibility belongs to L2): a factor by the rock at the surface of each sheet cell (sediments 3, arc volcanics 1.5, basalts and metamorphics 0.7, cratonic basement 0.5), young sediments weaker still (up to 10× when fresh, e-folding over 30 Myr), log-averaged over the markers by volume, scaled to a geometric mean of 1 on land and log-interpolated to each point. Fault damage is not yet included.
- *No hillslope diffusion:* at 0.001–0.01 m²/yr, D t / dx² is 0.4–4% in 10 Myr at 500 m. Hillslopes are L3's.
- *Depressions,* as at L1: one its catchment's sediment can fill within the run (the run's incision plus background denudation of 20 m/Myr) becomes a flat the water crosses; the others stay lakes.
- *L1's scale kept:* the run's net change, averaged over each L1 cell, is interpolated (as the relief's base) and removed, so L2 changes only the band it owns.
- *Outputs:* the 500 m surface, flow directions, weighted drainage area, Strahler order of the channels (1 km² and more), main-stem length upstream, lakes. Views: terrain with rivers (5 km² and more), change, drainage area, stream order; the surface is drawn bilinear between points.
- *Measured* (2 Myr):

  | Selection | Points | Time | Deepest incision | L1-scale change removed | Flow directions changed | Hollows filled / lakes | Hack's exponent |
  |---|---|---|---|---|---|---|---|
  | Test landmass, 28 coastal basins | 43,105 | 1.6 s | 131 m | 27 m RMS (2.4 m left) | 33% | 35 / 0 | 0.67 |
  | Genesis, 54 coastal basins | 216,751 | 3.0 s | 123 m | 8 m RMS | 23% | 3,065 (4,670 km²) / 0 | 0.65 |
  | Genesis, 7 mountain basins (3.2 km at the tap) | 1,156,859 | 7.8 s (setup 2.8 s, 0.45 s a step), +200 MB | 766 m | 69 m RMS | 12% | 17,823 (28,400 km²) / 46 | 0.57 |

  Every point drains to an outlet, and a run is deterministic (both tested). Hack's exponent comes out at 0.57–0.67, against about 0.5–0.6 for Earth's rivers: a little high on lowland selections, not yet understood. The network is still changing over 2 Myr (12–33% of flow directions), so the run length matters; it is a setting, and recorded with the session.

**Preview L2 relief** in the L2 panel first burns in the selection's outline (with progress; tiles are kept for the L1 run, so later selections reuse them), then draws the selected basins (the rest faded) from this function, per pixel at the view's own scale, leaving out octaves shorter than 2 pixels; scroll or pinch to zoom, drag to pan, tap for the elevation. The outline is drawn as a line, the land outside faded. The panel reports the land area inside the 500 m outline against L1's, the tiles, how often the domain widened, any L1 cells left with guesses, and the time spent evaluating tiles (and on how many workers) and stitching; and L2's planned grid (0.5 km; no time step yet, as nothing runs), the pixel size, the octaves drawn and the times. The session's checks include the relief at the tap and the outline's area, so a rerun verifies both.

Also:
- River network refinement constrained by L1's drainage graph.
- Strata tilted and folded near faults, so dip comes from fault kind and proximity.
- Water table from precipitation and formation permeability.
- Glacial carving where L1 recorded ice.
- Karst potential: carbonate formations (within marine shelf units) above the water table.

### L3: landscape (~50 m)

- Hillslope diffusion; landslide events on steep slopes over weak contacts; river terraces from recorded base-level changes.
- **Caves:** anisotropic shortest-path networks through soluble formations, biased along fractures, bedding planes and past water-table levels (cf. Paris et al. 2021).
- **Springs:** where terrain cuts a permeable-over-impermeable contact, where caves reach the surface, and along fault zones.

### L4: site (~5 m)

- Cliff bands from differential erosion at formation contacts (a resistant caprock over weak rock).
- Talus aprons, stream meanders, floodplains.
- Vegetation instances by point processes conditioned on biome, slope, water and soil.

### L5: detail (~1 m and below), generated lazily

- Deterministic per tile from hash(seed, level, tile), with overlapping halos so tiles agree at edges.
- Small cliffs, boulders, stream beds; caves dressed as signed distance fields along cave-graph edges, since a heightfield cannot represent overhangs.
- Storage note: a 100 km × 100 km region at 1 m is 10¹⁰ cells (about 40 GB of float32 heights), so this level must be computed on demand.

### Cross-cutting conventions

- Every handoff is band-limited; categorical data travels as geometry (contacts, traces) or proportions.
- Every level emits its own event catalogue for the levels below.
- Tile seams: pin drainage outlets and fluxes at tile edges; compute halos.

---

## 7. Known limitations and next steps

- **Resampling between equal-resolution grids** is the root of most remaining artifacts: about 5% of continental cells flicker in displayed rock per step, about 10–13 pinholes per 1,000 continental cells, and ragged plate boundaries. A finer body grid or moving particles would remove this family of problems at roughly double the cost.
- **The 75 km thickness cap is now the main continental volume sink.** With collisions keeping all consumed crust, excess above the cap is deleted. Over 200–600 Myr in one test it removed about 280k cell-km, against a net volume change of about −64k cell-km (−17%). Land still levels off because the area thermostat and accretion compensate, but routing the excess sideways (more collapse) instead of deleting it would be cleaner.
- **Erosion base-levels continents more than Earth's do.** Long-run land is about half the continental crust, versus Earth's three-quarters.
- **Fault traces are fragmented** into segments a few hundred kilometres long, and oceanic faults are not extracted.
- **Thrusting is simplified.** Collisions push a metamorphic unit rather than duplicating a neighbour's stack, and formation order is by age, so thrust sheets (older over younger) come out in the wrong order. Layers carry no dip.
- **Two fault episodes per parcel;** older history is overwritten in long-lived boundary zones.
- **Formation identity is inferred** at handoff from rock type, origin plate and age. Recording an explicit depositional-episode ID when each unit forms would be cleaner.
- **Small islands and peninsulas vanish** in the smoothed handoff (anything under about three cells). L1 should restore island chains from hotspot and arc records.
- **Seafloor-age streaks** come from nearest-cell resampling as plates rotate.
- **The marker core (section 3.14)** removes resampling from the material but is 3.5–4 times slower, settles at less continental crust (the lattice core's continents lean on resampled arc crust), and shows older seafloor. The view still samples moving markers, so the grid flickers. Next, per `L0_MESH_SPEC.md`: speed (warm starts, workers), then local forces and a torque balance (Stage 2).
- **L1:** see section 5.7.

---

## References and related work

- Cortial, Peytavie, Galin, Guérin (2019), *Procedural Tectonic Planets*.
- Cordonnier et al. (2016), large-scale terrain generation from tectonic uplift and fluvial erosion; (2017) coupled ecosystem and erosion; (2023) glacial erosion.
- Génevaux et al. (2013), terrain generation using procedural models based on hydrology.
- Paris et al. (2021), synthesising geologically coherent cave networks.
- England & McKenzie (1982), the thin viscous sheet model of continental deformation; Houseman & England (1986) on the Argand number.
- Ahnert (1970), denudation rate against relief in mid-latitude drainage basins.
- Braun & Willett (2013), the FastScape implicit stream-power solver.
- Barnes, Lehman & Mulla (2014), Priority-Flood depression filling; Zhou, Sun & Fu (2016), its faster variant used in L1.
- Hack (1957), the length–area scaling of river networks used as a check on L1.
- PlaTec / WorldEngine and Tectonics.js (game-oriented plate simulators); pyGPlates (plate reconstructions).
- Wilson cycle: J. Tuzo Wilson's observation that oceans open roughly along older sutures.
