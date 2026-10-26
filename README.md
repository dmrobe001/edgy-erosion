# Procedural Earth-like worlds from geological history: L0 plate tectonics, L1 landmasses

This repository holds the first two levels of a multi-level procedural generator for Earth-like game worlds, plus sketches of the levels below them. L0 is a planet-scale plate tectonics simulation. L1 takes one landmass at the present and replays the last 10–30 Myr of its history at finer resolution, with its crust deforming as a thin viscous sheet driven by L0's plates, so collisions crumple it into mountain belts and plateaus; it then settles the landmass's drainage basins on a 10 km grid, the tiles for L2. Both run in a browser (`plates.html`, a single self-contained page) and in Node.js for headless testing.

The goal is diversity that comes from **material and history** rather than from noise layered on hand-made rules. Each level simulates or procedurally generates the structure that is plausible at its own scale. It hands down only what it can actually resolve, plus a record of what happened, so finer levels inherit causes as well as shapes.

Contents: the big picture; quick start; L0 in detail; the L1 handoff format; L1 landmass; sketches of the levels below; known limitations and next steps. `DEVLOG.md` explains how the current version came about, including what didn't work.

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

---

## 2. Quick start

**In a browser:** open `plates.html`. It needs no server and no build step; it loads only Google Fonts. Press **Run**.

- **Views:** Elevation, Rock type, Fault history, L1 preview, Plates, Seafloor age, Crust thickness.
- **Projections:** a draggable globe (with optional camera spin) or an equirectangular map.
- **Tap any spot:** a panel under the map shows that cell's rock column, fault episodes, and plate motion.
- **Extract handoff:** pauses the run, traces the fault network, builds formations and domains, and overlays fault traces.
- **Download JSON:** saves the handoff. This uses the Claude artifact runtime's downloads capability (`window.claude.use("downloads")`); opened as a plain local file, the button stays hidden.
- **L1 landmass:** run L0 for a while (1–2 Gyr gives an interesting history; it keeps a saved state every 10 Myr), press **Pick a landmass**, tap land, then **Run L1**. L0 replays the window twice (tracing the crust, then driving the sheet) and ends exactly where it started. The L1 map shows terrain, crustal thickening, strain rate, surface rock and drainage basins, with rivers, zoom and pan, over L0's own elevation. A mid-sized continent takes well under a minute; a supercontinent a few minutes.

**Rebuild the page** from its parts after editing `src/`:

```sh
./build.sh        # concatenates src/head.html + src/sim.js + src/l1.js + src/ui.js into plates.html
```

**Headless tests** (Node.js 18+, no dependencies):

```sh
node tests/run_all.js              # quick pass over everything (about two minutes)
node tests/test_l1.js --n=96 --h=50 --win=20 --t0=1000 --land=3000   # L1 checks at the browser's defaults
node tests/test_budget.js --t=2400 # full-length crust budget run
node tests/bench.js --t=300        # ms per step at each resolution
```

### Repository layout

```
plates.html            the built single-file page (what the artifact serves)
build.sh               rebuilds plates.html from src/
src/head.html          markup, styles, explanatory notes; ends with an opening <script> tag
src/sim.js             the simulation core: a self-contained SIM module, usable in browser and Node
src/l1.js              L1: landmass pick, replay, thin viscous sheet with markers, flexure, crust budget, drainage pass (L1 module)
src/ui.js              rendering, controls, inspector, handoff panel, L1 panel
tests/                 headless tests and benchmarks (see section 3.13)
history/               simulation cores from earlier milestones, referenced by DEVLOG.md
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
| `bench.js` | ms per step, extraction and export cost at each resolution |

Representative numbers from development, in Node.js on a single core:

- **Cost per step:** about 14 ms at n = 64 and 22–25 ms at n = 96 after warm-up. The first tens of steps are slower; short benchmarks overstate cost. Extraction takes 0.2–0.5 s; the export is 2–7 MB.
- **Budget:** with a 50% target, land levels off at roughly 22–34% after the first few hundred Myr, continental crust holds at 46–50%, and peaks stay at 7–8 km.
- **Plate count and seafloor:** plates stay between 9 and 13; mean seafloor age wanders between about 70 and 220 Myr (Earth: about 60), so seafloor recycles somewhat slowly.
- **Events:** roughly 20–30 rifts and 20–30 sutures per 800 Myr, plus hotspot islands.
- **Rock flicker:** about 5% of continental cells change displayed rock per step, about 90% of it from resampling.

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
- **Pick.** At the present, a tap on land selects its landmass: the land cells connected to the tapped one, plus a ring of submerged continental crust around them. The plate carrying most of the land is the reference plate.
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
- **Viscosity.** Power law with exponent 3: η = η₀ · strength · (ε̇/ε̇₀)^(1/3 − 1), floored at 2% of ε̇₀. Old cratons (basement older than 500 Myr) are ten times stronger.
- **Parameters.** Two dimensionless numbers set the behaviour:
  - the **Argand number** Ar (default 3) sets η₀ at the reference strain rate ε̇₀ = 50 km/Myr ÷ ℓ. Higher means weaker crust that spreads into broader plateaus.
  - the **coupling length** ℓ = √(4η₀/κ) (default 500 km) sets how far a plate's drag reaches into the sheet.
- **Numerics.**
  - Staggered (Marker-and-Cell, MAC) grid, 50 km cells by default.
  - The operator is the gradient of a discrete dissipation function, so it is symmetric positive definite and conjugate gradients apply (Jacobi preconditioner, relative tolerance 10⁻⁴, at most 200 iterations), with two Picard iterations on the viscosity per solve.
  - Edges of the crust are traction-free; cells hold crust where the markers cover at least a quarter of them.
- **Markers.** There are 6 × 6 markers per L0 parcel, one at a random point (a fixed hash) in each of 6 × 6 sub-cells, so marker rows don't alias with the grid or bunch into lines under shortening. Each carries its parcel's identity, plate, rock column (as fractions of thickness), strength, the thickness and elevation it started with (L0's elevation smoothed over three passes, sampled bilinearly), and its own thickness and area.
  - Markers move with the sheet (midpoint rule, substeps limited to half a cell per substep).
  - They thicken by the local convergence: S → S·e^(−div·dt), A → A·e^(div·dt). Each marker's volume S·A is therefore exactly conserved.
  - Grid fields are bilinear weighted averages of marker values.
- **Volcanism.** L0 records the crust its arcs and hotspots add. Its collisions and collapse are not passed on, because the sheet makes those itself. Each marker thickens by what its L0 cell gained, as a new top unit of arc rock.
- **Elevation.** As in the first version: the inherited elevation plus the thickness change, minus a flexural root, e = e₀ + ΔS − 0.84·F[ΔS]. F is a Gaussian of width 1.18 × the 40 km flexural length, and the plate yields beyond 3 km from local balance. The inherited state counts as equilibrium.
- **Erosion during the sheet** is a smooth stand-in for rivers, which a 50 km grid cannot draw. Denudation is 0.5 km/Myr × (e / 5 km)^1.4: 0.5 km/Myr in ranges at 5 km and about 20 m/Myr on lowlands at 0.5 km. These are the magnitudes of Ahnert's (1970) relation between denudation and relief. Eroded rock is stripped from the top of the rock columns and leaves the system. A first version ran stream power on the sheet grid; its 50 km gullies were aligned with the grid, and the drainage pass inherited them as parallel stripe-shaped basins.
- **Budget.** A volume budget tracks initial crust, arcs and hotspots, erosion and crust carried off the grid. It closes to about 10⁻¹⁴ of the initial volume; this is tested.

### 5.3 The drainage pass

With tectonics frozen at the end, the sheet's surface is refined to 10 km cells: bilinear, plus fixed fractal relief of ±(0.15 + 0.15·e) km in the band the sheet cannot resolve (wavelengths from 150 km down), seeded from the world seed. Then 5 Myr of stream-power erosion runs (E = K A^m S, m = 0.5, n = 1, K = 0.1, implicit solver of Braun and Willett 2013, Priority-Flood+ε depression filling of Barnes et al. 2014 in the Zhou et al. 2016 variant, D8 receivers). As in L1's first version, each cell's node (used for distances, hence slopes) is displaced within the cell by up to 0.4 cell widths by a fixed hash, so neighbour directions vary from cell to cell and rivers are not locked to the grid's 0° and 45° directions. That is long enough for the network to organise, but not to change the landscape. Every land cell is then labelled by the outlet it drains to: a point on the coast, or the lowest cell of an inland sink. The map colours basins of at least 10,000 km²; smaller coastal basins are grey.

### 5.4 Measured behaviour

Measured on seed "genesis" in Node.js, unless stated:

- **Sheet solver against exact solutions.** A Newtonian sheet with no buoyancy has its left half dragged at +U and its right half at −U. The exact profile is U(1 − e^(−|x|/λ)). The solver matches it to 0.06% in a narrow strip with free sides (λ = √(3/4)·ℓ) and 0.03% in the middle of a wide sheet (λ = ℓ).
- **Stream-power solver.** It reaches the analytic steady state S·A^m = U/K at every cell to 5 × 10⁻¹⁵.
- **Hack's law** in the drainage pass comes out at 0.53–0.55 on the landmass below and 0.54 in the test setting, with no tuning (Earth's rivers: roughly 0.5–0.6).
- **Lattice alignment in the drainage pass** (landmass below, channels draining 5,000 km² or more): node jitter cuts the share of channel steps in straight runs of four or more identical steps from 52% to 33%. The share of 6-step reaches within 5° of a grid direction falls from 25% to 18% (uniform would be 22%).
- **A landmass at 1 Gyr** (medium grid; 41.5 M km² of land on plate 44, three lobes, two plates at the window's start): pass 1 replays 27 Myr in 2.1 s and traces 4,224 parcels. The sheet is 223 × 185 cells of 50 km with 152,928 markers, at 1.65 s per L0 step (22 steps). The replay and the sheet both leave L0 bit-identical. After 20 Myr: the most thickening is 13 km, in parallel bands at two collision fronts; the highest point is 6.9 km; the strain rate falls from about 0.3/Myr during the collision to 0.02/Myr once L0 has welded the plates. Erosion removed 71 M km³ (1.7 km averaged over the land), arcs and hotspots added 1.5 M km³, and the budget closes to 9 × 10⁻¹⁶. The drainage pass (1115 × 925 cells) finds 5,203 basins, of which about 380 are 10,000 km² or larger; the largest is 1.4 M km². It takes 2.2 s.
- **Solver tolerance.** On the same landmass, capping conjugate gradients at 200 iterations leaves a residual of about 5 × 10⁻³ and costs 1.66 s per step. A fully converged solve (10⁻⁶, four Picard iterations) costs 14.8 s per step. The final surfaces differ by 8 m on average (at most 0.35 km), against a mean thickness change of 1.5 km. A 30 km grid (5.4 s per step, never converging within the cap) gives the same peaks to within 0.15 km and the same most-thickening (12.5 km). The 50 km default is a cost choice, and the answer is not sensitive to it at that level.
- **The thickening bands are contacts, not grid.** At the strongest collision front, two bands about 300 km apart appear at sheet grids of 70, 50 and 30 km, and with or without marker jitter. What follows the grid is their width, about 1.5 cells each time: the shortening localises more narrowly than the grid resolves, as power-law crust tends to. Tracing the markers back shows the bands are two separate contacts. One is the suture between the two plates' crust; the other is a promontory of one margin that hit first. Their spacing comes from the shape of the coastlines L0 brought together.
- **A thickness spike** in an early run (one marker at 431 km for one step) came from passing on L0's orogenic collapse as volcanism: L0 turns an ocean cell into 28 km of continent in one step, and that volume landed on the few markers in the cell. Only arcs and hotspots are passed on now, as thickness per marker.
- **Cost in Chromium:** a 110 M km² supercontinent (365 × 366 sheet cells, 381k markers) took 5.9 s per step, about 2.5 minutes in all, with the drainage pass on 1825 × 1830 cells taking 8.8 s. Small landmasses take seconds.

### 5.5 What it does not do yet

- **Coupling to L0 is one-way.** L0 does not see the sheet's deformation.
- **Sutures stop the crumpling.** When L0 welds two plates, their markers get the same velocity, so a collision crumples the crust only until L0 welds it. That is L0's rule, and L1 inherits it.
- **Crust created during the window is missing.** Arcs that become continent and shelves that emerge have no markers.
- **Ocean is outside the sheet.** Continental crust is the only material.
- **Erosion during the sheet** is a smooth function of elevation, and eroded rock leaves the system instead of filling basins.
- **The drainage pass** has uniform rain and erodibility.
- **Cost.** The conjugate-gradient solve is the cost: a multigrid preconditioner, or a Web Worker, would make large landmasses quicker.

## 6. Sketches of the levels below

These are design sketches, not implementations (L1 above is partly built; its sketch lists what remains). Each level should follow the same pattern as L0: deterministic given the seed and its parent's output; owns one band of wavelengths; hands down band-limited fields plus history.

### L1: landmass: what remains

Built so far: section 5. The current plan, organised around the L1→L2 handoff (each item either feeds it, decides whether something belongs in it, or is parked):

**A. L1 work, in order**
1. **Finite strain per marker.** A 2×2 deformation gradient on each marker, updated with the sheet's velocity gradient: total horizontal shortening, its direction and vertical-axis rotation. This is what L2 needs to build thrust sheets and folds in crumple zones.
2. **Closed basins.** A water balance on the depressions the drainage pass finds, with a crude aridity field (latitude belts plus a rain shadow from the sheet's topography and a prevailing wind): a hollow overflows only if inflow beats evaporation over its lake. This decides basins that drain inland versus to the sea.
3. **Erodibility test (a decision gate).** Compare basin maps (adjusted Rand index) with uniform K against K driven by variables that actually predict rock strength:
   - rock type;
   - consolidation with age (young sediment soft, old lithified: ×10–100);
   - fault damage from L0's fault records (crushed zones erode fast, which also gives fault-guided valleys);
   - metamorphic grade.

   Add a threshold hillslope angle as a stand-in for landsliding, so strong contrasts don't build walls. If grounded contrasts of ×100 or more still don't move the major basins, detailed erodibility belongs to L2. Contrasts from type labels alone stay modest (×3–10): within-class variation exceeds the difference between classes, and contacts are blocky at L0's resolution.
4. **Glacial overprint (optional).** Gated by a world setting: whether there was recent ice, and how long ago it left. Only what is regional and permanent: the glaciated-region mask with its age, river diversions, large troughs and lake basins. It shares the lake logic of item 2 and gets the same basin-map test. No glacial buzzsaw: it is indistinguishable from older mountains.

**B. The L1→L2 handoff**

5. **Schema.** Passes are ordered by the age of the event they represent, oldest first. Every imprint carries its event class and age, plus a table of lifetimes by erasure mechanism (hillslope blurring, σ = √(2Dt); knickpoints sweeping up the network; lakes filling). Each basin tile carries:
   - its outlet and drainage graph;
   - whether it drains inland or to the sea;
   - rock stacks;
   - finite strain and vergence in crumple zones;
   - its glaciation record, if any.
6. **Export** and a basin-tile view in the page.

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
  - crust created during the window has no markers;
  - the coupling is one-way.

- **Inputs:** smoothed elevation, formation stacks and domains, fault zones, events, plate velocities (for present-day uplift).
- **Terrain:** uplift field from present plate motion and active faults; stream-power erosion (river incision ∝ drainage area^m × slope^n) with an O(n) implicit solver, as in Braun & Willett's FastScape. Run only long enough to relax wavelengths between L0's and L1's resolution.
- **Climate:** latitude bands and prevailing winds; orographic precipitation with rain shadows (e.g. the Smith & Barstad linear model); temperature from latitude and lapse rate; biomes by Köppen class.
- **Geology:** refine fault zones into strands. Strike sets orientation, slip sets spacing, and the across-to-along ratio sets style; oblique zones split into parallel thrust and strike-slip faults. Re-segment rifted coasts with transform offsets guided by fault records. Glaciation episodes enter as discrete events.
- **Outputs:** drainage graph with outlets pinned at tile edges; precipitation and temperature fields; formation contact surfaces; fault strands; event catalogue of glaciations, major floods and large landslides.
- **Tools:** Landlab (Python, from CSDMS, the Community Surface Dynamics Modeling System); fastscapelib; goSPL for a global version.

### L2: region (~500 m)

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
- **L1:** see section 5.5.

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
