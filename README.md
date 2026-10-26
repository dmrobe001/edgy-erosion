# Procedural Earth-like worlds from geological history: L0 plate tectonics, L1 continents

This repository holds the first two levels of a multi-level procedural generator for Earth-like game worlds, plus sketches of the levels below them. L0 is a planet-scale plate tectonics simulation. L1 (first version) takes one continent from L0, refines it about 8–16 times, and runs it forward for tens of millions of years with river erosion and flexural isostasy, while L0 keeps running underneath and supplies the tectonics. Both run in a browser (`plates.html`, a single self-contained page) and in Node.js for headless testing.

The goal is diversity that comes from **material and history** rather than from noise layered on hand-made rules. Each level simulates or procedurally generates the structure that is plausible at its own scale. It hands down only what it can actually resolve, plus a record of what happened, so finer levels inherit causes as well as shapes.

Contents: the big picture; quick start; L0 in detail; the L1 handoff format; L1 continent (first version); sketches of the levels below; known limitations and next steps. `DEVLOG.md` explains how the current version came about, including what didn't work.

---

## 1. The big picture

### Levels as conditional generators

Each level takes the coarser level's output as a condition and produces finer structure consistent with it. A level does not need to reach steady state. It only needs to replace the interpolation it inherits with structure that is plausible at its own resolution. Discrete, procedurally defined events are welcome wherever they produce richness that would be expensive to simulate physically.

| Level | Cell size | What it owns | Representative methods |
|---|---|---|---|
| L0 planet | ~80–160 km | Plates, crust type and age, rock columns, fault history, events | Rigid plates on a sphere (this repository) |
| L1 continent | ~5–10 km | Drainage basins, climate, biomes, fault strands, glaciation | Stream-power erosion, orographic rainfall |
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
- **L1 continent:** run L0 for a while (1–2 Gyr gives an interesting history), press **Pick a continent**, tap one, then **Run L1**. L0 advances by the window while L1 evolves the continent; the L1 map shows terrain with rivers, erosion, tectonic input and drainage basins, with zoom and pan. A 100 Myr window at ~13 km cells takes about 40 s in a desktop browser.

**Rebuild the page** from its parts after editing `src/`:

```sh
./build.sh        # concatenates src/head.html + src/sim.js + src/l1.js + src/ui.js into plates.html
```

**Headless tests** (Node.js 18+, no dependencies):

```sh
node tests/run_all.js              # quick pass over everything (about two minutes)
node tests/test_l1.js --n=96 --r=8 --win=100   # L1 checks at the browser's default resolution
node tests/test_budget.js --t=2400 # full-length crust budget run
node tests/bench.js --t=300        # ms per step at each resolution
```

### Repository layout

```
plates.html            the built single-file page (what the artifact serves)
build.sh               rebuilds plates.html from src/
src/head.html          markup, styles, explanatory notes; ends with an opening <script> tag
src/sim.js             the simulation core: a self-contained SIM module, usable in browser and Node
src/l1.js              L1: continent selection, refined grid, coupling to L0, erosion, flexure, crust budget (L1 module)
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
| `test_ui_smoke.js` | runs the whole page script against a stubbed DOM: every view, both projections, inspector, extraction, export |
| `test_budget.js` | land, continental crust and peak elevation over time for three seeds; fails if land falls under 10% or crust drifts more than 8 points from target |
| `test_netrotation.js` | net rotation with removal on and off; removal must zero it and leave relative motion unchanged |
| `test_handoff.js` | fault, formation and domain statistics; schema consistency of the exported JSON |
| `test_flicker.js` | rock-type flicker per step, split into surface processes versus resampling |
| `test_pinholes.js` | single-cell ocean pinholes inside freshly absorbed plates versus elsewhere |
| `test_volume_budget.js` | continental volume sources and sinks, via a patched copy of the sim |
| `test_l1.js` | L1: stream-power steady state against the analytic solution; on a real continent, crust budget closure, drainage to the sea, determinism, peak bound; prints Hack's law and the L0–L1 divergence |
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

## 5. L1: continent (first version)

L1 follows the plan agreed when it was started: rather than forcing a refined continent to match L0's final elevation, it starts from L0's state at some time t0 and runs forward through a window (default 100 Myr) **with L0 running alongside**, taking only L0's tectonics. That keeps L1's history its own (knickpoints, relict surfaces) while staying tied to plate motion. The time scales motivate it: isostatic rebound is effectively instant (~10 kyr), a plate crosses an L1 cell in ~100 kyr, rivers respond to a change in uplift in ~1–10 Myr, and the uplift pattern itself reorganises over ~10–50 Myr.

### 5.1 Domain, grid and lineage

- **Domain:** the connected continental crust of the picked cell's plate, in that plate's body frame, plus every cell of the same plate within three L0 cells (the ocean margin, which is base level).
- **Grid:** the same equal-angle cube-sphere, each L0 cell split r × r (r = 4, 8 or 16; ~13 km at r = 8 on the medium L0 grid). Equal-angle subdivision nests exactly, so each fine cell has exactly one parent and positions, areas and neighbours are analytic. Fine cell index = slot × r² + jj × r + ii, where slot numbers the domain's L0 cells. Neighbours are 8-connected; across cube faces they come from the fine-grid lookup.
- **Lineage:** each L0 parcel (slot) is followed as L0 runs. Through a weld it is mapped into the surviving plate's frame. A split is caught in the step it happens: the new plate holds exactly the parcels that moved, at the same body indices. A parcel whose body cell is empty afterwards was consumed; its fine cells leave the domain. Where neighbouring parcels end up on different plates (judged on a majority-filtered plate map, so single flips from L0's resampling noise don't count), the cells along the break become sea.

### 5.2 Each step

One L1 step per L0 step (Δt ≈ 0.9 Myr at the medium grid):

1. **L0 steps** with a recorder on: the thickness each fixed cell gains from arcs, collisions, hotspots and orogenic collapse (not diffusion or erosion, which L1 replaces). These are mapped onto the domain's parcels, smoothed once (volume-conserving) against resampling noise, and interpolated bilinearly onto fine cells as a thickening rate. Thickening that lands on (static) ocean cells is booked, not applied.
2. **Lower-crustal flow:** thick crust spreads sideways, as a diffusion of thickness with coefficient 300 km²/Myr × (T / 33 km)³ (channel flux grows with the cube of channel thickness). It runs on a coarser nested grid (~26 km), so it relaxes collision piles without filling river valleys.
3. **Thickness cap** at L0's 75 km.
4. **Elevation** from thickness by flexural isostasy applied to changes since t0 (the inherited state is taken as in equilibrium): e = e₀ + ΔT − 0.84·F[ΔT], with F a Gaussian of width 1.18 × flexural length (40 km), which matches the half-power point of an elastic plate's response 1/(1 + (kα)⁴). Long wavelengths follow L0's Airy relation; narrow loads are held up by the plate, but by no more than 3 km away from local balance (the plate yields).
5. **Routing:** sea cells (ocean and drowned continent) are outlets; a land edge of the domain is a closed boundary (drainage divide). Depressions are filled with Priority-Flood+ε (Barnes et al. 2014) in the Zhou et al. (2016) variant, receivers are steepest descent (D8) on the filled surface, and land the flood can't reach (cut off by consumed parcels) becomes an inland basin draining to its lowest cell.
6. **Erosion:** detachment-limited stream power, E = K A^m S with m = 0.5, n = 1, solved implicitly down the receiver tree (Braun & Willett 2013). Units are km, Myr and km², so K = 0.1 corresponds to about 1 × 10⁻⁷ yr⁻¹ in SI units. Eroded rock is exported to the sea.

Initial elevation is L0's elevation smoothed over the domain (three passes), interpolated bilinearly, plus fractal relief in the band L1 owns and L0 can't resolve (wavelengths from 1.5 L0 cells down to 2 L1 cells; 50 m in lowlands, rising with elevation), fixed in the body frame and seeded from the world seed. Without it, D8 routing on a bilinear surface draws parallel streams.

### 5.3 Bookkeeping and measured behaviour

The crust budget (continental crust in the domain) tracks the initial volume, tectonic input applied, rock eroded, volume lost to the cap, parcels consumed, cells drowned along rifts, and the small residue of interpolating crustal flow back from the coarse grid. It closes to about 10⁻¹⁴ of the initial volume.

Measured (seed "genesis", Node.js unless stated):

- **Solver:** on a planar grid with uniform uplift, the steady state satisfies S·A^m = U/K at every cell to 5 × 10⁻¹⁵, so the slope–area exponent is exactly −m.
- **Networks:** Hack's law (length ∝ area^h) comes out at h ≈ 0.53–0.54 with no tuning (Earth's rivers: roughly 0.5–0.6).
- **Choosing K.** Over 100 Myr on one continent (n = 64, r = 4), mean land elevation and how far L1's coarse-grained crust drifts from L0's own:

  | K | L1 mean elevation, start → 100 Myr | L0 over the same parcels | mean / rms thickness difference (L1 − L0) |
  |---|---|---|---|
  | 0.03 | 2.10 → 2.47 km | 2.11 → 2.20 km | +1.3 / 3.9 km |
  | 0.1 | 2.10 → 2.20 km | 2.11 → 2.20 km | +0.2 / 3.5 km |
  | 0.3 | 2.10 → 1.40 km | 2.11 → 2.20 km | −2.9 / 5.4 km |
  | 1 | 2.10 → 0.44 km | 2.11 → 2.20 km | −6.4 / 10.0 km |

  K = 0.1 is the default: L1's coarse-grained crust then tracks L0's to within about 0.2 km on average, and mean continental erosion comes out at 25–35 m/Myr, the same order as Earth's average continental denudation. Larger K carves more relief faster but drains the continent toward sea level well ahead of L0. This is a calibration choice, not a derivation.
- **Cost:** about 250 ms per coupled step at ~465k cells (r = 8 on the medium grid) in Chromium, so a 100 Myr window takes about 40 s. In Node at the same size: routing ~80 ms (the flood itself ~40 ms), flexure ~35 ms, crustal flow ~70 ms, L0's own step 70–80 ms.

## 6. Sketches of the levels below

These are design sketches, not implementations (L1 above is partly built; its sketch lists what remains). Each level should follow the same pattern as L0: deterministic given the seed and its parent's output; owns one band of wavelengths; hands down band-limited fields plus history.

### L1: continent (~5–10 km cells): what remains

Built so far: section 5. Still to do, roughly in order: erodibility from the exposed formation (and fault-zone weakening), with the exposed-formation map; deposition of eroded rock (basins, fans, shelves) as new L1 units; climate (winds, orographic rain, temperature) including drift in latitude over the window; fault strands; restoring island chains; growing the domain when crust accretes. The original sketch:

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
- **L1:** erodibility and rain are uniform; eroded rock leaves the system instead of filling basins; ocean cells are static, so tectonic input that lands on them is booked but not applied; the domain never grows (accreted crust is ignored); horizontal shortening is treated as vertical thickening. Coastal cells incise straight to sea level, which leaves a one-cell scarp along most coasts. Lower-crustal flow, the strength limit and the flexure kernel are simple stand-ins. A 100 Myr window is slower than the 30 s target (about 40 s at ~13 km cells); a Web Worker or the WebGPU solver from the `Moving-on` tag would help.

---

## References and related work

- Cortial, Peytavie, Galin, Guérin (2019), *Procedural Tectonic Planets*.
- Cordonnier et al. (2016), large-scale terrain generation from tectonic uplift and fluvial erosion; (2017) coupled ecosystem and erosion; (2023) glacial erosion.
- Génevaux et al. (2013), terrain generation using procedural models based on hydrology.
- Paris et al. (2021), synthesising geologically coherent cave networks.
- Braun & Willett (2013), the FastScape implicit stream-power solver.
- Barnes, Lehman & Mulla (2014), Priority-Flood depression filling; Zhou, Sun & Fu (2016), its faster variant used in L1.
- Hack (1957), the length–area scaling of river networks used as a check on L1.
- PlaTec / WorldEngine and Tectonics.js (game-oriented plate simulators); pyGPlates (plate reconstructions).
- Wilson cycle: J. Tuzo Wilson's observation that oceans open roughly along older sutures.
