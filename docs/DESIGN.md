# Hierarchical Terrain Generator — Design Document

Oct 24, 2026 · @Nick

## Purpose and scope

The generator builds a planet once from coarse tectonic simulation in about 10 seconds in the browser, then refines any region on demand, level by level, down to streams, rockfall and trees. The target is a wide variety of plausible local maps whose character is explained by material and history, not physically accurate steady states.

**Goals**

- Every local map is consistent with every coarser map above it: zooming out never contradicts what the player saw zoomed in.
- The same world seed and tile address always produce the same terrain, on any machine, in any exploration order.
- Variety comes from rock type, climate and event history inherited from coarser levels, with noise only as unresolved subgrid detail.
- The global stage completes in about 10 s; each refinement tile completes in under about 1 s so zooming feels continuous.

**Non-goals**

- Quantitative geophysical accuracy. Process models are chosen for the morphology they produce, not their calibration.
- Solving mantle flow (Stokes equations). Plate motion is prescribed kinematically.
- Reorganization across tile seams at fine levels, such as a river migrating across a tile edge. The design freezes this deliberately (see Boundary handling).

## Core design principles

Four rules govern every level; most later decisions follow from them.

**1. Consistency invariant (downscaling).** Each child level must reproduce its parent when coarse-grained, and may add only detail the parent did not determine. The restriction operator R is defined explicitly per field and checked after every child run:

```latex
\left| R(\text{child}) - \text{parent} \right| \le \varepsilon_{\text{field}}
```

For elevation, R is the block average over each parent cell. For water, R is the total discharge leaving through each parent-cell edge. This is multigrid restriction and prolongation, or the renormalization group run backwards.

**2. Material and history over noise.** Variety comes from rock type (granite domes, shale badlands, limestone karst), climate, and events such as glaciation, base-level drops and river capture. Noise is injected only as unresolved subgrid detail in each level's initial condition, and the level's processes then erode it into structure. Noise in, structure out.

**3. Relaxation time sets how much history a level runs.** Each process has a response time τ. A child starts from its parent's final state and runs for a few of its own τ, never replaying the full history. A past event survives as a relic, and is carried in tile state, only when the time since it is shorter than τ of the features it shaped.

| Feature | Approximate response time τ |
| --- | --- |
| Continental relief, mountain belts | 10–100 million years |
| Regional drainage network | 1–10 million years |
| Isostatic rebound after unloading | about 10,000 years |
| Hillslopes | 1,000–100,000 years |
| Stream channels, talus | years to centuries |

These τ values are order-of-magnitude design guides, not calibrated constants.

**4. Determinism.** Terrain is a pure function of world seed and tile address. Exploration order, hardware, and cache state must never change the result.

## Scale hierarchy

Four levels refine by about 32× in space each, so every tile is roughly 1024 × 1024 cells. Each level uses the cheapest process model that produces the right morphology at its scale.

| Level | Tile span | Cell size | Time simulated | Process model |
| --- | --- | --- | --- | --- |
| L0 Planet | whole sphere | about 30–50 km | about 500 million years | Rigid plates rotating about Euler poles; boundaries classified as convergent, divergent or transform by relative velocity; crustal thickening; isostasy for elevation; latitude climate bands |
| L1 Region | about 1,000 km | about 1 km | about 10 million years | Stream power erosion with uplift and hillslope diffusion, solved implicitly (FastScape-style); drainage graph extraction; orographic precipitation |
| L2 Landscape | about 32 km | about 32 m | about 100,000 years | Threshold landsliding (sandpile-style cellular automaton at the angle of repose), diffusion, alluvial fan deposition, glacial carving where cold, lithology-dependent erodibility |
| L3 Site | about 1 km | about 1 m | no time stepping | Statistical placement conditioned on inherited fields: channel width and meanders from discharge, talus below cliffs, boulders, trees from a spatial point process |

**L0.** Plates are Voronoi regions on an icosahedral sphere grid. Each crust parcel carries a material ID (oceanic, craton, accreted terrane, arc volcanic, sedimentary basin) and an age, which seed all later lithology. Plate-boundary history is logged as events with timestamps.

**L1.** Stream power erosion: elevation z changes with uplift U, erosion scaling with drainage area A and slope S, and diffusion D.

```latex
\frac{\partial z}{\partial t} = U - K A^{m} S^{n} + D \nabla^{2} z
```

Erodibility K varies with material ID and climate, which is where most regional variety comes from. This is the strongest candidate for GPU acceleration.

**L2.** Local processes act on the L1 surface plus subgrid noise. Glacial carving runs only if the event log contains a glaciation within the relaxation time of valley forms.

**L3.** No dynamics. Each feature type is a stamp or point process whose density and parameters are functions of slope, aspect, soil depth, discharge, precipitation and temperature. L3 may need splitting into a 10 m and a 1 m level if tile generation exceeds budget.

## Tile state and data structures

A tile is an immutable record: its address, raster fields, a river graph, and an event log. Children read their parent's tile and the shared corner and edge records; they never read sibling tiles.

```typescript
interface TileAddress { level: number; face: number; i: number; j: number }  // quadtree on a cube-sphere face

interface Tile {
  addr: TileAddress;
  seed: bigint;                 // hash(worldSeed, addr), see Determinism
  fields: {
    elevation: Float32Array;    // m
    material: Uint8Array;       // rock / crust type ID
    uplift: Float32Array;       // m per year
    precip: Float32Array;       // mm per year
    temperature: Float32Array;  // mean annual, °C
    regolith: Float32Array;     // loose soil and rock depth, m
    drainageArea: Float32Array; // m², includes area upstream of the tile
  };
  rivers: RiverGraph;
  events: TectonicOrClimateEvent[];
}

interface RiverGraph {
  nodes: { x: number; y: number; z: number }[];
  edges: { from: number; to: number; discharge: number; order: number }[];
  crossings: { edgeId: string; t: number; discharge: number; inflow: boolean }[]; // where rivers cross tile borders
}

interface TectonicOrClimateEvent {
  kind: "collision" | "rifting" | "glaciation" | "baseLevelDrop" | "riverCapture";
  time: number;                 // years before present
  region: Polygon;
  params: Record<string, number>;
}
```

Three choices here are deliberate:

- **drainageArea includes upstream area outside the tile.** Stream power erosion needs the full catchment; without it, rivers near an inflow edge erode as trickles.
- **The river graph is topology, not a raster.** It is what children must honor; their rasters are free to vary around it.
- **The event log is pruned per level** to events younger than that level's longest relaxation time, so children only inherit history that can still be visible.

## Boundary handling

Seams are decided by shared corner and edge records generated from the parent alone, never by whichever neighbor tile was generated first. This makes the world independent of exploration order.

**Dimension-ordered generation.** Within a level, entities are generated in order of dimension. Each depends only on the parent tile and its own seed:

1. **Corners (0D),** shared by four tiles: elevation and material from the parent, plus a seeded perturbation.
2. **Edges (1D),** shared by two tiles: an elevation profile interpolated from the parent between its two corners, plus seeded 1D detail. The edge also chooses where each parent river crosses it, anywhere within the stretch of edge that the coarse river cell spans, and records the crossing's discharge and drainage area.
3. **Interiors (2D):** solved with edge profiles as fixed (Dirichlet) boundary values and crossings as fixed river entry and exit points.

**Erosion boundary conditions in interiors.**

- Edges hold parent-derived elevations, so they do not act as drains. Only designated outflow crossings are sinks.
- Inflow crossings inject their upstream drainage area and discharge as point sources.
- The pinned river is pre-carved as a channel between its entry and exit before erosion runs, and flow directions that would let a tributary capture it are forbidden.
- A halo of a few parent cells around the tile is computed but discarded, so interior processes never see the hard boundary directly.

**Blending.** Smooth fields such as regolith depth and temperature may be crossfaded across halos. Channels, ridgelines and fault traces are never blended, because averaging two channels yields two half-channels.

**Accepted limitations.**

- Rivers cannot migrate across, or capture tributaries across, tile edges at levels below the one that placed them.
- A child-scale channel running parallel to an edge has no clean representation; it is either nudged into the interior or allowed one crossing pair.
- Features spanning a seam at L2 (landslides) use a staged pass: stage k of a tile may read stage k−1 of its eight neighbors, the same pattern chunk-based voxel games use for trees and structures.

## Determinism and seeding

Every random draw comes from a counter-based generator keyed by a hash of the world seed and the entity's address, so any tile, edge or corner can be regenerated alone and identically.

- **Seeds:** tile seed = hash(worldSeed, level, face, i, j); edge seed = hash(worldSeed, level, sorted pair of adjacent tile addresses); corner seed likewise from its four tiles. Sorting makes the edge seed identical from either side.
- **Generator:** a counter-based random number generator such as Philox, or a stateless integer hash, indexed by (seed, cell, draw number). Unlike a stateful generator, results do not depend on how threads are scheduled.
- **GPU floating point:** parallel reductions sum in hardware-dependent order, so cumulative quantities can differ across GPUs. Flow accumulation, discharge and anything compared against a threshold (landslide triggers, river initiation) use deterministic ordering or fixed-point integers.
- **Golden tests:** a fixed set of seeds and tile addresses has stored hashes of every output field. Any change to the hash is either a deliberate, versioned model change or a bug.
- **Versioning:** each world records the generator version. A model change creates a new version, never silently alters existing worlds.

## Platform, compute budget, and caching

Generation runs in Web Workers with WebGPU compute shaders where available and a WebAssembly CPU fallback. L0 gets the 10 s budget once; every finer tile must finish in about 1 s because tiles generate on demand while the player zooms.

| Level | Target time per tile | Main cost | Notes |
| --- | --- | --- | --- |
| L0 Planet | about 10 s, once per world | plate advection and isostasy over about 500 steps | Light on GPU; could run on CPU |
| L1 Region | under 1 s | implicit stream power solve, flow routing | Flow routing is a graph traversal, the hardest part to parallelize |
| L2 Landscape | under 1 s | landslide cellular automaton, diffusion | Embarrassingly parallel on GPU except avalanche cascades |
| L3 Site | under 0.3 s | point-process sampling, stamping | Cheap; bounded by mesh and texture upload |

These times are design targets, not measurements; milestone 2 benchmarks L1 first because it is the main risk.

**Memory.** A 1024 × 1024 tile holds 6 float32 fields at 4 MB each plus a 1 MB material field, about 25 MB uncompressed. Keep only the tiles on the current zoom path and their neighbors resident.

**Caching.** Finished tiles are stored in IndexedDB keyed by (generator version, world seed, tile address). Because generation is deterministic, the cache is purely a speed-up and can be cleared at any time. Neighbor tiles are prefetched in the direction the camera is moving.

## Validation and testing

Three properties are tested automatically: consistency across levels, determinism, and variety. Plausibility is judged against geomorphic statistics rather than by eye alone.

| Property | Test | Pass criterion |
| --- | --- | --- |
| Consistency | Apply restriction R to each child and compare with parent | Within ε per field; river discharge conserved at every crossing |
| Seam continuity | Compare shared edge values from both adjacent tiles | Bitwise identical |
| Order independence | Generate a tile neighborhood in two random orders | Bitwise identical outputs |
| Determinism | Golden seeds on at least two GPU vendors and the CPU fallback | Stored hashes match |
| Variety | Generate many local tiles with the same parent class but different seeds and histories | Spread in summary statistics, not just in noise |
| Plausibility | Compare with known landscape statistics | Hack's law (main channel length grows as roughly drainage area to the 0.6 power); slope–area relation with a concavity near 0.5; Horton stream-order ratios in observed ranges |

The plausibility statistics are standard results in geomorphology and give falsifiable checks on whether the drainage looks real. A landscape that fails them reads as artificial even when players cannot say why.

## Milestones

Five phases, each closed by an automated test rather than a date. The region benchmark comes second because L1 performance is the riskiest assumption in the design.

&#91;embedded content: roadmap · 5 phases, each with a pass/fail gate\]

Each gate reuses a test from Validation and testing, so a phase is done when its test passes on the golden seeds.

## Risks and open questions

The largest risk is L1 flow routing on the GPU within 1 s; the largest design question is how much history each level should carry.

**Risks**

- **L1 performance.** Implicit stream power needs flow receivers in topological order, which is sequential by nature. Fallback: parallel priority-flood depression filling plus fewer, larger time steps, or a coarser L1 grid.
- **Seam artifacts at L2.** Pinned edges may produce visible straight lines in landslide or fan deposits. Mitigation: wider discarded halos and the staged-neighbor pass.
- **Cross-GPU determinism.** Driver and vendor differences in floating point may break golden tests even with fixed ordering. Mitigation: fixed-point arithmetic for all thresholded quantities.
- **Generic look at L3.** Point processes conditioned on too few fields will look uniform. Mitigation: condition on history (fire, landslide age) as well as present fields.

**Open questions**

- [ ] Sphere grid at L0: icosahedral, cube-sphere, or a flat world with wraparound?
- [ ] Refinement factor: 32× per level with four levels, or 16× with five?
- [ ] Is L3 one level or two (10 m and 1 m)?
- [ ] Should players be able to see time (scrub history) or only the present?
- [ ] Which events beyond glaciation are worth logging, given their relaxation times?

## References

Cited from memory and not re-verified for this doc; check titles and years before relying on them.

- Braun, J. & Willett, S. D. (2013). A very efficient O(n), implicit and parallel method to solve the stream power equation. *Geomorphology*. The FastScape algorithm used at L1.
- Cordonnier, G. et al. (2016). Large scale terrain generation from tectonic uplift and fluvial erosion. *Computer Graphics Forum*. Stream power adapted for terrain generation.
- Cortial, Y. et al. (2019). Procedural tectonic planets. *Computer Graphics Forum*. Kinematic plate model close to L0.
- Génevaux, J.-D. et al. (2013). Terrain generation using procedural models based on hydrology. *ACM Transactions on Graphics*. River-network-first terrain construction.
- Bak, P., Tang, C. & Wiesenfeld, K. (1987). Self-organized criticality. *Physical Review Letters*. The sandpile model behind L2 landsliding.
- Salmon, J. K. et al. (2011). Parallel random numbers: as easy as 1, 2, 3. *Proceedings of SC11*. Counter-based generators including Philox.
- Hack, J. T. (1957). Studies of longitudinal stream profiles in Virginia and Maryland. *USGS Professional Paper 294-B*. Origin of Hack's law.
