# What the code's structures are

The project kept running into the question "is a plate an object or a process?" This file answers it for every structure in the code, using one fixed set of kinds. Each kind gets its identity from a different place, so it needs a different data structure and a different rule for crossing between levels.

## The kinds

| Kind | What it is | Where its identity comes from | Stored as |
|---|---|---|---|
| **Space** | The grid or frame things live on | Its geometry | Meshes, projections |
| **Field** | A value at every point | None: values at positions (Eulerian) | Grids |
| **Material** | Matter that moves and keeps its history | The same matter over time (Lagrangian) | Parcels, markers, rock columns |
| **Record** | History written into material or a field by processes | Whatever carries it | Rock units, fault episodes, time bins, logs |
| **Process** | A rate law: state in, change out | Its law and parameters | Functions, constants |
| **Source** | A forcing from outside the model's own dynamics | Given, not derived | Lists, tables |
| **Event** | A process seen through a window coarser than its duration: a threshold crossing, a stochastic trigger, or a topology change | Time, place, kind | Log entries |
| **Pattern** | A coherent shape the state passes through: picked out by a detection rule, followed by a linking rule | The two rules; names are conventions on branches of a lineage graph | Ideally recomputed and tracked, not stored |

Three working rules follow from the table:

1. **Material is the source of truth; fields are views of it.** Where both exist, a field is gathered from material, edited, and written back.
2. **Patterns should not cross levels.** What passes from a coarse level to a finer one should be fields, material, records and sources. Each level detects its own patterns at its own resolution.
3. **The erasure time of a feature belongs to other processes.** A feature persists until the processes acting at its scale remove its imprint. Its relic is whatever record it left in the material.

A pattern can be made a **primitive**: stored with an identity, because the model uses it as a cause. That is a legitimate modelling choice, but it should be labelled as one. Section "Findings" lists where the code does this.

## L0 (`src/sim.js`)

| Structure | Kind | Notes |
|---|---|---|
| `Cube` (`pos`, `nb`, `area`, `lookup`) | Space | The fixed (mantle) frame. Plates use the same cube as their body frame. |
| Plate body grids `P.T`, `P.B`, `P.S`, `P.K`, `P.col` | Material | Crust parcels in the plate's frame: thickness, birth time, sediment, continental flag, rock column. A parcel is a body cell. **Source of truth for crust.** |
| `P.pid` | Material identity | Parcel tags carried through splits, welds and hole fills while L1 is tracking. |
| `wT`, `wB`, `wS`, `wK`, `wCol`, `wR`, `wRA` | Field | Fixed-frame view of the plate grids. `gather()` builds it, the surface processes edit it, `scatter()` writes it back. |
| `owner`, `bodyIdx` | Field | Which plate, and which parcel, sits on each fixed cell: the join between field and material. |
| `elev` | Field | Derived from thickness, age and sediment by isostasy (`elevOf`). |
| `ColumnPool` (rock, born, thickness, plate per unit) | Material and Record | Stratigraphy: what built each column, when, on which plate. The long-term memory of L0. |
| `FaultPool` (two episodes per column) | Record | Deformation history carried by the parcel, in the plate frame. A relic that outlives the boundary that wrote it. |
| `slab`, `slabT` (10 Myr bins) | Record (field with a time axis) | Oceanic area consumed at trenches, in the mantle frame. L1 turns it into dynamic topography. |
| `trench`, `arcAdd`/`arcPlate`, `collAdd`/`collPlate`, `dT` | Process (per-step fluxes) | Where this step's subduction, arc growth, collision thickening and erosion act. |
| `C` constants (`KE`, `KD`, `ARC`, `COLL`, `SPLIT_TAU`, ...) | Process parameters | |
| `step()` stages: rotation, coverage and `resolve()`, gap filling (`newOcean`), surface processes, `syncUnsampled` | Process | `resolve()` is the subduction rule: which parcel stays on top where two plates overlap. `newOcean` is the only process that creates material. |
| Plate `P` (`id`, `Q`, `W`, `alive`, `parent`, `mergedInto`) | **Pattern made primitive** | In nature a plate is a region of lithosphere moving rigidly together, which emerges from mantle convection. Here it is stored because it causes motion. Its lineage is explicit (`parent` on split, `mergedInto` on weld) but is not kept as a graph. |
| Angular velocity random walk (`plateEvents`) | Process (stochastic forcing) | Plate speeds relax toward each plate's own preferred speed. |
| `split(P)` | Event (topology change) | Triggered by a rate (`SPLIT_TAU`) and a random cutting plane, not by stress. The side with `s ≤ 0` keeps the old id. That is the naming rule, and it is arbitrary. |
| `merge(a, b)` | Event (topology change) | Triggered when a collision counter passes `MERGE`. The larger plate keeps its id. |
| Collision zone (`pairs`, `pairLoc`, `collLogged`) | Pattern | Keyed by the plate pair. A decaying counter, located by its last contact point only. Its identity dies when either plate merges or splits (`purge`). The orogen itself survives as material: `OROGEN` rock and fault records. |
| Continent | Pattern (not stored) | Connected parcels with `K = 1`. Nothing tracks it. That is consistent with the table above. |
| Craton | Material | A rock type (`CRATON`) on columns, not a region. |
| `hotspots` (`v`, `i`, `wasUp`, `last`) | Source | Fixed in the mantle frame and immortal: 3–6 placed at the start. No plume head, drift, waning or death. |
| `events` log | Event | Mixes two kinds: causal events that change the state (rift, suture, plate lost) and observations of threshold crossings (island emerges, collision begins). Capped at 400. |
| `mergeLog` | Record | Which absorbed parcel filled each cell of the surviving plate: provenance for L1. |
| Formations and domains (handoff, README §3.11) | Pattern | Detected once at handoff by union-find on rock columns. No tracking needed, since they are only used at one instant. |

## L1 (`src/l1.js`)

| Structure | Kind | Notes |
|---|---|---|
| `landmassAt` (land, buffer, enclosed seas, connected crust) | Pattern detection | The landmass is detected **at the present** only. |
| `Replay` frames (`id`, `Q`, `parent` per frame) | Record | Parcel identities through the window. The landmass's identity going back in time is **its material**: "its crust is every parcel that ends up there". A pattern detected now, carried backwards by matter. This is the cleanest identity rule in the project. |
| Reference plate followed through its split ancestry | Pattern lineage | The map frame is a plate's lineage branch, followed back through `parent`. |
| `Projection` | Space | Azimuthal equidistant map about the landmass. |
| Sheet markers (position, thickness, area, deformation gradient F, parcel, rock column) | Material | Source of truth for the sheet's crust. Strain is exact per marker. |
| Sheet grid: thickness, reference thickness, inherited elevation, strength, viscosity, velocity `u`/`v`, finite strain | Field | Gathered from markers each substep. The velocity is solved for (process), then markers move. |
| Sheet momentum balance (conjugate-gradient solve) | Process | |
| Flexural root and elevation | Field (derived) | |
| `dynTopo` | Process applied to a Record | Turns L0's slab bins into a deflection field. |
| `addVolcanism` (from L0's `tecRec`) | Source, imported | L0's arcs and hotspots arrive as a thickening flux on the markers. |
| `erode` (Ahnert-type denudation) | Process | A smooth stand-in for rivers on the sheet grid. |
| `history`, `budget` | Record (diagnostics) | Aggregates per step; not used by any process. |
| `climate` (precipitation, evaporation, runoff) | Source (field) | Latitude tables plus orographic adjustment, at the present only. |
| `basins()`: elevation `e`, area `A`, discharge `Q`, receivers `rec`, uplift `U` | Field | The drainage pass (process) acting for a few Myr. |
| `basin` labels, lakes, closed basins, inland seas, `exitOf` | Pattern | Detected once at the end. A basin's identity is its outlet. No lineage: there is only one instant. |
| `handoff` tiles | Fields and Material for L2 | Fields and rock, which is what the working rules ask for. **Also** basin labels, which L2 then burns in (see Findings). |
| `basinSelection` | Pattern (chosen by the user) | A union of basin patterns. |

## L2 (`src/l2.js`)

| Structure | Kind | Notes |
|---|---|---|
| Point grid at 500 m, tiles, box | Space | |
| `reliefFrom`: base, added octaves, trench | Field with no process | The base is an interpolation of L1's elevation (a view). The added octaves are noise at every scale from 1 to 16 km with **no process that made them**. The trench stamps L1's trunk rivers (a pattern) into the field. |
| Burn-in labels (`boundary`, `refine`, `tileGraph`) | Pattern carried across levels | L1's basin membership enforced at 500 m: an L1 pattern imposed as a constraint on L2's field. |
| `erodibility` (`lnK` from the sheet's rock types and ages) | Material → Process parameter | Rock deciding a rate. The intended direction for material to influence processes. |
| `runPass`: uplift `U`, runoff weight, stream power, routing | Process | Forced by L1's fields (sources from L2's point of view). |
| `z`, `A`, `rec`, `dz` | Field | |
| Lakes, filled depressions | Pattern | Detected at the end of the run. |
| Strahler order, Hack's exponent | Pattern statistics | Diagnostics of the network. |
| L1-scale drift removal | Consistency constraint | Not a kind: the rule that keeps L2's large scales equal to L1's. |

## UI and sessions (`src/ui.js`, `tools/session.js`)

| Structure | Kind | Notes |
|---|---|---|
| Selections (landmass cell, tapped basin, radius) | Pattern chosen by the user | Logged as actions. |
| Session export, action log, snapshots | Record | Of the user's actions, for reproducing a run, not of the world. |

## Findings

Ordered from most to least consequential for the work ahead.

1. **L2's relief is a field with no generating process.** The added octaves are features without a cause, stamped at every scale at once. This is the ontological form of the grit diagnosed earlier. The fix in that language: put variation into material (erodibility) and sources (uplift), and let a process make the relief.
2. **Patterns cross from L1 to L2.** L1's basin labels are burnt into L2, and its trunk rivers are cut into L2's relief. Both impose an L1 pattern on L2's field rather than letting L2 detect its own. The nested rerun with a margin, discussed earlier, would pass fields, material and sources only, and leave basin outlines at 500 m to L2.
3. **Plates are a pattern made primitive, with incomplete lineage.** That is justified in a kinematic model. But:
   - splits happen at a fixed rate on a random plane rather than from stress;
   - which side keeps the id is arbitrary;
   - lineage is split between `parent` and `mergedInto` instead of one graph.

   Writing the lineage as a graph (nodes are plate ids over time intervals; edges are continue, split and weld) would make questions like "what did this plate descend from?" answerable, and would make the naming rule explicit.
4. **Hotspots are immortal, fixed sources.** Nature's hotspots have a life cycle: a flood-basalt birth, drift of 1–2 cm/yr, pulsing, and death after 50–150 Myr. As a source this is cheap to add, since nothing in the model has to derive it.
5. **The event log mixes causes and observations.** Rifts, sutures and lost plates change the state. "Island emerges" and "collision begins" only report threshold crossings. Separating them would make the log usable as a record (for example, ages of sutures) rather than only as narration.
6. **Collision zones lose identity when plates change.** The counter is keyed by plate pair and purged on a split or weld. The material relic (orogen rock, fault records) survives correctly, so this matters only if a collision zone needs to be followed as a pattern, for example to drive L1's sheet directly.
7. **Material records are rich, but L2 reads little of them.** Columns and fault episodes carry deformation history. L2's erodibility uses only rock type and the age of young sediments. Fault damage is already on the to-do list.
8. **No process is represented as events with an intensity.** Volcanism, erosion and collapse are all continuous fluxes. Event-like processes (eruptions, landslides, glacial advances) would need an intensity field, a magnitude distribution and an event log. This only matters once a level resolves their imprints.

## Kinds and the scale map

For the length–time scale map (the "Landscape Stommel diagram" artifact), the kinds sit as follows:

- A **pattern** or **material** feature is a box. Its top edge, the erasure time, is set by the **process** curves that cross it.
- An **event** is a process whose duration lies below the window's time step.
- A **source** is whatever the window takes as given from outside.
- A **record** is how a feature above the window (too slow for the run) reaches the window: as inherited material.
