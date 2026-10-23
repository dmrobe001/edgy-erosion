# The model

This page states what the simulation computes, in enough detail to re-implement it.
Every rule below is in `src/`; file names are given with each process.

Status labels:

- **[standard]**: a published, widely used landscape-evolution formulation. Cited.
- **[adapted]**: a standard idea re-expressed for this grid or combined in a way of our own. The adaptation is described.
- **[heuristic]**: chosen to behave sensibly and look right. Not derived. Treat its parameters as knobs, not measurements.

Units: metres, years, pascals, kilograms.

## State and grid

The domain is a rectangle of `cols × rows` pointy-top hexagons in "odd-r" offset layout
(`src/hexgrid.js`). Neighbouring centres are exactly `dx` apart, so the centres form a
triangular lattice and every cell has six equidistant neighbours.

| quantity | value |
|---|---|
| cell area | (√3/2) dx² |
| shared edge length | dx/√3 |
| row spacing | (√3/2) dx |
| domain extent | ≈ dx·(cols + ½) by (√3/2)·dx·rows |

Per-cell state (`src/simulation.js`):

```
rock[i]       bedrock surface elevation (m)
soil[i]       loose material on top: regolith, alluvium, talus (m)
z[i]          rock[i] + soil[i]
upliftCum[i]  total uplift so far (m), so strata can be looked up in the rock's own frame
```

A snapshot of these arrays, plus the parameter object, the time and the random number
generator state, reproduces a run exactly.

## Boundary conditions

Each grid edge (north, south, east, west) is either **open** or **closed**. Open edge
cells are fixed base level: they do not uplift, and any material that reaches them leaves
the domain. Closed edges have no flux across them. If `seaLevel` is set, every cell below
the current sea level is also an outlet, and rivers deposit sediment in the sea up to sea
level. Sea level can vary in time (see Coasts). If neither an open edge nor any sea is
left, for example because uplift has raised the whole sea floor, the lowest cell becomes
the sink of an internally drained basin.

## One time step

Operator splitting, in this order (`Simulation.step`):

1. **Uplift.** `rock += U(x, y)·dt`. Patterns: uniform, dome (Gaussian), tilt, or fault
   (a step: one side rises, the other doesn't). [standard]
2. **Weathering** (soil production). Bare rock turns into soil at
   `P = P0(rock) · exp(−soil/h*) · (1 + frostBoost·frost) · precip^moistureExponent · (1 + productionBoost·veg)`.
   The exponential "soil production function" is from Heimsath et al. 1997 (Nature 388).
   We integrate it exactly over the step: `h(t+dt) = h*·ln(exp(h/h*) + P·dt/h*)`.
   [standard]; the climate and vegetation multipliers are [heuristic].
3. **Vegetation.** `veg = Vmax · (1 − exp(−soil/hv))`, where `Vmax = P/(P + P½)` when
   `cover: 'auto'`, or a fixed number. Bare rock has no vegetation. [heuristic]
4. **Flow routing** (`src/flow.js`).
   - Priority-flood depression filling with an ε-gradient over lakes (Barnes, Lehman & Mulla 2014, Computers & Geosciences 62). [standard]
   - Each cell's receiver is chosen among its lower neighbours on the filled surface. With `randomReceivers` it is drawn at random, with probability ∝ drop^p, and re-drawn every step. This is the hex analogue of "Rho8" (Fairfield & Leymarie 1991, Water Resources Research 27). It removes the strong bias towards the six lattice directions; see the comparison in "Known limitations". [adapted]
   - Drainage area A is accumulated upstream-first, weighted by precipitation/(1 m/yr). With `routing: 'mfd'` area is shared among all lower neighbours ∝ drop^p (Freeman 1991, Computers & Geosciences 17). [standard]
   - The order in which cells leave the priority queue is a valid "stack" (every cell comes after its receiver), so no separate topological sort is needed.
5. **Rivers** (`src/processes/fluvial.js`).
   - Incision follows the stream power law `E = K·A^m·S^n`, with `K = K_soil·exp(−α·veg)·(1−c) + K_rock·w(x,y)·c` and `c = exp(−soil/H*)`. Here `w` is the rock-strength variation field. The soil term removes only soil and the rock term removes only rock, as in SPACE (Shobe, Tucker & Barnhart 2017, Geoscientific Model Development 10). [standard]
   - The time integration is implicit, downstream-first (Braun & Willett 2013, Geomorphology 180–181), so it is stable for any dt. For n ≠ 1, Newton's method solves each cell's implicit equation. [standard]
   - Simplification: we solve for the total lowering with the blended K and then split it between soil and rock in proportion. If the soil runs out within the step, the soil share is truncated. [adapted]
   - Sediment is routed upstream-first as a volume flux Qs. Each cell deposits at rate `D = G·Qs/A` (the ξ–q model: Davy & Lague 2009, JGR 114; dimensionless form from Yuan et al. 2019, JGR Earth Surface 124). We solve this per cell in closed form, explicitly rather than with Yuan's implicit iteration. `G = 0` is pure detachment-limited incision. [adapted]
   - Lakes trap sediment up to their spill level. Lake beds do not incise.
6. **Creep** (`src/processes/hillslope.js`). Soil flux per unit width is
   `q = D·S/(1 − (S/Sc)²)·(1 − exp(−soil/h0))`, the nonlinear law of Roering, Kirchner & Dietrich 1999 (Water Resources Research 35), with a depth dependence. Only soil moves, so a cell cannot export more soil than it has. The divergent factor is capped at `maxNonlinearity`. The run uses explicit sub-steps, with `dt_sub ≤ dx²/(6·D·fmax)`. On the hex lattice the flux across one edge changes the two cells' heights by `D·f·Δz·(2/3)/dx²·dt`, which is the standard isotropic hex Laplacian. [standard]
7. **Rockfall** (`src/processes/rockfall.js`). Each cell's drop to its lowest neighbour is treated as a rock face, with height Δz and angle β = atan(Δz/dx). The face is unstable when Δz exceeds the Culmann critical height
   `Hc = 4C/(ρg) · sin β · cos φ / (1 − cos(β − φ))`
   (Culmann 1866; used in landscape evolution by Densmore, Ellis & Anderson 1998, JGR 103). [standard] C and φ are thickness-weighted means over the layers the face cuts, so soft rock at the foot of a cliff weakens it (undercutting). [adapted] An unstable face fails at random with probability `1 − exp(−dt·(1 + frostFactor·frost)/τ)`. The wedge above the failure plane, at angle (β+φ)/2, is removed, capped at half the drop. It is dumped as loose debris on the lower cell. [heuristic]
7b. **Coasts** (`src/processes/coast.js`), when a sea exists. See the next section.
8. **Avalanching.** Soil or talus steeper than the critical slope Sc (raised by vegetation) slides to lower neighbours until it is at Sc, conserving volume. This is a "thermal erosion" relaxation (Musgrave, Kolb & Mace 1989, SIGGRAPH). It turns rockfall debris into talus slopes. [standard in terrain generation; heuristic as geomorphology]
9. **Boundaries.** Soil on open-edge cells is removed and counted as exported.

## Coasts

Sea level is `boundary.seaLevel` plus either a sinusoid (`seaLevelAmplitude`,
`seaLevelPeriod`, `seaLevelPhase`) or a piecewise-linear curve `seaLevelCurve`. Every
cell below the current sea level is sea. Waves act on a coastal band: sea shallower than
3 × `closureDepth`, plus land cells touching the sea.

1. **Wave exposure (0 to 1).** For each direction in the wave climate, a ray runs from
   the cell across open water until it reaches land or leaves the grid. A ray that leaves
   the grid is treated as open ocean.
   - Wave growth with open water upwind (fetch): `1 − exp(−fetch/fetchScale)`.
   - Damping by water shallower than the breaking depth (`waveHeight`/0.78):
     `exp(−shallow length/dissipationLength)`. This is what limits platform width
     (Trenhaile 2000, Marine Geology 166; Anderson, Densmore & Ellis 1999, Basin Research 11).
   - Rays run along the six lattice directions. Rays along those directions land on cell
     centres, so each cell's fetch follows from its upwind neighbour's: one O(n) sweep per
     direction instead of marching rays. The price is 60° angular resolution. A wave
     direction in between is shared between the two nearest lattice directions.
   - A share `refraction` of the energy is spread equally over all six directions. This
     stands in for waves bending round headlands. Without it, waves from a single
     direction leave long fingers of land whose sides are never attacked. [adapted / heuristic]
2. **Cliff retreat.** Every land cell touching the sea is attacked from the side. This
   includes a low shore on a gentle coast, which becomes a cliff as it retreats into
   rising ground. A notch at its foot deepens at
   `R = waveErodibility(rock at sea level) · weakness · exposure · exp(−soil/beachShielding)` (m/yr).
   The beach-shielding exponential follows SCAPE (Walkden & Hall 2005, Coastal
   Engineering 52). When the notch reaches dx, the column collapses to low-tide level.
   A share `debrisRetained` of the collapsed rock stays as loose sediment and the rest
   leaves as fines. A cliff retreats at most one cell per step, so keep dt below
   dx / (fastest retreat rate). [adapted: notch-and-collapse is our grid version of lateral retreat]
3. **Platform lowering.** Submerged cells lower at `platformRatio · R · f(depth)`. Here
   f = 1 within the tidal range and decays as `exp(−(depth − tide/2)/(depthScale·waveHeight))`
   below it: wave erosion decays with depth (Sunamura 1992, *Geomorphology of Rocky
   Coasts*). [standard form; the coefficients are heuristic]
4. **Attrition.** Shore sediment is ground to fines at
   `attrition · exposure · f(depth)` (m/yr). Without it, cliff-fall debris would bury the
   cliff foot for ever. With it, beaches persist mainly where waves are weak or sediment
   supply is high. [heuristic]
5. **Shore sediment transport.** Soil in the band diffuses with
   `D = marineDiffusivity · exposure · exp(−depth/closureDepth)`, down past the closure
   depth and up to the run-up limit. This builds beaches, shorefaces and delta fronts, as
   depth-decaying marine diffusion does in basin models (e.g. Yuan et al. 2019). Explicit
   substeps are capped at `maxSubsteps`. If the cap binds, D is reduced to the stable
   value. [standard in basin models]

Wave-driven sand transport along the shore is **not** modelled, so there are no spits,
barrier islands or cuspate forelands. Those need a shoreline model with sub-cell
position and a much faster time step, such as the Coastline Evolution Model of Ashton,
Murray & Arnoult 2001 (Nature 414). Refraction is only the crude share described above.
There are no tides beyond widening the erosion band: no tidal flats, inlets or marshes.

## Lithology

`src/lithology.js`. Rock type depends on position and on elevation *in the rock's frame*
(`z − upliftCum`), so strata travel upward with uplift. The stratigraphic column is a
list of layers below a top surface:

```
columnTop(x, y) = top − tan(dip)·s(x, y) + undulation(x, y)
```

where `s` is the distance along the dip direction from the domain centre. Layer
thicknesses are true thicknesses, perpendicular to the bedding; a layer of thickness t
dipping at angle δ spans t/cos δ vertically and t/sin δ across the map. The column can
repeat indefinitely, both up and down, or end in a basement rock. Vertical bodies (`disc` plugs, `dike`
walls) override the strata below their `top`. A smooth random field, `strengthVariation`,
multiplies erodibility and divides cohesion. It stands in for patchy fracture density.

Each rock has `erodibility` (stream power K), `cohesion` and `frictionAngle` (Culmann),
`weathering` (P0), `waveErodibility` (sea-cliff retreat rate) and a display colour. The defaults in `src/params.js` are
order-of-magnitude values for landscape-scale rock masses. Rock-mass cohesion is far below
intact-sample strength.

## Checks that are run (`npm test`)

| claim | test |
|---|---|
| neighbours are exactly dx apart and symmetric | `hexgrid.test.js` |
| routing leaves no pits; every cell drains to an outlet; area is conserved; the stack is topological | `flow.test.js` |
| at steady state, implicit stream power satisfies K·A^m·S^n = U to 1e-6 (n = 1 and n = 2) | `processes.test.js` |
| hex-lattice creep is isotropic: a bump's variance grows at 2Dt in both x and y, to within 1% | `processes.test.js` |
| avalanching conserves volume and leaves no soil slope above Sc | `processes.test.js` |
| Culmann matches its closed form; weak cliffs fail and strong ones stand | `processes.test.js` |
| volume budget closes: Δvolume = uplifted − exported, to 1e-8 relative, for full runs (including coastal presets) | `simulation.test.js` |
| with wave damping off, a cliff retreats at exactly its wave erodibility, to within one cell | `coast.test.js` |
| a widening platform slows retreat; a beach shields a cliff; a lagoon is less exposed than the open coast | `coast.test.js` |
| marine terraces: with only uplift and waves acting, shoreline angles sit within 6 m of (highstand − tide/2 + U·age), for the three youngest terraces. A 10% error in uplift rate or sea-level period fails the test | `coast.test.js` |
| same parameters give bit-identical output; snapshots resume exactly | `simulation.test.js` |

These are all checks of internal consistency and numerics. **None of them shows that the
landscapes are geologically realistic.** That is a matter of judgement and of calibrating
against real topography, which has not been done.

## Known limitations and scale dependence

- **Grid scale.** Rockfall only "switches on" when cliffs are many times `dx·tan φ` high
  compared with what the rock can hold up. On coarse grids, cliffs stay smeared into
  slopes. The critical soil slope, the creep stability limit and the rockfall
  probability all depend on `dx`. Stream power is grid-independent in principle, but
  small channels are resolved only down to one cell. Changing `dx` therefore changes the
  landscape, not just its resolution. Treat `dx` as a physical choice.
- **Lattice anisotropy.** With single-direction routing, channels lock onto the six lattice
  directions and domes erode into hexagonal stars. Random receivers remove most of this,
  but channels on uniform planar slopes still run parallel. Pinned base-level edges with
  fast uplift also produce boundary escarpments.
- **Time step.** Operator splitting is first-order. Rockfall is stochastic per step, so
  when dt is much larger than the failure timescale, every unstable face fails every step.
  Results converge as dt shrinks, but they are not dt-independent.
- **Not modelled:** explicit water depth or flood hydraulics; groundwater; chemical
  dissolution (karst); glaciers; wind; isostasy; lateral channel migration and meanders;
  deep-seated landslides; sediment grain size; porosity (rock and soil volumes are
  interchangeable 1:1).

## A short glossary

- **escarpment, scarp**: a long, steep slope or cliff separating two gentler surfaces. From French *escarpe*, the steep inner bank of a fortification ditch, and Italian *scarpa*, the sloping foot of a wall.
- **mesa / butte**: flat-topped hills capped by resistant rock. Mesas are wider than they are tall, buttes narrower. Spanish *mesa*, "table"; French *butte*, "mound".
- **cuesta**: an asymmetric ridge on gently dipping strata, with a steep scarp face and a gentle dip slope. Spanish, "slope".
- **talus / scree**: the apron of rockfall debris at the foot of a cliff, at its angle of repose.
- **regolith**: all loose material above solid bedrock. "Soil" in this code means regolith.
- **base level**: the lowest elevation to which rivers can erode, here the open edges or the sea.
- **knickpoint**: a sharp step in a river profile that migrates upstream as the river adjusts.
