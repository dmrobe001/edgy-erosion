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

**The processes do not assume hexagons.** They are written against a mesh contract
(`src/mesh.js`):
- per cell: area and centre position;
- per neighbour slot: neighbour index, centre-to-centre distance, and the length of the
  shared face. The number of neighbours may vary, and slots may be empty.

Diffusive fluxes use the finite-volume rule. A flux q (per unit face length) across the
face between i and j moves volume `q·face·dt`, which changes the two heights by
`∓q·face·dt/area`. On the regular hex grid this is exactly the isotropic hex Laplacian
used before. Two parts still need the regular hex lattice, and they say so if given
another mesh: the wave-exposure sweep, and zooming in.

`GraphMesh.irregular()` builds a deliberately messy test mesh: jittered centres, 2–6
neighbours, gaps between slots, and perturbed faces and areas. It checks conservation
and routing on arbitrary meshes. It is not a true Voronoi tessellation, so it does not
check convergence to the continuum equations. The grid can sit anywhere in world
coordinates (`grid.origin`).

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
   - Incision follows the stream power law `E = K·A^m·S^n`, where S is the slope to the receiver over the actual centre distance. Rock has `K_rock·w(x,y)`, where `w` is the rock-strength variation field. Loose sediment has `K_soil·exp(−α·veg)`. [standard]
   - **Strip, then cut.** Within a step, the river first entrains any sediment on its bed at the soil rate, and cuts rock only for the rest of the step. The stripping time is taken in proportion to the depth stripped. [adapted]
   - Why not the SPACE cover factor `exp(−soil/H*)` (Shobe, Tucker & Barnhart 2017, Geoscientific Model Development 10), which an earlier version used? Each step leaves a sediment layer whose thickness scales with dt, from weathering and deposition. A cover factor read from that layer made incision depend on dt without converging: quartering dt changed 150 kyr of the alpine preset by about 40 m. With strip-then-cut, the time spent stripping is a fraction of the step that does not depend on dt, and results converge (see Checks).
   - Each phase is implicit, downstream-first (Braun & Willett 2013, Geomorphology 180–181), so it is stable for any dt. For n ≠ 1, Newton's method solves each cell's implicit equation. [standard]
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

## Earthquakes

Faults listed in `earthquakes.faults` move in sudden slips instead of smooth uplift
(`src/processes/earthquakes.js`). Each fault has:
- a straight surface trace in fractions of the reference frame, which may run off the
  domain (the domain is then a window onto a longer fault);
- a dip and a long-term slip rate.

The side to the left of the trace's direction is raised, as for the `fault` uplift pattern.

**Which events, and how often.**
- **Magnitudes** follow a truncated Gutenberg–Richter law between `minMagnitude` (5.5:
  smaller events rarely break the surface) and `maxMagnitude`. Each magnitude unit up is
  10^b times rarer (Gutenberg & Richter 1944, *Bull. Seismol. Soc. Am.* 34).
- **Rupture length, average slip and down-dip width** come from the regressions of
  Wells & Coppersmith (1994, *Bull. Seismol. Soc. Am.* 84), all slip types. A rupture
  longer than the fault is cut to the fault.
- **Slip along a rupture** has an elliptical profile with the regression's mean.
- **The event rate** comes from a potency balance (moment balance without the shear
  modulus): rate = slipRate × faultLength / E[slip × length], with the expectation taken
  over the magnitude distribution. Events are a Poisson process with this rate, so on
  average the fault slips at its slip rate while any one spot sees a staircase.

**What an event does.**
- **Displacement.** Vertical throw = slip × sin(dip). The raised side rises by
  `upliftShare` of it and the other side drops by the rest. Both decay away from the
  trace as exp(−distance / W), with W the rupture width. [adapted: a stand-in for the
  elastic dislocation solution (Okada 1985); horizontal motion is ignored]
- **Strata.** Rock carries its layers with it, so the fault offsets the strata in the
  geology view and the cross-section.
- **Shaking.** Ground within about 10^(0.5M − 2.5) km of the rupture shakes for that step.
  Cliffs' Culmann critical height drops by `coseismicWeakening` × shaking, and unstable
  faces fail at once instead of stochastically: earthquake-triggered rockfall.
  [adapted: a one-parameter attenuation, not a ground-motion model]

`sim.earthquake(k, M)` slips fault k on demand (the browser's **Earthquake!** button).
`sim.quakes` holds the catalogue. The `earthquake-scarp` preset is a normal fault with
M 5.5–7.3 events building a range front over an alluvial plain. Its throw is measured
against the basin, which sediment keeps level.

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

## Zooming in (nested patches)

`zoomIn(parent, { center, factor, cols, rows })` in `src/nest.js` builds a finer
simulation of a patch of a running one. Its spacing is the parent's dx/factor. This is
backmapping, as in coarse-grained polymer simulation: the parent is frozen, the patch is
rebuilt at fine resolution with detail added, and the patch then runs its own physics.

The child inherits from the frozen parent:
- **Topography, soil and cumulative uplift**, by linear (barycentric) interpolation on
  the parent's triangle of centres. Elevation is then made mean-preserving: the child's
  average over every parent cell equals the parent's value, by repeatedly adding the
  interpolated misfit (16 passes, residual of millimetres). On top of that goes
  zero-mean roughness: noise at the parent's spacing, scaled to the local parent relief.
  [adapted]
- **Parameters, seed and reference frame.** Strata planes, plugs, fault lines, weak
  zones and rainfall therefore line up exactly with the parent's. Fractional positions
  in the parameters refer to the `frame` rectangle, which is the parent's domain.
- **Time and sea level.**
- **Edges.** Every edge cell is an outlet whose elevation follows the parent's last rate
  of change (`boundaryRate`).
- **Rivers entering the patch.** Parent cells outside the patch that drain into it pass
  their drainage area (`inflowArea`, shared as in the parent's multiple-flow routing) and
  sediment flux (`inflowSediment`) to the child cell where they enter.

Coupling is one-way: nothing the child does changes the parent. `compareWithParent`
reports three things over the parent cells the patch covers:
- the mean elevation, against the parent extrapolated at its last rate;
- the mean rate of change;
- the drainage area of the largest river leaving the patch.

**What the checks found.** Alpine preset, patch on the main river, 150 kyr from the same
starting state (patch-mean elevation 113.3 m):

| run | patch-mean elevation |
|---|---|
| parent, extrapolated at its last rate | 110.6 m |
| parent run on, dt 2000 yr (its own step) | 111.1 m |
| parent run on, dt 500 yr | 98.1 m |
| child, factor 1 (same resolution), dt 500 yr | 99.8 m |
| child, factor 4 (dx 7.5 m), dt 500 yr | 125.0 m |

How to read it:
- **The frozen-rate assumption holds** over this span (110.6 against 111.1 m).
- **The handoff is sound.** The same-resolution child matches the parent run at the same
  time step to within 1.7 m.
- **A time-step effect remains** (111.1 against 98.1 m). It is much smaller than before the
  river fix, and it shrinks with dt, but it is not negligible.
- **Refining changes the answer.** The factor-4 patch stands about 27 m higher. Finer
  cells resolve smaller drainage areas, which need steeper slopes to keep pace with
  uplift. This is the known grid-resolution dependence of stream-power models, not an
  error in the handoff.
- **The child reorganises its own drainage.** Its largest exiting river can lose
  discharge as divides migrate: in `examples/zoom.js` it went from 1.07 to 0.84 km² over
  150 kyr, against 1.13 km² in the parent.

Building the factor-1 control is what exposed the time-step defect described under Rivers.

## The planet

`new Planet(params)` (`src/planet.js`) models a whole planet on an icosahedral geodesic
mesh (`src/sphere.js`).

**The mesh.**
- **Construction.** An icosahedron's triangles are each cut into f × f smaller ones and
  pushed out onto the sphere. That gives 10f² + 2 cells: 12 pentagons, the rest hexagons.
  At f = 40 there are 16,002 cells about 190 km apart.
- **Geometry.**
  - Cells are the spherical Voronoi polygons of the vertices.
  - `dist` is the great-circle distance between centres.
  - `face` is the arc between the circumcentres on either side of an edge.
  - `area` is the spherical polygon area, which sums to 4πR² to round-off.
- **Map coordinates.** `x, y` are equirectangular, used only for drawing. Geometry uses
  the unit vectors in `pos`.
- **Accuracy.** With true geometry, the finite-volume processes are consistent on the
  sphere. The discrete Laplacian converges to the spherical-harmonic eigenvalues
  −l(l+1)/R² at second order in RMS (0.03% error at f = 16 for l = 1); the maximum
  error is first order, near the pentagons.

**Plates** are rigid and kinematic: each rotates about its Euler pole at a prescribed
rate. Plates start as a noisy, randomly weighted growth from scattered seeds, and their
speeds are drawn between `minSpeed` and `maxSpeed`.

The central data-structure choice: **each plate keeps its crust in its own frame**, on its
own copy of the mesh. Each plate-frame cell holds crust thickness, continental fraction,
formation time, soil and last-seen elevation. The world map is only a view of them: after
each small rotation (substeps move no point more than half a cell), every world cell looks
back along each plate's rotation to the plate cell now above it, by nearest neighbour.
Because rigid motion is never interpolated, coastlines stay sharp however far plates
travel.

Each world cell then falls into one of three cases:
- **Covered by one plate**: that plate owns it.
- **A gap** (no plate covers it): sea-floor spreading. New oceanic crust, 7 km thick, is
  added to the plate that was there.
  - The view opens gaps a whole cell at a time, so the new cell is dated from when it
    passes the gap: signed distance along the motion ÷ speed.
  - Without that correction, ages near ridges came out several Myr too old.
- **An overlap**: convergence.
  - **Which plate overrides.** Each plate pair keeps the polarity it took at first contact
    (the more buoyant column then). This is a deliberate choice: deciding cell by cell made
    plate edges interfinger, because ages of similar ocean floors alternate. A column meant
    to sink that floats more than `polarityFlip` (2 km) higher flips the polarity, as when
    a continent reaches a trench.
  - **Continents** beneath the overriding plate are too buoyant to sink and are scraped
    off whole (collision).
  - **Ocean floor** subducts, except for `sedimentAccretion` (half) of its sediment, which
    is scraped off as an accretionary wedge.
  - **Arcs.** `arcFraction` of the subducted volume returns as new continental crust on
    the overriding plate. All additions are spread over about `arcWidth` (250 km) behind
    the margin.
  - **Cleanup.** Plate-frame cells next to a consumed one that no world cell sampled are
    consumed too, so sinking plates leave no stray patches.

Plate behaviour beyond the boundary rules:
- **Collisions lock.** While continents collide along more than `lockingLength`, the two
  plates' relative rotation decays with timescale `collisionLocking`. Larger plates change
  less. [adapted: a kinematic stand-in for the forces that stall collisions] Only columns
  with more than 50% continental material count. Arc slivers counting is a bug the tests
  caught: island arcs alone were locking spreading plates.
- **Plate reorganisation.** Optionally, every `reorganizeEvery` years all plates draw new
  poles. Old sutures can then rift, a crude Wilson cycle.
- **Crustal flow.** Crust thicker than `flowThreshold` (45 km) spreads sideways within its
  plate, with diffusivity rising as (1 + (excess / 10 km)²) [adapted]. This keeps
  collision plateaus Tibet-like rather than spiky.

**Isostasy.** Every column floats (Airy, no flexure).
- **Freeboard.** A column of thickness h and density ρ has dry freeboard
  e = h(1 − ρ/ρm) − K. Below sea level the water's weight deepens it to e·ρm/(ρm − ρw).
- **Calibration.** K is set so that 35 km of continental crust stands at
  `referenceElevation` (600 m).
- **Density** mixes oceanic (2900) and continental (2800 kg/m³) by the continental fraction.
- **Ocean-floor cooling.** Oceanic lithosphere sinks with age as
  ridgeDepth + 350 m·√(age/Myr), until it reaches the isostatic depth of cold crust
  (about 5.6 km, after about 70 Myr) (Parsons & Sclater 1977, *J. Geophys. Res.* 82).
- **Trenches.** The edge of a plate sinking under its neighbour is a further
  `trenchDepth` deeper.
- **Rebound.** Erosion is compensated: removing 1 m of continental rock lowers the
  surface by only 1 − ρc/ρm ≈ 0.15 m.

**Sea level** is solved every step so the ocean's water volume stays constant (bisection
on the hypsometry). Young buoyant sea floor, sediment fill and growing continents push it
up. The water's own load is not re-balanced when sea level moves.

**Surface processes** run as an ordinary `Simulation` on the sphere mesh:
- **What runs.** Rivers, creep and weathering, with zonal rainfall (wet tropics, dry
  subtropics, wet mid-latitudes). Rockfall and waves are off.
- **Effective coefficients.** Erodibilities and diffusivities are effective values for
  ~150-km cells. At that size the cells hide all the relief that does the eroding, so the
  values are far from landscape ones.
- **Sediment.** Sediment a river carries past its first sea cell stays in that column as a
  delta, placed by isostasy the next step. Nothing leaves the planet.
- **Write-back.** Erosion and deposition change crust thickness in the plate frames.
- **Rates follow the rock.** `planet.upliftRate` and `planet.dzdt` compare each column
  with itself a step earlier, in its plate's frame. Measured at fixed map points instead,
  a mountain range sliding past would read as uplift on one flank and collapse on the other.

**Bookkeeping.**
- **Crust volume.** In the plate frames it changes by exactly created + arcs − subducted
  + surface change. The tests hold this to 1e-10.
- **Seams.** Where neighbouring world cells sample the same plate cell (a few per cent of
  cells, along seams where the rotated mesh does not line up), the map shows that column
  twice and skips another for that step. The volume itself is unaffected.

### From planet to landscape

`zoomFromPlanet(planet, { lon, lat, factor | dx, cols, rows })` lays a flat hex patch on
the planet with an azimuthal equidistant projection about (lon, lat). That projection keeps
true distances and directions from the centre. The patch then works like `zoomIn`:
- **Same machinery.** Both build a mapping from child cells to weighted parent cells, then
  share one refine step: mean-preserving elevation, roughness, edges following the
  parent's rate, and entering rivers.
- **Fields carried over.** The planet's tectonic uplift rate and rainfall become the
  patch's per-cell `uplift` and `precipitation` fields (`Simulation.setFields`).
- **Interpolation.** A Gaussian kernel 0.4 planet spacings wide, followed by 48
  mean-preserving passes. Linear interpolation on the planet's triangles left visible
  facets 100 km across.
- **Erosion parameters** are rescaled for the patch's dx: log-linear in log(dx) between
  the landscape calibration at 30 m and the planet's at 150 km
  (`regionalErosionParams`). [adapted: a rule of thumb, not a derivation]
- **Further down.** `zoomIn` refines the patch; uplift and rainfall fields pass down too.

The chain is planet (~190 km) → region (~5 km) → landscape (~600 m) → … In one run on the
`earth-like` planet, a patch centred on an active margin was compared with its frozen
planet over 2 Myr:
- The patch's mean followed the planet's extrapolation from −1597 m to +1043 m, against
  +1679 m.
- The landscape patch inside it stayed within 2–4% of its own parent's mean.

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
| results converge as dt shrinks: changes over the same span with dt, dt/2, dt/4, dt/8 differ by successively smaller amounts | `simulation.test.js` |
| the hex grid and an irregular mesh (2–6 neighbours) satisfy the mesh contract; all processes run on the irregular mesh with water draining everywhere and the volume budget closed; creep conserves volume with unequal cell areas | `mesh.test.js` |
| zoom: interpolation is exact for linear fields; the child preserves parent cell means to < 5 cm; strata, plugs and weak zones line up across scales; the largest river leaves the patch with the parent's discharge to within 10% (without inflows it carries about half); edges follow the parent's rate exactly | `nest.test.js` |
| zoom: a same-resolution patch matches the parent run on at the same dt to within 1% of relief, and its volume budget (uplift and edge motion + inflow − export) closes | `nest.test.js` |
| sphere: 10f²+2 cells with 12 pentagons; areas tile 4πR² to 1e-12; the Laplacian of l = 1 and l = 2 harmonics converges at second order (error < 0.5% at f = 16); greedy nearest-cell search equals brute force; interpolation reproduces linear fields; the azimuthal projection keeps distances from its centre; surface processes on the sphere close the volume budget | `sphere.test.js` |
| earthquakes: maximum-likelihood b (Aki 1965) recovers the Gutenberg–Richter b to 3%; over 1 Myr (~5600 events) the mean throw matches slip rate × sin(dip) × time to 5% (measured 1.3%); a fault offsets the strata; a fresh scarp under linear creep matches the error-function profile H/2·(1 + erf(d / 2√(Dt))) (Hanks et al. 1984) to 1% of its height (measured 0.04%); shaking fells cliffs that stand without it; runs are deterministic through snapshots | `earthquakes.test.js` |
| planet: isostasy gives the reference continent, the ridge depth, √age subsidence and 0.15 m of rebound per metre eroded; crust volume closes to 1e-10 of the total; a single rotating plate makes and destroys no crust and keeps its continent's area within 3% while turning by ωt; two plates spreading for 40 Myr give an age–distance slope matching the 4 cm/yr half rate to 10% (measured 4.07), and a trench where they converge; deterministic | `planet.test.js` |
| planet → patch: block means kept (< 0.5 m), sea level equal, uplift and rainfall fields carried over | `planet.test.js` |

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
- **Time step.** Operator splitting is first-order, and results converge as dt shrinks,
  but they are not dt-independent. On the alpine preset, the mean elevation change over
  150 kyr was −0.9, −6.5, −11.1 and −14.4 m at dt = 2000, 1000, 500 and 250 yr. Coarse
  steps therefore bias the result by several metres per 100 kyr in fast-eroding
  landscapes. Rockfall is stochastic per step: when dt is much larger than the failure
  timescale, every unstable face fails once per step.
- **Zoomed patches** are one-way coupled and their parent is frozen at its last rates.
  Refining a patch changes its relief and erosion rate (resolution dependence), so a
  refined patch drifts from what its parent predicts: about +27 m over 150 kyr in the
  alpine example above.
- **Planet.**
  - **Plates are kinematic.** They obey no forces; only collisions slow them. No new
    plates are born except through reorganisation.
  - **Isostasy is local.** Airy, with no flexure, so a mountain belt's foreland basin is
    missing.
  - **Nothing recycles continental crust** except subducted sediment, while arcs and
    sediment wedges add to it. With the defaults, continental area grows from 35% to about
    45% and land shrinks from about 29% to about 20% over 200 Myr. Sea level rises a few
    hundred metres as the basins fill.
  - **Calibration.** The planet's erosion coefficients are tuned only so that continents
    neither vanish nor grow without limit over ~10⁸ yr. They are not calibrated against
    real denudation rates.
  - **Cell size.** At 150–200 km cells, a "mountain range" is one or two cells wide.
- **Earthquakes** move the ground vertically only. There are no strike-slip offsets of
  channels, no aftershocks or clustering, and faults are fixed straight lines. A planet
  patch does not yet get faults from the plate boundaries crossing it: its tectonics
  arrive as a smooth uplift field.
- **Not modelled:** explicit water depth or flood hydraulics; groundwater; chemical
  dissolution (karst); glaciers; wind; isostasy at landscape scale; lateral channel migration and meanders;
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
