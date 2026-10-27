# Development log

How L0 and L1 got to their current form: what was tried, what broke, how it was diagnosed, and what fixed it. Snapshots of the simulation core at each milestone are in `history/`. Numbers are from headless Node.js runs unless stated otherwise.

A recurring lesson: **measure before tuning.** Several problems had a different cause than the first guess, and two "obvious" fixes would not have worked. Most diagnoses came from instrumenting a copy of the simulation and budgeting a quantity (crust area, crust volume, flicker events) by cause.

---

## v1: flat torus prototype (`history/v1_flat_torus_sim.js`)

**What:** a 2D grid wrapping in both directions, plates as rigid translations with sub-cell offsets that shift one whole cell at a time (the PlaTec approach), Airy isostasy, square-root-of-age seafloor, hotspots, diffusion and erosion, random rifts and collision-counter sutures.

**Worked:** fast (2.5–9 ms per step), and supercontinent cycles emerged without scripting.

**Problems and fixes:**
- **Continental crust decayed** from 30% to 12–16% over 1.5 Gyr because collisions destroy area. Keeping more of the loser's thickness, letting thick arcs accrete as continent, and spreading very thick crust onto adjacent seafloor ("orogenic collapse") stabilised it.
- **Seafloor was far too young:** 55% under 20 Myr, against Earth's mean age of about 60 Myr. New crust was being assigned to a random neighbouring plate, which made ragged boundaries that churned crust every step. Assigning gap crust to the plate whose trailing edge it is, plus slower plates and fewer of them, brought the mean age to 60–140 Myr.
- **Never fixed:** comb-like streaks in seafloor from whole-cell jumps, inherent to the approach.

## v2: sphere with rotating plates (`history/v2_sphere_rotating_plates_sim.js`)

**What:** cube-sphere grid with constant-time lookup; each plate stores crust in its own body frame and rotates about an Euler pole; each fixed cell asks candidate plates whether they cover it. Plates gained rock types (surface unit plus formation time).

**Problems and fixes:**
- **Continents grew to 90% of the planet.** Root cause: nearest-cell resampling between equal-resolution rotated grids is not one-to-one. When one fixed cell subducted a body cell, another fixed cell that also sampled it saw a gap and refilled it as fresh crust. The fix was to mark subducted cells with negative thickness and *restore* them if still sampled elsewhere that step.
- **Still creeping upward.** Instrumenting gains by mechanism showed hole-filling copied continental crust into previously oceanic cells about twice as often as the reverse. Preferring a source whose crust type matches what was there last step reduced the bias.
- **A related bug:** a fixed cell covered early in a step could keep pointing at a body cell cleared later the same step, then resurrect it with clamped thickness. Fixed by restoring any such cell after coverage.
- **Remaining drift:** handled with a **crust-budget thermostat** that speeds margin growth below target and founders drowned margins above it. It was deliberately labelled a fudge in the page.
- **Rock-type noise:** arc rock was tagged on every subduction contact, which painted the ocean with arc speckle. Tagging only on substantial build-up fixed it.

## v3: rock columns, fault records, net rotation (`history/v3_rock_columns_faults_netrotation_sim.js`)

**What:** every crust parcel points into a pooled record holding a rock stack (push on deposition, pop on erosion, reconcile thickness changes at depth into the basement) and an optional fault record (two episodes of across and along slip, strike stored in the plate frame). Records are cloned when crust is copied and garbage-collected when orphaned. Collisions thicken in proportion to the across-boundary share of motion. Handoff extraction traces faults into polylines and exports JSON.

**Design decision that worked:** attaching rock and fault history to the *parcel* rather than to fixed cells, so history moves, splits and welds with crust for free. Per-cell arrays on every plate would have cost about 40–80 MB; the pool keeps it proportional to live crust.

**Extraction tuning:** grouping fault cells by identical kind gave about 1,200 fragments averaging 2.7 points. Grouping by strike compatibility alone, classifying from summed slip, and requiring 60 km of slip gave about 260 faults averaging 800 km.

**Net rotation ("flying ice cube"):** the user noticed a coherent global drift. Measurement showed net lithospheric rotation of 0.1 °/Myr from the start (ten random rotation vectors don't sum to zero), growing to 0.3–0.45 °/Myr, comparable to the plates' own motion. A second, related creep: an unconstrained 3D random walk of Ω drifts outward in magnitude, so plate speeds crept toward the cap.

- **Fix for speed:** relax each plate's speed toward its own preferred value.
- **Fix for net rotation:** each step, compute net rotation with the exact discrete area-weighted tensor and subtract it from every plate. That zeroes it to machine precision (about 10⁻¹⁷) and leaves relative motion untouched, which was verified.
- **Kept as an option:** a toggle, because Earth does have a small net rotation relative to hotspots.

## v4: formations, domains, L1 preview (`history/v4_formations_domains_sim.js`)

**Problem:** displayed rock types flickered from step to step (9% of continental cells changed per step; 24% of neighbouring continental cells disagreed), and coasts were pixelated. Both would turn into checkerboards one level down.

**Diagnosis by decomposition:** 56% of flips came from surface processes on the same parcel, almost all thin clastic veneers crossing the 50 m display threshold. Sediment was routed to the single steepest neighbour, which changes step to step. The other 44% came from resampling as plates rotate.

**What worked:**
- Multiple-flow deposition, spreading sediment over all lower neighbours by drop.
- Displaying the dominant rock in the top kilometre instead of the topmost unit.
- Together these cut surface-process flicker from about 5% to under 1% of cells per step.
- **The real fix was at the handoff, not L0:** pass band-limited fields. Rock units are correlated into **formations** across neighbouring columns (same rock and origin plate, compatible age, order preserved), and columns are grouped into **domains** that interpolation never crosses. Formation thickness, not rock labels, is what gets smoothed and interpolated. That gives smooth contacts, pinch-outs, and rock boundaries that follow terrain at L1.
- Absorbing sliver domains under 8 cells reduced about 2,100 raw domains to about 320.

**What was not done, and why:** interpolating when sampling plate grids. Interpolating then writing back acts as numerical diffusion. Removing resampling noise properly needs a finer body grid or particles.

## v5: stale cells and pinholes (`history/v5_stale_cell_fixes_sim.js`)

**Problem:** after sutures, continents sprouted single-cell ocean "pinholes".

**The first guess was wrong.** It was "new ridge crust created in holes", and a fix for that changed nothing measurable. Profiling the pinholes showed 85% were seafloor older than 60 Myr, only 12% touched a plate boundary, and half vanished within three steps. They were *stale* plate-grid cells: cells no fixed cell sampled in a step miss that step's surface processes, so old seafloor survived inside continents that grew around it, and blinked into view when rotation brought it under a fixed cell. Welding made it worse, because pushing absorbed cells to their nearest targets is many-to-one and left holes, while the surviving plate's grid held stale leftovers in the region it had been subducting under.

**Fixes:**
- Welding *pulls* crust into every surviving-plate cell on the absorbed plate's territory, overwriting stale leftovers.
- Each step, unsampled cells surrounded by sampled neighbours take their crust type and mean thickness.
- A gap where the departing plate still has crust on three or more grid neighbours is treated as a hole, not as new seafloor.

**Results:**
- Pinholes fell from 16–22 to about 10 per 1,000 continental cells, and from 35–42 to about 14 inside freshly absorbed plates.
- Most pinholes still present 10 steps later disappeared (41 of 452 before, 3 of 146 after).
- Bonus: continental area stopped drifting (26–29% over 2.4 Gyr, previously 28–39%), so much of the earlier "budget" drift had been staleness too.

## v6 (current, `src/sim.js`): drowning continents

**Problem:** land above sea shrank steadily. The user proposed starting with more land.

**That alone would not have worked.** Starting at 40% continental crust still took land from 40% to 5% over 1.6 Gyr. Continental *area* held (the thermostat sees to that) but mean continental *thickness* fell from 41 to 23 km. The continents were drowning, not shrinking.

**Volume budget by mechanism:** over 400 Myr, continent–continent collisions destroyed about 148k cell-km of crust and kept only about 53k as thickening. Erosion was nearly volume-neutral (about 95% redeposited on land). Collisions deleting most of the consumed crust contradicted the geology the sim was meant to represent: buoyant continental crust underthrusts and stacks up rather than vanishing.

**Fixes:**
- Collisions keep all consumed crust, spread inland for head-on motion and at the boundary for sideways motion.
- Continental slices pushed onto seafloor accrete as continent; oceanic thickness is capped at 30 km. The first attempt produced 19–75 km "elevations" from oceanic cells receiving collision crust.
- A **Continental crust** slider sets both the starting fraction and the thermostat target (a higher start alone would just be pulled back down).

**Results:**
- With a 50% target, land levels off at about 22–34% instead of trending to zero.
- Peaks stay at 7–8 km, and continental crust holds within about 4 points of target.

**Known follow-up:** the volume test now shows the 75 km thickness cap as the dominant sink (about 280k cell-km deleted over 400 Myr), with continental volume still down about 17% over that window. Routing excess sideways instead of deleting it is the obvious next step.

## L1 v1 (`src/l1.js`): one continent, coupled to L0

**Design decision:** instead of holding a refined continent to L0's final elevation, start from L0 at t0 and run both forward, L0 supplying only tectonic thickening and L1 replacing L0's surface processes. Consistency between levels became a crust-volume budget rather than an elevation constraint. Rationale and time scales are in README section 5.

**L0 changes, all pure recording:** an optional per-step tectonic-thickening record (arcs, collisions, hotspots, collapse), a `parent` on split-off plates and `mergedInto` on welded ones, and elevation pulled out as a pure function `elevOf` shared by both levels. A run with the recorder on is bit-identical to the previous commit's simulation (checked by hashing elevation and plate ownership at 500 Myr).

**What went wrong first, and the fixes:**
- **22 km peaks along consumed parcels.** Land cells on the domain edge were treated as fixed base level, so they never eroded while collision forcing piled crust onto them. A land edge is now a closed boundary (a drainage divide).
- **Spikes up to 55 km after flexure.** Elevation came from the total thickness anomaly, with ocean cells entered as a "negative crust" load, which boosted land next to deep ocean as soon as coasts sharpened. Flexure now acts only on changes since t0; the inherited state counts as already in equilibrium.
- **40 km piles on narrow arc peninsulas.** Linear elastic support held up ~140 km of accumulated L0 thickening on a strip two cells wide, surrounded by static ocean. Three changes: thickness-dependent lower-crustal flow (∝ T³), a 3 km strength limit on departures from local balance, and applying the 75 km cap *before* computing elevation rather than after (the cap had been running after elevation within the step). Peaks now stay under ~10 km.
- **Calibrating K:** see the README table. K = 1 (low end of field estimates for m = 0.5) flattened the continent to near sea level within 75 Myr and left L1's crust 4–6 km thinner than L0's. K = 0.1 tracks L0 and gives tens of m/Myr of mean continental erosion. This anchors L1 to L0's own crude erosion rate, which is itself a choice; it is exposed as a slider.
- **Parallel rivers.** D8 routing on a bilinear surface draws parallel streams: the inherited slope drops ~13 m per cell, against 1 m of symmetry-breaking jitter. Fractal relief in L1's own band (from 1.5 L0 cells down to 2 L1 cells) gave dendritic networks; Hack's exponent stayed at 0.53–0.54.
- **A checkerboard of tiny seas.** A parcel momentarily absent from its plate was searched for on every descendant plate, including ones split off long before t0, whose grids had refilled the same body index with unrelated crust. Parcels jumped plates and every jump opened a "rift". Splits are now caught in the step they happen. (A majority filter tried first could not fix a checkerboard: it just inverted it.)
- **Terraced shading when zoomed.** Pixel-space differences over nearest-sampled cells step at every cell edge. Shading is now per cell, from a least-squares plane through the land neighbours.

- **A ring of islands along the shelf edge** (reported from the page: hundreds of metres to 2 km tall, with no tectonic input). Two guesses were wrong: the fractal noise (±100 m at most) and arc input (none there). Instrumenting perimeter cells showed lower-crustal flow had added 12–26 km of crust to them: it diffused *total* thickness, so the thin inherited rifted margin read as a gradient to flatten, crust was pumped outward, piled against the static ocean (no flux across), and the narrow strip was then held up flexurally. Flow now diffuses the thickness change since t0, as flexure already did; emerged shelf cells fell from 2,244 to 139, all but 2 of them under real tectonic input.
- **Lattice-aligned rivers.** Node jitter (above) cut the share of channel length in straight runs of 4+ identical steps from 53% to 29%.

- **Suture lookahead.** Three problems on the way. (1) Tracing a welding continent back to t0 through split ancestry alone failed when the crust had changed plates through an earlier weld (plate 20 welded into 14, then 22 split off 14 carrying former plate-20 crust); L0 now optionally records merge provenance and contacts are traced back through both splits and welds. (2) Two partner continents on one plate became two pieces that "welded" at +1 Myr with no links; same-plate continents are now one piece. (3) One frame per piece missed the real weld when the partner rifted first and only its new plate welded on; pieces now keep a frame per plate they occupy, and the map draws each rigid part where it is.

- **Continents that affect each other without welding** (reported from the page: square edges at odd angles, land that should be colliding but wasn't, grey patches eating land). Causes: each continent lived in its own plate's frame and only L0, at 100 km, said what happened where plates met; consumed L0 parcels removed whole 100 km blocks of L1 crust; pieces could overlap with no interaction; and only continents that would *weld* were included. Redesign, in tested layers: L0 parcel identities (recording only); a relevance pass on the lookahead; pieces brought in with a 15 Myr warm-up and dropped once separated; a world-space layer that resolves overlaps and contacts at L1 resolution; world-space drawing over L0's own elevation. What went wrong on the way:
  - Chaining relevance back through time ("within reach of whatever is relevant next step") compounds the reach every step: a 451-cell continent pulled in 11,133 parcels from ten plates. Relevance is now just "within reach of the target through contiguous land at this step"; 800 parcels from three plates.
  - Float32 storage of times put the last frame's time a hair below itself, so everything relevant to the end was dropped on the final step.
  - The world-space layer first indexed and searched every cell (510 ms per step). Working parcel-first (one lookup per parcel to find which parcels are near another part, or near another plate in L0) brought it to ~85 ms.
  - A dark rectangle inside the continent was a hole: crust consumed while the colliding plate overrode it in L0, later refilled in L0 by its hole filling but never in L1. L1 now refills holes too; a first version churned at moving collision fronts (refilled, overrun next step), so a hole must persist 5 Myr before it is refilled.
  - Input concentrated in one L0 cell interpolated into a square pyramid; it is now spread twice before interpolation.

**Speed:** the first version took ~430 ms per coupled step at 483k cells. Per-cell closures in the coarse↔fine transfers (flexure, crustal flow) were replaced by precomputed stencils, crustal-flow edge coefficients are computed once per step, and Priority-Flood uses the Zhou et al. (2016) variant (bit-identical filled surface, receivers and areas compared with the plain algorithm). Now ~250 ms in Chromium.

## L1 v2 (`src/l1.js`): one landmass, a thin viscous sheet, drainage basins

**Design change.** v1 conflated a continent with a landmass, and spent most of its effort tracking which crust mattered as rigid pieces met. The new unit is one landmass at the present, bounded by its coastline at the end of L1. The reasoning behind the change:
- Rivers mature within a few Myr, so fine erosion over a long window bought little.
- What L0 cannot do is deform crust where plates meet.

So L1 now:
1. picks a landmass at the present;
2. replays its last ~20 Myr from a saved L0 state, with the crust as a thin viscous sheet driven by L0's plates;
3. freezes tectonics and runs 5 Myr of fine erosion only to settle drainage basins, the tiles for L2.

The world-space overlap layer, lookahead relevance, pieces and refined cube-sphere grid of v1 are gone; picking at the end with replay answers "which crust" directly.

**L0 changes, all recording only (bit-identical runs, checked by hash):**
- `snapshot()` and `restore()`, with the generator's state exposed;
- `trackParcels()` and `untrackParcels()`;
- the tectonic recorder now records only arcs and hotspots (see below).

**Why L0's collision ranges seemed to fade** (measured before starting): parcels that rise above 4 km hold 4.4–4.7 km for 30 Myr. Peaks above 6 km relax to about 5.9 km within 5 Myr and to 5.1 km by 30 Myr. About 20% of the highest parcels are consumed or relabelled within that time, and diffusion contributes little. The "flash" of white was the top band of the colour map (white from 5.5 km), not missing physics. L0 was left alone.

**What went wrong on the way:**
- **Max thickness 113 km and 15 km peaks in the first run.** Grid thickness came from deposited volume ÷ area, which underestimates thickness where markers thin out at the crust's edge, while the reference thickness was an average. The mismatch gave spurious buoyancy at every edge. Markers now carry their own thickness and area, both strained by exp(∓div·dt), so each marker conserves its volume exactly and grid fields are plain averages.
- **A one-step spike to 431 km** (31 km elevation). L0's orogenic collapse turns an ocean cell into 28 km of continent in one step. Passed on as "volcanism" and shared by volume among the few markers in a sheet cell, it multiplied their thickness. Collisions and collapse are now not passed on (the sheet makes its own), and arcs and hotspots thicken each marker by what its L0 cell gained.
- **Conjugate gradients never converging** at 30 km cells (residual 10⁻² after 400 iterations; 9–10 s per step). At 50 km, 200 iterations land within tens of metres of a converged solve at a tenth of the cost (README 5.4). A multigrid preconditioner is the real fix.
- **Holes and a gap in the crust at t0.** Pinholes were parcels L0 relabelled during the window. The gap was crust overrun in the collision: absent at the end, so not traced. That crust is exactly what a crumple zone is made of. Consumed crust that touches the landmass during the window is now included, and pinholes are filled. After the fix, the two collision fronts showed parallel thickening bands of 5–14 km.
- **Grid-aligned gullies.** Stream power on the 50 km sheet grid cut 50 km gullies along the grid axes, and the drainage pass inherited them as parallel stripe-shaped basins. During the sheet, erosion is now a smooth function of elevation (calibrated to Ahnert's magnitudes). A first linear version (elevation ÷ 10 Myr) stripped 2–4 km from lowlands, because isostatic rebound keeps a plain's elevation up while it erodes. The fine drainage pass now draws all the rivers, with stronger fixed relief to organise around, and the basins come out dendritic.
- **Lattice-aligned rivers again, in the drainage pass.** v1's node jitter was ported to the planar grid. Straight runs of 4+ identical steps went from 52% to 33% of channel steps (v1: 53% to 29%), and Hack's exponent stayed at 0.53–0.55.
- **Are the thickening bands real?** They were about 300 km apart, close to both a few sheet cells and a few L0 cells. Running the same landmass at sheet grids of 70, 50 and 30 km, and with jittered markers, gave the same spacing (290–330 km by autocorrelation across the bands). Only the band width scaled with the grid (about 1.5 cells). Tracing markers back to t0 showed one band on the suture and the other at a promontory of the indenting margin, so the spacing is inherited coastline geometry. Marker seeding is now stratified-random anyway: under shortening, regular marker rows bunched into visible lines.
- **Closed basins, three attempts.** A water balance on every depression the filling found turned the dry belt into 800 tiny closed basins, one per hollow in the fine relief. Filling hollows with the sediment their catchments shed during the 5 Myr pass removed those, but also filled a 600 m deep, 700 km test bowl in a desert. That pass has tectonics frozen, while real closed basins persist because tectonics keeps them subsiding. The rule that held: only depressions of the resolved surface count (the sheet's surface smoothed over one cell, since nothing narrower than about three cells is resolved), and other hollows are filled as stand-in texture. The same bowl test also showed the first climate was too dry inland and too evaporative in the tropics: a 1,500 km rain-out length ignored rain recycled by vegetation, and evaporation followed temperature alone, though it peaks under clear subtropical skies.
- **Erodibility at L1: tested and left to L2.** The decision gate (README section 6, item 3) found rock type, age and fault damage, even at a 1,600× range, moving the basin map far less than re-seeding the synthetic fine relief. The suspicion that the relief was too loud did not hold up: it continues the sheet's spectrum.
- **Inland drainage.** Four steps to a result that holds up by aridity band:
  1. Channel losses proportional to channel width (√Q) could not both kill small desert rivers and let large ones cross deserts. A loss per km independent of size, as flood pulses soaking into beds give, does both.
  2. Classifying whole coastal basins by whether any water reaches the sea counted desert hillslopes as draining to the sea when wetter headwaters fed the same basin.
  3. Classifying per cell, by whether the cell's own flow path keeps carrying water to the coast, separates areic desert. Arid land then drains 88–89% inland on two landmasses, semi-arid 19–25%, humid none.
  4. With present tectonic rates in the pass, the sediment-fill rule is consistent again (it failed while tectonics was frozen) and replaces the resolved-depression rule.
- **Dynamic topography.**
  - **Zero mean.** The first field only pulled down, so every continent sank by about 200 m on average; it now averages zero, the return flow of the slabs lifting the rest.
  - **Kernel check.** Against the designed Gaussian it came out too wide (0.69 instead of 0.61 of the depth at one width) for two reasons. The test's grid (417 km cells) was coarser than its 380 km kernel, so the blur grid now has at least 48 cells a face. And full-weight neighbour averaging checkerboards over the few steps a narrow kernel needs; lazy half-steps fixed it (0.603 against 0.607).
  - **Robustness.** It barely made basins less contingent, which I had predicted. Measuring that properly needed a second metric: outlet shift per land cell, since label agreement penalises every moved divide.
- **The L2 handoff.**
  - **Tiles.** Basins under 10,000 km² first became tiles one by one: 152 tiles on a small landmass, many of them 100 km² specks. Grouping minor basins breadth first left small closed basins boxed in by major ones. These now join the neighbouring tile they share the most boundary with; any union of whole basins is still closed to water.
  - **Which divides are real.** A first test counted a divide as surviving a re-seed only if the exact 10 km edge stayed a divide. That kept 12–40%, with high crests faring worst, because the synthetic relief's amplitude grows with elevation and jitters crests by a cell. Distance to the nearest divide afterwards flattered every place with dense divides, lowlands full of tiny basins and escarpments full of parallel streams alike: rugged random cells were within 50 km of one 77–92% of the time. The test that held compares point pairs 30 km either side of each divide with random pairs as far apart. Persistence then rises steadily with local relief, from 45% (random 25%) to 90% (random 62%), with no sharp threshold. Three other predictors did worse: crest height, the convexity of the sheet's smooth surface, and divides of the drainage routed over that smooth surface alone, which were worse than chance because its near-flat lowlands are full of arbitrary divides. The firm/contingent flag (relief above the synthetic relief's peak-to-peak amplitude: 77% against 52%) is a convenience; the table travels with it.

- **What counts as the landmass.** Refined land sometimes ended in a straight edge against a lump of unrefined land. Measured on three 1 Gyr worlds before changing anything, the pick missed four things:
  - land touching only at a corner: 18–68 such contacts per world;
  - islands one or two sea cells off the coast;
  - enclosed seas: 13–55 cells on the largest landmass, left out of the pick because they are water;
  - land on oceanic crust: 1–3% of the land, never carried by the sheet, which seeded continental crust only.

  The pick now joins corner contacts, chains islands within a buffer of the coast (the buffer growing around each island that joins, so it is single-linkage clustering at the buffer distance), takes in everything the region encloses, and marks land on oceanic crust for the sheet to carry.

  Distance was first measured between cell centres, so a one-cell strait (about 100 km of water) needed a 210 km buffer to bridge. It is now coast to coast, the centre distance less a cell.

  After the change, all but 16 of 4,889 land cells end under the sheet's crust (seed b, 300 km buffer). The 16 carry markers, but the sheet moved the coast a cell away from L0's.

---

## Things that turned out not to matter, or were left alone

- **The first L1 fault-trace algorithm** (strict kind matching) mostly needed looser grouping, not a different tracer.
- **The camera spin** was not the cause of the perceived global drift; the drift was real (see v3).
- **Raising the arc-accretion threshold and lowering collapse probability** did little to the early crust-budget creep; the causes were the resampling bugs above. The thermostat remains, but now has much less to compensate for.

- **Square mountains and clipped coasts.** On a 1 Gyr landmass (seed genesis, 150 km buffer, 20 Myr window), the high ground of the L1 map came out in rectangles, and the coast had corners cut out of it.
  - **Rectangles.** I took the sheet's elevation apart, term by term, on the north of that landmass. The inherited surface was smooth. Crustal thickening from strain had the rectangles, up to 8 km of it, and so did strength. Cratons are 10 times stronger, and craton membership was taken parcel by parcel, so strength changed in steps along L0's 100 km cell edges. The weak crust between cratons took the shortening along those staircases. The craton share is now smoothed on L0's grid before markers sample it. Peak elevation on that landmass fell from 6.3 to 4.8 km, the same convergence spread over wider belts.
  - **Erosion and volcanism were blocky too, but that was not the visible cause.** Erosion took the same depth from every marker in a 50 km cell, and markers in one cell came to differ by up to 13 km in thickness. Volcanism added crust by L0 cell. Both are now interpolated bilinearly to each marker. Fixing these first left the rectangles as they were.
  - **Clipped coasts.** Three things, measured on the 10 km grid against L0's interpolated coast.
    - The sheet carried only continental crust, so the edge of its crust cut the coast where L0's land runs into an oceanic cell. It now also carries a coastal rim of oceanic cells on the land's own plate.
    - Crust still oceanic at the window's start was not seeded, so accreted crust was missing. Every parcel of the object is seeded now.
    - The inherited surface was smoothed across the coast, sinking coastal lowland toward the abyss next to it. Land is now smoothed with land and sea with sea.
    
    Land on L0's pick with no sheet crust fell from 728 to 492 fine cells. All 7 L0 cells left are parcels the sheet did carry and moved elsewhere. Of the 8,200 fine cells where L1 has crust under the sea but L0 has land, 2,300 are drowned by the mantle pull, which L1 adds and L0 does not have. The rest are drowned by erosion and spreading. Both are L1's dynamics, not gaps in the pick.
  - **No more L0 behind L1.** Drawing L0's surface outside the sheet's crust showed two shorelines that disagree for those reasons. The refined map now draws plain ocean there. The shelf fades into it across the edge cell, because a hard edge showed the crust's 50 km cells as a staircase in the sea.

- **Crosshatch on the refined map.** Zoomed in, mountains showed a grid of horizontal and vertical streaks at 10 km. The map drew each cell of the drainage pass as a flat block, and the hillshade lit only the blocks' edges. Land is now drawn bilinear between the pass's land cells (land only, so the coast does not move), as the L0 → L1 handoff interpolates L0's surface. The coast and rivers are still drawn cell by cell, and show 10 km steps when zoomed in.

- **Selecting basins for L2.** A buffer about the tapped basin, as at L1, would drag in a large desert basin behind a cluster of small coastal ones. The rule is instead: the tapped basin whole, whatever its size; the basins wholly within a radius of the tap; and what the selection encloses. A synthetic map with a basin inside a ring basin tests the enclosure rule, since the test landmasses had none to enclose. dx and dt now show in every panel: L0 104 km / 0.91 Myr at n = 96; the sheet 50 km / L0's step over its substeps; the drainage pass 10 km / 1 Myr.

- **L2's relief, a function of position.** Built so that basin boundaries can later be burnt in strip by strip. What the tests turned up:
  - *Leaking into L1's band.* Averaged over L1's 10 km cells, the first version differed from L1's surface by 53 m RMS on the genesis landmass, half of L1's own cell-to-cell relief. Taken apart: 29 m from the bicubic base, because an interpolant's cell mean is not its value at the cell centre, and 45 m from the added relief, because noise of lattice spacing 16 km has power at wavelengths L1 resolves. Correcting the base to preserve cell means and subtracting the added relief's own interpolated cell means brought it to 12 m (4 m on the test landmass). The high-pass removes most of the 16 km octave, so the test checks the amplitude law on the octaves of 4 km and less, where the band is L2's alone.
  - *Burnt-in rivers.* A trench whose sides rise at L1's valley-side slope flattened wedges kilometres wide; sides that rise by 3 RMS of the added relief per km fixed that. Burnt along L1's reaches, rivers were straight canals at 0° and 45°; warping the coordinates by a smooth 1.5 km noise lets them wander within L1's cell. The last reach dived toward L1's −1 km sea floor; it is now held at sea level. A V-shaped trench still let a 500 m grid dam rivers by up to 58 m, as the grid may have no cell on the river's line; a floor 1.2 km wide removed every dam.
  - *Plaid.* At high zoom, octaves sharing a lattice's axes added up to a grid of lines; each octave is now turned by its own angle.

- **Burning in basin outlines at 500 m.** Membership is where a point's water leaves the land, found by priority flood; the reference is the same flood over the selection's whole box. Labels come from canonical tiles, so neighbouring selections agree at their seams. What didn't work, in order:
  - *Halo from the steepest slope.* The first rule sized a tile's halo by the added relief's RMS over L1's steepest drop to a neighbour. On the coast it worked once every tile had a halo of at least two cells and a two-cell band was refined (mean displacement from 2.5 km for L1's outline to 0.02 km). On a gentle inland divide it stayed at two cells, and the outline was off by 1.9 km. The steepest drop is usually the wall of a small valley L1 carved. The regional slope (L1's surface blurred over 20 km) changed nothing either.
  - *Closing the halo's edge near divides.* The halo's edge took L1's labels, and the divide runs through the halo too, so water running along a gentle divide took a coin-flip label there. Seeding only the edge where L1's label is sure made things worse (coast 0.87 km): closed in, the water spilled somewhere arbitrary.
  - *Halos that agree with each other.* Widening a tile's halo until its labels stopped changing settled most inland tiles at 4 cells, because halos of 2 and 4 cells sent the same water to the same uncertain edge. Agreement proved nothing.
  - *What worked:* tracking where each label came from. A core cell's label is a guess only if its water reached the halo's edge near an L1 divide; the halo doubles until no label is a guess (to 24 cells). Every case tested now matches full routing exactly. The price is wide halos in flat country: 20 s for an 18,000 km² lowland basin.

- **A burn-in that never started.** On a landmass with a 982,800-cell drainage grid, listing the tiles still to route took 112 ms, longer than the 60 ms slice of a frame. The loop listed them, found its time gone and routed nothing, frame after frame ("0 tiles routed"). It now lists them once per phase and routes at least one tile per frame. The drawing loop had the same flaw: on a frame already over its budget it drew no row. The UI smoke test now burns in and draws under a clock that jumps 100 ms at every reading, which fails on the old code. That selection (7 basins, 301,500 km² at L1; 293,900 km² at 500 m) takes 300 tiles, 21.5 s in Node and 102 s with the CPU throttled four times.
- **Stitching instead of halos.** Each tile is now evaluated once and reduced to a graph between its edge points and outlets (Barnes's parallel priority-flood, adapted from levels to labels); the tiles' graphs are joined and flooded, and each tile flooded again from its edges. The first version matched full routing to within 0.015 km, but not exactly: all 246 mismatched points sat on tile edges at their own height, reached at that level from several lower neighbours. Routing the whole domain gives such a point its lowest neighbour's label; the stitched graph had picked a side by its own order. Edge points now take their lowest neighbour's label and the tiles are flooded again until nothing changes; two of the three cases became exact. The stitched domain depends on the selection, so two neighbouring selections now disagree on 2 of 287,000 shared points, where the halo version agreed exactly. Growing the domain more generously took fewer rounds but more tiles and cost more.
- **Web Workers.** Pass A went to workers built from the module's own source (`L2module` is a function so its source can be handed over; L1's heap became a class for the same reason). The first measurement showed no gain: each worker spends about a second building the relief from about 40 MB of copied inputs before its first tile, as long as the main thread took for all the tiles. Starting the workers when a basin is tapped hides that; then 3 workers halve the tile time. The stitching (about 2 s on the largest selection) is now the larger part, and stays on the main thread.
- **The L2 pass.** Landscape evolution at 500 m inside the burnt-in outline, L1's drainage pass one level down, with erodibility from rock now that the experiment of section 6 handed it to L2. What changed while building it:
  - *Lakes everywhere.* The first run left 7% of the area as standing water: the added relief's closed hollows, filled by the priority flood but never by sediment. L1's rule (a hollow its catchment's sediment can fill in the run becomes a flat) fills nearly all of them; 46 lakes remain on the mountain selection.
  - *Memory.* Arrays over the selection's whole box held 3.4 million points for 217,000 selected ones; they now hold only the selected points and the outlets beside them (200 MB for 1.2 million points).
  - *Crosshatch, again.* Drawn point by point, the run's surface showed the 500 m grid in the hillshade, as the drainage pass's had at L1; it is now drawn bilinear between points.
  - *Strahler order* is counted on channels only (1 km² and more); counted on every point, every hillslope pixel was a stream.

