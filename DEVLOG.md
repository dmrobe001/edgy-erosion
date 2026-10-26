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

---

## Things that turned out not to matter, or were left alone

- **The first L1 fault-trace algorithm** (strict kind matching) mostly needed looser grouping, not a different tracer.
- **The camera spin** was not the cause of the perceived global drift; the drift was real (see v3).
- **Raising the arc-accretion threshold and lowering collapse probability** did little to the early crust-budget creep; the causes were the resampling bugs above. The thermostat remains, but now has much less to compensate for.
