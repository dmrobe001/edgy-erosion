# Development log

How L0 got to its current form: what was tried, what broke, how it was diagnosed, and what fixed it. Snapshots of the simulation core at each milestone are in `history/`. Numbers are from headless Node.js runs unless stated otherwise.

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

---

## Things that turned out not to matter, or were left alone

- **The first L1 fault-trace algorithm** (strict kind matching) mostly needed looser grouping, not a different tracer.
- **The camera spin** was not the cause of the perceived global drift; the drift was real (see v3).
- **Raising the arc-accretion threshold and lowering collapse probability** did little to the early crust-budget creep; the causes were the resampling bugs above. The thermostat remains, but now has much less to compensate for.
