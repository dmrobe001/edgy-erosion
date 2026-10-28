# L0 on a marker mesh: specification

Status: revision 2. Stage 1 is built (`src/mesh.js`; README §3.14; results and deviations in §11.1). Later stages are proposals, written to be argued with before code exists. This document specifies a replacement for the L0 core in `src/sim.js`.

How claims are marked:

- **[derived]**: follows from geometry, conservation or a stated equation.
- **[measured]**: an Earth range from the literature, with the source named.
- **[calibrated]**: a value chosen to make the model's output match an Earth statistic, not derived.
- **[rule]**: a modelling rule standing in for physics we do not resolve.
- **[conjecture]**: expected behaviour, to be tested.

Units are km, Myr, km² and km³, with forces in N/m (per length of boundary) or N/m² (per area), and densities in kg/m³. Useful identities: 1 cm/yr = 10 km/Myr; 1 mm/yr = 1 km/Myr.

## 0. Decisions so far

| Question | Decision |
|---|---|
| Rigid plates, or a deforming sheet? | **A deforming sheet.** Rigid plates with a torque balance are an intermediate stage, and later a limit the sheet must reproduce (§5.3). |
| How plates weld | **Emergent.** Interface markers carry coupling damage that heals while relative motion across them is slow (§7.2). No counters, no timers. |
| Ocean–continent contact | Slab pull on the oceanic side; interface traction and trench suction on the continental side (§4.5). |
| Resolution | h₀ = 100 km to start. |
| Spin-up | **None.** The initial condition is shown as it is. Slabs start at zero length, so the first 10–20 Myr are a visible transient (§4.3). |
| Thickness and elevation | Every marker carries total crustal thickness and elevation explicitly, including below sea level (§3.1). |
| Plate interfaces | Markers sit exactly on interfaces and have one **side** per adjoining plate, each with its own properties: the split-node technique (§3.2). |
| Extremes | No speed caps. Typical behaviour comes from mid-range parameters; extremes have to come from configurations (§2, principle 7). |

---

## 1. Why

The current core stores each plate as a full lattice in its own frame and resamples every plate onto a fixed lattice every step. That has cost us in four ways:

1. **Resampling artefacts.** Rotating one lattice against another produces duplicates and pinholes. Three patches deal with them (README §3.3, `syncUnsampled`), and about 90% of displayed rock flicker comes from resampling (README §3.13).
2. **Plates are stored primitives with an arbitrary lineage.** Splits come from a timer and a random plane; welds from a counter (ONTOLOGY.md, findings 3 and 6).
3. **Plate motion has no cause.** Angular velocities are a random walk.
4. **Exogenous inputs have no process.** Hotspots are fixed and immortal.

The replacement keeps what works: rock columns, fault records, isostasy, ocean age from birth time, conservation of continental crust, and deterministic replay.

## 2. Principles

1. **Material is the source of truth.** Markers carry matter and history; fields are views of them (ONTOLOGY.md, rule 1).
2. **Plates are a pattern:** regions of the sheet bounded by interfaces, tracked by lineage.
3. **The triangulation is disposable.** It is used only to compute and can be rebuilt at any time, because no history lives on it.
4. **Every exogenous input has a process**, with a magnitude and a correlation time taken from Earth.
5. **Every quantity that should be conserved has a ledger**, tested every step.
6. **Determinism.** A run is a pure function of (seed, options, number of steps).
7. **Extremes come from configurations, not parameters.** Parameters are set to mid-range Earth values and calibrated against typical statistics. Rare behaviour (a plate at 15–20 cm/yr, a back-arc basin opening, a polarity flip) has to arise from rare arrangements of forces (two slabs pulling one plate, a plume head thinning its base). Nothing is capped except a safety limit, at 40 cm/yr, which logs an anomaly if reached.

---

## 3. State

### 3.1 Marker (interior)

One marker is one parcel of lithosphere, at spacing about h₀ (about 50,000 markers). Structure-of-arrays storage; a slot is reused after deletion, but an **id** never is.

| Field | Type | Meaning |
|---|---|---|
| `id` | uint32 | Permanent identity. |
| `x, y, z` | float64 | Position on the unit sphere, mantle frame. |
| `plate` | int16 | Plate this step (derived, §3.5; stored for speed). |
| `A` | float64 | Carried area, km². Changes by strain and by the ledgered operations. |
| `T` | float32 | **Total crustal thickness**, km. Stored explicitly. Crustal volume is V = A × T. |
| `e` | float32 | **Elevation**, km, relative to sea level, negative below. Stored explicitly; recomputed by isostasy, flexure and dynamic topography each step (§4.1). |
| `K` | uint8 | 1 if continental. |
| `B` | float32 | Birth time of the crust, Myr. |
| `S` | float32 | Sediment thickness, km. |
| `b` | float32 | Column buoyancy relative to the asthenosphere, kg/m² (§4.1). Derived each step; stored because forces and polarity read it. |
| `col` | int32 | Rock column in the shared pool, as now. |
| `flt` | int32 | Fault record, as now. |
| `D` | float32 | Damage, 0 to 1 (§7.1). |
| `parent` | uint32 | Id this marker was split from, if any. |

**Invariant [derived]:** crustal volume changes only through ledgered operations (§8). Strain conserves V: it changes A by e^(div·dt) and T by e^(−div·dt), exactly as L1's markers do.

### 3.2 Interface markers (split nodes)

An interface marker sits on the line where plates meet. It has one position, a few interface properties, and **one side per adjoining plate**: two on an ordinary boundary, three at a triple junction. Each side carries that plate's local properties.

This is the split-node technique used for faults in finite-element models (Melosh & Raefsky 1981). Each side is a separate vertex for its own plate's triangles, so velocity can jump across the interface while each plate's mesh stays continuous.

**Per side:**

| Field | Meaning |
|---|---|
| `plate` | The adjoining plate. |
| `T`, `e`, `K`, `B`, `S`, `b`, `col`, `flt` | That plate's material at the interface, its "frontier values". Refreshed when an interior marker is consumed into the side (§3.4). |
| velocity | That plate's velocity at this point, from the solve. |

**Per interface marker:**

| Field | Meaning |
|---|---|
| `kind` | ridge, trench, transform, collision or locked. Derived from relative velocity and polarity (§4.6). |
| `polarity` | Which side subducts (trenches only). Has hysteresis (§4.6). |
| `Dc` | Coupling damage, 0 to 1. Sets how much the interface resists slip (§7.2). |
| `slab` | The slab hanging from this point of the trench, a record: a list of segments {time entered, age at entry, crust thickness, length}. Total length L_slab, and attached or detached. |
| `quietFor` | How long relative motion has stayed below v_c (for healing, §7.2). |

**Mesh invariant [rule]:** no triangle edge may join interior markers of two different plates. Wherever one appears, an interface marker is inserted at its midpoint. So every boundary runs through interface markers.

### 3.3 Mesh

- **Built from** the Delaunay triangulation of all marker positions, interface markers counted once. On a sphere, this is the 3D convex hull of the points [derived].
- **Maintained by** local Lawson edge flips after each step's motion. Full rebuild at start, on restore, and as a fallback.
- **Each triangle belongs to one plate.** At an interface marker, it uses that plate's side.
- **Derived per step:** Voronoi areas, neighbour lists, cotangent weights for diffusion, and the finite-element matrices for the sheet (§5).

The mismatch m = A_vor − A (carried area against the space a marker occupies) is kept as a **diagnostic**. Gaps and overlaps are now handled at interfaces, so m should stay small everywhere. A large |m| flags a bookkeeping error (invariant I1).

### 3.4 What happens at each kind of interface [rule]

| Kind | Sides | Material bookkeeping |
|---|---|---|
| **Ridge** (sides separating) | The interface marker sits midway between where its two sides would go. | When the accumulated opening reaches h₀/2, a new interior oceanic marker is created on **each** side at the side's position: K = 0, B = now, T = 7 km, a fresh basalt column, damage 0. The side then resets to the axis. Spreading is symmetric by construction; asymmetric spreading would need a rule for ridge migration (open question 3). |
| **Trench** (sides converging, one subducts) | The interface marker moves with the overriding side, plus a rollback velocity (§4.5). | A subducting-side interior marker that comes within h₀/2 of the interface is **consumed**. Its area and volume are written to the interface's slab record and to the ledger; its properties become the subducting side's frontier values. |
| **Transform** (sliding past) | Sides slide along the interface. | No area change. Interface markers are re-spaced along the line (inserted or merged) to keep spacing about h₀. |
| **Collision** (converging, neither side can subduct) | Both sides held together by coupling. | Nothing is consumed. Convergence is taken up by sheet deformation, so crust thickens on both sides, and volume is conserved exactly [derived]. |
| **Locked** (healed) | Coupling damage below the weld threshold (§7.2). | The sides merge into one interior marker. The plates weld. |

### 3.5 Plates and lineage

- **A plate is** a connected set of triangles, with interfaces as its boundary. Union-find, each step.
- **Matching:** a component takes the id of the previous plate it shares the most markers with, provided it holds at least half of that plate's markers. Otherwise it gets a new id, with a `split` edge from its source. Two previous plates in one component give a `weld` edge; the larger keeps its id. **[rule]:** this is the naming convention.
- **Per plate (derived):** member list, area, continental area, centroid. A best-fit rigid rotation ω and an accumulated orientation Q, kept so fault strikes and L1's reference frames work as now.
- **The lineage graph** (nodes are plate id over a time interval; edges are continue, split, weld) is append-only and part of the snapshot.

### 3.6 Mantle state

| Item | Content |
|---|---|
| Flow modes | Real spherical-harmonic coefficients of a poloidal potential Φ and a toroidal potential Ψ, degrees 1–4. Basal mantle velocity v_m = ∇Φ + r̂ × ∇Ψ. |
| Plumes | {position, birth time, phase, flux, drift state, death time} (§6.2). |
| Slab record (global) | As now: consumed oceanic area, weighted by age, in 10 Myr bins on a fixed mantle grid, for dynamic topography. Filled from the interface slab records. |

### 3.7 Rendering grid (temporary compatibility)

The cube-sphere grid stays, as a view only: markers are interpolated onto it each step (nearest marker for categorical data, barycentric for continuous data), producing today's arrays (`owner`, `wK`, `wT`, `elev`, `wCol`, `wR`, `wRA`). The UI, the handoff export and L1 keep working during migration. Nothing is written back.

---

## 4. Forces

All forces are applied **locally**: per marker, per triangle, or per interface marker. Torques about any axis, and stresses from an irregular trench, follow from where the forces act; neither is computed as a lumped quantity.

### 4.1 Buoyancy of a column [derived]

Each column's buoyancy relative to the asthenosphere, per unit area:

> b = Δρ_c · T − ½ · α · ρ_m · ΔT · h_th(age)

| Term | Meaning | Value |
|---|---|---|
| Δρ_c · T | Light crust | Δρ_c ≈ 400 kg/m³ for oceanic crust (2,900 against 3,300), ≈ 500 for continental |
| ½ α ρ_m ΔT · h_th | Cold, dense thermal boundary layer | α = 3×10⁻⁵ /K, ρ_m = 3,300 kg/m³, ΔT = 1,300 K, so ½ α ρ_m ΔT ≈ 64 kg/m³ |
| h_th = 2.32 √(κ · age) | Half-space cooling thickness | κ = 10⁻⁶ m²/s |

**Consequences:**

- **Oceanic lithosphere with 7 km of crust becomes denser than the asthenosphere at about 11 Myr** [derived from the values above]. That agrees with the estimate of about 10 Myr in Cloos (1993) [measured]. Younger ocean resists subduction.
- **Thick oceanic crust (plateaus, hotspot tracks) and continental crust are buoyant**, and resist subduction.
- **Continental lithosphere's thermal term is taken as cancelled** by its chemically depleted, lighter mantle root (Jordan's "isopycnic" cratons [measured, as a hypothesis widely used]). So b > 0 for continents.

**Elevation** is set from the same column: isostasy as now, plus flexure, plus dynamic topography from the slab record and the mantle flow's upwelling. It is stored in `e`.

### 4.2 Gravitational potential energy (GPE) [derived]

A body force per area, f = −∇_s GPE, where GPE is the column's gravitational potential energy relative to a reference column. The same term gives **ridge push** in the oceans (GPE rises with age away from ridges, integrated about 2–3×10¹² N/m at 100 Myr [measured: Turcotte & Schubert, *Geodynamics*]) and **gravitational spreading** of thickened continents (the Φ term of L1's sheet). One force, not two.

### 4.3 Slab pull, per interface marker [derived, with one rule]

Each trench interface marker carries its own slab record, so the pull varies along the trench with the age and thickness of what has gone down **at that point**. It acts on the subducting side, perpendicular to the trench:

> F_sp = g · Σ over slab segments of L_seg · (½ α ρ_m ΔT · h_th(age_seg) + Δρ_ecl · T_seg)

| Term | Meaning |
|---|---|
| ½ α ρ_m ΔT · h_th(age_seg) | Thermal deficit of the segment, using its age when it entered the trench |
| Δρ_ecl · T_seg | Below about 60–100 km depth, oceanic crust turns to eclogite, about 200–300 kg/m³ **denser** than mantle [measured]. The crust adds to the pull once deep |
| L_seg | Length of the segment; total L_slab is capped near 700 km (the 660 km discontinuity) [rule]. The cap's material passes into the global slab record |

**Magnitude [derived from these values]:** for a 700 km slab, F_sp ≈ 3×10¹³ N/m at 10 Myr age, ≈ 4×10¹³ at 30 Myr, ≈ 7×10¹³ at 100 Myr. This matches the literature order of 10¹³–10¹⁴ N/m for the full negative buoyancy.

**Consequences:**

- Along a straight trench, an older segment pulls harder than a younger one, so the pull exerts a torque about a point on the trench as well as pulling toward it, and stresses the plate behind it. Your point, built in.
- With no spin-up, every slab starts at zero length. Pull grows as plates converge: at 5 cm/yr a slab reaches 700 km in about 14 Myr [derived]. The first 10–20 Myr are driven by GPE and mantle flow alone, and that transient is visible on purpose.

### 4.4 Slab resistance (the force the first estimate left out) [derived form, calibrated coefficient]

The first revision's sanity check gave plates at about 30 cm/yr. That came from **omitting the main resisting force**, not from taking every effect at its fast end: a slab can only sink as fast as the mantle around it lets it, and it has to bend at the trench. Most of the slab's negative buoyancy is supported by this viscous resistance, not transmitted to the plate [measured: Conrad & Hager 1999 on bending; Conrad & Lithgow-Bertelloni 2002 on how much pull is transmitted].

Per trench interface marker, opposing the subducting side's normal velocity v_n into the trench:

> F_res = C_s · v_n,  with C_s = c_geom · η_um · (1 + h_th / R_bend)³

| Term | Meaning |
|---|---|
| η_um | Upper-mantle viscosity, about 10²⁰–10²¹ Pa s [measured range] |
| c_geom | A dimensionless shape factor for a sinking sheet, of order 10–100 [calibrated] |
| (h_th / R_bend)³ term | Bending resistance: stiffer for thick (old) slabs bending at radius R_bend ≈ 200 km [measured form: Conrad & Hager 1999] |

**Balance [derived]:** a slab alone sinks at about v ≈ F_sp / C_s. With F_sp = 4×10¹³ N/m and C_s ≈ 10²²–10²³ Pa s (η_um = 10²¹ Pa s, c_geom of 10–100, a bending factor of 1–5), v ≈ 0.4–4×10⁻⁹ m/s, which is **about 1–13 cm/yr**. Earth's range comes out with mid-range parameters, and fast plates can come from configurations: a plate pulled by two slabs (India's proposed double subduction [measured: Jagoutz et al. 2015]), or a plume head thinning and lubricating its base (§6.2; also proposed for India [measured: Kumar et al. 2007]). c_geom is **[calibrated]** once, against Earth's root-mean-square (RMS) plate speed of about 4–5 cm/yr, with all other parameters at mid-range values.

### 4.5 The overriding plate [rule with measured motivation]

At an ocean–continent or ocean–ocean trench, the overriding side feels two forces:

1. **Interface traction**, equal and opposite to what the subducting side feels from the interface (Newton's third law). Its normal part is compressive, and grows with convergence speed and coupling (1 − Dc). Strong coupling squeezes the overriding plate, as in the Andes.
2. **Trench suction.** A slab sinking faster than the plate converges rolls back, and the mantle flow it drives pulls the overriding plate toward the trench:
   - Rollback speed: v_rb = max(0, v_sink − v_n), where v_sink = F_sp / C_s.
   - Suction: F_suc = c_suc · C_s · v_rb, toward the trench, on the overriding side, with c_suc **[calibrated]**, order 0.1–0.3.
   - The interface marker moves with the overriding side plus v_rb.

**Testable consequence [conjecture]:** with old slabs, suction exceeds compression, so the overriding plate stretches near the trench, damage grows (§7.1), and **back-arc basins** open (the Mariana type). With young, buoyant slabs, compression wins (the Andean, or Chilean, type). That is the classification of Uyeda & Kanamori (1979) [measured], and it should come out of the model rather than be put in.

### 4.6 Which side subducts: polarity [rule]

Decided per contiguous convergent segment of interface, not per marker, so trenches don't flicker:

| Sides in contact | Outcome |
|---|---|
| Ocean against continent | The ocean side subducts. |
| Ocean against ocean | The side whose frontier is **denser on average along the segment** (lower b, length-weighted) subducts. |
| Continent against continent | Neither: a **collision**. |
| Buoyant ocean arriving (b > 0: young ocean, a plateau, a hotspot track) at an active trench | Subduction continues while the attached slab's pull exceeds the arriving crust's buoyancy resistance; otherwise the segment **jams** and becomes a collision. |

**Hysteresis:** an established polarity reverses only if the mean buoyancy contrast has favoured the other side by more than δb **[calibrated]** continuously for τ_flip ≈ 5 Myr. Polarity reversals are real but rare (the Solomon arc after the Ontong Java plateau jammed its trench [measured]).

**Slab break-off [rule]:** a slab attached to a collision or jammed segment detaches after τ_bo ≈ 5–15 Myr [measured: break-off typically within about 10 Myr of collision]. Its pull vanishes, its material passes into the global slab record, and the region rebounds isostatically.

### 4.7 Basal drag and mantle traction

Per area, on each triangle:

> f = −C_b · k · (v − v_m)

| Term | Value |
|---|---|
| C_b = η_asthenosphere / h_asthenosphere | With η ≈ 10¹⁹–10²⁰ Pa s and h ≈ 100–200 km: C_b ≈ 10¹⁴–10¹⁵ Pa s/m [measured range] |
| k | 1 under ocean; k_keel ≈ 3–10 under thick continental roots [conjecture: the literature disagrees]; reduced above plume heads (§6.2) |

### 4.8 Interface coupling

Across each interface, a traction resists relative velocity Δv:

> τ = c_int · (1 − Dc) · Δv, capped at a frictional limit S_int · (1 − Dc)

c_int and S_int depend on kind: small for ridges and transforms, moderate for trenches, large for collisions **[calibrated]**: chosen so continental collision slows convergence by a factor of 2–4 within 10–20 Myr [measured: India slowed from about 15 to about 4 cm/yr around 50–40 Ma].

---

## 5. Motion

### 5.1 Stage 2: rigid plates and a torque balance

The forces of §4 are summed as torques per plate, r × F at each marker and interface marker, and solved for one ω per plate. The system is linear in the ω, a symmetric positive-definite block system of size 3P × 3P [derived]; at most about 90 × 90, solved directly. It's cheap, and it gets the forces right before deformation is added.

### 5.2 Stage 3: the sheet on the sphere

L1's thin viscous sheet (England & McKenzie 1982), moved onto the sphere and the marker mesh:

- **Unknown:** a tangential velocity at every vertex. Interface sides are separate vertices, so velocity jumps across interfaces.
- **Equation:** depth-integrated viscous stress, in balance with GPE (§4.2), basal drag and mantle traction (§4.7), slab pull and resistance (§4.3–4.4), suction (§4.5), and interface coupling (§4.8).
- **Discretization:** linear finite elements on the triangles. The operator is the gradient of a dissipation, so the system is symmetric positive definite (as L1's). Conjugate gradients with Picard iteration for the nonlinearity, as L1.
- **Viscosity:** η = η₀ · s · (1 − D) + η_min, with a power-law strain-rate dependence (n = 3, as L1).
  - s is the strength factor: cratons strong; ocean strengthening with √age.
  - η₀ high enough that undamaged lithosphere moves almost rigidly.

  Plates are then **regions of near-rigid motion** separated by weak interfaces and damaged zones, which is what the pattern definition of §3.5 detects.
- **Strain:** markers move with the vertex velocities, and strain changes A and T at constant V (as L1's markers). Crumple zones thicken by deformation; refinement (§3.3) adds markers there.

### 5.3 The rigid limit, as a test [derived]

With η₀ very large and interfaces weak, the sheet must reproduce the torque balance of §5.1: each plate's best-fit rotation from the sheet should match its ω from the balance. This tests the sheet solver against something simpler and independent.

---

## 6. Exogenous processes

### 6.1 Mantle flow [rule with measured parameters]

- **Form:** v_m = ∇Φ + r̂ × ∇Ψ, degrees 1–4.
- **Spectrum:** most power at degree 2 [measured: two antipodal large low-shear-velocity provinces (LLSVPs), under Africa and the Pacific].
- **Time:** each coefficient follows an Ornstein–Uhlenbeck process, dc = −c dt/τ_m + σ √(2/τ_m) dW:
  - τ_m ≈ 50–200 Myr [measured order: mantle overturn];
  - σ set so the RMS of v_m is 1–3 cm/yr [calibrated].
- **Coupling:** through basal traction (§4.7). Upwelling (the divergence of v_m) adds dynamic topography. An optional mild two-way coupling (slabs feeding downwelling into Φ) is off by default.

### 6.2 Plumes [rule with measured parameters]

| Phase | Rule | Earth range |
|---|---|---|
| Birth | Poisson process, rate λ_p; position drawn with probability rising toward the edges of the degree-2 upwellings | About 1 per 5–10 Myr, keeping 10–20 deep plumes alive [measured order]. Plumes from LLSVP margins: Burke & Torsvik (2004) |
| Head arrival | After about 10 Myr of rise: a large igneous province of 10⁶–10⁷ km³ (log-uniform) over 1 Myr, in a disc 1,000–2,000 km across; a swell of about 1 km; damage added; basal drag reduced locally for 10–20 Myr (thinned, lubricated base) | Deccan, Siberian Traps, Paraná–Etendeka [measured] |
| Tail | Flux 0.01–0.2 km³/yr; log-flux pulsing as an Ornstein–Uhlenbeck process with τ ≈ 10 Myr | Hawaii's flux has varied severalfold along the chain [measured] |
| Drift | Moves with the mantle flow at 0.3–0.5 of v_m | 0–2 cm/yr, briefly about 4 [measured: Tarduno 2003] |
| Death | Lifetime lognormal, median about 80 Myr | 50–130+ Myr [measured] |

A plume may "switch on" abruptly from L0's point of view: the head's arrival is an event.

---

## 7. Damage, cracking and welding

### 7.1 Damage in the plate [rule]

Depth-integrated strength S ≈ 10¹²–10¹³ N/m [measured order]: for ocean it rises with √age; for continent with thermal age, highest in cratons. Effective strength is S(1 − D).

> dD/dt = (1/τ_d) · max(0, σ_II / (S(1 − D)) − 1) · (1 − D) + plume term − D / τ_heal

- σ_II is the second invariant of the deviatoric depth-integrated stress, from the sheet.
- τ_heal is a few hundred Myr [conjecture]. Slow healing keeps old sutures weak, so rifts can reuse them: the inheritance of the Wilson cycle.

### 7.2 Coupling damage at interfaces, and welding [rule]

Each interface marker's coupling damage Dc evolves with the slip across it:

> dDc/dt = (|Δv| / v_c) · (1 − Dc) / τ_slip − [|Δv| < v_c] · Dc / τ_lock

- **Slip damages** the interface toward full weakness (Dc → 1).
- **Quiet heals it:** while relative motion stays below v_c (about 0.5 cm/yr [calibrated]), Dc decays with time constant τ_lock. Physically this stands for the shear zone cooling, recrystallizing and re-cementing so the two sides grip each other; the model does not resolve which.
- **Welding:** when Dc falls below D_weld (about 0.1) along a contiguous stretch of interface, those interface markers merge their sides into ordinary interior markers. Their plate damage D is set to the interface's history (high), so **a suture is born weak** and heals only over τ_heal.

**Consequences [conjecture]:** a continental collision that is still converging (India today, about 4 cm/yr) does not weld, because its slip keeps damaging the interface. A collision whose convergence has been stopped by the coupling forces welds after a few τ_lock. Welding emerges from forces rather than from a timer.

### 7.3 Cracking [rule]

When a marker's tensile normal stress exceeds its effective strength:

1. Grow a path both ways along mesh edges by Dijkstra's shortest-path algorithm, with edge cost = length × max(ε, S(1 − D) − σ_nn) × e^ξ. σ_nn is the tension across the edge; ξ is deterministic noise from a hash of the edge's marker ids and the seed, standard deviation about 0.3.
2. Stop when the path reaches an existing interface at both ends.
3. Accept it if its mean overstress is positive. The path's markers become interface markers (kind ridge) with Dc = 1, and the plate becomes two at the next component pass.

Compression doesn't crack plates. It raises damage, and the sheet thickens the crust.

---

## 8. One step, in order

The time step dt is limited so the fastest vertex moves at most 0.5 h₀ (a Courant condition) [rule]. dt varies with the fastest speed, and this is deterministic.

1. Advance the mantle (§6.1) and the plumes (§6.2).
2. Update each column: buoyancy b, GPE, elevation e (§4.1–4.2).
3. Set interface kinds and polarities (§4.6); break off slabs that are due (§4.6).
4. Solve for velocity: torque balance (Stage 2) or sheet (Stage 3+) (§5).
5. Move markers and interface sides; strain markers (A and T at constant V). Restore the mesh (§3.3).
6. Interface bookkeeping (§3.4): ridge accretion, consumption into slab records, transform re-spacing.
7. Damage and coupling damage (§7.1–7.2); welds; cracks (§7.3).
8. Refine and coarsen markers (§3.3, split and merge rules below).
9. Surface processes, ported unchanged from today's rules in Stage 1: diffusion (cotangent Laplacian), erosion with sediment routing, edifice decay, conversions, collapse, plume volcanism.
10. Plate components and lineage (§3.5).
11. Records and events: fault records per marker as now; `tecRec` splatted for L1; the event log split into **causal** events (rift, weld, plate lost, plume birth, head arrival, plume death, slab break-off, polarity flip) and **observed** ones (island emerges, collision begins).

**Refinement and coarsening.**
- **Split** a marker along its principal stretch direction when A_vor > 2 h₀², or where strain rate or thickness gradient exceeds a threshold. Children share A and V equally, copy the column and fault record, and record `parent`.
- **Merge** two neighbouring markers on one plate when both are small (A < 0.5 h₀²), their thicknesses differ by less than 2 km, and their columns pass the formation-linking rule (README §3.11). The merged marker keeps the older id; the merge is logged for L1.

---

## 9. Invariants and ledgers (tested every step)

| # | Invariant | Tolerance |
|---|---|---|
| I1 | Area: \|Σ A − 4πR²\| / 4πR² < 10⁻³; the largest \|m\| / h₀² bounded | As stated |
| I2 | Continental volume: change in Σ A·T over continental markers = arcs + accretion + plumes + conversions − erosion to ocean − conversions to ocean | 10⁻⁹ relative |
| I3 | Oceanic area: created at ridges − consumed at trenches − converted = change in oceanic area | 10⁻⁹ relative |
| I4 | Mesh: every triangle positively oriented; Delaunay after flips; every triangle on one plate; no edge joins interior markers of two plates | Exact |
| I5 | Identity: ids unique, never reused; lineage graphs acyclic | Exact |
| I6 | Determinism: identical state hash at step 500 across two runs, and after restore and replay | Exact |
| I7 | Forces: at the end of each solve, the net force and torque residuals on each plate are below tolerance | 10⁻⁶ relative |

---

## 10. What L1 needs, and how the new core supplies it

| L1 or UI use | Today | New core |
|---|---|---|
| `snapshot()` / `restore()`, exact replay | Typed-array copy | Same, plus interface markers and lineage graphs; the mesh is rebuilt on restore |
| Parcel identity (`trackParcels`, `P.pid`) | plate × N + body cell | Marker `id` plus the split and merge log |
| Reference frame (`refQ`, `refW`, split ancestry) | Plate `Q`, `W`, `parent` | Plate best-fit ω, accumulated Q, lineage graph |
| Seeding the sheet (`owner`, `bodyIdx`, `wK`, `wCol`, `pool`) | Body cells | Markers directly: L1's markers become subdivisions of L0's |
| `tecRec` | Per fixed cell per step | Splatted to the rendering grid |
| `slab`, `slabT` | Mantle-frame bins | Unchanged; filled from interface slab records |
| `extractFaults` | On columns | On markers, same format |
| `elev`, `cellKm`, `dt`, `t`, `steps`, `events`, `hotspots` | Lattice arrays and lists | Rendering grid, h₀, plume list, split event log |
| `buildL1`, `exportHandoff` | On the lattice | On the rendering grid during migration, on markers afterwards |

Old sessions keep reproducing: `tools/session.js` finds the commit whose sources match, and the lattice core is archived in `history/` as `v6_lattice_sim.js`.

---

## 11. Staged plan

The lattice core stays selectable in the page until Stage 1 passes, so both can run on one seed.

| Stage | Builds | Exit checks |
|---|---|---|
| **1. Markers, mesh, interfaces; today's kinematics** | Markers with explicit T and e; split-node interface markers; convex hull and flips; ridge accretion and consumption at interfaces; refinement; components and lineage; rendering grid. Velocities still a random walk, splits still by timer, welds still by counter, hotspots as now | I1–I6. Over 800 Myr on three seeds: land and continental fraction, peak elevation, plate count, rift and suture counts within today's ranges (README §3.13). Resampling flicker about 0. Cost within 2× of today |
| **2. Local forces, torque balance** | §4 and §5.1. Random walk removed | RMS speed 3–6 cm/yr after calibrating c_geom alone. Plates with more trench boundary move faster, rank correlation > 0.5 [measured: Forsyth & Uyeda 1975]. Collision slows convergence 2–4× within 20 Myr. Net rotation < 0.4°/Myr [measured: estimates 0.1–0.4°/Myr]. Speed distribution has a tail: report how often any plate exceeds 10 and 15 cm/yr, and for how long |
| **3. Sheet deformation** | §5.2 | Rigid limit matches Stage 2 (§5.3). Collision zones thicken; volume ledger closes. Plates detected as patterns agree with interface-bounded components |
| **4. Mantle and plumes** | §6. Fixed hotspots removed | Plume lifetimes, drift, large igneous province frequency within range. Hotspot tracks bend and age along their length |
| **5. Damage, cracking, emergent welding and polarity** | §4.6 hysteresis and break-off; §7. Split timer and merge counter removed | Rift rate 15–60 per Gyr without being set. Rifts start on old sutures more often than random paths through the same plate. Back-arc extension where slabs are old, compression where young (Uyeda & Kanamori). Time from continental contact to weld reported; India-like collisions still converging at 50 Myr are allowed |
| **6. Crust budget without the thermostat** | Arc production against sediment subduction and subduction erosion | Continental fraction 30–45% for 1 Gyr. Earth's continental volume has been roughly steady, with addition and loss each about 3 km³/yr [measured: Scholl & von Huene 2009]. If it fails, the thermostat stays and is labelled a fudge |
| **7. L1 on markers** | L1 seeds from L0's markers and reads the lineage graph | test_l1 passes in what it checks |

### 11.1 Stage 1 as built

**Where the build departs from this spec**, and why:

| Spec | Built | Why |
|---|---|---|
| Every marker carries its area A (§3.1) | Only continental markers carry area. Oceanic markers take their area from the geometry each step: their share of the triangles around them, plus a share of the interface markers beside them | Carried oceanic area drifted: the opening at ridges was stored on interface markers, which come and go. Oceanic crust is made and destroyed by definition, so only continental volume needs a ledger (DEVLOG) |
| Ridges emit a marker from each side once it has opened h₀ (§3.4) | A new oceanic marker is placed wherever a plate's nearest marker is more than 1.3 h₀ from an interface marker that is not at a convergent boundary | No state stored on short-lived interface markers |
| Interface sides carry the frontier's material (§3.2) | Sides carry plates and positions only; frontier values are read from the nearest interior marker of each plate. Interface markers hold no material and no area | Enough for Stage 1's rules. Material on sides comes with the sheet (Stage 3), where velocity jumps across interfaces |
| Polarity by segment-averaged buoyancy (§4.6, Stage 5) | Each connected stretch of convergent boundary takes the majority verdict of the lattice core's rule (continent, then thicker, then younger) | Neighbouring interface markers that disagreed ran into each other |
| Refinement and coarsening (§8) | Not built | Rigid plates keep marker spacing; the marker count holds steady |

**Exit checks** (`tests/test_mesh.js`, `tools/compare_cores.js`; three seeds, n = 96, 800 Myr):

| Check | Result |
|---|---|
| I1 area | Continental carried area against the space it covers: within about 1–1.5% throughout; stable |
| I2 continental volume ledger | Closes to about 10⁻⁹ every step |
| I4 mesh, I5 identity | Valid at every check; no edge joins two plates' interior markers; ids unique |
| I6 determinism | Two runs, and a snapshot restored and replayed, end in the same state |
| Land fraction | 0.29–0.31 at 800 Myr (lattice 0.17–0.23); README range 22–34% |
| Continental fraction | 0.35–0.36 (lattice 0.39–0.40, target 0.40): **below**. The lattice core counts about 3.5 times as many oceanic contacts as the marker core consumes, each depositing arc crust; seafloor turnover is the same on both. The lattice's continents lean on resampling |
| Rifts, sutures | 16–28 and 15–30 per 800 Myr (lattice 22–35 and 20–34); the seeds vary more than the cores |
| Plate count | 7–13 over the runs (lattice 7–13) |
| Mean ocean-floor age | 236–268 Myr at 800 Myr (lattice 124–206): **older** than the README range of 70–220 |
| Displayed-rock flicker | 0.55% of continental markers per step: resampling flicker gone from the material (the lattice's 5% is 90% resampling). The grid view still samples moving markers: 7–9% of cells per step |
| Cost | 100–126 ms per step against 29–34: **about 3.5–4 times**, not within 2× |

### 11.2 Stage 2 as built

Plates are driven by the forces of §4 and a torque balance per plate (§5.1), in `MeshWorld.solveVelocities` (`src/mesh.js`); the constants are in `PH` there, each labelled measured, derived, rule, calibrated or conjecture. The random walk remains selectable (`kin: "walk"`, and the page's Stage 1 option). Checks: `tools/stage2_checks.js`, and the force checks in `tests/test_mesh.js`.

**What the balance contains:**

- **Unknowns:** one angular velocity per plate, in the mantle's frame (no net rotation is removed). Stokes flow, so no inertia: each step's velocities follow from that step's forces.
- **Resisting:** basal drag per marker area, C_b = 3×10¹⁴ Pa s/m, five times larger under continental crust; slab resistance C_s·v_n per trench interface marker, on the subducting plate's velocity into the trench; interface coupling on the relative velocity of the two sides (3×10²⁰ Pa s at ridges and transforms, 10²¹ at trenches, 10²² where continent meets continent [conjecture: orders of magnitude, not fitted]), up to a frictional limit (below).
- **Driving:** slab pull at each trench interface marker, toward the overriding plate; the gradient of GPE as a body force on every triangle; trench suction on the overriding plate, c_suc = 0.2, using the previous step's velocity into the trench.
- **GPE** comes from floated columns (water, crust, a lithospheric mantle whose excess density falls linearly to zero at its base, asthenosphere), relative to a ridge [derived]. Ridge push at 100 Myr comes out 3.88×10¹² N/m, against Turcotte & Schubert's 3.9×10¹². The continental root's excess mass is set so 33 km of crust sits at sea level, as the elevation rule has it; a continent at sea level then has the GPE of 100 Myr old ocean.

**Where the build departs from this spec**, and why:

| Spec | Built | Why |
|---|---|---|
| Each trench interface marker carries its slab record (§4.3) | The record hangs from the overriding plate's interior marker nearest each consumed marker. It pulls at the nearest trench interface marker of its own plate within two rings | Interface markers come and go (§11.1); the overriding plate's markers ride with the trench. A first version pulled only where the record's marker was a trench marker's nearest frontier: 15% of the pull went unapplied at n = 48 and 25% at n = 96, one source of resolution dependence |
| Slab pull counts the thermal deficit and eclogite (§4.3) | Also counts oceanic crust above lEcl = 80 km as buoyant. Continental crust going down forms a buoyant plug at the top of the slab, counted up to lPlug = 150 km long [rule; ultra-high-pressure rocks show continental crust reaching 100–200 km, measured] | Without crust buoyancy, young ocean would not resist. Consumption already moves continental crust's volume into the orogen, so the slab cannot keep all of it: the first version counted all of it, and buoyant records swamped the pull (Forsyth–Uyeda ρ = 0.08) |
| Break-off after collision, Stage 5 (§4.6) | Built now as a rule: once no ocean has gone down behind a record for a step (and it is not behind an active oceanic trench), its oceanic part detaches with e-folding time 10 Myr. The plug of continental crust is exhumed on the same time scale once no more continental crust goes down | Without it, a slab behind a collision kept pulling for ever |
| Polarity by length-weighted mean buoyancy (§4.6) | Each interface marker compares its two frontier markers' b (§4.1) for ocean against ocean; each stretch of boundary takes its markers' majority | The Stage 1 vote, with b in place of age. No hysteresis yet (Stage 5) |
| Coupling capped at a frictional limit, reduced by damage (§4.8) | Linear up to the limit (secant coefficient by Picard iteration, from the previous step's velocities; stops when coefficients change by under 2%). Limits [measured orders]: 3×10¹¹ N/m at ridges and transforms (about 20 MPa over 15 km), 2×10¹² at trenches (about 15 MPa over 100–150 km of megathrust), 6×10¹² in collisions (the force Tibet transmits). No damage | Without the limit, coupling took the largest share of the driving torque (38–49%), so the fitted coefficient barely moved speeds. Damage is Stage 5 |
| Welds by collision counter (Stage 1, until Stage 5) | The counter counts collided area in cells of the n = 96 grid (where its threshold of 180 was set), not markers consumed | Counting markers welded plates after a quarter of the collided area at n = 96 compared with n = 48: about 7 plates against 13, and with them a third of the difference in speed |
| Mantle traction (§4.7, §6.1) | v_m = 0 | Mantle flow is Stage 4. The mantle is a fixed frame that resists |
| c_geom of order 10–100 (§4.4), fitted to RMS speed | 6, provisionally | At n = 48, before welding and slab attachment were made independent of resolution, 6 gave RMS 4.3–4.9 cm/yr. At n = 96 (the chosen resolution) the mean speed does not respond to it: 2.0–3.2 cm/yr for c_geom = 3, 4 and 6 |
| Boundary length per interface marker | Half the distance between its two nearest interface neighbours | Straightens the zigzag of edge midpoints: a great circle comes out 9% long (17% with the first estimate) |

**Exit checks.** Three seeds, n = 96, 100–600 Myr (the first 100 Myr is the start's transient: slabs start at zero length); `tools/stage2_checks.js --n=96 --t=600`:

| Check | Result |
|---|---|
| RMS speed 3–6 cm/yr after calibrating c_geom alone | 2.72, 3.15, 2.57 cm/yr. **Fails, and not for want of c_geom**: 3 and 4 give 2.0–3.2 too. Speeds cycle between about 1.5 and 6 cm/yr over 100–200 Myr, with spells near zero (below 0.2 cm/yr for up to 100 Myr on one seed) |
| Plates with more trench move faster, rank correlation > 0.5 (Forsyth & Uyeda 1975) | Pooled ρ = 0.32, 0.26, 0.47; per-sample medians 0.40–0.43. **Fails.** With 13 plates (n = 48 before the welding fix) it was 0.51–0.63 without the frictional limit and 0.36–0.47 with it |
| Collision slows convergence 2–4× within 20 Myr | **Fails the other way: collisions stop.** Of 105–120 contacts per seed, about 45% slow 2× or more and most of those (about 40% of all) stop or reverse; 3–7% fall in 2–4×. Median slowdown 1.6–2.2× at contact, 2.5× to "stopped" once collision is half the pair's boundary |
| Net rotation < 0.4°/Myr | Medians 0.08–0.09°/Myr; maxima 0.25–0.31. **Passes** |
| Speed tail: how often any plate exceeds 10 and 15 cm/yr | Plates larger than 1% of the sphere faster than 10 cm/yr 11–16% of steps; faster than 15 cm/yr 0.2–2.4%; longest spell above 10 cm/yr 16–22 Myr. Earth's best-documented extreme is India, 15–18 cm/yr for about 15 Myr |
| Plate count (reported) | 5.6–5.9 plates larger than 1% of the sphere on average; Earth has about 10 |
| I7 torque balance | Residual below 10⁻¹⁵ for each solve; the frictional limit is met to within the Picard tolerance (2%) |
| Uniform ocean, no slabs | No plate moves (4×10⁻¹⁶ cm/yr): the GPE discretization exerts no torque from a uniform field |
| Cost | About 210 ms a step at n = 96 headless, against 121 with the random walk: faster plates need 6.5 substeps a step against 2.7 |

**What balances the driving torque** (`budget` in `solveVelocities`: each term's torque along each plate's rotation axis, summed over plates; n = 96, one seed, 150–350 Myr). Driving: slab pull 77%, GPE 12%, suction 11%. Taken up by basal drag 49%, slab resistance 29%, interface coupling 22%. Before the frictional limit, coupling took 38–49% and absorbed any change in c_geom.

**Reading the failures [conjecture, to test next].** Under forces, colliding plates stay pressed together, so Stage 1's weld counter keeps filling, and the world welds down to about 6 large plates. Large plates carry more basal drag per kilometre of trench, which flattens the speed–trench relation, and when their slabs break off nothing starts a new trench: the world nearly stops until a timed split happens to open. Earth avoids this through what Stage 2 does not have: subduction starting where it is not already happening (at passive margins, at transforms, after a polarity flip) and mantle flow that keeps plates moving. The checks that fail are those the weld and split rules (Stage 1 stand-ins, to be replaced at Stage 5) and the missing mantle (Stage 4) control. Two cheap tests of this reading:
1. A weld rule that needs convergence to have stopped, not just collided area.
2. Subduction initiation at old passive margins when GPE contrast exceeds a threshold.

Each should raise the plate count, and with it speed and ρ, if the reading is right. Both were tested, with mantle flow added first: §11.3.

### 11.3 Mantle flow, and two rules tested against the Stage 2 failures

§11.2 ended with a conjecture: the failing checks are controlled by the Stage 1 weld and split rules and by the missing mantle. Mantle flow (Stage 4's §6.1, brought forward) was built first, then the two proposed rules, each behind a switch, and the four combinations compared. `tools/stage2_checks.js --n=96 --t=600 --ph.<name>=<value>` reruns any of them.

**Mantle flow, as built (§6.1):**

- **Form:** v_m = ∇_s Φ + r̂ × ∇_s Ψ, with Φ and Ψ sums of real spherical harmonics of degree 1–4, written as homogeneous harmonic polynomials. On the unit sphere a degree-l harmonic f has surface gradient ∇f − l f r. The degree-1 toroidal fields are left out: they are rigid rotations of the whole mantle, a choice of frame.
- **Scaling:** each of the 45 basis fields is scaled to unit RMS speed. They are mutually orthogonal, so the flow's mean square speed is the sum of its coefficients' squares (tested).
- **Spectrum [rule]:** amplitude by degree 0.5, 1, 0.5, 0.3 for l = 1–4 (the degree-2 peak); toroidal at 0.3 of poloidal [conjecture]. Overall RMS 2 cm/yr, the middle of the spec's 1–3 [not fitted].
- **Time:** each coefficient is an Ornstein–Uhlenbeck process with τ_m = 100 Myr, started from its stationary distribution: the mantle was already convecting, and only slabs start from zero. It draws from its own random stream.
- **Grid and coupling:** tabulated on the fixed grid every 2 Myr. It drives the plates through the basal drag of §4.7: torque R c (r̂ × v_m) per marker. Dynamic topography from its upwelling is not built.

**Rule 1: welding once a collision stops [rule].** A plate pair whose collision counter has crossed Stage 1's threshold (collided area) is latched. It welds once its convergence (length-weighted over their convergent boundary; none counts as zero) has stayed under 1 cm/yr for 10 Myr while the plates still touch. A collision still converging (India) stays a boundary however much crust it has consumed. A first version required the whole shared boundary to slip at under 0.5 cm/yr. With mantle drift and friction-limited coupling no long boundary between rigid plates gets that still: no plate welded in 500 Myr.

**Rule 2: subduction starting at old passive margins [conjecture].**
- **Where:** oceanic markers at least 80 Myr old beside continental crust of the same plate are margin.
- **How often:** per plate, at 3.3×10⁻³ per Myr per 1,000 km of such margin.
- **What happens:** the ocean connected to a chosen margin marker becomes a new plate. Each continental marker beside it gets a seed slab of that ocean, 150 km long, so the new trench pulls from its first step. This stands in for the foundering that starts real subduction.

**Fixes found on the way.** These change the default model too:

| Problem | Fix |
|---|---|
| The area budget broke (3.8%) once trenches formed along continents | An interface marker's surrounding area went to all its interior neighbours, but continental markers carry their own area, so their shares fell out. It now goes to oceanic neighbours only (within about 1% again) |
| Tiny plates ran at 35–66 cm/yr | A slab record pulled at its nearest trench marker even when another plate had made that slab. Each record now remembers the plate that fed it, followed through splits and welds |
| A plate of one marker "spun" at 10⁵–10⁶ cm/yr | Drag on a point cannot resist spin about that point. Each marker now resists spin as a disc of its area, and plates smaller than 10 cells of the n = 96 grid (about 90,000 km², or 10 markers on a coarser grid) join their longest-bordering neighbour, logged as "absorbed" |
| Microplates spun under off-centre pull at up to 80 cm/yr | Slab resistance acted only on motion into the trench. A slab dragged along the trench resists too: edgewise drag of a disc in Stokes flow is 2/3 of broadside drag [derived], applied to the sinking part only (sideways motion does not bend the slab) |

**The four combinations.** n = 96, c_geom = 12, three seeds, 100–600 Myr; ranges are across seeds:

| | Mantle only | + rule 1 | + rule 2 | + both |
|---|---|---|---|---|
| Plates larger than 1% (Earth about 13) | 7–8 | **11–12** | 17–20 | 19–21 |
| Subduction starts / welds / timed rifts per Gyr | 0 / 60–64 / 54–56 | 0 / 16–20 / 12–20 | 110–162 / 20–26 / 0 | 124–160 / 8–10 / 0 |
| RMS speed (cm/yr) | 2.7–3.1 | 3.6–4.0 | 3.2–3.3 | 3.1–3.3 |
| Forsyth–Uyeda ρ | 0.25–0.47 | 0.44–0.56 | 0.58–0.75 | 0.65–0.69 |
| A plate larger than 1% faster than 10 cm/yr | 2–15% of steps | 0–1% | 37–48% | 33–47% |
| Net rotation, median (°/Myr) | 0.07 | 0.07–0.08 | 0.04–0.05 | 0.03–0.04 |

- **Mantle flow alone does not cure the Stage 2 failures.** A large plate integrates a degree-2 flow over a large area and feels little net torque. With mantle flow and the Stage 1 counter, the world churns: welds and timed rifts both run at about 60 per Gyr, and about 8 plates survive.
- **Rule 1 brings the plate count to Earth's** and raises speed and the speed–trench correlation.
- **Rule 2 is right in kind but over-produces.** It keeps the world reorganising and gives the strongest speed–trench correlation. But it adds small oceanic plates faster than welds remove them, and with full slabs and little drag they run near the free sinking speed (about 14 cm/yr at c_geom = 12): the tail fattens and net rotation falls below Earth's estimates. At a quarter of the rate, starts fell only to 62–78 per Gyr (older margins accumulate), with 16–17 plates and the same tail. So the problem is in the rule's form, not just its rate. It stays available as an option (`initiation: true`; the page's experimental core option) and is off by default.

**Default now:** mantle flow and rule 1, rule 2 off, c_geom refitted to 10. Record: `tools/stage2_checks.js --n=96 --t=600`, three seeds, 100–600 Myr:

| Check | Result |
|---|---|
| RMS speed 4–5 cm/yr (fitted) | 4.46, 3.81, 4.37 |
| Forsyth–Uyeda ρ > 0.5 | 0.46, 0.38, 0.51 (per-sample medians 0.42–0.50): **marginal, below on two seeds**. At c_geom = 12 it is 0.44–0.56 |
| Collision slows convergence 2–4× within 20 Myr | Once collision is half the pair's boundary: 2–4× in 30% of cases, 2× or more in 40%, stopped in 10%; median 1.3–1.8×. **Closer, still short of typical** |
| Net rotation < 0.4°/Myr | Medians 0.087–0.091, maxima 0.16–0.17 |
| Speed tail | A plate larger than 1% faster than 10 cm/yr 0–33% of steps, faster than 15 cm/yr at most 1.3%; longest spell above 10 cm/yr 0–64 Myr |
| Plates larger than 1% | 10.6–11.9 (Earth about 13) |
| Welds and timed rifts per Gyr | 14–20 and 14–22 |
| Cost | About 310 ms a step at n = 96 headless, against 150 with the random walk: 5.3 substeps a step against 2.5, and about 55 ms for the force balance |

**What remains [conjecture].** The plates are still split by Stage 1's timer, and welded by a rule. New subduction starts only where the timer happens to split a plate. The next real gains should come from Stage 5: cracks that follow damage instead of the timer, and initiation where forces exceed strength, rather than either rule here.

### 11.4 Stage 5, first part: stress, heterogeneous strength with memory, and cracks

Plates crack where the stress inside them exceeds their strength; the timed splits are gone. Built in `MeshWorld.stressStep`, `failureStep` and `crack` (`src/mesh.js`); constants in `PH`.

**Stress recovery [derived form].**
- **The plate as a sheet.** Rigid plates have no strain, but they carry stress. Each plate is treated as a thin elastic sheet loaded by every force of the balance, with the solved velocities: basal and mantle drag per marker, GPE per triangle, and slab pull, resistance, suction and coupling at its boundary. This is the standard way intraplate stress is modelled (Richardson 1992; Coblentz and others).
- **The loads are in equilibrium by construction.** On a sphere the only rigid motions of a shell are rotations, so zero net torque, which the balance gives each plate, is exactly the condition for its tangential loads to balance. Assembled force by force, each plate's loads sum to 0.05–0.3% of their total magnitude in net torque. That residue comes from moving boundary forces onto frontier markers and from the coupling's Picard tolerance; it is removed exactly before the solve.
- **On a fixed grid.** The solve uses a grid of 8,000 Fibonacci nodes (about 250 km apart), triangulated once. Each node takes the plate of its nearest marker; a grid triangle belongs to a plate when its three nodes do.
  - **Loads:** a marker's load is shared among its plate's nodes around it by barycentric weight.
  - **Element:** plane-stress linear triangles, Poisson's ratio 0.25.
  - **Solver:** conjugate gradients with a Jacobi preconditioner, to 10⁻³; a spring to ground of 10⁻⁵ removes the rotations.
- **Stiffness tracks strength [rule].** Element stiffness is in proportion to effective strength, limited to 1:50. An isolated weak patch is shielded by the stronger plate around it; a weak band crossing the plate still carries the plate's load. With uniform stiffness, the stress resultants (N/m) would not depend on the elastic constants at all.
- **Interpolation:** stress goes back to markers from the nodes of their own plate, projected onto each marker's tangent plane. Within a grid spacing of a plate's edge, a marker takes the nearest fully surrounded node's stress: the edge is not resolved.
- **Measured** (no damage, n = 48, 60 Myr): σ_II median 5×10¹² N/m, 90th percentile 2×10¹³, 99th 4×10¹³. Earth's intraplate stress is 10¹²–10¹³ N/m (10–100 MPa over about 100 km), so the median is Earth's and the tail sits where slab pull is transmitted.

**Strength, heterogeneous, with memory.**

| Part | Rule | Kind |
|---|---|---|
| Ocean | 1.5×10¹³ N/m at 100 Myr, rising with √age (capped at 150 Myr) | measured order |
| Continent | 6×10¹² N/m at 33 km of crust, in proportion to 33/T within 0.5–1.5 (thick crust is weaker); twice in cratons (columns whose base unit is cratonic basement) | measured order |
| Heterogeneity | Lognormal, σ = 0.3, from a smooth random field with features about 500 km across, evaluated in the plate's own frame when a marker is first used, then kept: it travels with the material | conjecture |
| Damage D | Effective strength S·(1 − 0.75 D): fully damaged rock keeps a quarter (faults keep their friction) | conjecture |
| Damage growth | dD/dt = max(0, σ_II/(3 S_eff) − 1)(1 − D)/τ_d − D/τ_heal, τ_d = 10 Myr, τ_heal = 300 Myr. Compression and shear are compared with three times the tensile strength | measured order for the factor; conjecture for the times |
| Past failures | A crack's path takes D = 1 and its sides 0.7; a suture 0.8 beside the vanished boundary; the starting continents' old orogenic belts start at 0.5 | rule |
| Scale | All strengths times 6, fitted at n = 96 to a rift rate of 15–60 per Gyr | calibrated |

**Cracks (§7.3, on the stress grid).**
- **Seeds:** nodes whose every grid triangle belongs to their plate, where the greatest tension exceeds the node's strength (the area-weighted mean over the markers it stands for).
- **Paths:** Dijkstra over the plate's nodes, edge cost = length × max(5% of strength, strength − tension across the edge). The cheapest path to the plate's edge is taken on each side of the seed, along the direction of least tension, and the two are joined where they meet; a forked crack drops the shared stub.
- **Acceptance:** the mean overstress along the path must be positive, and each part must be at least 60 cells of the n = 96 grid (about 0.55 million km²).
- **Split:** each marker follows its node; the damage band is one grid spacing wide.
- **Limit:** at most two cracks per stress solve. The stress is solved every 5 Myr.

**Where the build departs from §7, and why:**

| Spec | Built | Why |
|---|---|---|
| σ from the sheet (Stage 3) | Recovered from the rigid plates' loads on a coarse elastic sheet | Gives stress without deformation; the sheet of Stage 3 will replace it |
| Paths along marker edges, noise per edge | Paths on the stress grid; heterogeneity as a smooth field | Marker-level paths and per-marker noise made cracks depend on resolution: the weakest point and the cheapest route both improve with more markers |
| Damage from σ_II against the strength | Against three times the tensile strength, with residual strength | Against tensile strength, compression near trenches damaged everything; with no residual strength, fully damaged ground cracked under any stress |
| Seeds at markers | Seeds at resolved nodes, off plate edges | Elastic stress under concentrated boundary loads is unbounded on refinement; within a grid spacing of an edge it is not resolved |

**How it got here** (details in DEVLOG):
- **Cracking ran away at first:** 85–300 rifts per Gyr, half the surface damaged, every check failing. The causes were, in turn:
  - damage feeding on itself under uniform stiffness;
  - compression compared with tensile strength;
  - zero residual strength;
  - point-load stress at plate edges;
  - seeds on edge markers blocking interior ones;
  - forked paths rejected;
  - marker-scale paths and noise.
- **Each fix was physical or about resolution, not a retuning.** Afterwards, damage covers 3–6% of the surface and rift rate responds smoothly to strength.

**Record.** n = 96, strength scale 6, three seeds, 100–450 Myr:

| Check | Result |
|---|---|
| Rift rate 15–60 per Gyr (fitted) | 74, 60, 63: the upper edge. At scale 5: 29–91 |
| Rifts start on old weaknesses more often than random (§11 Stage 5 check) | Damage along crack paths before they cracked 0.38–0.53, against 0.06–0.07 for their plates as a whole: **5–8 times** |
| Damage localized | Mean 0.03–0.04; 3–4% of the surface above 0.5 |
| Plates larger than 1% | 8.8–11.8 (Earth about 13) |
| RMS speed | 3.8–4.4 cm/yr |
| Forsyth–Uyeda ρ | 0.37–0.43 (per-sample medians 0.38–0.47): still below 0.5 |
| Collision slowdown | Median 1.2–1.6× once dominant; 2–4× in 12–25% |
| Net rotation | Medians 0.06–0.08°/Myr |
| Cost | About 550–800 ms a step at n = 96 (the stress solve, about 1 s every 5 Myr, and its loads and strength passes) |

**Open [conjecture unless stated]:**
1. **The fitted strengths are 6 times the literature's.** Old ocean comes out at 9×10¹³ N/m, continental crust at 3.6×10¹³ (7×10¹³ in cratons). I read this as the stresses being too high rather than the rocks being too strong: in this force model most of the slab pull net of resistance is carried through the plate, while estimates for Earth put most of a slab's weight into sinking the slab itself (Conrad & Lithgow-Bertelloni 2002; Schellart 2004).
2. **Some resolution dependence remains.** With the same constants, n = 48 rifts less often. One source is identified: the small-plate absorption threshold is set in markers below n = 96 (0.37 million km² at n = 48, 0.09 at n = 96), so more small plates survive at n = 96, each with edges to crack.
3. **Elasticity is not yielding.** Stress above strength is not capped and does not redistribute except through damage-softened stiffness at the next solve. A secant-plastic iteration (as for the coupling's friction limit) would cap it, at several solves per step.
4. **Cost.** The solve takes 400–1,300 iterations of about 1 ms each. Multigrid or the GPU would make it cheap enough to run every step and finer.

Statistics reported at every stage, against Earth:

| Statistic | Earth |
|---|---|
| Plate count and sizes | About 50 when microplates are counted, with a power-law tail of small plates [measured: Bird 2003] |
| Mean ocean-floor age | About 60 Myr; oldest in place about 180 Myr |
| Ocean area against age | Declines roughly linearly [measured: Parsons 1982] |
| Plate speeds | RMS about 4–5 cm/yr; India about 15–18 cm/yr for about 15 Myr is the best-documented extreme |

For comparison, today's core gives a mean ocean-floor age wandering between 70 and 220 Myr (README §3.13): it recycles seafloor too slowly.

---

## 12. Cost estimate [conjecture until measured]

Per step at about 50,000 markers, JavaScript, one core:

| Part | Estimate |
|---|---|
| Move markers, flips, Voronoi areas | 6–16 ms |
| Interface bookkeeping, refinement | 3–6 ms |
| Components and lineage | 2–3 ms |
| Torque balance (Stage 2) | 1–3 ms |
| Sheet solve (Stage 3+), Picard plus conjugate gradients, about 100k unknowns | 50–200 ms |
| Damage, cracking | 2–5 ms |
| Surface processes | 5–10 ms |
| Rendering grid | 3–6 ms |

Without the sheet: about 25–50 ms per step, against today's 22–25. With the sheet, it dominates. Options, in order: warm-start conjugate gradients from last step's velocity (velocities change slowly, so few iterations should be needed [conjecture]); solve the sheet every few steps and move with the last velocity between; Web Workers.

---

## 13. Open questions

1. **Subduction initiation.** New trenches currently appear only where plates already converge (on cracks that start closing). On Earth, old passive margins sometimes start subducting. One route: damage plus the buoyancy contrast at old ocean–continent margins. Not yet specified.
2. **Continental subduction and exhumation.** Ignored: continental crust never goes down an interface.
3. **Ridge migration and asymmetric spreading.** Ridges are symmetric by construction (§3.4).
4. **Parameter count.** The calibrated set is c_geom, c_suc, c_int and S_int per kind, δb, v_c, τ_lock, and σ of the mantle modes. Each has one target statistic named above. If two of them only matter through a ratio (as slab pull and drag did in revision 1), the spec should say so and calibrate the ratio.
