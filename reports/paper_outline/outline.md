# Paper outline — overcharged active cavities (draft, 2026-10-10; λ scan plan and venue added 2026-10-12)

Working title options
- *Flow-induced ejection sets the charge of active nematic cavities*
- *Why an active nematic cavity holds one charge in its bulk*
- *Topological overcharging of active nematics by confined flow*

## One-paragraph pitch
A cavity whose wall anchoring winds k times must contain charge k. A passive nematic puts all 2k +½ defects at the wall. An active one keeps 2k−2 pinned at the wall and moves +1 into the bulk. We show why: the anchoring has a (2k−2)-fold symmetry, the wall layer has a zero mode that ties each defect's tilt to its position, and the active vortex imposes on the layer a common, flow-aligned tilt (arccos(1/λ)) that a ring of 2k cannot adopt. Two defects are ejected. The mechanism works without activity (an imposed vortex does it) and appears in both MPCD and a continuum active nematic. Elongating the cavity beyond aspect ≈ 1.3 breaks the single coherent vortex; the layer then leaks, and it leaks only because the wall +½ are self-propelled (−½ layers do not leak).

## Structure (main text, ~4 figures + methods; or a longer PRE/Soft Matter version with 7–8 figures)

| # | Section | Key claims | Evidence status | Figures available | Production needed |
|-|--------|-----------------|---------|---------|-------|
| 1 | Setup and observation | Passive: 2k at the wall. Active circles: 2k−2 pinned + bulk +1, single vortex, tilt ≈ 62° | Established (circles: R/ℓa sweep, k = 2…10, 3 seeds; 10 seeds running) | maps_circles, phase_circles, flow_tilt | circle seeds s04–s10 (server, running) |
| 2 | Anchoring symmetry, tilt–position lock as a zero mode | (2k−2)-fold homeotropic sites; tilt and registry one variable; splay ≠ bend pins the commensurate ring; passive registry set by the initial condition | Established (lock 0.99; IC tests; modulated winding; flat-wall derivation) | fig_lock, anisotropy_pinning | none |
| 3 | Flow-induced ejection | A common tilt costs a 2k ring 200–500 kBT vs 155–200 to eject two; passive cavity + imposed vortex ejects (needs a threshold and flow alignment); continuum Q-tensor (imposed flow; active with backflow) reproduces it | Established | tilted_structures, fig_ejection, shear_tilt (b) | optional: passive swirl seeds (3 → 5) |
| 4 | The tilt angle | Flow alignment in the no-slip shear: τ = arccos(1/λ) at leading order (zero mode + Leslie angle + factor 2 for a +½); slip test; λ dependence; active continuum 61.7–62.1° vs MPCD ≈ 62° | Leading order established; the last 2–5° open (order-dependence compensation, splay ≠ bend; λ = 2 channel runs pending) | shear_tilt, fig_flatwall | **λ scan of the active MPCD circle (see 'λ scan plan')**; channel runs done |
| 5 | Hydrodynamics / effective boundary (could go to SI) | Bulk active stress drives the flow, the layer brakes it; layer = effective planar winding-1 boundary; onset R/ℓa ≈ 0.4–0.6 | Supported | dirstokes, effective_bc | small high-k cavities for the onset k-ordering (server data, to reduce) |
| 6 | Shapes | One transition at aspect 1.25–1.37, scale-free, ellipses = ovals; layer loss escape-driven, not a cascade; restoring force lost only with activity AND incoherent flow; −½ layers hold → self-propulsion; continuum model reproduces the transition | Transition established (3 seeds); mechanism supported (mirror: 1 size, 2 seeds; continuum knockout running) | phase_shapes, fig_shapes_mech, fig_restoring, fig_cont_ellipse | mirror shape sweep (chained), seeds s04–s10 for both windings (planned) |
| 6b | Wall law for the overcharging (added 2026-10-13) | Bulk charge per unit wall length saturates at aspect ≈ 1.5 and is then set by the winding density q = k/P alone: Q_b ≈ c(q)·P (∝ √area at fixed q), c ≈ 0.012 at low q, released fraction 0.61 → 0.11 from q 0.02 to 0.07 (crowding); ellipses = ovals = periodic channel (channel ≈ 8 % above the cavity plateau, the limit of the law); circles stay at +1 (topological) | Established (675 shape-sweep runs + 6 channels) | fig_perlength | channel at more q values; long channels (running) |

## Possible follow-up paper: active channels with winding anchoring (added 2026-10-13)
- Flat limit of the arclength anchoring (eMode 4, branch channelwind): periodic channel with opposite windings on the two walls, charge 2qL per period. Walls must sit ≥ 1 cell inside the domain.
- Found so far: same wall law as elongated cavities (≈ 40 % of the charge in the bulk at q 0.031, any width 20–60); vortex row with wavelength ≈ 1–1.5 W; homeotropic control = coherent unidirectional flow (known regime).
- Questions: vortex lattice in long channels; charge per vortex; does winding change the dancing-disclination regime (Shendruk et al., arXiv:1703.01531)? homeotropic width scan to locate laminar/dancing/turbulent in MPCD; wavy walls (eMode 4 uses the local normal).

## λ scan plan (added 2026-10-12)
Purpose: turn the tilt law into a tested prediction across λ (needed for PRL and above), and run the causal test of the tilt protection in elongated cavities. λ = `tumble` (species key, default 2); MPCD flow-aligning for λ > 1, tumbling below.

**A. Tilt law in active circles (paper section 4).**
- Cavity: circle, area 2500 (R/ℓa ≈ 1.5, as fineflow), k = 6, natural activity 0.07.
- λ ∈ {0.8 (tumbling control), 1.2, 1.5, 2, 3, 4, 6}, 3 seeds each (21 runs); subset k = 4 and 8 at λ ∈ {1.5, 3, 6}, 2 seeds (12 runs).
- 3000 t.u., defects every 2 t.u., flow every 20 t.u.; about 33 runs, roughly 3 h on the office machine.
- Measure: wall-defect tilt τ (mean, distribution, sign vs vortex sense); wall count (2k−2?), bulk charge (+1?); vortex coherence; v_rms and the active length (λ also changes the turbulence).
- Compare with arccos(1/λ) and with the particle-averaged λ_eff = λ(1+S4)/(2S2) (measured S2, S4 per run). Law: τ = arccos(1/λ) → 34° (λ 1.2), 48° (1.5), 60° (2), 71° (3), 76° (4), 80° (6). Caution: the PASSIVE imposed-swirl runs (runs/swirl_tumble) gave 57–65° at λ 1.2 (law fails there), 57–68° at 2, 75–83° at 4, and a residual 28–36° tilt even when tumbling (0.5). So the low-λ end is a genuine test, and agreement at 2 and 4 alone is not enough. Tumbling: no locked tilt expected; does the overcharged state survive?
- Figure: τ vs λ (MPCD points, both predictions, active continuum points at 2–3 λ values), inset tilt distributions.

**B. Causal test of the tilt protection in elongated cavities (paper section 6).**
- Running (2026-10-12, runs/lambda): aspect 2.0, λ ∈ {1.2, 4}, 8 seeds; circle controls 2 seeds each.
- Extend to λ ∈ {3, 6} if the first pair is clear.
- Prediction (inward jet = 0.40·cos τ·g(depth)): λ 4 → a 2–3x weaker resident push and a much lower depleting-escape rate. λ 1.2 ≈ λ 2. Evaluate by rate vs measured cos τ (and flow rms), not vs λ alone.
- Figure: depleting-escape rate vs ⟨cos τ⟩ of residents across λ and aspect (one master curve if the mechanism holds).

**C. Continuum cross-check (optional).** Active Q-tensor model (shear_tilt/active_bf.py) at λ ∈ {1.2, 3, 4}: tilt in the disk, and ellipse transition shift with λ.

Decision point: if A follows the law across λ and B collapses on cos τ, aim for PRL (circle mechanism and tilt law) with the escape physics as a companion PRE/Soft Matter paper, or everything as one PRX/Nat. Commun. paper.

## Venue (discussion 2026-10-12, to settle with Tyler)
- **PRL:** circle overcharging mechanism + tilt law (with the λ scan) + continuum check + onset criterion. Realistic.
- **PRX / Nat. Commun.:** one paper unifying overcharging, the tilt law and the escape physics (jet blocking by the wall, lobe erasure, untilting runaway, landing on Q = 0 lines — links to Head et al., Nat. Phys. 2024, and the 'topological kicks' paper). Needs the causal λ test, a trigger story, and production statistics.
- **Nat. Phys.:** only with an experimental realisation or a concrete experimental proposal (patterned anchoring of winding k > 1 in wells).
- Novelty against Klein et al. 2025 has to be the mechanism and the law, not the observation.

## Supplementary material
- Methods: MPCD (collision 21, activity, eMode 2 anchoring, imposed swirl, modulated winding), defect tracking and the head/tail convention, Frank constants and γ1, the Q-tensor models.
- Fluctuations and pinning (Hessian theory vs data), mobility (wall 0.02, bulk 0.0054).
- Half-integer windings; negative layer energetics; modulated winding.
- Point-defect models and why they fail (Stokeslet Langevin); the director-field Stokes fits.
- Tip weakness, interval statistics, coherence conditioning, escape/refill kinetics, bulk carriers.

## Claims to keep modest
- The last few degrees of the tilt (compensation picture not yet nailed down).
- Why MPCD's alignment departs from the particle-averaged law at λ = 2 (method question).
- The self-propulsion reading of the restoring-force loss until the continuum knockout and the mirror sweep are in.
- Onset k-ordering (prediction, untested).

## Risks / referee questions to pre-empt
- Small N for shapes → production seeds (planned, ~9–10 server days).
- MPCD-specific artefacts → the continuum model reproduces ejection, tilt and the shape transition.
- Grid pinning → registry follows the anchoring under modulated winding; passive registry set by the initial condition.
- Literature: Klein et al. 2025 (overcharged circles as a setup); our new parts: the active turbulent regime with exact 2k−2 + bulk +1, the mechanism, the tilt, the shapes.

## Order of work (proposal)
0. (2026-10-12) λ scan: B is running tonight; A next on the office machine (about 3 h), C optional. See 'λ scan plan'.
1. Finish the tilt (channel runs) and the self-propulsion test (continuum knockout; MPCD version with localized activity — needs an anchmod change).
2. Mirror shape sweep (chained), then seeds s04–s10 for both windings.
3. Draft sections 1–4 (circle story) while the shape production runs; add 6 at the end.
