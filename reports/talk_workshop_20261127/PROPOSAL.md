# Workshop talk, Friday 27 Nov 2026 — proposal (slides + script)

15 min talk + 5 min questions. 11 slides, about 1 350 scripted words (≈ 10½ min at 130 words/min). Pauses for the movies and pointing at
figures typically add 20–30 %, which brings it to the ≈ 14 min in the timing table; rehearse once with a timer.
Style taken from the previous deck (14h00_PLT-2551_Loewe.pptx): white slides, orange title bar with a full-sentence claim as the title,
large figures/movies, small citation bottom-right; title slide with the UC Chile + ANID logos and a large simulation image on the right.

Story in one line: **topology fixes how much charge the cavity holds; geometry, activated by the coherent flow, decides where it goes.**

Numbers follow the 22 Oct revision, not the PRL draft where they differ (see "Things to keep out of the talk" at the end).

---

## Timing overview

| # | Slide title (claim) | Visual | Time |
|---|---|---|---|
| 1 | Geometry, not topology, places the charge of active nematic cavities | Title; Fig. 1(d) snapshot (vorticity + tilted wall defects) | 0:20 |
| 2 | Wall anchoring fixes the total charge — but not where it sits | Schematic: winding anchoring k, passive ring of 2k +½ (Fig. 1a passive); prior work (Norton, Klein) | 1:15 |
| 3 | We simulate active nematics in cavities with imposed winding | Fig. 1(b) anchoring + model box (MPCD, ζ, ℓa, k 2–10, 675 shapes) | 0:50 |
| 4 | An active cavity keeps exactly 2k−2 defects at the wall and sends +1 into the bulk | **Movie S1** + Fig. 1(a) density maps + Fig. 1(c) | 1:50 |
| 5 | The anchoring geometry ties each wall defect's tilt to its position | Fig. 1(b) lock schematic; 2k−2 ring vs 2k ring cartoon | 1:45 |
| 6 | A coherent vortex imposes a common tilt, and the 2k ring must eject two | **Movie S2** (passive + imposed vortex) + Fig. 2(c) | 1:45 |
| 7 | The tilt is the flow-aligned angle in the no-slip shear: τ = arccos(1/λ) | Fig. 3(a) and 3(b) | 1:30 |
| 8 | Elongating the cavity breaks the vortex, and the wall layer leaks | **Movie S3** + Fig. 4(a), 4(b) | 1:30 |
| 9 | Untilted wall defects propel themselves out; −½ layers, without jets, do not leak | **Movie S4** + Fig. 4(e), 4(g) | 1:45 |
| 10 | *(work in progress)* Crowded walls leak ten times less: their defects resist the same inward flow | New figure (flow-at-core result, 23 Nov) | 1:00 |
| 11 | Conclusions | 4 bullets + Fig. 1(d) thumbnail; acknowledgements strip (FONDECYT 11250894, NIAMH-MPCD, Tyler) | 0:45 |
| | | **Total** | **14:15** |

Slide 10 is optional: drop it to reach ~13:15 with no rush, or keep it if the audience likes open problems. The result is 4 days old (see the caveats under slide 10).

---

## Slide-by-slide script

### 1 · Title (0:20)
**Slide:** "Geometry, not topology, places the charge of active nematic cavities". Benjamin Loewe (PUC), Tyler N. Shendruk (Edinburgh). Right half: Fig. 1(d), the cavity with its vortex and tilted wall defects.

> Thank you. I'll talk about what happens to topological charge when you confine an active nematic in a cavity whose walls force a winding. This is work with Tyler Shendruk in Edinburgh.

### 2 · The question (1:15)
**Slide:** left, the anchoring: the director turns k times along the wall, so the cavity must hold charge k. Middle, the passive answer: a ring of 2k +½ defects at the wall (Fig. 1a, passive panel). Right, two small boxes: "Norton et al.: turbulent disks screen the boundary charge in a thin layer" and "Klein et al.: up to q = 5, below pair creation, defects braid near the wall". Footer: "Topology says how much, not where".

> If the anchoring on the wall turns the director k times around, topology tells you that the cavity contains a total charge k. It does not tell you where that charge sits.
> For a passive nematic the answer is simple. All of it stays at the wall, as a ring of 2k half-integer defects, and the interior is defect-free.
> Active nematics are different. They create and annihilate pairs all the time. Previous work found that the boundary charge is screened: in turbulent disks it stays in a thin boundary layer, and below the pair-creation threshold the extra defects stay close to the wall.
> We pushed to large windings, up to k equal to ten, in the regime of continual pair creation. The question is simple: where does the charge go?

### 3 · Setup (0:50)
**Slide:** Fig. 1(b) (anchoring with winding k, homeotropic sites marked). Box: "NIAMH-MPCD: 2D extensile active nematic; activity ζ = 0.07, active length ℓa = 22 cells; flow-aligning, λ = 2; no-slip walls, strong anchoring θ = kφ. Circles R = 0.2–2.9 ℓa, k = 2–10; ellipses and Cassini ovals, aspect 1.1–2.5 (675 cavities); continuum Q-tensor model with backflow as a check."

> We use multi-particle collision dynamics, a particle-based solver for nematohydrodynamics, with extensile active dipoles. The relevant scale is the active length, about twenty cells here. The walls are no-slip and impose strong anchoring that winds k times. We vary the winding, the size from a fifth to three active lengths, and the shape, from circles to elongated ellipses and ovals. Everything I show is checked against a continuum Q-tensor model with backflow.

### 4 · Observation (1:50)
**Slide:** left, **Movie S1** (k = 6 circle: 10 wall +½, bulk +1, one vortex). Right top, Fig. 1(a): passive 12 at the wall vs active. Right bottom, Fig. 1(c): bulk charge = 1 and wall count = 2k−2 over the shaded range. Callout: "Exactly 2k−2 at the wall, exactly +1 in the bulk: k = 2–10, R = 0.8–2.2 ℓa (wall occupancy 0.997, bulk charge 1.01)".

> Here is what an active cavity does, for k equal to six. A passive cavity would have twelve defects at the wall. The active one keeps ten, sharply localised, and these sites almost never empty. In the bulk, pairs are created and annihilated all the time, but the net charge in the bulk is exactly plus one.
> You can also see that a single vortex fills the cavity, and the wall defects are all tilted in the same sense, following the vortex.
> This is not specific to k equal to six. Over windings from two to ten and sizes from about one to two active lengths, the wall holds two k minus two defects and the bulk plus one. Topology is satisfied, but the boundary does not screen the charge: exactly one unit always moves into the bulk. Too small and the vortex cannot form, too big and the cavity becomes turbulent.
> So why two k minus two, and why plus one?

### 5 · Geometry (1:45)
**Slide:** Fig. 1(b) enlarged: 2k−2 homeotropic sites, wall defects with their tails. Equation: τⱼ = (2k−2)φⱼ + c ("lock", obeyed with mean resultant length 0.99). Two cartoons: a ring of 2k−2 sliding rigidly (common tilt, no elastic cost) vs a ring of 2k whose tilts wind twice around (cannot tilt as a whole).

> The answer is in the geometry of the anchoring. Relative to the local wall normal, the anchored director turns by (k minus one) times the angle. So it is exactly homeotropic at two k minus two points around the wall.
> Now take a wall defect. Its orientation, the direction of its tail, is fixed by where it sits: the tilt from the normal is (two k minus two) times its position, plus a constant. The simulations follow this lock almost perfectly. So position and tilt are one degree of freedom.
> That has a consequence. A ring of two k minus two defects can sit on the homeotropic sites, untilted, and it can tilt as a whole simply by sliding all its defects together along the wall. A ring of two k cannot: with that spacing, the lock makes the tilts wind twice around the ring. There is no way to give all of them the same tilt.
> In a passive cavity this does not matter, because nothing asks for a common tilt, and the two k ring is perfectly stable.

### 6 · The flow does it (1:45)
**Slide:** left, **Movie S2** (passive MPCD, imposed vortex, the ring of 12 ejects two). Right, Fig. 2(c): wall count 12 → 10 when the tilt coherence reaches ≈ 1; no ejection for weak flow or a tumbling nematic. Small text: "Forcing a common tilt on the 2k ring costs ≈ 2× the energy of moving two +½ into the bulk; active continuum model reproduces the state".

> What asks for a common tilt is the flow. The coherent vortex shears the wall layer in one direction everywhere, and flow alignment tilts every wall defect the same way. For the ring of two k, forcing a common tilt costs roughly twice the elastic energy needed to move two defects into the bulk. So the ring ejects two. What is left is the largest ring that can tilt as a whole, and the bulk receives plus one.
> To show that it is the flow and not activity itself, we take a passive cavity and impose a vortex by hand. Once the flow is strong enough to give all wall defects the same tilt, the ring of twelve ejects exactly two, here in the movie. With weak flow, or with a tumbling nematic that has no preferred angle, nothing happens.
> And an active continuum model, started from the passive ring, organises itself into exactly the same state.

### 7 · The tilt angle (1:30)
**Slide:** left, Fig. 3(a) (flat wall: sliding the screening defect rotates the far field; in shear the far field sits at the Leslie angle; τ = 2θ_L). Right, Fig. 3(b): measured tilt vs λ on arccos(1/λ). Equation large: τ = arccos(1/λ).

> Which tilt does the flow select? Think of a flat no-slip wall with the winding anchoring along it. Sliding a wall defect along the wall rotates the director far from the wall, at no energy cost. In shear, the far field wants to sit at the flow-aligning, Leslie, angle. So the layer slides until the far field is aligned, and the defect tilt is twice the Leslie angle, which is arccos of one over lambda.
> The prediction does not contain the shear rate or the elastic constants. Only the ratio of strain rate to vorticity enters, and the no-slip wall fixes it.
> The active simulations follow it within three degrees for lambda from two to six. When we remove the shear next to the wall in the continuum model, the tilt is lost, and confining all shear to the wall band leaves it unchanged.

### 8 · Beyond circles (1:30)
**Slide:** left, **Movie S3** (aspect-2 ellipse: two vortices, the layer leaks, bulk +2 to +3). Right, Fig. 4(a) (tilt coherence collapses at aspect 1.25–1.37, every size, ellipses and ovals) and Fig. 4(b) (vacant fraction; crowded walls leak less).

> What if the cavity is not a circle? When we elongate it, the single vortex breaks, at an aspect ratio of about 1.3, for ellipses and ovals and at every size. Without the vortex, there is no common tilt, and the wall layer starts to lose defects, one at a time. The bulk now holds two, three or more.
> The released charge is not random. In elongated cavities it follows a wall law: the bulk charge per unit wall length depends only on the winding density, not on the shape or size, and a straight periodic channel is its flat limit.
> Stronger flow alignment does not save the layer, which tells us the protection really came from the coherent vortex.

### 9 · Why the layer leaks (1:45)
**Slide:** left, **Movie S4** (one escape, zoomed: the defect untilts, its head turns inward, it leaves). Right top, Fig. 4(e): escapers untilt over ≈ 40 t.u. before leaving; residents do not. Right bottom, Fig. 4(g): wall defects lost vs aspect, +k (solid) vs −k (dotted). Callout: "−k: wall −½, no jets, escape rate > 20× lower".

> Why do the defects leave? A plus one-half defect is self-propelled along its head. In the circle the tilt turns this jet along the wall, and the no-slip wall absorbs most of the inward part. Following single escapes, we see that a defect first untilts over about forty time units. Its head turns into the bulk, and then it runs away, against an elastic barrier that has not changed.
> The cleanest test is a mirror control. If we reverse the winding, the same geometry puts minus one-half defects at the wall, two more instead of two fewer. Minus one-half defects are not self-propelled. In exactly the same elongated cavities, the vortex breaks just as before, but the minus one-half layer does not leak: its escape rate is more than twenty times lower.
> So losing the vortex is not enough. The leak needs the wall defects' own jets.

### 10 · Work in progress: crowded walls (1:00) — optional
**Slide:** left, Fig. 4(b) highlight: vacant fraction drops with wall density. Right, new two-panel figure (to make): (i) inward flow at wall cores, sparse vs crowded (same or larger when crowded); (ii) response of the core to that flow at the same depth (slope 0.61 vs 0.09) next to the escape fraction (0.050 vs 0.0045). Tag: "preliminary".

> One open question. Densely packed walls leak much less: in aspect-two ellipses, about ten times less per defect. We first thought crowding weakens the inward push. It doesn't. Measuring the flow at each wall defect, crowded defects get the same inward flow as sparse ones, or more, but they barely follow it: at the same distance from the wall, their response is about seven times weaker. So crowding strengthens the hold on each defect rather than weakening the drive. Interestingly, the continuum model does not reproduce this regime: there, crowded layers either freeze or leak freely.

*Caveats if asked:* 5 seeds per density, one shape and activity; MPCD flow noise is controlled (equal noise fraction at the compared depth); "hold" = motion not carried by the measured flow, elastic in the natural reading, but not yet separated from self-propulsion.

### 11 · Conclusions (0:45)
**Slide:**
- Topology fixes how much charge; geometry decides where: the wall keeps the largest ring that can tilt as a whole, 2|k−1|.
- The coherent vortex activates the geometry by imposing a common tilt, τ = arccos(1/λ).
- Without the vortex (elongated cavities) the wall defects' jets turn inward and the layer leaks, to a wall law; −½ layers do not.
- Testable: microtubule–kinesin or cell layers in wells with patterned anchoring (2k−2 tilted defects in disks, leakage in elongated wells, none for negative winding).
Strip: acknowledgements (FONDECYT 11250894; NIAMH-MPCD, Shendruk lab).

> To summarise. Topology fixes how much charge an active cavity holds, but geometry decides where it goes. The wall keeps the largest ring that can tilt as a whole, two k minus two, and the coherent vortex is what makes that constraint bite, at the flow-aligned tilt. Without the vortex, the wall defects' jets turn inward and the layer leaks, while a layer of minus one-half defects does not. All of this should be testable in active nematics in wells with patterned anchoring. Thank you.

---

## Likely questions (5 min) — short answers

1. **Why exactly +1, not more?** Ejecting two makes the ring commensurate (2k−2 can tilt rigidly); ejecting more brings no tilt advantage and costs elastic energy. Above ≈ 2.2 ℓa the vortex is no longer coherent and the count drifts.
2. **Is this MPCD-specific?** The continuum Q-tensor model with backflow reproduces the ejection, the 2k−2 state, the tilt and the shape transition. It does not reproduce the details of leakage in crowded walls.
3. **What sets the tilt in crowded walls?** It falls below arccos(1/λ) in crowded walls (≈ 42° at a defect spacing of 7 cells, λ = 2), in active and passive layers alike. The fall is ordered by the core depth relative to the alignment length; the torque responsible is not identified yet.
4. **Pinning / Frank anisotropy?** With one Frank constant (K ≈ 44 kBT), the flow fixes the tilt and the lock fixes the positions. There is no static pinning; the splay–bend split in the fluctuations comes from the flow coupling of the algorithm.
5. **Experiments?** Microtubule nematics in wells: their wall defects are usually −½ and +½ are ejected, but that is without imposed winding. Patterned anchoring in wells would test 2k−2 vs 2k and the mirror.
6. **Activity dependence?** We fix ζ = 0.07 (between onset and the regime where MPCD compressibility matters) and vary the cavity size in units of ℓa instead.

---

## Things to keep out of the talk (withdrawn on 22 Oct, still in the PRL draft)
- "Pinned by unequal splay and bend, barrier 80–110 kBT" → say: flow fixes the tilt, the lock fixes the positions.
- "Elasticity alone tilts the layer by ≈ 35°" → do not mention.
- K₁ ≈ 42, K₃ ≈ 23 kBT → if asked, one constant K ≈ 44 kBT in the static response.
- Absolute ejection energies (375–475 vs 180–200 kBT) were computed with K₁ ≠ K₃ → the talk says only "roughly twice" (the ratio survives the one-constant revision).
- The mirror ratio "≈ 30×" → the 10-seed value is 27× at aspect ≥ 1.4 (22× with the looser criterion): "more than twenty times".
- The "topological current" channel movie (X_…mp4) was withdrawn as a tracking artefact; do not use.

## Material
- Figures: reports/paper/prl/figs/fig1–4 (PNG/PDF); panels would be re-cut per slide with larger fonts.
- Movies: reports/paper/si_movies/S1–S4 (7–15 s each, loop them).
- To make: slide 2 schematic (winding anchoring + passive ring), slide 5 ring cartoon (2k−2 vs 2k), slide 10 figure (from coreflow_rows.npz).
