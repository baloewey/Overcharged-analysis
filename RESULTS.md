# Results map: overcharged active nematic cavities

**Question.** An active nematic (NIAMH-MPCD, ζ = 0.07) sits in a cavity whose homeotropic wall anchoring is forced to wind k times (`customAnch`). Where does the excess charge k go? In circles it splits into a **+1 bulk** and a **k − 1 wall layer**. Is that split set by topology (k alone) or by the cavity's geometry?

**Short answer so far.** It is geometric. Circles keep the split over a wide window of size and density. Elongated cavities lose it at the same size and density, and that holds for ovals and ellipses, convex or not. The wall layer is a pinned lattice that leaks defects through empty sites. Crowding the wall (high k) pushes the structure back towards the topological picture.

**Read first:** `results/8_reports/master_document.pdf` (or `reports/master/master.pdf`), the full write-up as a paper skeleton.

**How to use this file.** Each section answers one question and points to the main figures in `results/` (copies, built by `tools/build_results.py`; nothing in `runs/` was moved). Each folder has a `FILES.md` that says what every file is and where it was copied from. Dates and the full reasoning are in the dated log of `RESEARCH.md`; open code issues are in `TODO.md`.

> **Being updated:** circles are going from 3 to 10 seeds per point (server, started 5 Oct; +k ≈ 9–10 Oct, −k ≈ 13–14 Oct). Sections 3, 4 and 7 will be refreshed then; rerun `python3 tools/build_results.py` after re-plotting.

Symbols: ℓ_a = active length (22.2 cells at ζ = 0.07); ρ = 2k/P = linear density of imposed defects; R/ℓ_a = cavity size; Q(r) or Q(c) = charge enclosed inside radius r (or inside a level c of an oval); knee = the +1 plateau of Q before the jump at the wall.

---

## 1. Is the method right? → `results/1_validation/`

**Answer: yes, within ρ ≤ 0.15.**
- Custom anchoring imposes the winding correctly: passive cavities reach alignment 0.96–0.99 with exactly 2|k| defects, all at the wall, bulk 0. The negative winding is an exact mirror.
- **Resolution limit ρ ≈ 0.15** (imposed defects ≥ 6–7 cells apart). Above it the anchoring fails in passive and active runs alike, so it is a limit of the method. Everything below uses ρ ≤ 0.15 (`resolution_limit.txt`).
- New wall shapes validated: arclength anchoring on ellipses (eMode 2), Cassini ovals (a = 0 reproduces the circle; passive ovals keep all defects at the wall even with a 13-cell neck): `cassini_walls_passive_tests.png`, `cassini_passive_pinched_oval.png`.
- Defect tracker rewritten; with walls it gives the same defect lists as the old one.
- Seed bug fixed (office 29 Sep, server 5 Oct). Runs before the fix are statistically independent but can't be reproduced exactly; runs since are byte-reproducible.

**Open:** checkpoint restarts crash (likely a field mismatch, found 5 Oct; `TODO.md`). Not used by our workflow.

## 2. The length unit → `results/2_active_length/`

**Answer: ℓ_a(0.07) = 22.2 ± 0.3 cells**, with no finite-size effect between 100² and 200² boxes (`finite_size_L200_vs_L100.txt`, `active_length_vs_activity.png`).
- Used only to normalise R and ρ. The ℓ_a(ζ) scaling is covered by published work and stays out of the paper (decision 29 Sep).
- The fluid is quite compressible at every activity used (rms Mach 0.66–1.0, `compressibility.txt`).

## 3. Circles: where is there a +1 bulk and a wall layer? → `results/3_circles/`

**Answer: from R ≈ 0.75 to 2 ℓ_a, at every density up to the resolution limit.** In that window the bulk carries +1 and the wall holds 2k − 2 pinned +½ (charge k − 1).
- **Main figure:** `layer_diagram_circles.png` (1112 runs, 3 seeds). The layer degrades above ≈ 2.5 ℓ_a and in very small cavities (< 0.5–0.6 ℓ_a).
- **Fixed size R = 1.43 ℓ_a, k = 2–10, 3 seeds:** bulk = +1.00 ± 0.02 at every k (`fixed_radius_scan/rscan_seedavg.txt`); 2k − 2 wall sites in all 27 runs.
- **Passive cavities:** bulk 0, all 2k defects at the wall. The +1 bulk is an active effect.
- **The effect doesn't fade at low density** (H2 not supported): +1 holds down to ρ = 0.02 at 1.43 ℓ_a.
- **How to measure the bulk charge:** six criteria agree for R ≈ 0.8–2 ℓ_a. The knee of Q(r) is the most robust; a fixed-depth cut fails at small and large R. We use the knee, with the density minimum and defect mobility as cross-checks (`bulk_charge_criteria/bulk_vs_R_by_criterion.png`).
- Wall defects sit on fixed angular sites; bulk defects move. That is the basis of the site rule (`defect_geography/`).

**Confidence:** high (3 seeds everywhere, 10 coming).
**Open:** large circles at high k (k = 15–20, R 2.0–2.6 ℓ_a) show a tilted plateau at Q ≈ 1.3–2.1: a genuine "less-charged layer", or the large-R degradation? Check with 10 seeds.

## 4. Negative winding (mirror runs) → `results/4_negative_winding/`

**Answer: the bulk stays positive.** With winding −|k| the wall carries −½ defects, but the bulk is not −1. It grows with |k|: 0.09 (|k| = 2), 0.47 (6), 0.83 (10), 1.0 from |k| ≈ 15 (`mirror_bulk_vs_k.png`, 570 runs).
- +1 is the attractor for both signs, but activity breaks the mirror symmetry: the wall has to absorb two extra −½ to free a +1 bulk of self-propelled +½.
- Why it falls short: the −½ wall lattice already has the period of a +1 bulk (2|k| + 2 sites), but 10–30% of its sites are empty at low density. Counting wall defects predicts the bulk charge (correlation 0.96; section 7).

**Confidence:** good (3 seeds; 10 coming).

## 5. Does activity enter only through ℓ_a? → `results/5_activity/`

**Answer: no, but ζ = 0.07 is a clean working point** (decided 1 Oct). At R = 1.43 ℓ_a, k = 2–10, 3 seeds (`all_activities_summary.txt`):

| ζ | what happens |
|---|---|
| 0.05 | turbulence onset: wall sites often empty, bulk 1.2–1.8, vortex fails at k = 9–10 |
| 0.06 | bulk slightly above 1 (1.05–1.15) |
| **0.07** | **bulk 1.01–1.09, flattest plateau** |
| 0.09 | bulk +1, but the vortex empties the centre (60–70% of central cells < 10 particles) |
| 0.11 | the layer–bulk gap closes for k ≥ 7 (layer width is fixed in cells, bulk density grows ∝ 1/ℓ_a²) |

- Even at ζ = 0.07 the centre is ~30% depleted by the vortex (`density_check.txt`). State this in the methods.

## 6. Shape: do elongated cavities keep the layer? (main result, H1) → `results/6_shape/`

**Answer: no. Elongation destroys the layer at the same size and density, so the split is geometric.**
- **Main figures:**
  - `accumulated_charge_circles_vs_ovals.png`: circles have a +1 plateau and a jump at the wall; ovals rise smoothly with no plateau.
  - `layer_area_fstar_vs_aspect.png`: the fraction of the cavity outside the Q = +1 level (f*) falls on **one curve against the aspect ratio** for all oval sizes, and **ellipses follow it too**. Circles stay at 0.20–0.32; above aspect ≈ 1.5 the "layer" would be bigger than the bulk (f* > ½).
  - `layer_depth_dstar.png`: size-matched control at half-length 1.53 ℓ_a: the depth of the +1 level doubles with elongation (0.23 → 0.50 ℓ_a vs 0.20 in circles).
- **Sharp crossover** in the small ovals (3 seeds): circle-like at aspect 1.17, partial at 1.29, lost from ≈ 1.4, while still convex (`small_ovals_1.53la/`).
- **Not a size effect:** the ovals that lose the layer sit inside the circles' layer region in every size measure (equivalent radius, half-length, inscribed radius, 2A/P) at matched density (`matching_ovals_to_circles.txt`).
- **Not the neck:** a single convex oval with the shape of one lobe loses the layer too. Pinched ovals show no restored layer per lobe (`tables/pinch_lobe_charge.txt`).
- **Crowded walls push back:** in the k scan the extra bulk charge stays below ~+1 and shrinks at high k; in the pinched oval k = 12–20 partly restores the layer, with a total bulk ≈ +1 (`k_scan/`, `high_k_pinched_oval.png`).
- Ovals carry two vortices, one per lobe, already at b/a = 1.6 (`lobes_and_vortices/`, 1 seed, tentative).

**Confidence:** good for the small ovals and the size series (3 seeds); the largest ovals and ellipses have 2–3 seeds.
**Next:** the oval sweep (`8_reports/oval_sweep_draft.md`, seed count to decide).

## 7. Mechanism: flow and wall-layer dynamics → `results/7_flow_and_wall_dynamics/`

**Answer: the wall layer is a pinned lattice. It never slides; it leaks through empty sites.**
- **Circulation** in circles collapses on R/ℓ_a and breaks down above ≈ 2.5 ℓ_a, the size where the layer is lost (`flow_mean_circulation.png`).
- **Wall +½ lean downstream** by 60–65° (30–40° at ρ ≈ 0.15); the tilt vanishes where the layer is lost (`wall_defect_tilt.png`).
- **Wall lattice** (circles, `wall_fourier_circles_static.png`, `_dynamic.png`): period 2k − 2 (+k) or 2|k| + 2 (−k), pinned, no sliding anywhere. +k: fully occupied from 0.75 to 2.3 ℓ_a; −k: 70–90% occupied (section 4).
- **Ovals vs circles** (`wall_dynamics_ovals_vs_circles.png`): wall +½ creep **against** the wall flow at a rate set by the **empty sites** (rank correlation −0.82 over all runs), not by how well the lattice is pinned. No drift in full circles; up to ~80% of the wall flow speed in elongated ovals; it stops when high k fills the wall.
- **Wall defects escape 10–20× faster** in elongated ovals than in circles of the same size, and the layer empties by undocking (≤ 10% of wall visits end in annihilation): `exchange_rates.txt`, `wall_defect_labels_circles_ovals.png`.

**Confidence:** circles good (3 seeds, 10 coming); ovals 1–3 seeds.
**Open:** why elongation creates vacancies (candidate: the two vortices meet along the wall and push defects off).

## 8. Documents → `results/8_reports/`

| Document | What |
|---|---|
| **`master_document.pdf`** | **Master document / paper skeleton: every result, grouped by shape (passive → distributions → charge → lattice → flow → survival), every method in detail; source `reports/master/`** |
| `project_report.pdf` | Full report, all results with figures (21 pages, 5 Oct) |
| `collaborator_brief/overcharged_brief.html` | Brief for Tyler (v7); published at https://claude.ai/artifact/8yAVtVFycdqht5aJMvHKwG |
| `defect_geography/defect_geography.html` | Spatial analysis of defects in circles (30 Sep) |
| `upstream_bug_report.md` | Bugs found in NIAMH-MPCD (seed, tracker, Makefile), for upstream |
| `oval_sweep_draft.md` | Planned oval sweep, waiting for a decision |
| `cassini_patch_plan.md` | How the Cassini walls were implemented |

---

## Where everything else is

| Path | What |
|---|---|
| `RESEARCH.md` | Current status (top), hypotheses H1–H8, dated results log (newest first) |
| `TODO.md` | Code and engineering issues with urgency flags |
| `runs/` | Raw simulation data **and** the analysis scripts, by run family: `cavity/` (circles; `server_main/` = main sweep copy), `cassini/` (ovals and ellipses: `P1–P4` pilots, `W1–W4` weekend set), `activeLength/`, validation runs |
| `reports/` | Sources of the report and the brief |
| `mpcd/` | Production solver (= repo `main`) |
| `niamh-overcharged/` | Git repo (private GitHub); `main` = production, `cassini` = oval walls |
| `tools/build_results.py` | Rebuilds `results/` from `runs/` and `reports/` |

Run-family names in `runs/cassini/`: P1 ovals at half-length 2.16 and 2.88 ℓ_a (1 seed) · P2 pinched ovals · P3 small ovals 1.53 ℓ_a (3 seeds) + high k · P4 lobe-shaped oval · W1 size series 2.16 ℓ_a (3 seeds) · W2 k scan · W3 ellipses · W4 seeds 2–3 for 2.88 ℓ_a.
