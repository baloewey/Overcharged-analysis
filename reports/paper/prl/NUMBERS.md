# PRL draft: origin of every number (audit 2026-10-15)

Every number in `main.tex` (text and captions) with its primary source: the analysis output that produced it, the conditions, and the date it was last checked.
Paths: `T/` = `reports/master/theory/`, `C/` = `runs/cavity/`. **Bold** = changed by the 2026-10-15 audit.
Circle-sweep numbers (marked ⟳) were re-checked against the 10-seed tables on 2026-10-16 (second pass): all unchanged within s.e.; fixes applied (k = 2–10 wording, spacing range, mirror band).

| Where (main.tex) | Number as quoted | Source | Conditions |
|---|---|---|---|
| Intro, Observation, Fig. 1 caption | **0.8 ≲ R/ℓa ≲ 2.2** (full layer + bulk +1) ⟳ | `C/analysis_weekend/criteria_all.csv` (q_slow, q_dens, q_knee); RESEARCH.md 2026-10-02 20:30 and status block | ζ 0.07, λ 2, +k, ρ ≤ 0.15; layer melts at 2.87 ℓa (bulk 1.8–1.9) |
| Model | ℓa = 22.2 cells (**measured here**; method of Kozhukhov 2024) | `runs/activeLength/L200/finite_size_results.txt` (22.21 ± 0.26) | ζ 0.07, periodic box |
| Model | K1 ≈ 42, K3 ≈ 23 kBT | `reports/master/data/frank_constants.txt` (42.1/22.6, 41.9/23.5) | passive MPCD |
| Model | **R = 0.2–2.9 ℓa**, k = 2–10 ⟳ | `C/server_main/main_act0.07/*/input.json` (0.23–2.87) | sweep has k = 1–20; text uses 2–10 |
| Model | aspect 1.1–2.5, 3 sizes, 675 runs | `T/energetics/perlength.txt`; `runs/shapes_server/manifest_shape_sweep.tsv` | ovals 1.12–2.46, ellipses 1.2–2.5; R_eq 1/1.5/2 ℓa; 3 seeds |
| Observation | **61–66°** wall-defect tilt at low wall density ⟳ | `T/tiltcrowd/tilt_crowding.txt` (coherent frames, C > 0.9, 10 seeds) | λ 2, k 2–10, R 0.8–2.2 ℓa: 64.6° (ρ ≈ 0.04), 63.0 (0.065), 59.2 (0.085), 54.7 (0.10), 49.4 (0.12), 44.2 (0.14); collapses on the spacing 2πR/(2k−2). Do NOT use flow_tilt.csv tilt_flow/tilt_coh at high ρ (sign flips and scatter pull them to 38°) |
| Anchoring | mean resultant 0.99 | `T/fluct/pinning/tilt_registry_active.txt` | k 3–10, ~147k frames |
| Anchoring | barrier 80–110 kBT (k 6) | `T/fluct/anisotropy_pinning.txt` (79.1–112.4) | core 0.5–2, first order in ΔK |
| Ejection | **375–475 kBT at k 6 (200–510, k 4–8)**; **180–200** (155–200) without flow | `T/ejection/ejection_energetics.txt` ("ring frustration"; "ejection cost without flow") | one-constant + anisotropic, cores 1.5/2.5; frustration = 2k ring forced from its natural winding tilts to a common τ |
| Fig. 2(a) | energies above the natural 2k ring | same file; `prl/figs/fig2_mechanism.py` | k 6 |
| Ejection | **tilt coherence 0.90–0.98**; tumbling never ejects | `T/provenance/swirl_summary.txt`; `T/reversals/swirl_tumble.txt` | passive k 6, imposed vortex; A 0.01 (λ 2: 0.90/0.96), λ 4: 0.96–0.98; λ 0.5: 12 wall, C 0.21 |
| Ejection | "reproduces the **thresholded ejection**" | `T/shear_tilt/scan.txt`; RESEARCH.md 2026-10-09 | continuum threshold U 0.2–0.3 (no map U ↔ A): qualitative |
| Tilt angle | ≈35° from elasticity alone (**defect-particle model**) | `T/defectmodel/scan_a2.txt` | A2 5, R/ℓa 1.0, half drive (median 35°; other conditions 9–31°) |
| Tilt angle | **29° at λ 2, 22° at λ 4 (above threshold)** | `T/shear_tilt/scan.txt` (95.9−66.7, 100.8−78.9) | U = 1; at U 0.3 λ 4 changes 11° |
| Tilt angle | **63.0, 71.7, 74.8, 77.9° at k 6**; k 4, 8 at λ 3, 6 within 3° | `T/energetics/tiltlaw.txt` | runs/lambda_circle, area 2500, 3 seeds |
| Tilt angle | breakdown below λ ≈ 2 | `T/energetics/tiltlaw.txt` | λ 1.5: wall 8.9/10, bulk +1.57, C 0.75 |
| Leaking | aspect 1.25–1.37, every size | `T/shapes/phase_shapes.txt` | both families, 3 sizes |
| Leaking | **continuum: similar transition at aspect 1.4–1.6** | RESEARCH.md 2026-10-09 (continuum ellipses) | k 6, λ 2, two activities |
| Leaking | **30–45 %** inward jet within 2 cells | `T/energetics/jetmap_ff12.txt` (0.30–0.33), `jetmap_oval.txt` (0.43–0.47) | circle + ε 1.4–2.0 / ovals |
| Leaking | untilting ≈40 t.u. | `T/energetics/etsg_reg.txt`, `selfflow_ff12.txt` | aspect ≥ 1.4 |
| Shielding | **ideal ring 1.6 kBT, measured full layers ≈17 (16) kBT** inside r < 0.6R | `T/bulkdyn/shielding_layer.txt` (circle_ideal 1.6; circle_meas 15.9 median); `vacancy_field.txt` (full 17) | one-constant K 32.5 kBT; layer part of V |
| Shielding | **14–27 kBT** ideal ellipse (without/with winding flux) | `shielding_layer.txt` (14.0); `T/bulkdyn/review_ellipse.txt` (26.8) | aspect 2, 10 ideal wall defects |
| Shielding | ≈27 kBT, 7 evenly spaced | `shielding_layer.txt` row ellipse_even7 (27.5; **row added**) | aspect 2 |
| Shielding, Fig. 4(f) | **180–230 kBT** (median, inside 0.6 level) measured ellipse | `shielding_layer.txt` (175.6); `review_ellipse.txt` (226.7 with flux) | 16 measured layers, 5–9 wall defects; 270 inside the 0.8 level |
| Wall exchange | ≈140 kBT vacancy/extra field; force points at hole in every frame (28/28) | `T/bulkdyn/vacancy_field.txt` (143/140; ⟨cos⟩ +1.00) | k 6, area 2500 |
| Wall exchange | **4–5× sooner, median 22 vs 99 t.u.**; goose another defect (60 %), anywhere | `vacancy_field.txt`; per size `T/bulkdyn/ddg.txt` (20/84, 32/188, 77/338) | shape-sweep circles R_eq 1–1.5, all k |
| Wall law | c ≈ 0.012 at q ≈ 0.03; **0.61 at q 0.02, 0.10 at q 0.07**; plateau above aspect 1.5 | `T/provenance/perlength_qbins.txt` | aspect ≥ 1.5, 63 runs per bin |
| Wall law | channel releases the same share **for W ≳ 18** | `T/energetics/perlength.py` docstring (0.39–0.46 at W 20–60 vs 0.41 at q 0.031); `T/channel/channel_narrow.txt` (W 12–14: 0.0–0.32) | L 192, m 6 |
| Wall law | circles hold +1 **up to R ≈ 2 ℓa** ⟳ | `perlength.txt` (+0.99/+1.02/+1.03 at R_eq 1/1.5/2); criteria_all.csv (≈ +1.8 at 2.9 ℓa) | |
| Mirror (TODO) | **13.5 of 14 in the circle (R ≈ 1.25 ℓa), 13.0 at aspect 2** | `T/provenance/neglayer_summary.txt` (13.39–13.65; 12.86–13.12); main sweep k −6 R 1.23 ℓa, same 4-cell band: 13.47 ± 0.01 (10 seeds) | wall = within 4 cells (mirror.csv uses 3 cells: 13.1 at 1.23 ℓa; falls with R: 12.5 at 1.43, 11.2 at 1.72) |
| Mirror (TODO) | preliminary: stronger restoring force | `T/selfprop/restoring_mirror.txt` (−1.8 to −2.2 vs +k −0.0 to −0.8) | R_eq 1.27, 2 seeds |
| Intro | Klein et al.: up to ten excess +½ | literature (Klein2026) | not checkable from project files |
