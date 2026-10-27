# Draft changes and decisions for the PRL and SI (prepared 27 Oct 2026; nothing applied)

Line numbers refer to `main_v2.tex` and `si.tex` as of 27 Oct. Each item gives the current text, the proposed text and the reason. Items marked **decide** need your choice.

---

## A. Text changes (proposed)

### A1. Frank constants (main_v2 l.72, Model paragraph)

**Now:** "measured Frank constants $K_1\approx42\,\kT$ (splay) and $K_3\approx23\,\kT$ (bend)"

**Proposed:** "Frank constant $K\approx44\,\kT$, isotropic in the static response (the splay–bend split of the director fluctuations, $K_1/K_3\approx2$, is produced by the director–flow coupling of the algorithm and is absent with the coupling off; SM)"

**Why:** isotropic static response (field healing-length test, 22 Oct). The fluctuation anisotropy grows with the coupling: $K_1/K_3$ = 1.00, 1.24, 1.98 and 4.3 at $\chi_\mathrm{HI}$ = 0, 0.25, 0.5 and 1. The Jan-2024 code shows the same 1.94–1.97 (27 Oct).

### A2. Pinning by splay ≠ bend (main_v2 l.98–99, Lock paragraph)

**Now:** "With unequal splay and bend constants the commensurate ring is pinned to the sites (barrier 80--110\,$\kT$ at $k=6$), and the incommensurate one is not."

**Proposed:** "With one Frank constant the commensurate ring has no static preference for any registry; the flow fixes the tilt and the lock, Eq.~(\ref{eq:lock}), fixes the positions: the defects sit where their tilt is the flow-selected one (positions follow the tilt with slope 1.00 over 1454 runs; SM)."

**Why:** the 80–110 kBT barrier was withdrawn on 22 Oct (site selection replaces pinning).

### A3. Fig. 2(a) caption, "lines: $K_1/K_3=1.82$" (l.107)

**Proposed:** replot with one constant, $K=44$. The one-constant energies scale exactly by 44/32.5 (`ejection/ejection_K44.txt`). Change the caption to "lines: one constant $K=44\,\kT$; bands: core sizes".

**Decide:** whether to replot now or keep the anisotropic band as a sensitivity range.

### A4. Tilt paragraph opening (l.130)

**Now:** "Elasticity alone would tilt the layer by only $\approx35^\circ$ (a defect-particle model with the computed $K_1\neq K_3$ pinning)."

**Proposed:** "Elasticity alone does not select a tilt: with one Frank constant a common tilt is a zero mode of the layer."

**Why:** the 35° came from the withdrawn splay–bend pinning.

### A5. Mirror control escape ratio (l.190)

**Now:** "its escape rate is $\approx30\times$ lower"

**Result (item 3, 27 Oct; `crowd_protect/mirror10/mirror10_summary.txt`; 10 seeds, validated tracker):** at aspect ≥ 1.4 the +k layer escapes 26× more often than the −k layer (strict criterion: resident, then deeper than 10 cells; 22× with the 6-cell criterion). Circles 4.7×, aspect 1.1–1.3 12×.

**Proposed:** keep "$\approx30\times$ lower", or write "$\approx25\times$ lower (10 seeds)". The −½ wall sits 0.65–0.69 below 2|k|+2 at every aspect.

### A6. Open questions (l.202–203)

**Now:** "what pins crowded layers below the flow-aligned tilt"

**Proposed:** "what lowers the tilt of crowded layers below the flow-aligned value (the drop follows the defect depth in units of the alignment length, $d/\ell_E$, but its torque is not identified)"

**Why:** "pins" refers to the withdrawn mechanism.

### A7. Escape and crowding in the main text (Beyond circles, l.163 ff.)

The main text never quotes the old ×100. If the crowding result goes in, proposed sentence:

"How much leaks is set by the wall density: at aspect $\geq1.4$ the escape probability per wall defect falls 18-fold from $\rho=0.04$ to 0.14 (escapes from full layers 11-fold), while entries are not slowed; near the wall a resident leaves with probability $\propto e^{-U(\rho)/T}$, with $U$ the elastic pull of the layer and $T$ set by the coherence (SM)."

**Decide:** whether crowding enters the PRL (it is a strong, clean effect), or stays in the SM.

### A8. SI S7, boundary-layer correction (si l.398–403, Eq. s_tiltcorr; l.459)

**Proposed:**
- Keep the derivation of $\ell_E$.
- Replace the ansatz sentence and Eq. s_tiltcorr with: "A torque passed through the band between wall and defects would over-rotate the layer. The natural ansatz (only the vorticity acts in an anchored band) predicts $\approx14^\circ$ at $\lambda=2$ (with $\gamma_1\to\chi_\mathrm{HI}\gamma_1$ and $K=44$), a forward rotation of the band and a registry spring, and all three fail direct tests (band rotated backwards by 1–2°, overtilt 2–4°, spring $\approx6\times$ too stiff; SM Sec.~X). The overtilt is carried by the rotated neighbourhood of each defect."
- In l.459 replace "$\approx16^\circ$" with "$\approx14^\circ$ (with $\chi_\mathrm{HI}$)" and drop (iii) "splay $\neq$ bend adds +2 to +3°".

### A9. SI methods, new paragraph "Tracking" (after the observables section)

"**Tracking.** Defects are linked frame to frame (2 t.u.) per sign, by global assignment (Hungarian) against the smaller of the distances to each track's last and predicted positions (prediction: last displacement), with a 9-cell gate and up to two missed frames. Nearest-neighbour linking without prediction drops fast defects (crowded escapers move at 0.7 cells/t.u.) and undercounts escapes by up to 20× in crowded layers. Validation:
- **Synthetic trajectories:** escape probability recovered within 4%, identity switches below 1% of links.
- **Exact closure:** the tracked net charge crossing any closed curve reproduces the change of the enclosed charge, which is topologically exact, within 1.5 charges.
- **Frame rate:** results unchanged when every second frame is dropped.
- **Gate:** results unchanged between 6- and 9-cell gates."

### A10. Channels (channel paper, not the PRL)

Withdraw the mid-line topological current; see `reports/master/theory/channel/CURRENT_WITHDRAWN.md`. The SI movie note at l.1053 ("X_topological_current_W24.mp4 belongs to the channel paper") should be dropped.

### A11. SI claims built on the old 10 t.u. chain (si l.579, 581, 835) — must be revised (validated recount, `crowd_protect/recheck/recheck4_summary.txt`)

| SI line | Current claim | Validated recount |
|---|---|---|
| l.579 | Escapes "3--7 times more frequent in elongated cavities" | ×4–30 per wall defect, depending on density (circles 0.002–0.005; elongated 0.009–0.16 per 20 t.u.) |
| l.581 | "Refill is not starved: a bulk $+\frac12$ reaches the wall band equally often in every shape (ratio 0.86--1.07)" | Refill per bulk $+\frac12$ per 20 t.u. is 0.006–0.019 in circles and 0.055–0.12 in elongated cavities (5–10× higher). In elongated cavities it falls about 2× with crowding (exits fall about 18×). |
| l.835 | "Per wall defect, the tips leak 2--2.5× faster than the flanks at aspect ≥ 1.4" | By region (\|x\|/A), tips and flanks leak at similar rates (sparse elongated: tip 0.145, middle 0.163, flank 0.175). Only the top quartile of wall curvature shows a ×1.4 excess (lagged per-event analysis). |

Proposed: rewrite these three items after you confirm. The old numbers are not usable.

---

## B. Figure decisions (with my recommendation)

### B1. Fig. 4(b), candidates in `figs/fig4b_alternatives.pdf` (corrected escape counts)

**Recommend D, escape probability vs aspect, one curve per density.** With the corrected counts it shows the whole story in one panel:
- every density sits at about $3\times10^{-3}$ at aspect 1;
- the curves split about 18-fold once coherence breaks (aspect ≈ 1.3);
- coherence switches the leak on, and density sets its size.

**Alternative: E**, the vacancy map with the $C=0.85$ contour, if (a) and (b) should merge.

### B2. Fig. 1(c), candidates in `figs/fig1c_alternatives.pdf`

**Recommend B:** keep both quantities but on one partition, $Q_b=k-N_\mathrm{wall}/2$ and the vacant fraction, with no knee. The charge then adds up by construction and the caveats disappear.

**Alternative: A**, the normalised vacant fraction alone, for a single panel. Either way the caption should say that beyond $R\approx2.2\,\la$ the layer becomes diffuse (deficit 4–12 cells from the wall).

### B3. SI panels for the TODO box (`figs/si_alternatives/`)

From `fig_si_measure_alternatives`, recommend:
- (a) classification snapshot, circle;
- (c) depth histograms with the 4- and 6-cell cuts;
- (d) enclosed-charge profiles;
- (e) tilt definition.

That makes one SI figure "How we measure".

From `fig_si_layer_alternatives`, recommend:
- (a–d) the four kymographs;
- (e) static fraction vs R (circles);
- (f) static fraction vs aspect per ρ.

That makes one SI figure "The layer outside the small-circle limit", plus a paragraph on the diffuse layer (large circles, elongated cavities; crowding restores it).

---

## C. Status notes for the paper

- **K = 44 everywhere:** steady-state and escape energetics scale ×1.3–1.4 and no conclusion changes (`energetics/steadyE_*_K44.txt`, `crowd_protect/etsg_K44_summary.txt`).
- **No activity shift since Jan 2024:** the activity calibration $\ell_a(0.07)=22.2$ stands (`theory/oldcode/oldcode_test.txt`).
- **Crowding report** (`reports/crowding_report`, 11 pp) has the corrected numbers and the tracker section.
