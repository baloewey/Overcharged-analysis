# Research notes — topological defects in cavities with custom anchoring

Engineering and code fixes live in `TODO.md`. This file holds the questions we're asking, the hypotheses, the planned runs, and a dated log of results.

---

## ▶ Current status (rewritten 2026-10-01 17:20 office time; top block 2026-10-05) — read this first after any interruption

Previous (long) status block archived in `reports/status_archive_2026-10-01.md`; details of every item are in the dated results-log entries below.

### Where we stand (main findings so far)
- **Circles, ζ = 0.07 (working activity, decided 2026-10-01):** bulk charge +1, wall layer k − 1 (2k − 2 pinned +½), for k = 2…10 at R/ℓ_a = 1.43, 3 seeds. Clean up to R/ℓ_a ≈ 2.2; at 2.9 the plateau degrades (bulk ≈ 1.8–1.9, thick sliding layer).
- **Best observable: accumulated charge Q(c)** (charge inside a level; circle c = r). Circles: plateau at +1 just inside the wall. Edge/bulk criteria that agree on circles: pinned sites, mobility, tilt (wall +½ tilted ≈ 63° to the normal, sign = vortex sense, nearly k-independent).
- **H3:** activity does not enter only through ℓ_a. Onset effects at ζ = 0.05 (unstable layer), compressibility at ζ ≥ 0.09 (vortex empties the centre; 60–70% of central cells < 10 particles), gap closes at 0.11. ζ = 0.07 balances both.
- **H1 — size control passed (2026-10-01 20:20, layer depth d*, `runs/cassini/fig_layer_depth.png`): at half-length 1.53 ℓ_a ovals have d*/ℓ_a 0.23 → 0.50 for b/a 2.5 → 1.3 vs 0.20 for circles of the same size; geometric effect confirmed.** Earlier first look (ovals = Cassini, branch `cassini`): elongated cavities (b/a 1.6 → 1.02, half-length 2.2–2.9 ℓ_a) have **no layer + bulk structure**: Q(c) rises smoothly, no plateau, few −½; only ~4–5.5 of 12 imposed +½ behave like wall defects. Circles of the same half-length keep the plateau. All ovals fall on one curve. Two vortices per oval already at b/a = 1.6 (convex). Tentative: 1 seed, size matched by half-length.

### ▶ 2026-10-11 ~14:05 — energetics with resolved flow (runs/fineflow, flow every 2 t.u.; runs/actlayout)
- fineflow (28 depleting escapes): elastic force -10.1/-11.3/-12.4 kBT/cell at 40/20/0 t.u. before (control -11.3); cost to 6 cells 69/73/65 (control 69); flow at the escaper inward
  +0.002/+0.053/+0.234 (control +0.063), |tangential| 0.09/0.12/0.24 (control 0.14) -> same picture with resolved flow: ordinary elastic barrier, local inward (and tangential) flow burst ~4x.
- actlayout (32 events): force -13.8/-11.6/-12.4 (control -12.8); cost 73/67/60 (control 72); inward flow +0.006/+0.027/+0.112 (control +0.030).
- Lead/lag on fineflow (tiltprotect/lead_lag.py) was killed with the session (exit 137); relaunched detached with the fine-flow time series (energetics_ts.py fineflow).

### ▶ 2026-10-11 ~13:45 — ACTIVITY LAYOUT (no imposed flow): near-wall activity alone sustains the coherent state, even at eps 2 (runs/actlayout; continuum bandns_*)
- MPCD, area 2500, k 6, wall +1/2 (coherence): circle: all 9.9 (0.96, single vortex) | band 4 10.2 (0.94, circ 0.87) | band 8 10.0 (0.97) | bulk-only 7.9 (0.51, P(two) 0.80).
  eps 2.0: all 7.1 (0.44, P(two) 0.41) | band 4 10.6 (0.85, circ 0.67, P(two) 0) | band 8 8.9 (0.84, circ 0.59) | bulk-only 6.0 (0.44). Bulk charge band-only 0.33-0.42 (towards passive).
- Continuum (no swirl): eps 1.0 band 4/8: 10, single vortex; eps 2.0 band 4: 10.7, single vortex (coh 0.76); band 8: 7.65, P(two) 0.41.
- -> Near-wall activity alone sustains the single coherent vortex, the coherent tilted layer and the full 2k-2 wall, even in an eps-2 ellipse that is fully broken with activity
  everywhere; bulk activity breaks the coherence and empties the layer. Strong support for the user's reading (tilted active layer protects the coherent near-wall circulation that
  holds it; untilting = loss of protection). The elongation instability is a bulk-activity effect.

### ▶ 2026-10-11 ~13:30 — aligned time series of depleting escapes (theory/energetics/energetics_ts.py -> fig_ets_existing.png); more statistics queued
- 71 depleting escapes (area-2500 swirl + forced-vortex runs), aligned on departure (-40..+4 t.u.): cos(tau) rises first (0.25 at -40 -> 0.75 at 0; residents 0.2-0.4); inward flow at the
  escaper stays at resident level (0.03-0.05) until ~-10, then 0.14 at departure and 0.19 just after; elastic force stays restoring (-11..-12 kBT/cell, slightly stronger near departure);
  cost to move 6 cells in drops only in the last ~10 t.u. (76 -> 66 kBT). Flow frames every 20 t.u. here (ramp shape smoothed).
- Statistics (user): shape-sweep energetics generalised to any ellipse/circle (energetics_general.py; 270 runs, up to 10 events each, energy only — the sweep's flow is every 100 t.u.)
  queued after runs/actlayout (queue_general.sh -> etsg_shapes.npy); 8 more fine-flow runs (runs/fineflow2: circle, eps 1.4, 1.6, 2.0 x 2 seeds) chained after runs/fineflow.
  Server not needed (office overnight); server stays on circle seeds + mirror sweep.

### ▶ 2026-10-11 ~12:40 — energetics of ejection from measured configurations (theory/energetics/energetics.py)
- Method: depleting escapes in area-2500 runs; at departure and 20/40 t.u. before (and a resident control): all measured defects fixed, probe moved inward 0-6 cells,
  anisotropic Frank energy with eMode-2 anchoring (ellipse_binding Ellipse, h 0.5, natural relaxed cores); measured flow at the probe (flow frames every 20 t.u.).
- Imposed-swirl active ellipses (52 events): elastic force -10.2/-10.6/-12.1 kBT/cell (40/20/0 t.u. before) vs control -11.2; cost to 6 cells 73/74/71 vs 73 kBT -> NO elastic
  softening before ejection; flow at the escaper inward +0.025/+0.052/+0.149 vs control +0.056 -> a local inward flow burst at departure (3x).
- Forced-vortex runs (~10 events each): bulk-only (wall inactive) inward flow at departure +0.059 vs control +0.025 (burst present without wall activity, weaker);
  band-only +0.19 vs +0.034; cost to 6 cells at departure 51-56 vs 73-74 kBT (possibly a consequence: probe already displaced).
- Reading: ejection = local flow event against an ordinary elastic barrier; the collective 'stage' is not in the elastic landscape of the configuration but in the frequency/intensity
  of near-wall inward flow bursts (cavity flow state). Queued: energetics + lead/lag on runs/fineflow (flow every 2 t.u.) and runs/actlayout.

### ▶ 2026-10-11 ~12:30 — user's reading: untilting weakens the layer, tilted active wall defects protect it; Tests A and B (theory/tiltprotect/)
- Reading (user): band-only activity keeps tilted, coherent layers (coh 0.92) that hold; bulk-only or incoherent layers deplete; residents untilt before depleting. Protection cannot be
  the reduced inward push (bulk-only wall defects have no push yet deplete most; tilted wall layer brakes the tangential flow) -> collective: tilted layer = spiral effective boundary
  -> active stress just inside drives near-wall circulation -> keeps the layer tilted, coherent and pinned (effective-boundary feedback loop).
- Test B (site_coherence.py, shape sweep; segments = sites, windows 50 t.u.): P(segment empty next window) vs local tilt coherence with neighbours: circle-like 0.017 (coh < 0) -> 0.0012
  (> 0.9); elongated 0.047 -> 0.0065-0.011; vs jitter (pinning): circle-like 0.0001 (< 0.5 cells) -> 0.022 (2-4); elongated 0.003 -> 0.049 (2-4), 0.081 (> 4). Joint: within every
  (coherence, jitter) cell elongated 2-15x (44x in the coherent pinned cell); circle-like rates on the elongated mix predict 0.0043 vs observed 0.0151 (circle 0.0017) -> ~19% of the gap.
  -> protection by coherent tilt and pinning holds site by site, but the class gap is mostly cavity-level.
- Test A (lead_lag.py): at 20 t.u. resolution the tilt-flow cross-correlations peak at lag 0 (unresolved). Launched runs/fineflow (natural activity, area 2500, k 6: circle, eps 1.6,
  eps 2.0 x2; flow AND defects every 2 t.u., union-find tracker, 3000 t.u.); lags +-60 t.u. in 6 t.u. steps. Running also: runs/actlayout (activity layouts, no swirl), continuum bandns_*.

### ▶ 2026-10-11 ~11:55 — CLEAN MPCD forced-vortex test: bulk activity removes wall defects; wall self-propulsion not needed (runs/forced_local)
- Wall +1/2 out of 10, constant swirl 0.06, area 2500, k 6: passive 10.0 / 10.2 (eps 1.5 / 2.0); wall band only active 9.6-9.7 / 8.5-9.2; bulk only active 7.8-8.0 / 6.1-6.3;
  active everywhere (swirl_ellipse_strong) 7.5 / 7.3. Bulk charge: passive ~0.9, band-only 1.0-1.4, bulk-only 1.6 / 2.3, everywhere 1.7. -> as in the continuum: the loss needs BULK
  activity; wall-defect activity is not needed and even helps hold the layer (bulk-only loses more than everywhere).
- Depleting-escape precursors in these runs (11 events per set, indicative): head-first build-up also with bulk-only activity (cos tau 0.10/0.32/0.43/0.51 at 40/20/10/0 t.u. before,
  controls 0.29-0.53) -> the head turning is a signature of being pulled out, NOT the cause; no tangential-flow drop and LOW wall-normal flow (0.03-0.04) in bulk-only runs;
  wall-band-only runs have high v_n (0.16, generated by the active wall defects) but little depletion -> the run-level v_n correlation is not causal.
- Revised reading: no local trigger found (kick-out = artefact; elastic/flow perturbations by bulk +1/2 not decisive; head-first and v_n = correlates). Working hypothesis: a STEADY-STATE
  SHIFT — bulk activity in an elongated cavity favours a larger bulk charge (continuum mirror also raises it); escapes implement the shift (cavity-level propensity, Poisson-like
  within a run with a restoring tendency). To test: bulk charge vs bulk activity with a passive wall band; partial uftrack analysis tonight.

### ▶ 2026-10-11 ~11:40 — DEPLETING escapes (theory/kickout/depleting.py -> depleting.txt)
- Depleting = resident (in the band >= 40 t.u.), layer at or below 2k-2 (smoothed), no swap, vacancy persists 10-50 t.u. after. ~11-12% of all escapes in both classes;
  rate 0.043 (circle-like) vs 0.183 (elongated) per wall defect per 1000 t.u. (x4.3).
- Precursors at the escape frame: cos(tau) 0.75 for depleting escapes vs 0.31-0.36 other escapes and 0.31-0.46 resident controls -> residents leave HEAD-FIRST (self-propulsion
  pointing into the bulk); local flow |v_n| 0.07-0.08 vs 0.054-0.067 and |v_t| 0.15 vs 0.17-0.20 (escapes in places with more wall-normal and less tangential flow);
  nearest bulk +1/2 heading away (-0.06..-0.14) vs towards for controls.
- Working mechanism (user's combination): tangential near-wall flow keeps wall defects tilted by flow alignment (self-propulsion mostly along the wall); where it weakens a defect
  untilts, heads into the bulk and propels itself out; bulk +1/2 heading at it would push it back. Elongated cavities: more normal, less tangential near-wall flow (run-level R^2 0.91).
  PRECURSOR CONFIRMED (depleting2.txt): cos(tau) of the escaper 40/20/10 t.u. before and at the escape: circle-like 0.48/0.58/0.65/0.75 (controls 0.43-0.46),
  elongated 0.29/0.49/0.63/0.75 (controls 0.24-0.31); non-depleting escapes lower than controls. -> residents untilt over 20-40 t.u. (head turns into the bulk) and then leave.
  Untilted = the flow-free resting state (tails at the homeotropic sites, cos tau = 1): flow-induced tilt redirects self-propulsion along the wall and protects the layer.
  Consistent: passive and -1/2 layers hold; continuum band-only activity + swirl holds. Not explained: continuum bulk-only activity still loses defects (different route).
  MISSING LINK, first evidence (theory/kickout/depleting_flow.py on the imposed-swirl active runs, flow every 20 t.u., 73 depleting escapes): local |v_t| at the site
  40/20/10/0 t.u. before: 0.24/0.21/0.17/0.13 (controls 0.18 flat) while cos(tau) 0.26/0.51/0.66/0.77; |v_n| 0.10-0.12 (controls 0.08-0.09) -> the flow along the wall
  drops, the defect untilts, then leaves. Small sample; full test on runs/uftrack (natural active flow, flow every 20 t.u., union-find tracker) tonight.

### ▶ 2026-10-11 ~11:00 — cavity-level predictors: wall-band flow explains the elongation effect (theory/kickout/runlevel.py -> runlevel.txt)
- 660 runs, log(escape rate per wall defect): elongation alone x4.9 (R^2 0.13); + rms wall-NORMAL flow at depth 2-8 (vn): elongation x0.3, vn +1.77/s.d., R^2 0.49; tangential flow
  in the band -1.46/s.d. (stabilising); k -1.41/s.d. (denser layers leak less). All predictors without elongation R^2 0.91; adding elongation: x0.9 (none left).
  -> the cavity-level property behind the leak is the near-wall flow: more flow towards/away from the wall, less along it, in elongated cavities. Caveat: escaping defects
  stir the flow themselves (reverse causality possible).
- User's point (11:00): focus on DEPLETING escapes (resident defect, layer at or below 2k-2, no swap, vacancy persists); steady-state escapes (revolving door, swaps, visitors
  leaving a crowded layer) are different dynamics. Analysis: theory/kickout/depleting.py.

### ▶ 2026-10-11 ~10:45 — approaching +1/2: elastic and flow models vs data; escape-factor decomposition (theory/kickout/elastic_kick*, flow_kick*, approach_range*, decomp*)
- Elastic barrier (anchored disk, k 6, 10-ring + bulk +1; barrier to depth 5, reference approacher far: 33 kBT one-constant / 37 anisotropic): approacher in front at 3-4 cells -9..-21;
  5 cells ~0; 6-14 cells +4..+11 (repulsion blocks the path); SIDE (depth 6, half a site along the wall) -19..-29 kBT (crowding: the layer + approacher ~ 11 -> one is cheap to expel).
- Flow (director-texture Stokes, zeta/eta 0.5): approacher in front adds inward velocity -0.12 (3 cells) ... -0.02 (14 cells) at the wall defect (pushes it back towards the wall,
  ~1/r); from the side +0.01. Approacher orientation fixed by the static texture (caveat). Wall defect's own self-propulsion in the static texture: +0.19 inward.
- Data (self excluded, region-matched controls): other bulk +1/2 DEPLETED around wall defects before escapes (ratio 0.3-0.8 within 8 cells, 1.1-1.2 beyond 10), mostly in front
  (depth 4-8, +-6 along: 0.2-0.6); side positions (9-12 along, depth 4-6) slightly enriched in circles (1.2-1.4), not elongated (~0.95); nearby ones point away (heading cosine
  r < 5: escapes 0.16 / -0.02 vs controls 0.56 / 0.41) -> bulk +1/2 heading at a wall defect go with it STAYING (consistent with the flow push-back), no kick-out.
- Decomposition (5.4M samples, pooled logistic, class gap 3.6x odds): depth +24%, orientation 0, nearest bulk +1/2 +1%, bulk +1/2 count 0, region -1%, coherence -38%, size 0;
  all local together +3%. Cavity vacancies +82%, neighbour gaps +65% (all +61%) — but vacancies are an outcome of escapes and do not raise the hazard within a run (avalanche
  analysis) -> escape propensity is a cavity-level property, not driven by the local environment of individual defects. Next: local flow at the wall defect (uftrack runs, flow every 20 t.u.)
  and the +1 pair's positions (identities with the union-find tracker).

### ▶ 2026-10-11 ~10:20 — loss channels; the 'kick-out' signal was an artefact (theory/kickout/channels.py, kick_identity.py)
- Loss channels (per wall defect per 1000 t.u.): escape 0.2-0.5 (circle-like) / 1.5-1.7 (elongated); annihilation with a bulk -1/2 at the wall 0.02-0.05 / 0.13-0.43 (8-25% of escape,
  growing with size); untracked 0.5-1.6 (legacy-tracker losses, a caveat). Swaps (a bulk +1/2 takes the vacated site within 20 t.u.) 13-43% of escapes.
- channels.py reported bulk +1/2 within 3 cells in the 20 t.u. before 29-62% of escapes vs 8-16% control ('kick-out', cf. Kozhukhov, Loewe, Thijssen, Shendruk arXiv:2503.19664).
  WITHDRAWN: with identities (kick_identity.py, merged split defects, gap bridging), the closest 'bulk +1/2' was the escaping defect itself at an earlier frame. Excluding it,
  another bulk +1/2 within 3 cells before an escape: ~0%, same as control. Instead: 21-61% (circle-like) / 22-39% (elongated) of escapers were themselves in the bulk within the
  previous 40 t.u. (random wall defect 6-16%): a revolving door of recent arrivals. Bulk -1/2 essentially never approach (0-3%). Identity classes (dynamical/persistent) not
  resolvable with the legacy tracker (median identity life ~3 frames); redo with union-find-tracker data (circle seeds s04-s10, mirror sweep) if needed.

### ▶ 2026-10-11 ~09:55 — clean MPCD forced-vortex test launched (runs/forced_local, bin/mpcd_wallact2.out)
- Constant imposed swirl 0.06 (as swirl_ellipse_strong, whose active-everywhere runs are the reference: wall 7.2-7.5 at eps 1.5, 7.3 at eps 2.0), area 2500, k 6:
  activity only in the bulk (wallActOff +4) and only in the wall band (wallActOff -4), eps 1.5 / 2.0, 2 seeds each; passive reference at the same swirl (1 per eps). 10 runs, ~2 h.
  Continuum analogue: band-only keeps 10, bulk-only loses 2. Analysis: overnight/counts.py forced_local + selfprop/restoring_swirl.samples2.

### ▶ 2026-10-11 morning — localized activity: wall-defect self-propulsion is NOT what removes the layer (ellipse_runs/summarize2.py; runs/wallact/compare.py)
- Continuum, matched flow: uniform zeta 0.75 (v_rms 0.65-1.12, matching the 4-cell knockout's 0.67-1.19): vacancies 0/0/1/1/1.6/3.1 vs knockout 0/0/1/2/2.6/4 -> at equal flow the
  knockout loses MORE defects; narrow knockout (2.5 cells): 0/0/1/1/2/3 (~ matched flow). No sign that wall activity removes defects.
- Continuum, forced vortex (swirl 0.1) + localized activity, eps 1.0/1.6/2.0: passive 10/10/10; wall band only active 10/10/10; cores only active 10/10/10; bulk only active 10/8/8;
  active everywhere 10.5/10/8 -> defect loss needs BULK activity (with elongation); wall/core activity is neither necessary nor sufficient.
- Continuum size test (R_eq 47.7, k 9): eps 1.0/1.2 circle state; eps 1.4 single vortex (coh 0.98) but 1.7 vacancies; eps >= 1.6 broken (two vortices); vacancies 1.7 -> 4.6,
  bulk charge -> +3.2: transition aspect about the same as at R 31.8, more vacancies at larger size (as in MPCD).
- MPCD knockout (wallActOff 4, paired seeds): CIRCLE loses its state (wall 9.8-10 -> 5.3-5.6, bulk +0.8 -> +2.5, circulation 0.86 -> 0.0-0.08, two vortices 70%, pull-back -2.4 -> -0.6):
  near-wall activity is needed for the coherent single vortex in MPCD -> the knockout is far from clean. Ellipses eps 1.6/2.0: wall 8.7 -> 7.6, 8.25 -> 6.9 (more loss), bulk +2 -> +2.5..2.8,
  pull-back slightly stronger (-0.6 -> -0.8..-1.1). -> knockout does not restore the layer; not interpretable as a self-propulsion test.
- Net: the self-propulsion hypothesis (of 2026-10-09 14:30) is not supported; in the continuum the loss is tied to bulk activity under elongation. Not via the vortex number: forced single vortex + bulk-only activity still loses 2 (P(two) 0-0.08), consistent with MPCD (leakiness the same in
  one- and two-vortex intervals). Open: the MPCD mirror (-1/2 hold,
  bulk charge falling). Clean MPCD analogue proposed: imposed swirl + activity only in the bulk (wallActOff +4) vs only in the wall band (wallActOff -4), bin/mpcd_wallact2.out.

### ▶ 2026-10-10 ~21:45 — the knockout is not clean; controls set up (user's objection)
- Continuum knockout (activity off within 4 cells, zeta 1.0) vs base: interior v_rms x0.71-0.80, circulation 0.92 -> 0.81 (circle), two-vortex fraction at eps 2.0 1.0 -> 0.0
  -> the knockout changes the flow, so 'knockout does not save the layer' is confounded (in the continuum, and likely in MPCD).
- Controls launched/queued (continuum, ellipse_runs/; summary summarize2.py now prints v_rms): (1) matched flow: uniform zeta 0.75 (match_*), compare layer at equal interior flow;
  (2) narrow knockout d = 2.5 (narrow_*; wall-defect cores at ~1.4 cells), queued; (3) user's idea, forced vortex + localized activity (swirl 0.1 within 0.6 B fixes the flow):
  passive / wall band only active (--bandonly 4) / cores only active (--actmode core: S < 0.6 within 4 cells, dilated) / bulk only active / active everywhere, eps 1.0/1.6/2.0
  (sw_*, 15 runs, swirl_queue.sh, max 12 concurrent). Short tests (eps 2.0, 300 t.u.): passive+swirl ejects to 10 (circ 0.74); band-only and core-only + swirl keep 10.
- Planned checks for every knockout (MPCD too): interior v_rms, circulation, P(two), bulk defect density, shear at the defect depth, tilt coherence (clean if within ~10%), and the
  self-propulsion signature (inward drift vs cos tau) must vanish.
- MPCD: wallActOff < 0 = activity only within |d| of the wall (branch anchmod, bin/mpcd_wallact2.out, off = bit-identical; swirl included) for the forced-vortex version.
  Core-only activity in MPCD possible via the cell S (computed every step) but needs threshold tests -> after the continuum result.

### ▶ 2026-10-10 evening — continuum knockout/mirror: self-propulsion NOT needed in the continuum; bulk Leslie angle at lam 2 follows lam_eff + 1.8 deg
- Continuum ellipses, zeta 1.0 (ellipse_runs/summarize2.py): wall band inactive (zeta = 0 within 4 cells): vacancies 0/0/1.0/2.0/2.6/4.0 (eps 1.0/1.2/1.4/1.6/2.0/2.5) vs base
  0/0/1.0/1.36/1.64/4.0 -> knocking out the wall defects' self-propulsion does not save the layer. Mirror k -6: wall -1/2 14.0 (eps <= 1.2; = 2|k|+2, the negative-layer prediction
  realised), 15-16 above the transition (gains -1/2); bulk charge +1 -> +1.5..+2 in ALL variants above the transition.
  [CORRECTED 2026-10-11: not 'goes with the two-vortex flow' — the forced-vortex runs lose defects with a single imposed vortex (P(two) 0-0.08); the two-vortex
  mechanism stays discarded, as in MPCD. What holds: bulk activity in an elongated cavity removes wall defects, independent of the vortex number.]
  CONFLICT with MPCD mirror P6 (wall -1/2 nearly full, positive bulk charge falling 0.6-0.7 -> 0.46 with elongation; strong pull-back). The continuum has little noise and one
  constant: its loss is a change of steady state, not escapes against noise -> MPCD wall-band knockout launched (runs/wallact, below) to decide for MPCD.
- Channel lam 2 (runs/leslie/leslie.txt): theta_L 32.8 (U 5, S 0.56), 29.9 (U 10, 0.88), 31.1 (U 20, 0.95); particle-averaged 30.8 / 28.4 / 29.2 -> MPCD bulk alignment follows the
  lam_eff S-dependence (non-monotonic, minimum near S 0.88) with a constant offset +1.8 deg; NOT bare lam. Cavity tilt vs 2 theta_L(bulk): 62.2 vs 65.6, 62.0 vs 59.8, 64.8 vs 62.2
  -> tilt = 2 theta_L(bulk) within +-3 deg; residual -3.4 / +2.1 / +2.6 grows with S (compensation picture).
- Launched: MPCD wall-band knockout (runs/wallact, bin/mpcd_wallact.out, wallActOff 4): ellipse eps 2.0 and 1.6 (Req 1.5, k 7) + circle (Req 1.5, k 6), s01-s02 = same inputs and
  seeds as shapes_server (paired). Continuum size test (ellipse_runs/big_*): R_eq 47.7 (1.5x), k 9, zeta 1.0, eps 1.0-2.5.

### ▶ 2026-10-10 ~11:10 — self-propulsion test launched (continuum) and MPCD localized activity implemented
- Continuum (ellipse_active.py, zeta 1.0, eps 1.0-2.5): (i) --nowall 4: activity off within 4 cells of the wall, ramped to full by 6 (bulk still active, flow still driven);
  (ii) --k -6 (mirror: wall -1/2, not self-propelled). 12 runs, ~3 h; summary ellipse_runs/summarize2.py. Prediction if self-propulsion of wall +1/2 removes the hold:
  (i) and (ii) keep their wall layer in elongated cavities while the base k 6 loses 1.4-4 defects.
- MPCD (user: localized activity also possible in MPCD): new input `wallActOff` on branch anchmod (uncommitted): activity scaled by clip((depth - d)/2, 0, 1) per cell before the
  collision, depth = min over inverted 4-fold BCs of -W/|grad W| (A = 1/aInv; exact near a circular wall). Binary bin/mpcd_wallact.out (clean build). Validation
  (runs/wallact_validation): wallActOff 0 bit-identical to mpcd_swirl.out; wallActOff 100 -> flow 0.067-0.081 vs passive 0.064-0.067 vs active 0.23 (the exact centre cell keeps
  activity because grad W = 0 there; irrelevant in real use). docs/InputGuide.md updated. Ready for an MPCD knockout run in ellipses once the continuum result is in.
- Paper outline drafted: reports/paper_outline/outline.md.

### ▶ 2026-10-10 morning — continuum ellipses reproduce the shape transition; bulk Leslie angle in MPCD depends on S (ellipse_runs/summarize.py; runs/leslie/leslie.txt)
- ACTIVE continuum model in ellipses (one-constant, k 6, lam 2, same area as R 31.83, noise 0.005, 12000 t.u.): eps 1.0/1.2 circle state at both zeta (10 wall, bulk +1, coh 1.0,
  tilt 61-62); eps 1.4: circle state at zeta 0.3 (tilt 60, circ 0.86), BROKEN at zeta 1.0 (9 wall, two vortices, coh 0.58); eps >= 1.6 broken at both (two vortices, coh 0.48-0.54,
  circ 0.18-0.31, wall 9 -> 6, bulk +1.5 -> +3, more vacancies at larger eps and zeta). Transition eps ~1.4-1.6 (MPCD 1.25-1.37), lower at stronger activity.
  -> the shape story (loss of the single coherent vortex -> layer loses defects to the bulk) holds outside MPCD.
- Bulk Leslie angle (Poiseuille channel, no anchoring): lam 1.2: 17.4 (U 5, S 0.56), 14.2 (U 20, S 0.95), U 10 half not aligned (~15-16 in the other half); lam 1.5 U 10: 22.3;
  lam 2 U 10: 29.9. Bare: 16.8 / 24.1 / 30.0; particle-averaged lam_eff: 18.9 / 14.4 / 12.0 (U 10) / 21.5 / 28.4. -> near tumbling MPCD alignment DOES depend on S roughly as the
  particle-averaged law; at lam 2 the measurement sits at the bare value but the predictions differ by only 1.6. The S-independence of the cavity tilt is then most likely a
  compensation (bulk alignment weakening with S vs a correction growing with S, e.g. splay != bend), not bare-lam alignment. Running: lam 2 channel at U 5 and 20 (runs/leslie/jobs2.txt)
  -> cavity correction = tilt - 2 theta_L(bulk) at each S.

### ▶ 2026-10-09 22:45 local — why does the MPCD tilt follow the bare lambda? Solver reading + bulk Leslie-angle runs (runs/leslie/)
- Solver (mpc.c ~4700, lc.c): per step LCcollision (each particle: Metropolis with independent uniform trial angles, 5 + 0.1 U S iterations, started from its previous angle;
  no acceptance -> keeps its old angle; distribution centred on the current cell director DIR from localPROP) THEN jefferysTorque (every particle rotated by
  shearSusc * dt * Larson-Jeffery rate with lam), then streaming. DIR of the next step = eigenvector of Q of the rotated particles. On paper this gives the moment law
  lam_eff = lam (1 + S4)/(2 S2) -> possibilities: (a) MPCD bulk alignment follows lam_eff and in the cavity a correction growing with S (K ~ S^2: splay != bend, boundary layer)
  cancels the S dependence by coincidence; (b) something in the algorithm gives the bare lam.
- Direct test: Poiseuille channel (grav 0.001, walls at y 4 and 64, MUN = MUT = 1 = no anchoring, 64 x 68 periodic in x), passive; aligned angle vs shear sign per row.
  lam 1.2 at U 5/10/20 (bare theta_L 16.8; lam_eff predicts 18.4 / 11.0 / ~13: discriminating), lam 1.5 and 2.0 at U 10. Analysis runs/leslie/leslie.py. Test: no-slip profile OK.

### ▶ 2026-10-09 22:15 local — active continuum model on ELLIPSES launched (theory/shear_tilt/ellipse_active.py, ellipse_runs/)
- Q-tensor + Stokes + backflow + extensile activity in ellipses of the reference circle's area (R_eq 31.83), arclength-uniform anchoring (eMode 2 analogue), k 6, lam 2, noise 0.005;
  eps 1.0/1.2/1.4/1.6/2.0/2.5 x zeta 0.3/1.0, 12000 t.u. (12 single-core runs, ~2-3 h). Check eps 1.0 (short): 10 wall + bulk +1, single vortex (circ 0.91), coherent tilt 62.
  Question: does the continuum model also lose coherence and leak wall defects above aspect ~1.3 (MPCD transition 1.25-1.37)? Summary: ellipse_runs/summarize.py.

### ▶ SUMMARY for the next session (written 2026-10-09 ~16:45 local) — read this first
- **Tilt angle:** flow alignment in the no-slip shear at the wall; leading order arccos(1/lam) = 60 deg (lam 2) vs measured ~62 (active) / 61-65 (passive swirl). The tilt is
  independent of the order parameter (S 0.56-0.94) and density (10/20/40) -> the particle-averaged lam_eff and the core-order correction are WITHDRAWN (correction 16:00).
  Residual 2-5 deg: splay != bend (+2..3 in the disk model); boundary-layer ansatz does not transfer to MPCD (gamma1 ~ 60 would give +16).
- **Active continuum model** (Q-tensor + backflow + activity, no imposed flow) self-organises: vortex -> ejection 12 -> 10 + bulk +1 -> tilt 61.7-62.1 deg at lam 2 (MPCD ~62).
- **Shapes:** wall defects lose their hold (restoring force) only with activity AND incoherent flow; passive ellipses at coherence 0.34 hold; active ellipses at 0.35-0.66 do not, even
  under 10x imposed swirl; non-self-propelled wall -1/2 (mirror) hold in incoherent active ellipses -> self-propulsion pushes +1/2 off the wall. Depth/orientation explain ~21%.
- **Numbers:** gamma1 ~ 60 (splay 62 +- 2, bend 56 +- 2); free-defect mobility 0.0054 (bulk) vs 0.02 (wall).
- **Server:** circle seeds s04-s10 running (169/2170 at 17:20 UTC Oct 9); mirror shape sweep (675) chained after it (chain_after_seeds10.sh, pid 1012571).
  Planned (to discuss): seeds s04-s10 for both windings in shapes.
- **Documents:** reports/theory_report/theory_report.pdf (for Tyler, 12 p., corrected) and master.pdf (100 p.) are current. Office machine: all my runs finished.

### ⚠ 2026-10-07 — SIGN ERROR in all Stokes/feedback work (read first)
- The defect angle in defects.dat is the direction of ∇·Q (lc.c `topoAngleLocal`, Vromans–Giomi) = the **tail** of a +½. Bulk +½ move against it (−0.56…−0.70 cells / 2 t.u., 66–72% of steps, ~60σ; `reports/master/theory/fluct/pinning/head_convention.py`) → extensile; head, self-propulsion and active force are along −angle.
- All +½ point forces (stokes_circle.wl, stokes_ellipse_static.wl, stokes_frames/solve.wl + fit.py, stokes_tilt_modes.wl, stage2) used +angle. Corrected: wall layer forces OPPOSE the circulation (non-negative fit gives it 0); bulk +½ drive it (corr +0.5 circles) but explain only 6–21% of the variance and do not predict the ellipse vortex number; feedback is negative (wall damps all tilt patterns, one-vortex most).
- Wall +½: tails at the wall, heads into the bulk leaning **upstream** (fits the observed upstream drift). "Tilt = vortex sense / lean downstream" statements refer to the tail.
- Master doc: correction box at the top of §hydro; stage-2 note. Next: Stokes with the full active force −ζ_eff ∇·(Q) from the measured director field (no defect idealisation), decomposed into wall band vs bulk.

### ▶ 2026-10-09 ~16:40 local — gamma1 final; density 40 (runs/rotvisc/rot_visc_modes.txt; theory/tiltdiag/tiltdiag_d40_U10.txt)
- gamma1 from director-mode decay, 4 seeds x 5000 t.u., modes m 2-4: splay 61.9 +- 2.3, bend 56.0 +- 2.0 kBT tau/cell^2 (12 estimates each) -> gamma1 ~ 60 (splay ~ bare; bend lower by backflow).
  Consistent with the free-defect mobility 0.0054 via mu ~ 4/(pi gamma1 ln(l/a)) for ln(l/a) ~ 4. The flat-wall boundary-layer ansatz with gamma1 60 (shearSusc 0.5) would add ~16 deg
  in MPCD (measured residual 2-5) -> its prefactor does not transfer (as stated).
- Density 40: tilt 62.6 (coh 0.92), S profile identical to density 20 -> tilt independent of density (10/20/40: 62.4 / 61.5-62.2 / 62.6).

### ▶ 2026-10-09 ~15:25 local — ACTIVE continuum model reproduces the active state (shear_tilt/scan_active2.txt); free-defect mobility (mobility/pair_mobility.txt)
- Q-tensor + Stokes with backflow and extensile activity (one-constant, no imposed flow), circle k 6 R 31.83, start from the passive 12-ring:
  zeta 0.1: vortex (|circ| 0.92), no ejection (12), coh 0.47, tilt 68; zeta 0.3 / 1.0 (lam 2): vortex + ejection 12 -> 10 + bulk +1, tilt 61.7 / 62.1 (MPCD active ~62),
  independent of activity; lam 4: 76.5; lam 1.2: 37.1. Passive + backflow + body force 0.03/0.1: ejection, tilt 65.5/65.8 (lam 2), 77.6 (lam 4).
  -> the full chain (activity -> vortex -> ejection -> 2k-2 + bulk +1 -> tilt ~ arccos(1/lam) + 2) self-organises in a continuum active nematic. (v1 failure was too low zeta.)
- Free +-1/2 mobility in the bulk (quench_L96_s1, 54 annihilations): r^2 linear in t_a - t, s = 1.11 cells^2/t.u. -> mu = 0.0054 (local slopes 1.30/1.09/1.11 for r 3-5/5-8/8-12),
  4x below the wall-defect value 0.02. Seed 2: 48 annihilations, s = 1.12 -> mu = 0.0055 (agrees).

### ▶ 2026-10-09 ~16:00 local — CORRECTION: the tilt does not follow the order-dependent alignment (theory/tiltdiag/tiltdiag_*.txt)
- Particle-level moments (base_s1/s2, trajOut): S4 = Maier-Saupe within 0.01 in bulk and layer; in cores slightly higher (S2 0.63-0.66, S4 0.27-0.30) -> lam_eff 1.98-2.01 (MS 1.92-1.95).
- Order varied with mfpot: measured tilt 62.2 (U 5, bulk S 0.56, core 0.40), 61.5-62.4 (U 10, S 0.87, core 0.65; also density 10), 64.8 (U 20, S 0.94, core 0.72).
  Jeffery-moment predictions: bulk-S 61.5 / 56.7 / 58.2, core-S 68.5 / 58-59 / 57.2 -> wrong trend (core version anti-correlated). The tilt is nearly independent of S
  (0.56-0.94) and density; best described by the BARE arccos(1/lam) = 60 deg + 2-5 deg. -> WITHDRAWN: the reduced order near the cores as the source of the residual, and
  the moment lam_eff as the leading term (the earlier "56.5 -> 60 -> 62-63" stack). Residual candidates: K1 != K3 (+2..3 in the model), boundary-layer correction.
  Why the moment law fails in MPCD: open (the LC collision resampling after the Jeffery step may reset the distribution differently than the moment calculation assumes).
- Wall density: only the half cell at the wall is depleted (~78%); no density effect on the tilt.
- Free-defect mobility (partial quench s1, 52 annihilations): r^2 linear in t_a - t (slope 1.07, stable over r 3-12) -> mu_bulk ~ 0.005, 4x below the wall value 0.02 (final pending).

### ▶ 2026-10-09 ~14:45 local — tilt diagnostics, first two runs (theory/tiltdiag/tiltdiag.py -> tiltdiag_<run>.txt)
- d10_U10 (density 10): bulk S 0.874, tilt 62.4, coh 0.87, wall 10.2. d20_U20 (mfpot 20): bulk S 0.943, tilt 64.8, coh 0.94. Particles per cell: only the half cell at the wall depleted
  (~80%), uniform beyond 0.5 cells -> wall density fluctuations a small effect. S near cores 0.67-0.68 (d10), 0.71-0.72 (U20). Higher S -> slightly larger tilt, as the Jeffery law
  predicts (lam_eff/lam minimum near S 0.85). Baselines with particle-level S4 pending.

### ▶ 2026-10-09 — PLANNED (user, to discuss): production seeds for shapes
- After the circle seeds s04-s10 and the 3-seed mirror shape sweep: 7 more seeds (s04-s10) for BOTH windings in ellipses and ovals (aim: answer small-N reviewer complaints).
  Mechanics: gen_shape_sweep.py --nseeds 10 / gen_shape_mirror.py --nseeds 10 write only the missing seed dirs (per-point seed numbering), same binary.
  Cost ~2 days per 3-seed sweep -> ~9-10 days at 96 (trim to 3 densities: ~6 days). Check disk first (~60-70 GB per 675 runs compressed?; 344 GB free on 2026-10-09).
  Open: full vs trimmed densities; order relative to the mirror first pass.

### ▶ 2026-10-09 ~17:40 UTC — MIRROR SHAPE SWEEP prepared and chained on the server (user decision: full set, after the circle seeds)
- Server ~/Overcharged/runs/shapes_mirror/: 675 runs = the positive shape sweep with k -> -k (same shapes, R_eq, |k|, steps, outputs; 3 seeds; seeds 81000000 + 100 i + s);
  generator gen_shape_mirror.py (uses gen_shape_sweep.py; local copy runs/shapes_server/), jobs_shape_mirror.txt (largest first), manifest_shape_mirror.tsv (k negative),
  binary ~/Overcharged/niamh-overcharged-cassini/mpcd.out (same as the positive sweep), run_queue.sh copied. 5.56e12 particle-steps (~2 days at 96).
- Validation (server, ~/Overcharged/runs/shapes_mirror_test/, passive 4000 steps): circle k -6 -> 12 -1/2, ellipse eps 2 k -3 -> 6, oval b/a 2 (eMode 3) k -4 -> 8; no +1/2.
  First check of Cassini ovals with negative winding.
- Chain: chain_after_seeds10.sh (pid 1012571, detached) waits for run_queue.sh jobs_seeds10.txt (pid 956445) to exit, then runs the mirror queue at 96; log chain.log,
  queue log queue_shape_mirror.log. If the seeds10 queue is stopped by hand, the mirror sweep starts anyway — kill pid 1012571 first if that is not wanted.

### ▶ 2026-10-09 ~14:50 — wall -1/2 keep their hold in active ellipses (selfprop/restoring_mirror.py -> restoring_mirror.txt; P6 mirror shapes, 2 seeds, A 2500)
- Pull-back at 2.5-3 cells: active k -6: circle -2.2, eps 1.2 -2.1, eps 1.5 -2.0, eps 2.0 -1.8..-2.0 (k -12 stronger); passive references -2.6..-3.3. Positive winding (+1/2) in active
  ellipses: -0.0..-0.8. -> non-self-propelled wall defects keep their hold in elongated active cavities; supports self-propulsion as the outward push (hypothesis of 14:30).
  Coherence checked (overnight/P6_mirror_flow.txt): mirror ellipses eps 1.5/2.0 have incoherent flow (circ/vrms -0.17..+0.10, P(two) 0.46-0.78), like positive elongated
  cavities, yet wall -1/2 13.0-13.4 of 14 and strong pull-back -> activity + incoherent flow removes the hold only for self-propelled (+1/2) wall defects.
  Remaining caveats: one size (R_eq 1.27 l_a), 2 seeds (the 675-run mirror sweep is chained on the server).
- Proposed (user, to run after the server circle seeds): full mirror of the 675-run shape sweep (k -> -k, 3 seeds), ~2 days on the server; trimmed option 405 runs (3 densities).

### ▶ 2026-10-09 ~14:30 — restoring force needs activity AND incoherence to collapse; combined tilt calculation matches (selfprop/restoring_coh.txt, restoring_swirl.txt, restoring_active_2tu.txt; shear_tilt/disk_mpcd.txt)
- Pull-back force -U'/kT_eff at 2-3 cells (2 t.u. steps): PASSIVE circles (swirl, coh 0.90-0.96) -1.4..-3.5; PASSIVE ellipses (swirl, coh 0.34-0.72) -2.3..-3.4 (even eps 2 at coh 0.34);
  ACTIVE circles (coh 0.93-0.98) -1.5..-3.5; ACTIVE elongated (shape sweep, coh 0.38-0.67; swirl ellipses up to A 0.1, coh 0.35-0.66) -0.0..-0.8 (high k 10-11: -1.4..-2.1).
  Shape sweep by interval coherence: elongated intervals with coh > 0.7 still weak (well <= 0.45 vs 1.8-2.5 circle-like). -> the hold is lost only with activity AND incoherent flow.
  Reading (hypothesis): self-propulsion pushes wall +1/2 off the wall in active cavities; only the coherent pinned ring of a circle binds enough to resist; passive cavities have no
  self-propulsion and hold at any coherence; active ellipses cannot reach the coherence (>0.9) that would hold them, even at 10x imposed swirl.
- **One combined tilt calculation (disk Q-tensor + measured no-slip flow + K1/K3 1.82 + MPCD Jeffery alignment law lam_c(S) = lam(1+S4(S))/2 point by point):** after ejection
  lam 2: 62.8-65.9 deg (MPCD 61-64), lam 4: 78.8-80.4 (MPCD 75-83), lam 1.2: 39.8 ejected / 58-68 not ejected at coherence ~0.5 (MPCD ~57 at coherence 0.54, only partly ejected).
  Replaces the additive estimate; lam 2 and 4 inside the measured ranges.
- Session crash ~14:00: detached jobs unaffected (tiltdiag 1/6 done, rotvisc running, scan_active2 running, driver_oct09.sh waiting).

### ▶ 2026-10-09 afternoon (4) — the restoring force collapses in elongated cavities (theory/selfprop/restoring.py -> restoring.txt); K1 != K3 in the tilt (shear_tilt/flatwall_aniso*, disk_aniso*)
- Fokker-Planck from 10 t.u. depth steps of 5.4M wall +1/2 samples: effective force -U'/kT_eff. Wall side identical (d 0.25: +5.7..+6.3 circle-like, +6.4..+7.0 elongated).
  Bulk side collapses: at d 2.75 circle-like -0.94..-1.22, elongated -0.07..-0.22; mean drift back to the wall there -0.54..-0.64 vs -0.08..-0.17 cells/10 t.u.; noise similar.
  Well (minimum -> 3.75 cells): circle-like 1.6-2.5 kT_eff, elongated 0.36-0.73 (tips 0.42, middles 0.52, flanks 0.71). -> weaker holding, not a stronger push.
- K1 != K3 (cubic L3 term, gradient-checked; K1/K3 = 1.82 at the same mean K): FLAT WALL: breaks the zero mode; the local anisotropy prefers a large tilt (rest 79 deg; reversed ratio 8 deg)
  and adds +6..+8 deg at lam 2 (73-74 vs 65-68) and +4..+5 at lam 4 (83-84 vs 79-80) under MPCD-like shear. Missing there: the global cost of the bulk +1 spiral, which pulls
  the tilt back (pins the active circle at tilt 0 without flow). DISK model with anisotropy (both effects; disk_aniso.txt): K1/K3 1.82 adds only +2..+3 deg (lam 2: 66.6/69.3 vs 64.8/66.7; lam 4: 80.1/81.0 vs 78.2/78.9; lam 1.2: 48.9 vs 46.3) —
  the bulk-spiral cost cancels most of the local +6..+8. Stacked estimate for MPCD lam 2: 56.5 (arccos(1/lam_eff)) -> ~60 (core order) -> ~62-63 (K1 != K3) vs measured 61-64;
  lam 4: ~77.5 vs 83; lam 1.2 (near tumbling) still fails (~38 vs 57). Additivity assumed.

### ▶ 2026-10-09 afternoon (3) — deeper wall defects and self-propulsion? (theory/selfprop/selfprop.py -> selfprop.txt, selfprop_events.npz)
- User's idea: elongated-cavity wall defects drift further from the wall -> stronger self-advection -> escapes. 5.4M wall +1/2 samples (675 runs).
- Depth: elongated 1.43-1.45 vs circle-like 1.35-1.38 (mean +0.08 cells); deep tail doubled (3-4 cells: 4.1% vs 2.0%). P(escape in 20 t.u.) rises steeply with depth
  (circle 0.0003 -> 0.073, elongated 0.0010 -> 0.125 from 0-1 to 3-4 cells); at equal depth elongated 2-4x higher.
- Orientation: inward head component cos(tau) LOWER in elongated (0.41 vs 0.51-0.55; 18% with head towards the wall vs 4%); P(escape) nearly flat in cos(tau). Opposite to the prediction.
- Self-propulsion signature: inward drift in 10 t.u. rises with cos(tau) in elongated (-0.02 -> +0.09), absent in circles (pinned); overall inward drift 0.04-0.07 vs 0.01-0.03.
- Reweighting circle-like P(esc | depth, cos tau) onto the elongated distribution explains 21% of the escape excess (upper bound: escaping defects count as deep).
  -> the mechanism operates but ~80% of the excess is weaker holding / stronger outward forcing at equal depth and orientation. Next: drift-vs-depth (effective restoring force) per class.

### ▶ 2026-10-09 ~12:50 — tilt excess: defect cores and wall density (theory/shear_tilt/core_S.py -> core_S.txt; runs/tiltdiag/); activity scan v1 failed
- User's question: can MPCD density fluctuations near the wall shift the tilt? Tested together with the core idea (alignment set where S is low).
- **Existing data (passive swirl, director every 100 t.u., ~15 frames):** cells within 1-1.5 of a wall-defect core have S 0.65-0.68 (bulk 0.88, layer away from cores 0.83);
  Maier-Saupe S4 -> pooled lam_eff 1.99-2.01 (lam 2), 3.96-3.99 (lam 4), 1.22 (lam 1.2) -> arccos(1/lam_eff) 60.0-60.2 / 75.5 / 34.8 deg vs measured 61-64 / 83 / 57.
  Closes most of the lam-2 gap (bulk-S value 56.5), little at lam 4 and 1.2. Caveat: Maier-Saupe S4 is wrong in core cells (mixed orientations) -> measure S4 from particles.
- **Runs launched (runs/tiltdiag/, bin/mpcd_swirl.out, 12:44, ~2.5-5 h):** passive k6 R31.83 lam 2 swirl 0.01; base_s1/s2_traj (densOut 10 t.u., director 2 t.u., trajOut = particle orientations every 50 t.u.
  -> S2, S4, density vs depth and core distance); d10_U10, d40_U10 (density per cell 10/40); d20_U5, d20_U20 (mfpot 5/20). Prediction if S-dependent alignment sets the excess:
  tilt rises with lower S (U 5, density 10), falls with higher S (U 20, density 40).
- **Activity + backflow Q-tensor scan v1 (scan_active.txt):** no active flow for zeta <= 0.03 (vrms ~1e-3, 12-ring survives) — activity range set from sqrt(K/zeta) ignoring the viscosity ratio;
  lam 4 diverged (NaN). Passive + backflow with body force 0.003/0.01: circulation 0.9, tilt 74-81 (coherence 0.31-0.55), no ejection (backflow raises the threshold?).
  v2 launched (scan_active2.py): zeta 0.1-1 (lam 2, 1.2), lam 4 at dt 0.005, stronger body force.

### ▶ 2026-10-09 ~12:40 — gamma1 runs launched; first-pass gamma1; strong-swirl ellipses done
- **gamma1 (rotational viscosity), first pass from the existing Frank boxes** (runs/Kcalib/rot_visc.py -> rot_visc.txt): decay of equilibrium director modes, gamma_eff = K q^2/Gamma_q:
  splay (q || y, director || x) ~ 40-70, bend ~ 25-70, mean ~ 50 (kBT tau/cell^2); noisy (each mode only 20-60 decay times). >> the 2-6 needed by the flat-wall tilt-correction ansatz:
  with gamma1 ~ 50 and shearSusc 0.5 (halves the flow torque) the ansatz predicts ~70-72 deg vs measured 61-64 -> **the correction prefactor does not transfer to MPCD**; leading term
  arccos(1/lam_eff) unaffected. Older results use the measured wall-defect mobility mu = 0.02, not gamma1 (stability onset, Langevin model, fluctuations): unaffected.
  The model-to-MPCD unit mapping (U = 1 <-> A ~ 0.015) assumed gamma1 ~ 20; with ~50 it shifts the other way (MPCD threshold maps to larger U) — recheck when gamma1 is final.
- **Runs launched** (runs/rotvisc/, production mpcd.out, run_queue.sh, 6 parallel, ~3.5 h): modes_L64_s1-4 (passive aligned 64^2 PBC, 5000 t.u., director every 2 t.u.) for gamma_splay/bend with
  ~4x the statistics; quench_L96_s1-2 (passive 96^2 PBC from random orientations, defects every 1 t.u., 1500 t.u.) for the FREE-defect mobility (pair annihilation) vs the wall value 0.02.
  Analysis: runs/Kcalib/rot_visc.py on the modes runs; quench: pair approach dr/dt vs r (to write).
- **Strong imposed swirl in ACTIVE ellipses (runs/swirl_ellipse_strong, A 0.06/0.1 = 6-10x before; overnight/swirl_ellipse_strong_summary.txt):** still NOT repaired: wall 7.2-7.5 (eps 1.5),
  7.3 (eps 2.0); bulk 1.7-1.95; tilt coherence 0.63-0.65 (eps 1.5), 0.35-0.37 (eps 2.0); circ/vrms 0.48-0.66. Activity keeps the coherence below circle values even at 10x forcing.
  Matched-coherence comparison: PASSIVE eps 1.5 with imposed swirl reached coherence 0.62-0.77 and wall 10 (+1 bulk); ACTIVE eps 1.5 at coherence 0.63-0.65 stays at 7.2-7.5 ->
  at similar coherence the passive layer refills and the active one does not: activity removes wall defects by more than lowering the flow coherence (escape kinetics).

### ▶ 2026-10-09 afternoon (2) — ANALYTIC TILT ANGLE: flat no-slip wall in shear (theory/shear_tilt/flatwall*.py/txt/json, mpcd_profile.py/txt, profile_dir.txt)
- **Result:** tau = arccos(1/lam_eff) + d sqrt(gamma1 s / K) (lam_eff^2 - 1)^(-1/4)   [rad for the correction; d = defect depth, s = wall shear rate].
  * Leading order (low Ericksen number): far from the wall the director aligns at the Leslie angle theta_L (cos 2theta_L = 1/lam); the screened layer is a zero mode
    (one-constant: registry shift <-> uniform rotation; the lock) so it follows rigidly; a +1/2 in a field rotated by theta_L has its axis rotated by 2 theta_L. Flat-wall
    Q-tensor model (periodic, one site spacing, simple shear, free far boundary): at s -> 0 tilt -> arccos(1/lam) (61.8 vs 60 at lam 2; 49.5/48.2 at 1.5; 76.9/75.5 at 4); far field
    2alpha = 2theta_L within 1 deg; registry x 360 = tilt (lock).
  * Correction (boundary layer): the anchored band of thickness d between wall and defects cannot align; its vorticity torque gamma1 s d / 2 is transmitted to the layer and
    over-rotates the director by d/(2 l_E sqrt(lam^2-1)), l_E = sqrt(K/(gamma1 s sqrt(lam^2-1))) (alignment length). Tested on 35 flat-wall runs (lam 1.1-4, s 0.001-0.1,
    a_s 10-40), measured defect depth, no free parameter: measured/predicted median 0.97, corr 0.98. Explains the growth near tumbling ((lam^2-1)^(-1/4)) and the weak sqrt(s)
    dependence on flow strength.
- **MPCD flow alignment:** per-particle Jeffery (larsonRotRate; shearSusc scales vorticity and strain equally -> angle unchanged). 2D moment: lam_eff = lam (1 + S4)/(2 S2);
  Maier-Saupe S4 at the measured S2 0.83-0.88 -> lam_eff = 1.80-1.83 (lam 2), 3.61-3.66 (lam 4), 1.08 (lam 1.2). Leading order: 56-57 deg (lam 2), 74 (lam 4), 23 (lam 1.2).
- **MPCD check (passive swirl, director profiles):** tilt = 2 alpha at the defect depth (1.75 cells): 65.5/61.6 vs measured 63.7/61.0 (lam 2); alpha in the layer 26-33 deg vs theta_L 28.
  Matching the measured 61-64 deg needs gamma1/K = 0.07-0.18 (gamma1 ~ 2-6 kBT tau/cell^2 at K 32.5) -> then lam 4: 77-79 (meas 83), lam 1.2: 32-37 (meas 57, coherence 0.54; near tumbling,
  linearisation fails). **Needed: an independent MPCD measurement of gamma1** (the earlier guess ~20 from the defect mobility would overshoot).
- Disk Q-model with the measured no-slip flow: bulk director at theta_L (29.1/36.7/16.5 vs 30/37.8/16.8); tilt = 2 x local Leslie angle with S-reduced alignment in the layer (46.0/64.7/77.0 vs 46.3/66.7/78.9).

### ▶ 2026-10-09 afternoon — tilt from the no-slip shear: continuum Q-tensor model (theory/shear_tilt/; shear_tilt.py, scan.py, scan.txt, shear_tilt.png)
- User's question: does the tilt theory include the shear from the no-slip wall? It did not (tilt theory purely elastic; Stokes solves had no-slip but no torque on the director). New model: passive 2D Beris–Edwards Q-tensor (one-constant LdG, core ~1 cell, S 0.876) in the anchored disk (k 6, R 31.83), imposed flow = measured passive-swirl v_θ(r) shape × U, one-way (no backflow), flow alignment λ_c = λ·S (so cos 2θ_L = 1/λ as in MPCD).
- **Reproduces the passive MPCD swirl results with no fitting:** 12 → 10 at the wall + bulk +1 above a flow threshold (λ 2: between U 0.2 and 0.3; coherence below it rises 0.07 → 0.56, then 1.00 — same pattern as swirl_threshold); tilt after ejection 65–67° at λ 2 for U 0.3–2.0 (independent of flow strength, MPCD 57–68°), 78–79° at λ 4 (MPCD 75–83°), 46° at λ 1.2 (needs U 1; MPCD 57–65°, 43° after the strong reversal); λ 0.5 never ejects (MPCD: never). Starting from the ejected 10-ring below threshold, the pair returns to the wall (12).
- **No-slip test (flow profile variants):** shear only in the wall band (rigid bulk) gives the same tilt (67.7 vs 66.7 at λ 2; 79 vs 78 at λ 4). No shear in the wall band (rigid rotation there, slip) changes the angle to 96° (λ 2) / 101° (λ 4) and raises the ejection threshold (λ 2, U 0.3: no ejection). → **The tilt angle is set by the no-slip shear in the wall band** through flow alignment (strain vs vorticity of simple shear); bulk shear helps reach the threshold but does not set the angle. Rough unit mapping (γ1 ~ 20 from μ 0.02): U = 1 ↔ MPCD swirl A ~ 0.015 (factor of a few) — threshold consistent with MPCD 0.003–0.01.
- Caveats: one-constant (no pinning; passive slides anyway), no backflow, LdG core size and λ mapping assumed. Next: add backflow / K1≠K3; active version (activity instead of imposed flow); the slope 2θ_L vs model at small λ.

### ▶ 2026-10-09 midday — do wall escapes cascade? NO in elongated cavities (theory/avalanche/avalanche.{py,txt}, vacancy_binding.{py,txt})
- User's question: does an escape (at the tips) destabilise the layer and help the next escape? 46k escapes, 675 runs.
- (1) Temporal: elongated escapes are Poisson-independent (g 0.85–0.95 at lag 10–60 t.u., 1.00 beyond: slight refractoriness). Circle-like escapes do cluster (g 1.1–1.7 at lag < 100, strongest at R_eq 1.0–1.5) — consistent with the two-at-a-time flow-induced ejection.
- (2) Spatial (pairs < 100 t.u. vs 500–2000 t.u.): same-site re-escape depleted (0.65–0.70, site empty until refilled); neighbour excess 1.2–1.5 in circle-like, only 1.04–1.09 in elongated.
- (3) Hazard per wall defect vs current vacancies, within-run: elongated flat or falling (0.89–0.98 with more vacancies; R_eq 1.0: 2.05 → 0.74) = restoring, not runaway; circle-like rising (1.2–1.8) — but common cause (incoherent-flow episodes raise both) not separated. Pooled curves rise steeply in both: cross-run confound. (3b) A tip escape does not raise non-tip escapes in the next 100 t.u. (0.98, 27.5k obs).
- (4) Frank model (anchored disk, k 6, 10-ring + bulk +1): a vacancy makes its two neighbours MORE bound (+8–10 kBT on 69–77 for a pull of 6 cells) and farther sites less (−5…−9), whatever the escaped +½'s position; one-constant and anisotropic agree. No elastic local cascade.
- → The elongated deficit is a higher, steady, independent escape rate per site (set by the flow/shape state), not an avalanche. Mild collective escape exists in circles, where it fits pair ejection by coherent flow.

### ▶ 2026-10-09 late morning — escape/refill kinetics, bulk carriers, coherence conditioning (theory/shapes/kinetics_carriers*.{py,csv,txt}, interval_coupling_conditioned.txt)
- **Kinetics (675 runs, second half, wall +½ / bulk +½ followed 20 t.u.):** escape per wall defect is 3–7× higher in elongated cavities (aspect ≥ 1.25) at every R_eq (g_esc 2.5–3.2 vs 0.4–1.0 per 1000 t.u.). At onset escape is tip-dominated (ε 1.2: tip 0.014 vs flank 0.003 per 20 t.u.; ε 1.4: 0.041 vs 0.012); by ε ≥ 2 the flanks catch up (0.033/0.029). A bulk +½ reaches the wall band equally often in all shapes (p_ref per 20 t.u. ratio elongated/circle 0.86–1.07). Arrivals ≫ escapes (≈30 vs 3 per 1000 t.u. in circles): most arrivals are transient visits. Net incorporation per vacancy (= escapes/vacancy in steady state) is ~0.5–0.6× the circle value. Across runs vacancies follow escape (corr 0.88 circle-like, 0.62 elongated). → **The deficit is mainly escape-driven**; refill is not starved of bulk +½, but a vacancy captures an arriving +½ less efficiently.
- **Bulk carriers:** the excess bulk charge (net 1.0–2.2 vs 0.7–1.1; corr 0.99 with vacancies, by charge conservation) is carried by **free +½** (only 1–5% of bulk +½ have a −½ within 4 cells; fewer in elongated). Elongated cavities have slightly more bulk +½ (×1.14) and markedly **fewer bulk −½** (0.11/0.86/2.39 vs 0.44/1.93/3.91 at R_eq 1/1.5/2): less bulk pair creation. The unpaired +½ are **not trapped at the tips**: mean |x|/A 0.39–0.42 in every shape (uniform in an ellipse ≈ 0.42), same depth profile as in circles.
- **Coherence conditioning (interval_coupling_conditioned.txt):** pooled over sizes, elongated vacancies fall from 2.74 (coh < 0.5) to 1.15 (coh > 0.95), but this is **mostly a size effect**. At fixed R_eq, coherent elongated intervals (coh > 0.85) still have 0.54 / 1.85 / 3.68 vacancies at R_eq 1 / 1.5 / 2 (circle-like: −0.01 / 0.12 / 0.54). Within-run (vacancies minus the run mean): elongated −0.10 at coh > 0.85 vs +0.04 at coh < 0.5; circle-like +0.73 at coh < 0.5 vs −0.03 at coh > 0.85. → In elongated cavities, **transient** circle-like coherence (100 t.u. intervals, 2–14% of the time) does not refill the layer, while in circles the wall count follows coherence strongly. Not excluded: refill needs **sustained** coherence (circle refill lags coherence by ~100 t.u.) — that is what runs/swirl_ellipse_strong tests.

### ▶ 2026-10-09 morning — overnight results (reports/master/theory/overnight/*_summary.txt, theory/shapes/interval_coupling.txt)
- **Causal test, ellipses with an imposed single vortex (swirl_ellipse):** PASSIVE eps 1.5: wall 12 → 10, bulk +1 (0.74–0.97), tilt coherence 0.62–0.77, one vortex → coherent flow gives the circle structure in an ellipse (curvature alone does not prevent a 2k−2 layer); eps 2.0: 12 → 10–11, coherence lower (0.23–0.44), two vortices return after the reversal. ACTIVE eps 1.5/2.0: imposed vortex does NOT repair the layer (wall 7–8, bulk 1.5–1.7, P(two) 0.36–0.56) even at A 0.03 with circ/v_rms 0.57–0.66.
- **Interval coupling (450 elongated runs, 22,500 intervals):** wall count independent of the instantaneous flow (|corr| ≤ 0.05 with circulation, ≤ 0.15 with tilt coherence at any lag; vacancies 2.43 one-vortex vs 2.13 two-vortex). Control (180 circle-like runs): wall count follows the flow (corr +0.25 with circulation at lag 0, +0.24 with coherence at lag +1; two-vortex intervals z −0.84).
- → NOT decisive (corrected after the user's objection, 2026-10-09): the imposed vortex never made the active ellipses circle-like coherent (tilt coherence 0.62–0.64 vs 0.97 in circles; two vortices still 36–40%), so activity disrupting the coherent flow remains a viable explanation of the deficit. The interval data argue against it only weakly (elongated intervals may never reach circle-like coherence). Sharper tests running: runs/swirl_ellipse_strong (active eps 1.5/2.0, constant swirl 0.06 and 0.1) and theory/shapes/interval_coupling_cond.py (vacancies binned by tilt coherence).
- **Ejection threshold (swirl_threshold, passive circle, constant swirl, 2 seeds):** no ejection at A ≤ 0.003 (tilt coherence 0.2–0.5); A 0.01: both seeds 12 → 10, bulk +0.75–0.93, coherence 0.90–0.96. Threshold 0.003–0.01 for a constant swirl; the earlier A 0.003 ejection had a reversal (transient helps).
- **Negative layer (swirl_neg, passive k −6):** wall stays 12 −½ for λ 0.5, 2 (both A); λ 4: 12 → 13 (bulk +0.5). Energetics favour gaining two but pair creation is needed and is rare without activity → kinetically blocked, not refuted. Fair test: active, or a seeded bulk pair.

### ▶ Night of 2026-10-08 → 09 (user asleep) — jobs and plan
- Running (office, bin/mpcd_swirl.out): runs/swirl_neg (4, passive k −6 imposed swirl: does the −½ layer gain two?), runs/swirl_ellipse (8: CAUSAL TEST — imposed single vortex in passive/active ellipses eps 1.5/2.0; does the layer return to 2k−2 with a common tilt?), theory/shapes/interval_coupling.py (detached; per-interval wall count vs flow coherence).
- Queued: runs/swirl_threshold (10 passive circle runs, constant swirl 0.0005–0.01, 2 seeds; ejection threshold), queue.sh starts after swirl_neg.
- Detached analysis driver (survives a session crash): reports/master/theory/overnight/driver.sh → swirl_neg_summary.txt, swirl_ellipse_summary.txt, swirl_threshold_summary.txt (counts.py: wall ±½, bulk charge, circulation/v_rms, P(two), tilt coherence per window); log driver.log. interval_coupling → theory/shapes/interval_coupling.txt.
- Then (needs the session): add shapes + causal test + ejection statistics to the master PDF, restructure §28 into the paper storyline.
- Phase diagram (theory/shapes/phase_diagram.py, phase_shapes.png): coherence/circulation/P(two) collapse vs aspect alone, transition aspect 1.25–1.37 for ellipses and ovals at all R_eq (A/ℓ_a 1.2–2.3: no length sets it); vacancies above it ≈ 1/3/4.5 at R_eq 1.0/1.5/2.0. The 2026-10-08 22:30 session crash (user's connection) killed interval_coupling (exit 137, no OOM); restarted detached.

### 2026-10-08 night — shape sweep analysed with the circle framework (theory/shapes/; data copied to runs/shapes_server, 675 runs)
- shape_sweep_analysis.py → shape_runs.csv, shape_escape_curvature.csv. Medians (all R_eq, k, rho): circle vac 0.04, Q_b 0.77, coh 0.94, circ 0.68, P(two) 0, esc 0.0008; ellipse 1.2 / ovals 1.12–1.22 circle-like (vac ≤ 0.17, coh 0.92–0.94, circ 0.60–0.66, P(two) 0); oval 1.29 intermediate (vac 0.95, coh 0.78, circ 0.39, P(two) 0.12); aspect ≥ 1.4 broken (vac 1.5–2.9, Q_b 1.4–2.0, coh 0.39–0.63, circ 0.04–0.30, P(two) 0.32–0.84, esc 0.012–0.020). Ellipses ≈ ovals at matched aspect → transition between aspect 1.2 and 1.4. Lock holds everywhere (0.91–0.97).
- Size: aspect ≥ 1.6 vacancies ≈ 1 / 3 / 4.4 and Q_b 1.1 / 1.8 / 2.5 at R_eq 1.0 / 1.5 / 2.0 ℓ_a (circles 0, 0.7–0.85).
- Escapes vs local curvature (pooled): rise at κR ≳ 1.25 (ellipse 0.019–0.032, oval 0.013–0.018) but non-circular shapes leak 5–10× more than circles even at κR ≈ 1 → global + local.
- Lock-produced tilt binding (lock_binding.py, natural cores, reliable): one vortex weakens tips 1–2 kBT vs flanks (eps 1.5: 16.4/17.4; eps 2: 14.7/16.6); two vortices weaken the flanks (meeting points) by 3–4 kBT (13.9, 12.2) — matches the eps 1.5 flank weakening in two-vortex intervals; does not explain weak tips in two-vortex intervals; bulk-pair placement near the tips strengthens tip binding (confound).
- Working hypothesis: elongation → loss of coherent circulation (coherence 0.94 → ~0.45) → no common tilt imposed → the layer is no longer held at 2k−2 and loses more defects with size. To test: phase diagram (aspect, R_eq); per-interval wall count vs flow coherence; director-field Stokes on the shapes.

### 2026-10-08 night — per-site misfit test (theory/ellipse_binding/site_misfit.py): the misfit mechanism is ruled out
- Local lock holds at every site, tips included: tail tilt − 2ψ_rel(s) concentration 0.92–0.99, offset ≈ 0 (ψ_rel = k 2πs/P − φ_n(s)). No forced core rotation: defects take position-set tilts (eps ≥ 1.5: tips 70–71°, middles 83–97°, flanks 68–76°). Deviation from a common 63°: rms 17–24° in every region; its forced-rotation cost 2–6 kBT, larger at the flanks → predicts no tip weakness → curvature-misfit mechanism ruled out.
- What differs at the tips: displacement from the nearest homeotropic point 2.8–3.1 cells (flanks 1.75–2.2), depth 2.1–2.2 (flanks 1.75–1.95).
- Noise check: fluid speed 1–5 cells from the wall not larger at tips (0.31/0.47, 0.23/0.26, 0.25/0.23 tips/flanks at eps 1.2/1.5/2.0); wall-defect rms step 0.83/0.82, 0.93/0.86, 0.96/0.95 → shallower tip wells are weaker binding, not stronger noise.
- Running: lock_binding.py — binding with lock-produced tilts (sites on the real homeotropic points of the ellipse, displaced to 63° by the lock, natural cores; reliable mode).

### 2026-10-08 night — measured binding in ellipses; analytic estimates (theory/ellipse_binding/measured_binding.py)
- Measured effective potential of wall +½ (U = −ln P(depth), units of the active effective temperature; W3 ellipses 2 seeds, circle 3 seeds): well depth U(6) − U_min: circle 3.08; eps 1.2 tips 3.26 / flanks 3.02; eps 1.5 tips 1.97–2.19 / flanks 2.64–3.20; eps 2.0 tips 1.90–2.07 / flanks 2.38–2.41. Tips shallower by ~0.5–1 → e^0.5–e^1 = 1.6–2.7× escape, matching the measured 2–2.5×; tip defects sit deeper (mode 1.75 vs 1.25). Vortex state minor (tips equal in both; eps 1.5 flanks weaker with two vortices). Static elasticity predicted only ~8% weaker tips. Caveat: equal active noise at tips and flanks assumed, not checked.
- Analytic (answer to the user): (1) one-constant energetics exact (images / conformal map); (2) forced core rotation c costs (K/2) C c², C = 2π/arccosh(d/a) (disk above a grounded line): 12-ring common tilt ≈ 120 kBT (a 1.5) – 400 (a 2.5) vs relaxed 350–480 → same order, explains the core-size dependence; (3) ellipse lock slope d tilt/ds = 2(k 2π/P − κ(s)): tip/flank 1.00 (circle), 0.75 (eps 1.5), 0.52 (eps 2) → a flow-imposed common tilt cannot be reached by a uniform registry shift; tip defects are left under-tilted → residual forced rotation → weaker binding at the tips, for one or two vortices, growing with eps — matches the measurement. Next: per-site residual rotation for a best-fit registry shift at 63°, binding loss from the capacitance formula, vs the measured tip/flank well depths.

### 2026-10-08 night — negative layer and ellipse binding (theory/ejection/ejection_neg.py, theory/ellipse_binding/binding.py)
- Negative winding k −6 (R 28.21): passive 12 −½ at the wall; active 13–14 (+ bulk ~0.6). Energetics: the 14-ring (2|k−1|, commensurate) has uniform natural tilts (±60° ≡ one arm into the bulk), the 12-ring's wind (50, 30, 10, …). Without flow, gaining two −½ costs +64…+181 kBT + 4 core energies (created ±½ pair). Under a common tilt the 12-ring costs +350…+480; the 14-ring tilts by registry (lock slope −14/3) almost free → 14-ring + bulk +1 cheaper by 196…392 kBT − 4 E_core. (Bugs fixed on the way: B tilted via cores instead of registry; numerical slope broken at the ±60° wrap.) Direct test running: runs/swirl_neg (passive k −6, imposed swirl, λ 2 A 0.01/0.02, λ 4 A 0.01, λ 0.5 A 0.01).
- Ellipse binding (user: which wall defect is weakest with one/two-vortex tilts?). Field = direct defect field + correction fixed by the eMode-2 anchoring on the wall; rigid cores; anisotropic; 10 wall +½ at depth 1.6 + bulk pair ±0.45A; binding = energy to move a defect 1.5 cells inward; h 0.35, mean over 4 grid offsets. Natural (untilted): 11–13 kBT ± 0.2–1, tips only ~1 kBT weaker than flanks (eps 1.5: 11.0 vs 12.1; eps 2: 11.2 vs 11.6) → static elasticity does not explain the 2–2.5× tip escapes. Imposed tilts (rigid core rotation, smallest rotation from natural): circle check FAILS (equivalent sites sd 3.5–4 kBT), binding drops to ~0 ± 5 everywhere, contradicting the stable tilted active circles → method not trustworthy for tilted states (first attempts: π-branch error, then ±8 kBT noise). Needs the tilt to come from a flow torque in the relaxation, or a measurement of binding from the simulations.

### ★ 2026-10-08 night — FLOW-INDUCED EJECTION: why the active wall holds 2k−2 (user's mechanism; master §29.14, §28)
- Mechanism (user): vortex → flow → flow tilts the wall defects → tilted 2k ring too expensive → two leave.
- Energetics (theory/ejection/: tilted_structures.py, ejection_energetics.py; rigid cores, relaxed; k 4/6/8; one-constant and anisotropic; core 1.5/2.5): ejection without flow +155…+200 kBT (validation: unrelaxed one-constant +183.7 vs image formula +184.0); forcing a common tilt on the 2k ring (natural tilts wind: 150, 90, 30… at k 6) costs +200…+500 kBT, ~independent of τ; the 2k−2 ring tilts by registry almost free (< 40); under a common tilt E_B − E_A = −20…−330 kBT (sign robust).
- Direct test: PASSIVE imposed-swirl runs eject: λ 2 A 0.003 12→11→10 (bulk 0→0.5→0.8), A 0.01/0.02 12→10 (bulk 0.9–1.0); λ 4 both A 12→10 (bulk 1.0); λ 1.2 partial (A 0.01: 12→11, bulk 0.5); λ 0.5 (tumbling) never. The ±57–83° 'passive tilts' were measured on the 10-ring after ejection.
- Reading: 2k−2 + bulk +1 follows from flow + flow alignment; activity supplies the flow. Predictions: bulk-charge onset tied to coherent flow onset; smaller λ keeps more defects at the wall under the same flow.

### 2026-10-08 evening — flow-alignment test and site-arrangement sweep (master §29.10, §29.11)
- Tilt vs tumble λ (runs/swirl_tumble, passive imposed swirl; swirl_tumble.txt): λ 0.5 (tumbling) ≈ 28–36°; 1.2: 57–65° (43° after the flip at A 0.01); 2: 57–68°; 4: 75–83°. Monotonic in λ → the angle is set by flow alignment at the wall; 2θ_L = arccos(1/λ) fails at 1.2 (34°) and the tumbling case still tilts → agreement at 2, 4 coincidental.
- Site-arrangement sweep (runs/anchmod_sweep, 14 runs, layer_observables.txt): m 1 A 1/2/3, m 2 A 1/1.5/2 → bulk charge 0.75–0.87 (baseline 0.76–0.83), wall 10 in 85–97% (88–93%), tilt 61–63° coh 0.96–0.97 (63–64°, 0.97), no reversals, no escapes (one run 0.05/1000), circ 0.59–0.75 (0.62–0.75), never two vortices. → layer robust to site arrangement; ellipse instability from curvature/elongated flow (circles never form a second vortex, so this sweep alone doesn't separate curvature from the second vortex).
- **Half-integer k (runs/halfk, analyse_halfk.txt; master §29.13):** active k 2.5/3.5/5.5 = integer behaviour: wall 3/5/9 (= 2k−2, odd) in 85–93% of frames, registry on the odd homeotropic set (lock D 0.00 ± 0.09, conc 0.99–1.00), tilt 63–65° coh 1.00, single vortex 0.66–0.74, no reversals. k 1.5 (1 site): no layer (wall 0.04), ~3 bulk +½ (Q ≈ 1.25), weak erratic flow 0.08–0.12, 6–7 reversals → layer needs ≥ 2 sites (smoothing depth ≈ 2R for N = 1). k 0.5: single bulk +½ at ~0.6R, strong steady vortex 0.42–0.49. Passive: depth law (k 0.5 near centre; 1.5–3.5 deep rings straddle the 4-cell band; 5.5 all 11 at wall, frustrated tilts).

### 2026-10-08 — imposed swirl results (swirl_response.txt; master §29.11)
- Active: A 0.005 does not reverse the vortex (circ 0.79 → 0.40, tilt stays +64): bistable. A 0.01 / 0.02: circulation zero crossing 375 / 295 t.u. after the flip, tilt 678 / 384 (lag ~90 at the stronger drive, as natural reversals); registry jump −7.2 / −7.3 cells (lock −6.8 / −7.0), transit 344 / 42 t.u. → stick–slip: pinned until the flow at the layer reverses, then jumps to the mirror registry.
- Passive: imposed flow tilts the layer to ±57–68°, nearly independent of the flow (0.02–0.09 at the layer), and the layer slides with the flow (40–44 cells; not pinned). → wording (user, 2026-10-08): in the cavities the tilt IS an active effect (the only flow is active flow), but its angle appears set by the layer's response to the presence of flow. Candidate: 2θ_L with θ_L = ½ arccos(1/λ) = 30° (λ = tumble = 2). Test: runs/swirl_tumble (passive, tumble 1.2 → ≈ 34°, 4 → ≈ 76°, 0.5 tumbling → no locked tilt).

### 2026-10-08 — wall-layer response to flow: reversals, build-up, imposed swirl
- Reversals (act 0.05, H3pilot k 9/10 + seeds; theory/reversals/registry_jumps.py): layer weak at this activity (tilt ±10–40°, vacancies, many noise reversals); jump timing/ratio unreliable; but the registry shift at each reversal follows the tilt lock event by event (e.g. +1.47 vs +1.42, +1.38 vs +1.41, −0.56 vs −0.54 cells). Act 0.07 never reverses.
- Build-up at act 0.07: production output starts after the warm-up (tilt already 61–65° for k 4/6); k 10 catches it: tilt 9° → 59° over ~800 t.u., layer moves ~1–1.4 cells at ~0.002 cells/t.u. vs fluid 0.03–0.17 → β ≈ 0.01–0.03. Dedicated no-warm-up runs (runs/tilt_buildup, buildup.py): dominated by layer formation (tilt swings ±140°, full-site phase slips) → inconclusive.
- **Imposed swirl (user's suggestion):** branch anchmod + body force: inputs swirlAcc, swirlRadius, swirlX/Y, swirlFlipTime, swirlPeriod (production steps only; mpc.c after acc_all). Binary bin/mpcd_swirl.out (mpcd_anchmod.out untouched for the running/queued jobs). Validation (runs/swirl_validation): swirl off = production bit for bit; on (A 1e-3, 0.6R, flip at 120): circulation +0.03 at 120, crosses zero ~70 t.u. after the flip. Runs: runs/swirl_test (k 6 R 31.83 eMode 2; passive A 3e-3/1e-2/2e-2, active 5e-3/1e-2/2e-2; within 0.6R; flip at production 1500; flow every 10 t.u.). Analysis: theory/reversals/swirl_response.py.

### 2026-10-08 — is the second vortex the driver of the layer instability in ellipses? (user's hypothesis; theory/vortexlayer/)
- Tracked escapes (wall +½ followed 20 t.u. by nearest-neighbour chains, escaped if deeper than 6 cells; vortex_vs_layer_tracked.txt; first version with single-frame escapes had too few events, discarded): eps 1.2 (always one vortex) 0.005 per wall defect per 20 t.u.; eps 1.5 one vortex 0.072 / two 0.066; eps 2.0 0.070 / 0.077 (n = 1500–2700 events per state). Vacancies 2.31/2.31, 2.70/2.73. Escapes larger at the tips (|x| > A/2: 0.08–0.10) than at the flanks where two vortices meet (0.04–0.05, eps 1.5). Lagged corr(two-vortex, escapes) |r| ≤ 0.11 at −200…+200 t.u.
- → No instantaneous link: the layer is ~14× leakier at eps ≥ 1.5 in both vortex states, escapes concentrate at the high-curvature tips. Elongation/tip curvature looks like the common cause of both. Caveat: slow healing could hide a time-integrated effect; the circle site-arrangement sweep (no tips) separates it.

### 2026-10-08 — annular-drive model: closure fails; wall defects are not advected
- Drive = f_eff(r) from the measured mean swirl (k 6, R 31.83; ≈ c⟨−∇·Q⟩_φ in the ring r/R 0.7–0.85, weaker wall brake), × sin τ / sin 63°, no Stokeslets (annular_profile_k06.json, scan_annular.txt): tilt does not lock at 63° (median 26–29°, broad, slips ≤ 130°/1000 t.u., 20–40 reversals); more drive → less tilt; half drive 37–55°.
- Measured (Rscan_R31.83_seeds/s02, k 6/10): fluid at layer depth 0.13–0.16 cells/t.u. (signed), wall +½ tangential velocity 0.0005 ± 0.0010 / 0.0004 ± 0.0005 → advection coefficient 0.002–0.004 ≈ 0. Tilted layer = static shifted registry, held by the anchored texture; the model's free advection makes it slide (pinning holds ≤ μF ≈ 0.03).
- Next: measure the wall-defect advection/mobility from the registry jumps at circulation reversals; then torque balance with a static drag.

### 2026-10-08 — flow alignment (theory/flowalign/)
- MPCD: Jeffery rotation (Larson) λ = tumble = 2, scaled by shearSusc 0.5 → flow-aligning, Leslie angle 30°. Not in any model.
- Residual director (measured − one-constant field of measured defects), signed by circulation: −1° at wall → −25° at 6–10 cells (k 6, 10; −10° k 2); grows with local shear (5° → 10–12°), sign follows shear. At wall-defect depth only 4–5° → no direct extra tilt.
- Bulk pitch: elastic 122–123° (= 90 + τ/2) vs measured 110° (2–6 cells from wall), 74° (k 6) / 87° (k 10) in 0.3–0.7R → flow alignment erases the elastic spiral in the bulk.
- Measured azimuthal active force (−∇·Q, signed): net drive in a ring r/R 0.6–0.9 (just behind the layer, ~ smoothing zone), brake r/R > 0.9, slightly negative 0.2–0.5R. → effective-boundary drive confined to an annulus; particle model's 1/r spiral force overestimates and misplaces the drive (fits locking needing ~half drive). Next model: measured annular drive + flow alignment of the bulk texture.

### 2026-10-08 — tilt size: splay ≠ bend beyond first order (theory/aniso2/)
- Particle-model tilt not locked (wanders, layer slips 2.4 sites/1500 t.u., many reversals) vs MPCD ±63° locked → need a non-cos τ term.
- First order in K1 − K3 is exactly ∝ cos τ (rigid τ/2 director shift): −66/−44/−22 kBT (k 4/6/8), higher harmonics 0. Bulk +1 adds ≈ −46 cos τ.
- Relaxed anisotropic Frank energy (relax_frank2.py: correction field ψ = 0 on the wall, rigid cores, energy outside cores; K1 = K3 check flat to 0.8 kBT cos 2τ; first soft-core version failed its check, discarded): k 6 R 31.83: E ≈ −(31–47) cos τ + (4.6–7.6) cos 2τ − (1–3) cos 3τ kBT, cos 2τ sign robust (core 1.5/2.5, h 0.5/0.35, bulk point/pair, cores fixed/free).
- Torque sin τ [D − |A1| + 4A2 cos τ] → supercritical onset, cos τ* = (|A1| − D)/(4A2), sliding above |A1| + 4A2.
- Particle model + A2 (scan_a2.txt): A2 20 → locked 64° (R/ℓ_a 1.0, half drive; no slip/flips), 67° (0.6, full drive); A2 10 → ~49°; A2 5 (computed) → ~35°. Full calibrated drive at R/ℓ_a ≥ 0.8 slides (missing wall-band brake).

### 2026-10-08 — defect model + spiral body force (user's request); half-integer k queued
- defect_langevin.py: spiral flow from the effective boundary (c = 90° + τ/2, force (ζ/η)S sin τ / r on r_c < r < R − 0.34ℓ, (ζ/η)S = 0.48, only when the wall holds 2k−2). Free model (k 6, R/ℓ_a 1.43): 2k−2 + bulk state present 27–30% at v0 0.15, circulation 0.29–0.38 then (0.007 without), but unstable (avalanches), no consistent tilt.
- Structured (2k−2 held radially, anisotropy pinning B 90/130 kBT, bulk pair kept 6 cells from wall; scan_structured.txt): k 6 onset between R/ℓ_a 0.6 and 0.8; circ 0.40/0.65/1.0/1.6 at 0.8/1.0/1.43/2.0, tilt·sign(circ) +24…+29°; k 4 onset 0.8–1.0. No spiral: circ ≤ 0.014, no tilt rule. vs MPCD: circulation right size, tilt smaller (30–40° vs 63°), onset later (0.7 vs 0.4), k-ordering as predicted. Missing: what keeps the +1 in the bulk (passive ground state = 2k at wall).
- Half-integer k (odd homeotropic points; code anchoring is sign-insensitive, production binary): runs/halfk, k 0.5/1.5/2.5/3.5/5.5, R 31.83, 2 active + 1 passive each; queue.sh starts them after the anchmod test.

### 2026-10-08 — symmetry argument for the fixed sites (discussed with the user; master §28)
- The anchoring θ_wall = kφ is (k−1)φ from the normal → homeotropic at φ_n = 2πn/(2k−2), tangential halfway; rotating the cavity by β gives kφ − (k−1)β → only a (2k−2)-fold symmetry (k = 1 is the only fully symmetric case).
- Active sites = homeotropic points (tails = radial splay side at the wall); tilt τ displaces them by τ/(2k−2).
- One-constant elasticity has a hidden internal symmetry (uniform director rotation is free) → rotate defects by β + add (k−1)β to the director = zero mode. K1 ≠ K3 breaks it → homeotropic points become minima.
- Commensurability: active 2k−2 matches the pattern (80–110 kBT, pinned from any IC, only chirality spontaneous); passive 2k incommensurate (< 1 kBT, frozen SSB, creeps).
- Test: modulated winding (next entry) — sites should follow the uneven homeotropic points; passive should become pinned (modulation breaks the hidden symmetry even with one constant).
- **Why 2k passive defects cannot be untilted:** lock tilt_j = 2ψ_rel(φ_j) = 2(k−1)φ_j + const (mod 2π) → zero tilt only on the 2k−2 homeotropic points; 2k defects cannot all sit there. Evenly spaced N: tilt step 2(k−1)·2π/N, zero only if N | 2(k−1); N = 2k → step −2π/k (k 6: −60°, measured +41, −20, −91, …), tilts wind −4π around the ring → frustrated, registry energy averages out. Interpretation: ejecting exactly two +½ leaves the largest commensurate ring (2k−2), which can be fully untilted → the +1 bulk also relieves the wall's tilt frustration.

### 2026-10-08 — modulated winding test (user's idea), branch `anchmod`
- Code: niamh-overcharged branch `anchmod` (off main, uncommitted working-tree edits): eMode 2 anchoring angle k·2πs/P + anchAmp·sin(anchMode·2πs/P + anchPhase) (SRDclss.h AAMP/AMOD/APH, read.c, lc.c, InputGuide). Binary `bin/mpcd_anchmod.out` (clean build). Validation (runs/anchmod_validation/v3): anchAmp 0 = production bit for bit; A = 1, m = 1 wall director follows 6φ + sin φ (9° rms vs 41° from 6φ). Not in checkpoints. Production mpcd/ untouched.
- **PRELIMINARY (2026-10-07 23:35, ~350 t.u. of production; runs still going, machine loaded by PhaseField jobs; analyse_anchmod_prelim.py → anchmod_results_PRELIMINARY.txt):** active (9 runs): wall sites within 0.22–0.34 cells rms of the predicted modulated homeotropic points shifted by the tilt (2.2–4.5 cells from the best even lattice); gaps 17–25 (A 1 m 1), 14.5–30 (A 1 m 2); same-chirality runs agree within 0.14–0.23 cells across ICs, opposite chirality differ by 7.5 cells (= ±τ/(2k−2) each). Passive (3): not yet IC-independent (2.8–7.8 cells apart, one near-coincident pair) — judge at the end.
- Correction: on 2026-10-07 ~23:10 I told the user the anchmod runs were done and halfk had started; wrong (bad pgrep pattern). halfk/queue.sh still waiting.
- **FINAL (anchmod_results.txt):** active 9 runs: rms to predicted sites 0.21–0.27 cells (even lattice 2.2–4.5); same-chirality runs agree 0.06–0.13 cells, opposite 7.6 cells; tilt 63–64°. Passive 3 runs: NOT pinned (site sets 3.3–8.5 cells apart). My prediction was wrong: A sin(mφ) couples via Σ e^{imφ_j} = 0 for 12 evenly spaced defects unless m is a multiple of 12; m = 1, 2 only distort the gaps (13.9–20.3 vs 16.7) and leave the rotation free.
- **Geometry vs topology inside the circle (user, 2026-10-07):** moving the pinned sites by the modulated winding at fixed curvature and k separates site arrangement from curvature (cf. ellipse escapes 10–30× faster). Sweep queued with flow output: runs/anchmod_sweep (k 6, R 31.83, m 1 A 0/1/2/3, m 2 A 1/1.5/2, 2 seeds, random start; queue.sh starts after anchmod_test + halfk). Scoring: reports/master/theory/anchmod/layer_observables.py (bulk charge, wall count, |tilt|, coherence, reversals, escape rate, circulation, two-vortex fraction). Baseline A = 0 (8 runs): Q_bulk 0.77–0.83 (P(Q=1) 0.61–0.68; deficit = bulk defects passing the 4–6 cell band), wall 10 in 90–93% of frames, |tilt| 63–64°, coherence 0.97, no reversals/escapes in the second half, |circ| 0.62–0.75, never two vortices.
- Runs (office, started 2026-10-07 ~22:00): runs/anchmod_test — active k 6 R 31.83 A 0.5 m 1, A 1 m 1, A 1 m 2 (3 ICs each), passive A 1 m 1 (3 ICs). Analysis: reports/master/theory/anchmod/analyse_anchmod.py (prediction: sites where (k−1)φ + A sin mφ = τ/2 mod π; passive should become IC-independent).

### 2026-10-08 — passive registry = frozen spontaneous symmetry breaking (user's hypothesis confirmed)
- All production runs start with oDist = 2 (director along y). 8 passive k 6 R 31.83 runs with other initial directors (runs/pinning_ic): arg Z_12 = −2.62 (x), −1.44 (45°), −0.07/−0.23/+0.69 (y), +0.01/+1.20/+0.48 (random) ≈ 2(β0 − 90°); creep 0.3–0.5 rad / 3000 t.u.; Rayleigh p ≈ 0.1. The earlier cross-radius passive agreement was the shared IC. rand_s2 trapped an extra ±½ pair.
- **Active registry IS pinned** (runs/pinning_ic_active, ic_test_active.py): 6 runs (random ×3, x, 45°, y) + 2 production seeds all end at arg Z_10 = ±1.17 rad (= ±tilt 63–64°, one per circulation sense; concentration 1.00, no creep); lock-corrected D = +0.03 ± 0.07 in all 8 (Rbar 1.00, p ≈ 3e-4). Site fixed by physics; only chirality is spontaneous.

### ▶ THEORY STATUS (2026-10-08) — read master §28 "Where the theory stands" first
- Defect Langevin model (theory/defectmodel/: unit-disk Green's function, slaved orientations, 2k +1/2): no stable 2k−2 + bulk state (all stay below v0 ≈ 0.1, most leave above 0.15), circulation ≤ 0.02 → needs orientation-dependent wall binding + spiral bulk stress. Negative result, documented §29c.
- Effective-boundary stability (theory/stability/effbc_stability.py): loop registry → spiral boundary (c = 90° + τ/2) → bulk circulation → layer advection, vs anisotropy pinning (rate ∝ 1/r_s²). Calibrated: onset R/ℓ_a ≈ 0.47–0.50 (k 6), 0.4 (k 8), 0.6–0.75 (k 4) for r_c = 0.25 ℓ_a (0.30–0.63 for r_c 0–0.4 ℓ_a); measured 0.39–0.44. Predicted ordering (onset ↓ with k) untested: main sweep has onsets only for k 2–4 (0.46, 0.39, 0.40); defect-only tilt coherence too noisy at small N → need circulation of small high-k cavities (ρ 0.16–0.5; reduce flows on the server).

### 2026-10-07 night — director-field Stokes, effective planar boundary (master §29c)
- **Director-field Stokes (theory/dirstokes/):** f = −∇·Q from measured director (interval-averaged), wall band (<4 cells) vs bulk, FEM no-slip. One positive coefficient: R² 0.62–0.72 circles, 0.48–0.65 ellipses; circulation corr 0.61–0.85; c = 0.42–0.67 (vs ζ_eff/η ≈ 0.24). **Bulk drives** (alone R² 0.47–0.73, c > 0); **wall band brakes** (circulation corr −0.56…−0.85, share 9–17%). Ellipses 1.5/2.0: bulk reproduces the two-vortex part (corr 0.87–0.91), vortex number right 72–82%. → the two vortices come from the bulk active stress.
- **Effective boundary behind the layer (theory/effbc/):** with lc.c anchoring θ = kφ and the measured registry, the bulk +1 sees θ = φ + c, c = 90° (planar) untilted, 90° + τ/2 tilted (121.5° at 63°); patchiness < 0.1 rad beyond 0.25ℓ behind the ring (0.34ℓ from the wall). Spiral at 121.5° drives the same circulation sense as the tilt → elastic feedback loop.
- **Spacing criterion vs data (151 points):** bulk-charge existence set by R/ℓ_a ≥ 0.76 (99%, AUC 0.997), not by spacing (AUC 0.59); excess charge at R/ℓ_a ≳ 2.2 also size-set; spacing matters for pinning only on the crowded side (ℓ ≲ 0.4 ℓ_a). k = 2 (smoothing depth ≈ R) anomalous in the one point available (q_dens 0.44).

### ▶ Wednesday 2026-10-07 — server
- Shape sweep (`~/Overcharged/runs/shapes`, 675 runs) **finished 2026-10-07 18:44, 0 failures** → analyse (3 seeds per shape).
- Circle 10-seed queue **resumed 2026-10-07 19:18** (`runs/cavity/jobs_seeds10.txt`, 2170 jobs, 96 parallel, log `queue_seeds10.log`); none had finished before the pause (the 192 FAIL lines are the runs killed at the pause; they restart from scratch).
- SSH from Claude Code needs `SSH_AUTH_SOCK=/run/user/1000/gcr/ssh` when the session is started from a remote login (key is passphrase-protected).

### ▶ Tuesday 2026-10-06 evening — overnight jobs (office machine) — read this first
- User's list for tonight: (1) stage 2 on s02 + larger disorder ensembles; (2) active wandering vs passive; (3) stage 1 + mobility into master PDF and RESEARCH.md; (4) random walks with thermal + fluid noise for the normal spread; plus K1 ≠ K3 pinning.
- DONE (night): `reports/master/theory/fluct/stage2/` — 8 wolframscript processes `run_limited.wl 2_c00…07` (enlarged ensembles + s02 snapshots, ~5 h), then `queue3.sh` starts `3_c00…03` (realistic disorder) when the old e12 run ends. Logs `log2_c*.txt`, `log3_c*.txt`; results `M_2_c*.json`, `M_3_c*.json`. Stop: kill the wolframscript PIDs (`pgrep -f run_limited`), not pkill -f.
- DONE: `runs/fluidnoise/passive_R27_k06` (passive MPCD, coarseOut every 1 t.u.) → `theory/fluct/langevin/fluid_noise.py` then `langevin.py`.
- Done tonight: active wandering, exact MC, anisotropy pinning, master subsection `sections/29b_fluct.tex` (stage-2 paragraph to finish when the runs end).

### ▶ Monday 2026-10-05 afternoon — read this first
- **Master document (paper skeleton): `reports/master/master.pdf`** (built 2026-10-05 from `reports/master/master.tex` + `sections/*.tex`; new figure scripts in `reports/master/scripts/`: distributions, passive circles, Kaplan–Meier survival + docking, figure collection). New results in it: passive depth law ⟨s⟩ ≈ 0.56 + 0.094/ρ (all 468 passive circles); docking probability 4–10% (clean circles) vs 46–73% (elongated); one vortex → two at the same aspect as the layer loss (small ovals). Known code issues listed in its Appendix B.
- **Results map: `RESULTS.md`** (by question: validation, ℓ_a, circles, mirror, activity, shape, mechanism, documents) with copies of all figures/tables/movies/reports in `results/` (built by `tools/build_results.py`, copy only; rerun after new analyses). Nothing in `runs/` or `reports/` was moved.
- Power cut ~13:00 (office rebooted; nothing was running, no data lost; server unaffected). Office: nothing running.
- Done today (log entries dated 2026-10-05): flow/enstrophy/tilt in circles; wall Fourier for all circles (+k, −k); edge/bulk criteria over the whole sweep (knee = primary criterion); **wall dynamics in ovals/ellipses** (`runs/cassini/fig_wall_dynamics.png`): no sliding-ordered layer anywhere; wall +½ drift upstream, set by the vacancy fraction of the wall lattice (Spearman −0.82), not by pinning loss.
- Report updated 2026-10-05 ~16:00: new section "Monday" (`reports/project_report/monday_section.tex`; report.pdf now 21 pages). Brief updated to v7 (same link).
- **Server production binary fixed 2026-10-05 ~16:30** (seed honoured; = repo `main`; reproducibility test passed; details in TODO.md). Old binary kept in `~/Overcharged/mpcd_legacy_2026-09-28/`.
- **RUNNING (server, launched 2026-10-05 18:21 UTC): shape sweep, first pass** (user decision: ovals before the extra circle seeds, after finding Klein et al. 2025; see "Related literature"). `~/Overcharged/runs/shapes/` (generator `gen_shape_sweep.py`, also in `runs/cassini/`): circles (15 points), Cassini ovals b/a 3.0–1.1 (9 shapes, 135 points), ellipses ε 1.2–2.5 (5 aspects, 75 points) × R_eq = √(A/π) = 1.0, 1.5, 2.0 ℓ_a × ρ ≈ 0.04, 0.06, 0.09, 0.12, 0.14 (ρ_eff 0.036–0.146), 3 seeds = 675 runs, 5.6e12 particle-steps ≈ 2.5 days → ETA ≈ 8 Oct. Cassini binary (`~/Overcharged/niamh-overcharged-cassini/mpcd.out`, e63155b; seed honoured, checked in synopsis). Seeds 80000000 + 100·i_point + s (`manifest_shape_sweep.tsv`); `--nseeds 10` later adds only missing seeds (10 seeds total ≈ 8.5 days). Log `queue_shape_sweep.log`. Relaunch: `cd ~/Overcharged/runs/shapes && setsid -f nohup ./run_queue.sh jobs_shape_sweep.txt 96 >> queue_shape_sweep.log 2>&1 < /dev/null`.
- **DONE (office, 6 Oct morning): negative-winding shapes** `runs/cassini/P6_mirror_shapes/` (20/20, 0 failures); analysed (results log 2026-10-06, master document section "Negative winding in ellipses"). **DONE: Frank-constant calibration** `runs/Kcalib/` (3/3).
- **PAUSED: circle seeds s04–s10** (2170 runs, `~/Overcharged/runs/cavity/jobs_seeds10.txt`, stopped 2026-10-05 ~18:10 UTC, 0 finished; partial runs restart from scratch). Resume after the shape sweep with the same command as before: `cd ~/Overcharged/runs/cavity && setsid -f nohup ./run_queue.sh jobs_seeds10.txt 96 >> queue_seeds10.log 2>&1 < /dev/null`.
- **How to stop a queue** (the hint in `run_queue.sh` is incomplete): kill the `xargs` worker pool first (`ps -C xargs -o pid,cmd`, then `kill <pid>`), then `pkill -u baloewe -x mpcd.out`. Killing only `run_queue.sh` leaves xargs starting new jobs; never use `pkill -f <pattern>` inside an ssh command (it matches the ssh shell itself).
- Then: the oval sweep (seed count to decide; user wants 10 per point); share brief v7 with Tyler.

### Monday 2026-10-05 morning: weekend results
- **All weekend work done**: office set 52/52 (finished 4 Oct 01:50), server sweep 1304/1304 (3 Oct), all analysed. Report updated: `reports/project_report/report.pdf` (14 pages, new section "Weekend"); brief **v6** at the artifact link (ellipses, layer diagram, mirror runs).
- **H1:** elongation destroys the layer in **ovals and ellipses** (aspect ≳ 1.4–1.5) at matched size and density; f* collapses on aspect ratio with seeds; circles keep the layer for R ≈ 0.75–2 ℓ_a at all ρ ≤ 0.15 (`runs/cavity/analysis_weekend/layer_diagram_circles.png`), and the ovals sit inside that region in every size measure → shape effect.
- **k dependence:** circles: Q at fixed depth independent of k (traditional layer). Ovals: bounded excess (≤ ~+1) at moderate k, shrinking at high k (k = 12: plateau and wall lattice reappear) → crowded wall restores the topological structure.
- **H5 mirror:** bulk stays positive, reaches +1 only for |k| ≳ 13–15; activity breaks the mirror symmetry.
- **H7:** wall defects escape 10–20× faster in elongated cavities; the layer empties by undocking.
- New tools: `wall_fourier.py` (S(m), f_static), `charge_knee.py` (plateau value; smoothing 0.045), `exchange_rates.py`, `main_sweep_analysis.py`, `mirror_analysis.py`, `gen_ellipse.py`, `layer_area.py`.
- **Waiting on the user:** approve the oval sweep and its open choices (`reports/oval_sweep_draft.md`; validated binary ready on the server); share brief v7 with Tyler. Open analysis items: tilted plateau at Q ≈ 1.3–2.1 in large high-k circles (genuine less-charged layer?); plateau-width measure superseded by knee + f* + wall Fourier.

### Weekend plan 2026-10-02 → 05 (done; kept for reference)
- **Office queue** `runs/cassini/queue_weekend.log` (52 runs, 8 parallel, branch binary e63155b; jobs `jobs_weekend.txt`, started 16:22):
  W1 size collapse: half-length 48, b/a 1.3/1.6/2.0/2.5, k = 6, 3 seeds, 10k + 100k (`W1_size/`);
  W2 k scan: half-length 34, b/a 2.0 and 1.6, k = 2/4/8/10/12, 3 seeds, 10k + 100k (`W2_kscan/`; k = 6 = P3);
  W3 ellipses: area 2500, aspect ε = 1.2/1.5/2.0, k = 6, eMode 2, 2 seeds (`W3_ellipse/`, `gen_ellipse.py`);
  W4 extra seeds s02, s03 for H64 b/a 1.3 and 1.1 (10k + 40k, as P1; `W4_H64seeds/`).
- **Server:** sweep **complete** 3 Oct 12:30 UTC (0 failures), idle; new binary validated. **Oval sweep drafted, not launched:** `reports/oval_sweep_draft.md` (324 runs, ~1.2 days; generator `runs/cassini/gen_oval_sweep.py`); open choices listed there for the user. Branch `cassini` binary built in `~/Overcharged/niamh-overcharged-cassini` (e63155b; production binary untouched). User decision: **prepare an oval sweep (validate the new binary, draft job list) but do NOT launch until the user approves (Monday)**.
- **Analysis tasks for Claude:** main sweep circles (fetched to `runs/cavity/server_main/`: 834 active + 278 passive, 2026-10-02): Q(c), d*, f*, site rule vs R/ℓ_a and ρ (H2, H4); mirror runs (H5) when batch 3 ends; analyses of W1–W4 as they finish (f* collapse with seeds, k dependence, ellipses on the f* curve); plateau-width rework; exchange rates (H7); Monday: update report + brief.
- Watchers (background, this session): W1, W2, W3+W4 completion, server sweep completion.
- **User's ideas (2026-10-02 16:30), added to the weekend analysis:**
  1. *Traditional vs less-charged layer:* read the bulk charge from the knee of Q(c) (minimum of dQ/dc before the wall jump), not from Q = 1; layer charge = k − Q_knee; no knee = no layer. Sharpest test: **dq_b/dk** from the W2 k scan (traditional layer: dq_b/dk = 0, dq_layer/dk = 1; less-charged layer: 0 < dq_b/dk < 1). Recompute d*, f* at the knee.
  2. *Apples with apples:* match linear defect density ρ = 2k/P (small ovals at k = 6: 0.060–0.069 vs circle R 31.8: 0.060 — roughly matched) and activity number; for ovals test candidate lengths (√(A/π), half-length, inscribed/lobe radius, A/P) for which collapses the data; compare at equal √(A/π)/ℓ_a by interpolating between circle sizes.
  3. *Heat maps* of oval defect positions (time-averaged, wall overlay); statistic: site contrast (peak/background of the wall density) vs aspect.
  4. *Fourier along the wall arclength* (near-wall band): static structure (modes of the time-averaged wall density; pinned lattice → peaks at q = 2π(2k−2)/P), instantaneous S(q) = ⟨|Σ e^{iq s_j}|²⟩/N (peak survives sliding if ordered), and F(q,t) at the peak q (decay time = layer mobility). Pinned: both peaks; sliding-ordered: instantaneous only (suspected at R = 2.87 ℓ_a); no layer: neither.

### Overnight 2026-10-01 → read this first
- **H1 settled for k = 6 at ζ = 0.07: the charge split is geometric.** Small ovals at half-length 1.53 ℓ_a (3 seeds): circle-like at aspect 1.17 (b/a 2.5), partial at 1.29 (2.0), lost at ≥ 1.4 (b/a ≤ 1.6, still convex). Layer depth d*/ℓ_a 0.23 → 0.50 vs 0.20 for circles (`runs/cassini/fig_layer_depth.png`).
- Pinch-off: lobes don't decouple; a single lobe-shaped oval loses the layer too → shape (aspect), not the neck.
- High k (15–20) in the pinched oval partly restores the layer, total bulk ≈ +1 (≈ +½ per lobe) → crowded walls favour the topological picture. 1 seed, near the resolution limit.
- ζ = 0.09, 3 seeds: q_b 1.01–1.05, clean bulk charge (centre still depleted).
- Server batches 2+3 stuck at 544/1304 since ~17:00 (largest mirror runs, ~26 h each at 96-way); expect bunches of completions; ETA likely 4 Oct.
- Report: `reports/project_report/report.pdf`; brief v5 (with the corrected pinch-off bullet) at the artifact link.

### Running now
| Where | What | ETA |
|---|---|---|
| server | batches 2+3 (positive winding low ρ done; mirror −k running), 544/1304, 0 failures | ≈ 3–4 Oct (check rate) |
| office | Cassini P3: small ovals, half-length 1.53 ℓ_a, b/a 2.5/2.0/1.6/1.3, k = 6, 3 seeds (the H1 size control); high k = 12/15/18/20 on H48 b/a 1.02 | tonight |
| office | Cassini P1 H64 b/a 1.6, P2 H64 b/a 1.05 (last two) | ~18:00 |
| office | ζ = 0.09 seeds 2–3 (10/18) | ~19–20:00 |
| office | ~~overnight analysis~~ done 20:18 `runs/night_2026-10-01.sh` (log `runs/night_2026-10-01.log`, outputs `runs/night_2026-10-01/`): waits for all office runs, then accumulated-charge figure with the small ovals, cassini_charge / lobe_profiles / edge_criteria / bc_labels on P3 (+ circles 1.43, 1.72 for reference), renders + lobe rotation, ζ = 0.09 three-seed analysis, server status | when the runs finish (tonight) |

### Priorities (order agreed with the user 2026-10-01)
1. **H1 evidence (core of paper 1):** small-oval control (running) → `fig_accumulated_charge` update; seeds for the H48 ovals; same-area comparison. Ovals go into paper 1 if confirmed (user, 2026-10-01; to discuss with Tyler); ellipses = office exploration.
2. **Main sweep analysis** (server data, circles): Q(c) plateau vs R/ℓ_a and ρ (H2, H4); mirror runs (H5) when batch 3 is done.
3. **Mechanism:** where the BC defects go (label tracking), exchange rates, lobe vortices.
4. Wrap up H3 (ζ = 0.09 seeds) and the density/compressibility note for the SI.
5. Maintenance: server rebuild with the seed fix after the sweeps; checkpoint-restart crash (TODO); upstream bug report; merge `cassini` into `mpcd/` when validated.
6. Communication: update the collaborator brief (v5) with the H1 figure.
Parked: low-Mach (higher kBT) control — discuss with collaborator; channel. Dropped (user, 2026-10-01): b < a (two separate ovals) — just two independent cavities, nothing to learn. Instead, for the pinch-off limit (b/a → 1⁺, already run at 1.01–1.02): compare with a single lobe-sized cavity at winding ≈ k/2 (two decoupled lobes would give bulk ≈ +2 in total).

### If a machine went down
- Server: see `queue_chain.sh` / `run_queue.sh` (skip DONE runs): `cd ~/Overcharged/runs/cavity && setsid -f nohup ./run_queue.sh jobs_batch23.txt 96 > queue_batch23.log 2>&1 < /dev/null`.
- Office: rerun the queue of the affected set (`run_queue.sh` / `runs/cassini/run_queue_cassini.sh <jobs file> N`); DONE runs are skipped, unfinished ones restart.

## Related literature and novelty (checked 2026-10-05; web search, not exhaustive)
- **Norton, Baskaran, Opathalage, Langeslay, Fraden, Baskaran, Hagan, PRE 97, 012702 (2018)** (`1708.05773v2.pdf` in the project root): continuum (Beris–Edwards-type) model in disks, boundary charge −1, 0, +1. Circulating state (two co-rotating +½ inside, −½ at the wall) for every topology; in the turbulent state the boundary charge sits in a thin "topological boundary layer" around a neutral interior; hydrodynamic stability of wall-defect orientations. Message: **active dynamics are insensitive to topology.** Overlaps our circle results at |k| ≤ 1.
- **Klein, Soto Franco, Sabbir, Deutsch, Kliegman, Selinger, Mitchell, Beller, PNAS (2025), arXiv:2503.10880**: **disks with steady-winding anchoring of q = 1/2 … 5 (n = 1–10 excess +½) — our boundary condition**, but **below the turbulence threshold** (no pair creation): braided periodic orbits (golden/silver braids), flow gyres = 4|q−1|; "defects tend to self-screen excess topological charge by localizing near boundaries". Other shapes: epicycloids (cardioid, nephroid, trefoiloid) with tangential anchoring, cusps pin −½ (also arXiv:2401.04363, cardioid). **No ellipses/ovals, no turbulent regime.**
- Hardoüin et al., Nat. Commun. (2022), arXiv:2012.02740: active boundary layers in experiments; −½ accumulate at (slip) walls, topological polarisation. Hardoüin et al., Commun. Phys. (2019): channels, wall-created defects. Opathalage et al., PNAS (2019): circular wells, circulating state → turbulence. Theillard et al., Soft Matter (2017): active suspensions in circles and racetracks. Head et al., arXiv:2607.10234 (Jul 2026): 3D cylinders, defect delocalisation (no overlap).
- Own prior work: Loewe & Shendruk, Japan poster (`Poster_Japan_v4 (1).pdf`, old code): overcharged circles ±15, "topological Stern layer", "dynamical planar anchoring", tilt/shielding, flow speed.
- **Novelty after the search:** overcharged circles are no longer new as a setup (Klein 2025), but the **active-turbulent regime** (pair creation; exact 2k − 2 pinned wall lattice + bulk +1, mirror asymmetry with vacancies) is not covered; **no 2D active-nematic study of ellipses/Cassini ovals with imposed winding found** → the shape result is the main novelty ("topology is screened, geometry is not"). Klein's gyre-counting (defects vs vortices) is a natural theoretical hook for our one→two-vortex mechanism. Competition risk: the Beller group already combines winding BCs with non-circular shapes.
- To do: a citation search on Norton 2018 and Klein 2025 (Google Scholar "cited by") before writing.

## Goal
Active nematic (NIAMH-MPCD) in cavities whose homeotropic wall anchoring has an imposed winding k (`customAnch`). Known result: the excess charge collects at the wall, leaving an effectively planar-anchored bulk with net +1. Main question: **is this charge-layer effect topological or geometric?** Plan: validate on circles, then move to other shapes.

## Collaborator brief (web page)
- Private page "Overcharged Cavities": https://claude.ai/artifact/8yAVtVFycdqht5aJMvHKwG (v7, 2026-10-05: + wall-layer dynamics section with `fig_wall_dynamics.png`, mirror vacancies, criteria choice, seed fix on the server, 10-seed runs in "Running now"; v5 2026-10-01 ovals, v6 2026-10-04 weekend). Source: `reports/overcharged_brief.html`. To update from a later session, publish that file with this URL (keeps the link). Share it from the page's Share menu (the user does this).

## Spatial-analysis report (2026-09-30)
- Private page "Defect Geography in Cavities": https://claude.ai/artifact/FyTEbtvqPCy66sfz1vz3xp — source `reports/defect_geography.html`, figures `runs/cavity/analysis_spatial/report/` from `runs/cavity/spatial_report_figs.py` (all positions in r/R; g shown as ++ and −− curves, +− separate, as the user asked). Fig. 7: bulk g_++ collapses in r/R (∝ r/ℓ_a) at ζ = 0.05, 0.07, 0.09; g_−− only up to ~0.5 R at ζ = 0.05 (3-seed update 14:10).

## Code repository (2026-09-30)
- Private GitHub `baloewey/NIAMH-MPCD-overcharged` (`git@github.com:baloewey/NIAMH-MPCD-overcharged.git`), local clone `~/Physics/Overcharged/niamh-overcharged/`. `main` = the production code in `mpcd/`: upstream **2a38960** (PR #160, 2025-10-01) + 3 commits (seed fix, union-find tracker, customAnch/eMode 0/1/2); source identical to `mpcd/`, builds, same binary size. Remote `upstream` = Shendruk-Lab (fetch only, push disabled). Feature work (e.g. Cassini) goes on branches here.
- Note: `mpcd/` is based on 2a38960, not on a2cf7b2; the patches in `reports/upstream_patches/` were rebased onto a2cf7b2 for the upstream report. The Makefile fix (patch 0003) is **not** in production/`main`.

## Upstream bug report (2026-09-29)
- Private page "NIAMH-MPCD Issue Report": https://claude.ai/artifact/4jeoHSW715f71EBdfs3nqA — source `reports/upstream_bug_report.md` (paste-ready for GitHub issues) → `reports/upstream_bug_report.html` via `tools/md_to_artifact.py`.
- Patches on top of upstream master a2cf7b2 (2026-08-03): `reports/upstream_patches/fixes/` (0001 seed, 0002 tracker, 0003 Makefile) and `reports/upstream_patches/customanch/` (customAnch + eMode 0/1/2). All apply cleanly (`git am`) and build. Upstream (github.com/Shendruk-Lab/NIAMH-MPCD, public) is active with reviewed PRs (#170/#171, Aug 2026). Not pushed anywhere; the user decides how to send it.

## Publication plan (2026-09-28)
- **Paper 1 = circle + ellipse.** Aim to write sooner rather than later.
  - Circle: ℓ_a calibration, main sweep (ρ, R/ℓ_a), mirror (negative winding), H2 (low density), H4 (circulation), H3 (activity).
  - Ellipse: eMode 0 vs 1 (and arclength-uniform eMode 2 if needed to separate geometry from anchoring distribution).
- **Methods scope (user decision 2026-09-29):** the paper states that it uses nematic MPCD (cite the existing method papers) and does **not** report the ℓ_a(ζ) scaling, which two published papers already cover. ℓ_a enters only as the measured calibration ℓ_a(0.07) = 22.2 ± 0.3 used to normalise R and ρ. The exponent and compressibility analyses stay internal.
- **Future, separate project: parallelisation, needed for a jump to 3D** (user, 2026-09-29). 2D sweeps don't need it (running one job per thread already saturates memory bandwidth: 28.7 MPUPS with 96 jobs). In 3D a 100³ box at ρ̄ = 20 is 2×10⁷ particles ≈ 30 s/step serially → infeasible. Earlier attempt (Tim) scaled badly because of cache/memory stalls. Likely causes seen here: array-of-structs particles reached through per-cell linked lists in random memory order, one global RNG state, false sharing, load imbalance. Suggested path before any rewrite: sort particles by cell → per-thread xoshiro streams (jump function) → OpenMP over cells (dynamic schedule) → structure-of-arrays only if needed.
- **Channel (two walls, PBC):** only if needed. If the ellipse results are clear and compelling, either drop it or keep it for the first reply to referees. If run: both relative senses eventually, **opposite sense (ellipse limit) first**.
- **Plan change (user, 2026-10-01), pending the small-oval control and discussion with Tyler (collaborator):** if the geometric effect is confirmed, **the ovals (Cassini) go into paper 1** (more elegant and original). **Future server runs will be on ovals** (needs the `cassini` branch binary on the server: build it from the repo after the current sweeps, which also brings the seed fix). Ellipses become an office-machine exploration while the server runs ovals; add them to the paper if interesting. Question for the ellipses (user): can the layer be destroyed in convex cavities? Note: Cassini ovals with b/a ≥ √2 are convex too, and b/a = 1.6 already shows no layer (b/a 2.0, 2.5 running in P3), so ellipses test a second convex family (different curvature profile, no flattened sides).
- **Decision (user, 2026-10-01): working activity ζ = 0.07**, chosen to balance two effects: below it, onset effects (unstable wall layer, few dynamical defects; ζ = 0.05); above it, compressibility (centrifugal density stratification in the vortex, depleted centre; ζ ≥ 0.09). Present the results as the MPCD result: density stratification at finite Mach number happens in any particle-based solver, so state it, don't try to remove it. How much of the stratification analysis goes into the SI is to be decided with the collaborator. No low-Mach (higher kBT) control planned for now.
- **To discuss with the collaborator: a low-Mach control by raising kBT.** It would reduce the Mach number (c_s² ∝ kBT) and so the stratification, but it changes everything else too: viscosity, the nematic parameters and ℓ_a, so the whole calibration (ℓ_a, resolution limit, R/ℓ_a) would have to be redone. The user considers it a pain; worth raising only if a referee or the SI discussion asks how much the density stratification matters.


## Definitions
| Symbol | Meaning |
|---|---|
| ℓ_a | active length = √(A / ⟨N_defects⟩) in a large periodic box (`runs/activeLength`) |
| k | anchoring winding = `customAnch`; forces 2k +½ defects |
| ρ | bare linear defect density = 2k / (2πR) = k / (πR) |
| R/ℓ_a | activity number |
| ρ·ℓ_a | normalized density (maybe; decide once there's data) |
| Q(r) | enclosed charge = director winding on a circle of radius r (from `dirSOut`) |

---

## Hypotheses / open questions

### H1 — Topological or geometric? (main question)
Compare circles with other shapes at equal k (same topology) and at equal wall curvature or perimeter (same geometry).
- **Ellipse test design (user's plan + refinement, 2026-09-28):**
  1. **eMode 0 vs eMode 1 on the same ellipse** (user's test): same k, same shape, different distribution of the imposed winding along the wall (tip/flank density ratio (a/b)³ vs a/b). If the result is the same → consistent with a topological effect. If it differs → something beyond k matters (local anchoring distribution).
     **Constraint (2026-09-29):** only feasible for a/b ≲ 2 and mean ρ ≲ 0.05 (the local density limit at the tips; see ellipse notes).
  2. **Circle vs ellipse with arclength-uniform anchoring** (needs `eMode=2`, see TODO): same k, same uniform local density, different shape → isolates the geometry (curvature) itself. Needed because eMode 0's density follows the curvature, so a 0-vs-1 difference could still come from the shape acting through the anchoring.
  - With eMode 0/1/2 on one ellipse the non-uniformity steps (a/b)³ → a/b → 1; a smooth trend in the observables against that ratio would point to local density as the controlling factor.

### H2 — The effect disappears at low density (proposed 2026-09-28)
Forced wall defects sit ~1/ρ apart. If 1/ρ ≫ ℓ_a they no longer interact through the active length scale and may come off into the bulk. Then bulk charge > 1 and the charge layer fades.
- Prediction: crossover near **ρ·ℓ_a ≈ 1** → ρ ≈ 0.045 at act 0.07 (ℓ_a = 22.2 ± 0.3, L=200).
- Test: fixed-R scan, R = 100/π (R/ℓ_a ≈ 1.43), k = 2…10 → ρ = 0.02…0.10, ρ·ℓ_a = 0.44…2.2 (`runs/cavity/Rscan_R31.83`, done).
- If confirmed: extend the main sweep below ρ = 0.05 and confirm with fixed-ρ points (e.g. ρ = 0.03, k = 2, 3).
- **Result 2026-09-28 (R-scan): NOT supported at R/ℓ_a = 1.5** — bulk net charge ≈ +1 down to ρ = 0.02 (ρℓ_a = 0.43) once the edge/bulk cut is placed at the wall-layer boundary (≈ R − 4). May still apply at other R/ℓ_a.

### Reference: density effects paper (user is a co-author)
Kozhukhov, Loewe & Shendruk, "Mitigating density fluctuations in particle-based active nematic simulations", Commun. Phys. (2024), https://www.nature.com/articles/s42005-024-01734-2 — ℓ_d ∝ α^(−½) for all four activity formulations; activity-induced density fluctuations (feedback: dense regions → stronger activity) grow at high α (≳ α†); modulated formulations halve the maximum cell density. Same parameters as here (200², ρ = 20, U = 10, γ_R = 0.01, λ = 2, χ = 0.5, δt = 0.1). **Activity values are NOT directly comparable** to the current code (unknown edit shifted turbulence onset from ≈ 0.02 to higher activity; see TODO) → cutoffs must come from measured density fluctuations / Mach number, and results should be normalised by the measured ℓ_a.

### H3 — Activity enters only through ℓ_a
Results should collapse in terms of R/ℓ_a and ρ·ℓ_a, independent of the bare activity. The user has doubts.
- Test: rerun a subset of the cavity sweep at act 0.05 and 0.09.
- **If it holds:** activity becomes a control knob. Changing activity changes ℓ_a, so ρ·ℓ_a and R/ℓ_a can be tuned at *fixed geometry* (fixed R, k). This separates geometry from density cleanly and helps H1: a shape can be fixed while ρ·ℓ_a is scanned.
- Needs ℓ_a(act) from the active-length sweep, plus a check that it doesn't depend on box size.

### H4 — Bulk charge drops to +1 when the nematic flow starts to circulate (proposed 2026-09-28)
The reduction of the bulk net charge to +1 coincides with activity pushing the **fluid flow** into a net circulation (vortex), which also rotates the defects' orientation.
- **Order parameter (user's definition):** the absolute value of the time- and space-averaged azimuthal flow of each run, averaged over all runs (seeds) of a parameter set. Normalised by the mean flow speed so it's dimensionless and comparable across activities:
  - **Φ_persist = |⟨⟨v_φ⟩_x⟩_t| / ⟨|v|⟩** — the user's definition: persistent sense of rotation.
  - **Φ_inst = ⟨|⟨v_φ⟩_x|⟩_t / ⟨|v|⟩** — companion: instantaneous coherent circulation. It doesn't cancel when the vortex reverses (the sense of rotation is spontaneous and can flip within a run).
  - Φ_persist ≈ Φ_inst → steady vortex; Φ_inst ≫ Φ_persist → reversing vortex (count reversals too).
- Test: does the rise of Φ (vs R/ℓ_a, ρ·ℓ_a) coincide with the bulk charge dropping to +1 (P(q_bulk = 1) → 1, ⟨q_bulk⟩ → 1)? Secondary: do bulk-charge switches line up with vortex reversals within runs? Do the bulk +½ defects point along φ̂?
- Data: `flowOut` on **all seeds** (main sweep every 1000 steps, R-scan every 500; frames are interval averages). Tool: `runs/cavity/circulation.py`.
- Reference point: confined active nematics show a circulating two-(+½) state at R ~ ℓ_a before turbulence at larger R.

### H5 — Mirror runs with negative winding (customAnch = −k)
Total charge −k → at least 2|k| **−½** defects at the wall. −½ defects are not self-propelled, so a bulk +1 here can't come from wall defects swimming off. The excess over the natural +1 of a homeotropic circle is −(|k|+1), not k−1, so the two sets aren't mirror-symmetric even topologically. Same ρ definition with |k|: ρ = 2|k|/(2πR).

### Ellipse design notes (for H1; recorded 2026-09-28)
- **The imposed density isn't uniform along the wall.** Arclength per unit ellipse parameter t is √(a² sin²t + b² cos²t): b at the tips, a on the flanks. Tip/flank ratio of the local bare density ρ(s):
  - normal-angle wrapping (plain `customAnch`): (a/b)³ (follows the curvature κ);
  - parameter-uniform (`eMode=1`, `Hack_correction.nb`): a/b;
  - arclength-uniform (not implemented; see TODO): 1.
  → At high eccentricity, density and curvature both peak at the tips (a geometry/density confound). Compare `eMode=1` vs an arclength-uniform mode at the same a/b to separate them. Define ρ locally, ρ(s), or report the mean 2|k|/perimeter plus the tip/flank ratio.
- **Edge/bulk contour.** Circles no longer follow the wall. Use curves at fixed distance d inside the wall; they develop cusps once d > tip radius of curvature b²/a. For very elongated ellipses the tip wall layer is ill-defined → rely on the charge audit (`defect_audit.py`, needs no contour) and report Q on scaled/level-set contours with care.
- **Density limit applies LOCALLY (2026-09-29, from the circle resolution scan: anchoring fails where ρ(s) ≳ 0.15).** Local density ρ(s) = (1/π) dψ/ds. Peak/mean ratio and the resulting maximum mean density:

  | a/b | peak/mean eMode 0 | eMode 1 | arclength | max mean ρ eMode 0 | eMode 1 | arclength |
  |---|---|---|---|---|---|---|
  | 1.5 | 1.9 | 1.3 | 1 | 0.079 | 0.12 | 0.15 |
  | 2 | 3.1 | 1.5 | 1 | 0.049 | 0.097 | 0.15 |
  | 3 | 6.4 | 2.1 | 1 | 0.024 | 0.071 | 0.15 |
  | 4 | 10.9 | 2.7 | 1 | 0.014 | 0.055 | 0.15 |

  - eMode 0 at the tip: ρ_tip = k / (π r_tip), with r_tip = b²/a — the tip acts like a small circle of the tip radius of curvature.
  - **Conclusion (user): parameter-uniform (eMode 1) and arclength-uniform (eMode 2) are the only viable wrappings in general.** eMode 0 is usable only in a narrow window, **a/b ≲ 2 and mean ρ ≲ 0.05**. The eMode 0 vs 1 topology test (H1 step 1) must be done there (large perimeters → the expensive corner). Arclength-uniform is the only mode whose limit doesn't tighten with eccentricity → effectively required for high eccentricity.
- **Resolution:** require tip radius of curvature b²/a ≥ 5 cells (analogue of R ≥ 5), and check the local-density limit above for the chosen wrapping. Expect tracker errors (`miss`/`unpr` in the audit) to concentrate at the tips → audit every ellipse run.

### Channel geometry: infinite-eccentricity limit (proposed 2026-09-28)
Two vertical walls at x = 0 and x = W with custom homeotropic anchoring; periodic in y with period L_y. Vary the separation W.
- **Anchoring must depend on position along the wall**, ψ(y) = ψ₀ + 2π k y / L_y. The normal-angle hack gives no winding on a flat wall. Same mechanism as the arclength-uniform `eMode=2` → implement once (see TODO).
- **Net charge per period q = k_R − k_L** (both windings in the +y sense):
  - **opposite sense** (k_L = −k_R): q = 2k_R, same-sign defects on both walls. **The ellipse limit**: an ellipse of winding k ↔ k_R = −k_L = k/2 per period.
  - **same sense** (k_L = k_R): q = 0, +½ on one wall and −½ on the other. Charge-neutral control: can the wall layers' excess charges annihilate across the channel as W shrinks?
- Constraints: winding per period per wall ∈ ℤ/2 (director continuity across the PBC); ρ = 2|k_wall| / L_y; L_y ≥ several ℓ_a to avoid periodic-image effects → at fixed L_y, ρ comes in steps of 1/L_y. Parameters: W/ℓ_a, ρ (or ρℓ_a), relative sense.
- Observables: charge in strips at fixed x (edge/bulk by distance from the walls); **net channel flow ⟨v_y⟩** as the H4 analogue of circulation; the charge audit per period.
- Reference: confined active nematic channels go unidirectional flow → "dancing disclinations" → turbulence as W/ℓ_a grows (Shendruk et al., likely this same code). A baseline for untwisted anchoring.

### H6 — Wall-layer formation and stability (proposed 2026-09-29)
Define the wall layer from the **passive** run of the same (R, k): the activity-free position of the imposed defects. The inner edge d_p(R, k) comes from the passive radial defect distribution (e.g. 95% of passive defects outside it). At small k the passive defects may sit well away from the wall (k = 1: r ≈ 0.65R) → "no layer" rule needed.
- Observation (R-scan): the active layer is diffuse and thicker at low k, sharp at high k → interpret as a **forming / less stable layer** at low k.
- Observables vs the passive reference: layer width (radial spread, active vs passive), occupancy N_layer(t) and its fluctuations, exchange rate bulk ↔ layer.
- Data: `main_passive/` (batch 1).

### H7 — Docking dynamics at the layer edge (proposed 2026-09-29)
Bulk +½ defects that reach the inner edge of the layer circulate along it, "looking for a gap" to insert (user's observation).
- Observables: azimuthal speed of bulk +½ vs distance from the layer edge; docking events (time spent circulating before insertion, where it happens); azimuthal gaps between neighbouring layer defects vs the imposed spacing 2πR/(2k) → do insertions pick the largest gaps?
- Needs **defect trajectories** (frame-to-frame linking; tool to write). defectsOut raised to every 20 steps in batch 1 so linking works (displacement ≪ spacing) over most of the grid.

---

### H8 — Figure-8 cavity: separate circulation without changing topology (proposed 2026-09-30, user)
- **Conditional on H4:** worth doing if flow circulation turns out to matter for forming the wall layer (or the +1 bulk).
- **Idea:** take an elliptical cavity and pinch it vertically at the middle into a figure-8 (peanut) shape. It stays one simply connected cavity with one wall, so the imposed winding k and the topology are unchanged. But each lobe should host its own vortex: two circulation cells, probably counter-rotating, instead of one.
- **What it separates:** topology (same k, same total charge constraint) from flow organisation (one vortex vs two). Possible outcomes:
  - if the layer and bulk charge follow the topology: one +1 bulk for the whole cavity, wherever it ends up;
  - if they follow the circulation: each lobe behaves like its own cavity, e.g. +1 per lobe (bulk total +2, layer k − 2), with the neck acting as a barrier.
- **Control parameter:** the neck width w relative to the lobe size and ℓ_a. A wide neck should give one vortex, a narrow neck two; scanning w finds the crossover and whether the charge distribution switches with it.
- **Shape (checked in the source 2026-09-30):** NIAMH-MPCD walls are not a free implicit function. `calcW()` / `surf_func()` / `normal()` (`bc.c`, `mtools.c`) only handle W = Σ|A_i(x_i − Q_i)|^{P_i} − R^{P_4}, plus a wavy term B₀cos(B₁θ) (θ = polar angle for P = 2) and a non-4-fold branch (`ROTSYMM`, cos/sin of mφ/4). A Cassini oval ((x−a)² + y²)((x+a)² + y²) = b⁴ has an x²y² cross term, so it needs a small patch: a new shape flag in `calcW`, `surf_func` and `normal` (gradient: 4(x² + y² − a²)x, 4(x² + y² + a²)y). **No-patch alternative:** the existing wavy wall with P = 2, B = [B₀, 2, 0] gives r² = R² − B₀cos2θ, a peanut family from the circle (B₀ = 0) to two touching circles (|B₀| = R²); `normalWavy()` includes the cos term correctly. Concave neck: Cassini for b/a < √2 (waist/half-length < 0.58), wavy for |B₀| > R²/3 (waist/half-length < 0.71). Cassini is the nicer family: the neck closes smoothly at b = a (lemniscate), and for b < a it splits into two separate ovals, i.e. a genuine topology change, which could be a later control. Both have an explicit polar form r(θ) (Cassini: r² = a²cos2θ + √(b⁴ − a⁴sin²2θ)), so the eMode 2 arclength table is a numerical integral of √(r² + r′²).
- **Implementation needs:** a non-convex wall shape. Options: an implicit-surface BC if NIAMH-MPCD supports one (e.g. a Cassini oval or a higher-order superellipse term), or several overlapping walls, which needs care with the anchoring at the junctions. The arclength anchoring (eMode 2) needs an arclength table for the general shape. The site-based edge rule (`edge_sites.py`) and the charge audit work for any shape; `circulation.py` needs a per-lobe version (circulation about each lobe centre).
- **Confound: the pinch also changes the geometry.** The neck makes the wall concave (negative curvature), which circles and ellipses don't have, so a change in the charge distribution could come from curvature rather than from the two vortices. Design the neck-width scan around the one-vortex → two-vortex transition: if the charge distribution jumps at the transition while the curvature changes smoothly, it's the flow. Include a pinched cavity that still has one vortex as the geometric control.
- **Measure the rotation sense of each lobe** (co- vs counter-rotating) as a function of neck width, and whether the charge distribution follows it. Prior work, to verify before citing: Wioland, Woodhouse, Dunkel & Goldstein, Nature Physics (2016), "Ferromagnetic and antiferromagnetic order in bacterial vortex lattices": coupled chambers of bacterial suspension (polar active fluid, no nematic defects) where neighbouring vortices lock co- or counter-rotating depending on the gap. Suggests the lobes' relative sense may flip with neck width.
- Other design notes: each lobe needs R_lobe/ℓ_a ≳ 1.4 (steady vortex), so about twice the area and cost of the current R-scan runs; watch whether the neck traps defects or holds charge belonging to neither lobe (site rule + charge audit show this).
- **Patch plan:** `reports/cassini_patch_plan.md` (2026-09-30). To be implemented on a **separate copy of the code** (branch in `niamh-overcharged`, see Code repository), never in the production `mpcd/` tree; user preference for major changes.
- **Lobe coordinates (2026-10-01):** w = z² − a² maps the Cassini family to polar coordinates: |w| = d1·d2 = c² (levels; c = 0 the foci, c = a the lemniscate, c = b the wall), arg w = φ (orthogonal rectangular hyperbolas through the foci). Focus-to-wall curves z(c; φ) = ±√(a² + c² e^{iφ}), c ∈ [0, b] (+: right lobe). For c < a each level is a pair of lobe ovals, so the region c < a splits the cavity into lobes without a cut. Use for per-lobe enclosed charge Q(c) and lobe-centred profiles (level spacing in physical distance ∝ 1/|∇(d1 d2)|). Explorer: `lobeCoordsExplorer` in `runs/cassini/cassini_oval_manipulate.wl`.
- After paper 1 (circle + ellipse); fits a follow-up on geometry vs topology.

## Planned runs

| Run set | Status | Notes |
|---|---|---|
| ℓ_a(act), 100² PBC, act 0.050–0.090 | running, ends ~18:15 2026-09-28 | 30k steps, 1 seed |
| ℓ_a finite-size check, 200² PBC, act 0.05/0.07/0.09 | queued after the L100 sweep | 25k steps, overnight |
| Fixed-R scan R = 31.83, k = 2…10 (ρ = 0.02–0.10), act 0.07 | queued after the L100 sweep | 10k warm-up + 40k steps, dirSOut/500 |
| Main circular sweep: ρ = 0.05…0.5 (step 0.05), k = 1…20, R ≥ 5, act 0.07 | design agreed, not generated | 3 seeds first (10 later); 10k warm-up + 100k; first seed also writes flow/enstrophy + dirSOut at a low rate; ~1100 core-hours for 3 seeds (¾ from ρ = 0.05). Maybe start at ρ ≥ 0.10 |
| Activity-dependence subset (H3) | after the main sweep's first results | act 0.05, 0.09 |
| Mirror sweep, customAnch = −k, ρ ≥ 0.10 (H5) | generator ready (`gen_sweep.py --sign -1` → `km05` dirs, seeds +5000000), not generated | same grid/outputs as the main sweep, ~220 core-hours for 3 seeds |

---

## Results log

### 2026-10-06 — Wall-lattice fluctuations, mobility and pinning (`reports/master/theory/fluct/`: measure_passive.py, fluct_theory.wl, fig_fluct.py, mobility.py, pinning/, stage2/)
- **Harmonic elastic theory (Hessian + equipartition) vs passive data.** Internal tangential rms: measured/theory = 1.35 ± 4% for all circles (N = 4–20, R 31.83 and 27.28) and passive ellipses (eps 1.2–2.0) → one constant K_eff = 17.8 kBT (cf. K3 = 23, K1 = 42). Normal rms underpredicted by ~0.3–0.4 cells (core/detection). Figure `fluct.png`.
- **Sliding (rigid rotation) is a zero mode of the one-constant theory** (circle exactly; ellipse eigenvalue 1e-7–1e-4 → rms 84/10/2.7 cells at eps 1.2/1.5/2.0), but measured collective rms is 0.3–1 cell everywhere.
- **Mobility** (model-free, mu = Gamma_m <|u_m|^2>/kBT, passive circles): mu ≈ 0.03 (0.025–0.048) cells²/(kBT t.u.).
- **Pinning is not the grid** (galInv = 1, random shift): lattice registry with the anchoring pattern is reproducible across seeds and radii (active: Rbar of arg Z_{2k-2} 0.43–0.64 vs random ~0.2; passive Z_{2k}: 0.65–0.97 vs ~0.3). **Correction 2026-10-08: the passive agreement is inherited from the common initial condition (see 2026-10-08 entry).**
- **Tilt and position are locked by the anchoring:** tilt_j = (2k−2) phi_j + const. Passive: collective tilt vs lattice rotation slope 0.87–0.92 × (2k−2), corr 0.95–0.97 (k ≥ 5); per-defect mean tilt steps by −2·(2π/N) along a 2k ring. Active (N = 2k−2, ~150k frames, k = 3–10): D = arg Z − tau concentrated with Rbar 0.99, mean 0.00; within-run slope 1.04–1.15 (theory 1). Sites = positions where heads point straight at the wall; a coherent tilt tau requires a lattice rotation tau/(2k−2) (k = 6, 63°: 3.5 cells at R = 32). Pinning = cost of tilting a wall +½ (zero in one-constant Frank theory; nonzero with K1 ≠ K3 or orientation-dependent wall energy).
- **Exact potential + detector (mc_exact.py):** Metropolis on the full disk energy, wall offset 0.61, positions rounded to the cell grid (both trackers write the centre of mass of the flagged core cells; usual 2×2 cluster → integer coordinates in 93–98% of cases, rest at 1/6 or 1/2; same for legacy and default tracker — checked 2026-10-08): with **K = K3 = 23 kBT** the tangential rms matches within 2–4% (N = 8–20), depth within 0.05. K_eff = 17.8 of the harmonic fit was detector rounding + anharmonicity; tangential displacements cost bend. Normal: residual in quadrature ≈ 0.15 (N 20), 0.35 (N 12), 0.5 (N 4–8) cells → random-walk test with fluid noise (langevin/).
- **Splay ≠ bend pins the registry (anisotropy_pinning.py):** first order in Δ = (K1 − K3)/2 = 9.45 kBT on the exact one-constant field (smooth core cutoff). Active ring (2k−2 + central +1): minimum exactly at zero tilt (heads at the wall), barrier 80–110 kBT (k 6), predicted collective rms 0.45–0.51 cells vs 0.29–0.36 measured (pinned active circles, R 31.8). Passive 2k ring: barrier only 0.1–1 kBT, core-dependent → passive confinement unexplained at first order. One-constant models (Norton's Q-tensor) would have free sliding.
- **Active wandering (active_wander.py, sigma_gap = sd(gap)/√2 = stage-2 sigma_t):** active circles k 6, R/ℓa 0.4–2.2: 0.6–1.5 cells (passive level), count 10 in 75–93% of frames; R/ℓa 2.9 melts (≈ 9). Active ellipses: eps 1.2 1.3 (count 9.8); eps 1.5 6.3–6.5 (7.7 ± 1.1); eps 2.0 8.1–8.4 (7.3 ± 1.1), both seeds; free sliding. Realistic stage-2 disorder: sigma_t ≈ 1 (circles, 1.2), 6–8 with 2–3 vacancies (1.5–2).
- **Random walks (langevin/):** fluid noise measured in a passive run with instantaneous cell velocities (runs/fluidnoise/passive_R27_k06): rms 0.042–0.049 cells/t.u., tau 0.9–2.6, isotropic. Langevin (full potential, K3, mu 0.03, thermal + OU fluid advection, rounding): fluid ×1 adds 2–5%; closing the normal gap needs ×3, which overshoots the tangential (N 12: 1.36–1.89 vs 1.15). Normal excess is direction-specific (detector near the wall), not noise. Matched mode times → mu ≈ 0.02.
- **Stage 2, measured snapshots s01:** eps 1.5 and 2.0: top mode one-vortex in 40/40 snapshots each, while the measured flow is two-vortex in 72% / 55% of intervals; gain ratio only slightly lower in two-vortex intervals (1.22 vs 1.29; 1.16 vs 1.20); measured tilt pattern barely projects on the top mode (0.16–0.17 vs random 0.13–0.15). Overnight: enlarged ensembles (sigma 0.5–6, 1–3 vacancies, 60 samples, eps 1–2.5; stage2/cfg_2_*), realistic disorder (sigma 6/8 × 2/3 vacancies; cfg_3_*, queued), s02 snapshots.
- **Stage 2 final (2026-10-06 night, stage2b_summary.txt, 4845 configs):** P(top mode two-vortex) = 0 up to eps 1.5 for every class (sigma_t 0.5–8, 0–3 vacancies, 60–80 samples); eps 2.0: 1/60 only in the two most disordered classes; eps 2.5: 8–15% with realistic disorder (sigma 6–8, 2–3 vac), 2–5% milder. Median gain ratio 1.11–1.32. Measured snapshots s01 + s02 (eps 1.2/1.5/2.0): 240/240 one-vortex top mode while the flow is two-vortex in 68–72% (1.5) and 48–55% (2.0) of intervals. → wall feedback never selects two vortices at eps ≤ 2; the bulk does.
- **Stage 2 (imperfect lattices, Stokes feedback modes):** positional disorder (sigma_t 1, 2, 4 cells) and 1–2 vacancies never make the two-vortex mode the top mode (0/20 in every case, eps 1–2); median gain ratio one/two 1.08–1.28. Measured-snapshot comparison pending.

### 2026-10-06 — Wall tilt patterns in ellipses: coherence and Stokes feedback modes (`reports/master/scripts/wall_tilt_coherence.py`, `reports/master/theory/stokes_tilt_modes.wl`)
- Measured tilt coherence of the wall layer: 1.00 (circle), 0.98–0.99 (ε 1.2), 0.43 (ε 1.5), 0.34 (ε 2.0); the tilt agrees with its own half's rotation for 70–79% of wall defects.
- Static split layer (left +63°, right −63°): two counter-rotating vortices, no net rotation. Per-half rotation relative to a coherent layer: 0.68 (circle), 0.73 (1.2), 0.94 (2.0).
- Feedback operator M_ij (tangential velocity at i from a tangential force at j): the uniform (one-vortex) mode has the top gain for every ε up to 2.5. The dipolar (two-vortex) mode is second, and the margin shrinks from 13% (circle) to 5% (ε 2.5).
- Conclusion: the wall hydrodynamics allow two vortices but never select them; elongation makes the two nearly degenerate, so a small bulk push selects the two-vortex state. Consistent with the vortex leading the tilt in circles. Added to master §hydro (Tables tiltsplit, tiltmodes).

### 2026-10-06 — Hydrodynamics in ellipses (`reports/master/theory/stokes_ellipse_static.wl`, `stokes_frames/prep_shape.py`, `solve.wl`)
- **Static:** a coherently tilted wall layer (10 +½, evenly spaced, depth 1.5, tilt 63° all the same way) drives a single vortex in every ellipse (area 2500): halves same sense, whole-cavity rotation 0.050 (circle) → 0.048 / 0.043 / 0.033 f/η (ε 1.2 / 1.5 / 2.0).
- **Measured flows (W3, seed 1):** halves same sense 100% / 28% / 45% of intervals, whole-cavity rotation 0.57 / 0.26 / 0.09 (ε 1.2 / 1.5 / 2.0).
- **Frame-by-frame fit on measured defects:** R² all / wall only / bulk only 0.72/0.70/0.21, 0.45/0.39/0.06, 0.34/0.23/0.11; f/η per wall defect 10.0 / 9.5 / 7.6; half-rotation correlation model–measured 0.25–0.43 (ε 1.2, single steady vortex: little variance), **0.94 (ε 1.5), 0.85 (ε 2.0)**; model same-sense fraction 100 / 40 / 55%.
- **Reading:** the wall defects' tilts follow the local vortex and their forces reproduce and reinforce the two-vortex dynamics, but a coherent layer would drive one vortex → the two-vortex state originates in the bulk; the wall layer's share of the flow falls with elongation (70% → 23%). Master document: hydrodynamics subsection, paragraph "Ellipses".

### 2026-10-06 — Hydrodynamics: Stokes model of defect forces (user's idea; `reports/master/theory/stokes_circle.wl`, `stokes_frames/`)
- Static Stokes flows (Mathematica FEM, no-slip disk R = 31.83, smeared point forces σ = 1.5): +½ = force along its orientation, −½ = three inward forces at 120° (tripole). **No-slip screening:** a tangential force at depth 1–2 cells drives 5–7× less circulation than one deep inside; a wall layer of N = 2…18 at the predicted depths (5 → 1.5 cells) gives only 3× more circulation for 9× more defects. Tripole flows ≈ 0.
- **Measured circulation is k-independent** (fixed-radius scan, k = 2–10: |⟨v_θ⟩| 0.61–0.72, same profile, peak ≈ 1 at r ≈ 0.6–0.8 R).
- **Frame-by-frame regression** (measured defects of each 50-t.u. flow interval → Stokes flows of wall +½, bulk +½, −½; regress measured cell velocities; k = 2, 6, 10, 40 intervals each): wall +½ alone explain 67% (k6) and 76% (k10) of the flow field, bulk +½ 8%, tripoles nothing; circulation model/measured 0.52/0.55 (k6), 0.64/0.68 (k10); fails at k = 2 (R² 0.25, circulation 0.25 vs 0.61). Fitted force per wall defect f/η = 20 / 11.5 / 9.5 (k 2/6/10) → crowding weakens the push; bulk coefficient negative.
- **Force scale from the code:** collOp 21 active stress ζ_eff ≈ 4 (per-particle kick ζ N_cell/2 ≈ 0.7, dipole ≈ 3.5 per cell), MPCD viscosity η ≈ 17 (2D MPCD-AT formulas, ρ = 20, Δt = 0.1) → F = π ζ_eff S ρ_cut ≈ 11 ρ_cut → fitted f/η ⇔ ρ_cut ≈ 31 / 18 / 15 cells: physically plausible, shrinking with crowding. (Corrects the earlier "√(K/ζ) ≈ ℓ_a": with ζ_eff, ℓ_a ≈ 8 √(K/ζ_eff).)
- **Timing:** fluctuation cross-correlation (steady vortices, 1170 intervals) peaks at lag 0 to −50 with a slight asymmetry favouring circulation leading (±150: 0.36 vs 0.22). At vortex reversals (ζ = 0.05, k = 8–10, 11 matched events): **circulation flips first in 9/11**, median lead ≈ 90 t.u. (mean 50 ± 30).
- **Reading:** two-way coupling. The vortex sets the wall tilt (follower in time); the tilted wall layer then sustains the circulation (engine in magnitude), with screening + weaker crowded forces keeping it k-independent. Reversals start in the bulk.

### 2026-10-06 — Negative winding in ellipses: the shape effect is mainly active (`runs/cassini/P6_mirror_shapes/`, `reports/master/data/mirror_shapes_wall.tsv`, `shielding_mirror_*.tsv`, figure `reports/master/figs/mirror_ellipses.png`)
- Circle and ellipses ε 1.2/1.5/2.0, area 2500, eMode 2, k = −6 and −12, 2 seeds; passive references: exactly 2|k| −½, no +½ (eMode 2 with negative winding validated).
- **The −½ wall layer survives elongation:** occupancy 0.93–0.99; escape rate of wall −½ rises only ≈3× circle → ε 2 (k = −6: 2.3 → 6.4–7.4 ×10⁻⁴), vs ≈20× for wall +½ with positive winding. Lattice order along the wall drops at k = −6 (pinning 0.56–0.77 → 0.22–0.32) but defects stay at the wall.
- **Shielding weakens modestly:** P(annihilation | bulk +½ meets wall −½) k = −6: 0.28 (circle) → 0.31 (ε 2); k = −12: 0.14 → 0.21; bulk −½ targets 0.51–0.59. Positive bulk at k = −6: ≈0.70 → 0.47; at k = −12 no shape effect (0.79–0.94).
- **Conclusion:** geometry weakens the elastic shielding for both signs, modestly (as the theory says); the dramatic layer loss with +k needs self-propelled wall defects carried off by the two-vortex flow → the shape effect is mainly active, enabled by geometry.
- Also: in the −k circle the wall −½ creep along the wall (median 11–12 cells per 500 time units, 65% beyond half a spacing), unlike +½ in circles.
- Frank constant size check: 96×96 box K1 = 41.5, K3 = 23.0 k_BT (64×64: 42, 23).

### 2026-10-05 evening — Frank constants of the MPCD nematic (`runs/Kcalib/`, `frank_constant.py`, `reports/master/data/frank_constants.txt`)
- Passive defect-free periodic 64×64 boxes (density 20, mfpot U = 10 k_BT, collOp 21, act 0, aligned along x, 5k + 20k steps, dirSOut 20), 2 seeds; equilibrium director fluctuations, A k_BT/⟨|θ_q|²⟩ = K1 q_y² + K3 q_x², |q| < 0.8 with a noise constant.
- **K1 (splay) = 42.0 k_BT, K3 (bend) = 23.0 k_BT** (seeds agree within 1 k_BT), one-constant K ≈ 32.5 k_BT → πK ≈ 100 k_BT. ⟨S⟩ = 0.876. 96×96 box running (size check).
- Theory energies in k_BT: approach barrier ≈ 4, pull 4→1 cell ≈ 30 (free pair 70), repulsion between wall defects ≈ 60, leftover field in ε = 2 ellipse ≈ 6 (inner 60%) / 33 (inside 0.8 level), circle 0.3.
- ~~Consistency: √(K/ζ) = 21.5 cells vs ℓ_a = 22.2~~ — coincidence with the bare ζ; corrected 2026-10-06: collOp 21 active stress ζ_eff ≈ 4 → √(K/ζ_eff) ≈ 2.8, ℓ_a ≈ 8 √(K/ζ_eff).
- Literature: 3D nematic MPCD isotropic, K ≈ 113 U per unit length (Shendruk & Yeomans 2015); 2D here anisotropic.

### 2026-10-05 evening — Elastic shielding, analytic (user's poster free energy; `reports/master/theory/shielding_disk.wl`, `shielding_disk_results.txt`, `approach_energy.png`)
- Model: one-constant nematic, strong anchoring = Neumann condition for the field ψ conjugate to the director angle → **same-sign images** at R²/z̄ (wall repulsive). E/(πK) = −Σ_{i≠j} q_i q_j ln|z_i − z_j| − Σ_{i,j} q_i q_j ln|R² − z_i z̄_j|; reproduces the poster's U_ij and U_W (poster writes −E).
- **Ring identity:** N evenly spaced wall charges at r_w act on a bulk defect only through ln|1 − (r/r_w)^N e^{iNφ}| + ln|1 − (r r_w/R²)^N e^{iNφ}|: zero on angular average (2D shell theorem, images included) → **the bulk does not feel the winding k** (exact in a circle); corrugation ∝ (r/r_w)^N. Interior variation over r < 0.7R: 0.019 πK (k = −6, R = 31.8).
- **Passive depth law derived:** ring minimum r* = R[(N−1)/(3N−1)]^{1/2N} (checked: dE/dr = 0, d²E/dr² > 0); depth → (ln 3/4π)·spacing = 0.087·spacing. **Against all 310 passive circles (ρ ≤ 0.15): measured = 0.984 × predicted + 0.61 cells, r = 0.996**; offset 0.61 (IQR 0.60–0.62) identical for both signs (effective wall ~half a cell inside). With a central bulk charge Q_b: r* = R[(N−1+2Q_b/q)/(3N−1+2Q_b/q)]^{1/2N}.
- **Approach of a bulk +½ to the wall −½ layer** (k = −6, R = 31.8, 14 sites at the predicted depth): energy released 8 → 4 cells: free pair 0.35 πK, towards a wall −½ −0.03 (a small barrier ≈ 0.04 πK at 4–5 cells); 4 → 1 cell: free 0.69, wall 0.30 (43%); between two wall −½: repulsive all the way (+0.6 πK at contact). → screened, short-range attraction reachable only along narrow channels; matches the shielding test (suppression in the approach; annihilation ratio 0.41).
- Ellipse: first numerical attempt (method of fundamental solutions, `ellipse_shielding.py`) failed its own disk check (ring depth 1.97 vs exact 1.50, circle tip ratio 1.40) → discarded. **Redone by conformal mapping** (`shielding_ellipse.wl`; w = m^{1/4} sn((2K/π) arcsin(z/c), m), Neumann Green's function + one harmonic correction H with ∂_nH = 2π/P − |w'|): map exact to 1e−13, ε → 1 reproduces the disk depth (1.33 vs 1.328). Passive equilibrium of 12 +½ (area 2500): **evenly spaced in arclength** (flank gap/tip gap 1.00–1.01), not the conductor distribution; leftover field from the layer on a bulk +½: inside the 0.6 level 0.003 (circle) → 0.008 / 0.021 / 0.059 πK (ε 1.2 / 1.5 / 2.0), weak compared with near-wall approach energies (0.3–0.7 πK). → elastic shielding explains the circle; in ellipses it is imperfect but weak, so the flow is the likely driver of the layer loss.
- **Correction (2026-10-05):** the earlier "passive tip excess ≈ √ε" (count of wall defects per length at tips vs flanks) was a discreteness artefact: every passive ellipse has exactly 6 of 12 defects in the tip regions. Gap-based measure (`reports/master/scripts/wall_gaps.py`, `data/wall_gaps.txt`): passive flank gap/tip gap 0.97–0.98 (uniform), active 0.77–0.87 (tips depleted; the active conclusion stands). Master document corrected (Sec. 16.1) and theory section added (Sec. 19).

### 2026-10-05 evening — Shielding test in circles (user's argument; `reports/master/scripts/shielding_test.py`, `reports/master/data/shielding_both.tsv`, 114 runs per sign, ρ 0.05–0.12, R 1.0–2.3 ℓ_a, 3 seeds)
- Encounter = first approach of a mobile defect within 4 cells of an opposite-sign defect; annihilation = both tracks end within 10 time units, ≤ 4 cells apart. Wall targets (s < 4) vs bulk targets (s > 6), same runs.
- **−k:** bulk +½ meeting a wall −½ annihilate in 0.218 ± 0.002 of encounters vs 0.534 ± 0.002 with bulk −½ (per run 0.22 vs 0.54). Suppression mostly in the approach: 28% of wall encounters get within 2 cells (bulk 57%); once close 0.46 vs 0.57. → screened attraction (hydrodynamic contribution not excluded).
- **+k:** bulk −½ rarely reach the wall +½ (2,635 vs 78,007 encounters with bulk +½); 0.31 vs 0.52.
- → The edge is shielded in circles for both signs. Next: the same on the −k ellipses (P6, running overnight).

### 2026-10-05 15:30 — Wall dynamics in ovals and ellipses: upstream drift into vacancies (`runs/cassini/wall_fourier_dyn.py` → `wall_fourier_dyn.csv`; circles k = 6 on the server → `wall_fourier_dyn_circles_k06.csv`; figure `runs/cassini/fig_wall_dynamics.{png,pdf}`)
- Same static + dynamic measures as the circle analysis below, with positions by arclength along the real wall; 81 active oval/ellipse runs (P1–P4, W1–W4), 0 failures. Static numbers reproduce the 2–3 Oct values exactly.
- **No sliding-ordered layer in ovals either.** Where the lattice exists it is pinned (lattice phase speed ≈ 0). The lattice-phase sliding speed is noise where S/N ≲ 0.4 (weak order, tiny and reversing wall flow), so it is replaced by a **per-defect drift**: wall +½ matched frame to frame (≤ 1.5 cells), ds/dt along the wall, signed with the near-wall tangential flow v_w (cells 1–5 from the wall), block standard errors.
- Near-wall flow turns with the bulk vortex (sign agreement 100% with one strong vortex; 60–90% in elongated cavities, where |v_w| drops from 0.2–0.35 to 0.03–0.07).
- **Wall +½ drift against the wall flow, and the drift follows the vacancies of the wall lattice, not the loss of pinning:** drift/|v_w| vs vacancy fraction 1 − N_wall/(2k − 2): Spearman −0.82 over all runs, −0.93 among pinned lattices (f_static ≥ 0.75); vs f_static only +0.18. Examples: circles k = 6 drift ≈ 0 (|ratio| < 1%) for R = 0.7–1.7 ℓ_a, −2% at 2.15, −18% at 2.87 ℓ_a (occupancy 0.78); oval H48 b/a 2.0 is pinned (f 0.84) but 41% empty and drifts at −69%; near-circular ovals/ellipse (aspect 1.18–1.2, full lattice) −1 to −9%; elongated ovals −10 to −80%. Reading: defects hop upstream into empty sites (vacancies move downstream).
- k scan (half-length 1.53 ℓ_a): drift large at k = 2–4 (−25% at aspect 1.29, −66 to −69% at 1.51), ≈ 0 from k = 6 (1.29) and k = 8 (1.51) on, as the wall fills (occupancy → 0.9–0.97). Pinched H48 b/a 1.02, k 12–20: drift ≈ 0. Crowded walls stop the drift as they restore the lattice.
- Caveats: the drift ratio is noisy for tiny circles without a vortex (R < 0.4 ℓ_a) and for 1-seed P1/P2 cavities; drift includes defects entering/leaving the band.

### 2026-10-05 12:35 — Edge/bulk criteria compared over the whole circle sweep (`runs/cavity/criteria_all.py` → `analysis_weekend/criteria_all.csv`; `plots_criteria.py` → `analysis_weekend/criteria_maps.pdf`; 930 runs, ρ ≤ 0.15, ±k)
- Criteria: knee of Q(r), fixed layer r < 11R/12, pinned-site rule, density minimum, mobility (slow), slow + tilted (+k only).
- **+k:** all criteria agree on q_b ≈ 1 for R ≈ 0.8–2 ℓ_a at every ρ. The knee is the most robust (≈ 1 from ~0.6 to ~2.4 ℓ_a) but becomes undefined (no plateau) at large R and low ρ. Fixed-layer cut fails at small R (layer thicker than R/12 → q_b ≫ 1 below ~1 ℓ_a) and large R. Density-minimum and mobility give q_b → 0 below ~0.6–0.8 ℓ_a (no bulk defects). Above ~2–2.5 ℓ_a all criteria except the knee rise above 1 (defects leave the layer), earliest at high ρ. Slow + tilted overcounts (1.3–1.6) at high ρ, where the wall tilt is smaller (30–40°, flow entry), so the tilt cut misclassifies wall defects there.
- **−|k|:** knee q_b rises with ρ and R (0.3–0.5 at ρ ≤ 0.05; → 0.9–1 only at ρ ≳ 0.10 and R ≳ 1.5 ℓ_a); site rule and mobility saturate at ~0.8. Consistent with the vacant −½ wall lattice (wall Fourier entry below). Fixed-layer cut is wrong for −k at R ≲ 1.5 ℓ_a (strongly negative).
- → For the paper: knee as the primary bulk-charge criterion, density minimum / mobility as cross-checks; don't use a fixed-depth cut.

### 2026-10-05 — Flow, enstrophy and wall-defect tilt in circles (930 runs, ρ ≤ 0.15; `runs/cavity/flow_reduce.py` on the server → `flow_tilt.py` → `plots_flow.py`, output `runs/cavity/analysis_flow/flow_maps.pdf`)
- Mean flow |⟨v_θ⟩| (per-frame |spatial mean|, time-averaged) collapses on R/ℓ_a for all ρ: rises to ~0.35 at R ≈ 0.6 ℓ_a, ~0.9–1.0 at 2–2.5 ℓ_a (Mach ~ 1). +k slightly faster than −k at R ≈ 1–2 ℓ_a. +k collapses at R ≳ 2.5 ℓ_a (circulation fraction 0.9 → 0.2–0.4, reversals): the same R where the layer is lost. Circulation fraction ≈ 0.8–0.9 between R ≈ 0.5 and 2.5 ℓ_a; below 0.4 ℓ_a no coherent vortex (k = 2: reversals, agree ≈ 0.5).
- Enstrophy ⟨ω²/2⟩ (interior cells) peaks ~0.01–0.015 at R ≈ 0.4–0.5 ℓ_a and decays to 0.004 at 2.9 ℓ_a; +k above −k by ~20–30% at R < 2 ℓ_a.
- Tilt Δ = ψ − φ (code orientation − polar angle), signed with sign⟨v_θ⟩: wall +½ lean downstream by ~60–65° at ρ ≤ 0.06, decreasing to ~30–40° at ρ = 0.14–0.15 (crowded wall); Δ drops to ~0 where the layer is lost (R ≳ 2.5 ℓ_a) and at R < 0.4 ℓ_a. Wall −½ (three-fold) lean downstream by ~8–15°, peak ~25° at R ≈ 0.5 ℓ_a. The sign agrees with the circulation in ~100% of frames inside the layer window.

### 2026-10-05 — Wall Fourier for all circles, +k and mirror, static + dynamic (`runs/cavity/wall_fourier_all.py`, `plots_fourier.py`, output `runs/cavity/analysis_fourier/fourier_maps.pdf`)
- Lattice period m* = m_lat (2k − 2 for +k, **2|k| + 2 for −k**, i.e. the period of a +1 bulk) for every run with R ≳ 0.5–0.75 ℓ_a, both windings.
- **Mirror runs: the wall −½ lattice has the +1-bulk period but is under-occupied** (N_wall/m_lat 0.7–0.9 at low ρ, → 1 at high ρ / large |k|): the vacancies are why the mirror bulk is below +1. The bulk charge from the wall count (k − sign·N_wall/2) tracks q_slow (r = 0.96) and q_dens (median diff 0.05).
- +k: occupancy 1.0 from 0.75 to 2.3 ℓ_a, drops to 0.6–0.85 at R ≳ 2.5 ℓ_a (layer loss); instantaneous order S/N 0.8–0.9 (low ρ) vs 0.5–0.8 (high ρ); f_static ≈ 0.95–1 (pinned). Mirror: S/N only 0.4–0.6, f_static 0.5–0.9 rising with R.
- Dynamics: no steady sliding anywhere (lattice speed < 2–3% of |⟨v_θ⟩| for R > 0.6 ℓ_a; only tiny cavities rotate). Lattice memory F(m*, 1000) ≈ 1 for +k in the layer window; mirror 0.5–0.9. Rearrangement time τ_fl ~ 10 time units for +k at low ρ, 100–1000 at high ρ (+k) and ~60–100 for −k. The "sliding-ordered" class never appears in circles: the layer is pinned or vacant/melting.

### 2026-10-04 — W3 ellipses and W4 extra seeds; final f* figure (`runs/cassini/W3_ellipse/`, `W4_H64seeds/`, `fig_layer_area_all.{png,pdf}`)
- Ellipse support added to `lobe_profiles.profile` (levels = scaled copies of the wall, c = √((x/ε)² + y²), area inside c = π ε c²).
- **Ellipses (area 2500, k = 6, eMode 2, 2 seeds): ε 1.2: knee 0.98–1.07, f* 0.31, pinned lattice (f_static 0.98) — circle-like; ε 1.5: no knee, f* 0.51–0.52, f_static 0.33–0.38; ε 2.0: no knee, f* 0.57, f_static 0.06–0.15.** On the Cassini curve up to aspect 1.5; slightly below at 2.0 (0.57 vs 0.60–0.67). → **User's question "can the layer be destroyed in convex cavities?": yes — elongation (aspect ≳ 1.4–1.5) is the generic control for both convex families.**
- W4 (H64, 3 seeds now): b/a 1.3: f* 0.53–0.60, d*/ℓ_a 0.78–1.00 (peak near b/a 1.3 holds; seed 1 was the top of the range); b/a 1.1: f* 0.66–0.74, d*/ℓ_a 0.56–0.62. At the largest size f* sits slightly lower than for smaller ovals at the same aspect (~0.56 vs 0.63–0.66 at aspect 1.94): collapse good but not perfect at 2.88 ℓ_a.
- Weekend office queue finished 2026-10-04 01:50 (52/52).

### 2026-10-03 — W2 k scan (half-length 34, b/a 2.0 and 1.6, k = 2…12, 3 seeds; `runs/cassini/W2_kscan/kscan_summary.txt`, `wall_fourier.txt`)
- **Idea 1 answered with fixed-depth Q** (works without a knee): circle R 31.8 (k = 2–10): Q(c/b = 0.85) = 0.94–0.99, independent of k (slope −0.002) → traditional layer, all extra winding goes to the wall.
- **Ovals: bounded, non-monotonic excess.** b/a 1.6: Q(0.85) = 1.30 (k2), 2.06 (k4), 1.96 (k6), 1.63 (k8), 1.42 (k10), 1.31 (k12); knees reappear at k = 10 (1/3, 1.54) and k = 12 (3/3, 1.40). b/a 2.0: 1.12, 1.54, 1.51, 1.49, 1.36, 1.24; knee at k = 12 (3/3, 1.31). Neither "traditional" (constant) nor "all winding into the bulk" (slope ~1): the bulk carries at most ~+1 extra at moderate k (ρ ≈ 0.04–0.07), and the excess shrinks as the wall gets crowded (k = 12, ρ ≈ 0.13).
- Wall Fourier (1 seed): b/a 1.6 lattice reappears with k (f_static 0.03, 0.21, 0.11, 0.39, 0.70; S 1.2 → 5.3 for k 2 → 12); b/a 2.0 partly pinned throughout (0.61–0.96), occupancy → 2k − 2 at high k.
- Same pattern as the high-k pinched oval and the mirror runs: **a crowded wall restores the traditional (topological) structure**; elongation matters most at moderate imposed density.

### 2026-10-03 — W1 size collapse with seeds (half-length 48 = 2.16 ℓ_a, k = 6, 3 seeds; `runs/cassini/W1_size/`, `fig_layer_area_seeds.{png,pdf}`)
- **f* collapse confirmed with seeds:** f* at half-length 34 / 48: b/a 2.5 0.29 / 0.28 · 2.0 0.39 / 0.42–0.45 · 1.6 0.53 / 0.49–0.54 · 1.3 0.63 / 0.65–0.67 (seed spread ≤ 0.03). Aspect ratio is the control parameter at fixed k.
- Knee: none at half-length 48 for any shape, including b/a 2.5 (H34: knee 1.05). Wall Fourier: b/a 2.5 still a pinned 10-site lattice (S 7.1–7.3, f_static 0.97); 2.0 partly pinned (S 3.0–3.9, f 0.81–0.88, ~6 of 10 sites occupied); 1.6, 1.3 no lattice (S 1.5–2.7, f 0–0.24, one seed 0.67). At √(A/π) ≈ 2 ℓ_a the bulk side (plateau) degrades first while the wall lattice survives — as for circles at 2.0–2.5 ℓ_a.

### 2026-10-03 10:30 — H5 mirror runs (winding −|k|, 570 runs; `runs/cavity/mirror_analysis.py`, `analysis_weekend/mirror.csv`, `mirror_summary.txt`, `mirror_Q_check.png`)
- Server sweep complete 2026-10-03 12:30 UTC (1304/1304, 0 failures); all data fetched (`runs/cavity/server_main/`: 834 + 278 positive, 570 + 190 mirror).
- **Under negative winding the bulk charge is positive, not −1**: Q(r) rises to a positive peak, then plunges to −|k| at the wall. Peak (bulk) charge vs |k| (ρ ≤ 0.15): 0.00 (1), 0.09 (2), 0.27 (4), 0.47 (6), 0.68 (8), 0.83 (10), 0.92 (12), 1.00–1.02 (15–20). Wall −½ per frame ≈ 2|k| + 2·(bulk charge) (e.g. |k| = 10: 20.8; |k| = 20: 41.9 ≈ 2|k| + 2).
- At fixed R/ℓ_a the bulk charge rises with |k| (R 1.0–1.5 ℓ_a: 0.42 → 0.64 → 0.78 → 0.84 → 0.93 for |k| 3–4 → 13–16); at fixed |k| it peaks at R ≈ 1.5–2 ℓ_a and falls for larger cavities.
- **Reading: +1 (Euler characteristic) is the attractor for both winding signs, but mirror symmetry is broken by activity**: with +k the +1 bulk appears at every k ≥ 2 (made of self-propelled +½); with −k the wall must absorb two extra −½ (their +½ partners form the bulk +1), which happens gradually and fully only for |k| ≳ 13–15. Consistent with the high-k pinched oval: a crowded wall favours the topological picture.
- **Server binary validation (2026-10-03):** branch `cassini` build (`~/Overcharged/niamh-overcharged-cassini/mpcd.out`, e63155b) — circle k = 3 passive: alignment 0.971 (production 0.972, office 0.971), Cassini a = 0: 0.970, audit 100%, seed honoured (`runs/cassini/server_validation/`). The server production binary predates the tracker patch (rejects "defectTracker"): the whole sweep ran with the legacy tracker, as intended.

### 2026-10-02 21:15 — wall ↔ bulk exchange rates, H7 (`runs/cavity/exchange_rates.py`, `analysis_weekend/exchange_rates.txt`)
- Gap-closed +½ tracks; wall = within 3 cells, bulk = beyond 6 (hysteresis); escape rate = visits ending in the bulk per wall defect per time unit; residence from fully observed visits.
- Escape rate (×10⁻³ per time unit) / mean residence: circle 1.43 ℓ_a 0.35–0.48 / 880–950; 2.15 ℓ_a 0.28 / 2050; 2.87 ℓ_a **4.5** / 220; oval b/a 2.5 0.68–0.77 / 730–970; 2.0 2.7–3.4 / 200–250; 1.6 **5.8–6.2** / 140–150; 1.3 5.5–5.8 / 150; lobe reference (k = 3) 7.6 / 125.
- **Wall defects in elongated ovals escape 10–20× faster than in circles of the same size**, rising monotonically with elongation; the large circle reaches oval-like rates (degraded layer). ≤ 10% of wall visits end in annihilation: the layer empties by undocking.

### 2026-10-02 20:30 — main sweep (circles), layer diagram (`runs/cavity/main_sweep_analysis.py`, `analysis_weekend/main_sweep.csv`, `main_sweep_summary.txt`, `layer_diagram_circles.{png,pdf}`)
- 1112 runs analysed (834 active ζ = 0.07, k = 1–20, ρ = 0.03–0.50, 3 seeds; 278 passive), 0 failures: knee of Q(r), d*, f*, site rule, wall Fourier per run.
- **Knee finder fixed** (`charge_knee.py`, smoothing 9 bins = 0.045 in c/b): with 5 bins it caught the flat shoulder of the cell-quantised layer steps (e.g. k = 5 at 0.80–0.90 ℓ_a gave a spurious knee at Q ≈ 1.9; true plateau at +1.0). Recomputed for all runs; known cases unchanged (circles ≈ 1.0, oval b/a 2.5 1.03, elongated ovals: no knee).
- **Circles (ρ ≤ 0.15): +1 bulk with a pinned wall layer from R ≈ 0.75 to 2.0–2.3 ℓ_a for all ρ** (knee in 89–100% of runs, Q_knee ≈ 1.0, f_static ≥ 0.78); degrades above ~2.5 ℓ_a (knee in 15%) and in very small cavities (< 0.5–0.6 ℓ_a: Q_knee → k or none, no room for a bulk). Beyond the resolution limit (ρ ≥ 0.2) mostly no knee.
- **Inscribed-radius test passed:** at ρ ≈ 0.06, circles with R = 0.72 and 0.96 ℓ_a: knee 1.03–1.07, f_static 0.97–0.98, d* 0.16–0.17 ℓ_a. The ovals that lose the layer have inscribed radii 0.79–1.01 ℓ_a and √(A/π) 1.15–1.8 ℓ_a, **inside the circle layer region** in every size measure → the loss is a shape effect.
- Possible genuine "less-charged layer": large circles at high ρ (k = 15–20, ρ 0.10–0.12, R 2.0–2.6 ℓ_a): tilted plateau at Q ≈ 1.3–2.1 with all 2k − 2 sites pinned (f_static 0.65–0.95). To check (near the resolution limit and the large-R degradation).
- **Figure:** layer diagram in (ρ, R/ℓ_a) with the ovals overlaid at (2k/P, √(A/π)/ℓ_a): ovals sit in the middle of the circle layer region but have no knee (b/a 2.5 has one).

### 2026-10-02 17:30 — weekend analysis, part 1: wall Fourier, knee, matching table (`runs/cavity/wall_fourier.py`, `runs/cassini/charge_knee.py`, `runs/cavity/analysis_weekend/`)
- **Wall Fourier (idea 4)** (+½ within 4 cells of the wall, arclength along the wall; m* = peak of S(m), f_static = |⟨ρ_m⟩|²/⟨|ρ_m|²⟩): passive circle k5: m* 10 = 2k, S 8.1, f 0.94; circles k6 at 1.43 / 2.15 ℓ_a: m* 10 = 2k − 2, S 9.0, f 0.99 (pinned lattice); 2.87 ℓ_a: S 7.1, f 0.87 (partly sliding); ovals H34 b/a 2.5: S 8.1, f 0.99 (pinned) · 2.0: 6.7, 0.97 · 1.6, 1.3: **S 2.1, f 0.33** (weak order, not pinned); pinched k6: S 1.7, f 0.59; pinched k18: m* 30 (2k − 2 = 34), S 4.8, f 0.74 (partial lattice). Decay time of the fluctuating part only meaningful when f is small (b/a 1.6: ~190, 1.3: ~370 time units). Maps: `wall_fourier_test_maps.png`. (Bug fixed: off-by-one in the peak index; band 4 cells.)
- **Knee of Q(c) (idea 1)** (local slope minimum after the −½ core, flatter than both sides, reached by a rise ≥ 0.5): circles up to 2.15 ℓ_a: Q_knee 0.97–1.09 (k = 3 and 6); 2.87 ℓ_a: knee in 1 of 3 seeds; oval b/a 2.5: 1.05–1.08; **all ovals b/a ≤ 2.0 (all sizes, k up to 20): no knee**; lobe reference: weak knee (0.77) in 1 of 3. No "less-charged layer" (knee at +2, +3) found so far.
- **Matching (idea 2)** (`matching_table.txt`): oval H34 b/a 2.5 ≡ circle R 1.43 ℓ_a in ρ (0.060) and all lengths (√(A/π) 1.42, half-length 1.53, inscribed 1.30, 2A/P 1.41) — and behaves identically. Ovals b/a 1.6, 1.3 matched in ρ (0.066–0.069); their √(A/π) 1.27/1.15, 2A/P 1.23/1.07 are within the clean-circle range (≥ 1.07); only the **inscribed radius (1.01 / 0.79 ℓ_a) is below all circles compared so far**. Test: main-sweep circles at ρ ≈ 0.06 with R = 10.6 / 15.9 / 21.2 (0.48 / 0.72 / 0.96 ℓ_a, k = 2/3/4).
- Main sweep analysis running (`main_sweep_analysis.py` → `analysis_weekend/main_sweep.csv`, 1112 runs).

### 2026-10-01 22:40 — layer area fraction f* (user's idea; `runs/cassini/layer_area.py`, `layer_area_all.txt`, `fig_layer_area.{png,pdf}`)
- f* = A(c > c*) / A_total (fraction of the cavity outside the level where Q = +1); w_eff = A(c > c*) / P (equivalent strip width).
- Circles (k = 6, 1.43–2.15 ℓ_a): f* 0.20–0.26; k = 3 circles (1.08–1.43 ℓ_a) 0.26–0.32; 2.87 ℓ_a 0.20–0.27. w_eff/ℓ_a 0.17–0.25 for every clean layer (0.29–0.39 at 2.87).
- Small ovals (1.53 ℓ_a): f* 0.29 (b/a 2.5) → 0.39 (2.0) → 0.53 (1.6) → 0.63 (1.3); larger ovals 0.49–0.71; lobe-shaped oval 0.47; pinched oval k 6 → 20: 0.71 → 0.58.
- **Threshold f* = 1/2** ("layer" larger than the bulk it surrounds): all circles ≤ 0.32, ovals that lost the layer 0.5–0.71, partial cases (b/a 2.0, lobe reference) 0.39–0.47 → a continuous crossover; 1/2 is the defensible bound beyond which it is not a layer.
- **Against aspect ratio, f* of all ovals collapses onto one curve regardless of size** (1.53, 2.16, 2.88 ℓ_a), crossing 1/2 at aspect ≈ 1.5; the lobe reference falls on it. w_eff does not collapse (grows with size). Interpretation: for a true layer f* ∝ width/size (falls with size, as in circles 1.43 → 2.15 ℓ_a), whereas a charge distribution spread through the cavity gives a size-independent f*. **The size scaling of f* itself diagnoses "layer vs no layer"**, and aspect ratio is the control parameter.

### 2026-10-01 21:50 — lobe-shaped reference: the lobe shape alone loses the layer (`runs/cassini/P4_lobe_ref/`)
- Single convex Cassini oval with one lobe's area and aspect (b/a 1.8, half-length 30.44, area 2171, aspect 1.38), k = 3, 3 seeds, 10k + 100k steps: **3 sites (circle: 2k − 2 = 4), q_b ≈ 1.9 of 3, d*/ℓ_a 0.31–0.32** (k = 3 circles of the same area 0.18–0.22), Q(c) smooth (1.29 / 1.66 / 2.07 at c/b 0.80 / 0.85 / 0.90), labels at the wall after 1000: 0.16–0.28.
- → The pinched-oval lobes behave like single cavities of their shape; **no neck coupling needed**. Consistent with the small-oval crossover in aspect (1.17 circle-like, 1.29 partial, ≥ 1.38 lost). The control parameter looks like the aspect ratio / elongation, not the neck.
- Final project report with all of today: `reports/project_report/report.pdf` (11 pages).

### 2026-10-01 20:30 — overnight analysis: small ovals (P3, 3 seeds), high k, ζ = 0.09 (`runs/night_2026-10-01/`)
- **Small ovals, half-length 1.53 ℓ_a, k = 6, 3 seeds — sharp crossover with elongation:**
  b/a 2.5 (aspect 1.17): **circle-like** — 10 sites (= 2k − 2), q_b 1.07–1.11, P(+1) 0.80–0.85, labels at the wall after 1000: 0.82–0.86, mobility/tilt agree with the site rule 97–98%, d*/ℓ_a 0.22–0.23.
  b/a 2.0 (aspect 1.29): partial — 8–9 sites, q_b 1.6–2.1, labels 0.60–0.65, criteria agree 81–86%, d*/ℓ_a 0.29–0.31.
  b/a 1.6, 1.3: **layer lost** — q_b 3.3–3.6, labels 0.37–0.50, criteria agree 60–73%, d*/ℓ_a 0.41–0.50.
- **High k on the pinched oval (H48 b/a 1.02, k = 6 → 20; imposed density 0.05 → 0.17):** the layer partly re-forms. Total Q(c) at c/b 0.85 / 0.90: k6 1.76 / 2.48 · k12 1.43 / 2.00 · k15 1.15 / 1.51 · k18 1.06 / 1.36 · k20 1.03 / 1.27, then a sharp jump within the last ~5% (layer). Each lobe carries ~+0.5–0.65 before the jump at k ≥ 15 (Q_L ≈ Q_R), i.e. **a total bulk ≈ +1 for the whole cavity, not +1 per lobe**. d*/ℓ_a 0.38 → 0.36 → 0.33 → 0.30 → 0.30 (circle 0.20); labels at the wall after 1000: 0.76 (k12), 0.84 (k15, k18), 0.63 (k20, near the resolution limit). A crowded wall (high imposed density) stabilises the layer and pushes the oval back towards the topological picture (total +1). k = 18–20 is at the edge of the anchoring resolution.
- **ζ = 0.09, 3 seeds:** 2k − 2 sites in all runs, q_b 1.01–1.05 (SE ≤ 0.03), P(+1) 0.77–0.92, Φ 0.88–0.99 — clean like 0.07 in the bulk charge (centre density depletion still applies).
- Server: batches 2+3 at 544/1304 since ~17:00 (no new completions in 3 h, 96 running): the largest mirror runs (R up to 64, ~26 h each at 96-way) are going through; completions will come in bunches. ETA probably 4 Oct rather than 3 Oct.

### 2026-10-01 20:20 — lobe-sized references for the pinch-off (`runs/cavity/server_k03/`, `runs/cassini/layer_depth_k03_circles.txt`; P4 queued)
- Lobe of the H64 b/a 1.01 oval: area 2171, aspect 1.39 (≈ circle R 26.3 = 1.18 ℓ_a). **k = 3 circles (half the winding, server sweep, 3 seeds): R = 23.9 (1.08 ℓ_a) d* = 4.1 cells (0.18–0.19 ℓ_a), R = 31.8 (1.43 ℓ_a) 4.5–5.0 cells (0.20–0.22 ℓ_a); both with a +1 plateau (slope ≈ 3) and the jump to 3 at the wall.** The pinched-oval lobes cross +1 ~10 cells from the wall with no plateau (slope ≈ 9) → a lobe does not behave like a circle of its size.
- Queued 20:16 (office, `runs/cassini/P4_lobe_ref/H30.4_ba1.8_k03/s0{1,2,3}`, `queue_lobe_ref.log`): single convex oval b/a 1.8, half-length 30.44 (area 2171, aspect 1.38 = lobe), k = 3, 3 seeds, 10k + 100k steps. Layer → neck coupling destroys it in the pinched oval; no layer → lobe shape suffices. ETA ~22:00.

### 2026-10-01 20:20 — layer depth d* from Q(c) (user's idea; `runs/cassini/layer_depth.py`, `layer_depth_k06.txt`, `fig_layer_depth.{png,pdf}`)
- c* = outermost upward crossing of Q(c) = +1; d* = arclength-weighted mean wall distance along the level c = c* (circle: R − c*). Normalised by ℓ_a (better than by size: d*/R is ~0.10–0.14 in all circles).
- **Circles** (3 seeds): d*/ℓ_a = 0.19–0.20 (R 1.43 ℓ_a), 0.20–0.21 (1.72), 0.22–0.27 (2.15), 0.31–0.42 (2.87); d* ≈ 4.5–6 cells up to 2.15 ℓ_a.
- **Size control passed: small ovals at half-length 1.53 ℓ_a (3 seeds, P3; b/a 2.5 seed 3 still running):** d*/ℓ_a = 0.23 (b/a 2.5), 0.29–0.31 (2.0), 0.41 (1.6), 0.49–0.50 (1.3); seeds agree within ~0.02. Same size as the clean circles, d* doubles with elongation, back to the circle value as the oval rounds → **geometric effect confirmed** (convex b/a 2.0, 1.6 included).
- Larger ovals (1 seed): half-length 2.16 ℓ_a 0.38–0.67, 2.88 ℓ_a 0.46–1.00. **Non-monotonic in the neck:** maximum at b/a ≈ 1.3 ((a/b)² ≈ 0.6), decreasing towards the lemniscate (b/a 1.1 → 1.01: 0.58 → 0.46 at 2.88 ℓ_a). **Correction (user's question, same evening; `runs/cassini/pinch_lobe_charge.txt`): this is NOT a restored layer.** Per-lobe Q_L, Q_R (c < a) for b/a 1.01–1.1 pass through +1 smoothly at c/b ≈ 0.86–0.89 with slope dQ/d(c/b) ≈ 8–10 (circle plateau ≈ 3, circle layer jump ≈ 50–100); the whole cavity passes +2 without a plateau; Q_L ≈ Q_R; curves barely change from 1.1 to 1.01. The d* drop is because each lobe is smaller. Lobes of ~1–1.4 ℓ_a (size of a clean circle) do not behave like circles even with a 13-cell neck → lobe shape or neck coupling still matters. Reference needed: a single lobe-shaped cavity at winding ≈ k/2.
- Plateau width (|Q − 1| < 0.15 band) is unreliable for ovals (negative values from band selection when Q just passes through 1); circles ≈ 2.4–2.9 cells. To be reworked; d* is the robust number.

### 2026-10-01 18:00 — where the wall (BC) defects go: label tracking (`runs/cavity/bc_labels.py`, `analysis_spatial/bc_labels_k06.{png,txt}`)
- Labels at reference times every 100 time units: circles = +½ on pinned sites; ovals = slow & tilted (window = P10–P90 of this run's slow near-wall tilt: 0.57–0.75 to 1.66–1.78 rad), s < 6. Label follows its gap-closed track; at a paired annihilation it passes to the birth partner of the −½. (A cycle guard never triggered except once; slowness was the O(N²) event pairing, now windowed.)
- **Circles** (3 seeds each): labels identified per reference ≈ 2k − 2 = 10 (7.5–10 at 2.87); fraction still within 3 cells of the wall after 1000 time units: **0.87–0.94 (1.43 ℓ_a), 0.84–0.88 (2.15), 0.41–0.59 (2.87)**; median wall distance of tracked labels stays 1.7–2.3; 2–17% (1.43) up to 31–45% (2.87) left (s > 6) and came back.
- **Ovals** (1 seed each): only 2.4–5 labels per reference; fraction at the wall after 200 / 1000: b/a 1.6: 0.35 / 0.09 · 1.3: 0.27 / 0.11 · 1.1 (H48): 0.43 / 0.35 · 1.02: 0.52 / 0.39 · 1.1 (H64): 0.13 / 0.06. Median wall distance of tracked labels rises to 4–8 cells within ~100–200 time units. **Wall defects in ovals leave the wall within ~100–200 time units and wander into the bulk; in circles they stay for ≫ 1000.** Tighter necks (1.1, 1.02) hold them somewhat longer.
- Caveat: label chains are lost more often in ovals (tracked fraction after 1000: 0.2–0.8 vs ≥ 0.75 in circles); a lost label counts as "not at the wall" in the at-wall fraction, while the median distance uses tracked labels only — both show the same trend.

### 2026-10-01 18:20 — gap-closing defect tracker (`runs/cavity/defect_tracks.py`)
- Frame-to-frame optimal assignment per sign + gap closing (up to 4 missing frames, radius rlink + 1.5·(g−1)·dt) + pairing check (creation/annihilation must be ±½ pairs: track start/end with an opposite-sign start/end within 6 cells and 5 frames).
- With the old 3-cell link radius ~90% of track breaks were jumps of 3–6 cells per frame (normal motion at ~1 cell/time unit, dt = 2, cell-quantised positions), not missing detections: median track 1 frame, only ~20% of starts/ends paired.
- **Settings rlink = 8, gap = 4:** circle 2.15 ℓ_a: 16,800 → 1,821 tracks; paired starts/ends 0.92/0.94 (circle 1.43: 0.95/0.96; ovals b/a 1.6, 1.02: 0.90/0.95, 0.89/0.95). Gap closing adds little once rlink is right.
- Only 1–3 tracks span > 90% of a 10,000-time-unit circle run (all at s ≈ 1.8, wall defects): the 2k − 2 wall defects are themselves replaced over time (undocking/docking or annihilation with a −½ near the wall), so following "the BC defects" needs the label-passing rule (label goes to the surviving partner at annihilation).

### 2026-10-01 17:45 — defect density in Cassini lobe coordinates (`runs/cassini/lobe_profiles.py`, `lobe_profiles_k06.png/.csv`)
- Shells between Cassini levels c (c = √(d1 d2); 0 foci, a lemniscate, b wall), area from a fine grid; n_±(c) = defects per shell area per frame; Q(c) = accumulated charge inside level c. For a = 0 this is the circle's radial profile.
- **Circles (k = 6, 1.43 / 2.15 / 2.87 ℓ_a):** +½ peak at c/b ≈ 0.97 (layer), dip just inside (gap), Q(c) has a knee/plateau at ≈ 1 around c/b ≈ 0.85–0.92, and a negative dip towards the centre (−½ excess: Q min −0.25 / −0.8 / −1.1).
- **Active ovals (H48 b/a 1.6, 1.1, 1.02; H64 1.1):** n_+ rises steadily towards the wall with no peak–gap structure; Q(c) rises smoothly, no plateau, never negative (few −½: e.g. 2.1 per frame at b/a 1.6 vs ~4.7 in the circle of the same size). The excess charge is spread over the whole cavity rather than split into layer + bulk. Passive oval: all +½ at c/b > 0.82 as expected.
- **Summary figure:** `runs/cassini/fig_accumulated_charge.{png,pdf}` (`fig_accumulated_charge.py`; picks up the P3 small ovals automatically when done). Circles: negative dip (−½ core), plateau at +1 for c/b ≈ 0.85–0.93, jump to k within the last few % (layer). Ovals: monotonic, never negative, no plateau, crossing +1 at c/b ≈ 0.78; **all four ovals (b/a 1.6–1.02, half-length 2.16 and 2.88 ℓ_a) fall on nearly one curve**, so the neck width barely matters once the shape is elongated. User (2026-10-01): this is the plot to show the qualitative circle/oval difference; a diffusive layer is hardly a layer.
- Caveat: near the wall the c-shells have non-uniform physical width (wider at the neck), which would smear a thin layer; the wall-distance profiles (16:45 entry) show the same absence of a layer, so the conclusion does not rest on c alone. Per-lobe Q at c = a is only meaningful when a/b is well below 1 (at b/a = 1.02 the lemniscate is almost the wall).

### 2026-10-01 17:20 — cut-free edge criteria (mobility, orientation), tested on circles, applied to ovals (`runs/cavity/edge_criteria.py`, `analysis_spatial/edge_criteria_circles.png`, `edge_criteria_apply.txt`)
- Per +½: mobility D = net displacement over 20 time units along its linked track; orientation Δ = defect angle (defects.dat col. 4) − outward wall normal.
- **Circles (site rule as reference):** wall defects are slow (median D 1.0–1.4 vs 10–11 cells in the bulk) and **tilted at |Δ| ≈ 1.1 rad (≈ 63°) to the normal, with the sign of the vortex sense** (6/6 runs: Δ = +1.1 ↔ v_φ > 0). Bulk +½ cluster near |Δ| ≈ π − 0.5. Each criterion alone agrees with the site rule on 97–99% of +½ at R/ℓ_a ≤ 2.15 (0.89 at 2.87).
- Edge-like = slow (D < 5) and/or tilted (0.5 < |Δ| < 1.8, s < 6), counted per frame (occurrences / frames used; normalising by all +½ biases towards edge because bulk tracks break more often — fixed). Bulk charge q_b = k − (edge-like)/2, k = 6:
  circle 1.43 ℓ_a: tilted 1.07, slow 1.05, both 1.06 (site rule 1.05) · 2.15: 1.11 / 1.25 / 1.24 · 2.87: 1.93 / 2.71 / 2.66 (site 1.85; part of the layer slides: tilted but mobile)
  ovals half-length 2.16 ℓ_a, b/a 1.6 / 1.1 / 1.02: tilted 3.42 / 3.58 / 3.26, slow 3.66 / 3.44 / 3.10, both 3.85 / 3.87 / 3.57 · half-length 2.88, b/a 1.1: 4.03 / 4.25 / 4.45
- **Result: with any of the cut-free criteria, the ovals have q_b ≈ 3.1–3.9 (edge holding only ~4–5.5 of the 12 imposed +½) against ≈ 1.1–1.25 for circles of the same half-length.** The orientation criterion (most generous, counts sliding wall defects) still gives 3.3–3.6. Holds already at b/a = 1.6 (convex). Pending: the small-oval control (P3, half-length 1.53 ℓ_a), seeds, and the same-area comparison.

### 2026-10-01 16:45 — H1 first look: bulk vs edge charge, ovals vs circles of the same size (k = 6, ζ = 0.07; `runs/cassini/size_vs_shape_k06.txt`, `cassini_charge.py`)
- In the active ovals the site rule finds only 3–7 sites (circles: 2k − 2): **the wall layer largely detaches**. Defect profile vs distance from the oval wall: passive: all 12 +½ at s = 1–4; active: only ~4–5 +½ within s < 3, the rest spread into the interior with no gap; near-wall speed 0.34–0.40 (passive jitter 0.30).
- Size control from the server sweep (circles k = 6, 3 seeds each, 100k steps; `runs/cavity/server_k06/`, fetched 2026-10-01): near-wall +½ (s < 3) and enclosed charge Q(s > 5):
  circle R/ℓ_a 1.43: 9.7, 0.99 · 1.72: 9.6, 0.99 · 2.15: 8.8, 1.06 · 2.87: 6.2, **1.80**
  ovals, half-length 2.16 ℓ_a (b/a 1.6, 1.3, 1.1, 1.02): 4.6, 4.4, 5.0, 5.2 near-wall; Q(s > 5) **2.6, 2.7, 2.5, 2.3** · half-length 2.88 (b/a 1.1): 3.4, **3.3**
- **Correction (same day):** in the ovals Q(s) has **no plateau** (falls steadily, ~4.5 at s = 2 → 1.5 at s = 10), so the quoted oval "bulk" values are just Q at s = 5 and cut-dependent; the near-wall count uses s < 3 borrowed from passive runs/circles. Robust statement: at every s the ovals have fewer defects near the wall than circles of the same size, i.e. **no separate edge layer**; report full Q(s) curves (s/ℓ_a) rather than a single bulk number. Testing cut-free criteria (mobility, orientation) on circles first.
- **Both size and shape matter.** Circles keep the +1 bulk / (k − 1) layer split up to R/ℓ_a ≈ 2.2 and start to lose it at 2.9 (bulk ≈ 1.8). Ovals of the same half-length lose it much more (bulk ≈ 2.3–2.7 at 2.16, 3.3 at 2.88), already at b/a = 1.6 where the wall is convex. So the charge split is **geometric**, not fixed by topology (all b/a > 1 shapes are simply connected). Possible mechanism: the two vortices of an elongated cavity meet along the wall (stagnation regions) and push layer defects into the bulk; a single circular vortex shears along the wall uniformly.
- Caveats: ovals 1 seed, 40k steps (circles 3 seeds, 100k); "same size" taken as half-length = R (oval area is smaller: H48 b/a 1.6 area 5000 vs circle R = 47.7 area 7160); ρ = 2k/P is similar (~0.04–0.05). Needs seeds, a smaller-oval series (half-length ≲ 1.7 ℓ_a, where circles are clean), and rounder ovals (b/a 2–3).
- Also for H2/H4 / main sweep: circles at R/ℓ_a ≈ 2.9 already lose the clean +1 bulk at ζ = 0.07.

### 2026-10-01 16:20 — Cassini lobe rotation, 5 finished active runs (`runs/cassini/lobe_rotation.py`, movies in `runs/cassini/movies/`)
- Mean azimuthal velocity about each half's centroid (last 80% of 4,000 time units), co-rotation fraction, reversals:
  H48 b/a 1.6 (waist 64): L +0.32 / R −0.26, |v_φ|/|v| 0.67/0.54, co-rotating 15% of frames, 3/3 reversals, whole-cavity rotation +0.05.
  H48 b/a 1.3 (49): +0.18 / −0.17, 0.48/0.44, co-rot 34%, 8/5.  H48 b/a 1.1 (30): −0.04 / −0.12, 0.54/0.60, co-rot 75%, 7/4.
  H48 b/a 1.02 (13.5): +0.19 / −0.27, 0.64/0.72, co-rot 14%, 2/0.  H64 b/a 1.1 (40): −0.23 / −0.25, 0.49/0.54, co-rot 88%, 4/4.
- **The convex "control" (b/a = 1.6) already has two counter-rotating vortices**, not one: at these sizes the elongation (aspect ~1.5) is enough. A single-vortex reference needs a rounder shape (b/a ≳ 2–2.5) or a smaller cavity.
- Relative sense is not monotonic in neck width: counter (1.6), mostly counter (1.3), mostly co (1.1, both sizes), counter (1.02). One seed each and frequent reversals, so tentative; needs seeds and longer runs before any claim.

### 2026-10-01 15:45 — density check in cavities (k = 8 at the H3 radii; `runs/cavity/density_check/`, `analysis_spatial/density_check.txt`)
- Per-cell occupancy, interior cells (mean 20, Poisson SD 4.47): SD 6.6 / 9.0 / 9.7 at ζ = 0.07 / 0.09 / 0.11 (1.5 / 2.0 / 2.2 × Poisson); P(n < 10) = 4.5% / 15% / 17%; P(n < 5) = 0.13% / 1.4% / 1.9%; empty cells 0 / 7×10⁻⁵ / 7×10⁻⁵.
- **Strong radial stratification: the vortex centrifuges particles to the wall.** ⟨n⟩ near the wall (s < 5) 25 / 28 / 29, at the centre (s > 20) **14 / 9 / 8**; there P(n < 10) = 15% / 60% / 70%. The gap region (s = 3–5) is over-dense (25–28), not depleted, so low occupancy does **not** explain the closing gap at ζ = 0.11.
- −½ cores sit in depleted cells (⟨n⟩ = 15 / 11 / 11), +½ at ~21–23.
- Interpretation: compressibility (Mach number). Consistent with the 2026-09-29 box result (rms Ma 0.66 → 1.0 from ζ = 0.05 to 0.09): with v ~ 0.4–0.5 and c_s = 1, centrifugal compression ~ exp(∫v²/r dr) gives order-one density contrast. collOp 21 is a density-modulated dipole, so the active stress itself varies with this stratification.
- **Consequences:** (i) ζ ≥ 0.09 cavities have a strongly depleted centre (most central cells < 10 particles), affecting the −½ region and the bulk charge there; ζ = 0.09 looked clean in q_b but is already compromised in the centre. (ii) Even ζ = 0.07 has a 30% central depletion: state it in the methods; check whether H4/bulk-charge results depend on it (e.g. compare with a lower-Mach run: higher kBT or lower activity at the same R/ℓ_a). (iii) Supports ζ = 0.07 as the working activity and argues against using ζ ≥ 0.09 for H3.

### 2026-10-01 15:00 — H3 at ζ = 0.06 and 0.11 (R = 1.434 ℓ_a, 3 seeds × k = 2…10; `runs/cavity/H3high/`, `analysis_spatial/H3_all_activities.txt`, `enclosed_charge.txt`)
- ℓ_a (200², 2 seeds): 0.06 → 24.33 ± 0.09 (seeds 23.99, 24.39), 0.11 → 19.06 ± 0.09 (19.11, 18.80); R = 34.89, 27.34. Seed-to-seed spread 1.2–1.5σ, so single-seed ℓ_a errors are optimistic.
- Site rule, mean over seeds (k = 2…10): **ζ = 0.06**: q_b = 1.05–1.15 (k = 9: 1.06), P(+1) 0.67–0.83, Φ 0.85–0.95, 1.3–1.7 pairs/frame (between 0.05 and 0.07). **ζ = 0.07**: q_b 1.01–1.09, P(+1) 0.89–0.96. **ζ = 0.11**: k ≤ 6: q_b 1.0–1.2, P(+1) 0.6–0.8; **k ≥ 7: edge/bulk split ill-defined** (both rules erratic: q_b 0.5–0.9 site rule, 0.2–1.2 density cut), although circulation is strongest (Φ 0.98) and pairs most frequent (2.3/frame).
- **Site-rule fix:** at ζ = 0.11, k ≥ 7 the auto capture radius ran to the 8-cell window edge (no dip) and captured bulk defects (q_b ≈ −0.5). Now capped at 0.4 × median site spacing (`edge_sites.py --dfrac`, default 0.4); unchanged at ζ = 0.07 (cap never binds, max |Δq_b| = 0).
- **Cause: the layer–bulk gap does not scale with ℓ_a.** Radial +½ density at k = 8: the layer peak (s ≈ 1.2 cells) is the same at all ζ; the gap at s ≈ 3–4.5 is deep at ζ = 0.07 (5× below bulk), shallow at 0.06, almost absent at 0.11 (bulk +½ at full speed by s ≈ 4). Layer width is set in cells, bulk density ∝ 1/ℓ_a², so at high activity bulk defects press against the layer.
- **Cut-free diagnostic: enclosed charge Q(s)** (net charge further than s from the wall). ζ = 0.07: Q ≈ 1.0 flat for s = 3–5 (1.08 → 0.92), the "+1 bulk" is a plateau. ζ = 0.06: plateau at ~1.1–1.25. ζ = 0.11: no plateau, Q falls 1.25 → 0.75 over s = 3 → 5, crossing 1 only at s ≈ 4. Inside, Q falls below 0 towards the centre at all ζ (−½ concentrate at the centre, +½ ring at s ≈ 5–10).
- **ζ = 0.09 added** (H3pilot, R = 29.63, 1 seed; `analysis_spatial/act009_summary.txt`): clean like 0.07: 2k − 2 sites, edge occupancy 2k − 2 within 1%, q_b = 0.99–1.08, P(+1) 0.79–0.92, Φ 0.88–0.99, ~2.1 pairs/frame. Gap at s ≈ 3.8–4.2 is 1.2 vs bulk 4–5 (×10⁻³), i.e. ~4× (0.07: ~5×, 0.11: ~1.5×). Enclosed-charge slope between s = 3 and 5 at k = 8 (plateau flatness, per cell): 0.06: 0.13, **0.07: 0.07**, 0.09: 0.15, 0.11: 0.32 — flattest at 0.07, steepening monotonically above it.
- **H3 reading:** activity does not enter only through ℓ_a. At the onset (0.05) the layer is unstable; at high activity (0.11) the layer is intact but the gap that separates it from the bulk closes, because the layer has a fixed width in cells. The clean "+1 bulk" regime is where the gap is deep (ζ ≈ 0.07–0.09 at this R). Test: ζ = 0.11 at a larger R in cells (e.g. R = 31.8 like 0.07, or R/ℓ_a = 2) to see whether a plateau reappears; and Q(s) plateau flatness as the observable for the main sweep.

### 2026-10-01 14:00 — first active Cassini runs (P1, H = 48 done: b/a = 1.1, 1.3; k = 6, ζ = 0.07, 1 seed)
- Rendering: `runs/cassini/render_cassini.py` → `runs/cassini/movies/*.mp4` and `*_last.png` (vorticity background, director, +½ with orientation tail, −½ triangles). Lobe rotation measured as the mean azimuthal velocity about each half's centroid (the area-mean vorticity of a no-slip region is ~0 by Stokes, so it can't give the sense).
- **Each lobe hosts its own coherent vortex** (|v_φ|/|v| ≈ 0.45–0.6 in each lobe; whole-cavity v_φ about the centre ≈ 0), but the **sense reverses over time** (several sign changes in 4,000 time units).
- b/a = 1.1 (narrow neck, waist 30 cells): lobes **co-rotating in 75% of frames**, with an opposite-vorticity band in the neck. b/a = 1.3 (waist 49): mostly **counter-rotating** (same sign in 34% of frames; left mostly CCW, right mostly CW). Tentative (1 seed each); compare with Wioland et al. 2016 (gap-dependent co-/counter-rotation in coupled bacterial vortices, to verify).
- Imposed defects stay at the wall (12 +½ for k = 6). Charge per lobe / site rule not yet analysed.

### 2026-10-01 — Cassini walls implemented and validated (branch `cassini`, commits 362ca8b, e63155b; `runs/cassini/T1`, `T2`)
- Implementation as in `reports/cassini_patch_plan.md`: `"shape": 1`, `"cassiniA"`, W = d1·d2 − R², normal = ∇(d1·d2), secant crossing times, eMode 3 (arclength by polar angle), setup guards, synopsis + checkpoint. Separate fix commit: checkpoints now store customAnch/eMode (they didn't), domain walls get CA/EMODE defaults (were uninitialised after realloc). Checkpoint restarts segfault in both old and new builds (pre-existing; TODO).
- Passive k = 3 tests (`cassini_look.py`, figure `runs/cassini/tests_T1T2.png`): a = 0 vs circle R = 15: alignment 0.970/0.961 (eMode 0/3) vs 0.971, accessible volume 706.8 vs 706.9, audit 100%, 6 +½ per frame. b/a = 1.6 and 1.2 (concave neck), eMode 3: alignment 0.976, audit 100%, volume within 0.1% of the oval area. eMode 0 at b/a = 1.2 works too (0.964, audit 99.7%) but bunches the imposed defects towards the tips; eMode 3 spreads them along the wall including the neck.

### 2026-09-30 14:10 — act 0.05, 3 seeds (`runs/cavity/analysis_spatial/act005_seedavg.txt`, `act005_sites.csv`; `analysis_lunch/circ_act005_3seeds.csv`)
- Seeds 2–3 (`H3pilot_seeds/act0.05/`, fixed seeds, legacy tracker) done 13:46, no failures. R = 37.51 (R/ℓ_a = 1.43), k = 2…10.
- **The pilot's picture holds with 3 seeds.** Site rule: 2k − 2 sites for k ≤ 8 in every seed, but **occupancy below 2k − 2** (e.g. k = 6: 8.99 ± 0.19 of 10; k = 8: 12.4 of 14). **q_b = 1.16–1.81**, rising with k up to k = 8, and P(q_b = +1) only 0.36–0.68 (vs 0.9 at ζ = 0.07). Circulation Φ_pers = 0.83–0.92 for k ≤ 8.
- **k = 9, 10: circulation breaks down in all seeds** (Φ_pers = 0.37 ± 0.21, 0.51 ± 0.10), and site detection picks up spurious interior maxima (16–27 sites), so q_b there (≈ 1.15) is not reliable.
- Dynamical pairs 0.6–1.3 per frame, no pair in 22–50% of frames (ζ = 0.07: ~1.9, ~5%).
- Bulk pair correlations (report Fig. 7, now 27 runs at ζ = 0.05): **+½/+½ collapses** with ζ = 0.07 and 0.09 in r/R; **−½/−½ collapses only up to ~0.5 R** — beyond that ζ = 0.05 stays at g ≈ 0.75–0.9 (narrow errors), i.e. far −½ pairs are rarer. (The 1-seed version looked like a full collapse; corrected.)
- Reading: at the onset, the +½ arrangement still scales with ℓ_a, the −½ one differs at large separation, and the wall layer is unstable (vacant sites, excess bulk charge), and at high k the vortex itself fails. Relevant for H4: at ζ = 0.05, k = 9–10 are cases where circulation switches off at fixed R/ℓ_a.

### 2026-09-30 — pinned-site edge classification + seed-pooled g(r | s) (`runs/cavity/edge_sites.py`, `analysis_spatial/rscan_sites*`, `ref_sites*`)
- **Rule:** sites = local maxima (> 0.02 per cell² per frame) of the time-averaged, σ = 1-smoothed +½ density map; no radial cut, no circle geometry. A +½ is "edge" if within d_site of a site; d_site = minimum of that run's distance-to-nearest-site histogram (1.5–8 cells). Fixed d_site = 2–3 is too tight (passive k=5: only 7.8–9.6 of 10 captured).
- **R-scan, act 0.07, 3 seeds** (`rscan_sites_seedavg.txt`): **sites = 2k − 2 in all 27 runs**; d_site shrinks with k (6.2 → 3.1 cells, ≈ half the site spacing); N_edge = 2k − 2 to within 1%; **q_b = 1.01–1.05** (k = 10: 1.09 ± 0.06, one seed with a transient), **P(q_b = +1) = 0.89–0.96**, equal to or above the density cut; the two rules agree on 96–99.7% of +½ defects (lowest at k = 2).
- **Passive:** k ≥ 4: 2k sites, all defects edge, q_b = 0. k = 2, 3: static defects sit *away* from the wall (k = 2: 5 sites at s ≈ 5.7) — the rule finds pinned defects, not wall defects; fine for q_b (0) but "pinned" ≠ "at the wall" in passive runs.
- **act 0.09:** 2k − 2 sites, q_b = 0.99–1.08, like act 0.07. **act 0.05 (1 seed):** 2k − 2 sites (except k = 9: 27, flow reversals), but **sites are often vacant** (N_edge 5–10% below 2k − 2) and q_b = 1.1–1.6: at the onset the layer is less stable and exchanges defects with the bulk. Waiting for seeds 2–3.
- **Pooled g(r | zone), all k and seeds** (`rscan_sites_gr_zones.png`): bulk (outer and centre alike): **+½/−½ peak g ≈ 3.2 at r ≈ 2.5**, tail > 1 out to r ≈ 20–25; **+½/+½ and −½/−½ have the same correlation hole** (g ≈ 0 for r < 3, 0.5 at r ≈ 6, 1 at r ≈ 15–20, i.e. ~ℓ_a); the hole is slightly longer in the centre. Edge ring: +½/+½ lattice peak at r ≈ 11 (site spacing, pooled over k), +½/−½ ≈ 1 (layer defects don't bind −½). Per-k maps (`rscan_sites_gr_pooled.png`) are consistent across k; the red diagonal (partner in the far-side layer) is still an artefact of the thin layer peak in the expectation (not fixed by finer θ sampling; next: finer radial bins for n_b near the wall).

### 2026-09-30 — spatial distribution of defects, edge/bulk from pinning (`runs/cavity/analysis_spatial/`, `defect_spatial.py`)
- User's idea: layer defects are stationary → sharp peaks in the time-averaged density; bulk defects move → flat/noisy. Runs: passive k=5, ζ=0.07 k=2/5/10 (seed 2), ζ=0.05 and 0.09 k=5 (R = 1.43 ℓ_a).
- **Confirmed.** Time-averaged maps: layer +½ sit at **fixed angular sites** (discrete spots, no ring smearing, although the bulk flow circulates). Site count from the angle histogram of +½ with s < 4.2: **2k − 2 in every active run** (2, 8, 18 for k = 2, 5, 10; also at ζ = 0.05 and 0.09), **2k in the passive run** (10).
- Radial density n_+(s): sharp peak at s ≈ 1–2, a gap at s ≈ 3–5, then a flat noisy bulk plateau (~10⁻³ per cell² per frame); n_− rises towards the centre and is ≈ 0 for s < 5. Speed: layer +½ move at 0.2–0.35 cells/time, the same as the passive defects (0.3 = tracker jitter floor, passive defects don't move), bulk +½ at 0.8–1.0, peaking at s ≈ 8–12 and decreasing towards the centre. The speed step and the density gap both sit at s ≈ 3–5, i.e. at d_ρ.
- **Proposed edge definition:** a defect is "edge" if it sits on one of the fixed sites of the time-averaged map (pinned), bulk otherwise. Needs no radial cut and carries over to ellipses (no circle geometry). To do: per-defect classification (distance to nearest site, or per-track mobility), compare with d_ρ on all R-scan runs.
- Pair correlations g_ab(r) (reference = pairs between frames i and i + nf/2, keeps the one-body profile): only meaningful for mobile (bulk) defects; for the static layer the lag reference contains self-pairs (passive g shows the site lattice). Bulk, s > 4.2, preliminary (1 seed, few −½): +½/+½ depleted below r ≈ 10–15, −½/−½ depleted below r ≈ 5, +½/−½ enhanced (g ≈ 2–3.5) at r ≲ 3–5 (fresh or annihilating pairs), strongest at ζ = 0.05. Needs pooling over seeds and k.
- **Ring-resolved versions (user's suggestion, rings of width 1):** (1) same-ring angular g(s, Δφ) (`*_ring_gphi.png`): clean for the layer (peaks at the site spacing Δφ = 2π/(2k−2)), but same-ring bulk pairs are too rare (−½/−½: no ring with ≥ 30 pairs). (2) **g(r | s)**: reference defect binned by ring, partner anywhere, expectation = measured radial density n_b integrated over the annulus around the reference (`*_ring_gr.png`); within-frame only, no self-pairs, radial profile divided out. Readable at 1 seed: +½/+½ exclusion zone r ≲ 5–8 at every ring, layer lattice stripes for s < 4; +½/−½ and −½/+½ (same pairs, consistency check) enhanced at r ≲ 10 for bulk references; −½/−½ depleted at r ≲ 5–10 (noisy). The red diagonal at the largest r (partner in the far-side layer) is a discretisation artefact of the thin layer peak in the expectation.

### 2026-09-30 — lunch runs: R-scan seeds 2–3, H3 pilot, passive references (`runs/cavity/analysis_lunch/`)
- All legacy tracker, audit q = k in ≥ 98% of frames. Seed 1 (`Rscan_R31.83/`) has defectsOut/50 → its docking rates are not comparable (linking warning); seeds 2–3 have defectsOut/20.
- **R-scan R = 31.83 (R/ℓ_a = 1.43), act 0.07, 3 seeds** (`rscan_seedavg.txt`): with the density edge d_ρ ≈ 4 cells, **q_b = 1.00 ± 0.02 for every k = 2…10**, P(q_b = 1) ≈ 0.9, and the layer holds N_L = 2k − 2 +½ defects (charge k − 1). Seed-to-seed scatter is small. With the passive edge d_p (which shrinks from 8.3 to 2.3 cells with k), q_b runs 0.44 → 1.20: the passive layer is thinner than the active one, so d_p is the wrong cut for k ≳ 6 (and too wide at k = 2). **Use d_ρ.**
- **Passive references:** k ≥ 3: all 2k +½ defects in the wall layer, q_b = 0, no fluctuations. k = 2: no layer (defects in the bulk; q_b(ρ) = 1.95).
- **H3 pilot (1 seed, R = 1.434 ℓ_a):** act 0.09 (R = 29.6) matches act 0.07: q_b = 1.00–1.06, thin layer (fl ≈ 0.03). **act 0.05 (R = 37.5) does not:** q_b(ρ) = 0.83–1.87, rising with k, wider/fluctuating layer (fl ≈ 0.1, d_ρ 5–8 cells), and circulation breaks down at k = 9, 10 (Φ_pers = 0.02 / 0.60). Tentatively **against H3** (activity doesn't enter only through ℓ_a), but: 1 seed, and the density edge is less well defined at act 0.05. Needs 3 seeds before quoting.
- **Dynamical defects vs activity at fixed R/ℓ_a = 1.43** (`pairs_by_activity.txt`): mean ±½ pairs per frame ≈ 1.2 at act 0.05 (0 pairs in 15–45% of frames), ≈ 1.9 at act 0.07 (0 pairs in ~5%), ≈ 2.1 at act 0.09 (~4%). A bulk defect count πR²/ℓ_a² ≈ 6.5 would be expected at every activity. The user's reading (2026-09-30): act 0.05 is at the lower edge of the active-turbulent regime (barely active); a student's run also showed very few dynamical defects there. Consistent with the bulk data: ⟨N⟩ drops fast below act 0.06, τ_int is large at 0.05 (69 vs ~25), and ℓ_a(0.05) differs between L = 100 and 200 (27.5 vs 26.2). **So the act 0.05 difference is more likely proximity to the onset than a failure of ℓ_a scaling:** H3 should be tested only well inside the turbulent regime (act ≳ 0.06), with act 0.05 as the onset limit. Seeds 2–3 at act 0.05 launched 2026-09-30 11:32 (office, `H3pilot_seeds/act0.05/s0{2,3}/`, queue `queue_h3_seeds.log`, fixed seeds 81500+k / 82500+k, legacy tracker, same inputs otherwise).
- circulation.py with the fixed `--wall 2.5` cut underestimates the layer (see next step 6). Φ ≈ 0.9–0.98 for all act 0.07/0.09 runs.

### 2026-09-28 — passive baseline (act = 0), R = 15 (`runs/cavityValidation/passive`)
- The custom anchoring is imposed correctly: alignment 0.96–0.99; `customAnch` = 0 and 1 give identical results.
- The tracker is exact: always 2k +½ defects, net charge = k in 100% of frames, no spurious wall defects.
- k = 1: two +½ defects at r ≈ 9–10. **k ≥ 2: all defects in a wall layer ~2–3 cells thick, bulk charge = 0** (Q(r) = 0 for r ≤ 11).

### 2026-09-28 — passive mirror check, customAnch = −2, −3, R = 15 (`runs/cavityValidation/passive/km2, km3`)
- Exact mirror of k = +2, +3: Q(R−1) = −2.00 / −3.00, alignment 0.99, the tracker finds 2|k| −½ defects (no +½), net charge = k in 100% of frames.
- Passive bulk charge 0 again (all defects in the wall layer) → **any positive/negative asymmetry in active runs is an active effect** (H5).

### 2026-09-28 — active reproduction, k = 5, ρ = 0.1, R = 15.9, act 0.07 (`runs/cavityValidation/active/k5_rho0.10`)
- **Bulk (r < R − 2.5) net charge +1 in 94% of frames**, carried by two +½ defects at r ≈ 9 (IQR 7.8–10.5); the other +4 sit at the wall. Stable (0.99 in both halves of the run).
- Q(r) rises smoothly to ≈ 1 by r ≈ 13, then jumps to ≈ 5 at the wall → sharp edge/bulk split at r ≈ R − 2 to R − 3.
- **Compared with passive (bulk 0): the +1 bulk comes from activity.**
- Candidate observable: P(bulk charge = +1) against (R/ℓ_a, ρ·ℓ_a).

### 2026-09-29 — arclength-uniform anchoring (eMode 2) validated (`runs/emode2Validation/`)
- Circle R = 15, k = 3: eMode 2 ≡ eMode 0 (alignment 0.984 vs 0.985, 6 defects, Q = 3, audit 100%).
- Ellipse a = 20, b = 10 (tip radius 5), k = 4, passive: fraction of the director winding in the tip regions |x| > a/√2, measured 0.35 cells inside the wall (total winding 3.93–4.00): eMode 0 0.694 (pred. 0.709, normal-angle), eMode 1 0.478 (pred. 0.500, parameter), **eMode 2 0.384 (pred. 0.396, arclength)**. eMode 0 fails the audit in 14% of frames (tip density ≈ 0.25 > 0.15) — independent confirmation of the local resolution limit.

### 2026-09-29 — compressibility vs activity (`runs/activeLength/compressibility/compressibility_results.txt`)
- Density SD 6.1 → 9.8 (Poisson 4.47) and rms Mach 0.66 → 1.0 (isothermal c_s = 1) from act 0.05 to 0.09: **large compressibility at every activity used, rising smoothly, no threshold.**
- Ma_rms = 1.00 exactly at 0.085, where ℓ_a(ζ) flattens: suggestive (depends on the choice of c_s), not an independent cutoff.
- ℓ_a ∝ ζ^(−½) up to ζ = 0.08 (α = 0.485 ± 0.031, both boxes); full range α = 0.39, single power law rejected. Cavity results use the measured ℓ_a(0.07) = 22.2 ± 0.3, so the exponent doesn't enter them; for H3, measure ℓ_a at each activity used.

### 2026-09-29 — resolution limit of the custom anchoring (`runs/cavity/resolution_scan/resolution_results.txt`)
- Set by the imposed spacing 1/ρ, not by R; **passive fails like active** → a method limit (cell-scale anchoring).
- ρ ≤ 0.10 clean; **ρ = 0.15 good**; ρ = 0.20 degraded (alignment ~0.78, up to 20% frames missing imposed defects); **ρ ≥ 0.25 fails** (audit ≤ 30%).
- → Physically meaningful range **ρ ≲ 0.15** (spacing ≳ 6–7 cells). Batch 1's ρ ≥ 0.25 part documents the limit only.
- Figure: `runs/activeLength/active_length_vs_activity.png/.pdf` (ℓ_a vs activity, both boxes, fit) made the same day.

### 2026-09-29 04:14 — ℓ_a finite-size check, L = 200 vs 100 (`runs/activeLength/L200/finite_size_results.txt`)
- ℓ_a(L=200) = 26.16 ± 0.23 (act 0.05), **22.21 ± 0.26 (0.07)**, 20.66 ± 0.12 (0.09). Differences from L=100: −4.8%, +4.0%, +2.3% (1.1–1.6σ, no consistent sign) → **no significant finite-size effect** (box/ℓ_a 4–5 vs 8–10).
- **Reference ℓ_a(0.07) = 22.2 ± 0.3** (L=200; weighted mean of both boxes 22.0 ± 0.2). ~4% above the L=100 value; no conclusion changes.
- α(0.05→0.09) = 0.40 ± 0.02 at L=200 (L=100 fit 0.38 ± 0.04): the departure from ζ^(−½) is real, not finite size. Local slope still flattens: 0.49 (0.05→0.07) vs 0.29 (0.07→0.09).

### 2026-09-28 19:20 — fixed-R scan R = 31.83 (R/ℓ_a = 1.5), k = 2…10, ρ = 0.02…0.10 (`runs/cavity/Rscan_R31.83/rscan_results.txt`)
- Audit clean (99–100% of frames q = k). **Strong steady vortex in every run** (Φ_persist = Φ_inst = 0.89–0.98). ~1.3–2.2 dynamical ±½ pairs on average (vs ~0 at R = 15.9).
- **Bulk charge ≈ 1.0–1.2 for all ρ** with bulk = r < R − 4; P(q_b = 1) ≈ 0.74–0.95. With the earlier cut R − 2.5 it looked like a low-ρ rise (1.55 at ρ = 0.02) → **artefact**: the wall layer is diffuse and thicker at low k (r ≈ 28–30 for k = 2) and sharp at high k (1–2 cells from the wall for k = 10).
- Structure: −½ near the centre, +½ ring at r ≈ 20–26 (Q(r) < 0 at intermediate r), Q ≈ 1 at r ≈ 28 for all k.
- → **H2 not supported here; H4 consistent but untested** (vortex and bulk +1 both "on" everywhere). Edge/bulk cut must be data-driven (density minimum / Q plateau), not a fixed d.

### 2026-09-28 — circulation tool test, k = 5, ρ = 0.1, R = 15.9, act 0.07, 3 short seeds (5k warm-up + 10k steps; scratch)
- Strong circulation: Φ_inst = 0.84 / 0.71 / 0.91; bulk charge mostly +1 (P = 0.74 / 0.74 / 0.96). Consistent with H4 (one parameter set only).
- Seed 2: Φ_persist = 0.16 vs Φ_inst = 0.71 → **the vortex reversed within the run**. Seed-averaging |⟨v_φ⟩| alone would understate circulation → report both.

### 2026-09-28 18:16 — ℓ_a, L = 100, full data (`runs/activeLength/L100/active_length_results.txt`)
- ℓ_a = 27.5 ± 1.1 (act 0.050), 25.3 (0.055), 22.9 (0.060, 0.065), 21.4 (0.070–0.080), 20.6 (0.085), 20.2 ± 0.4 (0.090).
- Reference for the cavity work at the time: ℓ_a(0.07) = 21.36 ± 0.47 (box 100) — superseded by the L=200 value 22.2 ± 0.3.
- Power law: α = 0.38 ± 0.04 (χ²/dof 12.1/7); α = ½ disfavoured (20.6/8). Flattens at act ≥ 0.07. Possible causes: finite size (box/ℓ_a ≈ 5), a minimum set by the defect core / cell size, or underestimated errors (one seed, n_eff 19–75).
- Counting: no split defects (merged = raw); neutral-frame filter keeps 62–78% and biases ⟨N⟩ down (8/9 activities), within errors.

### 2026-09-28 — ℓ_a preliminary (L = 100, first ~300 frames; superseded)
- ℓ_a ≈ 33 (act 0.05) → ≈ 19 (act 0.09); not yet converged. Box/ℓ_a ≈ 3–5 → finite-size check queued.
- The defect-splitting artefact isn't seen (merged = raw count). The neutral-frame filter is biased, as expected (edge-strip losses).

### 2026-10-11 — What drives the inward flow burst at a depleting escape? (selfflow.py)
Script: reports/master/theory/energetics/selfflow.py; figure fig_selfflow.png. Data: fineflow (flow + defects every 2 t.u.), 3 ellipse runs, depleting escapes (circle run: 0 events), 668 aligned rows.
- Decomposed the inward fluid velocity at the escaper into core (r<2) and ring (4<r<8). The ring stays flat or slightly OUTWARD (−0.01 to −0.04) through the event, so the burst is localized on the defect, not an external stream.
- Own-jet calibration: for bulk +1/2, (core − ring) flow projected on the head = 0.32 cells/t.u. (median; n 1033). (A first version measured defect velocity relative to the core fluid: −0.12. That doesn't measure self-propulsion, because the core fluid already contains the defect's own jet and the defect lags it. Discarded.)
- Prediction u_self·cos τ against the measured core − ring excess:

  | lag (t.u.) | prediction | measured | cos τ |
  |---|---|---|---|
  | −20 | 0.12 | 0.07 | 0.37 |
  | −10 | 0.19 | 0.13 | 0.59 |
  | −2 | 0.23 | 0.22 | 0.73 |
  | 0 | 0.24 | 0.25 | 0.75 |
  | +2 | 0.24 | 0.34 | 0.75 |

  Near the wall the jet runs about 0.05 weaker than the bulk value (consistent with no-slip damping at depth ≈ 2). It converges at departure. The defect's own inward velocity follows the core fluid.
- READING: the inward burst is the escaper's own active jet. Untilting rotates it inward, and it then carries the defect out. The order is untilt (from −40) → jet turns inward → departure. There's no evidence of an external flow pushing it. This resolves "untilting first, then flow": they are one mechanism. What remains open is what triggers the untilting. The ring flow and the elastic barrier are unchanged before the event (energetics), so it may be fluctuations of the tilt mode or local vortex incoherence.
- Tension with the bulk-only activity runs (wall inactive, layer still empties): there the wall defects have no own jet. The depletion route must then be different (bulk-flow driven). To check: run selfflow on actlayout bulk-only runs if they have fine output (they don't; they would need a rerun at fine output).

### 2026-10-11 — Escape movies (energetics/movie/escape_movie.py)
- escape_eps2.0_s2_ev4.mp4 (representative), ev12.mp4 and ev6.mp4 (atypical: head starts into the wall, cos τ = −0.97). Each covers −60 to +24 t.u. at 4 fps (frames every 2 t.u.). Shown: vorticity (smoothed σ = 1 cell), flow arrows, defect comets, escaper trail and head, inward normal, circles at r = 2/4/8, and insets of cos τ and the inward core/ring flow. Frames are kept in frames_ev4/ and frames_ev12/.
- Per-event survey (survey.py): single events are heterogeneous. Some escapers sit untilted (cos τ ≈ 0.9) for the whole 60 t.u. and others turn from head-into-wall within 20 t.u. The ensemble rise in cos τ is a mean over these. In ev4 the core jet appears only once cos τ > 0.5, from about −5 t.u. After that it reaches 0.3–0.7 while the ring flow turns outward (−0.12): the defect's own jet, as selfflow.py concluded.

### 2026-10-11 — Lead/lag of near-wall flow vs wall-defect tilt, fineflow (tiltprotect/lead_lag_fineflow.txt)
Frame spacing 2 t.u.; lags ±30 frames (±60 t.u.); lag > 0 means the flow leads.
- Magnitude pairs (|v_t| vs |sin τ|, |v_t| vs alignment): in all three ellipse runs the peak is at lag 0 to +3 frames, i.e. the flow leads by 0–6 t.u. at most. Peak 0.17–0.37, half-width about ±6 frames (±12 t.u.). In the circle (eps1.0) the correlation is weak and noisy (peak 0.2 at lag 0, negative side-lobes).
- The signed pair (v_t vs τ) is broad and flat (0.1–0.25 across ±60 t.u.) with no sharp peak. It reflects the slow persistent vortex sense and carries no timing information.
- READING: flow and tilt change essentially together, on the scale of the tilt relaxation. At most the flow leads by about one frame spacing. This fits the flow-alignment picture (tilt slaved to the local shear within a few t.u.). It does not show what triggers the untilting before a depleting escape. That needs an event-conditioned test (does local wall shear drop before cos τ rises?). Note that the ring inward flow in selfflow stays flat before departure.

### 2026-10-11 — Event-conditioned local wall flow before depleting escapes (energetics/shear_trigger.py, shear_trigger.txt)
Fineflow, 3 ellipse runs. v_t = flow along the wall at the escaper's wall point (depth 1–5, ±4 cells along).
- First version: each series signed by its own pre-event mean. Both the escaper and control curves decayed to zero, which is the regression-to-mean artefact. Discarded.
- Corrected version: signed by the instantaneous cavity-wide wall circulation. |global| is only about 0.04 in these elongated cavities.
  - Control wall +½ sites follow the global sense: +0.02 to +0.06.
  - Escaper sites show about 0 (±0.02) from −60 to −24 t.u., before the untilting is under way. From about −18 t.u. they turn slightly counter to the global sense, to −0.04 ± 0.02, at the same time as cos τ rises.
- READING (weak, ~1.5–2σ): escapers sit where the local wall flow does not follow the cavity circulation (a stagnation or counter-flow patch), so there is no shear to hold the tilt. This is consistent with "incoherent flow removes the tilt protection". The late counter-flow may partly be the escaper's own jet as it rotates. It needs the 12-run fineflow + fineflow2 set before any claim.

### 2026-10-11 (evening) — Shape-sweep energetics of depleting escapes (energetics/plot_etsg.py -> fig_etsg_shapes.png, etsg_shapes.txt)
Shape sweep (ellipses ε 1.2–2.5 + circles; R_eq 1.0/1.5/2.0; 3 seeds; defects every 2 t.u.); energy and tilt only (flow every 100 t.u.). Events at departure: circle 22, ε1.2 143, ε1.4 318, ε1.6 306, ε2.0 312, ε2.5 335.
- Untilting is universal. At −40 t.u. escapers have the same cos τ as residents in every shape (Δ ≈ 0). They then rise to about 0.75 at departure, while residents stay flat. The onset is around −40 t.u. at every aspect and size.
- Resident cos τ falls with aspect: 0.53 circle, 0.45 ε1.2, 0.41 ε1.4, 0.35 ε1.6, 0.32 ε2.0, 0.25 ε2.5. The wall layer becomes less ordered in tilt as the cavity elongates (cos τ can't tell a larger tilt from mixed or wall-pointing heads).
- Elastic restoring force: no softening before departure at any aspect. It stays at −11 to −14 kBT/cell and becomes more restoring as the defect moves out.
- Cost to move 6 cells in (E6): escapers are 4–5 kBT cheaper than residents at −40 t.u. for ε ≥ 1.4, and about 10 kBT cheaper at ε 2.5. The gap grows with aspect and holds at all sizes (size split in the txt). There is no such gap in the circle or at ε 1.2. So there is a mild elastic predisposition of escape sites in elongated cavities (about 10–15 % of the barrier). It is set long before the event, not a precursor. The control is matched in |x|/A (within 0.15), so this is not just tips vs flanks. It could be local layer spacing (a vacancy or gap next to the site).
- The last 10 t.u. drop in E6 (≈ 8–12 kBT) is the defect moving off the wall, as before.

### 2026-10-11 (evening) — Identities with the union-find tracker (kickout/uftrack_identity.py -> uftrack_identity.txt)
- runs/uftrack (R_eq 1.5; circle k 6, ε 1.6/2.0 k 7; 2 seeds; defects every 2 t.u.): the identity life is still a median of 2–3 frames (90th percentile 11 circle, 47 ellipse). The flicker is at the detection level (defects drop out and reappear frame to frame), not in the tracker. So the identity classification of nearby bulk +1/2 stays unresolvable, as with the legacy tracker. This line is dropped. Chain-following with gap bridging (energetics.events, selfflow) remains the working approach.
- Side note: 27–40 % of escapers were in the bulk within the previous 40 t.u. (random wall defect: 8–18 %). This is the revolving door again, which the depleting-escape filter excludes.

### 2026-10-11 (night) — 12-run fine-flow analyses (fineflow + fineflow2: circle, ε 1.4, 1.6, 2.0; selfflow_ff12, shear_trigger_ff12, ets_ff12)
- selfflow_ff12: the same picture with about 2x the events. Bulk +1/2 own jet 0.317 (n 3051). The core − ring excess at the escaper is 0.07/0.12/0.17/0.24 at −20/−10/−4/0 t.u., while 0.32·cos τ = 0.13/0.19/0.22/0.23. The ring stays at about −0.01 to −0.02 throughout.
- shear_trigger_ff12:
  - Controls: wall flow +0.00 to +0.03 along the circulation.
  - Escaper sites, −60 to −24 t.u.: −0.01 to +0.005. That is only 0.02–0.03 below the controls (~1.5σ), so a pre-existing stagnation patch stays weak.
  - From −18 t.u.: counter to the circulation, −0.05 ± 0.012 at −12…0 (4σ). This coincides with the untilting; it is likely the tangential part of the escaper's own jet as it rotates.

### 2026-10-11 (night) — Force model fit for wall +1/2 (energetics/forcefit.py -> forcefit_ff12.txt)
Model: v_in = α·(u_core − u_ring) + β·u_ring + μ_n·F0 + c, along the inward normal, per frame. Data: 1486 escaper samples (lags −40…+6 t.u.), 211 residents, 1107 bulk +1/2 (along the head). Errors from a bootstrap over events.
- Bulk +1/2: α = 0.63 ± 0.05, β = 0.54 ± 0.02 (R² 0.45). A free +1/2 moves at about 60 % of its own jet plus about half the background flow. Mean speed along the head 0.14 cells/t.u.
- Escapers: α = 0.85 ± 0.04, β = 0.97 ± 0.08, μ_n = 0.0022 ± 0.0008, c = +0.02 ± 0.01 (R² 0.31). With residents added: μ_n = 0.0016 ± 0.0008.
  - The normal mobility of a wall defect is 0.002 ± 0.001 cells²/(kBT t.u.), about 2.5x below the free-defect value (0.0054) and 10x below the along-wall value from layer fluctuations (0.02–0.03).
  - So the elastic pull, μ_n·F0 with F0 ≈ −12, is only about −0.03 cells/t.u.
  - Caveat: F0 varies little and is noisy (errors-in-variables bias pulls μ_n down), and it trades off against the constant. Treat it as an order-of-magnitude bound: μ_n ≲ 0.003.
- Budget by lag (escapers; measured v vs model = jet + ring + elastic):

  | lag (t.u.) | measured v | model | α·jet | β·ring | μ·F0 |
  |---|---|---|---|---|---|
  | −30 | +0.02 | +0.03 | +0.04 | 0 | −0.03 |
  | −20 | 0 | +0.03 | +0.06 | −0.02 | −0.03 |
  | −10 | +0.08 | +0.07 | +0.10 | −0.02 | −0.03 |
  | −4 | +0.16 | +0.12 | +0.15 | −0.02 | −0.03 |
  | 0 | +0.34 | +0.19 | +0.20 | −0.01 | −0.03 |
  | +4 | +0.26 | +0.28 | | | |

  The model tracks the run-up. It misses the burst at departure (0.34–0.40 measured vs 0.19–0.25), probably detection jumps as the defect crosses the band edge.
- Residents: v ≈ 0, own excess flow only 0.055 (vs 0.32·cos τ ≈ 0.12–0.15 expected from their tilt). The jet of a wall defect is suppressed by more than half at depth ≈ 2, on top of being turned along the wall by the tilt. The fit gives c = −0.06 ± 0.03, a weak hint of extra pinning (static friction).
- Jet modelled from the tilt alone (0.32 cos τ): R² 0.10 vs 0.30 with the measured jet. The jet near the wall is not just tilt × bulk strength; wall suppression (depth) matters.
- READING (force picture): the elastic restoring force is real (−12 kBT/cell) but, with the small normal mobility, gives only about −0.03 cells/t.u. of drift toward the wall. A resident's suppressed, tangentially turned jet gives about +0.05 inward, roughly balanced by the pull and pinning (net 0). Untilting plus moving out of the suppression zone raises the inward jet to 0.15–0.20 by −4…0 t.u., 5–7x the elastic pull, and the defect goes. The hold of the layer is thus mainly geometric (tilt plus wall suppression of the jet), with elasticity a weak backstop.

### 2026-10-11 (night) — Own-jet vs depth and tilt: the wall blocks the NORMAL component, not the strength (energetics/jetmap.py -> jetmap_ff12.txt, fig_jetmap_ff12.png)
38,565 +1/2 samples from the 12 fine-flow runs. Jet = (flow within 2 cells) − (4–8 cell ring, fluid cells only).
- Jet strength along the head is NOT reduced near the wall: 0.36–0.40 at depth 0–3 vs 0.32 in the bulk (×1.15–1.27). It is weaker only for heads pointing into the wall (cos τ < −0.6: 0.26–0.29).
- The inward component is blocked. g(d) = transmitted fraction of the normal component (inward jet vs head-projection jet·cos τ): 0.30 (d < 1), 0.33 (1–2), 0.58 (2–3), 0.86 (3–4), ≈ 1 beyond 4 cells. This is no-penetration at the wall; the jet is deflected along the wall.
- This corrects the earlier note "jet suppressed by more than half at depth ≈ 2" (forcefit entry): it is the normal projection that is suppressed, ×0.33.
- Closed law: inward jet = 0.40·cos τ·g(d).
  - It predicts the residents' measured inward jet exactly: 0.056 vs 0.055 at depth 2.1, cos τ 0.34.
  - Force fit with the closed jet instead of the measured one (no flow data except the ring): R² 0.22 (vs 0.30), α = 0.79 ± 0.07, μ_n = 0.0027 ± 0.0010. Run-up reproduced: −20/−10/−4 t.u. model 0.04/0.08/0.14 vs measured 0.00/0.08/0.16.
- READING: a hydrodynamic trap with a runaway.
  - In the band (d ≈ 2) only a third of the normal jet acts. Even an untilted defect gets about 0.4·0.75·0.33 ≈ 0.10 inward against about −0.03 elastic drift (μ_n·F0), so it creeps out.
  - Each cell of depth then raises g (0.33 → 0.58 → 0.86 → 1). The push grows 3x between 2 and 4 cells, so the escape accelerates. This matches the depth curve: flat until about −12 t.u., then accelerating.
  - Tilt sets the projection, the wall sets the transmission, and elasticity is a weak backstop.

### 2026-10-11 (night) — Vorticity around wall +1/2 and escapers (energetics/vortdef.py -> vortdef_ff12.txt, fig_vortdef_ff12.png)
User's observation from the movies: escapers end up at about zero vorticity once detached (riding their own flow). Before that they sit in signed vorticity, as if the wall suppressed the lobe on its side; the vorticity ramps out as they detach.
Method: 12 fine-flow runs. Vorticity smoothed σ = 1 cell. Core and sides signed by the cavity circulation (> 0 co-rotating). The +1/2 lobes (2.5 cells either side of the head axis) are measured raw in the defect frame, in units of the free bulk +1/2 lobe strength O0 = 0.118 1/t.u.
- All +1/2 by depth:
  - The wall-side lobe is ABSENT in the band: −0.43/−0.20/+0.03 of the free value at depth 0–1/1–2/2–3. It recovers to 0.35 (3–4), 0.48 (4–5), 0.73 (5–6), 0.86–0.89 (6–10). This is the same depth scale as the jet transmission g(d) (0.33/0.58/0.86/1).
  - The core vorticity is then the bulk-side lobe plus the boundary layer: −0.084 (counter-rotating) at depth 1–2, −0.057 at 2–3, about 0 by 5–6, +0.035 (co-rotating) in the bulk.
  - The bulk-side lobe is ×1.6–1.9 in the band. Caveat: the tilt sense is locked to the circulation, so the counter-rotating boundary-layer vorticity adds to this lobe; part of the enhancement is background.
- Escapers, aligned on departure:
  - The core goes from −0.044 (−40 t.u.) to −0.002 ± 0.010 at −10…−8 t.u. while the depth is still about 2.3. Residents at −10: −0.024 ± 0.015. It returns to −0.02…−0.03 at departure and after (depth 3.5–5.6, background level for that depth).
  - The wall-side lobe grows before the defect moves: 0.07 (−40) → 0.10 (−20) → 0.22 (−10) → 0.39 (−4) → 0.49 (0) → 0.8 (+2…+6). Residents stay at −0.05…+0.13. At −10 the escaper's wall lobe (0.22) is well above that of all +1/2 at the same depth (0.03), so the recovery starts with the untilting, before detachment. (Part of it is geometric: an untilted defect's lateral lobe sits farther from the wall.)
- READING: the user's picture holds and joins the jet result. The jet of a +1/2 is the flow between its two counter-rotating lobes. Near the wall the wall-side lobe is erased (no-slip), the dipole is one-sided, the defect sits in net counter-rotating vorticity, and its jet is deflected along the wall (g ≈ 0.3). Untilting swings the lobes so that the wall-side lobe moves off the wall and recovers. The dipole re-symmetrises (core vorticity → 0, "riding its own flow"), the jet straightens inward, and the defect leaves. The vorticity is the flow-field signature of the same transmission factor g(d).

### 2026-10-11 (night) — Q criterion around defects (Head, Doré, …, Shendruk, Nat. Phys. 20, 492 (2024)) (energetics/qcrit.py -> qcrit_ff12s25.txt, fig_qcrit_ff12s25.png; σ = 1 version *_s1)
Definitions: Q = det(∇u) (= ½(|Ω|² − |E|²) if incompressible); Qn = (|Ω|² − |E|²)/(|Ω|² + |E|²), from −1 (strain) through 0 (simple shear) to +1 (rotation). d_Q = distance to the nearest Q = 0 line. 12 fine-flow runs.
- Smoothing matters. At σ = 1 cell, MPCD noise makes Q = 0 lines everywhere (random points: P(d_Q < 1) ≈ 0.7). Use σ = 2.5 for d_Q.
- Head et al. reproduced (σ = 2.5, bulk, depth > 6):
  - +1/2 are closer to Q = 0 lines than random points: median d_Q 0.85 vs 1.41; P(d_Q > 3) suppressed ×0.14–0.31.
  - −1/2 sit in rotation-dominated cores, away from the lines: Qn +0.47, median d_Q 3.2, enriched ×3–6 at d_Q > 4.
- Wall band: +1/2 within 1 cell of Q = 0 for 0.48–0.72 of samples, vs random 0.30–0.40. Wall defects also sit preferentially on viscometric lines.
- Escapers vs residents (σ = 2.5):
  - Residents: Qn ≈ 0 ± 0.06, P(d_Q < 1) ≈ 0.64.
  - Escapers: Qn −0.10 (−40 t.u.) → +0.10…+0.25 (−30 to −12, rotation-dominated, ~2–3σ above residents) → 0 (−6) → −0.17 (departure, strain-dominated).
  - P(d_Q < 1) for escapers is 0.57–0.66 up to −10, then rises to 0.89 at departure and 0.95–0.99 after. The detaching escaper locks onto a Q = 0 line, more tightly than typical +1/2 at that depth (0.66–0.68).
  - σ = 1 gives the strain excursion (Qn −0.35…−0.40 from −14 to 0) but not the rotation phase. So the Qn time course is scale-dependent; the landing on Q = 0 at departure is seen at both scales.
- READING (links to vortdef): with its wall-side lobe erased, the escaper's core sits inside its remaining (bulk-side) lobe. That is rotation-dominated at the 2.5-cell scale and is the net counter-rotating vorticity seen before. As the dipole re-symmetrises during untilting, the core returns to simple shear and the detached defect rides a Q = 0 line. That is the free +1/2 state of Head et al. ("riding its own flow"). Escape = the wall +1/2 regaining its free-defect viscometric constraint. Caveats: scale dependence; residents do not show the rotation phase, so it is specific to the run-up, but the error bars are about 0.06.

### 2026-10-11 (night) — Tips vs flanks (energetics/tipcount.py; shape sweep, second half of each run)
Region by the ellipse angle parameter t = atan2(y/B, x/A): tip |cos t| > 0.866, flank < 0.5, mid otherwise. Rates are per resident wall defect in the region.

| aspect | escapes tip / mid / flank | rate per 1000 t.u., tip / mid / flank |
|---|---|---|
| circle | 4 / 3 / 4 | 0.006 each |
| 1.2 | 61 / 33 / 5 | 0.113 / 0.054 / 0.007 |
| 1.4 | 173 / 107 / 111 | 0.355 / 0.189 / 0.164 |
| 1.6 | 145 / 102 / 83 | 0.317 / 0.174 / 0.122 |
| 2.0 | 135 / 103 / 91 | 0.309 / 0.164 / 0.123 |
| 2.5 | 115 / 122 / 104 | 0.260 / 0.190 / 0.128 |

- At onset (1.2) escapes are almost only at the tips. From 1.4 on every region leaks: tips 2–2.5x the flank rate per defect, but only about 40 % of all escapes. The flank rate alone is 20x the circle rate.
- The fine-flow set (12 runs) has only 65 depleting escapes (24/31/10), too few to split flow quantities.
- Queued:
  - runs/fineflow3: 24 runs (ε 2.0 and 1.6 × 12 seeds, 6000 t.u., flow + defects every 2 t.u.), chained after runs/bulkfine (fineflow3/chain.sh, 20 in parallel).
  - Shape-sweep energetics rerun with the region attached: energetics_general.py shapes_reg -> etsg_shapes_reg.npy; plot with plot_etsg_reg.py.

### 2026-10-11 (night) — QUEUED: causal test of the tilt protection via flow alignment (runs/lambda)
- Idea: tilt follows λ (passive swirl: about 60° at λ = 2, about 80° at 4, similar at 1.2). Inward jet = 0.40 cos τ g(d) predicts a 2–3x weaker inward push for residents at λ = 4, so a much lower depletion rate in elongated cavities. λ = 1.2 is the control (similar tilt).
- 20 runs: ε 2.0 × λ ∈ {4, 1.2} × 8 seeds, plus circle × λ × 2 seeds. 3000 t.u., defects every 2 t.u., flow every 20 t.u. Chained after fineflow3 (lambda/chain.sh). λ = 2 reference: fineflow/fineflow2 ε 2.0 and the shape sweep.
- Judge it by whether the escape rate tracks the measured cos τ (and the flow rms), not by λ alone: λ also changes the active turbulence.

### 2026-10-11 — Cloud backup set up (OneDrive UC, rclone)
- rclone v1.75.2 in ~/bin/rclone (official build, checksum verified). Remote "onedrive:" = institutional OneDrive for Business (5 TiB). Login done once through an SSH tunnel on port 53682.
- backup/backup_overcharged.sh: `rclone copy` (never deletes remotely) of ~/Physics/Overcharged to onedrive:Backups/Overcharged/office.
  - Skips deleteme*, __pycache__, its own logs, and run directories with a live solver (pid, no DONE).
  - Exits if another copy is still running. Logs in backup/logs/.
- Cron: nightly at 03:00. First full upload (about 19.4 GiB, 17.6k files) started 22:55.
- Pending: a server backup (63 GB) needs rclone on the server plus its own token (tunnel login to 146.155.46.150).

### 2026-10-12 ~00:15 — Escapes WITHOUT wall activity (runs/bulkfine: wallActOff 4, 10 runs, flow + defects every 2 t.u.) vs activity everywhere (fineflow 12 runs) (fig_compare_bulk.png; *_bulk.txt)
- The masking works. In the passive band a +1/2 carries only 20–28 % of the free jet (0.07–0.10 vs 0.35; leakage from the active fluid beyond 4 cells), rising to full strength at 7–9 cells. Its vortex lobes are nearly absent (wall-side −0.03, bulk-side 0.2–0.4 of free, vs about 0 / 1.7 in the active case).
- The cavity state is different:
  - Wall-band circulation 0.016–0.03 (active 0.13).
  - The wall band is strain-dominated: Qn −0.5…−0.65 (active ≈ 0).
  - Residents sit at cos τ ≈ 0.2 with no flow-locked tilt.
- Depletion still happens: rates 0.2–1.0 per wall defect per 1000 t.u., small numbers (5–13 events per set). This includes the CIRCLE (0.2–0.5; the active circle has 0), consistent with the forced_local and actlayout result that bulk-only activity empties the layer at any shape.
- Escape kinematics are nearly identical in both cases (depth and inward speed vs lag overlap). Part of this is built into the event definition (resident < 4 cells for 40 t.u., then > 6 cells), so it is not informative by itself.
- What differs:
  - (i) NO untilting precursor. Bulk-only escapers stay at cos τ ≈ 0.2 until −20 t.u. and turn only DURING the motion (0.34 at −10, 0.57 at 0), vs a steady rise from −40 in the active case.
  - (ii) NO own jet in the band (core − ring ≈ 0 until −10). It appears only as the defect reaches the active zone (0.08 at departure, 0.14 at +4).
  - (iii) Background flow slightly INWARD (ring +0.01…+0.04), vs slightly outward (−0.01…−0.02) in the active case.
  - (iv) No counter-flow at the site, and no landing on Q = 0 (P(d_Q < 1) peaks at 0.75 vs 0.89–0.99).
- Force fit (bulk-only escapers): α·jet ≈ 0 and β·ring ≈ +0.02…+0.04 against μF0 ≈ −0.025. It needs a constant c = +0.05 ± 0.014 (active: +0.02 ± 0.01). In the bulk-only route the inward motion before −10 t.u. is NOT accounted for by local flow. Candidates: noise plus selection (conditioning on escape), and larger-scale bulk flow or bulk defects (to test).
- READING: two routes. (A) Active wall: tilt-protected layer, released by untilting, then the own jet carries the defect out (the jet/lobe/Q signatures). (B) Passive wall, active bulk: no protection and no jet. The wall is quiet and strain-dominated with weak circulation; defects drift out with background flow and fluctuations against a weak elastic hold (μ_n F0 ≈ −0.025), and the tilt follows the motion. Route B is how the bulk-only runs empty the layer even in circles. It is the route that the active wall's coherent circulation suppresses in active circles.

### 2026-10-12 ~00:10 — Paper outline: λ scan plan + venue discussion added (reports/paper_outline/outline.md/.pdf; previous version outline_2026-10-10.md)
- A: tilt law in active circles across λ ∈ {0.8, 1.2, 1.5, 2, 3, 4, 6} (k 6, 3 seeds; k 4/8 subset), about 33 runs, roughly 3 h office. B: causal escape-rate test (running, runs/lambda). C: optional continuum. Low-λ end is a real test (passive swirl at λ 1.2 gave 57–65° vs arccos 34°).
- Venue: PRL realistic (circle mechanism + tilt law across λ); PRX/Nat. Commun. if the escape physics is unified and the λ causal test works; Nat. Phys. only with experiments. To settle with Tyler.
- (00:20) Plan A queued: runs/lambda_circle (33 runs: k 6 × λ {0.8,1.2,1.5,2,3,4,6} × 3 seeds; k 4/8 × λ {1.5,3,6} × 2), chained after runs/lambda (lambda_circle/chain.sh, 28 parallel); director also every 20 t.u. (S fields for λ_eff).

### 2026-10-12 ~11:30 — λ TEST of the tilt protection (runs/lambda; energetics/lambda_an.py -> lambda_B.txt). λ = 2 reference: fineflow/fineflow2.
Per-run means, second half; k 6 (2k−2 = 10), area 2500.

| set | wall +1/2 | bulk charge | resident ⟨cos τ⟩ | depleting rate per 1000 t.u. | v_rms | wall v_t rms |
|---|---|---|---|---|---|---|
| ε 2.0, λ 1.2 (8 runs) | 7.12 | +2.45 | 0.238 ± 0.016 | 1.01 ± 0.13 (43 ev) | 0.46 | 0.055 |
| ε 2.0, λ 2 (4 runs) | 7.14 | +2.44 | 0.191 ± 0.025 | 0.65 ± 0.23 (14 ev) | 0.57 | 0.061 |
| ε 2.0, λ 4 (8 runs) | 7.65 | +2.18 | 0.143 ± 0.013 | 0.066 ± 0.046 (3 ev) | 0.70 | 0.063 |
| circle, λ 1.2 (2 runs) | 7.92 | +2.05 | 0.39 | 0.34 (4 ev) | 0.49 | 0.127 |
| circle, λ 2 (3 runs) | 9.98 | +1.01 | 0.46 | 0 | 0.74 | 0.36 |
| circle, λ 4 (2 runs) | 9.98 | +1.01 | 0.26 | 0 | 1.04 | 0.34 |

- (1) KINETICS follow the tilt, as predicted. In the ε 2 ellipse the depleting-escape rate falls 15x from λ 1.2 to 4 (1.01 → 0.65 → 0.066) as the residents' inward projection ⟨cos τ⟩ falls (0.24 → 0.19 → 0.14). The confound goes the other way: v_rms RISES with λ (0.46 → 0.70) and wall circulation is unchanged (0.055–0.063). So more turbulence with fewer escapes; the drop is not a weaker-flow effect. The response is strongly nonlinear (1.7x in ⟨cos τ⟩ gives 15x in rate), as expected for a threshold or rare-event escape.
- (2) BUT the steady-state wall count barely moves: 7.12 / 7.14 / 7.65, with no drift in time. Constant in 6 time bins, so not a slow transient. Gross turnover of the smoothed wall count halves at λ 4 (loss/gain 17.8/17.4 vs 33.8/33.8 per run at λ 2). The tilt protection controls the RATES (both loss and refill), not where the steady state sits. The depleted layer of elongated cavities (bulk charge about +2.2 to +2.5 instead of +1) is a steady-state property set elsewhere. This revives the "steady-state shift" hypothesis (2026-10-11 ~11:55): bulk activity in an elongated cavity favours a larger bulk charge; escapes are the pathway, and tilt and jets set how fast the pathway runs.
- (3) NEW: at λ 1.2 even the CIRCLE loses the overcharged state (wall 7.9, bulk +2.05, wall circulation 0.13 vs 0.36 at λ 2). Weak flow alignment breaks the coherent single vortex. λ 1.2 is NOT a neutral control. Plan A (runs/lambda_circle, λ 0.8–6, running) will map where the circle state survives.
- Implications for the story:
  - Section 6 should separate (a) what sets the steady-state wall population (cavity-level: coherence or bulk-charge preference) from (b) the escape pathway and its kinetics (tilt protection, jet blocking, untilting runaway). The λ test supports (b) causally but says (a) is not set by the tilt.
  - Next: explain (a) directly. Options: Frank energy of the measured configurations with 2k−2 vs fewer wall defects (+ bulk); the bulk charge as a function of bulk activity and coherence; λ dependence of the steady state at ε 1.4 (near the transition), where a small shift could matter.

### 2026-10-12 ~12:15 — Loss vs refill across λ (user: restoring events may also slow down → same steady state) (energetics/turnover.py -> turnover_B.txt; refill.py)
- All-channel turnover, ε 2.0, rates per 1000 t.u.:

  | λ | exits per wall defect | entries per vacancy | annihilations at the wall per wall defect |
  |---|---|---|---|
  | 1.2 | 3.57 | 5.33 | 0.52 |
  | 2 | 2.81 | 3.61 | 0.56 |
  | 4 | 0.46 | 0.98 | 0.15 |

  Depleting escapes are only 14–28 % of the exits.
- CONFIRMED (user's point): stronger flow alignment slows BOTH directions. Exits drop 6x and entries 3.7x from λ 2 to 4, so the balance and the wall count barely move (7.14 → 7.65). The tilt and jet blocking act on crossings of the wall band in both directions, not on the steady state.
- The crude balance n* = N0 k_entry/(k_exit + k_annih + k_entry) gives 5.2–6.1 vs measured 7.1–7.7. Gains are undercounted: pair creation in the band and untracked flicker are not counted. Use it for ratios only.
- Launched: steadyE.py A (Frank energy of measured steady-state configurations vs one refill or one ejection; circle, ε 1.4, 1.6, 2.0 at λ 2; ε 2.0 at λ 1.2 and 4; 4 cores) -> steadyE_A.txt.

### 2026-10-12 ~11:20 — HYPOTHESIS (user): self-limiting bulk charge
More bulk charge pushes the wall +1/2 closer to the wall. There less of their jet is transmitted (g ≈ 0.3 below depth 2 vs 0.58 at 2–3) and the elastic hold is stronger, so exits slow as the bulk charges up. The same repulsion pushes bulk +1/2 toward the wall, so entries rise. The steady state is where the two meet. Fits the earlier "Poisson-like with a restoring tendency" within runs.
- Test (energetics/hazard.py, running on fineflow 1/2/3 and lambda ε 2.0; 6 cores): within-run fluctuation of the bulk charge vs wall depth, exit hazard per wall defect, and entry hazard per vacancy. Elastic side: dE_eject vs n_bulk from steadyE_A.
- Expectation if true: depth decreases with dQ; exit hazard decreases with dQ; entry hazard increases with dQ.
- Note: office load about 56 on 32 threads (lambda_circle 28 runs + analyses), so expect slower runs.

### 2026-10-12 ~12:00 — Theory report: draft of the escape-dynamics section (reports/theory_report/sec_escape_draft.tex; not yet \input, PDF not rebuilt)
- Content:
  - depleting-escape definition;
  - intact elastic barrier (fig_etsg_shapes);
  - untilting → local own jet (fig_selfflow_ff12);
  - jet blocking g(d) and the closed law (fig_jetmap_ff12);
  - lobe erasure (fig_vortdef_ff12);
  - Q criterion / Head et al. (fig_qcrit_ff12s25);
  - force balance + budget table;
  - escape-in-one-paragraph box;
  - tips vs flanks (fig_etsg_reg);
  - passive-wall route (fig_compare_bulk);
  - λ test: rates vs steady state.
- Placeholder: "What sets the steady state?" (steadyE, hazard, lambda_circle).
- Test-compiled standalone in the scratchpad: 6 pages, no LaTeX errors (only the cross-ref to sec:shapes, defined in the main file).
- At finalisation:
  - Refresh the force-fit and Q numbers with the 36-run set (ff36, running).
  - Enlarge fonts in the multi-panel figures (compare_bulk, vortdef legends).
  - Update the abstract, the six-line box, the "Which ingredient removes the hold?" paragraph (the activity-layout result supersedes it) and the open-questions table.

### 2026-10-12 ~12:00 — What sets the steady state? (energetics/steadyE.py -> steadyE_A.txt; hazard.py -> hazard_ff36.txt)
**(1) The depleted layer is NOT an elastic minimum.** Frank energy of the measured steady-state configurations (16 frames per set):
- Moving one bulk +1/2 into the largest wall gap LOWERS the energy by 80–100 kBT in elongated cavities. Best mover, median: −82 to −99; P(< 0) = 0.75–1.00. λ 4: −40.
- In the circle, an 11th wall defect is also favourable (−28; P 0.81). Sanity check: passive cavities keep all 2k at the wall.
- Ejecting a wall +1/2 costs +51…+65 (cheapest), +77…+84 (mean) in every set.
- So elasticity everywhere wants a fuller wall. The overcharged bulk (+1 in circles, +2.2…+2.5 in elongated cavities) is held up by the active dynamics against a large elastic preference. The steady state is a dynamical balance of exit and entry fluxes, not an energy minimum. (Caveat: the passive point-defect model with relaxed cores does not include the active tilt.)

**(2) Self-limiting bulk charge (user's hypothesis): supported on the exit side.** Within-run fluctuations, 49 elongated runs, 1006 exits, 539 entries:

| dQ_bulk | wall depth − run mean (cells) | exit hazard per wall defect per 1000 t.u. |
|---|---|---|
| −1.5 | +0.22 | 3.6 |
| −1 | +0.12 | 3.7 |
| −0.5 | +0.05 | 3.0 |
| 0 | −0.01 | 2.3 |
| +0.5 | −0.06 | 2.1 |
| +1 | −0.10 | 2.1 |
| +1.5 | −0.16 | 1.3 |

- More bulk charge pushes the wall defects towards the wall, and they exit less (×1.8 from dQ −1 to +1).
- Elastic side (steadyE, within sets): ejection cost rises with the number of bulk +1/2 in 4 of 5 elongated sets (corr +0.29 to +0.67; e.g. ε 1.6: 80 → 112 kBT from 2 to 6 bulk +1/2), and refill becomes more favourable (corr down to −0.74). Weak statistics: 16 frames per set.
- BUT entries per vacancy also FALL with bulk charge (6.4 → 3.3; the −1.5 bin is an artefact of near-zero vacancies). Refill is not enhanced; it looks supply- or capture-limited. The restoring feedback acts through the exits.
- READING: the bulk charge limits itself mainly by suppressing exits. More bulk +1/2 compress the layer towards the wall, where the jet transmission g(d) is smallest and the ejection barrier higher, so fewer wall defects leave. Combined with the λ test, the steady state is where the exit flux (falling with bulk charge) meets an entry flux that is itself slow (capture-limited, and slowed by the same tilt and jet blocking). It is far above the elastic optimum.
- (12:15) User: refill falling with bulk charge is NEEDED for the bulk charge to rise, so it fits. The steady state is where the TOTAL exit and entry fluxes cross.
  Launched:
  - perdef.py ff36: per-defect escape model. Push P = 0.40 cos τ g(d) per wall defect and frame; logistic P(exit | push); does the push distribution explain the measured drop of the exit hazard with bulk charge? Depth and tilt trends per unit bulk charge give a saturation prediction; also total exit vs entry flux by dQ.
  - bulkmap.py: where the bulk charge sits in the shape sweep, by elliptic radius and |x|/A.
- (12:15) User: if the excess charge sits near the middle, it presses less on the tip wall defects → tips leak more. Launched bulkpress.py: push of the bulk charge on each wall +1/2 along the outward normal (2D defect interaction, no images), by aspect × region (tip/mid/flank), plus depth; complements bulkmap.py.

### 2026-10-12 ~12:30 — Per-defect model vs bulk charge; where the bulk charge sits (perdef_ff36_v1.txt, bulkmap.txt, fig_bulkmap.png)
**Per-defect jet push does NOT explain the feedback.**
- 49 elongated runs, 194,513 wall-defect frames. Logistic P(exit in 20 t.u.) on the push P = 0.40 cos τ g(d): b = +6.25 (the push predicts exits).
- With more bulk charge (dQ −1 → +1) depth falls 0.12 cells per unit, but cos τ RISES +0.037 per unit, so the mean push is flat or slightly up.
- Predicted exit hazard is flat (1.00 at all dQ), while measured falls 1.92 → 0.90 (relative to dQ = 0). The exit suppression with bulk charge comes from something outside the jet push (tilt × transmission).
- Candidates: the elastic hold or barrier growing with bulk charge (steadyE: ejection cost rises with n_bulk), or crowding of the bulk in front of the wall defect. Rerunning perdef with the bulk pressure (push of the bulk charge along the normal) as a covariate.
- Total fluxes (per frame): exits fall 0.041 → 0.023 and entries RISE 0.11 → 0.39 with dQ. Both are restoring, so the steady state is stable.
  - CORRECTION to the earlier note: entries per vacancy fall, but the total entry flux rises with bulk charge (vacancies grow faster).
  - Absolute fluxes do not balance (entries > exits at every dQ) because the event definitions are asymmetric (flicker, pair creation, annihilation); use trends only.

**Where the bulk charge sits** (shape sweep, defects > 4 cells from the wall, per unit normalised area):
- Circle: the bulk +1 sits as a ring just inside the band (net peak at ρ ≈ 0.8). The CENTRE carries net NEGATIVE charge (−0.9 per area at ρ < 0.3), i.e. −1/2 there.
- Elongated (1.4–2.5): net charge positive everywhere and rising towards the wall (+0.3 centre → +1.2–1.3 at ρ 0.8). Along the long axis the density per normalised area is nearly flat, rising slightly towards the tips (ε 2: +0.63 at |x|/A < 0.2, +0.95–0.99 at 0.8–0.9).
- So the excess charge is NOT concentrated in the middle. Per normalised area it is uniform to slightly tip-ward and hugs the wall. Mean |x|/A of bulk +1/2 is 0.42–0.45 at every aspect.
- Whether the tip wall defects nonetheless feel less pressure (geometry: the tip faces a narrow region) is measured directly by bulkpress.py (running).
- bulkpress.txt (shape sweep; push of the bulk charge on wall +1/2 along the outward normal, 2D defect interaction, no images):
  - Tips, mid and flanks feel the SAME push at every aspect (within 1–2 %: ε 2.0 tip 0.0543, mid 0.0534, flank 0.0537 πK/cell). The user's tip-pressure hypothesis is not supported: the tip weakness is not a pressure deficit. It is the lower tilt there (resident cos τ 0.46 vs 0.25–0.29, etsg_reg).
  - The push grows with aspect along with the bulk charge (0.030 circle → 0.055 at ε 2.5), yet wall defects sit DEEPER in elongated cavities (1.31 → 1.51 cells). Across shapes, bulk pressure does not set the depth. Within aspect, tips respond most to bulk-charge fluctuations (slope +0.017 vs +0.010 flanks).

### 2026-10-12 ~13:00 — Defect seeding: new solver option on branch `seeddir` (off anchmod, uncommitted) + generator
- Solver (niamh-overcharged, branch seeddir):
  - `oDist = 10` (DIRFILE): each particle's orientation is set from its cell's director angle, read from json `dirInitFile` (lines "ix iy theta"; unlisted cells along y), with ±0.01 rad jitter. 2D only.
  - Files: definitions.h (DIRFILE), SRDclss.h (dirInitFile), read.c (key), init.c (loadDirInit + orient case + load before setcoord), init.h.
  - Binary bin/mpcd_seeddir.out; builds cleanly. The anchmod working-tree edits came along onto the branch.
- Generator reports/master/theory/seeding/seed_director.py:
  - Field θ = Σ q arg(z − z_i) + ψ, relaxed with the binding.py Frank functional (K1 42 / K3 23.1) and the eMode-2 anchoring on the wall; cores relaxed.
  - Sampled at cell centres via cos 2θ / sin 2θ.
  - Options: --wall N (evenly in arclength, depth 1.6), --nplus/--nminus (inner ellipse) or --bulk "x,y,q;..."; total charge must equal k.
  - Writes dirinit.dat, dirinit_seed.txt, dirinit.png and sets oDist/dirInitFile in input.json.
- Validation in progress (runs/seeddir_test: ε 2.0, k 6, full = wall 10 + 2 bulk +1/2; depl = wall 4 + 8 bulk +1/2; 20 t.u., defects every step); checker check_seed.py.
- TODO: bit-identity check of the binary with oDist ≠ 10 against mpcd_wallact2.out.
- (13:35) Seeding VALIDATED:
  - The seeded topology survives exactly. First output frame (t 10): full = 12 defects, wall 10, bulk +1.0; depl = wall 4, bulk +4.0; total +6 always.
  - Positions: 10/12 (full) and 6/12 (depl) seeds within 2.5 cells at t 10, then decorrelating as the active bulk +1/2 move (depl: wall 3, bulk +4.5 by t 20).
  - Note: the first defect frame is at t = 10 even with warmUp 0, so check the frame cadence at t < 10 if needed.
  - Bit identity: with oDist 2, mpcd_seeddir.out reproduces mpcd_wallact2.out byte for byte (defects, flow, director; 100 steps, seed 777).
- (13:35) perdef rerun with the bulk pressure (perdef_ff36.txt):
  - Pressure + depth explain only a small part of the exit-hazard variation with dQ (predicted ratio 1.18 vs measured 1.92 at dQ −1; good on the high side: 0.93 vs 0.90).
  - Likely missing covariate: crowding of the layer itself. Low dQ means more wall defects; exits from a crowded layer are revolving-door visitors.
  - Next: add n_wall and the neighbour spacing per defect.
- (13:45) perdef + CROWDING (refit from perdef_ff36.npy):
  - Exit flag rate per wall defect vs wall count relative to 2k−2: 0.068 (full layer, 0 vacancies), 0.050 (−1), 0.037 (−2), 0.034 (−3, −4), 0.028 (−5).
  - With dQ the wall count goes from +0.1 to −5.1 relative to 2k−2.
  - Model jet + pressure + depth + crowding predicts the dQ dependence 1.64/1.41/1.19/1/0.86/0.75/0.66 vs measured 2.05/1.92/1.48/1/0.89/0.90/0.84; pressure + depth alone give 1.28 … 0.91. Crowding alone (with the jet): 1.44 … 0.64.
  - READING: the self-limitation is mainly a WALL-CROWDING effect. A sparser layer holds its remaining defects more strongly, consistent with the elastic model: a vacancy makes its neighbours more bound by +8–10 kBT (2026-10-10 avalanche analysis). As the bulk charges up, the wall empties, each remaining wall defect sits in a deeper well, and exits slow. "Elasticity becomes too strong" (user) holds, through the wall layer's own spacing rather than bulk pressure.
  - The jet push and crowding are roughly independent factors (jet coefficient unchanged at about 6.2).

### 2026-10-12 ~13:15 — 36-run escape signatures (fineflow 1+2+3) and the region split (forcefit/vortdef/qcrit_ff36, regsplit_ff36.txt)
- Confirms the 12-run picture with tighter errors:
  - Force fit: α = 0.87 ± 0.02, β = 0.98 ± 0.05, μ_n = 0.0016 ± 0.0004 (bulk +1/2 α 0.64, β 0.58).
  - Lobe depth profile: wall-side −0.13/0.10/0.35/0.48/0.60/0.78 at depth 1–2 … 6–8.
  - Q: escapers rotation-dominated at −20/−10 t.u. (Qn +0.10/+0.13 ± 0.02 vs residents −0.06/0.00); land on Q = 0 at departure (P(d_Q < 1) 0.89, 0.91 at +4).
- By region (tip 190 / mid 140 / flank 81 events):
  - SAME mechanism everywhere: jet coefficient 0.82–0.90; wall-side lobe recovers 0.04–0.10 → 0.18–0.20 → 0.40–0.58 → 0.67–0.77 (−20/−10/0/+4 t.u.); landing on Q = 0 at 0.87–0.91.
  - Difference: at the TIPS the elastic pull is dynamically absent (μ_n = 0.0005 ± 0.0005 vs 0.0024–0.0027 mid/flank) and background flow is carried fully (β 1.19 vs 0.9). Together with the less tilted tip layer, this is the tip weakness.
  - Flanks need a larger constant (+0.036): residual inward drift not captured.
- (13:20) QUEUED: runs/hysteresis (seeded starts, bin/mpcd_seeddir.out).
  - ε 2.0, k 6, area 2500: start wall 10 + 2 bulk +1/2 (bulk +1), wall 7 + 5 (+2.5), wall 4 + 8 (+4); 4 seeds each.
  - Circle controls: wall 12 + 0 (passive-like 2k) and wall 6 + 6 (+3); 2 seeds each.
  - 4000 t.u., warmUp 0, defects every 2 t.u., flow every 20 t.u.; 16 parallel; chained after lambda_circle (hysteresis/chain.sh).
  - Question: does every start converge to the same wall count (unique dynamical steady state, about 7.1 at ε 2) and do circles return to 10 + bulk +1?
- (13:35) Report draft updated to the 36-run numbers (force fit, lobes, Q, budget table, region split) and the crowding result (steady-state paragraph); fig_compare_bulk regenerated 2x2 with larger fonts. Waiting for lambda_circle (tilt section) before the full rebuild.

### 2026-10-12 ~14:40 — TILT LAW ACROSS λ in active circles (runs/lambda_circle, 28/33 done; energetics/tiltlaw.py -> tiltlaw.txt, fig_tiltlaw.png)
Area 2500, natural activity, 3 seeds at k 6 (2 at k 4).

| k | λ | mean \|τ\| | arccos(1/λ) | wall +1/2 vs 2k−2 | bulk charge | tilt coherence |
|---|---|---|---|---|---|---|
| 6 | 2 | 63.0 ± 0.4 | 60.0 | 10.04 / 10 | +0.98 | 0.99 |
| 6 | 3 | 71.7 ± 0.2 | 70.5 | 10.00 / 10 | +1.00 | 1.00 |
| 6 | 4 | 74.8 ± 0.1 | 75.5 | 9.98 / 10 | +1.02 | 1.00 |
| 6 | 6 | 77.9 ± 0.1 | 80.4 | 9.97 / 10 | +1.05 | 1.00 |
| 4 | 3 | 71.6 | 70.5 | 6.0 / 6 | +1.00 | 1.00 |
| 4 | 6 | 78.1 | 80.4 | 6.0 / 6 | +1.04 | 1.00 |

- For λ ≥ 2 the overcharged state holds (wall exactly 2k−2, bulk +1, coherence ≥ 0.99, circulation never flips) and the tilt follows arccos(1/λ) within 3° (+3° at λ 2, −2.5° at λ 6) at both k. The tilt law is confirmed as a prediction across λ.
- THRESHOLD: for λ ≤ 1.5 the circle state breaks.
  - k 6: λ 1.5 wall 8.9, bulk +1.6, coherence 0.76; λ 1.2 wall 8.2, +1.9, 0.59; λ 0.8 (tumbling) wall 7.5, +2.3, 0.29, circulation flips.
  - k 4, λ 1.5: 4.7/6, +1.6. k 8, λ 1.5: 12.2/14, +1.9, coherence 0.23 (1 run; 4 more running).
  - λ_c lies between 1.5 and 2 at k 4–8.
  - Mean |τ| at λ ≤ 1.5 (60–70°) is not the law: the layers are incoherent.
- Link to the escape mechanism (prediction, to test): a flow-aligned wall defect has inward projection cos τ = 1/λ. At λ 2 half its jet points into the bulk; at 1.2, 83 %. The circle layer may lose its protection below λ ≈ 2 for the same reason elongated layers leak: too much of the jet points inward.
  - Alternative ordering: low λ breaks the coherent vortex first, then the layer leaks.
  - To test: escape kinetics and turnover in the low-λ circles (lambda_circle, flow every 20 t.u.).
- Ordering test in low-λ circles (corr of tilt coherence(t) with wall count(t + lag), smoothed 20 t.u.): λ 1.2 rises from +0.22 at lag 0 to +0.42 at +60 t.u. (coherent intervals are followed by refilling, i.e. coherence leads); λ 1.5 broad and symmetric (+0.2 to +0.3); λ 0.8 weakly negative. No clean ordering; correlations are slow (≫ 60 t.u.). λ 1.2 drifts during the run (coherence 0.23 → 0.57, wall 7.4 → 8.1 from t 1500 to 3000), so it may not have reached steady state. The 'jets point inward below λ ≈ 2' link stays a hypothesis.
- (14:55) Theory report: tilt-law-across-λ paragraph + fig_tiltlaw inserted in sec:tilt of theory_report.tex (k 8 numbers to update when the last 5 runs finish).
- (15:10) Theory report REBUILT (20 pages, clean): new Sec. 7 escape dynamics (\input sec_escape_draft.tex); tilt law across λ in Sec. 5; activity-layout paragraph replaces 'Which ingredient (open)'; abstract, 7-line box, summary table and one-paragraph version updated; previous version kept as theory_report_2026-10-09.tex. Pending: k 8 numbers (last 5 runs).
- (15:10) lambda_circle complete: k 8 λ 3/6: 70.0/77.5° (arccos 70.5/80.4), wall 14.00/13.98 of 14, bulk +1.00/+1.06; λ 1.5: 12.3/14, +1.85, coherence 0.27. Report Sec. 5 updated with k 8 and rebuilt (final for Tyler pending the hysteresis result). Hysteresis runs started 15:04.

### 2026-10-12 ~15:30 — Ovals vs ellipses (user: no qualitative difference so far?)
- Known: one transition at aspect 1.25–1.37 for both, collapse on aspect, escapes rise with curvature in both (ellipses leak somewhat more at high κR), pinch irrelevant. Not known: the escape mechanism in ovals; all of this week's mechanism work used ellipses and circles.
- Launched:
  - (1) energetics/ovaldef.py on the shape sweep (225 ellipses + 405 ovals, same code, region by |x|/x_max): depleting escapes by region, untilting precursor (escaper vs control cos τ), resident tilt by region, exits vs wall crowding -> ovaldef.txt.
  - (3) runs/ovalfine: 12 Cassini runs (R_eq 1.5; ba 1.30, aspect 1.94, k 7; ba 1.60, aspect 1.51, k 6; 6 seeds; 4000 t.u.; flow + defects every 2 t.u.). Binary runs/cassini/bin/mpcd_cassini_e63155b.out. Chained after runs/hysteresis (ovalfine/chain.sh, 12 parallel).
  - Analyses planned: selfflow, jetmap, vortdef, qcrit. forcefit's F0 uses the ellipse elastic model, so its elastic term is invalid for ovals.
- Pending (2): Frank energetics for Cassini walls (needs a Cassini version of the binding grid).

### 2026-10-12 ~15:25 — OVALS vs ELLIPSES at the escape level (energetics/ovaldef.py -> ovaldef.txt; 225 ellipse + 405 oval runs, same code, region by |x|/x_max)
- Untilting precursor: identical. Escaper cos τ rises from 0.19–0.52 at −40 t.u. to 0.71–0.78 at departure at every aspect in both shapes; controls stay flat. Example at aspect ≈ 2: ellipse 0.24 → 0.73, oval 0.26 → 0.72.
- Tip weakness: the same in both. Residents at the tips are less tilted (cos τ 0.48–0.64 vs 0.12–0.44 on the flanks), and flanks leak least per wall defect.
  - Difference: in ovals the MID section leaks as much as or more than the tips (aspect 1.3: mid 0.24 vs tip 0.19; 1.7: 0.21 vs 0.18). In ellipses tips lead at aspect ≥ 1.6. Possibly the oval's flatter shoulders.
- Rates: comparable magnitudes (0.09–0.31 per wall defect per 1000 t.u. at aspect ≥ 1.4 in both). Onset in both: aspect 1.2 leaks almost only at tips and mid (ellipse 44/46/9, oval 99/123/17). Near-circular oval (aspect 1.1): rates ≈ 0.003–0.02.
- Crowding self-limitation: the same in both, WITHIN runs. Exit rate relative to the run mean vs wall count relative to the run mean (elongated):
  - ellipse −2: 0.53, −1: 0.59, 0: 0.95, +1: 1.43, +3: 1.89;
  - oval −2: 0.79, −1: 0.58, 0: 0.84, +1: 1.55, +3: 1.46.
  - (The pooled table in ovaldef.txt shows the opposite trend because it mixes sizes and k: larger cavities have both more vacancies and more exits. Use within-run only.)
- READING: no qualitative difference between ovals and ellipses at the escape level either (precursor, tip tilt weakness, crowding feedback). Quantitative: ovals spread their leak more towards the mid section. Flow-level signatures await runs/ovalfine.

### 2026-10-12 ~16:00 — HYSTERESIS TEST: the steady state is unique (runs/hysteresis, seeded starts, branch seeddir)
Wall +1/2 (mean over seeds) in time bins 0–20 / 20–100 / 100–300 / 300–700 / 700–1500 / 1500–2500 / 2500–4000 t.u.:

| start | n | wall +1/2 by time bin | late per seed | bulk charge early → late |
|---|---|---|---|---|
| ε 2.0, wall 10 | 4 | 9.7 / 5.0 / 6.8 / 7.3 / 7.2 / 7.3 / 7.16 | 7.05–7.24 | +1.1 → +2.43 |
| ε 2.0, wall 7 | 4 | 7.0 / 6.0 / 7.1 / 7.2 / 7.2 / 7.0 / 7.17 | 6.97–7.39 | +2.5 → +2.43 |
| ε 2.0, wall 4 | 4 | 4.3 / 5.9 / 6.9 / 6.9 / 7.0 / 7.3 / 7.24 | 7.11–7.36 | +3.8 → +2.39 |
| circle, wall 12 | 2 | 11.8 / 5.2 / 7.8 / 8.6 / 10.0 / 10.0 / 10.00 | | +0.1 → +1.00 |
| circle, wall 6 | 2 | 6.1 / 5.0 / 7.1 / 8.1 / 9.9 / 10.0 / 9.92 | | +3.0 → +1.04 |

- Every start converges to the same state: ε 2: 7.2 ± 0.1 wall defects, bulk +2.4, matching the unseeded fineflow value of 7.1. Circles converge to 2k−2 = 10 with bulk +1.0. No hysteresis.
- Relaxation: about 300–700 t.u. in the ellipse, about 1000 t.u. in the circle.
- The seeded passive-like circle (wall 12) ends at 10 + bulk +1: the active ejection, reproduced from a seeded start.
- Transient: ordered starts (wall 10 or 12) first crash to about 5 wall defects at 20–100 t.u. (onset of active turbulence from the aligned seed) and then recover.
- IMPLICATION: the depleted layer of elongated cavities is a unique dynamical steady state, not a history-dependent or kinetically trapped one. Together with steadyE (far above the elastic optimum) and the crowding feedback, the steady state is a flux balance set by the cavity's active state.
- (16:05) Report: hysteresis result added to the steady-state paragraph; summary table updated (unique steady state; ovals = ellipses at the escape level); open list renumbered. Rebuilt, 20 pages, clean.
- (16:15) Report page 1 fixed: the author/date line overflowed by 193 pt (one unbreakable line) and pushed the title block off-centre; split into two lines; box title 'six' → 'seven lines'. No overfull boxes left.
- (17:10) Server: user keeps the current order (circle seeds incl. negative winding → mirror shape sweep → later shape production s04–s10). Estimate: circle seeds 4–7 more days (893 neg. runs not started; neg. runs ~27.6 h each so far), mirror sweep ~2 days after.

### 2026-10-12 ~17:25 — Flow-level escape signatures in OVALS (runs/ovalfine: 12 Cassini runs, aspect 1.94 / 1.51, R_eq 1.5; *_oval.txt, qcrit_ovals25.txt)
Oval vs ellipse (36 fine-flow runs):
- Untilting and own jet: cos τ 0.23 → 0.75; core − ring excess 0.04 → 0.25 at departure vs u_jet cos τ 0.24 (bulk jet 0.320); ring flat (−0.02…0). Same as ellipses.
- Jet blocking g(d), oval vs ellipse: 0.47/0.43/0.60/0.83/0.94 vs 0.30/0.33/0.58/0.86/0.97 at depth 0–1/1–2/2–3/3–4/4–5. Slightly more transmission within 2 cells of the wall in ovals; same beyond. Jet strength near the wall ×1.2–1.3 of bulk, as in ellipses.
- Lobe erasure by depth, oval vs ellipse: −0.10/0.11/0.36/0.49/0.63/0.83 vs −0.13/0.10/0.35/0.48/0.60/0.78. Escaper wall-side lobe recovery: 0.01/0.08/0.20/0.40/0.55/0.83 (−40/−20/−10/−4/0/+4) vs 0.03/0.08/0.19/0.34/0.46/0.69.
- Q: rotation phase before escape Qn +0.10/+0.13 at −20/−10 t.u. (ellipse +0.10/+0.13); landing on Q = 0 lines 0.92 at departure, 0.89 at +4 (ellipse 0.89/0.91). Bulk −1/2 rotation-dominated (+0.14, d_Q 2.1).
- READING: the escape mechanism is the same in ovals and ellipses at every level checked (defect, flow, vorticity, Q). Only small quantitative differences (slightly weaker jet blocking at the very wall; more leakage from the mid section).
- (17:30) Report: oval paragraph added to Sec. 7; open item (vi) removed; rebuilt clean.

### 2026-10-12 ~17:25 — NEGATIVE-LAYER test with seeding (runs/neglayer, branch seeddir binary)
- Question: the continuum mirror (k −6) reaches the "negative-layer prediction" (14 = 2|k|+2 wall −1/2 + bulk +1); MPCD mirror runs stay near 12–13 with less bulk charge. Is MPCD kinetically blocked (pair creation needed), or is that its steady state?
- 14 runs, k −6, area 2500, eMode 2, 4000 t.u., defects every 2 t.u.:
  - circle active: start 12 −1/2 (bulk 0) and 14 −1/2 + 2 +1/2 (bulk +1), 3 seeds each; plus passive 1 seed each;
  - ε 2.0 active: same two starts, 3 seeds each.
- Seeded passive Frank energy: 14 + 2 is 370–445 kBT vs 12 + 0 296–332 (untilted, cores relaxed). The prediction state is not elastically favoured without tilt.
- 14 parallel, started 17:23.

### 2026-10-12 ~17:30 — Cassini Frank-energy model (reports/master/theory/cassini_binding/cassini_binding.py)
- CassiniAB(a, b) gives the attributes binding.Ellipse.energy needs: grid, inside mask (|z−a||z+a| < b²), wall from the polar parametrisation, arclength, outward normals from ∇(d1 d2), wall band. Anchoring k 2π s/P from (x_max, 0), as eMode 3. Includes contains(z).
- Validation: a = 0 reproduces the circle exactly (467.7 vs 467.7 kBT, same perimeter).
- steadyE.py generalised: Cassini for "shape": 1, EllipseAB(A, B) for any ellipse or circle. Ejection targets must be inside and more than 2 cells from the wall.
- Running: steadyE.py oval on ovalfine (aspect 1.51, 1.94) + shape-sweep ovals ba 2.50 (1.18) and 1.80 (1.38), R_eq 1.5, k 6; 8 cores.
- (17:30) shear_trigger_ff36 (36 fine-flow runs): the pre-existing "stagnation patch" is NOT confirmed.
  - Escape vs control sites from −60 to −36 t.u.: identical (−0.001…+0.006 vs 0.000…+0.006, ± 0.006).
  - The counter-flow at the escape site starts with the untilting (−24: −0.010; −12: −0.021; 0: −0.034 ± 0.006; controls ≈ 0), so it is consistent with the escaper's own jet.
  - No local flow precursor found: the untilting trigger is open, possibly a fluctuation of the tilt in incoherent flow.
- (17:40) QUEUED for tonight (user: limit to 8 cores / 16 threads; runs at 12 parallel, analyses ≤ 4):
  - (1) runs/lowlam_fine: circles k 6, area 2500, λ 1.5 / 1.75 / 2.0 × 4 seeds, flow + defects every 2 t.u. Escape mechanics at low λ: do circles leak through the same untilting/jet route, and in which order do coherence and layer fail? λ_c between 1.5 and 2.
  - (2) runs/lamact: circles λ 1.5 / 1.75 / 2.0 × activity 0.05 / 0.10 × 2 seeds. Does λ_c move with activity (cavity size relative to the active length) or is it set by alignment alone (cos τ = 1/λ)?
  - Chain: lowlam_fine/chain_night.sh, after neglayer (14 runs).
- (17:50) QUEUED (user): precursor.py on the shape sweep (ellipses + ovals, aspect ≥ 1.35): escaper vs control correlations with wall neighbours (gaps, tilt, tilt-sense agreement, layer coherence) and bulk defects (nearest ±1/2 distance, bulk count within 8, nearest +1/2 bearing/heading/orientation), lags −60…0 t.u.; 3 processes after steadyE oval.
- (17:58) QUEUED (user): rotdiff.py, is the untilting just rotational diffusion? (1) Kramers–Moyal drift/diffusion of cos τ from resident tracks vs escapers' drift at the same cos τ (excess = driven); (2) excursion profiles to cos τ ≥ 0.7, escapers vs returning resident excursions; (3) P(escape | excursion). Shape sweep elongated, 3 processes, after precursor.py. Caveat: frame-to-frame angle jitter is large (±0.1 in cos τ), so B includes measurement noise; drift is robust with the statistics.
- (17:35) steadyE_oval.txt (Cassini model; ovals R_eq 1.5):
  - Moving one bulk +1/2 into the largest wall gap lowers the Frank energy by 81–89 kBT at aspect 1.4–1.9 (P(< 0) 0.90–0.95). In the near-circular oval (aspect 1.2) the gain is 41 kBT (circle: 28).
  - Ejecting costs +45…+55 (cheapest), +63…+78 (mean).
  - Same as ellipses: the depleted oval layer is far above its elastic optimum.
- (18:05) QUEUED (user): lobes.py, per-lobe analysis of ovals (b/a 1.1, 1.3 pinched; 1.45, 1.6 convex; ellipses aspect ≥ 2 as transfer controls): bulk charge per lobe, P(Q_L = Q_R) and correlation, lobe circulation and co-rotation, and TRANSFERS between lobes (bulk ±1/2 and wall +1/2 crossing x = 0, per bulk defect per 1000 t.u.). After lamact + rotdiff, 4 processes. Test run (b/a 1.1, k 5): 100 +1/2 and 6 −1/2 neck crossings in 5000 t.u.; Q_L +1.36, Q_R +1.55.

### 2026-10-12 ~17:50 — FLAT CHANNEL LIMIT of the winding anchoring (user: ovals add nothing new; very elongated → vortex lattice?)
- Flat limit of the arclength rule: on a long straight wall the anchoring angle advances at a constant winding density q = k/P. Around a cavity the two long walls are traversed in opposite directions, so in a periodic channel θ = θ0 + 2πqx (bottom) and θ0 − 2πqx (top), with charge 2qL per period L.
- Solver: branch `channelwind` (off seeddir, uncommitted), new eMode 4 in oriBC (lc.c): anchoring angle = wall-normal angle + customAnch·2π·x/L_x (+ O[2]). customAnch = windings per period (multiple of 1/2). Bottom +m, top −m → charge 2m per period. InputGuide and SRDclss.h updated. Binary bin/mpcd_channelwind.out.
- Matching to the ellipse: ε 2, area 2500 has P ≈ 194, so q ≈ 0.031 and m ≈ 6 per wall over L ≈ 192.
- Validation running (runs/channel_test): L 96, W 40, m 3 per wall (charge 6 per period); passive and active (act 0.07), 300 t.u. Analysis: reports/master/theory/channel/channel_an.py (wall +1/2 per wall, bulk charge, wall tilt, mid-line vorticity sign changes and dominant wavenumber, snapshot).

### 2026-10-12 ~18:05 — Solver branches committed (user request), private repo niamh-overcharged
- anchmod 34c588c: imposed swirl, wallActOff, anchoring modulation (on main fa68326).
- seeddir 5ca1c79: oDist 10 / dirInitFile (on anchmod).
- channelwind 0d05ddc: eMode 4 planar winding (on seeddir).
- Each built cleanly at its commit. The final working tree is identical to the pre-commit backup. main stays at fa68326 (production). Not pushed.
- (18:15) Pushed anchmod, seeddir, channelwind to origin (private repo baloewey/NIAMH-MPCD-overcharged); upstream tracking set. main unchanged.
- (18:25) Analysis repo: git init in ~/Physics/Overcharged (branch main), whitelist .gitignore (code, notes, LaTeX, input.json/jobs/tsv, reports/results tables, report figs; no simulation data, no solver, no binaries); initial commit 2ff0df8 (4495 files, 58 MB, largest 1.9 MB; secret scan clean); local identity Benjamin Loewe <baloewe@uc.cl>. Waiting for the user to create an empty private GitHub repo, then push.
- (18:30) Analysis repo pushed: origin = git@github.com:baloewey/Overcharged-analysis.git (private), main tracks origin/main (126c3d9). Policy (user): commit only when authorized, ask periodically.
- (18:40) CHANNEL FIX: the top-wall anchoring failed because the wall sat on the domain edge (y = W). The solver shifts walls by the random grid shift, so the cut cells were often the empty halo row and the strong ghost anchoring never acted (traced with an instrumented scratch build: angle and sign were correct). The early bottom/top asymmetry, the total charge 0 and the −1/2 at the top were artefacts of this.
  - With the walls one cell inside the domain (domain height W + 2, walls at y = 1 and W + 1), the passive channel L 96, W 40, m 3 gives total charge +6.0, 6 +1/2 at EACH wall, none in the bulk; director rotation +5.93π (bottom row), −5.93π (top row), 0 (mid-line). Exactly as predicted; no code change needed.
  - RULE for eMode-4 channels: keep the walls ≥ 1 cell inside the domain. To add to the InputGuide (channelwind branch).
  - The passive channel is the flat analogue of the passive cavity (all 2·2m +1/2 at the walls). Active test running (runs/channel_test/active_margin, 400 t.u., flow every 5 t.u.).
- (18:55) ACTIVE CHANNEL (L 96, W 40, m 3, charge 6; channel_test/active_margin, 400 t.u.):
  - Starts passive-like (12 at the walls), sheds +1/2 into the bulk within ~100 t.u., then wall +1/2 7–8 (bottom 3.6, top 3.0), bulk charge +2.7 (45 % of the charge; aspect-2 ellipse: 2.4 of 6, 40 %), a few bulk −1/2.
  - Flow: a row of alternating vortices, ~8 mid-line vorticity sign changes per 96 cells, dominant wavelength ~48 cells (≈ W).
  - The channel looks like the limit of elongated cavities.
  - QUEUED after lamact (runs/channel, 8 parallel): L 192 with the aspect-2 ellipse winding density (m 6 per wall):
    - W 40 × 2 seeds; W 40 m 0 (plain homeotropic) × 2 as the no-winding control;
    - W 20 and W 60 with m 6 × 2 each (vortex spacing vs width);
    - 4000 t.u., flow + defects every 2 t.u.

### 2026-10-12 ~19:10 — Is the untilting just rotational diffusion? NO (energetics/rotdiff.py -> rotdiff.txt; Langevin null model)
Shape sweep, elongated (405 runs). x = cos τ, steps every 2 t.u.: 7.4 M resident step pairs, 59,900 escaper pairs.
- Residents: drift A(x) restores towards x ≈ 0.44 (+0.024/t.u. at 0.1, −0.020 at 0.9; relaxation ~20 t.u.). Diffusion B 0.005–0.03 per t.u.; frame jitter is included.
- Escapers (−40…0 t.u.): drift is opposite to the restoring drift at the same tilt. At x 0.6–0.8, escapers +0.005 vs residents −0.014 (excess +0.019, 39σ); 0.8–1.0: −0.005 vs −0.020 (+48σ); 0.2–0.4: +0.019 vs +0.008 (+14σ).
  - Caveat: conditioning on a future escape biases the drift even in a purely diffusive process (Doob h-transform), so the excess drift alone does not prove driving.
- Decisive test, excursion shape. Profiles aligned on the upward crossing of x = 0.7, at −40/−20/−10/−6/−2/0 t.u.:
  - Langevin model with the measured resident A(x), B(x): 0.41/0.32/0.35/0.39/0.53/0.80. It reproduces the residents' returning excursions almost exactly (data 0.41/0.31/0.32/0.39/0.51/0.79).
  - Escapers: 0.24/0.36/0.48/0.56/0.62/0.77. They start MORE tilted than residents (0.24 vs 0.41) and then ramp up steadily over ~40 t.u., well above the diffusive profile from −20 on (0.48 vs 0.35 at −10).
  - Diffusive excursions are short spikes (jitter plus fast fluctuations) that rise only in the last ~4 t.u.
- Only 0.7 % of resident excursions to x ≥ 0.7 end in an escape (692 of 101,252). Those sit deeper at the crossing (2.64 vs 1.32 cells), consistent with the depth-dependent jet transmission g(d).
- READING: the untilting before an escape is a slow, persistent, driven rotation, not rotational diffusion of the tilt. Something rotates the escaper over ~40 t.u. against the flow-aligning restoring torque, starting from a more tilted than average state.

### 2026-10-12 ~19:10 — Precursor correlations (energetics/precursor.py -> precursor.txt; 405 elongated runs, 2955 escapes)
Escaper vs control (same frame, matched |x|/x_max):
- LOCAL CROWDING precursor (before the untilting): the escaper's nearer wall neighbour is closer, gap_min 0.776 vs 0.800 of the mean spacing at −40 t.u. (−3.9σ). It grows to 0.763 vs 0.827 at departure (−10σ). The larger gap is the same at −40 (1.19 vs 1.23 → +2σ later).
- EMPTY FRONT precursor: at −40 t.u. fewer bulk defects within 8 cells (0.095 vs 0.135, −4.3σ) and the nearest bulk +1/2 farther (18.1 vs 16.9 cells, +5σ). The difference reverses later (−10…0: nearer, −5σ) as the escaper moves out.
- Neighbours: their cos τ falls (more tilted) and their tilt sense agrees less with the escaper's from −20 on (−3…−7σ). This follows the escaper's untilting, not before it.
- Layer coherence slightly higher around escapers (+3σ, 0.40 vs 0.38); small.
- Nearest bulk +1/2 heading towards the escaper at −20 (+5σ, 0.15 vs 0.05); orientation difference shifts sign between −40 and −10.
- READING: two weak but significant precursors present before the untilting starts: a crowded side (a close wall neighbour) and an empty bulk in front. Both are about 3–15 % effects at 4–5σ. They fit the crowding feedback (neighbour repulsion along the wall) and the absence of bulk +1/2 that would push back (earlier approach_range result). They make a site susceptible; they do not determine it (most crowded sites do not escape).

### 2026-10-12 ~19:20 — NEGATIVE-LAYER test result (runs/neglayer, seeded starts, k −6, area 2500)
Wall −1/2 in time bins 0–20 / 20–100 / 100–300 / 300–700 / 700–1500 / 1500–2500 / 2500–4000 t.u.:

| set | n | wall −1/2 by time bin | late per seed | late: bulk +1/2 / bulk −1/2 |
|---|---|---|---|---|
| circle active, start 12 | 3 | 12.0 / 12.0 / 12.6 / 13.3 / 13.5 / 13.6 / 13.51 | 13.44–13.62 | 2.3 / 0.8 |
| circle active, start 14 (+ bulk +1) | 3 | 14.0 / 13.9 / 13.4 / 13.4 / 13.8 / 13.6 / 13.48 | 13.39–13.65 | 2.4 / 1.0 |
| ε 2.0 active, start 12 | 3 | … → 12.98 | 12.86–13.05 | 1.9 / 1.0 |
| ε 2.0 active, start 14 | 3 | … → 13.09 | 13.04–13.12 | 2.1 / 1.1 |
| circle passive, start 12 | 1 | stays 12.0 | | none |
| circle passive, start 14 + 2 bulk +1/2 | 1 | relaxes to 12.0 by ~700 t.u. | | none |

- MPCD is NOT kinetically blocked. Active starts from 12 and from 14 converge to the same state: circle 13.5 wall −1/2 (bulk net ≈ +0.7); ε 2: 13.0 (bulk net ≈ +0.5). Unique steady state, between 2|k| = 12 and the continuum's 2|k|+2 = 14 (+1).
- Passive: 14 + bulk +1 relaxes back to 12 (the seeded +1/2 annihilate with wall −1/2). The negative-layer prediction state is not even elastically stable without activity, consistent with its higher seeded Frank energy (370–445 vs 296–332 kBT).
- So the MPCD/continuum mirror discrepancy is a genuine difference of the active steady state, not a kinetic trap. MPCD reaches a partial gain (≈ +1.5 of the +2) with a smaller bulk charge (+0.5…+0.7), consistent with the earlier MPCD mirror sweep (bulk +0.6–0.7). Candidate cause: MPCD fluctuations / pair creation and annihilation in the bulk versus the noiseless continuum. Open.

### 2026-10-12 ~21:15 — Low-λ circles with fine output (runs/lowlam_fine; tiltlaw_lowlam.txt, selfflow_lowlam.txt)
- Threshold region, k 6, area 2500, 4 seeds each:

| λ | wall +1/2 (of 10) | bulk charge | tilt coherence | tilt (deg) vs arccos(1/λ) | depleting escapes per wall defect per 1000 t.u. |
|---|---|---|---|---|---|
| 2.0 | 9.94 | +1.03 | 0.99 | 62.1 vs 60.0 | 0.034 (1 event) |
| 1.75 | 9.60 | +1.20 | 0.93 | 59.1 vs 55.2 | 0.21 (6) |
| 1.5 | 8.65 | +1.68 | 0.73 | 61.8 vs 48.2 (incoherent) | 0.39 (10) |

  The breakdown is a smooth crossover between λ ≈ 1.6 and 2, with escape rates at λ ≤ 1.75 comparable to elongated cavities (0.1–0.3).
- Escape mechanics in the leaky circles (λ 1.5 + 1.75; 16 depleting escapes, indicative): same route as in elongated cavities.
  - Escaper cos τ 0.34 → 0.83 at departure.
  - Core − ring excess 0.02 → 0.25 at departure vs own-jet estimate 0.255.
  - Ring flat.
- READING: weak flow alignment lets circles leak through the same untilt-and-jet route. Consistent with the link cos τ = 1/λ (more of the jet points inward at low λ). Ordering (coherence vs layer first) still unresolved; small numbers.

### 2026-10-12 ~22:20 — LONG CHANNEL: vortex lattice and overcharging (user idea; possible follow-up paper, maybe with wavy walls)
- Motivation: homeotropic active channels show activity-induced "dancing disclinations" (Shendruk et al., arXiv:1703.01531, Soft Matter 2017; not shown in MPCD as far as we know) and vortex lattices at intermediate widths. With the winding anchoring the channel also carries a net charge. Questions:
  - (1) Is there a vortex lattice, and how does its period depend on W?
  - (2) How is the imposed charge shared between wall and bulk, and do bulk +1/2 lock to the lattice?
  - (3) Does the winding change the dance (vs m 0)?
- Queued (runs/channel_long, chained after lamact, 2 parallel next to the L 192 scan): L 768, W 40, m 24 per wall (ellipse winding density, charge 48 per period) and m 0 (homeotropic control). 4000 t.u., defects every 2 t.u., flow every 10 t.u. ~650k particles each, ~20 h per run → tomorrow evening.
- Plan: choose further widths for long channels from the L 192 scan (W 20/40/60, m 6 and m 0), possibly on the server once free. Wavy walls: eMode 4 takes the angle from the local wall normal, so it should combine with NIAMH wavy walls (winding along x, fine for small amplitude).
- (22:30) Long channels ON HOLD (user will restart the office PC tomorrow, NVIDIA driver issue): channel_long/chain.sh killed; inputs ready in runs/channel_long (jobs.txt). Start after the reboot (./run_queue.sh jobs.txt 2), possibly with widths chosen from the L 192 scan. Restart window: after ~06:00 when runs/channel (L 192) has finished. After reboot: SSH agent /run/user/1000/gcr/ssh needs a desktop/keyring login.

### 2026-10-12 ~23:10 — λ threshold vs ACTIVITY (runs/lamact + lowlam_fine; circles k 6, area 2500)
**[RETRACTED 2026-10-13 (user): activities outside the working window are confounded. ζ 0.05 = onset effects (unstable wall layer, few dynamical defects); ζ ≥ 0.09 = compressibility (vortex empties the centre, rms Mach ~1); decision of 2026-10-01. The 'threshold moves with activity' reading below does not stand; only the ζ 0.07 rows are usable. Clean replacement: vary R/ℓa via cavity size at ζ 0.07 (circles R_eq 1.0/1.5/2.0 × λ 1.5/1.75/2.0).]**

| act | λ | wall +1/2 (of 10) | bulk charge | tilt coherence | P(circulation keeps its sign) |
|---|---|---|---|---|---|
| 0.05 | 1.5 | 7.8 | +2.1 | 0.38 | 0.51 |
| 0.05 | 1.75 | 7.7 | +2.2 | 0.37 | 0.67 |
| 0.05 | 2.0 | 8.9 | +1.5 | 0.78 | 0.93 |
| 0.07 | 1.5 | 8.65 | +1.7 | 0.73 | 0.91 |
| 0.07 | 1.75 | 9.6 | +1.2 | 0.93 | 1.00 |
| 0.07 | 2.0 | 9.94 | +1.03 | 0.99 | 1.00 |
| 0.10 | 1.5 | 9.95 | +1.03 | 0.96 | 1.00 |
| 0.10 | 1.75 | 10.0 | +0.98 | 0.98 | 1.00 |
| 0.10 | 2.0 | 9.94 | +1.03 | 0.98 | 1.00 |

- The λ threshold MOVES with activity: stronger activity keeps the overcharged circle state down to λ 1.5; weaker activity loses it even at λ 2.
- So it is not set by alignment geometry alone (cos τ = 1/λ). It combines with the cavity size relative to the active length (R/ℓa grows with activity). Weaker flow alignment needs a larger R/ℓa (stronger activity or a bigger cavity) for the coherent single vortex.
- Act 0.05 at λ 2 (area 2500) sits near the known circle onset (R/ℓa ≈ 0.4–0.6).
- Next: map the threshold in the (λ, R/ℓa) plane using ℓa measured at each activity, and check whether the shape-transition aspect shifts the same way.
- (23:25) QUEUED FOR AFTER THE RESTART: runs/lamsize (18 circles, ζ 0.07, R_eq 1.0/1.5/2.0 × λ 1.5/1.75/2.0 × 2 seeds; replaces the retracted activity scan) + runs/channel_long (2). One start script: runs/after_restart.sh (16 threads: 2 + 14). The long-channel widths may be revised from the L 192 scan before starting.

### 2026-10-13 morning — overnight results
**Channel width scan** (runs/channel, L 192, ellipse winding density m 6 per wall = charge 12; walls at y = 1 and W + 1; 4000 t.u., second half):

| W | wall +1/2 (bottom / top) | bulk charge | bulk share of the charge | dominant vorticity wavenumber per L (wavelength) |
|---|---|---|---|---|
| 20 | 7.0–7.4 / 7.0–7.3 | +4.9…+5.0 | 41 % | 6–7 (≈ 29 cells ≈ 1.4 W) |
| 40 | 6.5–7.2 / 6.6–7.5 | +4.7…+5.6 | 39–46 % | 3–5 (≈ 38–64) |
| 60 | 6.8–7.0 / 6.9–7.3 | +5.0…+5.2 | 42 % | 2–3 (≈ 64–96 ≈ W–1.5 W) |

- The overcharging fraction is ≈ 40–45 % at every width, the same as the aspect-2 ellipse (40 %).
- The vortex scale follows the width (wavelength ≈ 1–1.5 W), suggesting a W-set vortex row/lattice; ~19–24 mid-line vorticity sign changes per 192 cells.
- No net channel flow with winding (|centre v_x| ≤ 0.13, except W 20: ±0.4).
- Homeotropic control (m 0, W 40): essentially no defects (≈ one pair) and a strong unidirectional flowing state, centre v_x ≈ 2.6. This is far above MPCD's sound speed (c_s = 1, Mach ≈ 2.6), so treat it as SUSPECT (compressibility). No dancing state at this width/activity. The winding (defects) suppresses the net flow.

**Per-lobe oval analysis** (lobes.txt; 270 runs; ovals b/a 1.1/1.3 pinched, 1.45/1.6 convex; ellipses aspect 2/2.5 as controls):
- Bulk charge is shared symmetrically (Q_L ≈ Q_R, each ≈ half the total: +0.75 … +1.6 by size) but NOT locked: P(Q_L = Q_R) only 0.19–0.31, correlation ≈ 0 to −0.3 (slight anticorrelation in small cavities = exchange). No "+1 per lobe".
- Transfers between lobes are frequent: bulk +1/2 cross the neck 1.7–6.2 times per bulk defect per 1000 t.u. in pinched ovals vs 2.2–7.8 across x = 0 in ellipses. The narrowest neck (b/a 1.1) lowers the rate by only ~30–50 %. −1/2 cross rarely (0.1–1.1); wall +1/2 move along the wall through the neck 10–24 times per 1000 t.u.
- Rotation: pinched lobes tend to CO-rotate (b/a 1.1: 0.53–0.72) while convex elongated ovals mostly COUNTER-rotate (b/a 1.45/1.6: 0.13–0.40; ellipses 0.40–0.66). Per-lobe circulation is larger in ovals (0.19–0.33) than in ellipses (0.13–0.20). The neck changes the vortex pairing (cf. Wioland et al. 2016) but not the charge split. Possibly the one oval-specific feature; small effect, check seeds.

### 2026-10-13 ~10:00 — Bulk charge PER UNIT WALL LENGTH (energetics/perlength.py -> perlength.txt, fig_perlength.png; shape sweep 675 runs + channels)
- vs aspect (pooled over q = 0.019–0.072): circle 0.0052 (= +1/P, falls with size: 0.0071/0.0049/0.0037 at R_eq 1/1.5/2); ellipse 1.2: 0.0060; 1.4: 0.0092; 1.6: 0.0098; 2.0: 0.0102; 2.5: 0.0103. Ovals: same curve (1.29: 0.0081 … 1.5–2.5: 0.0095–0.0097).
  - Above aspect ~1.5 the per-length value PLATEAUS and is SIZE-INDEPENDENT (0.0096–0.0106 at all R_eq). Released fraction Q_b/k ≈ 0.29–0.31.
- vs winding density q = k/P (elongated, aspect ≥ 1.5):

| q | n | Q_b per length | released fraction Q_b/k | fraction of 2q kept at the wall |
|---|---|---|---|---|
| 0.015–0.025 | 63 | 0.0121 | 0.61 | 0.40 |
| 0.025–0.035 | 63 | 0.0122 | 0.41 | 0.60 |
| 0.035–0.045 | 24 | 0.0104 | 0.24 | 0.76 |
| 0.045–0.060 | 66 | 0.0089 | 0.18 | 0.82 |
| 0.060–0.080 | 99 | 0.0074 | 0.11 | 0.89 |

- At the channels' q (0.031) the elongated cavities give 0.0122 per length, released fraction 0.41, wall kept 0.60. That is identical to the channel (0.012–0.0145; 0.39–0.46; 0.55–0.62).
- READING: overcharging in elongated cavities is a LOCAL WALL PROPERTY set by the winding density q, independent of size, shape family (ellipse/oval/channel) and aspect once ≳ 1.5.
  - The bulk charge per wall length is ≈ 0.012 at low q and decreases at high q.
  - The released fraction falls steeply with q (0.61 → 0.11): crowded walls hold their defects. This is the crowding feedback, and it fits the earlier observation that high k partly restores the layer.
  - The circle stays topological (+1 total).
- (10:15) Scaling with cavity AREA (elongated, aspect ≥ 1.5):
  - At fixed winding density q, Q_b ∝ area^0.47 (fit with q bins; perimeter ∝ area^0.50). So Q_b ∝ P, i.e. ∝ √area: per-length constant (q 0.03: 0.0119/0.0122/0.0126), per-area falling (0.00117/0.00080/0.00062 at areas 1550/3487/6199).
  - At fixed k, faster: e.g. k 6: +1.47/+2.59/+3.60, ~area^0.65. This is because q = k/P falls with size and the released fraction rises at lower q.
  - Circles: +0.99/+1.02/+1.03, independent of area (topological).
- (10:05) Wall law added to the theory report (Sec. 7: 'A wall law for the overcharging' + 'The channel limit', fig_perlength; summary table) and to the paper outline (row 6b + 'Possible follow-up paper: active channels'). Channel at q 0.031: 0.0132 ± 0.0004 per length vs cavity plateau 0.0122 ± 0.0001 (aspect ≥ 1.5).
- (10:15) ADDED to after_restart.sh: runs/channel_narrow (16 runs: L 192, W 10/12/16/24 × winding m 6 (q 0.031) and homeotropic m 0 × 2 seeds; 4000 t.u., flow + defects every 2 t.u.). Runs after lamsize in the same 14 slots. Tests where the wall law Q_b/(LW) = 2c(q)/W breaks and the homeotropic laminar/dancing regimes.
- (10:25) Channel defect dynamics (runs/channel, second half):
  - WALL LAYER PINNED at the homeotropic anchoring points (spacing L/2m = 16 cells). Site contrast of wall +1/2 positions 1.6–2.3 (uniform = 1). Median |dx| along the wall 1.5–2 cells over 20 t.u. and 2–3 cells over 200 t.u. for defects that stay in the band (survivors only, n 6–30 at 200 t.u.; one W 20 seed 11). So the wall layer behaves as in the cavities: pinned sites with turnover by escapes and refills.
  - BULK DYNAMICAL DEFECTS:

| channel | bulk +1/2 | bulk −1/2 | pair creations per 1000 t.u. |
|---|---|---|---|
| W 20, winding | 10.6–10.8 | 0.7–1.0 | 20–26 |
| W 40, winding | 11.5–14.2 | 2.2–3.1 | 55–68 |
| W 60, winding | 17.0–18.1 | 7.1–7.7 | 95–98 |
| W 40, homeotropic | 1.5 | 1.5 | 10–12 |

  - The widest channel has many dynamical defects on top of the ~5 net charge: about 7.4 ± 1/2 pairs per 192 × 60.
  - At W 40 the winding raises pair creation 5–6× over the homeotropic channel (laminar). The imposed charge destabilises the laminar state and seeds active turbulence.
- (10:45) Theory report: new Supplementary material (appendix) S1 'Layer versus bulk defects near a wall' with fig_layer_profile (reports/master/theory/channel/layer_profile.py): no density gap, mobility crossover at 3–5 cells, weak site locking in channels; 23 pages, clean.

### 2026-10-13 ~11:00 — Bulk overcharge vs dynamical defect nucleation (energetics/nucleation.py -> nucleation.txt; hysteresis runs)
- ACROSS SHAPES (shape sweep, same size): overcharged elongated bulks have FEWER dynamical pairs and nucleate less per unit bulk area.

| R_eq | group | bulk charge | pair creations per bulk area per 1000 t.u. (×1e−3) | bulk −1/2 density (×1e−4) |
|---|---|---|---|---|
| 1 | circle | +0.99 | 7.6 | 5.1 |
| 1 | aspect ≥ 1.4 | +1.49 | 4.0 | 1.0 |
| 1.5 | circle | +1.03 | 11.2 | 7.8 |
| 1.5 | aspect ≥ 1.4 | +2.21 | 7.0 | 3.1 |
| 2 | circle | +1.03 | 16.6 | 8.1 |
| 2 | aspect ≥ 1.4 | +2.95 | 9.3 | 4.8 |

  Nucleation is ~40–50 % lower and the dynamical −1/2 density 2–5× lower. Confounded with shape and flow state (charge density and aspect are collinear; the regression cannot separate them).
- WITHIN RUNS (elongated): creation rate rises with the bulk-charge fluctuation (×0.82 at dQ −1.5, ×1.29 at +1, ×1.5 at +1.5…+2). Likely a common cause: active bursts drive both escapes and nucleation. Q_b also counts −1/2, so the coupling is partly mechanical. Not causal evidence.
- SEEDED STARTS (hysteresis): no usable window. The bulk charge converges within 20–100 t.u. (all starts at +3.0…+3.5, then +2.4), so the seeded difference is gone before nucleation can be compared.
- Channels: winding vs homeotropic (W 40) 5–6x more pairs, but that compares an active state with a laminar one; not a test.
- OPEN: needs a comparison at fixed geometry and flow with different bulk charge. Options:
  - (a) channels at fixed W with different winding density q (bulk charge per length 0.0074–0.012, released fraction 0.11–0.61);
  - (b) positive vs mirror (negative) winding at the same shape (bulk +2.4 vs +0.5…+0.7 at aspect 2; mirror sweep on the server after the circle seeds).
- (11:05) ADDED to after_restart.sh: runs/channel_qscan (L 192, W 40, m 3/9/12 per wall × 2 seeds; with runs/channel m 6 and m 0 → q scan 0, 0.016, 0.031, 0.047, 0.063). Test of bulk charge vs pair nucleation at fixed geometry. Runs after channel_narrow in the 14 slots.
- (11:15) ADDED to after_restart.sh (user): runs/channel_mixed, per-wall winding at L 192, W 40 (wall charge = customAnch on the bottom, −customAnch on the top): negative (−6/−6, total −12), mixed zero-net (+6/−6: +1/2 layer bottom, −1/2 top), unequal same sign (+6/+3, total +9), unequal opposite (+6/−3, total +3); 2 active seeds each + passive checks for negative and mixed (10 runs). After channel_qscan.

### 2026-10-13 ~11:30 — "Topological capacitor": mixed-winding channels and the topological current
- Idea (user):
  - Mixed channel (wall charges +Q bottom / −Q top) = capacitor-like. Passive candidates: an unscreened defect-free texture rotating uniformly along x (cost ∝ K q² W L) vs screened plates (+1/2 layer bottom, −1/2 layer top, field-free bulk; cost ∝ L). That implies a passive crossover width.
  - Active: released +1/2 (bottom) and −1/2 (top) can cross and annihilate with the opposite layer; the walls renucleate pairs. This would be a steady topological current (+y), the analogue of breakdown/discharge; with bulk pair creation possibly an avalanche (Townsend).
- QUEUED: runs/channel_mixedW (W 10/16/24/60, wall charges +6/−6, 2 active + 1 passive each; 12 runs) after channel_mixed in after_restart.sh.
- ANALYSIS: reports/master/theory/channel/channel_current.py (--group for summaries). Per run: layer occupancies (±1/2 at each wall), net charge profile across the width, current J_cross (mid-line crossings Σq(up − down)) and J_flux (Σq v_y over bulk defects), per unit length per 1000 t.u.; nucleation and annihilation rates and where (bottom band / bulk / top band).
- Tests:
  - all-positive W 40 active: J_cross −0.007 (≈ 0, as required by symmetry); nucleation 0.36, annihilation 0.52 per length per 1000 t.u., 67–75 % in the bulk.
  - homeotropic: J 0, nucleation and annihilation only in the bulk.
  - passive all-positive: all +1/2 at the walls, J 0.
  - The early buggy run (top wall effectively of opposite sense) is a passive mixed channel at W 40: SCREENED state (6 +1/2 bottom, 6 −1/2 top, no bulk charge), J 0. So at W 40 the passive capacitor is screened; the narrow widths will show whether the unscreened texture takes over.
- (11:40) Scope split (user): paper 1 = cavities (channel only as a wall-law check); channels = possible paper 2; topological capacitor/current = possible paper 3 (or part of 2). Recorded in the paper outline.
- (11:50) Seeds: keep 10 seeds for the circles (user). Shape production: decide 5–7 vs 10 seeds per point after the circle seeds finish, from their measured seed-to-seed variance (more seeds near the transition aspect 1.2–1.4 and at the extreme sizes).
- (16:08) Office PC rebooted (NVIDIA driver 595.99.02, GPU working; SSH agent OK). Started runs/after_restart.sh: channel_long (2) + lamsize → channel_narrow → channel_qscan → channel_mixed → channel_mixedW (14 slots).
- (18:50) reveal.js presentation of the theory report (for a possible meeting with Tyler): reports/presentation/index.html (reveal.js 4.6 vendored in reports/presentation/reveal/, figures + escape movie in figs/; 20 slides with speaker notes; works offline). Zip: reports/presentation.zip.
- (19:00) Presentation fix: figure heights were in vh (viewport units), which reveal.js scales a second time, so slides overflowed on screens larger than 1280×800. Switched to slide pixels (430 px single, 220 px stacked); all 20 slides checked at 1920×1080 and the densest at 1024×768. Zip rebuilt.
- (19:10) Presentation published as a private Artifact: https://claude.ai/artifact/4xDw6JqRYcgv6oboMWh1Ri (from reports/presentation/publish.html = index.html without the html/head/body skeleton; reveal.js files, 16 figures and the escape movie as supporting files). Share with Tyler via the page's Share menu.

### 2026-10-13 — Theory report: channel figure, abstract sentence, self-consistency pass (24 pp, clean build)
- Added fig:channel (channel_fig.py) and the wall-law/channel sentence to the abstract.
- Consistency fixes: mirror (−k) paragraph updated with the seeded neglayer result (not kinetically blocked; 13.5/13.0 wall −½, open item (iv) reworded, added to Supported); active tilt ≈62° (was 60°); story line (5) and Sec 1 "grows with size" → wall law; along-wall mobility 0.02 vs normal μ_n distinguished; "wall defects not advected" restricted to the resident layer; Sec 6 restoring-force collapse labelled as effective (elastic force unchanged, Sec 7); lowlam_fine λ 1.75 crossover added to the tilt-law section; λ 1.2 circle numbers reconciled (7.9–8.2); escape-section data sentence now says ff12 for jet/g(d), ff36 for vorticity/Q/force fit; draft/TODO comments removed; one-paragraph version gains the wall law; channelwind branch listed.
- Flagged, not changed: Sec 6 medians "1, 3, 4.5 vacancies / bulk 1.1, 1.8, 2.5" do not satisfy bulk = 1 + vacancies/2 (separate medians over k, density?); needs a check against the shape-sweep table.
- (later, user: "λ part hard to follow, reads as if something is wrong but doesn't matter") Rewrote the Sec 7 λ-test paragraph as prediction → confirmed (rates) → refill slowed by the same gate → wall count barely moves → "gate" interpretation; compared λ 2 vs 4 (λ 1.2 set aside as not clean: circle coherence fails there). Rewrote the Sec 5 low-λ "possible link" bullet: tilt alone changes rates (λ test), so the circle breakdown is more likely coherence loss with escapes as the route; ordering still unresolved. Abstract wording adjusted to match.
- (later, user: the arccos-correction part was confusing) Restructured Sec 5 after the flat-wall theory: "The mechanism in MPCD" (director rotation, tilt = 2× rotation) → tilt law across λ → "The last few degrees": residual +3.0/+1.2/−0.7/−2.5° at λ 2/3/4/6 (falls with λ, changes sign), three candidates with verdicts (boundary layer absent, +16° predicted; order-dependent alignment ±2–3°, not the trend; splay≠bend +2–3° at λ 2 and 4, explains λ 2 only) → bottom line → other caveats. Old splay≠bend table and "Not settled" list folded in. Abstract, Supported row and open item (iii) updated (boundary layer no longer credited). Note: energetics/tiltlaw.txt now holds the lamsize/lowlam-style rows; the λ 2–6 table lives in this file (2026-10-12 ~14:40).
- (later) fig_perlength redrawn (user: grey bands unclear): colour = q bin, channel as a diamond at an extra 'channel' x position in (a,b) and at q 0.031 in (c), elongated vs circle-like in (c); REPLOT=1 mode in perlength.py; caption rewritten.

### 2026-10-13 ~19:25 — Overnight (user going home)
- Chain (detached, survives the session): lamsize 12/18 at 19:20 → channel_narrow (16) → channel_qscan (6) → channel_mixed (10) → channel_mixedW (12), 14 slots; channel_long (2, L768) until ~19:00 Oct 14.
- runs/overnight_analysis.sh (detached, nice 19, log runs/overnight_analysis.log): per finished set → energetics/tiltlaw_lamsize.txt; channel/channel_narrow.txt, current_narrow.txt, channel_qscan.txt, current_qscan.txt, current_mixed.txt, current_mixedW.txt. tiltlaw.py now takes OUT env (so tiltlaw.txt is not overwritten again).
- On resume: read these outputs; lamsize → λ threshold vs R/ℓa; qscan → nucleation vs bulk charge; mixed → capacitor current.
- (later) fig_qcrit_ff36s25 panel (a) redrawn (user: mean Qn misleading): per-event fractions — within 1 cell of a Q=0 line and Qn>0, escaper vs resident, s.e. bands, random baseline; REPLOT=1 mode in qcrit.py. Text: mean Qn ~0 for residents is cancellation (median |Qn| 0.31); escapers at departure 89% within 1 cell but |Qn|≈0.43, strain side (64% Qn<0); run-up 'rotation-dominated phase' reworded as a shift (61–65% Qn>0 vs 42–51%).
  1-cell tolerance justified in text: tracker positions on the cell grid (20848 of ~22k samples at integer positions in fineflow/eps1.0_s1).

### 2026-10-13 ~22:45 — λ threshold vs cavity size (runs/lamsize, k 6, ζ 0.07, 2 seeds per cell; tiltlaw.py per R_eq)
| R_eq/ℓa | λ | wall (of 10) | bulk | tilt mean (arccos) | coherence | v_rms |
|---|---|---|---|---|---|---|
| 1.0 | 1.5 / 1.75 / 2 | 9.84 / 9.54 / 10.08 | +1.08 / +1.23 / +0.96 | 53.4 / 56.0 / 58.2 (48.2/55.2/60.0) | 0.66 / 0.70 / 0.97 | 0.44 / 0.48 / 0.58 |
| 1.5 | 1.5 / 1.75 / 2 | 8.93 / 9.97 / 9.96 | +1.54 / +1.02 / +1.02 | 61.7 / 60.7 / 64.5 | 0.90 / 0.99 / 1.00 | 0.63 / 0.75 / 0.87 |
| 2.0 | 1.5 / 1.75 / 2 | 9.02 / 9.64 / 9.72 | +1.49 / +1.18 / +1.15 | 61.1 / 62.2 / 64.6 | 0.98 / 0.99 / 1.00 | 0.82 / 0.88 / 0.97 |
- PRELIMINARY (2 seeds). Small cavity: loses COHERENCE at low λ but keeps the full wall. Large cavities: stay COHERENT (0.98 at λ 1.5, R_eq 2) but lose about one wall defect; R_eq 2 leaks slightly even at λ 2 (9.72, bulk +1.15).
- So coherence loss and layer loss are separable and go opposite ways with size → contradicts the Sec 5 wording "more likely the loss of coherent flow" (written today). Layer loss at fixed coherence fits the jet route (cos τ = 1/λ inward share) with stronger flow in larger cavities.
- Tilt at λ 2 grows with size (58.2 / 64.5 / 64.6) along with v_rms → the residual over arccos(1/λ) depends on flow strength (cf. boundary-layer term ∝ √s). Revisit "last few degrees".
- (later) qcrit panel (a) re-redrawn (user: colours/labels confusing, resident label wrong): two stacked sub-panels (top: within 1 cell of Q=0; bottom: Qn>0), colour = class (escaper black, control grey dashed, random dotted). 'Resident' control = random other wall +1/2 in the same frame and region, >3 cells away, NOT required to stay in the layer (same definition in energetics_general, vortdef, qcrit); now defined once at the start of Sec 7 and named 'other wall +1/2' in this figure.
- (later) 'resident' → 'control' throughout Sec 7, figure legends (plot_etsg, plot_etsg_reg, vortdef [REPLOT mode added], qcrit) and slides. Kept the genuine meaning only in the escape definition (reworded: 'a wall +1/2 that has been in the band ≥ 40 t.u.'); 'per resident wall defect' → 'per wall defect'; λ-test cos τ is the mean over ALL wall +1/2 (lambda_an), now said so; Sec 3 'resident layer' → 'wall layer as a whole'.
- (23:10) QUEUED lamsize batch 2 (user: 2 seeds too few; add circles smaller than ℓa): runs/lamsize/jobs2.txt, 67 runs — s3–s5 for R_eq 1/1.5/2 × λ 1.5/1.75/2 (27) and s1–s5 for R_eq 0.5/0.75 × λ 1.5/1.75/2/3 (40; λ 3 as a high-alignment reference for the small circles). Same settings (k 6, ζ 0.07, density 20, 30000 steps, domain ceil(2R)+2). Detached runs/lamsize/run_batch2.sh waits for the channel chain (PID 7052) then runs 14 in parallel (~5–6 h, set by the R_eq 2 runs), then writes energetics/tiltlaw_lamsize_byR.txt (all seeds, per R_eq).
- (later) fig_vortdef_ff36: legends of (a, b) moved below the panels, colour = location, style = escaper/control (user: legends collided with curves).
- (later) fig_vortdef_ff36 rebuilt with a sampling-geometry sketch (a), explicit y labels (b: circulation-signed ω; c,d: lobe vorticity / Ω0 = 0.12, defect frame; d core ω on its own right axis), caption spells out every definition.
- (later) Fig. 5 (dirstokes) redrawn: (b) bulk and wall-band predictions on the run's single fitted c (the old per-source scale flipped the wall band's sign and made it look as good as the bulk); wall band now small and anti-correlated as stated in the text. (a) axis = mean azimuthal velocity; caption spells out the model (f = -div Q from measured director, wall band <4 cells vs bulk, 2D no-slip Stokes, one c per run, 40 x 50 t.u. intervals).
- (later) Fig. 7 (shear_tilt): U defined (v_theta = U f(r), f = measured passive-swirl profile normalised to max 1 → U = peak azimuthal speed, model units); (c) relabelled: groups by λ, colour legend for the three flow profiles, title 'flow variants at U = 1 (above threshold)'; caption rewritten.
- (later) Fig. 7 (c) grouped by λ (shaded blocks, values on bars, arccos(1/λ) per group, legend below); Fig. 8 (b) relabelled: x = predicted excess tilt from Eq. (tilt) with each run's d, s, λ and γ1/K = 1 in model units; y = measured excess tilt; λ colour bar; caption spells out the 35 runs and that nothing is fitted.
- (23:50) TILT COHERENCE: two definitions were in use. Most scripts (shape sweep/phase diagram, ejection figure, Q-tensor models, kickout, halfk) use C = |<exp(i tau_j)>| per frame, time-averaged; tiltlaw.py (Fig. 9, lowlam_fine, lamsize) used |<sign tau>|. tiltlaw.py now reports C (and the sign version as an extra column); fig_tiltlaw and tiltlaw.txt regenerated (λ 2–6 table restored). Changed numbers: circles λ≥2 0.96–0.99 (was ≥0.99); k6 λ 1.5/1.2/0.8: 0.75/0.61/0.42 (was 0.76/0.59/0.29); k8 λ1.5 0.50 (0.27); lowlam_fine λ 2/1.75/1.5: 0.96/0.90/0.73 (0.99/0.93/0.73).
  lamsize with C: R_eq 1: 0.73/0.75/0.94 (λ 1.5/1.75/2); R_eq 1.5: 0.86/0.96/0.98; R_eq 2: 0.94/0.97/0.98 → the size trend (small loses coherence, large stays coherent but leaks) holds.
- Report: C defined once in Sec 2 (new paragraph 'Tilt coherence', with the tilt τ definition); Fig. 9 panels labelled (a–c), title 'ζ = 0.07, no imposed flow' replaces 'natural activity'; caption rewritten.
- (later) Fig. 3 (tilted_structures) titles name the elasticity (one-constant K1=K3=32.5; anisotropic K1=42, K3=23.1, K1/K3=1.82); Fig. 6 (effective_bc): 'patchiness' replaced by its definition (rms over the circle of θ−φ−c at distance δ inside the defect ring, δ in site spacings), axis labels and caption rewritten.

### 2026-10-14 ~00:20 — House figure style (user: capitalised labels, larger fonts)
- reports/figstyle/sitecustomize.py: monkeypatch loaded via PYTHONPATH (reports/figstyle/run.sh script.py). Capitalises the first ordinary word of axis labels, titles, suptitles, legend entries/titles and annotations (after "(a) "); leaves symbols/variables (k, s, q, λ, cos, R_eq, k02, $..$) alone; "rms" → "RMS". Scales every font size by 1.3 (FIGSTYLE_SCALE).
- reports/figstyle/regen_report_figs.sh: redraws all 23 theory-report figures and copies them to theory_report/figs and presentation/figs. fig_cont_ellipse had no script (inline on 10-09) → recovered from the transcript as theory_report/fig_cont_ellipse.py.
- Per-figure layout fixes for the larger text: shorter axis labels (definitions moved to captions), legends below panels (etsg, etsg_reg, qcrit, cont_ellipse, phase_shapes), panel titles added where missing (etsg_shapes, phase_shapes, maps_circles), bigger canvases (dirstokes, shear_tilt, channel, make_figs, vortdef). Report rebuilt (25 pp, clean); slides republished.
- (00:45) Vector figures: sitecustomize also writes a PDF next to every PNG (TrueType fonts, pdf.fonttype 42; collections > 1500 elements and meshes rasterised at ≥300 dpi in the PDF only). Report uses the PDFs (\DeclareGraphicsExtensions{.pdf,.png}; extensions dropped from \includegraphics): all 23 figures vector, 18–68 KB each (channel 418 KB), report PDF 1.4 MB, fonts TrueType/Type 1 only (no Type 3). Slides keep the PNGs.
- (01:00) Page-21 viewer error 'insufficient data for an image': matplotlib 3.6 writes low-colour rasterised images as palette PNGs with bit depth < 8 but omits BitsPerComponent from the predictor parameters (poppler fails; fixed upstream in 3.7). Patched PdfFile._writeImg in sitecustomize; fig_channel regenerated; all figure PDFs and the report scanned with mutool: no errors.
- (01:30) Final pass for undefined terms/quantities (report now 26 pp, clean). Added: MPCD expansion and units (cells, t.u., kBT, density 20); how ℓ_a is measured; λ as the Bretherton parameter (tumbling < 1); aspect ratio (overall length/width; ε for ellipses) and R_eq = √(area/π); "wall charge k−1"; ψ_rel, s, P, φ_n in the lock; lock "concentration" → mean resultant length; registry; μ units; swirl amplitude A; R² for the Stokes model; continuum circulation ⟨v_φ⟩/v_rms; K and γ1 in the Ericksen–Leslie equation; flat-wall wavenumber renamed q → κ (q is the winding density); shearSusc described in words; γ1 units (τ clashed with the tilt); S2, S4 defined and the MPCD mean-field strength renamed U → U_MF (U is the imposed flow in Fig. 7); depleting escape pointer; vacancies, two-vortex criterion, escape, tip/mid/flank regions, A (long semi-axis), Fokker–Planck force/potential units (kBT_eff); Q, Ω, E, Q_n and "viscometric" in the text; "best mover" rephrased; eMode 4 → "a new anchoring option"; Q_b and released fraction; missing comma in "Where things are"; date/"New" line updated.
- Figure: selfflow legend u_self → u_jet (matches text); selfflow.py got a REPLOT mode (from selfflow_ff12.txt), used by regen_report_figs.sh (saves ~11 min).
- (10:08 Oct 14) lamsize batch 2 started early (user OK to exceed 16 threads): 14 slots in parallel with the channel chain; summary → energetics/tiltlaw_lamsize_byR.txt when done.
- (10:40 Oct 14) PAPER SKELETONS (user: one PRL, one long format): reports/paper/prl/main.tex (revtex PRL, 3 pp now, ~1200 words of draft text: abstract, observation, anchoring geometry, ejection, tilt law, beyond circles + wall law; intro/discussion TODO; 4-figure plan as composites in source comments) and reports/paper/long/main.tex (revtex PRE, 9 figure slots, 10 sections + 3 appendices, lead paragraphs drafted, %% CONTENT lists per section with report sources). Shared reports/paper/refs.bib (entries from memory, marked TODO-VERIFY; Klein 2025 and Loewe 2024 need full refs). Both build clean.
- (11:00) refs.bib verified via Crossref/arXiv: Head 2024 Nat. Phys. 20, 492–500; Shendruk 2017 Soft Matter 13, 3853; Kozhukhov & Shendruk 2022 Sci. Adv. 8, eabo5788; Kozhukhov, Loewe & Shendruk 2024 Commun. Phys. 7, 251 (key Kozhukhov2024, replaces Loewe2024); Klein et al. arXiv:2503.10880 (preprint only); Doostmohammadi 2018 Nat. Commun. 9, 3246.
- (11:20) Paper drafts: authors Benjamin Loewe (Facultad de Física, PUC Chile, Santiago 7820436) and Tyler N. Shendruk (School of Physics and Astronomy, Univ. of Edinburgh, Peter Guthrie Tait Road, Edinburgh EH9 3FD, UK), as in arXiv:2503.19664.
- (11:40) PRL Supplemental Material skeleton: reports/paper/prl/si.tex (8 pp, clean): S1 methods, S2 observables, S3 circles, S4 lock/pinning, S5 energetics, S6 imposed vortex, S7 continuum models (passive Q, active Q, flat wall), S8 tilt mechanism and residuals, S9 hydrodynamics, S10 shapes, S11 escape, S12 steady state + wall law + channel, S13 layer vs bulk, S14 movies; 18 figure slots (report PDFs), Tables S1–S2 placeholders, S-numbering.
- (12:00) Added Norton, Baskaran, Opathalage, Langeslay, Fraden, Baskaran, Hagan, PRE 97, 012702 (2018) (verified via Crossref) to refs.bib; cited in the PRL intro and the long-format intro plan. Relevance: imposed net charge in disks; active flow screens the anchoring; an APS 2017 abstract of theirs reports net +1 kept at strong activity, +3/2 below a threshold → read the paper for overlap with our 2k−2 + bulk +1 result.
- (12:10) Correction: Norton 2018 had been read in full on 2026-10-05 (RESEARCH 'Related literature and novelty'); PRL intro rewritten from those notes (boundary charges −1/0/+1, topological boundary layer, insensitivity; Klein: up to 10 excess +½ below the turbulence threshold, braids) with the novelty 'topology is screened, geometry is not'. Klein's PNAS listing unverified (only arXiv found).
- (12:40) Literature update (user: add refs + cited-by search):
  - Klein et al. IS published: PNAS 123(28), e2516670123 (8 Jul 2026), doi 10.1073/pnas.2516670123 (key Klein2026). The 10-05 'PNAS' note was right; today's first searches missed it.
  - Added (Crossref-verified): Hardoüin et al. Nat. Commun. 13, 6675 (2022); Hardoüin et al. Commun. Phys. 2, 121 (2019); Opathalage et al. PNAS 116, 4788 (2019); Theillard, Alonso-Matilla & Saintillan, Soft Matter 13, 363 (2017). Cited in the PRL intro and the long-format intro plan.
  - Cited-by Norton 2018 (Semantic Scholar, 57 papers): none imposes |k| > 1 in the turbulent regime or uses non-circular cavities with winding anchoring. Worth reading for the paper: Houston 2025 Proc. R. Soc. A "The role of elastic anisotropy in active nematics" (our K1 ≠ K3 pinning); Schimming et al. 2023 PRE "Friction-mediated phase transition in confined active nematics" (disks); de Oliveira et al. 2023 J. Mol. Liq. (ring-shaped cavity); Mirantsev 2021 EPJE (circular region); Tang et al. 2024 (circular rotations); Vafa et al. 2023 PRE (periodic orbits on cones); Gulati et al. 2022 and Wagner et al. 2021/2023 (channels, for paper 2).
  - Cited-by Klein 2026: none indexed yet (published July 2026). No study of elliptical/oval cavities with winding anchoring found → shape novelty stands (as of 2026-10-14). Repeat the search before submission.
- (13:00) Mirror placeholders added: PRL 'Mirror control' paragraph + Fig. 4(d) plan; SI new Sec. S12 'Mirror control: negative winding' (S13–S15 renumbered: steady state, layer, movies); long format new Sec. 'Mirror control' with Fig. 10 slot. All pending the server mirror shape sweep.

### 2026-10-14 ~14:00 — Bulk overcharge vs dynamical defects: DOPING (mass-action) TEST (reports/master/theory/channel/doping.py → doping.txt, doping_all.txt, fig_doping_all.png)
- Hypothesis (user's question; my proposal): pairs created at a fixed rate per area, annihilating ∝ n+ n− → n+ n− = n0² and n− = (√(δ²+4n0²) − δ)/2 (semiconductor-like doping).
- Shape sweep (675 runs, per run): log n− anti-correlates with log ρ_Q (−0.45/−0.64/−0.58 at R_eq 1/1.5/2), slope log n− vs log n+ −0.82/−0.49/−0.17; n+n− not constant. Confounded (k, aspect, size).
- Channels: the q-scan (W 40, m 3–12) barely changes the BULK excess (1.1–1.8e-3 per cell², crowding caps the release) → poor lever; the width is the lever (excess 1.0 → 6.5e-3 from W 60 to 16, but geometry changes too). 4/6 q-scan runs partial.
- Result: −½ ARE suppressed by the excess, but NOT by mass action. Per run n− falls ~3× (0.74 → 0.24e-3) while the excess rises 6×; n+n− rises (0.7 → 2.0e-6) instead of staying constant. Locally (2-cell strips across the width, pooled): n− drops from 0.80 to ~0.40e-3 as the excess goes 0.4 → 0.85e-3 and then PLATEAUS at ~0.4 up to an excess of 4.2e-3 (slope log n− vs log n+ = 0.01).
- Pair creation per bulk area shows no systematic drop with the excess in channels (5–14e-3 per 1000 t.u.); annihilation exceeds bulk creation (−½ created in the wall bands, 15–75 % of creation in narrow channels).
- Robust new observation: strong LIKE-CHARGE EXCLUSION among bulk +½: within 6 cells, ++ pairs are 0.09–0.17 of a uniform distribution (−− 0.1–0.9; +− 1.0–2.2) in every winding channel → the excess +½ form a correlated, spaced-out liquid (Coulomb-gas-like?). Homeotropic m 0 (laminar flow, few defects): ++ 0.4–1.6, +− 6–7.
- Next: finish q-scan; mixed-winding channels (winding walls, zero net charge) as the control at equal walls; g(r) of +½ (screening length, any order?); same in cavities (circle sizes at fixed +1 vary the density cleanly).
- (14:40) PAIR CORRELATIONS (channel/pair_corr.py → pair_corr.txt, fig_pair_corr.png; g(r) = same-frame pairs / pairs between frames 200 t.u. apart, which removes inhomogeneous one-body densities; two normalisation bugs fixed on the way).
  RETRACTION of the 14:00 'like-charge exclusion' as an overcharge effect: the neutral periodic box (ζ 0.07, L 200, no walls) shows the same exclusion. g++ for r < 4 / 4–8 / 8–16 cells: box 0.10/0.44/0.83; circles (bulk +1) 0.09/0.28/0.68; ellipses 0.05/0.35/0.83; channels W 40–60 0.05/0.39/0.87; narrow 0.05/0.41/0.88. g−− and g+− also collapse on the box curve (g+− ≈ 2.3–2.7 at r < 4: pairs). The 14:00 numbers compared with a uniform gas, not with neutral turbulence.
  Small overcharge-specific differences at most: like neighbours slightly more spaced relative to a Poisson gas (nn ratio 1.23–1.28 vs 1.14 in the box); circles slightly stronger exclusion at 4–16 cells. Homeotropic channel: too few defects, noisy.
  → Short-range structure of the bulk defects is that of ordinary active turbulence; the overcharge changes the NUMBERS (−½ suppressed, plateau) not the local correlations. Not a paper result as it stands.
- (15:00) SAME SHAPE AND SIZE, DIFFERENT WINDING k (from nucleation.npy; user recalled this test — the 10-13 analysis had compared ACROSS shapes only):
  within each (family, aspect, R_eq), k varies the bulk charge density 1.4–3.3×. Elongated (27 groups): creation rate and −½ density RISE with the bulk charge density (median log-log slopes +0.83 and +0.63; positive in 27/27). Circle-like (8 groups): ~flat (−0.06, −0.11).
  → At fixed geometry, more bulk charge comes with MORE nucleation, not less. The across-shape reduction (40–50 % fewer creations in elongated cavities) is a shape/flow-state effect, not a charge effect. Caveat: k also changes the wall (winding density, crowding), so this is not a pure charge test either.
- Lifetimes: no tracked lifetime statistics existed. Mean bulk −½ lifetime from Little's law (density / creation rate per area): circles 67/71/47 t.u. at R_eq 1/1.5/2; elongated 22/45/52 t.u.; within a shape, weakly shorter with more bulk charge (slope −0.12, 18/27 negative). A track-based lifetime distribution (pair creation → annihilation, partner identity) still to do if needed.

### 2026-10-14 ~16:00 — Pair lifetimes (tracked) and 'duck-duck-goose' in circles (reports/master/theory/bulkdyn/: pair_lifetime.py → pair_lifetime.txt, fig_pair_lifetime.png; ddg.py → ddg.txt)
- PAIR LIFETIMES (kick_identity.track, gaps ≤ 3 frames bridged; pair births ≥ 4 cells from the wall). 73–81 % of detected pairs are FLICKER (−½ seen < 10 t.u.), in the neutral box too. Established pairs (≥ 10 t.u.): mean −½ lifetime 32–38 t.u. in every geometry (box 35, circles 37, ellipses 38, channels 32). The −½ almost never dies with its own birth partner (≤ 1 %); among identified deaths it annihilates with another +½ (box 25 %, circles 24 %, ellipses 47 %, channels 40 % of all; the rest LOST = track ends without a nearby +½ death, 44–71 %, tracker caveat). Partner exchange is generic active turbulence; the higher 'other +½' share in overcharged ellipses/channels may reflect the excess +½, but the lost fraction differs between groups, so not conclusive.
- DUCK-DUCK-GOOSE (circles, shape sweep, tracked +½; entry bulk > 6 → band < 4, exit band → > 6): CONFIRMED in circles. After an entry the next exit comes 4–5× sooner than chance (median 20/32/77 t.u. vs null 84/188/338 at R_eq 1/1.5/2; P(exit within 20 t.u.) 0.50/0.44/0.34 vs 0.15/0.11/0.06), and it is ANOTHER defect in 56/68/79 % of cases (the rest: the newcomer bounces back out). The goose is usually not the neighbour: median 2.2/3.2/3.7 site spacings from the entry site, only weakly enhanced within 1.5 sites (0.34/0.25/0.23 vs uniform 0.29/0.18/0.13) → a GLOBAL occupancy constraint (the ring holds 2k−2), not a local push. Elongated (ellipse 2.0): no coupling (14 vs 16 t.u., P 0.67 vs 0.62): there the layer is not count-locked.
  Orbiting: bulk +½ in the 4–8 cell ring travel ~2 site spacings per 40 t.u. along the wall (circles); those about to enter travel only ~0.5 in the 40 t.u. before entry → they circulate with the vortex and enter after slowing near a spot.
- (12:15) ELASTIC SHIELDING FIGURE (user: Japan-poster potential; master §29 had the theory: same-sign images, exact shell theorem in circles, conformal map for ellipses).
  Version 1 = poster construction: potential of one bulk +½ from the wall layer alone (layer charges + images + boundary winding), one-constant, Mathematica (bulkdyn/shield_wl/probe_maps.wl, conformal map from theory/shielding_ellipse.wl), layers from layers.json (measured frames of the fine-output runs + idealised full layers). Range of the LAYER's field (probe's own wall term removed): inside the 0.6 / 0.8 level, kBT:
     circle ideal 1.6 / 25; circle measured (9–11 wall) 16 [7–138] / 49; ellipse 2 ideal (10) 14 / 59; ellipse 2 measured (5–9 wall) 176 [54–369] / 270.
  → shielding is near-exact in the circle, ~10× weaker in an ideally filled ellipse (geometry; matches master Tab. ellth), and lost in the real depleted ellipse layer: the gaps leave wells of ~100–300 kBT reaching into the bulk (consistent with refill lowering the Frank energy by 80–100 kBT). The degradation is mainly the vacancies, not the shape.
  Version 2 = PMF −ln ρ+ from the simulations: circle shows the sharp 10-site ring with a flat, depleted interior; ellipse layer blurred, interior flatter (s.d. 0.43 vs 0.73).
  First attempt (shielding.py, numerical grid solver over the full measured frames incl. bulk defects) superseded: dominated by bulk-bulk interactions, not the layer. Figure: bulkdyn/fig_shielding_layer.png (plot_shielding.py).
- (12:30) Shielding added to the drafts: PRL paragraph + figure (to be merged into Fig. 4), SI new section 'Elastic shielding by the wall layer', long format content + Fig. 7b slot. fig_shielding_layer copied to theory_report/figs.
- (12:50) Shielding figure, user: depleted ellipse map 'saturated almost everywhere, low only near one sector'. Cause: display (cap at 100 kBT vs a 180 kBT range), not a bug; and my description ('wells at the vacancies') was wrong. Check: 7 EVENLY spaced wall +½ in the same ellipse → layer-field range 27 kBT (0.6 level), vs 176 for the measured 7. So the uneven spacing (gaps), not the count, kills the shielding: the gaps leave unscreened ½ charges whose 2D log field (~πK ln(ratio) ≈ 100 kBT per gap) TILTS the whole interior, downhill towards the largest gaps. Figure redrawn with a square-root colour scale to 300 kBT and contours every 20; paper text corrected (PRL, SI, long). Caveat noted: the other bulk charges, absent from the 'layer alone' construction, partly screen this field.
- (13:10) Per-wall-defect hold in depleted ellipse layers (user: do some defects get more shielding, others less?). bulkdyn/shield_wl/binding_per_defect.wl, binding_full.wl: energy cost to move each wall +½ 6 cells inward (one-constant, conformal map), 16 measured aspect-2 frames (5–9 wall defects).
  Layer alone: ideal full layer 38 kBT each (37–42); 7 evenly spaced 18; measured median 21, 10–90 % [0, 36], some even slightly unstable. With the measured bulk defects included: median 65 [47, 88]; none below 10 kBT.
  Dependence on the adjacent gap (with bulk): 55 (normal spacing) → 66 (1.2–1.8 sites) → 75 (1.8–2.5) → 68 (> 2.5) kBT: defects next to a gap are held MORE, not less (same as the anisotropic model: vacancy → neighbours +8–10 kBT). Within a frame the hold varies by ~45 kBT across the layer.
  → Depletion makes the layer uneven, not weak: the excess bulk +½ push the remaining wall defects onto the wall (layer alone 21 → with bulk 65), and gap neighbours are the best held. Consistent with crowding self-limitation (sparser layer holds its defects more strongly). The bulk feels an unscreened, tilted field; the wall defects are not elastically exposed.
- (13:40) Shielding, positive vs MIRROR (user): bulkdyn/plot_shielding_mirror.py → fig_shielding_mirror.png, shielding_mirror.txt; each panel minus its interior (0.8-level) minimum, own diverging scale. Mirror layers from runs/neglayer (k −6, area 2500, active; 24 frames each).
  Range of the layer's field felt by a bulk +½, inside the 0.6 level (median, kBT): k +6 circle ideal 1.6 / measured 16; ellipse ideal 14 / measured 176. k −6 circle ideal 0.3 / measured 116 [13–254]; ellipse ideal 4.9 / measured 167.
  → Ideal −½ rings shield even better than +½ rings (14 sites), but the MEASURED −½ layers shield poorly even in the circle: 12–15 wall −½ (not count-locked at 14; mean 13.5) and irregular. The +½ circle layer is the regular, count-locked one (pinned commensurate ring of 2k−2) and the only well-shielded measured case.
  User's point: the +k plots matter for duck-duck-goose: a vacancy leaves an unscreened charge whose 2D log field spans the whole cavity and draws a bulk +½ into the hole (and an extra wall defect pushes one out) → explains why the 'goose' is non-local (median 2–4 sites from the entry).
- (14:10) VACANCY vs FULL CIRCLE (bulkdyn/vacancy_field.py → fig_vacancy_field.png, vacancy_field.txt): measured frames of the k 6, area-2500 circles (fineflow ×3) with 9/10/11 wall +½; layer-alone field felt by a bulk +½ (analytic disk, one constant).
  Range in the inner 60 %: full (10, 703 frames) 17 kBT [IQR 13–21]; one vacancy (9, 28 frames) 143 [135–149]; one extra (11, 22 frames) 140. With a vacancy the interior force points at the gap in EVERY frame (<cos> = +1.00).
  Duck-duck-goose (shape-sweep circles R_eq 1–1.5, 30 runs): next exit 22 t.u. after an entry vs 104 at random times; another defect 60 %; goose distance nearly uniform around the ring (median 2.5 sites).
  → Mechanism for the non-local exchange: a vacancy is an unscreened half-charge whose log field (~140 kBT) spans the cavity and pulls a bulk +½ into the hole; an extra wall defect does the same in reverse, pushing a defect out anywhere on the ring.

### ▶ TO REVIEW WITH BENJAMIN (requested 2026-10-14 ~14:30) — shielding / vacancy maps (bulkdyn/fig_shielding_layer, fig_shielding_mirror, fig_vacancy_field)
User's observations:
 1. No potential drop near some of the layer defects (expected wells/peaks at every site).
 2. An asymmetry emerging in the circle maps.
 3. The ellipse tips look particularly well shielded.
Candidate explanations to check (not verified):
 - (1) Masking: maps are drawn only inside the 0.9 level (≈ 2.8 cells from the circle wall, more at the ellipse flanks); wall defects sit 1–4 cells deep, so the deeper ones lie inside the plotted region (local well/peak visible) and the shallow ones outside (invisible). Check: plot to 0.98 level and mark each defect's depth.
 - (1) Colour scale: TwoSlopeNorm anchored at the interior minimum clips near-wall structure differently per panel.
 - (2) Real layers are irregular (spacing, depth): any non-uniformity is an unscreened multipole whose log field tilts the interior (shell theorem needs equal charges at equal radius and spacing). Check: field of the measured layer vs the same count re-spaced evenly at the measured depths; and vs depth irregularity alone.
 - (3) Tips: uniform flux per arclength + higher curvature; the layer-alone construction omits bulk charges; check the per-defect field vs curvature and against the ideal-layer ellipse maps.
 - Also re-derive the probe-own-image subtraction for the ellipse (H term with Q including the probe).
Paper text using these maps carries a red TODO pointing here (PRL shielding/vacancy paragraphs).
- (15:20) REVIEW, circle checks done (bulkdyn/review_shielding.py → fig_review_shielding.png, review_shielding.txt):
  (1) 'No drop near some layer defects' = MASKING + self term. Plotted to 0.99 R, every wall +½ shows its local structure: a repulsive peak for the +½ probe (layer term), with attractive pockets in the gaps next to the wall. Shallow defects (depth 0.6–1.1 cells) lay outside the earlier 0.9-level mask, so their peaks were invisible; in the 'total' maps the probe's own wall repulsion (−q² ln(R²−r²), ~60 kBT at the edge) fills the gap pockets. For the +½ layer the per-site signal is a peak, not a drop; wells appear only for the −½ (mirror) layer.
  (2) Circle asymmetry = angular jitter of the layer, at the tracker's resolution. Full layers (703 frames): layer-term range inside 0.6 R 17 kBT [14–21]; same defects re-spaced evenly at their measured depths 2.4; measured angles at the mean depth 16.8; ideal 1.4. Angular s.d. only 1.7° (≈ 0.8 cell); an ideal ring with ±0.5 / ±1 cell random position noise already gives 7 / 13 kBT. So per-frame interior tilts of ~10–20 kBT in full circles are dominated by small displacements (partly grid quantisation of defect positions) and average out over time; the vacancy/extra fields (~140 kBT) are an order of magnitude above this noise floor.
  (3) Ellipse tips and the ellipse self-term subtraction: not yet checked (Mathematica, conformal map).
- (17:20) CURRENT DECOMPOSITION in mixed channels (channel/current_decomp.py → current_decomp.txt; user: no crossings → is the flux carried by bulk-created defects?).
  UNIT CORRECTION: J_flux (channel_current.py) and the new motion/event terms are Σ q v_y / L, i.e. the current INTEGRATED over the bulk width (32 cells at W 40). Per-line current = these / 32. Then motion+events agree with the mid-line crossings: mixed +6/−6: total/32 = +0.026 / +0.035 vs crossings +0.021 / +0.027; unequal opposite (+6/−3): +0.036 / +0.031 vs +0.023 / +0.039. So no contradiction: the flux is small and IS the crossing current.
  Significance: symmetric controls give per-line currents of −0.003 to +0.020 (crossings) / −0.003 to +0.013 (total/32) → mixed channels at +0.02–0.04 are at most 1.5–2× the noise; not established. Long/wide mixed runs needed (channel_mixedW tonight).
  Carriers of the motion term in the mixed channels: +½ that came from the + wall band moving up (+0.48–0.74, integrated units) dominate; pair-born +½ ~0; −½ moving down +0.24 / −0.06 (+0.26 / +0.43 in the unequal-opposite channels). Pair EVENTS oppose the current: creations −0.18 to −0.22 (pairs are born with the +½ BELOW the −½), annihilations +0.09 to +0.19. → No sign that bulk-created pairs carry the current; if anything the charge is carried by defects ejected from the walls, and pair creation polarises against it.
- (17:45) Theory report updated (28 pp, clean): λ×size window table (Sec 5), shielding + vacancy + duck-duck-goose + mirror shielding (2 figures), bulk-overcharge consequences (doping fails, nucleation at fixed shape, correlations/lifetimes generic), narrow + mixed channels with the current units correction; summary table and open item (ii) updated.
- (17:50) Session crash ~17:20 killed the restyled movie render only (runs unaffected: mixedW 6/12, long channels running). Rerendered detached: reports/master/theory/channel/movie_mixed0_Qb+6_Qt-6_act_s1.mp4 (t 3000–3600, 300 frames, 15 fps; vorticity smoothed 2 cells, flow arrows, +½ comets, −½ tripods, wall counts, bulk charge, cumulative mid-line crossings). Script channel_movie.py.
- (18:00) Control movie: reports/master/theory/channel/movie_W40_m6_L192_s1.mp4 (+6/+6, same window t 3000–3600). At t 3300 the control already shows +3.0 cumulative mid-line crossings (mixed: +2.0): over a few hundred t.u. the running crossing count fluctuates by several units in a symmetric channel, consistent with no established current.

### 2026-10-15 ~12:20 — Overnight results: mixed winding at several widths, long channels; server
- Office: all queues finished (mixedW 21:59, long channels done); machine idle. Backup 03:16 OK. Server: circle seeds finished 07:05 UTC (2170/2170); mirror shape sweep started automatically (~2 days → ~17 Oct).
- MIXED-WINDING CURRENT (current_decomp.py; per-line current = mid-line crossings per unit length per 1000 t.u.):
  mixed (+6/−6): W16 +0.059/+0.034, W24 +0.057/+0.046, W40 +0.021/+0.027, W60 +0.018/+0.026 → positive in 8/8 active runs; passive 0.
  symmetric controls (+6/+6, W16–60, 8 runs): −0.007 … +0.022, mean +0.006, s.d. ≈ 0.011.
  → A steady topological current from the + wall to the − wall IS present, clearest in narrow channels (W16–24: ~0.05, 3–5 s.d. above the controls), marginal at W40–60 (~0.02). Layers (mixed): + wall keeps 6.4–7.6 of 12 (W16: 10.4), − wall 10.9–12.7 −½.
  W10: tracker finds no defects at all (also passive) → channel too narrow for the analysis; excluded.
- LONG CHANNELS (L 768, W 40): winding m 24 (charge 48): wall 27.7 + 28.1, bulk +20.4 (42 %, the wall law holds at 4× the length); dominant vortex wavenumber 14 per L (wavelength ≈ 55 cells ≈ 1.4 W). Lattice persistence not yet analysed. Homeotropic m 0: laminar, centre |v_x| 2.6 (compressibility concern, as at L 192). Current 0 in both (crossings +0.001).
- (12:40) QUEUED runs/channel_current (26 runs, 26 parallel, 6000 t.u. each): mixed +6/−6 at W16, W24 (s3–s8) and W32 (s1–4); symmetric +6/+6 controls W16, W24 (s3–6), W32 (s1–2). Purpose: topological current vs width with statistics.
- (13:30) REVIEW, ellipse checks done (bulkdyn/review_ellipse.py → fig_review_ellipse.png, review_ellipse.txt; shield_wl/ellipse_review.wl):
  (3) Tips are NOT better shielded. Mean |force| from the layer term (kBT/cell) inside the 0.8 level: k +6 ideal tips 2.2 vs centre/flanks 2.0; measured 7.2 vs 6.0; k −6 ideal 1.1 vs 1.2; measured 5.3 vs 5.3. The impression came from comparing ranges over regions of different size and from the colour scale. Measured +k layers are actually DEPLETED at the tips (3.4 wall +½ in |x| > 0.6 A vs 4.95 ideal; gaps 1.30 vs 1.13 site spacings), consistent with tip escapes.
  Flux coefficient (uniform-winding term with the charges present, as used, vs the imposed winding k): changes the numbers but not the picture: +k ideal inside 0.6 level 14 → 27 kBT, measured 176 → 227; −k ideal 4.9 → 27.5, measured 167 → 172. With the layer alone the two are not equivalent (the missing charge is carried by bulk defects); the 'charges present' version is the self-consistent Neumann problem and stays the default; quote both as a bracket.
  → All three review points resolved: (1) masking + self term, (2) angular jitter at tracker resolution, (3) no tip effect. The red TODO on the PRL shielding/vacancy paragraphs can be lifted after Benjamin's look.
- (14:20) LONG-CHANNEL VORTEX ROW (channel/vortex_row.py → fig_vortex_row.png, vortex_row.txt; L 768, W 40, mid-band vorticity every 10 t.u.):
  winding (charge 48): 38.5 ± 3.6 sign domains (≈19 vortex pairs), broad structure-factor peak at λ ≈ 59 cells = 1.5 W (half width 7 of 13 modes), domain-length CV 0.55, pattern correlation 0.88/0.40/0.24/0.11 after 10/50/100/200 t.u. (decorrelation ≈ 60 t.u.), no net drift (+0.06 cells/t.u.). → NOT a lattice: a fluctuating, liquid-like vortex row with a preferred spacing ≈ 1.5 W. Bulk +½: 1.31 per vortex cell, spread uniformly within the cell (not at vortex centres).
  homeotropic: laminar flow (|v_x| 2.3) advecting a vortex pattern at −2.5 cells/t.u. (wavelength 1.3 W, decorrelation ≈ 50 t.u. in the co-moving frame). (A first version limited the drift search to ±60 cells and wrongly showed fast decorrelation; fixed with FFT cross-correlation over all shifts.)
- (15:00) CURRENT STATISTICS (channel/current_stats.py → current_stats.txt, fig_current_stats.png; existing runs, 2 seeds per width; new 6000-t.u. runs pending):
  The motion+events estimator (Σ q v_y + pair dipoles, / bulk width) is far less noisy than mid-line crossings (controls ±0.002 vs ±0.012). Per-line current (per unit length per 1000 t.u.), motion+events:
  W16 mixed +0.044 ± 0.012 vs control 0.000 ± 0.001; W24 +0.048 ± 0.002 vs +0.005 ± 0.002; W40 +0.031 ± 0.005 vs −0.001 ± 0.002; W60 +0.022 ± 0.003 vs +0.003 ± 0.003.
  → the current is present at every width, ~0.045–0.05 for W ≤ 24, falling to ~0.02–0.03 at W 40–60. Carried by +½ ejected from the + wall (narrow: essentially all of it; W ≥ 40: about half, the rest by other +½); −½ and pair events contribute little.
- (15:00) SI movie candidates rendering (channel/cavity_movie.py, current_movie.py): S1 active circle (overcharged, tilted layer, vortex), S2 passive circle under an imposed vortex (ejection), S3 aspect-2 ellipse (leaking), + existing escape movie (energetics/movie/escape_eps2.0_s2_ev4.mp4) and channel movies; current movie with accumulated-current graph (W24 mixed vs control).
