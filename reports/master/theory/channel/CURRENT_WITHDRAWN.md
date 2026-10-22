# The mid-line "topological current" in mixed-winding channels is a tracking artefact (2026-10-22)

**Claim withdrawn.** Mixed channels (+6/−6) carry a net topological current J ≈ 0.03–0.06 per unit length per 1000 t.u. that falls as 1/W; symmetric controls carry ≈ 0.
Sources: current_stats.py/.txt, current_decomp.py, current_movie.py, crowd_an.py "current crossings"; briefing_20261021 (channels); paper_outline "topological current".

**Why it is wrong.** The mid-line of a periodic channel is a closed loop. The net charge above it, Q_up(t), changes only when charge crosses it; pair creation and annihilation are neutral, and a pair born straddling the line counts correctly. Q_up is therefore an exact, tracking-free measure of the net charge carried across. Its slope is zero in every run:

| Width W | Current from Q_up | Previously reported (tracked) |
|---|---|---|
| 16–60, mixed and symmetric | \|J\| ≤ 0.001 | 0.011–0.064 in mixed channels |

Over each run Q_up only fluctuates within a band of about 4 charge units (crossing2/qup_current.txt).

**Where the old count went wrong.** It linked each defect to its nearest same-sign neighbour within 2.5 cells per 2 t.u. frame.
- 18% of bulk defects move more than 2.5 cells in a frame, and 1% move more than 13.6 cells.
- Their tracks break (median track length 4 t.u.). Fast trips, which are mostly downward +½, go uncounted.
- The uncounted crossings are not random, so the count is biased, not just noisy.

**The fix.** A predictive tracker (crossing2/tracker2.py) uses the last displacement as a prediction, a global Hungarian assignment, a 9-cell gate and up to 2 missed frames. It gives a median track length of 46 t.u. and reproduces Q_up(t): end value +1.0 vs +0.9, rms difference 2 (W24_mixed_s3, 6000 t.u.).

**What is real.** There is heavy traffic across the mid-line. In W24 mixed over 6000 t.u.:
- 537 +½ crossings up and 477 down;
- 101 −½ crossings up and 71 down.

The +½ carry a net +30 charge up. This is cancelled by −½ moving up (−15) and by pairs born across the line with the −½ above (−14). Figure: crossing2/fig_heights.png.

**To redo** for the channel paper: any quantity built on the old crossing count, including the J ∝ 1/W law, the decomposition into wall-ejected and other +½, and the crowd_an current columns.
