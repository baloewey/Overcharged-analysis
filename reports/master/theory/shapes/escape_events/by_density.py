#!/usr/bin/env python3
"""Within-run gap and tilt effects on the lagged escape outcome, split by the run's wall density rho = 2k/P (shape_runs.csv); within-run gap spread (CV)."""
import csv, numpy as np
Z = np.load("events.npz"); E = Z["E"]; run = E[:, 0].astype(int); tilt = np.abs(E[:, 2]); gap = 0.5 * (E[:, 3] + E[:, 4]); late = E[:, 9] > 0.5
meta = list(csv.DictReader(open("runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open("../shape_runs.csv"))}
rho_r = np.array([float(SR[m["run"]]["rho"]) for m in meta]); lev = np.argmin(np.abs(rho_r[:, None] - np.array([0.04, 0.06, 0.09, 0.12, 0.14])), 1)
print("OUTCOME lagged escape (in wall at +40 t.u., deeper than 6 cells at +60 t.u.). Within-run quartile contrast: p(top gap quartile)/p(bottom), pooled per level (MH-style).")
print("rho  | runs | p_late x1000 | median gap | within-run gap CV | p ratio top/bottom gap quartile | p ratio top/bottom |tilt| quartile")
for L, r0 in enumerate((0.04, 0.06, 0.09, 0.12, 0.14)):
    cv, nh, dh, nt, dt_, P = [], 0, 0, 0, 0, []
    for j in np.flatnonzero(lev == L):
        m = run == j
        if m.sum() < 100: continue
        g, t, e = gap[m], tilt[m], late[m]; cv.append(g.std() / g.mean()); P.append(e.mean())
        qg = np.quantile(g, [0.25, 0.75]); qt = np.quantile(t, [0.25, 0.75])
        nh += e[g >= qg[1]].mean() * m.sum(); dh += e[g <= qg[0]].mean() * m.sum()
        nt += e[t >= qt[1]].mean() * m.sum(); dt_ += e[t <= qt[0]].mean() * m.sum()
    sel = np.isin(run, np.flatnonzero(lev == L))
    print(f"{r0:.2f} | {len(cv):4d} | {1000 * late[sel].mean():8.2f}     | {np.median(gap[sel]):6.1f}     | {np.median(cv):.2f}              | {nh / max(dh, 1e-12):.2f}                            | {nt / max(dt_, 1e-12):.2f}")
