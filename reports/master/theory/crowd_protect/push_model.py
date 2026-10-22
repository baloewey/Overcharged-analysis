#!/usr/bin/env python3
"""Test 1b (2026-10-22): can the per-defect push model explain crowding protection? Push P = 0.40 cos tau g(d) (g = inward jet transmitted, jetmap_ff12),
logistic P(exit in 20 t.u.) = 1/(1 + exp(-(a + b P))), a -3.66, b +6.25 (energetics/perdef_ff36.txt, k 6 elongated runs). Applied to every resident wall +1/2
start in the shape sweep (shapes/escape_events/events.npz; tilt from the outward normal, same convention as perdef.py) and averaged per wall density;
compared with the measured 20 t.u. escape probability. Only ratios between densities are compared (the escape definitions differ)."""
import csv, numpy as np
S = "/home/baloewe/Physics/Overcharged/reports/master/theory/shapes"
GX = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 6]); GY = np.array([0.30, 0.33, 0.58, 0.86, 0.97, 1.0]); A_, B_ = -3.66, 6.25
E = np.load(S + "/escape_events/events.npz")["E"]; run = E[:, 0].astype(int); ct = np.cos(np.radians(E[:, 2])); d = E[:, 5]; esc = E[:, 7]
meta = list(csv.DictReader(open(S + "/escape_events/runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
rho = np.array([float(SR[m["run"]]["rho"]) for m in meta]); asp = np.array([float(SR[m["run"]]["aspect"]) for m in meta])
lev = np.argmin(np.abs(rho[:, None] - np.array([0.04, 0.06, 0.09, 0.12, 0.14])), 1)
P = 0.40 * ct * np.interp(d, GX, GY); pred = 1 / (1 + np.exp(-(A_ + B_ * P)))
L = ["PUSH MODEL vs CROWDING: per density level, resident wall +1/2 means and the model's predicted exit probability vs the measured one (20 t.u.)",
     "subset | rho | mean depth | <cos tau> | <push> | predicted p | measured p | pred/pred(0.04) | meas/meas(0.04)"]
for nm, sel in (("circles", asp < 1.05), ("aspect>=1.4", asp >= 1.4), ("all", asp > 0)):
    base = None
    for i, r0 in enumerate((0.04, 0.06, 0.09, 0.12, 0.14)):
        m = np.isin(run, np.flatnonzero((lev == i) & sel)); pp, mm = pred[m].mean(), esc[m].mean()
        if base is None: base = (pp, mm)
        L.append(f"{nm:11s} | {r0:.2f} | {d[m].mean():.2f} | {ct[m].mean():+.2f} | {P[m].mean():+.3f} | {pp:.4f} | {mm:.4f} | {pp / base[0]:.2f} | {mm / base[1]:.3f}")
open("push_model.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
