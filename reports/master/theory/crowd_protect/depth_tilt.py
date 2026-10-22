#!/usr/bin/env python3
"""Does depth set the tilt, defect by defect? (user list item 3, 2026-10-22). depth_memory.npz: for every resident wall +1/2 (in the wall for the previous
20 t.u.), the 20 t.u. window mean of depth and of cos(tilt) (tilt from the outward normal; smaller cos = more tilted). Within each run: Spearman rank correlation
between window depth and window |tilt| (via -cos), and the slope of mean |tilt| (deg) per cell of depth across within-run depth quintiles. Pooled per rho level
and shape class (median over runs). Prediction from the d sqrt(s) ordering: deeper -> more tilted (positive correlation). -> depth_tilt.txt"""
import csv, numpy as np
Z = np.load("depth_memory.npz", allow_pickle=True); SR = {r["run"]: r for r in csv.DictReader(open("../shapes/shape_runs.csv"))}
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); rk = lambda x: np.argsort(np.argsort(x))
L = ["DEPTH vs TILT within runs (20 t.u. window means; median over runs): Spearman rho(depth, |tilt|) | slope of |tilt| vs depth (deg/cell, quintile means) | runs"]
for cls, lo, hi in (("circles", 0.99, 1.01), ("aspect>=1.4", 1.39, 9)):
    for i, r0 in enumerate(LV):
        cs, sl = [], []
        for r, S in zip(Z["runs"], Z["S"]):
            s_ = SR[r]
            if not (lo <= float(s_["aspect"]) <= hi and np.argmin(abs(float(s_["rho"]) - LV)) == i) or len(S) < 200: continue
            d = S[:, 0]; tl = np.degrees(np.arccos(np.clip(S[:, 1], -1, 1))); ok = np.isfinite(d) & np.isfinite(tl); d, tl = d[ok], tl[ok]
            cs.append(np.corrcoef(rk(d), rk(tl))[0, 1]); q = np.quantile(d, np.linspace(0, 1, 6)); b = np.clip(np.digitize(d, q[1:-1]), 0, 4)
            md = [d[b == j].mean() for j in range(5)]; mt = [tl[b == j].mean() for j in range(5)]; sl.append(np.polyfit(md, mt, 1)[0])
        L.append(f"{cls:11s} rho {r0:.2f}: {np.median(cs):+.2f} | {np.median(sl):+5.1f} | {len(cs)}")
open("depth_tilt.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
