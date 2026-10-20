#!/usr/bin/env python3
"""SI S3 numbers (2026-10-20): active circles, +k, rho <= 0.15, k 2-10, act 0.07, lambda 2, 10 seeds (runs/cavity/analysis_weekend/main_sweep.csv: N_wall = +1/2 within
4 cells of the wall; criteria_all.csv: q_knee, has_knee, q_dens (bulk charge from the density-minimum edge)). Per R/l_a bin: runs, mean wall count minus (2k-2), fraction of
runs with a knee, mean knee bulk charge, mean q_dens; plus tilt coherence / circulation from the flow table if present. -> circles_size.txt"""
import csv, numpy as np
LA = 22.21; A = "/home/baloewe/Physics/Overcharged/runs/cavity/analysis_weekend/"
M = {r["run"]: r for r in csv.DictReader(open(A + "main_sweep.csv")) if r["set"] == "main_act0.07"}
C = {r["run"]: r for r in csv.DictReader(open(A + "criteria_all.csv")) if r["set"] == "pos"}
f = lambda x: float(x) if x not in ("", "nan", None) else np.nan
rows = []
for run, m in M.items():
    k = int(m["k"]); rho = f(m["rho"])
    if not (2 <= k <= 10 and rho <= 0.151): continue
    c = C.get(run, {}); rows.append(dict(k=k, R=f(m["R"]) / LA, dw=f(m["N_wall"]) - (2 * k - 2), knee=f(c.get("q_knee")), has=f(c.get("has_knee")), qd=f(c.get("q_dens"))))
bins = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.2), (2.2, 2.6), (2.6, 3.0)]
L = ["R/l_a bin | runs | wall +1/2 minus (2k-2): mean (sd) | fraction with a knee | knee bulk charge (mean over runs with a knee) | q_dens bulk charge"]
for lo, hi in bins:
    g = [r for r in rows if lo <= r["R"] < hi]
    if not g: continue
    dw = np.array([r["dw"] for r in g]); kn = np.array([r["knee"] for r in g]); hs = np.array([r["has"] for r in g]); qd = np.array([r["qd"] for r in g])
    L.append(f"{lo:.1f}-{hi:.1f} | {len(g):4d} | {np.nanmean(dw):+.2f} ({np.nanstd(dw):.2f}) | {np.nanmean(hs):.2f} | {np.nanmean(kn):.2f} | {np.nanmean(qd):.2f}")
w = [r for r in rows if 0.8 <= r["R"] <= 2.2]
L.append(f"\nwindow 0.8 <= R/l_a <= 2.2: runs {len(w)}; occupancy N_wall/(2k-2) {np.nanmean([1 + r['dw'] / (2 * r['k'] - 2) for r in w]):.3f}; knee bulk {np.nanmean([r['knee'] for r in w]):.2f} (knee in {np.nanmean([r['has'] for r in w]):.2f} of runs)")
for k in range(2, 11):
    g = sorted([r for r in rows if r["k"] == k], key=lambda r: r["R"])
    full = [r["R"] for r in g if abs(r["dw"]) < 0.25]
    L.append(f"k {k:2d}: R/l_a range in the sweep {g[0]['R']:.2f}-{g[-1]['R']:.2f}; wall count within 0.25 of 2k-2 for R/l_a {min(full) if full else float('nan'):.2f}-{max(full) if full else float('nan'):.2f} (runs)")
open("circles_size.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
