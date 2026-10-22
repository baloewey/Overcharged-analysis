#!/usr/bin/env python3
"""Summaries of test 1a (depth_memory.npz) and test 2 (wallflow.npz) by wall density level and shape class. -> analyse_12.txt"""
import csv, numpy as np
S = "/home/baloewe/Physics/Overcharged/reports/master/theory/shapes"; SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
RHO = (0.04, 0.06, 0.09, 0.12, 0.14); L = []
def levels(runs):
    rho = np.array([float(SR[r]["rho"]) for r in runs]); asp = np.array([float(SR[r]["aspect"]) for r in runs])
    return np.argmin(np.abs(rho[:, None] - np.array(RHO)), 1), asp
# ---- test 1a
Z = np.load("depth_memory.npz", allow_pickle=True); lev, asp = levels(Z["runs"]); lags = Z["LAGS"]
L.append("TEST 1a (i): autocorrelation of resident depth and cos tilt (normalised, pooled runs), lags in t.u.")
for nm, a in (("depth", Z["acd"]), ("cos tilt", Z["acc"])):
    for i, r0 in enumerate(RHO):
        m = lev == i; c = a[m].mean(0); c = c / c[0]; tc = lags[np.argmax(c < np.exp(-1))] if (c < np.exp(-1)).any() else np.nan
        L.append(f"  {nm:8s} rho {r0:.2f}: " + " ".join(f"{lags[q]:.0f}:{c[q]:.2f}" for q in (1, 2, 5, 10, 20, 30)) + f" | 1/e time {tc:.0f} t.u.")
L.append("\nTEST 1a (ii): window-averaged (20 t.u.) depth and cos tilt -> escape. Within-run MH odds ratio (upper vs lower half), prompt (+20 t.u.) and lagged (+40 in wall, +60 out)")
def mh(rows_by_run, xcol, ycol):
    num = den = 0.0
    for S_ in rows_by_run:
        if len(S_) < 40: continue
        x = S_[:, xcol]; y = S_[:, ycol] > 0.5; hi = x > np.median(x); n = len(S_)
        num += (y & hi).sum() * (~y & ~hi).sum() / n; den += (~y & hi).sum() * (y & ~hi).sum() / n
    return num / den if den else np.nan
for i, r0 in enumerate(list(RHO) + ["all"]):
    rs = [s for s, l in zip(Z["S"], lev) if (r0 == "all" or l == i) and len(s)]
    L.append(f"  rho {r0}: window depth prompt {mh(rs, 0, 5):.2f} lagged {mh(rs, 0, 4):.2f} | window cos prompt {mh(rs, 1, 5):.2f} lagged {mh(rs, 1, 4):.2f} | "
             f"instant depth prompt {mh(rs, 2, 5):.2f} lagged {mh(rs, 2, 4):.2f}")
# ---- test 2
W = np.load("wallflow.npz"); lev, asp = levels(W["runs"]); DB = W["DB"]; dc = 0.5 * (DB[:-1] + DB[1:]); prof = W["prof"]
names = ["mean inward v", "rms inward v", "mean |v_t|", "rms |v|", "signed v_t"]
for nm, sel in (("circles", asp < 1.05), ("aspect>=1.4", asp >= 1.4)):
    L.append(f"\nTEST 2: near-wall flow, {nm}: medians over runs (MPCD units); columns wall distance d = " + " ".join(f"{x:.2f}" for x in dc[[1, 3, 5, 7, 11, 15, 19, 23]]))
    for q in (1, 2, 3, 0):
        for i, r0 in enumerate(RHO):
            m = (lev == i) & sel; p = np.median(prof[m, q], 0)
            L.append(f"  {names[q]:13s} rho {r0:.2f}: " + " ".join(f"{p[j]:+.4f}" for j in (1, 3, 5, 7, 11, 15, 19, 23)))
    for i, r0 in enumerate(RHO):
        m = (lev == i) & sel; vt = np.median(prof[m, 2], 0); sh = (vt[np.argmin(np.abs(dc - 3.25))] - vt[np.argmin(np.abs(dc - 0.75))]) / 2.5
        L.append(f"  near-wall shear d|v_t|/dd (0.75 -> 3.25 cells) rho {r0:.2f}: {sh:.4f}; rms inward at 2-6 cells {np.median(prof[m, 1][:, (dc > 2) & (dc < 6)].mean(1)):.4f}")
open("analyse_12.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
