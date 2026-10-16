#!/usr/bin/env python3
"""Passive imposed-vortex crowding test (2026-10-15): wall +1/2 count, bulk charge, tilt coherence C and tilt of the coherent layer (C > 0.9) vs wall
spacing, second half of each run. Compare with the active circles at the same R (tiltcrowd/bulk_spiral.txt: 64.4 deg at k 3 -> 49 deg at k 20).
-> summary.txt"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = ["PASSIVE + IMPOSED VORTEX (R 31.83, swirl 0.01 within 0.6R, lambda 2), second half",
     "run | k | wall +1/2 | spacing 2 pi R/N_wall | bulk charge | median C | frac C>0.9 | tilt of coherent layer (deg) | tilt all frames"]
for run in sorted(glob.glob(H + "/k*_s*") + glob.glob(H + "/A03_k*")):
    if not os.path.exists(run + "/DONE"): continue
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
    nw, qb, C, tc, ta = [], [], [], [], []
    for f in fr:
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4); nw.append(w.sum()); qb.append(0.5 * np.sign(f[R - r >= 4, 2]).sum())
        if w.sum() < 3: continue
        z = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))); C.append(abs(z)); ta.append(abs(np.angle(z)))
        if abs(z) > 0.9: tc.append(abs(np.angle(z)))
    N = np.mean(nw); L.append(f"{os.path.basename(run):7s} | {k:2d} | {N:5.2f} | {2*np.pi*(R-1.5)/max(N,1):5.1f} | {np.mean(qb):+.2f} | {np.median(C):.2f} | {np.mean(np.array(C)>0.9):.2f} | "
                             f"{np.degrees(np.median(tc)) if len(tc) > 20 else float('nan'):5.1f} | {np.degrees(np.median(ta)):5.1f}")
print("\n".join(L)); open(os.path.join(H, "summary.txt"), "w").write("\n".join(L) + "\n")
