#!/usr/bin/env python3
"""Resting state of seeded passive rings (no flow): layer tilt (signed coherent mean, deg; 0 = radial arm at the wall = tails on the homeotropic sites,
180 = the other site type), coherence, wall +1/2 count, bulk charge, in time windows. -> summary.txt"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); W = [(0, 20), (20, 100), (100, 300), (300, 600), (600, 1000)]
L = ["RESTING STATE of seeded passive 2k-2 rings (R 31.83, no flow): per window t [0,20) [20,100) [100,300) [300,600) [600,1000): signed layer tilt (C) / wall +1/2 / bulk charge"]
for run in sorted(glob.glob(H + "/k*_t*")):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
    cells = []
    for lo, hi in W:
        zs, nw, qb = [], [], []
        for i in np.flatnonzero((t >= lo) & (t < hi)):
            f = fr[i]
            if not len(f): continue
            r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4); nw.append(w.sum()); qb.append(0.5 * np.sign(f[R - r >= 4, 2]).sum())
            if w.sum() >= 3: zs.append(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))))
        z = np.mean(zs) if zs else np.nan
        cells.append(f"{np.degrees(np.angle(z)):+6.1f} ({np.mean(np.abs(zs)) if zs else np.nan:.2f}) / {np.mean(nw) if nw else np.nan:4.1f} / {np.mean(qb) if qb else np.nan:+.1f}")
    L.append(f"{os.path.basename(run)} (k {k}, ring {2*k-2}, spacing {2*np.pi*R/(2*k-2):4.1f}): " + " | ".join(cells))
print("\n".join(L)); open(os.path.join(H, "summary.txt"), "w").write("\n".join(L) + "\n")
