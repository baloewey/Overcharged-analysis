#!/usr/bin/env python3
"""Validation of seeded runs: defects detected at early frames vs the seeded list (match within 2.5 cells, same sign), wall count and bulk charge over time.
Usage: check_seed.py RUNDIR [...]"""
import os, sys, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/energetics")
import energetics as En  # noqa
from scipy.spatial import cKDTree
for run in sys.argv[1:]:
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]
    S = np.loadtxt(run + "/dirinit_seed.txt"); W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
    t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t)
    print(f"== {os.path.basename(run)}: seeded {len(S)} defects (wall +1/2 {int(((S[:, 2] > 0) & (tree.query(S[:, :2])[0] < 4)).sum())}, bulk charge {S[tree.query(S[:, :2])[0] >= 4, 2].sum():+.1f})")
    for n in [0, 1, 2, 5, 10, 20, 50, 100, 150, len(fr) - 1]:
        if n >= len(fr): continue
        f = fr[n]
        if not len(f): print(f"  t {t[n]:6.1f}: no defects"); continue
        d = tree.query(f[:, :2])[0]; nw = int(((f[:, 2] > 0) & (d < 4)).sum()); qb = 0.5 * np.sign(f[d >= 4, 2]).sum()
        matched = sum(1 for x, y, q in S if np.any((np.hypot(f[:, 0] - x, f[:, 1] - y) < 2.5) & (np.sign(f[:, 2]) == np.sign(q))))
        print(f"  t {t[n]:6.1f}: {len(f):2d} defects, total charge {0.5 * np.sign(f[:, 2]).sum():+.1f}, wall +1/2 {nw:2d}, bulk {qb:+.1f}; seeds matched within 2.5 cells {matched}/{len(S)}")
