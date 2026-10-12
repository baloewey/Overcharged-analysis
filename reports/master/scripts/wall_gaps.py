#!/usr/bin/env python3
"""
Tip vs flank spacing of wall defects in ellipses (replaces the count-based tip/flank ratio, which is dominated by the discreteness of
12 defects). Per frame: wall-sign defects within 4 cells of the wall, arclength s of the nearest wall point, sorted; gaps between
neighbours (cyclic); each gap is assigned to tip (its midpoint's wall point has |x| > A/sqrt(2)) or flank. Density ratio
= <gap>_flank / <gap>_tip: 1 for uniform spacing, eps for the conductor (fully screening) distribution.
Run from runs/: python3 ../reports/master/scripts/wall_gaps.py
"""
import sys, json, numpy as np
sys.path.insert(0, "cavity")
from edge_criteria import read_with_angle
from wall_fourier import wall_polyline, project
runs = [(f"passive eps{e}", f"cassini/P5_passive_refs/ellipse_eps{e}_k06") for e in ("1.2", "1.5", "2.0")]
runs += [(f"active eps{e} s{s}", f"cassini/W3_ellipse/eps{e}_k06/s0{s}") for e in ("1.2", "1.5", "2.0") for s in (1, 2)]
print("run\tgap_tip\tgap_flank\tdensity_ratio(flank_gap/tip_gap)\tframes")
for lab, r in runs:
    bc = json.load(open(r + "/input.json"))["BC"][0]; sg = 1 if bc["customAnch"] > 0 else -1
    W, s, P = wall_polyline(bc); A = np.abs(W[:, 0]).max(); tipw = np.abs(W[:, 0]) > A / np.sqrt(2)
    t, fr = read_with_angle(r, np.array(bc["Q"][:2], float), 0.2)
    gt, gf = [], []
    for f in fr:
        p = f[np.sign(f[:, 2]) == sg, :2]
        if len(p) < 3: continue
        d, sa = project(W, s, p); sa = np.sort(sa[d < 4.0])
        if len(sa) < 3: continue
        gaps = np.diff(np.r_[sa, sa[0] + P]); mids = (sa + gaps / 2) % P
        tip = tipw[np.searchsorted(s, mids) % len(W)]
        gt += list(gaps[tip]); gf += list(gaps[~tip])
    print(f"{lab}\t{np.mean(gt):.2f}\t{np.mean(gf):.2f}\t{np.mean(gf) / np.mean(gt):.3f}\t{len(fr)}")
