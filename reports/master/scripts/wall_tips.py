#!/usr/bin/env python3
"""
Tip vs flank distribution of defects in ellipses (master document). Frames after 20%.
Tips = wall points with |x| > A/sqrt(2) (A = semi-major axis); wall defects = within 4 cells of the wall, assigned to the
nearest wall point. Reported: wall +1/2 per 100 cells of wall at the tips and on the flanks and their ratio (uniform
arclength = 1; conductor, i.e. fully screening, = eps); net charge per frame in the cavity regions |x| > A/sqrt(2) and the rest.
Run from runs/:  python3 ../reports/master/scripts/wall_tips.py > ../reports/master/data/wall_tips.txt
"""
import sys, json, numpy as np
sys.path.insert(0, "cavity")
from edge_criteria import read_with_angle
from wall_fourier import wall_polyline, project
BAND = 4.0
runs = [(f"passive eps{e}", f"cassini/P5_passive_refs/ellipse_eps{e}_k06") for e in ("1.2", "1.5", "2.0")]
runs += [(f"active eps{e} s{s}", f"cassini/W3_ellipse/eps{e}_k06/s0{s}") for e in ("1.2", "1.5", "2.0") for s in (1, 2)]
print("run\twall+1/2_tips_per100\twall+1/2_flanks_per100\tratio\tnetQ_tipregions\tnetQ_rest")
for lab, r in runs:
    bc = json.load(open(r + "/input.json"))["BC"][0]
    W, s, P = wall_polyline(bc); A = np.abs(W[:, 0]).max()
    t, fr = read_with_angle(r, np.array(bc["Q"][:2], float), 0.2); nf = len(fr)
    tipw = np.abs(W[:, 0]) > A / np.sqrt(2)
    seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T); Lt, Lf = seg[tipw].sum(), seg[~tipw].sum()
    wt = wf = qt = qf = 0.0
    for f in fr:
        xy = f[:, :2]; q = np.sign(f[:, 2]) * 0.5
        d, sa = project(W, s, xy); tip = tipw[np.searchsorted(s, sa) % len(W)]; wall = (d < BAND) & (q > 0)
        wt += np.sum(wall & tip); wf += np.sum(wall & ~tip)
        rt = np.abs(xy[:, 0]) > A / np.sqrt(2); qt += q[rt].sum(); qf += q[~rt].sum()
    dt, df = wt / nf / Lt * 100, wf / nf / Lf * 100
    print(f"{lab}\t{dt:.2f}\t{df:.2f}\t{dt / df:.2f}\t{qt / nf:.2f}\t{qf / nf:.2f}")
