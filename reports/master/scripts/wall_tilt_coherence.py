#!/usr/bin/env python3
"""Do the wall +1/2 of an ellipse all tilt the same way? Tilt = sin(psi - theta_n), psi = defect angle (defects.dat), theta_n = angle of
the outward wall normal at the nearest wall point. Wall = within 4 cells of the wall. Per defect frame (after 20%): coherence = |mean sign|
over wall defects (1: all the same way, 0: half and half); fraction tilted with the majority. Per flow interval: rotation sense of each half
(mean azimuthal velocity about the half's centroid) and the fraction of wall defects in that half whose tilt sign agrees with it.
Usage (from runs/): wall_tilt_coherence.py"""
import gzip, json, sys, numpy as np
sys.path.insert(0, "cavity")
from edge_criteria import read_with_angle
from wall_fourier import wall_polyline, project
runs = [("circle 1.43 la (k6)", "cavity/Rscan_R31.83_seeds/s02/k06", 50)]
runs += [(f"ellipse {e} s{s}", f"cassini/W3_ellipse/eps{e}_k06/s0{s}", 100) for e in ("1.2", "1.5", "2.0") for s in (1, 2)]
print(f"{'run':22s} {'coherence':>9} {'majority':>8} | {'tilt agrees with own half (by half rotation)':>44} | halves same sense")
for lab, r, dtf in runs:
    bc = json.load(open(r + "/input.json"))["BC"][0]; Q = np.array(bc["Q"][:2], float)
    W, s, P = wall_polyline(bc); nrm_ang = np.arctan2(*(np.gradient(W, axis=0)[:, ::-1] * [-1, 1]).T[::-1]) if False else None
    # outward normal of the polyline: rotate the tangent by -90 deg (polyline runs ccw)
    T = np.roll(W, -1, axis=0) - np.roll(W, 1, axis=0); n_ang = np.arctan2(-T[:, 0], T[:, 1])
    t, fr = read_with_angle(r, Q, 0.0); n0 = len(fr) // 5
    coh, maj, rec = [], [], []
    for j in range(n0, len(fr)):
        f = fr[j]; p = f[f[:, 2] > 0]
        if not len(p): continue
        d, sa = project(W, s, p[:, :2]); w = d < 4
        if w.sum() < 2: continue
        k = np.searchsorted(s, sa[w]) % len(W)
        sg = np.sign(np.sin(p[w, 3] - n_ang[k]))
        coh.append(abs(sg.mean())); maj.append(max(np.mean(sg > 0), np.mean(sg < 0)))
        rec.append((t[j], p[w, 0], sg))
    rows = [l.split() for l in gzip.open(r + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()]
    A = np.array(rows, float); ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]
    inside = (x / np.abs(W[:, 0]).max()) ** 2 + (y / np.abs(W[:, 1]).max()) ** 2 < 0.95
    def halfrot(m):
        cx, cy = x[m].mean(), y[m].mean(); dx, dy = x[m] - cx, y[m] - cy; rr = np.maximum(np.hypot(dx, dy), 0.5)
        return np.mean((dx * A[:, m, 5] - dy * A[:, m, 4]) / rr, axis=1)
    L, R = halfrot(inside & (x < 0)), halfrot(inside & (x > 0))
    agree = []
    for tt, xs, sg in rec:
        i = np.searchsorted(ts, tt)
        if i >= len(ts): continue
        own = np.where(xs < 0, np.sign(L[i]), np.sign(R[i])); agree.append(np.mean(sg == own))
    same = np.mean(np.sign(L[len(ts) // 5:]) == np.sign(R[len(ts) // 5:]))
    print(f"{lab:22s} {np.mean(coh):9.2f} {np.mean(maj):8.2f} | {np.mean(agree):44.2f} | {same:.0%}")
