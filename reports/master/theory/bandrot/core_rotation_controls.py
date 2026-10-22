#!/usr/bin/env python3
"""No-flow controls for core_rotation.py (2026-10-22): the estimator should give ~0 without flow. Continuum flat wall at rest (s = 0, a_s 8, 12, 20);
MPCD seeded passive layers without flow (runs/rest_pinning, director frames t <= 250 while the seeded layer is intact). -> core_rotation_controls.txt"""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore"); sys.path.insert(0, "."); import core_rotation as CR
from multiprocessing import Pool
def mpcd_early(run):
    import json, gzip, os
    from edge_criteria import read_with_angle
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
    with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
    td = np.unique(D[:, 0]); td = td[td <= 250]; t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
    n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij"); ins = np.hypot(X, Y) < R
    vals = [[], []]
    for tt in td:
        i = int(np.argmin(np.abs(t - tt)))
        if abs(t[i] - tt) > 1e-6: continue
        f = fr[i]; rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
        if w.sum() != N or abs(f[:, 2].sum() - k) > 1e-6: continue
        tau = np.angle(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))); sg = np.sign(tau) or 1.0
        d = D[D[:, 0] == tt]; th2 = np.zeros((n, m)); S = np.zeros((n, m)); ix, iy = d[:, 1].astype(int), d[:, 2].astype(int)
        th2[ix, iy] = 2 * np.arctan2(d[:, 5], d[:, 4]); S[ix, iy] = d[:, 7]; th2s, _, _ = CR.TS.fields(X, Y, list(f[:, 0] + 1j * f[:, 1]), list(f[:, 2]), R)
        for x0, y0 in zip(f[w, 0], f[w, 1]):
            valid = ins.copy()
            for x1, y1 in zip(f[:, 0], f[:, 1]):
                if (x1, y1) != (x0, y0): valid &= np.hypot(X - x1, Y - y1) > 1.5
            for j, v in enumerate(CR.per_defect(th2, S, th2s, X, Y, x0, y0, valid, 1.0)): vals[j].append(sg * v)
    return run.split("/")[-1], len(vals[0]), np.degrees(np.nanmedian(vals[0])) if vals[0] else np.nan, np.degrees(np.nanmedian(vals[1])) if vals[1] else np.nan
if __name__ == "__main__":
    import glob
    with Pool(6) as p:
        C = p.map(CR.continuum, [(a_s, 0.0, 2.0) for a_s in (8.0, 12.0, 20.0)])
        M = p.map(mpcd_early, sorted(glob.glob("/home/baloewe/Physics/Overcharged/runs/rest_pinning/k*_t*")))
    L = ["NO-FLOW CONTROLS for the core-rotation estimator (deg; ring 2-3 | ring 3-4.5). Signed by the (seeded/rest) tilt; tilt 0 runs signed +."]
    L.append("CONTINUUM flat wall at rest: a_s | tilt | depth | ring 2-3 | ring 3-4.5")
    for x in C: L.append(f"{x['job'][0]:4.0f} | {x['tilt']:6.1f} | {x['depth']:.2f} | {x['r23']:+6.1f} | {x['r345']:+6.1f}")
    L.append("MPCD seeded passive layers, no flow (t <= 250): run | samples | ring 2-3 | ring 3-4.5")
    for nm, n, a, b in M: L.append(f"{nm} | {n:4d} | {a:+6.1f} | {b:+6.1f}")
    open("core_rotation_controls.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
