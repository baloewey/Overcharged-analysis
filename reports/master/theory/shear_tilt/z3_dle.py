#!/usr/bin/env python3
"""d/l_E of the continuum active disk at zeta 3 (user, 2026-10-19): from the final fields (fields_act_3_k04/16_fields.npz): wall +1/2 depth d = R - r (median),
wall shear s = d<v_theta>/dr - <v_theta>/r (azimuthal average, 1-4 cells from the wall; same estimator spirit as MPCD e_band + w_band), and
l_E = sqrt(K/(gamma1 s sqrt(lambda^2 - 1))) with the model's K = 1.535, gamma1 = 1.53, lambda 2. Compared with MPCD (d/l_E 0.5-1.2; tiltcrowd) and the
continuum at zeta 0.3 (flat-wall/disk estimates 0.05-0.29). -> z3_dle.txt"""
import sys, json, numpy as np
import shear_tilt as M
K, G1, LAM = 1.535, 1.53, 2.0
L = []
for k in (4, 16):
    F = np.load(f"fields_act_3_k{k:02d}_fields.npz"); q1, q2, vx, vy = F["q1"], F["q2"], F["vx"], F["vy"]; R = float(F["R"])
    M.k = k; D = M.defects(q1, q2); wall = [(px, py) for px, py, q, a in D if q > 0 and np.hypot(px, py) > R - 6]
    d = np.array([R - np.hypot(px, py) for px, py in wall])
    r = np.hypot(M.X, M.Y); ph = np.arctan2(M.Y, M.X); vt = -vx * np.sin(ph) + vy * np.cos(ph); ins = M.inside
    edges = np.arange(R - 8, R + 0.01, 0.5); rc = 0.5 * (edges[1:] + edges[:-1]); prof = np.array([vt[ins & (r >= a) & (r < b)].mean() for a, b in zip(edges[:-1], edges[1:])])
    band = (R - rc >= 1) & (R - rc <= 4); sl = np.polyfit(rc[band], prof[band], 1)[0]; s = abs(sl - np.mean(prof[band] / rc[band]))
    lE = np.sqrt(K / (G1 * s * np.sqrt(LAM ** 2 - 1)))
    J = json.load(open(f"var_act_3_k{k:02d}_fields.json"))
    L.append(f"k {k:2d}: wall +1/2 {len(wall)} (json nwall {J['nwall']:.1f}); depth d median {np.median(d):.2f} (IQR {np.percentile(d,25):.2f}-{np.percentile(d,75):.2f}) cells; wall shear s {s:.4f}; "
             f"l_E {lE:.2f}; d/l_E {np.median(d)/lE:.2f}; tilt {J['tilt']:.1f} (first run {json.load(open(f'var_act_3_k{k:02d}.json'))['tilt']:.1f}); v_theta at 2 cells {np.interp(R-2, rc, prof):.3f}")
open("z3_dle.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
