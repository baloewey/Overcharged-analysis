#!/usr/bin/env python3
"""Independent test of the annulus theory (annulus_b2.py): the chi scan (runs/chiscan circ_k*_chi*, R 31.83, k 4/8/12, 3 seeds) changes MPCD's shear
susceptibility chi at fixed size. Theory: near-wall director deficit alpha - theta_L comes from chi != 1 plus elasticity; at chi = 1 it should vanish (up to
the interior pull). Profiles with curvature_profiles.one; annulus with the run's chi, D 0.73, lambda 2, y_L 1.25, y_c 8.25. -> chi_test.txt"""
import os, glob, numpy as np
from multiprocessing import Pool
import curvature_profiles as CP, annulus_b2 as AB
H = os.path.dirname(os.path.abspath(__file__)); YS = [1.25, 2.25, 3.25, 4.25, 5.25, 6.25]
runs = sorted(glob.glob("/home/baloewe/Physics/Overcharged/runs/chiscan/circ_k*_chi*_s*"))
with Pool(9) as p: res = p.map(CP.one, runs)
L = ["CHI-SCAN TEST of the annulus theory (R 31.83, lambda 2; alpha from the tangent, deg; mean over seeds)",
     "k | chi | runs | row | y: " + " ".join(f"{y:5.2f}" for y in YS)]
for k in (4, 8, 12):
    for chi in (0.25, 0.5, 1.0):
        g = [r for r, run in zip(res, runs) if r and "err" not in r and f"_k{k:02d}_" in run and f"_chi{chi}_" in run]
        if not g: continue
        v = np.nanmean([r["v"] for r in g], 0); am = np.nanmean([r["alpha"] for r in g], 0); R = g[0]["R"]
        s, _ = AB.solve(R, v, am, chi, 1.25, 8.25); s1, _ = AB.solve(R, v, am, 1.0, 1.25, 8.25)
        L.append(f"{k:2d} | {chi:4.2f} | {len(g)} | measured      " + " ".join(f"{np.interp(y, CP.YC, am):+5.1f}" for y in YS))
        L.append(f"   |      |   | alpha_J(chi)  " + " ".join(f"{AB.alphaJ(R, v, chi, y):+5.1f}" for y in YS))
        L.append(f"   |      |   | annulus chi   " + " ".join(f"{np.degrees(s.sol(R - y)[0]):+5.1f}" for y in YS) + ("" if s.success else "  (bvp FAIL)"))
        L.append(f"   |      |   | annulus chi=1 " + " ".join(f"{np.degrees(s1.sol(R - y)[0]):+5.1f}" for y in YS))
open(os.path.join(H, "chi_test.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
