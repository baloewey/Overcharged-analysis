#!/usr/bin/env python3
"""Fixed-spacing shear scan vs the d sqrt(s) ordering (2026-10-19). Passive imposed-vortex runs: runs/swirl_shear (k 6 and 13, swirl A 0.005-0.1, R 31.83) measured
exactly as depth_correction.py (defect_flow.one: depth d, coherent-layer tilt; wallshear.one: wall shear s = e_band + w_band). Reference fit tilt = tau0 + B d sqrt(s)
on the earlier passive crowding set (swirl_crowd, 14 runs, as depth_correction.py); prediction for the shear scan.
Question: at FIXED spacing, does changing the shear move the tilt along the same line (d/l_E ordering), or does spacing act separately?
-> shear_scan_dsqrts.txt"""
import os, sys, glob, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); ROOT = "/home/baloewe/Physics/Overcharged"
import depth_correction as DC
if __name__ == "__main__":
    ref = sorted(glob.glob(ROOT + "/runs/swirl_crowd/k0[468]_s*") + glob.glob(ROOT + "/runs/swirl_crowd/A03_k*_s*")); ref = [r for r in ref if "A03_k16" not in r]
    scan = sorted(glob.glob(ROOT + "/runs/swirl_shear/k*_A*"))
    with Pool(6) as p: R1 = [r for r in p.map(DC.one, ref) if r]; R2 = p.map(DC.one, scan)
    x = np.array([r["d"] * np.sqrt(r["s"]) for r in R1]); y = np.array([r["tilt"] for r in R1]); B, t0 = np.polyfit(x, y, 1)
    L = [f"Reference (swirl_crowd passive, {len(R1)} runs): tilt = {t0:.1f} + {B:.1f} d sqrt(s); rms {np.std(y - (t0 + B * x)):.1f} deg",
         "shear scan run | k | depth d | wall shear s | d sqrt(s) | tilt (coherent, deg) | predicted | residual"]
    for run, r in zip(scan, R2):
        if not r or not np.isfinite(r["tilt"]): L.append(f"{os.path.basename(run):11s} | no coherent layer / failed"); continue
        xx = r["d"] * np.sqrt(r["s"]); pr = t0 + B * xx
        L.append(f"{os.path.basename(run):11s} | {r['k']:2d} | {r['d']:.2f} | {r['s']:.3f} | {xx:.3f} | {r['tilt']:5.1f} | {pr:5.1f} | {r['tilt'] - pr:+5.1f}")
    open(os.path.join(H, "shear_scan_dsqrts.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
