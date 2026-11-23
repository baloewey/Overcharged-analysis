#!/usr/bin/env python3
"""Errors-in-variables check of the response slope (coreflow_response.txt). The MPCD flow at a core (disk r < 2, 2 t.u. average) carries thermal noise that
attenuates the OLS slope of v_n on u_n by var_signal / var_total, more when the signal is small (crowded). Noise variance per run from the flow itself:
u_n at a FIXED wall point (cells within r < 2 of points 2 cells inside the wall, every 40th wall point): noise is white between 2-t.u. frames (thermal
velocities decorrelate within a few collision times), the signal is not -> var_noise ~ var(u) - cov(u_t, u_t+2 t.u.) (upper bound if the signal also
decorrelates over 2 t.u.: then the corrected slope is an over-correction). -> noise_check.txt"""
import glob, json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coreflow as C
out = ["NOISE CHECK of the MPCD flow at wall points (r < 2 disk, 2 cells inside the wall): var_total, lag-1 (2 t.u.) covariance, noise var = difference"]
for k in ("05", "14"):
    vt, vn = [], []
    for run in sorted(glob.glob(C.R + f"/coreflow/ell2_k{k}_s*")):
        J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = np.array(bc["Q"][:2], float); W, _, _ = C.wall_polyline(bc); w = C.Wall(W)
        tf, UX, UY = C.mpcd_flow(run, *J["domain"]); UX, UY = UX[len(tf) // 10:], UY[len(tf) // 10:]
        X, Y = np.meshgrid(np.arange(J["domain"][0]) + 0.5 - Q[0], np.arange(J["domain"][1]) + 0.5 - Q[1], indexing="ij"); IN = w.inside(X, Y)
        for j in range(0, len(W), 40):
            p = W[j] + 2 * w.n[j]; c = (np.hypot(X - p[0], Y - p[1]) < 2) & IN
            u = (UX[:, c].mean(1) * w.n[j, 0] + UY[:, c].mean(1) * w.n[j, 1]); u2 = (u[:-1] + u[1:]) / 2      # same t +- 2 average as coreflow.py
            vt.append(u2.var()); a = u2 - u2.mean(); vn.append(u2.var() - np.mean(a[:-2] * a[2:]))         # lag 2 frames: the t+-2 averages do not overlap
    out.append(f"k{k}: var_total {np.mean(vt):.5f} (sd {np.sqrt(np.mean(vt)):.4f}) | noise var {np.mean(vn):.5f} (sd {np.sqrt(np.mean(vn)):.4f}) | noise fraction {np.mean(vn) / np.mean(vt):.2f}")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "noise_check.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out))
