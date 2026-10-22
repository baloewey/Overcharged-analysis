#!/usr/bin/env python3
"""Continuum version of locallock/local_lock.py (user, 2026-10-22): split the tilt of the continuum disks into twice the near-wall director angle
one cell inside the wall defects (ring average of the angle from the local tangent, signed by the circulation) and the excess.
Final fields of the active one-constant disks (fields_act_<zeta>_k04/16, zeta 3-8) and of the swirl controls (passive-like: zeta 3 + imposed swirl).
-> continuum_local_lock.txt"""
import glob, json, numpy as np, warnings
warnings.filterwarnings("ignore")
import shear_tilt as M
L = ["field | k | tilt | depth | 2 alpha_1 (local field) | excess | 2 alpha at depth 4-5 | circulation sign"]
for f in sorted(glob.glob("fields_act_*_k*.npz")) + sorted(glob.glob("fields_swirl_*_k*.npz")):
    k = 4 if "k04" in f else 16; M.k = k; F = np.load(f); q1, q2, vx, vy = F["q1"], F["q2"], F["vx"], F["vy"]; R = float(F["R"])
    D = M.defects(q1, q2); w = [(px, py, a) for px, py, q, a in D if q > 0 and np.hypot(px, py) > R - 6]
    if len(w) < 3: continue
    ph = np.arctan2(M.Y, M.X); r = np.hypot(M.X, M.Y); vt = -vx * np.sin(ph) + vy * np.cos(ph); circ = np.sign(np.mean(vt[M.inside & (r < R - 2)]))
    pw = np.array([np.arctan2(py, px) for px, py, _ in w]); tw = np.angle(np.mean(np.exp(1j * (np.array([a for _, _, a in w]) - pw))))
    dep = np.median([R - np.hypot(px, py) for px, py, _ in w])
    S2 = 2 * np.hypot(q1, q2); th2 = np.angle(q1 + 1j * q2); z = S2 * np.exp(1j * circ * (th2 - 2 * ph - np.pi))
    def alpha(u0, u1):
        m = M.inside & (R - r >= u0) & (R - r < u1); return np.degrees(np.angle(z[m].mean()) / 2)
    a1 = alpha(dep + 0.5, dep + 1.5); a45 = alpha(4, 5); tau = abs(np.degrees(tw))
    L.append(f"{f:32s} | {k:2d} | {tau:5.1f} | {dep:.2f} | {2*a1:5.1f} | {tau-2*a1:+5.1f} | {2*a45:5.1f} | {circ:+.0f}")
open("continuum_local_lock.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
