#!/usr/bin/env python3
"""Continuum band/bulk tilt vs MPCD's chi (run.sh). Per run, profile frames with t >= 500 and |circulation| > 0.3 (single coherent vortex):
  alpha(y) director angle from the tangent (deg, signed with the circulation) and uniformity, from the summed ring Z; v(y) azimuthal velocity;
  Lambda(y) = |v_r + v/r| / |v_r - v/r|; alpha_J(y) = local aligned angle with advection and chi: cos 2 alpha_J = [1 + (1/chi - 1)(1 - Lambda)] / lambda;
  plus layer tilt, coherence, wall +1/2, bulk +1/2 and -1/2 counts (pair creation = turbulence indicator), |circ|. theta_L = 30 deg (lambda 2).
-> summary.txt"""
import glob, json, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = 31.830989; LAM = 2.0; YC = 0.25 + 0.5 * np.arange(32); YS = (1.25, 2.25, 3.25, 4.25, 6.25, 8.25, 10.25)
L = ["CONTINUUM BAND / BULK TILT vs chi (circle R 31.83, lambda 2; theta_L = 30 deg, 2 theta_L = 60). Frames t >= 500 with |circ| > 0.3.",
     "run | frames | layer tilt (coh) | wall +1/2 | bulk +1/2, -1/2 | |circ| || alpha(y) at y = " + " ".join(f"{y:g}" for y in YS) + " || alpha_J(chi) at same y"]
for f in sorted(glob.glob(H + "/z*_k*_chi*_s*.json")):
    J = json.load(open(f)); chi = J["args"]["chi"]; ob = [o for o in J["obs"] if o["t"] >= 500 and abs(o.get("circ", 0)) > 0.3 and "pZr" in o]
    if not ob: L.append(f"{os.path.basename(f)[:-5]:22s} | no coherent frames"); continue
    Z = np.sum([np.array(o["pZr"]) + 1j * np.array(o["pZi"]) for o in ob], 0); S = np.sum([o["pS"] for o in ob], 0); N = np.sum([o["pN"] for o in ob], 0)
    V = np.sum([o["pV"] for o in ob], 0) / np.maximum(N, 1); al = np.degrees(0.5 * np.angle(Z)); un = np.abs(Z) / np.maximum(S, 1e-9)
    r = R - YC; dv = -np.gradient(V, YC); om = dv + V / r; e = 0.5 * (dv - V / r); Lam = np.abs(om) / np.maximum(np.abs(2 * e), 1e-12)
    c = (1 + (1 / chi - 1) * (1 - Lam)) / LAM; aJ = np.where(np.abs(c) <= 1, np.degrees(0.5 * np.arccos(np.clip(c, -1, 1))), np.nan)
    tl = [o["tilt"] for o in ob if "tilt" in o and o.get("coh", 0) > 0.5]; co = np.mean([o.get("coh", np.nan) for o in ob])
    g = lambda arr, y: np.interp(y, YC, arr)
    L.append(f"{os.path.basename(f)[:-5]:22s} | {len(ob):3d} | {np.mean(np.abs(tl)) if tl else np.nan:5.1f} ({co:.2f}) | {np.mean([o['nwall'] for o in ob]):4.1f} | "
             f"{np.mean([o['nbp'] for o in ob]):4.1f}, {np.mean([o['nbm'] for o in ob]):4.1f} | {np.mean([abs(o['circ']) for o in ob]):.2f} || "
             + " ".join(f"{g(al, y):+5.1f}({g(un, y):.2f})" for y in YS) + " || " + " ".join(f"{g(aJ, y):5.1f}" for y in YS))
open(os.path.join(H, "summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
