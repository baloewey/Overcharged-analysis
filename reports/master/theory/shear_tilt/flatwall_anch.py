#!/usr/bin/env python3
"""Flat wall with FINITE anchoring (user, 2026-10-16): the continuum models clamp Q at the wall (Dirichlet, infinite anchoring, S = S0); MPCD anchors
particle by particle. Here the wall is a Rapini-Papoular surface term: ghost rows copy the first interior row (no elastic flux through the wall) and the first
interior row relaxes towards the anchored Q_w with rate Gamma W / h, i.e. F_s = (W/2) int |Q - Q_w|^2 dx (extrapolation length ~ L/W cells; L = 1).
Scan: K1/K3 1.0 and 1.82, lambda 2, shear s 0.03, site spacing a_s 8/12/20/30, W 1, 0.3, 0.1 (and the Dirichlet reference from flatwall_depth.txt).
Recorded at t = 2000 / 4000 / 8000 (convergence). Also the director deviation from the anchoring on the first interior row (median |delta|).
-> flatwall_anch.json/.txt"""
import json, os, numpy as np
from multiprocessing import Pool
import flatwall_aniso as A, flatwall as F
np.seterr(all="ignore")
LAM = 2.0
def one(j):
    ratio, s, a_s, W = j; h = A.h; S0 = A.S0; H = 30.0; dt = 0.02
    L1, L3, K1, K3 = A.params(ratio); nx, ny = int(round(a_s / h)), int(round(H / h)) + 4
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 2) * h; X, Y = np.meshgrid(x, y, indexing="ij"); O = A.ops(nx, ny)
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2):   # Neumann at the wall (ghost rows copy row 2) and at the top
        q1[:, 0] = q1[:, 2]; q1[:, 1] = q1[:, 2]; q2[:, 0] = q2[:, 2]; q2[:, 1] = q2[:, 2]
        q1[:, -2:] = q1[:, -3:-2]; q2[:, -2:] = q2[:, -3:-2]; return q1, q2
    vx = s * np.maximum(Y, 0); Oxy = 0.5 * s; Exy = 0.5 * s; dx = O[0]; lam_c = LAM * S0; kap = W / h
    def surf(q1, q2):
        d1 = np.zeros_like(q1); d2 = np.zeros_like(q2); d1[:, 2] = -kap * (q1[:, 2] - qw1); d2[:, 2] = -kap * (q2[:, 2] - qw2); return d1, d2
    q1, q2 = bc(q1, q2); t = 0.0
    while t < 500:
        H1, H2 = A.Hfield(q1, q2, L1, L3, O); s1, s2 = surf(q1, q2); q1, q2 = bc(q1 + dt * (H1 + s1), q2 + dt * (H2 + s2)); t += dt
    out = []; t = 0.0
    for T in (2000.0, 4000.0, 8000.0):
        while t < T - 1e-9:
            H1, H2 = A.Hfield(q1, q2, L1, L3, O); s1, s2 = surf(q1, q2)
            d1 = -vx * dx(q1) + 2 * Oxy * q2 + H1 + s1; d2 = -vx * dx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2 + s2
            q1, q2 = bc(q1 + dt * d1, q2 + dt * d2); t += dt
        m = F.measure(q1[:, 1:-1], q2[:, 1:-1], X[:, 1:-1], Y[:, 1:-1], a_s)
        dev = np.degrees(np.abs(np.angle(np.exp(1j * (np.arctan2(q2[:, 2], q1[:, 2]) - 2 * thw))) / 2)); Sw = 2 * np.hypot(q1[:, 2], q2[:, 2])
        out.append((T, m.get("tilt", np.nan), m.get("depth", np.nan), float(np.median(dev)), float(np.median(Sw)), m.get("ndef", 0)))
    return dict(ratio=ratio, s=s, a_s=a_s, W=W, out=out)
if __name__ == "__main__":
    J = [(r, 0.03, a, W) for r in (1.0, 1.82) for W in (1.0, 0.3, 0.1) for a in (8.0, 12.0, 20.0, 30.0)]
    with Pool(8) as pool: R = pool.map(one, J, chunksize=1)
    json.dump(R, open("flatwall_anch.json", "w"), default=float)
    L = ["FLAT WALL, FINITE ANCHORING (Rapini-Papoular W; lambda 2, s 0.03): tilt (deg), depth, median wall-director deviation from the anchoring (deg), wall S, ndef; t = 2000 / 4000 / 8000"]
    for r in R:
        L.append(f"K1/K3 {r['ratio']:.2f} W {r['W']:.2f} a_s {r['a_s']:4.0f}: " + " | ".join(f"t{T:.0f}: tilt {tl:6.1f} d {d:.2f} dev {dv:4.1f} Sw {sw:.2f} n {nd}" for T, tl, d, dv, sw, nd in r["out"]))
    open("flatwall_anch.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
