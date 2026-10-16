#!/usr/bin/env python3
"""Flat-wall rest state check (2026-10-15): is the 'at rest' tilt of flatwall_aniso a converged elastic minimum? (It disagrees with the circle pinning
calculation fluct/anisotropy_pinning, which puts the commensurate layer at tau = 0.) No flow; K1/K3 = 1.82, 1.0, 0.55; a_s 10/20/40; two initial
registries of the interior texture (phase 0 and a quarter period); tilt and registry recorded at t = 500, 1500, 4500, 9000. Note the far boundary is
zero-gradient (free far-field angle), so the registry is not fixed by a bounded bulk as in the circle.  -> flatwall_rest.txt"""
import numpy as np
from multiprocessing import Pool
import flatwall_aniso as A, flatwall as F
np.seterr(all="ignore")
def one(j):
    ratio, a_s, ph = j; h = A.h; S0 = A.S0; H = 30.0; dt = 0.02
    L1, L3, K1, K3 = A.params(ratio); nx, ny = int(round(a_s / h)), int(round(H / h)) + 4
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 2) * h; X, Y = np.meshgrid(x, y, indexing="ij"); O = A.ops(nx, ny)
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * (X - ph * a_s) / a_s) + (1 - w) * np.pi * ph; q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2):
        q1[:, :2] = qw1[:, None]; q2[:, :2] = qw2[:, None]; q1[:, -2:] = q1[:, -3:-2]; q2[:, -2:] = q2[:, -3:-2]; return q1, q2
    q1, q2 = bc(q1, q2); out = []; t = 0.0
    for T in (500, 1500, 4500, 9000):
        while t < T - 1e-9:
            H1, H2 = A.Hfield(q1, q2, L1, L3, O); q1, q2 = bc(q1 + dt * H1, q2 + dt * H2); t += dt
        m = F.measure(q1[:, 1:-1], q2[:, 1:-1], X[:, 1:-1], Y[:, 1:-1], a_s); E = A.energy(q1, q2, L1, L3, O)
        out.append((T, m.get("tilt", np.nan), m.get("registry", np.nan) * 360, m.get("far_alpha", np.nan), E))
    return (ratio, a_s, ph, out)
if __name__ == "__main__":
    J = [(r, a, p) for r in (1.82, 1.0, 0.55) for a in (10.0, 20.0, 40.0) for p in (0.0, 0.25)]
    with Pool(6) as pool: R = pool.map(one, J, chunksize=1)
    L = ["FLAT WALL AT REST: tilt (deg), registry*360, far-field angle (deg), energy, at t = 500 / 1500 / 4500 / 9000"]
    for ratio, a_s, ph, out in R:
        L.append(f"K1/K3 {ratio:4.2f} a_s {a_s:4.0f} init phase {ph:.2f}: " + " | ".join(f"t{T:.0f}: {tl:6.1f} {rg:6.1f} {fa:6.1f} E {E:.3f}" for T, tl, rg, fa, E in out))
    open("flatwall_rest.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
