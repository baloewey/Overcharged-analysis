#!/usr/bin/env python3
"""Flat wall with a MELTED wall (user, 2026-10-16): Dirichlet anchoring of the director direction (strong, as MPCD: wall deviation 2-5 deg) but
reduced order at the wall, Q_w = (S_w/2)(cos 2 theta_w, sin 2 theta_w), S_w < S0 (MPCD: S 0.91 -> 0.79 between defects, 0.65 at cores when crowded).
Scan via env SW (wall order), otherwise as flatwall_depth. Original docstring:
Flat-wall depth test (user, 2026-10-16): does the tilt depend on the defect depth through d / l_E, l_E = sqrt(K/(gamma1 s sqrt(lambda^2-1))) the
alignment length? (MPCD: tilt linear in d sqrt(s), tiltcrowd/depth_correction.txt.) Q-tensor flat wall (flatwall_aniso.run machinery; K/gamma1 = 1 here),
lambda 2, K1/K3 1.0 and 1.82, shear s 0.003/0.01/0.03/0.1 (l_E 14.4/7.6/4.4/2.4 cells), site spacing a_s 8/12/20/30 (sets the depth).
Convergence: the sheared state is recorded at t = 3000, 6000 and 12000 after a 500 t.u. rest start. -> flatwall_depth.json/.txt"""
import json, os, numpy as np
from multiprocessing import Pool
import flatwall_aniso as A, flatwall as F
np.seterr(all="ignore")
LAM = 2.0
def one(j):
    ratio, s, a_s, SW = j; h = A.h; S0 = A.S0; H = 30.0; dt = float(os.environ.get('DT', 0.02))
    L1, L3, K1, K3 = A.params(ratio); nx, ny = int(round(a_s / h)), int(round(H / h)) + 4
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 2) * h; X, Y = np.meshgrid(x, y, indexing="ij"); O = A.ops(nx, ny)
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * SW * np.cos(2 * thw), 0.5 * SW * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2):
        q1[:, :2] = qw1[:, None]; q2[:, :2] = qw2[:, None]; q1[:, -2:] = q1[:, -3:-2]; q2[:, -2:] = q2[:, -3:-2]; return q1, q2
    vx = s * np.maximum(Y, 0); Oxy = 0.5 * s; Exy = 0.5 * s; dx = O[0]; lam_c = LAM * S0
    q1, q2 = bc(q1, q2); t = 0.0
    while t < 500:
        H1, H2 = A.Hfield(q1, q2, L1, L3, O); q1, q2 = bc(q1 + dt * H1, q2 + dt * H2); t += dt
    out = []; t = 0.0
    for T in [float(v) for v in os.environ.get('TIMES', '3000,6000,12000').split(',')]:
        while t < T - 1e-9:
            H1, H2 = A.Hfield(q1, q2, L1, L3, O)
            d1 = -vx * dx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * dx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
            q1, q2 = bc(q1 + dt * d1, q2 + dt * d2); t += dt
        m = F.measure(q1[:, 1:-1], q2[:, 1:-1], X[:, 1:-1], Y[:, 1:-1], a_s); out.append((T, m.get("tilt", np.nan), m.get("depth", np.nan), m.get("far_alpha", np.nan)))
    lE = np.sqrt(1.0 / (s * np.sqrt(LAM ** 2 - 1)))
    return dict(ratio=ratio, s=s, a_s=a_s, SW=SW, lE=lE, out=out)
if __name__ == "__main__":
    SS = [float(v) for v in os.environ.get('SS', '0.003,0.01,0.03,0.1').split(',')]; TAG = os.environ.get('TAG', '')
    SWS = [float(v) for v in os.environ.get('SWS', '0.876').split(',')]
    J = [(r, s, a, sw) for r in (1.0, 1.82) for s in SS for sw in SWS for a in (8.0, 12.0, 20.0, 30.0)]
    with Pool(8) as pool: R = pool.map(one, J, chunksize=1)
    json.dump(R, open(f"flatwall_depth{TAG}.json", "w"), default=float)
    L = ["FLAT WALL DEPTH TEST (lambda 2; K/gamma1 = 1): tilt (deg), depth d, d/l_E at t = 3000 / 6000 / 12000 (convergence)"]
    for r in R:
        L.append(f"K1/K3 {r['ratio']:.2f} S_w {r['SW']:.2f} s {r['s']:.3f} a_s {r['a_s']:4.0f} l_E {r['lE']:5.2f}: " + " | ".join(f"t{T:.0f}: tilt {tl:6.1f} d {d:.2f} d/lE {d / r['lE']:.2f}" for T, tl, d, fa in r["out"]))
    open(f"flatwall_depth{TAG}.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
