#!/usr/bin/env python3
"""Does the size of the melted core set the overtilt? (user, 2026-10-22). Flat no-slip wall in simple shear (flatwall.py model: one +1/2 per site spacing
a_s, anchoring rotating by pi per period, one-constant Q-tensor with L = Gamma = 1, S0 = 0.876, lambda 2, one-way flow). The bulk terms are scaled by b,
H = lap Q + b (Q - c tr(Q^2) Q), so the core (coherence) length xi = 1/sqrt(b) changes while K = 2 L S0^2 and S0 stay fixed: b = 4, 1, 0.25 -> xi = 0.5,
1, 2 cells. Shear s = 0.03 (l_E = 4.4 cells), a_s = 8 (crowded) and 20 (sparse). Tilt at t = 3000 and 6000 after a 300 t.u. rest start (convergence).
Overtilt = tilt - arccos(1/lambda) = tilt - 60 deg. -> flatwall_core.txt"""
import os, json, numpy as np
from multiprocessing import Pool
import flatwall as F
S0, h, c = F.S0, F.h, F.c
def run(job):
    b, a_s, s, lam = job; H = 40.0; dt = min(0.04, 0.04 / b)
    nx, ny = int(round(a_s / h)), int(round(H / h)) + 2; x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 1) * h; X, Y = np.meshgrid(x, y, indexing="ij")
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-Y / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2): q1[:, 0] = qw1; q2[:, 0] = qw2; q1[:, -1] = q1[:, -2]; q2[:, -1] = q2[:, -2]; return q1, q2
    lap = lambda f: (np.roll(f, 1, 0) + np.roll(f, -1, 0) - 2 * f) / h ** 2 + (np.roll(f, 1, 1) + np.roll(f, -1, 1) - 2 * f) / h ** 2
    gx = lambda f: (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
    vx = s * np.maximum(Y, 0); Exy = Oxy = 0.5 * s; lam_c = lam * S0
    def step(q1, q2, flow):
        tr = 2 * (q1 ** 2 + q2 ** 2); H1 = lap(q1) + b * (q1 - c * tr * q1); H2 = lap(q2) + b * (q2 - c * tr * q2)
        if flow: d1 = -vx * gx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * gx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
        else: d1, d2 = H1, H2
        return bc(q1 + dt * d1, q2 + dt * d2)
    q1, q2 = bc(q1, q2)
    for _ in range(int(300 / dt)): q1, q2 = step(q1, q2, False)
    out = []
    for T in (3000, 6000):
        for _ in range(int((T - (out[-1][0] if out else 0)) / dt)): q1, q2 = step(q1, q2, True)
        m = F.measure(q1, q2, X, Y, a_s); out.append((T, m.get("tilt", np.nan), m.get("depth", np.nan), m["ndef"]))
    S = 2 * np.hypot(q1, q2); core = (S < 0.5 * S0) & (Y > 0.3)
    return job, out, float(core.sum() * h * h)
if __name__ == "__main__":
    jobs = [(b, a_s, 0.03, 2.0) for b in (4.0, 1.0, 0.25) for a_s in (8.0, 20.0)]
    with Pool(6) as p: res = p.map(run, jobs)
    L = ["FLAT WALL, core size vs overtilt (lambda 2, s 0.03; tilt - 60 deg); xi = 1/sqrt(b) cells; columns: t 3000 | t 6000 (tilt, depth, defects); melted area (S < S0/2, cells^2)"]
    for (b, a_s, s, lam), out, area in res:
        L.append(f"xi {1/np.sqrt(b):.2f}  a_s {a_s:4.0f} | " + " | ".join(f"t {T}: tilt {tl:6.1f} (overtilt {tl-60:+5.1f}) depth {d:.2f} n {n}" for T, tl, d, n in out) + f" | melted area {area:.1f}")
    open("flatwall_core.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
