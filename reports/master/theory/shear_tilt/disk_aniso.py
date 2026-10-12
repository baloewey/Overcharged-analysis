#!/usr/bin/env python3
"""Disk Q-tensor model with the measured no-slip flow (shear_tilt.py) plus splay != bend (cubic L3 term of flatwall_aniso.py), so that both the local
(core/layer) and the global (bulk +1 spiral) effects of K1 != K3 are present. Mean K kept at 2 S^2. Ring of 12 relaxed without flow, then flow.
Scan: K1/K3 = 1 and 1.82 (MPCD), lam 1.2 / 2 / 4, U 0.3 / 1.0. -> disk_aniso.txt"""
import json, numpy as np
from multiprocessing import Pool
import shear_tilt as M, flatwall_aniso as A
np.seterr(all="ignore")
h = M.h
def dx(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
def dy(f): return (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
def H_aniso(q1, q2, L1, L3):
    tr = 2 * (q1 ** 2 + q2 ** 2); g1x, g1y, g2x, g2y = dx(q1), dy(q1), dx(q2), dy(q2)
    Gxx = g1x ** 2 + g2x ** 2; Gyy = g1y ** 2 + g2y ** 2; Gxy = g1x * g1y + g2x * g2y
    d1 = dx(q1 * g1x + q2 * g1y) + dy(q2 * g1x - q1 * g1y); d2 = dx(q1 * g2x + q2 * g2y) + dy(q2 * g2x - q1 * g2y)
    H1 = L1 * M.lap(q1) + q1 - M.c * tr * q1 - 0.5 * L3 * (Gxx - Gyy) + L3 * d1
    H2 = L1 * M.lap(q2) + q2 - M.c * tr * q2 - L3 * Gxy + L3 * d2
    return H1, H2
def step(q1, q2, dt, lam_c, fl, L1, L3):
    H1, H2 = H_aniso(q1, q2, L1, L3)
    if fl is None: d1, d2 = H1, H2
    else:
        vx, vy, Exx, Exy, Oxy = fl; g1x, g1y = dx(q1), dy(q1); g2x, g2y = dx(q2), dy(q2)
        d1 = -(vx * g1x + vy * g1y) + lam_c * Exx + 2 * Oxy * q2 + H1; d2 = -(vx * g2x + vy * g2y) + lam_c * Exy - 2 * Oxy * q1 + H2
    q1 = q1 + dt * d1; q2 = q2 + dt * d2; q1[M.ghost] = M.qb1[M.ghost]; q2[M.ghost] = M.qb2[M.ghost]
    out = ~(M.inside | M.ghost); q1[out] = 0; q2[out] = 0; return q1, q2
def run(j, T0=600.0, T1=3000.0, dt=0.02, every=250.0):
    ratio, lam, U = j; L1, L3, K1, K3 = A.params(ratio); ring = 12
    zs = list((M.R - 2.6) * np.exp(2j * np.pi * (np.arange(ring) + 0.5) / ring)); th2, _, _ = M.fields(M.X, M.Y, zs, [0.5] * ring, M.R)
    q1 = np.where(M.inside, 0.5 * M.S0 * np.cos(th2), M.qb1); q2 = np.where(M.inside, 0.5 * M.S0 * np.sin(th2), M.qb2)
    for _ in range(int(T0 / dt)): q1, q2 = step(q1, q2, dt, 0, None, L1, L3)
    fl = M.setup_flow(U, "noslip"); obs = []; t = 0.0
    while t < T1 - 1e-9:
        for _ in range(int(every / dt)): q1, q2 = step(q1, q2, dt, lam * M.S0, fl, L1, L3)
        t += every; o = M.observe(q1, q2, 1); o["t"] = t; obs.append(o)
    last = obs[-4:]; z = [np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]
    return dict(ratio=ratio, lam=lam, U=U, nwall=float(np.mean([o["nwall"] for o in last])), nbp=float(np.mean([o.get("nbulk_p", 0) for o in last])),
                coh=float(np.nanmean([o.get("coh", np.nan) for o in last])), tilt=float(np.degrees(np.angle(np.mean(z)))) if z else np.nan, obs=obs)
if __name__ == "__main__":
    JOBS = [(r, l, U) for r in (1.0, 1.82) for l in (1.2, 2.0, 4.0) for U in (0.3, 1.0)]
    with Pool(6) as p: R = p.map(run, JOBS, chunksize=1)
    json.dump(R, open("disk_aniso.json", "w"), default=float); L = []
    for r in R: L.append(f"K1/K3 {r['ratio']:4.2f} lam {r['lam']:3.1f} U {r['U']:.1f}: wall {r['nwall']:4.1f} bulk +{r['nbp']:.1f}  coh {r['coh']:.2f}  tilt {r['tilt']:+6.1f}")
    open("disk_aniso.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
