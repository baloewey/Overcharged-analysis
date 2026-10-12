#!/usr/bin/env python3
"""One calculation with all three tilt contributions: disk Q-tensor model with the measured no-slip flow, K1/K3 = 1.82 (disk_aniso), and MPCD's Jeffery alignment law:
the strain coupling lam_c(S) E with lam_c(S) = lam (1 + S4(S)) / 2 point by point (local S; S4 from a Maier-Saupe distribution), so that the director alignment
parameter is lam_eff(S) = lam (1 + S4)/(2 S) as in MPCD (reduced order near the cores -> stronger alignment, automatically). lam = MPCD tumble parameter.
Compare: arccos(1/lam_eff(S_bulk)) and the MPCD passive-swirl tilts (61-64 / 75-83 / ~57 deg at lam 2 / 4 / 1.2). -> disk_mpcd.txt"""
import json, numpy as np
from multiprocessing import Pool
from scipy.special import i0e, i1e
from scipy.optimize import brentq
import shear_tilt as M, flatwall_aniso as A, disk_aniso as D
np.seterr(all="ignore")
Sg = np.linspace(0.002, 0.995, 400); S4g = np.array([1 - 2 * s / brentq(lambda x: i1e(x) / i0e(x) - s, 1e-8, 5000) for s in Sg])
def lamc(S, lam): return lam * (1 + np.interp(np.clip(S, Sg[0], Sg[-1]), Sg, S4g)) / 2
def step(q1, q2, dt, lam, fl, L1, L3):
    H1, H2 = D.H_aniso(q1, q2, L1, L3)
    if fl is None: d1, d2 = H1, H2
    else:
        vx, vy, Exx, Exy, Oxy = fl; lc = lamc(2 * np.hypot(q1, q2), lam); g1x, g1y = D.dx(q1), D.dy(q1); g2x, g2y = D.dx(q2), D.dy(q2)
        d1 = -(vx * g1x + vy * g1y) + lc * Exx + 2 * Oxy * q2 + H1; d2 = -(vx * g2x + vy * g2y) + lc * Exy - 2 * Oxy * q1 + H2
    q1 = q1 + dt * d1; q2 = q2 + dt * d2; q1[M.ghost] = M.qb1[M.ghost]; q2[M.ghost] = M.qb2[M.ghost]
    out = ~(M.inside | M.ghost); q1[out] = 0; q2[out] = 0; return q1, q2
def run(j, T0=600.0, T1=3000.0, dt=0.02, every=250.0):
    ratio, lam, U = j; L1, L3, K1, K3 = A.params(ratio); ring = 12
    zs = list((M.R - 2.6) * np.exp(2j * np.pi * (np.arange(ring) + 0.5) / ring)); th2, _, _ = M.fields(M.X, M.Y, zs, [0.5] * ring, M.R)
    q1 = np.where(M.inside, 0.5 * M.S0 * np.cos(th2), M.qb1); q2 = np.where(M.inside, 0.5 * M.S0 * np.sin(th2), M.qb2)
    for _ in range(int(T0 / dt)): q1, q2 = step(q1, q2, dt, lam, None, L1, L3)
    fl = M.setup_flow(U, "noslip"); obs = []; t = 0.0
    while t < T1 - 1e-9:
        for _ in range(int(every / dt)): q1, q2 = step(q1, q2, dt, lam, fl, L1, L3)
        t += every; o = M.observe(q1, q2, 1); o["t"] = t; obs.append(o)
    last = obs[-4:]; z = [np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]
    Sb = np.mean(2 * np.hypot(q1, q2)[M.inside & (M.Rr < M.R - 8)])
    return dict(ratio=ratio, lam=lam, U=U, nwall=float(np.mean([o["nwall"] for o in last])), coh=float(np.nanmean([o.get("coh", np.nan) for o in last])),
                tilt=float(np.degrees(np.angle(np.mean(z)))) if z else np.nan, S_bulk=float(Sb), lam_eff_bulk=float(lamc(Sb, lam) / Sb), obs=obs)
if __name__ == "__main__":
    JOBS = [(r, l, U) for r in (1.82, 1.0) for l in (1.2, 2.0, 4.0) for U in (0.3, 1.0)]
    with Pool(4) as p: R = p.map(run, JOBS, chunksize=1)
    json.dump(R, open("disk_mpcd.json", "w"), default=float); L = []
    meas = {1.2: "~57 (weakly coherent)", 2.0: "61-64", 4.0: "75-83"}
    for r in R:
        lead = np.degrees(np.arccos(1 / r["lam_eff_bulk"])) if r["lam_eff_bulk"] > 1 else np.nan
        L.append(f"K1/K3 {r['ratio']:4.2f} MPCD lam {r['lam']:3.1f} U {r['U']:.1f}: wall {r['nwall']:4.1f} coh {r['coh']:.2f} tilt {r['tilt']:+6.1f}   "
                 f"[bulk S {r['S_bulk']:.3f}, lam_eff {r['lam_eff_bulk']:.2f}, arccos(1/lam_eff) {lead:5.1f}; MPCD measured {meas[r['lam']]}]")
    open("disk_mpcd.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
