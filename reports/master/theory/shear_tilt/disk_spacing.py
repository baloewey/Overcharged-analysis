#!/usr/bin/env python3
"""Continuum circular cavity at different wall-layer spacings (user, 2026-10-16): does the passive Q-tensor disk with the measured no-slip flow show the MPCD
tilt drop with crowding? R = 31.83 (as the MPCD k-scan), winding k = 4, 6, 8, 10, 13, 16, 20 (ring of 2k-2 + bulk pair seeded; spacing 2 pi (R-2.6)/(2k-2) =
30.6 ... 5.1 cells), lambda 2, imposed vortex U 0.3 and 1.0. Two models: (a) one-constant (shear_tilt.run); (b) K1/K3 1.82 + MPCD Jeffery law lam_c(S)
(disk_mpcd.step). Tilt = coherent mean of the wall-defect tilts, last 4 observations (of every 250 t.u. up to 3000). -> disk_spacing.json/.txt"""
import json, numpy as np
from multiprocessing import Pool
import shear_tilt as M
np.seterr(all="ignore")
import os
RATIO = float(os.environ.get("RATIO", 1.82)); TAG = os.environ.get("TAG", "")


def setk(k):
    M.k = k; M.qb1 = 0.5 * M.S0 * np.cos(2 * k * M.PH); M.qb2 = 0.5 * M.S0 * np.sin(2 * k * M.PH)


def one(j):
    model, k, U = j; setk(k); ring = 2 * k - 2
    import disk_mpcd as DM, flatwall_aniso as A
    if model == "oneK":
        obs = M.run(2.0, U, "noslip", ring, T0=300.0, T1=3000.0)
    else:
        L1, L3, K1, K3 = A.params(RATIO); dt = 0.02
        zs = list((M.R - 2.6) * np.exp(2j * np.pi * np.arange(ring) / ring)) + [0.45j * M.R, -0.45j * M.R]; th2, _, _ = M.fields(M.X, M.Y, zs, [0.5] * (ring + 2), M.R)
        q1 = np.where(M.inside, 0.5 * M.S0 * np.cos(th2), M.qb1); q2 = np.where(M.inside, 0.5 * M.S0 * np.sin(th2), M.qb2)
        for _ in range(int(300 / dt)): q1, q2 = DM.step(q1, q2, dt, 2.0, None, L1, L3)
        fl = M.setup_flow(U, "noslip"); obs = []; t = 0.0
        while t < 3000 - 1e-9:
            for _ in range(int(250 / dt)): q1, q2 = DM.step(q1, q2, dt, 2.0, fl, L1, L3)
            t += 250; o = M.observe(q1, q2, 1); o["t"] = t; obs.append(o)
    last = [o for o in obs[-4:] if "tilt" in o]
    z = [np.exp(1j * np.radians(o["tilt"])) for o in last]
    return dict(model=model, k=k, U=U, spacing=2 * np.pi * (M.R - 2.6) / ring, nwall=float(np.mean([o["nwall"] for o in obs[-4:]])),
                coh=float(np.nanmean([o.get("coh", np.nan) for o in last])) if last else np.nan, tilt=float(np.degrees(np.angle(np.mean(z)))) if z else np.nan,
                tilt_series=[o.get("tilt", np.nan) for o in obs])


if __name__ == "__main__":
    J = [(m, k, U) for m in os.environ.get("MODELS", "oneK,aniso_mpcd").split(",") for k in (4, 6, 8, 10, 13, 16, 20) for U in (0.3, 1.0)]
    with Pool(8) as p: R = p.map(one, J, chunksize=1)
    json.dump(R, open(f"disk_spacing{TAG}.json", "w"), default=float)
    L = ["CONTINUUM DISK, wall-layer spacing scan (R 31.83, lambda 2, imposed no-slip vortex): wall +1/2 (target 2k-2), coherence, tilt (deg; arccos(1/2) = 60)"]
    for r in R: L.append(f"{r['model']:10s} k {r['k']:2d} (spacing {r['spacing']:4.1f}) U {r['U']:.1f}: wall {r['nwall']:5.1f}/{2*r['k']-2} coh {r['coh']:.2f} tilt {r['tilt']:+6.1f}  | series " + " ".join(f"{t:.0f}" for t in r["tilt_series"][-6:]))
    open(f"disk_spacing{TAG}.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
