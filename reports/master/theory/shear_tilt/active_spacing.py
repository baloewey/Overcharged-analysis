#!/usr/bin/env python3
"""Active continuum (Q-tensor + backflow + extensile activity, no imposed flow; active_bf.run) in the disk at different wall-layer spacings (user, 2026-10-16):
does the self-organised overcharged state show the MPCD tilt drop with crowding? R 31.83, lambda 2, zeta 0.3 and 1.0 (the values that gave the overcharged
state at k 6), k = 4, 6, 8, 10, 13, 16, 20 with the 2k-2 ring + bulk pair seeded. Last 6 observations (every 250 t.u. to 4000). -> active_spacing.json/.txt"""
import json, numpy as np
from multiprocessing import Pool
import shear_tilt as M
np.seterr(all="ignore")


def one(j):
    k, zeta = j; M.k = k; M.qb1 = 0.5 * M.S0 * np.cos(2 * k * M.PH); M.qb2 = 0.5 * M.S0 * np.sin(2 * k * M.PH)
    import active_bf as B
    obs, _ = B.run(2.0, zeta, 0.0, ring=2 * k - 2, T0=300, T1=4000, dt=0.02, every=250, nflow=1); last = obs[-6:]
    zs = [np.exp(1j * np.radians(o["tilt"])) for o in last if "tilt" in o]
    f = lambda key: float(np.nanmean([o.get(key, np.nan) for o in last]))
    return dict(k=k, zeta=zeta, spacing=2 * np.pi * (M.R - 2.6) / (2 * k - 2), nwall=f("nwall"), coh=f("coh"), circ=f("circ") if any("circ" in o for o in last) else np.nan,
                tilt=float(np.degrees(np.angle(np.mean(zs)))) if zs else np.nan, series=[(o.get("nwall"), o.get("tilt", np.nan)) for o in obs])


if __name__ == "__main__":
    J = [(k, z) for z in (0.3, 1.0) for k in (4, 6, 8, 10, 13, 16, 20)]
    with Pool(7) as p: R = p.map(one, J, chunksize=1)
    json.dump(R, open("active_spacing.json", "w"), default=float)
    L = ["ACTIVE CONTINUUM DISK, wall-layer spacing scan (R 31.83, lambda 2, Q-tensor + backflow + activity, no imposed flow): wall +1/2 vs 2k-2, coherence, circulation, tilt (deg)"]
    for r in R: L.append(f"zeta {r['zeta']:.1f} k {r['k']:2d} (spacing {r['spacing']:4.1f}): wall {r['nwall']:5.1f}/{2*r['k']-2} coh {r['coh']:.2f} circ {r['circ']:+.2f} tilt {r['tilt']:+6.1f} | last tilts " + " ".join(f"{t:.0f}" for _, t in r["series"][-6:]))
    open("active_spacing.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
