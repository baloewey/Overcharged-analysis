#!/usr/bin/env python3
"""Tilt vs anchoring-site spacing a_s on a flat no-slip wall in shear (user, 2026-10-15: is the crowding drop of the MPCD tilt elastic?).
Hypothesis: with K1 != K3 the registry shift is pinned (tau -> rest value) with a spacing-independent barrier, while the flow's grip grows ~ s a_s^2,
so tau(a_s) = tau_L / (1 + (a*/a_s)^2)-like with a* ~ sqrt(dK / (eta s)); one-constant (no pinning) should show no such drop.
Scan: K1/K3 1.82 (MPCD) and 1.0 (control), lambda 2, shear s 0.003 / 0.01 / 0.03, a_s 6-40 (flatwall_aniso.run, H 30).  -> flatwall_spacing.json, .txt"""
import json, numpy as np
from multiprocessing import Pool
import flatwall_aniso as A
np.seterr(all="ignore")
JOBS = [(r, s, a) for r in (1.82, 1.0) for s in (0.003, 0.01, 0.03) for a in (6, 8, 10, 14, 20, 28, 40)]
def one(j):
    r, s, a = j
    try: o0, o1 = A.run(2.0, s, ratio=r, a_s=float(a)); return dict(ratio=r, s=s, a_s=a, rest=o0, **o1)
    except Exception as ex: return dict(ratio=r, s=s, a_s=a, err=repr(ex))
if __name__ == "__main__":
    with Pool(6) as p: R = p.map(one, JOBS, chunksize=1)
    json.dump(R, open("flatwall_spacing.json", "w"), default=float); L = ["FLAT WALL, lambda 2: tilt vs site spacing a_s (cells; h 0.5), arccos(1/2) = 60 deg"]
    for r in R:
        if "err" in r: L.append(f"K1/K3 {r['ratio']:.2f} s {r['s']} a_s {r['a_s']}: ERROR {r['err']}"); continue
        L.append(f"K1/K3 {r['ratio']:.2f} s {r['s']:.3f} a_s {r['a_s']:3d}: tilt {r.get('tilt', float('nan')):6.1f}  registry*360 {r.get('registry', float('nan')) * 360:6.1f}  2*far {2 * r['far_alpha']:6.1f}  depth {r.get('depth', float('nan')):.2f}"
                 f"  Er=s*a_s^2 {r['s'] * r['a_s'] ** 2:6.2f} | at rest: tilt {r['rest'].get('tilt', float('nan')):6.1f}")
    open("flatwall_spacing.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
