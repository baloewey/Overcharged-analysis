"""Scan for shear_tilt.py -> scan.json, scan.txt. Steady values = mean over the last 4 samples (t >= T1 - 750)."""
import json, os, numpy as np
from multiprocessing import Pool
import shear_tilt as M
JOBS = [("lam", l, U, "noslip", 12) for l in (0.5, 1.2, 2.0, 4.0) for U in (0.1, 0.3, 1.0)] + \
       [("U", 2.0, U, "noslip", 12) for U in (0.01, 0.03, 0.06, 0.2, 0.6, 2.0)] + \
       [("variant", l, U, v, ring) for l in (0.5, 2.0, 4.0) for U in (0.3, 1.0) for v in ("noslip", "slip", "wallonly") for ring in (12, 10)]
def one(j):
    tag, l, U, v, ring = j; obs = M.run(l, U, v, ring, T0=300, T1=2000, every=250); last = obs[-4:]
    f = lambda key: float(np.nanmean([o.get(key, np.nan) for o in last]))
    tl = np.exp(1j * np.radians([o.get("tilt", np.nan) for o in last]))
    return dict(tag=tag, lam=l, U=U, variant=v, ring=ring, nwall=f("nwall"), nbulk_p=f("nbulk_p"), nbulk_m=f("nbulk_m"), coh=f("coh"),
                tilt=float(np.degrees(np.angle(np.nanmean(tl)))), obs=obs)
if __name__ == "__main__":
    with Pool(6) as p: res = p.map(one, JOBS, chunksize=1)
    json.dump(res, open("scan.json", "w")); L = []
    for r in res:
        L.append(f"{r['tag']:8s} lam {r['lam']:3.1f} U {r['U']:5.2f} {r['variant']:8s} start {r['ring']:2d}: wall {r['nwall']:4.1f} bulk +{r['nbulk_p']:.1f}/-{r['nbulk_m']:.1f} "
                 f"coh {r['coh']:.2f} tilt {r['tilt']:+6.1f}")
    open("scan.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
