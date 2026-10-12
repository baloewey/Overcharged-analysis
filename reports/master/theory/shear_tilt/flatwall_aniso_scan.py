import json, numpy as np
from multiprocessing import Pool
import flatwall_aniso as A
np.seterr(all="ignore")
JOBS = [(r, l, s) for r in (1.0, 1.82, 1 / 1.82) for l in (1.2, 2.0, 4.0) for s in (0.001, 0.01, 0.03)]
def one(j):
    r, l, s = j; o0, o1 = A.run(l, s, ratio=r); return dict(ratio=r, lam=l, s=s, rest=o0, **o1)
if __name__ == "__main__":
    with Pool(7) as p: R = p.map(one, JOBS, chunksize=1)
    json.dump(R, open("flatwall_aniso_scan.json", "w"), default=float); L = []
    for r in R:
        t0 = np.degrees(np.arccos(1 / r["lam"]))
        L.append(f"K1/K3 {r['ratio']:4.2f} lam {r['lam']:3.1f} s {r['s']:.3f}: tilt {r.get('tilt', np.nan):6.1f} (arccos(1/lam) {t0:5.1f}, excess {r.get('tilt', np.nan) - t0:+5.1f})  "
                 f"2*far {2 * r['far_alpha']:6.1f}  depth {r.get('depth', np.nan):.2f}  registry*360 {360 * r.get('registry', np.nan):6.1f}  | at rest: tilt {r['rest'].get('tilt', np.nan):6.1f} "
                 f"registry*360 {360 * r['rest'].get('registry', np.nan):6.1f}")
    open("flatwall_aniso_scan.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
