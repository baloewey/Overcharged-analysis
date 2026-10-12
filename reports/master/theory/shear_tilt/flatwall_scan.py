import json, numpy as np
from multiprocessing import Pool
import flatwall as F
np.seterr(all="ignore")
JOBS = [(l, s, 20.0, 30.0) for l in (1.1, 1.2, 1.5, 2.0, 3.0, 4.0) for s in (0.001, 0.003, 0.01, 0.03, 0.1)] + [(2.0, s, a, 30.0) for a in (10.0, 40.0) for s in (0.003, 0.03)] + \
       [(2.0, 0.03, 20.0, 60.0)]
def one(j):
    l, s, a, H = j; o0, o1 = F.run(l, s, a_s=a, H=H, T0=1500, T1=3000); return dict(lam=l, s=s, a_s=a, H=H, **o1)
if __name__ == "__main__":
    with Pool(5) as p: R = p.map(one, JOBS, chunksize=1)
    json.dump(R, open("flatwall_scan.json", "w"))
    for r in R:
        th = np.degrees(np.arccos(1 / r["lam"])) / 2
        print(f"lam {r['lam']:3.1f} s {r['s']:.3f} a_s {r['a_s']:4.0f} H {r['H']:3.0f}: tilt {r.get('tilt', np.nan):6.1f}  2*far {2 * r['far_alpha']:6.1f}  2*theta_L {2 * th:5.1f}  "
              f"2*alpha(depth) {2 * r.get('alpha_at_depth', np.nan):6.1f}  depth {r.get('depth', np.nan):.2f}  registry*360 {360 * r.get('registry', np.nan):6.1f}  Er=s*a_s^2 {r['s'] * r['a_s'] ** 2:.2f}")
