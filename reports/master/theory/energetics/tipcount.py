#!/usr/bin/env python3
"""Where along the wall do depleting escapes happen? (user, 2026-10-11: tips vs flat parts). Site = escaper position at departure; ellipse angle parameter
t = atan2(y/B, x/A); TIP if |cos t| > 0.866 (within 30 deg of the long axis in parameter), FLANK if |cos t| < 0.5, MID otherwise. Normalised by the time-averaged number of resident
wall +1/2 in each region (second half) -> escapes per wall defect per 1000 t.u. Usage: tipcount.py TAG RUNDIR [...] -> tipcount_<TAG>.npy (per-event rows)"""
import os, sys, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
def reg(p, A, B):
    c = abs(np.cos(np.arctan2(p[1] / B, p[0] / A))); return 0 if c > 0.866 else (2 if c < 0.5 else 1)
def one(run):
    try:
        J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); ainv = bc.get("aInv", [1, 1, 1]); B = float(bc["R"]); A = B * float(ainv[0] or 1)
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); sg = 1 if k > 0 else -1; occ = np.zeros(3); nf = 0
        for f in fr[len(fr) // 2:]:
            if not len(f): continue
            d = tree.query(f[:, :2])[0]; nf += 1
            for p in f[(np.sign(f[:, 2]) == sg) & (d < 4), :2]: occ[reg(p, A, B)] += 1
        T = (t[-1] - t[len(t) // 2]) if len(t) > 2 else 0
        evs = [reg(e["chain"][e["n"]], A, B) for e in ev if t[e["n"]] >= t[len(t) // 2]]
        return (os.path.relpath(run, En.ROOT), A / B, A * B, k, occ / max(nf, 1), T, evs)
    except Exception as ex:
        return (run, None, str(ex))
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(24) as pool: res = pool.map(one, runs, chunksize=2)
    res = [r for r in res if r[1] is not None]; np.save(os.path.join(H, f"tipcount_{tag}.npy"), np.array(res, dtype=object), allow_pickle=True)
    for lo, hi in ((0.99, 1.01), (1.15, 1.25), (1.35, 1.45), (1.55, 1.65), (1.95, 2.05), (2.45, 2.55)):
        S = [r for r in res if lo <= r[1] <= hi]
        if not S: continue
        n = np.zeros(3); occT = np.zeros(3)
        for r in S:
            for g in r[6]: n[g] += 1
            occT += r[4] * r[5]
        rate = 1000 * n / np.maximum(occT, 1e-9)
        print(f"aspect {lo + 0.01:.1f}: {len(S)} runs; escapes tip/mid/flank {n.astype(int)}; resident occupancy {np.round(sum(r[4] for r in S) / len(S), 2)}; rate per wall defect per 1000 t.u. {np.round(rate, 3)}")
