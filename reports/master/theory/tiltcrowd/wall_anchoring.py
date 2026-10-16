#!/usr/bin/env python3
"""Effective anchoring in MPCD (user, 2026-10-16): how far does the cell director next to the wall deviate from the anchored direction theta_w = k phi,
between and at the wall defects, against the layer spacing? (The continuum clamps it: deviation 0, S = S0.) Main circle sweep at R = 31.83, k 3-20,
director fields every 100 t.u., second half. Cells binned by depth R - r (0-1, 1-2, 2-3, 3-4 cells; cell centres) and by the arc distance to the nearest
wall +1/2 (defect frame nearest in time): 'core' < 1.5 cells, 'between' > a/4 (a = spacing). Per bin: median |delta| (deg), <S>, and the fidelity
<S cos 2 delta>/<S>. -> wall_anchoring.txt, wall_anchoring.npz"""
import os, sys, glob, gzip, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa
import bulk_spiral as BS  # noqa  (frames reader)
DB = [(0, 1), (1, 2), (2, 3), (3, 4)]


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); a = 2 * np.pi * R / (2 * k - 2)
    D = BS.frames(run + "/directorfield.dat.gz", 8); td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]
    t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
    acc = {(b, z): [[], [], []] for b in range(4) for z in ("core", "between")}
    for tt in td:
        d = D[D[:, 0] == tt]; x = d[:, 1] + 0.5 - Q[0]; y = d[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ph = np.arctan2(y, x); dep = R - r
        m = (dep >= 0) & (dep < 4) & (d[:, 7] > 0)
        if not m.any(): continue
        f = fr[int(np.argmin(np.abs(t - tt)))]
        if not len(f): continue
        rw = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rw < 4)
        if w.sum() < 2: continue
        phw = np.arctan2(f[w, 1], f[w, 0]); dphi = np.abs(np.angle(np.exp(1j * (ph[m][:, None] - phw[None, :])))).min(1) * R
        delta = np.angle(np.exp(2j * (np.arctan2(d[m, 5], d[m, 4]) - k * ph[m]))) / 2; S = d[m, 7]
        for b, (lo, hi) in enumerate(DB):
            sel = (dep[m] >= lo) & (dep[m] < hi)
            for z, zs in (("core", dphi < 1.5), ("between", dphi > a / 4)):
                q = sel & zs
                if q.any(): acc[(b, z)][0] += list(np.abs(delta[q])); acc[(b, z)][1] += list(S[q]); acc[(b, z)][2] += list(S[q] * np.cos(2 * delta[q]))
    out = {}
    for key, (dl, S, SC) in acc.items():
        out[key] = (np.degrees(np.median(dl)) if dl else np.nan, np.mean(S) if S else np.nan, np.sum(SC) / np.sum(S) if S else np.nan, len(dl))
    return dict(k=k, a=a, out=out)


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz")
                  if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.05)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs)
    ks = sorted(set(x["k"] for x in res))
    L = ["MPCD WALL ANCHORING vs spacing (R 31.83, second half): median |deviation from k phi| (deg) / <S> / fidelity <S cos 2 delta>/<S>, by depth bin, BETWEEN defects (> a/4) and at CORES (< 1.5 cells)",
         "k | spacing | between: depth 0-1, 1-2, 2-3, 3-4 | cores: depth 0-1, 1-2"]
    for kk in ks:
        S_ = [x for x in res if x["k"] == kk]
        def agg(b, z, i): return np.nanmean([x["out"][(b, z)][i] for x in S_])
        L.append(f"{kk:2d} | {S_[0]['a']:5.1f} | " + "  ".join(f"{agg(b,'between',0):5.1f}/{agg(b,'between',1):.2f}/{agg(b,'between',2):.2f}" for b in range(4)) + " | " +
                 "  ".join(f"{agg(b,'core',0):5.1f}/{agg(b,'core',1):.2f}/{agg(b,'core',2):.2f}" for b in range(2)))
    np.savez(os.path.join(H, "wall_anchoring.npz"), res=np.array(res, dtype=object), allow_pickle=True)
    print("\n".join(L)); open(os.path.join(H, "wall_anchoring.txt"), "w").write("\n".join(L) + "\n")
