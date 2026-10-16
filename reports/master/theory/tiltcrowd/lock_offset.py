#!/usr/bin/env python3
"""Is the wall-defect orientation relaxed on its position (the lock), or twisted out of it in crowded layers? (user, 2026-10-15: the log repulsion assumes
orientations relaxed to their minimum; unrelaxed defects exert torques; maybe likewise with the wall.)
Lock: tau_j = (2k-2) phi_j + c with c fixed by the anchoring (circle). If the flow twists each core by delta beyond its relaxed orientation, c shifts by
delta, signed with the tilt. Per coherent frame (C > 0.9) of the 10-seed circle window: c_frame = arg<exp(i (psi_j - (2k-1) phi_j))>, signed by the sign of
the frame-mean tilt; also the spread of c_j within the frame. -> lock_offset.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); cs, sp, ts = [], [], []
    for f in fr:
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        if w.sum() < 3: continue
        phi = np.arctan2(f[w, 1], f[w, 0]); tau = f[w, 3] - phi; z = np.mean(np.exp(1j * tau))
        if abs(z) < 0.9: continue
        s = np.sign(np.angle(z)); zc = np.mean(np.exp(1j * (f[w, 3] - (2 * k - 1) * phi)))
        cs.append(s * np.angle(zc)); sp.append(abs(zc)); ts.append(abs(np.angle(z)))
    return dict(k=k, R=R, a=2 * np.pi * R / (2 * k - 2), c=np.degrees(np.median(cs)) if len(cs) > 20 else np.nan, res=np.median(sp) if sp else np.nan,
                tilt=np.degrees(np.median(ts)) if len(ts) > 20 else np.nan)


if __name__ == "__main__":
    runs = []
    for d in sorted(glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/DONE")):
        d = os.path.dirname(d); bc = json.load(open(d + "/input.json"))["BC"][0]
        if 2 <= bc["customAnch"] <= 10 and 0.8 <= bc["R"] / LA <= 2.2: runs.append(d)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs, chunksize=4)
    a = np.array([x["a"] for x in res]); c = np.array([x["c"] for x in res]); T = np.array([x["tilt"] for x in res]); rs = np.array([x["res"] for x in res])
    L = ["LOCK OFFSET (circle window, 10 seeds, coherent frames C > 0.9): c = arg<exp(i(psi - (2k-1)phi))> signed by the tilt sense (deg); resultant = spread of c within a frame",
         "spacing bin | runs | tilt (median) | c signed (median) | resultant of c"]
    for lo, hi in ((5, 8), (8, 10), (10, 12), (12, 15), (15, 20), (20, 30), (30, 60)):
        s = (a >= lo) & (a < hi) & np.isfinite(c)
        if s.any(): L.append(f"[{lo},{hi}) | {s.sum():3d} | {np.median(T[s]):5.1f} | {np.median(c[s]):+5.1f} | {np.median(rs[s]):.3f}")
    ok = np.isfinite(c) & np.isfinite(T); L.append(f"corr(c, tilt) {np.corrcoef(c[ok], T[ok])[0,1]:+.2f}; slope dc/dtilt {np.polyfit(T[ok], c[ok], 1)[0]:+.2f}")
    print("\n".join(L)); open(os.path.join(H, "lock_offset.txt"), "w").write("\n".join(L) + "\n")
