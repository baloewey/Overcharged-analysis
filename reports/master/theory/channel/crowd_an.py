#!/usr/bin/env python3
"""Crowded winding channels (runs/channel_crowd; user 2026-10-15: crowding is the natural next step). L 192; symmetric W 40 (wall charge m 18/24/30/36 per wall)
and mixed W 24 (m 3/12/18/24) and W 40 (m 12/18/24); 2 seeds each. Per run (second half): wall +1/2 and -1/2 within 4 cells of each wall, bulk net charge
(>= 4 cells), bulk charge per unit length, released fraction (sym: Q_b / total imposed charge), wall-layer tilt (coherent frames), and for mixed channels the
topological current (current_decomp.one, crossings + motion/events, per unit length per 1000 t.u.). Compared with the earlier m 6 channels (current_stats.txt).
-> crowd_an.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from current_decomp import one as cur_one  # noqa
R_ = os.environ.get("RDIR", "/home/baloewe/Physics/Overcharged/runs/channel_crowd")


def one(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; BC = J["BC"]; y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; W = y1 - y0
    cb, ct = BC[2].get("customAnch", 0), -BC[3].get("customAnch", 0); mixed = cb * ct < 0
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); acc = []; tl = []
    for f in fr:
        if not len(f): continue
        y = f[:, 1]; bot = y - y0 < 4; top = y1 - y < 4; bulk = ~bot & ~top; q = np.sign(f[:, 2])
        acc.append(((bot & (q > 0)).sum(), (bot & (q < 0)).sum(), (top & (q > 0)).sum(), (top & (q < 0)).sum(), 0.5 * q[bulk].sum()))
        w = bot & (q > 0)
        if w.sum() >= 3:
            z = np.mean(np.exp(1j * (f[w, 3] + np.pi / 2)))           # outward normal of the bottom wall is -y
            if abs(z) > 0.9: tl.append(abs(np.angle(z)))
    A = np.array(acc).mean(0) if acc else np.zeros(5); out = dict(run=os.path.basename(run), W=W, m=abs(cb), cb=cb, ct=ct, mixed=mixed, bp=A[0], bm=A[1], tp=A[2], tm=A[3], qb=A[4], L=L,
                                         tilt=np.degrees(np.median(tl)) if len(tl) > 20 else np.nan)
    if mixed:
        c = cur_one(run); out["cross"] = c["cross"]; out["mot"] = (sum(c["Jm"].values()) + c["ev_c"] + c["ev_a"]) / (W - 8)
    return out


if __name__ == "__main__":
    runs = sorted(r.rstrip("/") for r in glob.glob(R_ + "/" + os.environ.get("PAT", "W*_s*") + "/") if os.path.exists(r + "DONE"))
    with Pool(6) as pool: res = pool.map(one, runs)
    L = ["CROWDED WINDING CHANNELS (L 192, second half; wall = within 4 cells; per run, then mean of 2 seeds)",
         "run | W | wall charges (bottom/top) | bottom +1/2, -1/2 | top +1/2, -1/2 | bulk charge (per length) | released fraction | tilt bottom layer | current crossings / motion+events | charge closure: measured total - imposed"]
    for r in res:
        tot = r["cb"] + r["ct"]; rel = r["qb"] / tot if tot else np.nan
        L.append(f"{r['run']:20s} | {r['W']:.0f} | {r['cb']:+d}/{r['ct']:+d} | {r['bp']:5.1f} {r['bm']:4.1f} | {r['tp']:5.1f} {r['tm']:4.1f} | {r['qb']:+6.2f} ({r['qb']/r['L']:+.4f}) | {rel:5.2f} | {r['tilt']:5.1f} | "
                 + (f"{r['cross']:+.4f} / {r['mot']:+.4f}" if r["mixed"] else "-") + f" | {0.5 * (r['bp'] - r['bm'] + r['tp'] - r['tm']) + r['qb'] - tot:+.2f}")
    print("\n".join(L)); open(os.path.join(H, os.environ.get("OUT", "crowd_an.txt")), "w").write("\n".join(L) + "\n")
