#!/usr/bin/env python3
"""Do winding channels pump fluid along x? (user, 2026-10-16) Per run (periodic channel, second half): channel-averaged <v_x>(t) over the fluid cells of each
flow frame; time mean, rms, fraction of time > 0, sign changes per 1000 t.u., and the time-averaged profile <v_x>(y) (shown at the walls, quarter and middle).
Compared with the rms flow speed. Valid (correct-density) sets only. -> net_flow.txt"""
import os, sys, glob, json, gzip, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); R_ = "/home/baloewe/Physics/Overcharged/runs"
def opn(fn): return gzip.open(fn + ".gz", "rt") if os.path.exists(fn + ".gz") else open(fn)
def one(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; BC = J["BC"]; y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; W = y1 - y0
    if abs(J["species"][0]["pop"] / (L * (W + 2)) - 20) > 1: return None
    rows = np.loadtxt((l for l in opn(run + "/flowfield.dat") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    ts = np.unique(rows[:, 0]); ts = ts[ts >= ts[0] + 0.5 * (ts[-1] - ts[0])]; m = (rows[:, 2] >= y0) & (rows[:, 2] < y1)
    sel = np.isin(rows[:, 0], ts) & m; r = rows[sel]
    tt, inv = np.unique(r[:, 0], return_inverse=True); vx = np.bincount(inv, r[:, 4]) / np.bincount(inv); v2 = np.bincount(inv, r[:, 4] ** 2 + r[:, 5] ** 2) / np.bincount(inv)
    yb = (r[:, 2] - y0).astype(int); prof = np.bincount(yb, r[:, 4]) / np.maximum(np.bincount(yb), 1)
    dt = np.median(np.diff(tt)); sc = np.sum(np.diff(np.sign(vx)) != 0) / (len(vx) * dt) * 1000
    cb, ct = BC[2].get("customAnch", 0), -BC[3].get("customAnch", 0); act = J["species"][0]["act"]
    return dict(run=os.path.relpath(run, R_), W=W, cb=cb, ct=ct, act=act, mean=vx.mean(), rms=vx.std(), fpos=np.mean(vx > 0), sc=sc, vrms=np.sqrt(v2.mean()), prof=prof)
if __name__ == "__main__":
    runs = sorted(set(r.rstrip("/") for p in ("channel/W*_L192_*/", "channel_mixedW/W*/", "channel_mixed/*/", "channel_long/*/", "channel_current/W24_*/", "channel_narrow/W*/", "channel_qscan/*/", "channel_fix/W*/", "channel_curve/W*/")
                      for r in glob.glob(f"{R_}/{p}") if os.path.exists(r + "DONE")))
    with Pool(6) as pool: res = [x for x in pool.map(one, runs) if x]
    L = ["NET FLOW ALONG PERIODIC CHANNELS (second half; <v_x> = channel average per flow frame; v_rms = rms speed). Valid-density runs.",
         "run | W | wall charges | act | <v_x> mean | rms of <v_x>(t) | frac > 0 | sign changes /1000 t.u. | v_rms | <v_x>(y): wall+2, W/4, W/2, 3W/4, wall-2"]
    for x in res:
        p = x["prof"]; W = int(x["W"]); idx = [2, W // 4, W // 2, 3 * W // 4, W - 3]
        L.append(f"{x['run']:44s} | {W:3d} | {x['cb']:+d}/{x['ct']:+d} | {x['act']:.2f} | {x['mean']:+.4f} | {x['rms']:.4f} | {x['fpos']:.2f} | {x['sc']:6.1f} | {x['vrms']:.4f} | " + " ".join(f"{p[i]:+.3f}" for i in idx if i < len(p)))
    print("\n".join(L)); open(os.path.join(H, os.environ.get("OUT", "net_flow.txt")), "w").write("\n".join(L) + "\n")
