#!/usr/bin/env python3
"""Narrow mixed-winding channels: screened vs unscreened (user, 2026-10-16; runs/channel_capacitor). Walls +m (bottom) / -m (top), L 192, passive (and 2 active
seeds at m 6). Per run (second half, every frame): +1/2 and -1/2 within 4 cells of each wall (for W < 8 the bands overlap: a defect is assigned to the nearer wall),
bulk counts, fraction of frames with no defect at all, and the screening fraction S = (bottom +1/2 + top -1/2) / (4m) (1 = each wall fully screened by 2m defects).
Crossover W*: the width where S falls through 1/2 at each m; test W* ~ 1/m (m W* constant). -> capacitor.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); R_ = "/home/baloewe/Physics/Overcharged/runs/channel_capacitor"
def one(run):
    J = json.load(open(run + "/input.json")); BC = J["BC"]; y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; W = y1 - y0; m = abs(BC[2].get("customAnch", 0)); act = J["species"][0]["act"]
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); A = []; empty = 0
    for f in fr:
        if not len(f): empty += 1; A.append((0, 0, 0, 0, 0, 0, 0)); continue
        y = f[:, 1]; db, dt_ = y - y0, y1 - y; bot = (db < 4) & (db <= dt_); top = (dt_ < 4) & (dt_ < db); bulk = ~bot & ~top; q = np.sign(f[:, 2])
        lo = db < W / 2
        A.append(((bot & (q > 0)).sum(), (bot & (q < 0)).sum(), (top & (q > 0)).sum(), (top & (q < 0)).sum(), (bulk & (q > 0)).sum(), (bulk & (q < 0)).sum(), (lo & (q > 0)).sum() + (~lo & (q < 0)).sum()))
    A = np.array(A, float).mean(0); S = (A[0] + A[3]) / (4 * m); S2 = A[6] / (4 * m)
    return dict(run=os.path.basename(run), W=W, m=m, act=act, bp=A[0], bm=A[1], tp=A[2], tm=A[3], blp=A[4], blm=A[5], empty=empty / max(len(fr), 1), S=S, S2=S2)
if __name__ == "__main__":
    runs = sorted(r.rstrip("/") for r in glob.glob(R_ + "/m*/") if os.path.exists(r + "DONE"))
    runs += sorted(r.rstrip("/") for p in ("channel_narrow/W*_mixed*", "channel_mixedW/W*_mixed*", "channel_mixed/mixed0*passive*") for r in glob.glob(f"/home/baloewe/Physics/Overcharged/runs/{p}/") if os.path.exists(r + "DONE"))
    with Pool(6) as p: res = p.map(one, runs)
    L = ["NARROW MIXED CHANNELS (+m bottom / -m top, L 192; second half). S = (bottom +1/2 + top -1/2)/(4m): 1 = fully screened walls, 0 = no wall defects.",
         "S2 = (+1/2 in the bottom half + -1/2 in the top half)/(4m), band-free (screening defects of sparse windings sit 2-6 cells deep, beyond the 4-cell band).",
         "run | m | W | act | bottom +1/2 -1/2 | top +1/2 -1/2 | bulk +1/2 -1/2 | frames with no defect | S | S2 | m*W | W/site spacing (W 2m/L)"]
    for r in sorted(res, key=lambda r: (r["m"], r["act"], r["W"])):
        L.append(f"{r['run']:28s} | {r['m']:2d} | {r['W']:3.0f} | {r['act']:.2f} | {r['bp']:5.2f} {r['bm']:4.2f} | {r['tp']:4.2f} {r['tm']:5.2f} | {r['blp']:4.2f} {r['blm']:4.2f} | {r['empty']:.2f} | {r['S']:.2f} | {r['S2']:.2f} | {r['m']*r['W']:.0f} | {r['W']*2*r['m']/192:.2f}")
    open(os.path.join(H, "capacitor.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
