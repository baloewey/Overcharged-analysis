#!/usr/bin/env python3
"""Tracking-free mid-line current (2026-10-22). The mid-line of a periodic channel is a closed loop, so the net topological charge above it, Q_up(t)
(all +-1/2 with y > y_mid, the top wall layer included), changes only when charge crosses it; pair creation/annihilation is neutral (a pair born
straddling the line counts as a crossing, correctly). J_true = slope of Q_up(t) per unit length per 1000 t.u. (second half, least squares; s.e. from
10 blocks). Compared with the tracked current of current_stats.py (current_decomp.one 'cross'). Same runs as current_stats.py. -> qup_current.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, "..")
from edge_criteria import read_with_angle  # noqa
from current_stats import classify  # noqa
from current_decomp import one as tracked  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs"
def qup(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); t = np.array(t); q = np.array([0.5 * np.sign(f[f[:, 1] > ym, 2]).sum() if len(f) else np.nan for f in fr])
    ok = np.isfinite(q); t, q = t[ok], q[ok]; sl = np.polyfit(t, q, 1)[0]
    blocks = [np.polyfit(tb, qb, 1)[0] for tb, qb in zip(np.array_split(t, 10), np.array_split(q, 10))]
    x = tracked(run)
    return run, 1000 * sl / L, 1000 * np.std(blocks) / np.sqrt(10) / L, x["cross"], q.max() - q.min()
if __name__ == "__main__":
    cands = sorted(set(glob.glob(f"{R_}/channel_current/*/") + glob.glob(f"{R_}/channel_fix/W*/") + glob.glob(f"{R_}/channel_curve/W*/") + glob.glob(f"{R_}/channel_mixedW/*/") +
                       glob.glob(f"{R_}/channel_mixed/mixed0*/") + glob.glob(f"{R_}/channel_narrow/W*_m6_*/") + glob.glob(f"{R_}/channel/W*_m6_L192_*/")))
    runs = [r.rstrip("/") for r in cands if os.path.exists(r + "DONE")]; cls = {r: classify(r) for r in runs}; runs = [r for r in runs if cls[r] and cls[r][1] >= 14]
    with Pool(8) as p: res = p.map(qup, runs)
    G = {}
    for r, j, e, c, rg in res: G.setdefault(cls[r], []).append((j, c, rg))
    Lt = ["MID-LINE CURRENT: tracking-free (slope of Q_up) vs tracked crossings (current_stats.py), per unit length per 1000 t.u., mean +- s.e.m. over seeds; range of Q_up(t)",
          "W | type | n | J from Q_up | J tracked | Q_up range (max - min) over the second half"]
    for (typ, W) in sorted(G, key=lambda k: (k[1], k[0])):
        a = np.array(G[(typ, W)]); se = lambda v: v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
        Lt.append(f"{W:3d} | {typ:5s} | {len(a)} | {a[:, 0].mean():+.4f} +- {se(a[:, 0]):.4f} | {a[:, 1].mean():+.4f} +- {se(a[:, 1]):.4f} | {a[:, 2].mean():.1f}")
    open("qup_current.txt", "w").write("\n".join(Lt) + "\n"); print("\n".join(Lt))
