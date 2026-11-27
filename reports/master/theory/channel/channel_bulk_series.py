#!/usr/bin/env python3
"""Bulk charge of the symmetric winding channels across densities (user, 2026-11-26; for the Fig. 4(c) candidate). Defects only (no flow), same estimator as
channel_an.py: bulk = net charge of defects >= 4 cells from both walls, t in the second half of each run (channel_an used t 2000-4000 of 4000).
Runs: channel_qscan (m 3/9/12, W 40), channel_qcross (m 6/7/8/10/11, W 40), channel (W40_m6). Per wall winding m (customAnch), so 2m wall +1/2 per wall;
rho = 2m / L (+1/2 per unit wall length, as rho = 2k/P in cavities); Q_b/P with P = 2L. -> channel_bulk_series.txt, channel_bulk_series.npy"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
R = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
runs = sorted(glob.glob(R + "/channel_qscan/W*_m*_*/") + glob.glob(R + "/channel_qcross/W40_m*_s*/") + glob.glob(R + "/channel/W40_m6_L192_s*/"))
rows = []
for run in runs:
    run = run.rstrip("/")
    if not os.path.exists(run + "/DONE"): continue
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; m = J["BC"][2].get("customAnch", 0)
    if m <= 0 or J["BC"][3].get("customAnch", 0) != -m: continue                    # symmetric +m/+m only
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); sel = t >= t[0] + 0.5 * (t[-1] - t[0])
    qb, nw = [], []
    for f, s in zip(fr, sel):
        if not s or not len(f): continue
        near = (f[:, 1] - y0 < 4) | (y1 - f[:, 1] < 4); qb.append(0.5 * np.sign(f[~near, 2]).sum()); nw.append(((f[:, 2] > 0) & near).sum())
    rows.append((m, L, y1 - y0, 2 * m / L, np.mean(qb), np.mean(qb) / (2 * L), np.mean(nw) / (4 * m), os.path.relpath(run, R), t[-1]))
L_ = ["run | m per wall | L | W | rho = 2m/L | bulk charge | Q_b/P (P = 2L) | wall occupancy N_wall/(4m) | t_end"]
for r in rows: L_.append(f"{r[7]:32s} | {r[0]:2d} | {r[1]} | {r[2]} | {r[3]:.4f} | {r[4]:+.2f} | {r[5]:.4f} | {r[6]:.2f} | {r[8]:.0f}")
ms = sorted({r[0] for r in rows}); L_.append("\nper m (mean over runs): m | rho | Q_b/P | s.e. | runs")
out = []
for m in ms:
    g = [r for r in rows if r[0] == m]; v = np.array([r[5] for r in g]); se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
    out.append((m, g[0][3], v.mean(), se, len(v))); L_.append(f"{m:2d} | {g[0][3]:.4f} | {v.mean():.4f} | {se:.4f} | {len(v)}")
np.save(os.path.join(H, "channel_bulk_series.npy"), np.array(out)); open(os.path.join(H, "channel_bulk_series.txt"), "w").write("\n".join(L_) + "\n"); print("\n".join(L_))
