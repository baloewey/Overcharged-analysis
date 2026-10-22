#!/usr/bin/env python3
"""Diagnosis of the mid-line crossing count (user, 2026-10-22: 'fast-moving defects are playing tricks on the accounting').
(1) frame-to-frame displacement of each defect to its nearest same-sign defect in the next frame (periodic in x): distribution, bulk vs wall band.
(2) tracking-free check: Q_up(t) = net charge above the mid-line (all +-1/2, walls included). Pair creation/annihilation is neutral, so Q_up(t) - Q_up(0)
    is exactly the net charge carried up across the mid-line (plus detection noise). Compared with the old crossing count (crossing_movie.track)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, "..")
from edge_criteria import read_with_angle  # noqa
import crossing_movie as CM  # noqa
run = sys.argv[1]; J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); dt = t[1] - t[0]; D = {"bulk": [], "wall": []}
for a, b in zip(fr[:-1], fr[1:]):
    for sg in (1, -1):
        p = a[np.sign(a[:, 2]) == sg]; q = b[np.sign(b[:, 2]) == sg]
        if not len(p) or not len(q): continue
        for x, y in p[:, :2]:
            dx = (q[:, 0] - x + L / 2) % L - L / 2; d = np.hypot(dx, q[:, 1] - y).min(); D["wall" if min(y - y0, y1 - y) < 4 else "bulk"].append(d)
for k, v in D.items():
    v = np.array(v); print(f"{k}: n {len(v)}; nearest-match displacement per frame ({dt:.0f} t.u.) quantiles 50/90/99: {np.quantile(v, [.5, .9, .99]).round(2)}; fraction > 2.5 cells {np.mean(v > 2.5):.3f}")
qu = np.array([0.5 * np.sign(f[f[:, 1] > ym, 2]).sum() for f in fr]); ntot = np.array([0.5 * np.sign(f[:, 2]).sum() for f in fr])
_, _, _, _, cum, ev = CM.track(run)
print(f"total charge per frame: mean {ntot.mean():.2f}, frames off the mode {np.mean(ntot != np.median(ntot)):.3f}")
for T in (1000, 2000, 4000, t[-1]):
    i = int(np.argmin(np.abs(np.array(t) - T))); j = slice(max(i - 25, 0), i + 25)
    print(f"t {t[i]:6.0f}: Q_up - Q_up(0) = {np.mean(qu[j]) - np.mean(qu[:50]):+6.1f} (50-frame means) | old crossing count {cum[i, 0]:+6.1f}")
