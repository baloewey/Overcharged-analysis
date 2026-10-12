#!/usr/bin/env python3
"""Onset of the coherent wall tilt vs cavity size, per k (test of the effective-boundary instability, which predicts the onset to move to
smaller R/l_a as k grows). Main sweep, active 0.07, all densities (rho up to 0.50 -> small cavities), seeds s01-s03, defects only.
Per frame: wall +1/2 (within 6 cells of the wall, as flow_tilt.py), Delta = psi - phi (psi = recorded angle = tail direction);
coherence = |circular mean of exp(i Delta)| x sign-free mean angle: tilt_coh = time average of |circ. mean angle of Delta| (rad).
Onset: R/l_a where tilt_coh (seed mean) first exceeds 30 deg. Writes tilt_onset_k.txt."""
import glob, json, os, sys
import numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.2); vals = []
    for f in fr[::5]:
        p = f[f[:, 2] > 0]; r = np.hypot(p[:, 0], p[:, 1]); m = R - r < 6
        if m.sum() < 2: continue
        dl = p[m, 3] - np.arctan2(p[m, 1], p[m, 0]); vals.append(abs(np.angle(np.mean(np.exp(1j * dl)))))
    return run, R, np.mean(vals) if vals else np.nan
if __name__ == "__main__":
    runs = [os.path.dirname(f) for f in sorted(glob.glob(ROOT + "/cavity/server_main/main_act0.07/rho*/k*/s0*/DONE"))]
    runs = [r for r in runs if int(r.split("/")[-2][1:]) in (2, 3, 4, 5, 6, 8, 10)]
    with Pool(8) as pool: res = pool.map(one, runs)
    from collections import defaultdict
    D = defaultdict(lambda: defaultdict(list))
    for run, R, tc in res: D[int(run.split("/")[-2][1:])][round(R / LA, 3)].append(tc)
    L = []
    for k in sorted(D):
        xs = np.array(sorted(D[k])); y = np.degrees(np.array([np.nanmean(D[k][x]) for x in xs]))
        on = np.nan
        for i in range(1, len(xs)):
            if y[i - 1] < 30 <= y[i]: on = xs[i - 1] + (30 - y[i - 1]) * (xs[i] - xs[i - 1]) / (y[i] - y[i - 1]); break
        L.append(f"k={k:2d}: onset R/l_a {on:.2f}; tilt_coh (deg) by R/l_a: " + " ".join(f"{x:.2f}:{v:.0f}" for x, v in zip(xs, y)))
        print(L[-1])
    open(os.path.join(H, "tilt_onset_k.txt"), "w").write("\n".join(L) + "\n")
