#!/usr/bin/env python3
"""Pinning vs spontaneous symmetry breaking (user's question, 2026-10-07). All production runs start from the same director (oDist = 2,
aligned along y), so a shared registry could be inherited from the start. Eight passive k = 6, R = 31.83 runs (runs/pinning_ic) differ only
in the initial orientation (random x4, along x, at 45 deg, along y x2 with new seeds). Lattice phase: arg Z, Z = <(1/N) sum_j exp(i N phi_j)>
over frames with N = 12 wall +1/2 (phase of the 12-site lattice, defined mod 2 pi; site angle offset = arg Z / 12), in the second half of
the run, plus its time course (phase in 300 t.u. windows) and |Z| (lattice order). Reference: the original run (Rscan_R31.83_passive/k06).
Pinned: all runs share arg Z (up to the spread seen within runs). Spontaneous: random starts scatter over [0, 2 pi).
Writes ic_test.txt."""
import glob, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
N = 12
runs = [("reference (alignY, seed 90006)", ROOT + "/cavity/Rscan_R31.83_passive/k06")] + [(os.path.basename(d), d) for d in sorted(glob.glob(ROOT + "/pinning_ic/*"))]
phs = []
for lab, run in runs:
    if not os.path.exists(run + "/DONE"): say(f"{lab}: not finished"); continue
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.0); t = np.array(t)
    zs, ts = [], []
    for tt, g in zip(t, fr):
        p = g[g[:, 2] > 0]; m = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if m.sum() == N and len(p) == N: zs.append(np.exp(1j * N * np.arctan2(p[m, 1], p[m, 0])).mean()); ts.append(tt)
    zs, ts = np.array(zs), np.array(ts)
    if len(zs) < 50: say(f"{lab}: only {len(zs)} complete frames"); continue
    half = ts > ts.min() + (ts.max() - ts.min()) / 2
    Z = zs[half].mean(); phs.append(np.angle(Z))
    win = [np.angle(zs[(ts >= a) & (ts < a + 300)].mean()) for a in np.arange(ts.min(), ts.max(), 300) if ((ts >= a) & (ts < a + 300)).sum() > 10]
    say(f"{lab:32s}: |Z| {abs(Z):.2f}, arg Z {np.angle(Z):+.2f} rad (site offset {np.degrees(np.angle(Z)) / N:+.2f} deg); "
        f"time course (300 t.u. windows): " + " ".join(f"{w:+.2f}" for w in win))
if len(phs) > 2:
    phs = np.array(phs); Rb = abs(np.mean(np.exp(1j * phs)))
    say(f"\nphase agreement over {len(phs)} runs: Rbar {Rb:.2f} (random phases ~{1 / np.sqrt(len(phs)):.2f}); Rayleigh p ~ {np.exp(-len(phs) * Rb ** 2):.1e}; "
        f"circular spread {np.degrees(np.sqrt(-2 * np.log(max(Rb, 1e-9)))) / N:.2f} deg of site angle")
open(os.path.join(H, "ic_test.txt"), "w").write("\n".join(L) + "\n")
