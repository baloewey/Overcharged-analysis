#!/usr/bin/env python3
"""Active registry vs initial condition (k = 6, R = 31.83, act 0.07; runs/pinning_ic_active + the production seeds of the same cavity).
Per run: frames with 2k-2 = 10 wall +1/2 (within 4 cells): Z = (1/10) sum exp(i 10 phi_j), and the lock-corrected registry D = arg Z - tau
(tau = circular mean of tail angle - phi), which removes the tilt-driven shift; time-averaged over the second half. Pinned at the predicted
site: arg <exp(iD)> ~ 0 in every run regardless of the initial director. Writes ic_test_active.txt."""
import glob, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
runs = [(os.path.basename(d), d) for d in sorted(glob.glob(ROOT + "/pinning_ic_active/*"))] + \
       [(f"production {s} (alignY)", ROOT + f"/cavity/Rscan_R31.83_seeds/{s}/k06") for s in ("s01", "s02", "s03")]
DS, AZ = [], []
for lab, run in runs:
    if not os.path.exists(run + "/DONE"): say(f"{lab}: not finished"); continue
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.5)
    d, a, tl = [], [], []
    for f in fr:
        p = f[f[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); m = R - np.hypot(p[:, 0], p[:, 1]) < 4
        if m.sum() != 10: continue
        z = np.exp(1j * 10 * ph[m]).mean(); tau = np.angle(np.mean(np.exp(1j * (p[m, 3] - ph[m]))))
        a.append(z); d.append(np.exp(1j * (np.angle(z) - tau))); tl.append(tau)
    if len(a) < 20: say(f"{lab}: only {len(a)} frames with 10 wall defects"); continue
    Dm = np.mean(d); Am = np.mean(np.exp(1j * np.angle(a))); DS.append(np.angle(Dm)); AZ.append(np.angle(Am))
    say(f"{lab:26s}: {len(a)} frames; registry arg Z {np.angle(Am):+.2f} (concentration {abs(Am):.2f}); lock-corrected D {np.angle(Dm):+.2f} "
        f"(concentration {abs(Dm):.2f}); mean |tilt| {np.degrees(np.mean(np.abs(tl))):.0f} deg")
if len(DS) > 2:
    for nm, v in (("arg Z", AZ), ("D", DS)):
        v = np.array(v); Rb = abs(np.mean(np.exp(1j * v)))
        say(f"across {len(v)} runs: {nm} Rbar {Rb:.2f} (random ~{1 / np.sqrt(len(v)):.2f}), Rayleigh p ~ {np.exp(-len(v) * Rb ** 2):.1e}, mean {np.angle(np.mean(np.exp(1j * v))):+.2f}")
open(os.path.join(H, "ic_test_active.txt"), "w").write("\n".join(L) + "\n")
