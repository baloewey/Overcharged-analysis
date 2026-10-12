#!/usr/bin/env python3
"""Active circles (N = 2k-2 wall +1/2): one-constant theory says a uniform head tilt tau requires a rigid lattice rotation alpha = tau/(2k-2)
away from the zero-tilt registry, i.e. D = arg Z - tau = const for all frames and runs, Z = (1/N_w) sum_j exp(i (2k-2) phi_j),
tau = circular mean over wall +1/2 of (head angle - phi_j). Compare the concentration (Rbar over all frames of all runs) of
arg Z alone (registry), tau alone (tilt) and D (registry corrected for the tilt), and the per-frame correlation of arg Z with tau.
Frames with exactly 2k-2 wall +1/2 only. Writes tilt_registry_active.txt."""
import glob, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); out = []
def rb(x): return abs(np.mean(np.exp(1j * x)))
for k in range(3, 11):
    A, T, D, runs_used, cors = [], [], [], 0, []
    for f in sorted(glob.glob(f"{ROOT}/cavity/server_main/main_act0.07/rho*/k{k:02d}/s*/defects.dat.gz")):
        run = os.path.dirname(f); bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]
        t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.2)
        a, ta = [], []
        for g in fr[::5]:
            p = g[g[:, 2] > 0]; ph = np.arctan2(p[:, 1], p[:, 0]); m = R - np.hypot(p[:, 0], p[:, 1]) < 4
            if m.sum() != 2 * k - 2: continue
            z = np.exp(1j * (2 * k - 2) * ph[m]).mean()
            if abs(z) < 0.5: continue
            a.append(np.angle(z)); ta.append(np.angle(np.mean(np.exp(1j * (p[m, 3] - ph[m])))))
        if len(a) < 20: continue
        a, ta = np.array(a), np.array(ta); runs_used += 1
        A += list(a); T += list(ta); D += list(np.angle(np.exp(1j * (a - ta))))
        if np.std(np.angle(np.exp(1j * (ta - np.angle(np.mean(np.exp(1j * ta))))))) > 0.1:
            ac = np.angle(np.exp(1j * (a - np.angle(np.mean(np.exp(1j * a)))))); tc = np.angle(np.exp(1j * (ta - np.angle(np.mean(np.exp(1j * ta))))))
            cors.append((np.corrcoef(tc, ac)[0, 1], np.polyfit(tc, ac, 1)[0]))
    A, T, D = map(np.array, (A, T, D))
    line = (f"k={k:2d}: {runs_used} runs, {len(A)} frames; concentration Rbar: registry arg Z {rb(A):.2f}, tilt tau {rb(T):.2f}, "
            f"corrected D = arg Z - tau {rb(D):.2f} (mean {np.angle(np.mean(np.exp(1j * D))):+.2f}); within-run corr(arg Z, tau) median "
            f"{np.median([c[0] for c in cors]) if cors else np.nan:+.2f}, slope median {np.median([c[1] for c in cors]) if cors else np.nan:+.2f} (theory 1)")
    print(line); out.append(line)
open(os.path.join(H, "tilt_registry_active.txt"), "w").write("\n".join(out) + "\n")
