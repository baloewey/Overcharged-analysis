#!/usr/bin/env python3
"""One-constant theory: the director is fixed by the defect positions, so a rigid rotation of the wall lattice by alpha tilts every +1/2 head
by (2k-2) alpha relative to its wall normal (theta_rest at defect j shifts by sum_{i!=j} q_i alpha = (k - 1/2) alpha; head = 2 theta_rest;
normal shifts by alpha). Test on passive circles: per frame, collective angular displacement alpha (mean tracked arclength shift / R) vs
mean head tilt (head angle minus outward normal, circular mean, minus its time mean). Writes tilt_slaving.txt."""
import os, sys, json
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); FL = os.path.dirname(H)
src = open(os.path.join(FL, "measure_passive.py")).read().split("res = {}")[0]; ns = {}; exec(compile(src, "mp", "exec"), ns)
out = []
for k in range(2, 11):
    lab = f"circle_R31.83_k{k:02d}"; run = os.path.join(ns["ROOT"], ns["RUNS"][lab])
    bc = json.load(open(run + "/input.json"))["BC"][0]; N = 2 * k; R = bc["R"]; band = 9.0 if k <= 3 else 6.0 if k <= 5 else 4.0
    t, fr = ns["read_with_angle"](run, np.array(bc["Q"][:2]), 0.2)
    S, TL = [], []
    for f in fr:
        p = f[f[:, 2] > 0]; r = np.hypot(p[:, 0], p[:, 1]); m = R - r < band
        if m.sum() != N or len(p) != N: continue
        ph = np.arctan2(p[m, 1], p[m, 0]); o = np.argsort(ph % (2 * np.pi))
        S.append((ph % (2 * np.pi))[o]); TL.append((p[m, 3][o] - ph[o]))
    S = np.array(S); TL = np.array(TL)
    for n in range(1, len(S)):
        best = min(range(N), key=lambda sh: (((np.roll(S[n], -sh) - S[n - 1] + np.pi) % (2 * np.pi) - np.pi) ** 2).sum())
        S[n] = np.roll(S[n], -best); TL[n] = np.roll(TL[n], -best)
    U = np.unwrap(S, axis=0); al = (U - U.mean(0)).mean(1)                  # collective rotation (rad)
    TL = np.array(TL); mu = np.angle(np.mean(np.exp(1j * TL), axis=0))       # per-defect circular time mean
    dtl = np.angle(np.exp(1j * (TL - mu)))                                     # per-defect tilt fluctuation
    tl = dtl.mean(1)                                                           # uniform (collective) tilt fluctuation
    slope = np.polyfit(al, tl, 1)[0]; c = np.corrcoef(al, tl)[0, 1]
    line = (f"k={k:2d} N={N:2d}: {len(S)} frames, sd(alpha) {al.std():.4f} rad ({al.std() * R:.2f} cells), sd(mean tilt) {tl.std():.3f} rad; "
            f"slope d tilt/d alpha {slope:+.2f} (one-constant prediction {2 * k - 2}), correlation {c:+.2f}; mean tilt about normal "
            f"per defect: " + " ".join(f"{np.degrees(x):+.0f}" for x in mu) + f" deg; sd per-defect tilt {dtl.std():.2f} rad")
    print(line); out.append(line)
open(os.path.join(H, "tilt_slaving.txt"), "w").write("\n".join(out) + "\n")
