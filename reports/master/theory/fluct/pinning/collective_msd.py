#!/usr/bin/env python3
"""Pinned or spontaneously broken? Mean-square displacement of the collective (mean) arclength coordinate X(t) of the passive wall lattice.
Free (zero mode, spontaneous symmetry breaking with a frozen-in registry): MSD = 2 D_c t with D_c = mu kBT / N (rigid sliding of N defects,
mu = 0.02-0.03 cells^2/(kBT t.u.)); pinned: MSD saturates at 2 var(X). Same tracking as measure_passive.py. Writes collective_msd.txt/png."""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); FL = os.path.dirname(H)
src = open(os.path.join(FL, "measure_passive.py")).read().split("res = {}")[0]; ns = {}; exec(compile(src, "mp", "exec"), ns)
L = []; fig, ax = plt.subplots(figsize=(5.5, 4))
for lab in ["circle_R31.83_k04", "circle_R31.83_k06", "circle_R31.83_k08", "circle_R31.83_k10", "circle_R27.28_k06",
            "ellipse_eps1.5_k06", "ellipse_eps2.0_k06"]:
    run = os.path.join(ns["ROOT"], ns["RUNS"][lab]); bc = json.load(open(run + "/input.json"))["BC"][0]; N = 2 * bc["customAnch"]
    band = 6.0 if lab.endswith(("k04", "k05")) else 4.0
    t, fr = ns["read_with_angle"](run, np.array(bc["Q"][:2]), 0.0); W, s, P = ns["wall_polyline"](bc)
    T, S = [], []
    for tt, f in zip(t, fr):
        p = f[f[:, 2] > 0]; d, sa = ns["project"](W, s, p[:, :2]); m = d < band
        if m.sum() == N and len(p) == N: T.append(tt); S.append(np.sort(sa[m]))
    S = np.array(S); T = np.array(T)
    for n in range(1, len(S)):
        best = min(range(N), key=lambda sh: (((np.roll(S[n], -sh) - S[n - 1] + P / 2) % P - P / 2) ** 2).sum()); S[n] = np.roll(S[n], -best)
    X = (np.unwrap(S * 2 * np.pi / P, axis=0) * P / (2 * np.pi)).mean(1)
    dt = np.median(np.diff(T)); idx = np.round((T - T[0]) / dt).astype(int); Xf = np.full(idx[-1] + 1, np.nan); Xf[idx] = X
    lags = np.unique(np.round(np.logspace(0, np.log10(len(Xf) / 3), 25)).astype(int))
    msd = np.array([np.nanmean((Xf[l:] - Xf[:-l]) ** 2) for l in lags])
    free = 2 * (0.025 / N) * lags * dt
    late = msd[-5:].mean(); slope = np.polyfit(np.log(lags[-8:] * dt), np.log(msd[-8:]), 1)[0]
    L.append(f"{lab:20s} N={N:2d}: {len(S)} frames over {T[-1] - T[0]:.0f} t.u.; collective drift first->last {X[-1] - X[0]:+.2f} cells; "
             f"MSD at lag {lags[-1] * dt:.0f}: {late:.2f} cells^2 (free sliding would give {free[-1]:.1f}); late log-slope {slope:.2f} (1 = diffusive, 0 = pinned)")
    print(L[-1]); ax.loglog(lags * dt, msd, "o-", ms=3, label=lab)
ax.loglog(lags * dt, 2 * (0.025 / 12) * lags * dt, "k--", label="free sliding, N = 12")
ax.set_xlabel("lag (time units)"); ax.set_ylabel(r"MSD of the collective coordinate (cells$^2$)"); ax.legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(H, "collective_msd.png"), dpi=150)
open(os.path.join(H, "collective_msd.txt"), "w").write("\n".join(L) + "\n")
