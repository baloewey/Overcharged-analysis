#!/usr/bin/env python3
"""Defect mobility from passive wall-lattice dynamics (model-free): for each Fourier mode m of the tangential displacements (circles),
overdamped harmonic relaxation gives Gamma_m = mu k_m with k_m = kBT/<|u_m|^2>, so mu = Gamma_m <|u_m|^2> / kBT (kBT = 1).
Gamma_m from the 1/e time of the mode autocorrelation. Uses the tracking of measure_passive.py. Writes mobility.txt."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("mp", os.path.join(HERE, "measure_passive.py"))
src = open(os.path.join(HERE, "measure_passive.py")).read().split("res = {}")[0]      # functions only
ns = {}; exec(compile(src, "mp", "exec"), ns)
lines = []
for lab in ["circle_R31.83_k%02d" % k for k in range(4, 11)] + ["circle_R27.28_k06"]:
    run = os.path.join(ns["ROOT"], ns["RUNS"][lab])
    j = json.load(open(os.path.join(run, "input.json"))); bc = j["BC"][0]; N = 2 * bc["customAnch"]
    t, fr = ns["read_with_angle"](run, np.array(bc["Q"][:2]), 0.2)
    W, s, P = ns["wall_polyline"](bc)
    band = 6.0 if lab.endswith(("k04", "k05")) else 4.0
    S = []; T = []
    for tt, f in zip(t, fr):
        p = f[f[:, 2] > 0]; d, sa = ns["project"](W, s, p[:, :2]); m = d < band
        if m.sum() != N or len(p) != N: continue
        S.append(np.sort(sa[m])); T.append(tt)
    S = np.array(S); T = np.array(T); dt = np.median(np.diff(t))
    idx = np.round((T - T[0]) / dt).astype(int)                       # frame index (gaps where the count was wrong)
    for n in range(1, len(S)):
        best = min(range(N), key=lambda sh: (((np.roll(S[n], -sh) - S[n - 1] + P / 2) % P - P / 2) ** 2).sum())
        S[n] = np.roll(S[n], -best)
    U = np.unwrap(S * 2 * np.pi / P, axis=0) * P / (2 * np.pi); dU = U - U.mean(0)
    um0 = np.fft.fft(dU, axis=1) / np.sqrt(N)
    um = np.full((idx[-1] + 1, N), np.nan + 0j); um[idx] = um0
    out = []
    for m in range(1, N // 2 + 1):
        x = um[:, m]; var = np.nanmean(np.abs(x) ** 2)
        ac = np.array([np.nanmean(np.real(x[: len(x) - L] * np.conj(x[L:]))) / var for L in range(0, 200)])
        k = np.argmax(ac < np.exp(-1)) if (ac < np.exp(-1)).any() else None
        if not k: continue
        tau = dt * (k - 1 + (ac[k - 1] - np.exp(-1)) / (ac[k - 1] - ac[k]))
        out.append((m, var, tau, var / tau))
    mus = [o[3] for o in out]
    lines.append(f"{lab}: {len(S)} frames (dt {dt:g}); per mode m: <|u_m|^2>, tau_m, mu = " +
                 "; ".join(f"m{m}: {v:.2f}, {ta:.1f}, {mu:.3f}" for m, v, ta, mu in out) + f"  -> median mu {np.median(mus):.3f}")
    print(lines[-1])
open(os.path.join(HERE, "mobility.txt"), "w").write("\n".join(lines) + "\n")
