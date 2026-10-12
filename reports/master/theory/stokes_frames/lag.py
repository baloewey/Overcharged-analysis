#!/usr/bin/env python3
"""Who leads, wall tilt or circulation? Fixed-radius scan (R = 31.83), seeds 2-3, k = 2..10. Per flow interval (50 time units):
circulation C = spatial mean of v_theta (r < R-1); wall tilt T = mean over wall +1/2 (s < 4) and defect frames in the interval of
sin(psi - phi) (tangential component of the defect orientation, psi = defects.dat angle, phi = polar angle). Fluctuations (mean
removed per run) -> normalised cross-correlation corr(T(t), C(t + lag)), pooled over runs; lag > 0: tilt leads."""
import gzip, json, os, sys, numpy as np
sys.path.insert(0, ".")
from edge_criteria import read_with_angle
lags = range(-3, 4); acc = {l: [] for l in lags}; nrev = 0; nint = 0
for s in ("s02", "s03"):
    for k in range(2, 11):
        d = f"Rscan_R31.83_seeds/{s}/k{k:02d}"; bc = json.load(open(d + "/input.json"))["BC"][0]; R = bc["R"]; Q = bc["Q"]
        t, fr = read_with_angle(d, np.array(Q[:2], float), 0.0)
        rows = [l.split() for l in gzip.open(d + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()]
        A = np.array(rows, float); ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
        x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ok = (r < R - 1) & (r > 0.5)
        C = np.mean((x[ok] * A[:, ok, 5] - y[ok] * A[:, ok, 4]) / r[ok], axis=1)
        T = np.full(len(ts), np.nan)
        for i, tf in enumerate(ts):
            vals = []
            for j in np.where((t > tf - 50) & (t <= tf))[0]:
                f = fr[j]; p = f[f[:, 2] > 0]; s_ = R - np.hypot(p[:, 0], p[:, 1]); w = s_ < 4
                if w.any(): vals.append(np.mean(np.sin(p[w, 3] - np.arctan2(p[w, 1], p[w, 0]))))
            if vals: T[i] = np.mean(vals)
        n0 = len(ts) // 5; C, T = C[n0:], T[n0:]
        good = np.isfinite(T); C, T = C[good], T[good]
        nrev += np.sum(np.diff(np.sign(C)) != 0); nint += len(C)
        c, tt = C - C.mean(), T - T.mean()
        for l in lags:
            a, b = (tt[:len(tt) - l], c[l:]) if l >= 0 else (tt[-l:], c[:len(c) + l])
            acc[l].append(np.corrcoef(a, b)[0, 1])
print(f"{nint} flow intervals, {nrev} circulation sign changes")
print("lag (time units; > 0 = tilt leads):  " + "  ".join(f"{50 * l:+5d}" for l in lags))
print("corr(tilt(t), circ(t+lag)), median:  " + "  ".join(f"{np.nanmedian(acc[l]):+5.2f}" for l in lags))
print("                    mean +- s.e.:    " + "  ".join(f"{np.nanmean(acc[l]):+5.2f}±{np.nanstd(acc[l]) / np.sqrt(len(acc[l])):.2f}" for l in lags))
