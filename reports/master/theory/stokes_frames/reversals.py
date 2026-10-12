#!/usr/bin/env python3
"""At vortex reversals (activity 0.05, k = 9, 10, runs with flow every 50 time units): which flips first, circulation or wall tilt?
Circulation C per flow interval (value assigned to the interval centre, t - 25); zero crossing by linear interpolation.
Wall tilt T(t) = mean sin(psi - phi) over wall +1/2 (s < 4) per defect frame (2 time units), smoothed over 10 time units; its zero
crossing closest to the circulation crossing within +-200 time units. dt = t_tilt - t_circ (< 0: tilt flips first)."""
import glob, gzip, json, sys, numpy as np
sys.path.insert(0, ".")
from edge_criteria import read_with_angle
dts = []
runs = sorted(set(glob.glob("H3pilot/act0.05/k09") + glob.glob("H3pilot/act0.05/k10") + glob.glob("H3pilot_seeds/act0.05/s0*/k09") + glob.glob("H3pilot_seeds/act0.05/s0*/k10")))
for d in runs:
    bc = json.load(open(d + "/input.json"))["BC"][0]; R = bc["R"]; Q = bc["Q"]
    t, fr = read_with_angle(d, np.array(Q[:2], float), 0.0)
    op = gzip.open(d + "/flowfield.dat.gz", "rt") if glob.glob(d + "/flowfield.dat.gz") else open(d + "/flowfield.dat")
    rows = [l.split() for l in op if len(l.split()) == 7 and l.split()[0][0].isdigit()]
    A = np.array(rows, float); ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ok = (r < R - 1) & (r > 0.5)
    C = np.mean((x[ok] * A[:, ok, 5] - y[ok] * A[:, ok, 4]) / r[ok], axis=1); tc = ts - 25
    T = np.full(len(t), np.nan)
    for j, f in enumerate(fr):
        p = f[f[:, 2] > 0]; s_ = R - np.hypot(p[:, 0], p[:, 1]); w = s_ < 4
        if w.any(): T[j] = np.mean(np.sin(p[w, 3] - np.arctan2(p[w, 1], p[w, 0])))
    good = np.isfinite(T); tt, TT = t[good], T[good]
    TT = np.convolve(TT, np.ones(5) / 5, "same")
    tz = tt[:-1][np.sign(TT[:-1]) != np.sign(TT[1:])]
    nrev = 0
    for i in range(1, len(C) - 1):
        if np.sign(C[i]) != np.sign(C[i + 1]) and abs(C[i]) > 0.05 and abs(C[i + 1]) > 0.05:
            t0 = tc[i] + (tc[i + 1] - tc[i]) * C[i] / (C[i] - C[i + 1]); nrev += 1
            cand = tz[np.abs(tz - t0) < 200]
            if len(cand): dts.append(cand[np.argmin(np.abs(cand - t0))] - t0)
    print(f"{d}: {nrev} clear reversals")
dts = np.array(dts)
print(f"{len(dts)} reversals matched: t_tilt - t_circ median {np.median(dts):+.1f} time units, mean {dts.mean():+.1f} +- {dts.std() / np.sqrt(len(dts)):.1f};"
      f" tilt first in {np.mean(dts < 0):.0%}, circulation first in {np.mean(dts > 0):.0%}")
print("values:", np.round(np.sort(dts)).astype(int))
