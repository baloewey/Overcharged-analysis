#!/usr/bin/env python3
"""Error bars for the response test (2026-10-21): kresponse.py profiles split into 4 time blocks (t >= 200 after the first frame); fits with the bulk
angle free and with it fixed at the field angle (10 deg; valid when xi << half-width 19.5). Reports xi mean +- sd/sqrt(4) per run and
K1/K3 = (xi_planar/xi_homeo)^2 with propagated error. H 4: xi ~ half-width, so only the fixed-bulk fit is meaningful (the free fit trades bulk angle against xi).
-> kresponse_blocks.txt"""
import os, gzip, numpy as np
RD = "/home/baloewe/Physics/Overcharged/runs/Kresponse"; T0, NB = 200.0, 4
def profiles(run, geom):
    f = os.path.join(run, "directorfield.dat.gz"); rows = []
    with gzip.open(f, "rt") as fh:
        for line in fh:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            rows.append((float(p[0]), int(p[2]), np.arctan2(float(p[5]), float(p[4]))))
    a = np.array(rows); t = a[:, 0] - a[0, 0]; keep = t >= T0; a, t = a[keep], t[keep]
    edges = np.linspace(t.min(), t.max() + 1e-9, NB + 1); out = []
    for b in range(NB):
        s = (t >= edges[b]) & (t < edges[b + 1]); y = a[s, 1].astype(int); ang = a[s, 2]
        c2 = np.bincount(y, np.cos(2 * ang), 42); s2 = np.bincount(y, np.sin(2 * ang), 42); th = 0.5 * np.arctan2(s2, c2)
        if geom == "homeo": th = np.angle(np.exp(2j * (th - np.pi / 2))) / 2
        out.append(np.degrees(th))
    return out
def fit(th, fixed=None):
    y = np.arange(42); sel = (y >= 2) & (y <= 39); best = None
    for xi in np.arange(0.5, 40.0, 0.01):
        g = np.cosh((y[sel] - 20.5) / xi) / np.cosh(19.5 / xi)
        if fixed is None:
            X = np.column_stack([1 - g, g]); c, *_ = np.linalg.lstsq(X, th[sel], rcond=None); r = np.mean((X @ c - th[sel]) ** 2)
        else:
            yy = th[sel] - fixed * (1 - g); c, *_ = np.linalg.lstsq(g[:, None], yy, rcond=None); r = np.mean((g * c[0] - yy) ** 2)
        if best is None or r < best[0]: best = (r, xi)
    return best[1]
L = ["run | xi free (mean +- se over 4 blocks) | xi with bulk fixed at the field angle"]; R = {}
for coup in ("prod", "susc0"):
    for geom in ("planar", "homeo"):
        for H in (4, 8, 12, 16):
            P = profiles(f"{RD}/{coup}_{geom}_H{H}", geom); sgn = 1 if geom == "planar" else -1
            xf = np.array([fit(p) for p in P]); xx = np.array([fit(p, 10.0 * sgn) for p in P])
            R[(coup, geom, H)] = (xf, xx)
            L.append(f"{coup}_{geom}_H{H:<2d} | {xf.mean():5.2f} +- {xf.std(ddof=1)/2:4.2f} | {xx.mean():5.2f} +- {xx.std(ddof=1)/2:4.2f}   blocks {np.round(xx, 2)}")
L.append("\nK1/K3 (response) = (xi_splay/xi_bend)^2  [free fit | fixed-bulk fit]")
for coup in ("prod", "susc0"):
    for H in (4, 8, 12, 16):
        out = []
        for k in (0, 1):
            a, b = R[(coup, "planar", H)][k], R[(coup, "homeo", H)][k]
            r = (a.mean() / b.mean()) ** 2; e = r * 2 * np.hypot(a.std(ddof=1) / 2 / a.mean(), b.std(ddof=1) / 2 / b.mean()); out.append(f"{r:.2f} +- {e:.2f}")
        L.append(f"{coup} H {H:2d}: {out[0]} | {out[1]}")
open("kresponse_blocks.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
