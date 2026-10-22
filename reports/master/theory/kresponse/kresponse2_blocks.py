#!/usr/bin/env python3
"""Kresponse2 (W 80, H 4/6, 3 seeds; adapted from the Kresponse version). Error bars for the response test (2026-10-21): kresponse.py profiles split into 4 time blocks (t >= 200 after the first frame); fits with the bulk
angle free and with it fixed at the field angle (10 deg; valid when xi << half-width 19.5). Reports xi mean +- sd/sqrt(4) per run and
K1/K3 = (xi_planar/xi_homeo)^2 with propagated error. H 4: xi ~ half-width, so only the fixed-bulk fit is meaningful (the free fit trades bulk angle against xi).
-> kresponse_blocks.txt"""
import os, gzip, numpy as np
RD = "/home/baloewe/Physics/Overcharged/runs/Kresponse2"; T0, NB = 500.0, 4
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
        c2 = np.bincount(y, np.cos(2 * ang), 82); s2 = np.bincount(y, np.sin(2 * ang), 82); th = 0.5 * np.arctan2(s2, c2)
        if geom == "homeo": th = np.angle(np.exp(2j * (th - np.pi / 2))) / 2
        out.append(np.degrees(th))
    return out
def fit(th, fixed=None):
    y = np.arange(82); sel = (y >= 2) & (y <= 79); best = None
    for xi in np.arange(0.5, 60.0, 0.01):
        g = np.cosh((y[sel] - 40.5) / xi) / np.cosh(39.5 / xi)
        if fixed is None:
            X = np.column_stack([1 - g, g]); c, *_ = np.linalg.lstsq(X, th[sel], rcond=None); r = np.mean((X @ c - th[sel]) ** 2)
        else:
            yy = th[sel] - fixed * (1 - g); c, *_ = np.linalg.lstsq(g[:, None], yy, rcond=None); r = np.mean((g * c[0] - yy) ** 2)
        if best is None or r < best[0]: best = (r, xi)
    return best[1]
L = ["run | xi free (mean +- se over 4 blocks) | xi with bulk fixed at the field angle"]; R = {}
import glob
for coup in ("prod", "susc0"):
    for geom in ("planar", "homeo"):
        for H in (4, 6):
          for run in sorted(glob.glob(f"{RD}/{coup}_{geom}_H{H}_s*")):
            P = profiles(run, geom); H_ = H; sgn = 1 if geom == "planar" else -1
            xf = np.array([fit(p) for p in P]); xx = np.array([fit(p, 10.0 * sgn) for p in P])
            R.setdefault((coup, geom, H), []).append((xf, xx))
            L.append(f"{run.split('/')[-1]:22s} | {xf.mean():5.2f} +- {xf.std(ddof=1)/2:4.2f} | {xx.mean():5.2f} +- {xx.std(ddof=1)/2:4.2f}   blocks {np.round(xx, 2)}")
L.append("\nK1/K3 (response) = (xi_splay/xi_bend)^2  [free fit | fixed-bulk fit]")
L.append("(pooled over seeds x blocks; per-seed ratios in brackets)")
for coup in ("prod", "susc0"):
    for H in (4, 6):
        if (coup, "planar", H) not in R: continue
        out = []
        for k in (0, 1):
            a = np.concatenate([v[k] for v in R[(coup, "planar", H)]]); b = np.concatenate([v[k] for v in R[(coup, "homeo", H)]])
            r = (a.mean() / b.mean()) ** 2; e = r * 2 * np.hypot(a.std(ddof=1) / np.sqrt(len(a)) / a.mean(), b.std(ddof=1) / np.sqrt(len(b)) / b.mean())
            per = [(pa[k].mean() / pb[k].mean()) ** 2 for pa, pb in zip(R[(coup, "planar", H)], R[(coup, "homeo", H)])]
            out.append(f"{r:.2f} +- {e:.2f} [{', '.join(f'{x:.2f}' for x in per)}]")
        L.append(f"{coup} H {H:2d}: {out[0]} | {out[1]}")
open("kresponse2_blocks.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
