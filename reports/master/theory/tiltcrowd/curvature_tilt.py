#!/usr/bin/env python3
"""Curvature (1/R) correction of the tilt law at matched defect spacing (user, 2026-11-23). Flat-wall argument: tau = arccos(1/lambda) with strain = vorticity
at the no-slip wall. In a circle, the azimuthal flow at depth d below the wall has strain/vorticity |E|/W ~ (1 + d/R)/(1 - d/R) ~ 1 + 2d/R, so
cos(2 theta_L) = (W/|E|)/lambda and the tilt rises by dtau ~ (2d/R)/sqrt(lambda^2 - 1) rad (lambda 2: 66.2 deg per unit d/R). Crowding lowers the tilt at small
spacing, and at fixed k a small circle is also a crowded one, so the test is done at FIXED wall density rho = 2k/P (fixed spacing): R = k/(pi rho).
Data: defect_flow_main10.txt (main active sweep, lambda 2, per run: spacing a, depth d, |E|/W at the defect and 2 cells deeper, measured tilt).
Fits: (1) per rho, OLS tilt vs d/R across k (seeds pooled); (2) pooled with one offset per rho (fixed effects): tilt = alpha_rho + beta d/R, and
|E|/W (2 cells deeper) = a_rho + b (d+2)/R (prediction b = 2). Only rho with >= 3 k values and R/l_a in the single-vortex window 0.8-2.2 (l_a 22.21).
-> curvature_tilt.txt"""
import os, re, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21; LAM = 2.0
rows = []
for line in open(os.path.join(H, "defect_flow_main10.txt")):
    m = re.match(r"rho([\d.]+)/k(\d+)/s(\d+)\s*\|", line)
    if not m: continue
    f = [x.strip() for x in line.split("|")]; rho, k = float(m.group(1)), int(m.group(2)); a, d = float(f[2]), float(f[3])
    at, deep, tilt = f[4].split(), f[5].split(), float(f[6]); R = k / (np.pi * rho)
    rows.append(dict(rho=rho, k=k, a=a, d=d, R=R, r0=float(at[4]), r2=float(deep[0]), tilt=tilt))
X = [r for r in rows if 0.8 <= r["R"] / LA <= 2.2 and np.isfinite(r["tilt"])]
pred = np.degrees(2 / np.sqrt(LAM ** 2 - 1))
L = ["CURVATURE TEST OF THE TILT LAW AT FIXED WALL DENSITY (main active sweep, lambda 2; R = k/(pi rho); window 0.8 <= R/l_a <= 2.2)",
     f"prediction (no-slip azimuthal flow at depth d): dtau = {pred:.1f} deg x d/R; |E|/W at depth d' = 1 + 2 d'/R",
     "rho | spacing | k values (R/l_a) | runs | mean depth | tilt by k (mean over seeds) | slope dtilt/d(d/R) (deg) +- s.e. | |E|/W (2 cells deeper) by k"]
rhos = sorted({r["rho"] for r in X}); FE = []
for rho in rhos:
    g = [r for r in X if r["rho"] == rho]; ks = sorted({r["k"] for r in g})
    if len(ks) < 3: continue
    x = np.array([r["d"] / r["R"] for r in g]); y = np.array([r["tilt"] for r in g])
    A = np.c_[np.ones_like(x), x]; c, res, *_ = np.linalg.lstsq(A, y, rcond=None); s2 = np.sum((y - A @ c) ** 2) / max(len(y) - 2, 1)
    se = np.sqrt(s2 * np.linalg.inv(A.T @ A)[1, 1]); FE += g
    byk = " ".join(f"{k}:{np.mean([r['tilt'] for r in g if r['k'] == k]):.1f}" for k in ks)
    byr = " ".join(f"{k}:{np.mean([r['r2'] for r in g if r['k'] == k]):.2f}" for k in ks)
    L.append(f"{rho:.2f} | {np.mean([r['a'] for r in g]):5.1f} | {ks[0]}-{ks[-1]} ({min(r['R'] for r in g) / LA:.2f}-{max(r['R'] for r in g) / LA:.2f}) | {len(g):3d} | "
             f"{np.mean([r['d'] for r in g]):.2f} | {byk} | {c[1]:+7.1f} +- {se:5.1f} | {byr}")
# pooled fixed-effects fits
rl = sorted({r["rho"] for r in FE}); D = np.array([[1.0 * (r["rho"] == q) for q in rl] for r in FE])
for name, yk, xf, p in (("tilt", "tilt", lambda r: r["d"] / r["R"], pred), ("|E|/W 2 cells deeper", "r2", lambda r: (r["d"] + 2) / r["R"], 2.0),
                        ("|E|/W at the defect", "r0", lambda r: r["d"] / r["R"], 2.0)):
    x = np.array([xf(r) for r in FE]); y = np.array([r[yk] for r in FE]); A = np.c_[D, x]; c, *_ = np.linalg.lstsq(A, y, rcond=None)
    s2 = np.sum((y - A @ c) ** 2) / (len(y) - A.shape[1]); se = np.sqrt(s2 * np.linalg.inv(A.T @ A)[-1, -1])
    L.append(f"POOLED (one offset per rho, {len(rl)} densities, {len(FE)} runs): {name} slope = {c[-1]:+.2f} +- {se:.2f} per unit (predicted {p:+.2f}); "
             f"range of the regressor {x.min():.3f}-{x.max():.3f} -> predicted span {p * (x.max() - x.min()):.2f}, fitted span {c[-1] * (x.max() - x.min()):.2f}")
open(os.path.join(H, "curvature_tilt.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
