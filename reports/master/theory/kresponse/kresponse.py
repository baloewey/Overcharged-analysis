#!/usr/bin/env python3
"""Static response test of the Frank anisotropy (runs/Kresponse; 2026-10-21).
Passive channel, walls at y = 1 and 41 (fluid rows 1..40), periodic in x (L 64). Planar walls + field 10 deg from x: director angle theta(y) from x
heals from the wall value to the field-imposed bulk value; with the director ~x the gradient is along y = SPLAY. Homeotropic walls + field 10 deg
from y: theta measured from y; gradient along the director = BEND.
Profile: time- and x-averaged Q tensor per row (frames t >= T0), theta_row = 1/2 atan2(<sin 2a>, <cos 2a>) - axis. Fit
    theta(y) = theta_b + (theta_w - theta_b) cosh((y - y_c)/xi) / cosh(h/xi)
(theta_b bulk, theta_w apparent wall value, h = 19.5 half-width, y_c = 20.5 centre; numpy only: linear in theta_b, theta_w, grid in xi), rows 2..39 (first fluid row excluded). Same field in both geometries ->
K1/K3 (response) = (xi_planar / xi_homeo)^2. Coupling on (prod) vs off (susc0: control, expected 1).
Usage: kresponse.py [RUNS_DIR] -> kresponse.txt, kresponse_profiles.npz"""
import sys, os, gzip, numpy as np
RD = sys.argv[1] if len(sys.argv) > 1 else "/home/baloewe/Physics/Overcharged/runs/Kresponse"; T0 = 200.0
def profile(run, geom):
    f = os.path.join(run, "directorfield.dat"); op = gzip.open(f + ".gz", "rt") if os.path.exists(f + ".gz") else open(f)
    c2 = np.zeros(42); s2 = np.zeros(42); n = np.zeros(42); t0 = None
    with op as fh:
        for line in fh:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            t = float(p[0]); t0 = t if t0 is None else t0
            if t - t0 < T0: continue
            y = int(p[2]); a = np.arctan2(float(p[5]), float(p[4])); c2[y] += np.cos(2 * a); s2[y] += np.sin(2 * a); n[y] += 1
    ok = n > 0; th = np.full(42, np.nan); th[ok] = 0.5 * np.arctan2(s2[ok], c2[ok])
    if geom == "homeo": th = np.angle(np.exp(2j * (th - np.pi / 2))) / 2     # angle from y
    return np.degrees(th), n
def model(y, tb, tw, xi, yc, h=19.5): return tb + (tw - tb) * np.cosh((y - yc) / xi) / np.cosh(h / xi)
out, prof = ["run | bulk angle (deg) | apparent wall angle | xi (cells) | centre | rms residual (deg) | frames-rows"], {}
res = {}
for coup in ("prod", "susc0"):
    for geom in ("planar", "homeo"):
        for H in (4, 8, 12, 16):
            run = f"{RD}/{coup}_{geom}_H{H}"
            if not (os.path.exists(run + "/directorfield.dat") or os.path.exists(run + "/directorfield.dat.gz")): continue
            th, n = profile(run, geom); y = np.arange(42); sel = (y >= 2) & (y <= 39) & np.isfinite(th)
            best = None
            for xi in np.arange(0.5, 60.0, 0.02):
                g = np.cosh((y[sel] - 20.5) / xi) / np.cosh(19.5 / xi); X = np.column_stack([1 - g, g])
                c, *_ = np.linalg.lstsq(X, th[sel], rcond=None); r = np.sqrt(np.mean((X @ c - th[sel]) ** 2))
                if best is None or r < best[0]: best = (r, c[0], c[1], xi)
            r, p = best[0], [best[1], best[2], best[3], 20.5]
            res[(coup, geom, H)] = p[2]; prof[f"{coup}_{geom}_H{H}"] = th
            out.append(f"{coup}_{geom}_H{H:<2d} | {p[0]:6.2f} | {p[1]:6.2f} | {abs(p[2]):6.2f} | {p[3]:5.2f} | {r:.2f} | {int(n[20])}")
out.append("\nK1/K3 (response) = (xi_planar/xi_homeo)^2; field scaling check: xi * H should be constant")
for coup in ("prod", "susc0"):
    for H in (4, 8, 12, 16):
        a, b = res.get((coup, "planar", H)), res.get((coup, "homeo", H))
        if a is not None and b is not None: out.append(f"{coup} H {H:2d}: xi splay {abs(a):6.2f} (xi*H {abs(a)*H:6.1f}), xi bend {abs(b):6.2f} (xi*H {abs(b)*H:6.1f}) -> K1/K3 = {(a/b)**2:.2f}")
open("kresponse.txt", "w").write("\n".join(out) + "\n"); np.savez("kresponse_profiles.npz", **prof); print("\n".join(out))
