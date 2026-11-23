#!/usr/bin/env python3
"""Is the crowding drop of the tilt a threshold or a continuous fall? (user, 2026-11-23). Data: defect_flow_main10.txt (active circles, lambda 2, 10 seeds,
per run: spacing a, tilt), R/l_a 0.8-2.2. Models (least squares on run values): (H) hinge: T0 for a >= a*, T0 - b (a* - a) below; (L) linear in 1/a;
(E) T0 - B exp(-a/a0). Compare residuals (and BIC); bootstrap over runs for a*. Also the spacing binned means. -> crowd_onset.txt"""
import os, re, numpy as np
from scipy.optimize import curve_fit
H = os.path.dirname(os.path.abspath(__file__)); a, t, R = [], [], []
for line in open(os.path.join(H, "defect_flow_main10.txt")):
    m = re.match(r"rho([\d.]+)/k(\d+)/s(\d+)\s*\|", line)
    if not m: continue
    f = [x.strip() for x in line.split("|")]; k = int(m.group(2)); rr = k / (np.pi * float(m.group(1)))
    if 0.8 <= rr / 22.21 <= 2.2: a.append(float(f[2])); t.append(float(f[6])); R.append(rr)
a, t = np.array(a), np.array(t); ok = np.isfinite(t); a, t = a[ok], t[ok]; n = len(a)
hinge = lambda x, T0, b, s: T0 - b * np.maximum(s - x, 0)
lin = lambda x, T0, c: T0 - c / x
ex = lambda x, T0, B, a0: T0 - B * np.exp(-x / a0)
L = [f"CROWDING DROP: THRESHOLD OR CONTINUOUS? ({n} active circles, lambda 2, R/l_a 0.8-2.2)", "binned means (spacing a in cells: tilt +- s.e., runs):"]
for lo, hi in ((5, 7), (7, 9), (9, 11), (11, 13), (13, 15), (15, 17), (17, 20), (20, 24), (24, 30), (30, 40), (40, 70)):
    m = (a >= lo) & (a < hi)
    if m.sum() >= 3: L.append(f"  {lo:2d}-{hi:2d}: {t[m].mean():5.1f} +- {t[m].std() / np.sqrt(m.sum()):.1f} ({m.sum()})")
res = {}
for name, f, p0 in (("hinge", hinge, (64.5, 1.5, 18)), ("1/a", lin, (66, 100)), ("exp", ex, (65, 40, 6))):
    p, _ = curve_fit(f, a, t, p0=p0, maxfev=20000); rss = np.sum((t - f(a, *p)) ** 2); k = len(p); bic = n * np.log(rss / n) + k * np.log(n)
    res[name] = (p, rss, bic); L.append(f"{name:6s}: params {np.round(p, 2).tolist()} | rms {np.sqrt(rss / n):.2f} deg | BIC {bic:.1f}")
rng = np.random.default_rng(1); st = []
for _ in range(500):
    i = rng.integers(0, n, n)
    try: st.append(curve_fit(hinge, a[i], t[i], p0=res["hinge"][0], maxfev=20000)[0][2])
    except Exception: pass
st = np.array(st); L.append(f"hinge onset a* = {res['hinge'][0][2]:.1f} cells (bootstrap 68%: {np.percentile(st, 16):.1f}-{np.percentile(st, 84):.1f}); plateau {res['hinge'][0][0]:.1f} deg, slope {res['hinge'][0][1]:.2f} deg/cell")
L.append(f"candidate lengths (chi 0.5): l_a 22.2; 2 pi l_E with l_E 2.3-4.4 -> 14-28; core depth d 1.2-1.8; 2 pi d -> 8-11; registry stiffness collapse (crowd_elastic) between k 6 and 8 at R 31.8 -> spacing 20-14")
open(os.path.join(H, "crowd_onset.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
