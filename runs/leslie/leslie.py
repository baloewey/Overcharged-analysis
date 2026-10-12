#!/usr/bin/env python3
"""Bulk Leslie angle in the MPCD nematic (Poiseuille channel, walls at y = 4 and 64, no anchoring), production second half.
Per row y: director angle from the flow direction x (S-weighted nematic mean), S, and the local shear rate dv_x/dy (flowfield).
The aligned angle has the sign of the shear: theta = +theta_L (dv/dy > 0, lower half), -theta_L (upper half). Rows with |dv/dy| above half its maximum are pooled
(theta * sign(dv/dy)). Predictions: bare 0.5 arccos(1/lam); particle-averaged 0.5 arccos(1/lam_eff), lam_eff = lam (1 + S4)/(2 S2), Maier-Saupe S4 at the measured S."""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt")
from mpcd_profile import S4_of
def frames(fn, ncol, t0):
    cur, tc = [], None
    import gzip
    for line in (gzip.open(fn + ".gz", "rt") if not os.path.exists(fn) else open(fn)):
        p = line.split()
        if len(p) != ncol or not p[0][0].isdigit(): continue
        t = float(p[0])
        if t != tc:
            if cur and tc >= t0: yield tc, np.array(cur)
            cur, tc = [], t
        if t >= t0: cur.append([float(x) for x in p])
    if cur and tc >= t0: yield tc, np.array(cur)
L = []
for run in sorted(glob.glob("lam*_U*")):
    if not os.path.exists(run + "/DONE"): continue
    J = json.load(open(run + "/input.json")); lam = J["species"][0]["tumble"]; U = J["species"][0]["mfpot"]; T = J["simSteps"] * 0.1; t0 = T / 2
    z, Ss, n = {}, {}, {}
    for t, c in frames(run + "/directorfield.dat", 8, t0):
        y = c[:, 2].astype(int); S = c[:, 7]; th = np.arctan2(c[:, 5], c[:, 4]); ok = S > 0
        for yy in np.unique(y[ok]): m = ok & (y == yy); z[yy] = z.get(yy, 0) + np.sum(S[m] * np.exp(2j * th[m])); Ss[yy] = Ss.get(yy, 0) + S[m].sum(); n[yy] = n.get(yy, 0) + m.sum()
    v = {}
    for t, c in frames(run + "/flowfield.dat", 7, t0):
        y = c[:, 2].astype(int)
        for yy in np.unique(y): v.setdefault(yy, []).append(c[y == yy, 4].mean())
    ys = np.array(sorted(set(z) & set(v))); vm = np.array([np.mean(v[yy]) for yy in ys]); g = np.gradient(vm, ys)
    th = np.array([np.angle(z[yy]) / 2 for yy in ys]); S = np.array([Ss[yy] / n[yy] for yy in ys])
    sel = (np.abs(g) > 0.5 * np.abs(g).max()) & (ys > 6) & (ys < 62)
    thL = np.degrees(np.mean(th[sel] * np.sign(g[sel]))); Sb = S[sel].mean(); le = lam * (1 + S4_of(min(Sb, .995))) / (2 * Sb)
    bare = np.degrees(0.5 * np.arccos(1 / lam)); pav = np.degrees(0.5 * np.arccos(1 / le)) if le > 1 else float("nan")
    L.append(f"{run}: lam {lam} U {U}: measured |theta_L| {thL:5.1f} deg (rows {sel.sum()}, S {Sb:.3f}, max shear {np.abs(g).max():.4f}) | bare {bare:5.1f} | lam_eff {le:.2f} -> {pav:5.1f}")
    L.append("   rows y: theta*sign(shear) " + " ".join(f"{yy}:{np.degrees(a * np.sign(b)):+.0f}" for yy, a, b in zip(ys[::4], th[::4], g[::4])))
print("\n".join(L)); open("leslie.txt", "w").write("\n".join(L) + "\n")
