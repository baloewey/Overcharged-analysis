#!/usr/bin/env python3
"""Lobe rotation statistics for Cassini runs: mean azimuthal velocity about the centroid of each half
(x < Q_x, x > Q_x), after dropping the first 20% of flow frames. Usage: lobe_rotation.py RUN [RUN ...]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_cassini as rc
print(f"{'run':<22} {'b/a':>5} {'H':>3} {'waist':>5}  {'<v_phi> L':>9} {'R':>7}  {'|v_phi|/|v| L':>13} {'R':>5}  {'co-rot':>6}  {'reversals L':>11} {'R':>3}  {'whole':>6}")
for run in sys.argv[1:]:
    g = rc.geom(run); Lx, Ly = g["Lx"], g["Ly"]
    tf, F = rc.read_grid_frames(os.path.join(run, "flowfield.dat"), 7, Lx, Ly)
    xc, yc = np.meshgrid(np.arange(Lx) + 0.5, np.arange(Ly) + 0.5, indexing="ij"); mask = rc.inside(g, xc, yc)
    L = mask & (xc < g["Q"][0]); R = mask & (xc > g["Q"][0])
    def vphi(m, sel):
        c = np.array([xc[sel].mean(), yc[sel].mean()]); X, Y = xc[sel] - c[0], yc[sel] - c[1]; r = np.hypot(X, Y); ok = r > 1
        return np.mean((X[ok] * F[m, :, :, 1][sel][ok] - Y[ok] * F[m, :, :, 0][sel][ok]) / r[ok])
    s = range(len(tf) // 5, len(tf))
    vl = np.array([vphi(m, L) for m in s]); vr = np.array([vphi(m, R) for m in s]); vw = np.array([vphi(m, mask) for m in s])
    sp = np.array([np.mean(np.hypot(F[m, :, :, 0][mask], F[m, :, :, 1][mask])) for m in s])
    rev = lambda v: int(np.sum(np.diff(np.sign(v[np.abs(v) > 0.05])) != 0))
    ba = g["b"] / g["a"] if g["a"] else float("inf"); H = np.hypot(g["a"], g["b"]); waist = 2 * np.sqrt(max(g["b"] ** 2 - g["a"] ** 2, 0))
    print(f"{os.path.basename(os.path.normpath(run)):<22} {ba:5.2f} {H:3.0f} {waist:5.1f}  {vl.mean():+9.3f} {vr.mean():+7.3f}  {np.mean(np.abs(vl))/sp.mean():13.2f} {np.mean(np.abs(vr))/sp.mean():5.2f}  {np.mean(np.sign(vl)==np.sign(vr)):6.2f}  {rev(vl):11d} {rev(vr):3d}  {vw.mean():+6.3f}")
