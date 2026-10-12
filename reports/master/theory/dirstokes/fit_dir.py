#!/usr/bin/env python3
"""Measured flow vs Stokes flow of the full measured active force (-div Q), wall band + bulk (solve_dir.wl).
Per run: (a) one coefficient c for the total: v = c (u_w + u_b); expected c ~ zeta_eff/eta ~ 4/17 ~ 0.24 (Sec. frank) if the calibration
holds (no density factor). (b) separate c_w, c_b. Reported: R^2 (fraction of the measured velocity variance, over cells and intervals),
correlation of the cavity circulation per interval; ellipses: half-cavity rotations, one-vortex (symmetric) and two-vortex
(antisymmetric) parts, and the vortex number per interval (halves co-/counter-rotating) predicted by the total, the wall band alone and
the bulk alone (each with its own fitted c). Writes fit_dir.txt."""
import os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
def circ(x, y, v): r = np.hypot(x, y); return np.mean((x * v[..., 1] - y * v[..., 0]) / np.maximum(r, 0.5), axis=-1)
def halves(x, y, v):
    out = []
    for m in (x < 0, x > 0):
        cx, cy = x[m].mean(), y[m].mean(); dx, dy = x[m] - cx, y[m] - cy
        out.append(np.mean((dx * v[..., m, 1] - dy * v[..., m, 0]) / np.maximum(np.hypot(dx, dy), 0.5), axis=-1))
    return np.stack(out, -1)
def fit(Xs, y):
    c, *_ = np.linalg.lstsq(np.column_stack(Xs), y, rcond=None); r2 = 1 - np.sum((np.column_stack(Xs) @ c - y) ** 2) / np.sum((y - 0) ** 2); return c, r2
for tag in ("k02", "k06", "k10", "e12_s01", "e12_s02", "e15_s01", "e15_s02", "e20_s01", "e20_s02"):
    d = np.load(os.path.join(H, tag + ".npz")); x, y = d["x"], d["y"]; nint, nc = d["vx"].shape
    meas = np.stack([d["vx"], d["vy"]], -1)
    uw = np.loadtxt(os.path.join(H, tag + "_uw.csv"), delimiter=",").reshape(nint, nc, 2)
    ub = np.loadtxt(os.path.join(H, tag + "_ub.csv"), delimiter=",").reshape(nint, nc, 2)
    yv = meas.ravel()
    (c,), r2 = fit([(uw + ub).ravel()], yv); (cw, cb), r2s = fit([uw.ravel(), ub.ravel()], yv)
    (cwo,), r2w = fit([uw.ravel()], yv); (cbo,), r2b = fit([ub.ravel()], yv)
    cm = circ(x, y, meas); ct = circ(x, y, c * (uw + ub)); cwc = circ(x, y, uw); cbc = circ(x, y, ub)
    say(f"{tag:8s}: total c = {c:+.3f} (R2 {r2:.2f}); separate c_wall {cw:+.3f}, c_bulk {cb:+.3f} (R2 {r2s:.2f}); alone: wall R2 {r2w:.2f} (c {cwo:+.3f}), "
        f"bulk R2 {r2b:.2f} (c {cbo:+.3f}); circulation corr: total {np.corrcoef(ct, cm)[0, 1]:+.2f}, wall {np.corrcoef(cwc, cm)[0, 1]:+.2f}, bulk {np.corrcoef(cbc, cm)[0, 1]:+.2f}; "
        f"circulation share wall/bulk (|.| mean) {np.mean(np.abs(cwc)) / (np.mean(np.abs(cwc)) + np.mean(np.abs(cbc))):.2f}")
    if tag.startswith("e"):
        Hm = halves(x, y, meas); two = Hm[:, 0] * Hm[:, 1] < 0
        line = f"          measured two-vortex in {two.mean():.0%} of intervals;"
        for lab, U in (("total", uw + ub), ("wall", uw), ("bulk", ub)):
            Hp = halves(x, y, U); s_p, a_p = (Hp[:, 0] + Hp[:, 1]) / 2, (Hp[:, 0] - Hp[:, 1]) / 2; s_m, a_m = (Hm[:, 0] + Hm[:, 1]) / 2, (Hm[:, 0] - Hm[:, 1]) / 2
            pt = Hp[:, 0] * Hp[:, 1] < 0
            line += (f" {lab}: corr one-vortex part {np.corrcoef(s_p, s_m)[0, 1]:+.2f}, two-vortex part {np.corrcoef(a_p, a_m)[0, 1]:+.2f}, "
                     f"vortex number right {np.mean(pt == two):.0%};")
        say(line)
open(os.path.join(H, "fit_dir.txt"), "w").write("\n".join(L) + "\n")
