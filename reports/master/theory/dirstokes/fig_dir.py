#!/usr/bin/env python3
"""Figure for the director-field Stokes model: (a) circulation per interval, measured vs predicted (total, one fitted c per run), circles and
ellipses; (b) ellipses: two-vortex (antisymmetric) part of the half-cavity rotations, measured vs bulk-only and wall-only predictions;
(c) example interval (eps 1.5, two-vortex): measured flow and predicted flow (streamlines)."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(H, "fit_dir.py")).read().split("for tag in (")[0].split("H = os.path")[0])
def circ(x, y, v): r = np.hypot(x, y); return np.mean((x * v[..., 1] - y * v[..., 0]) / np.maximum(r, 0.5), axis=-1)
def halves(x, y, v):
    out = []
    for m in (x < 0, x > 0):
        cx, cy = x[m].mean(), y[m].mean(); dx, dy = x[m] - cx, y[m] - cy
        out.append(np.mean((dx * v[..., m, 1] - dy * v[..., m, 0]) / np.maximum(np.hypot(dx, dy), 0.5), axis=-1))
    return np.stack(out, -1)
def load(tag):
    d = np.load(os.path.join(H, tag + ".npz")); nint, nc = d["vx"].shape
    uw = np.loadtxt(os.path.join(H, tag + "_uw.csv"), delimiter=",").reshape(nint, nc, 2); ub = np.loadtxt(os.path.join(H, tag + "_ub.csv"), delimiter=",").reshape(nint, nc, 2)
    meas = np.stack([d["vx"], d["vy"]], -1); u = uw + ub; c = np.sum(u * meas) / np.sum(u * u)
    return d, meas, uw, ub, c
fig = plt.figure(figsize=(13, 4)); ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2)
for tag, col in (("k02", "C0"), ("k06", "C1"), ("k10", "C2"), ("e12_s01", "C4"), ("e15_s01", "C3"), ("e20_s01", "C5")):
    d, meas, uw, ub, c = load(tag); ax1.plot(c * circ(d["x"], d["y"], uw + ub), circ(d["x"], d["y"], meas), "o", ms=3, color=col, label=tag)
ax1.axline((0, 0), slope=1, color="k", lw=0.6); ax1.set_xlabel("predicted circulation (full active force)"); ax1.set_ylabel("measured circulation")
ax1.set_title("(a) circulation per interval", fontsize=10); ax1.legend(fontsize=7)
for tag in ("e15_s01", "e15_s02", "e20_s01", "e20_s02"):
    d, meas, uw, ub, c = load(tag); x, y = d["x"], d["y"]; Hm = halves(x, y, meas); am = (Hm[:, 0] - Hm[:, 1]) / 2
    for U, mk, lab in ((ub, "o", "bulk"), (uw, "x", "wall band")):
        Hp = halves(x, y, U); ap = (Hp[:, 0] - Hp[:, 1]) / 2; cc = np.sum(ap * am) / np.sum(ap * ap)
        ax2.plot(cc * ap, am, mk, ms=3, color="C3" if lab == "bulk" else "gray", label=lab if tag == "e15_s01" else None)
ax2.axline((0, 0), slope=1, color="k", lw=0.6); ax2.set_xlabel("predicted two-vortex part (own fitted scale)"); ax2.set_ylabel("measured two-vortex part")
ax2.set_title("(b) ellipses 1.5, 2.0: two-vortex part", fontsize=10); ax2.legend(fontsize=7)
d, meas, uw, ub, c = load("e15_s01"); x, y = d["x"], d["y"]; Hm = halves(x, y, meas); n = int(np.argmax(-(Hm[:, 0] * Hm[:, 1])))
ax3 = fig.add_subplot(1, 3, 3); A, B = float(d["A"]), float(d["B"])
from scipy.interpolate import griddata
gx, gy = np.meshgrid(np.linspace(-A, A, 80), np.linspace(-B, B, 50))
for U, col, lab, off in ((meas[n], "k", "measured", 0), (c * (uw[n] + ub[n]), "C3", "full active force", 0)):
    vx = griddata((x, y), U[:, 0], (gx, gy)); vy = griddata((x, y), U[:, 1], (gx, gy))
    ax3.streamplot(gx, gy, np.nan_to_num(vx), np.nan_to_num(vy), color=col, density=0.8, linewidth=0.7)
t = np.linspace(0, 2 * np.pi, 200); ax3.plot(A * np.cos(t), B * np.sin(t), "k-", lw=1); ax3.set_aspect("equal")
ax3.set_title(f"(c) eps 1.5, interval {n}: measured (black) vs model (red)", fontsize=9); ax3.axis("off")
fig.tight_layout(); fig.savefig(os.path.join(H, "dirstokes.png"), dpi=150); print("ok")
