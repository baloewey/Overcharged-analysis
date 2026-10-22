#!/usr/bin/env python3
"""Candidate SI panels illustrating how we measure (user TODO, 2026-10-22).
(a,b) classification snapshots: wall band (d < 4 cells, dark), buffer (4-6), bulk; wall +1/2 (dark red, tail drawn), bulk +1/2 (light red), -1/2 (blue tripods);
small circle and sparse aspect-2 ellipse (cases as in the layer sheet). (c) depth histograms of +1/2 (time-averaged, per frame per cell of depth) for the four
cases, with the 4- and 6-cell cuts. (d) enclosed charge Q(r) against distance from the wall (circles, k 6, three radii; crowd_protect/fig1c/bands_all.npz):
plateau 1 inside the layer, jump to k at the wall; the 4-cell cut. (e) tilt definition: tail direction vs outward normal, tau. (f) an escape: depth of one
tracked wall +1/2 against time, with the 4-cell band and the 6-cell escape depth (sparse ellipse)."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.spatial import cKDTree
from common import CASES, load, tracks, LA
T = "/home/baloewe/Physics/Overcharged/reports/master/theory"; RED, LRED, BLUE = "#b2182b", "#f4a582", "#2166ac"
fig = plt.figure(figsize=(17, 9.5)); g = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.28, left=0.05, right=0.99, top=0.93, bottom=0.07)
def snapshot(ax, D, n):
    W = D["W"]; xs = np.linspace(W[:, 0].min() - 1, W[:, 0].max() + 1, 300); ys = np.linspace(W[:, 1].min() - 1, W[:, 1].max() + 1, 300); X, Y = np.meshgrid(xs, ys)
    ins = Path(W).contains_points(np.c_[X.ravel(), Y.ravel()]).reshape(X.shape); d = D["tree"].query(np.c_[X.ravel(), Y.ravel()])[0].reshape(X.shape)
    band = np.where(ins, np.where(d < 4, 2, np.where(d < 6, 1, 0)), np.nan)
    ax.contourf(X, Y, band, levels=[-0.5, 0.5, 1.5, 2.5], colors=["white", "#f0f0f0", "#d9d9d9"]); ax.plot(np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]], "k-", lw=1)
    f = D["fr"][n]; dd, jj = D["tree"].query(f[:, :2])
    for (x, y, q, a), dv, j in zip(f, dd, jj):
        if q > 0:
            c = RED if dv < 4 else LRED; ax.plot([x, x + 3 * np.cos(a)], [y, y + 3 * np.sin(a)], color=c, lw=1.8); ax.plot(x, y, "o", ms=4, color=c, mec="w", mew=0.5)
        else:
            for k in range(3): b = a + 2 * np.pi * k / 3; ax.plot([x, x + 2.2 * np.cos(b)], [y, y + 2.2 * np.sin(b)], color=BLUE, lw=1.4)
    ax.set_aspect("equal"); ax.axis("off")
    ax.plot([], [], "s", color="#d9d9d9", ms=9, label="Wall band, $d<4$"); ax.plot([], [], "s", color="#f0f0f0", ms=9, mec="0.7", label="Buffer, $4\\leq d<6$")
    ax.plot([], [], "-", color=RED, lw=2, label="Wall $+\\frac{1}{2}$ (tail)"); ax.plot([], [], "-", color=LRED, lw=2, label="Bulk $+\\frac{1}{2}$"); ax.plot([], [], "-", color=BLUE, lw=2, label="$-\\frac{1}{2}$")
DS = [load(r, 0.5) for _, r in CASES]
ax = fig.add_subplot(g[0, 0]); snapshot(ax, DS[0], len(DS[0]["fr"]) // 2); ax.set_title(CASES[0][0], fontsize=11); ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8.5, frameon=False)
ax = fig.add_subplot(g[0, 1]); snapshot(ax, DS[2], len(DS[2]["fr"]) // 2); ax.set_title(CASES[2][0], fontsize=11)
ax = fig.add_subplot(g[0, 2]); E = np.arange(0, 30.01, 0.5); CC = ("#08306b", "#4a90c2", "#e6550d", "#4d1a02")
for (lab, _), D, c in zip(CASES, DS, CC):
    h = np.zeros(len(E) - 1)
    for f in D["fr"]:
        if len(f): p = f[f[:, 2] > 0]; h += np.histogram(D["tree"].query(p[:, :2])[0], E)[0]
    ax.semilogy(0.5 * (E[:-1] + E[1:]), np.maximum(h / len(D["fr"]) / 0.5, 1e-3), "-", color=c, lw=1.6, label=lab)
for x, l in ((4, "wall cut"), (6, "escape depth")): ax.axvline(x, color="0.4", ls="--", lw=0.9); ax.text(x + 0.3, 0.012, l, fontsize=9, rotation=90, va="bottom", color="0.3")
ax.set_xlabel("Distance from the wall $d$ (cells)"); ax.set_ylabel(r"$+\frac{1}{2}$ per frame per cell of depth"); ax.legend(fontsize=8, loc="upper right"); ax.set_xlim(0, 30); ax.set_ylim(1e-2, 10)
ax = fig.add_subplot(g[1, 0]); Z = np.load(T + "/crowd_protect/fig1c/bands_all.npz"); Ed = Z["edges"]; dc = Ed[1:]
for R0, c in ((1.43, "#08306b"), (2.15, "#4a90c2"), (2.87, "#9ecae1")):
    m = (Z["k"] == 6) & (np.abs(Z["R"] - R0) < 0.02); q = 0.5 * (Z["Hp"][m] - Z["Hm"][m]).mean(0); enc = 6 - np.cumsum(q)   # charge farther than d from the wall
    ax.plot(np.r_[0, dc[:-1]], np.r_[6, enc[:-1]], "-", color=c, lw=1.8, label=f"$R={R0}\\,\\ell_a$")
ax.axhline(1, color="0.5", ls=":", lw=0.9); ax.axvline(4, color="0.4", ls="--", lw=0.9); ax.set_xlim(0, 25); ax.set_xlabel("Distance from the wall $d$ (cells)"); ax.set_ylabel("Charge farther than $d$ from the wall")
ax.legend(fontsize=9, title="Circles, $k=6$", title_fontsize=9); ax.text(4.3, 4.5, "4-cell cut", fontsize=9, color="0.3")
ax = fig.add_subplot(g[1, 1]); th = np.linspace(-0.5, 0.5, 50); ax.plot(np.sin(th) * 6, 6 * np.cos(th) - 6, "k-", lw=2)
x0, y0 = 0.0, -1.6; tau = np.radians(55); ax.annotate("", xy=(0, 1.2), xytext=(0, -0.0), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.5)); ax.text(0.15, 1.0, "outward normal", fontsize=10)
ax.plot([x0, x0 + 2.2 * np.sin(tau)], [y0, y0 + 2.2 * np.cos(tau)], color=RED, lw=3); ax.plot(x0, y0, "o", color=RED, ms=9, mec="w")
ax.plot([x0, x0], [y0, y0 + 2.0], ":", color="0.4"); aa = np.linspace(0, tau, 30); ax.plot(x0 + 0.9 * np.sin(aa), y0 + 0.9 * np.cos(aa), "-", color="0.2"); ax.text(0.35, y0 + 1.05, r"$\tau$", fontsize=14)
ax.text(1.95, -0.25, "tail", fontsize=10, color=RED); ax.text(-1.5, -2.6, "head", fontsize=10, color=RED); ax.annotate("", xy=(x0 - 1.3 * np.sin(tau), y0 - 1.3 * np.cos(tau)), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=RED))
ax.text(-3.0, 1.9, "Tilt $\\tau$: tail vs outward normal\nCoherence $C=|\\langle e^{i\\tau_j}\\rangle|$ over wall defects", fontsize=10)
ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 2.6); ax.set_aspect("equal"); ax.axis("off"); ax.set_title("Tilt definition", fontsize=11)
ax = fig.add_subplot(g[1, 2]); tr = tracks(DS[2]); dt = DS[2]["t"][1] - DS[2]["t"][0]; best = None
for v in tr.values():
    a = np.array(v); d = a[:, 3]
    if len(a) > 60 and (d[:30] < 4).all() and d.max() > 10: best = a; break
ax.plot((best[:, 0] - best[0, 0]) * dt, best[:, 3], "-", color=RED, lw=1.6); ax.axhspan(0, 4, color="#d9d9d9", zorder=0); ax.axhline(6, color="0.4", ls="--", lw=0.9)
ax.text(1, 6.3, "escape depth (6 cells)", fontsize=9, color="0.3"); ax.text(0.98, 0.04, "wall band", fontsize=9, color="0.3", transform=ax.transAxes, ha="right")
ax.set_xlabel("Time (t.u.)"); ax.set_ylabel("Distance from the wall (cells)"); ax.set_title("An escape (tracked wall $+\\frac{1}{2}$, sparse ellipse)", fontsize=11)
for a, s in zip(fig.axes, "abcdef"): a.text(-0.06 if s in "abe" else -0.13, 1.06 if s == "a" else 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
for e in ("pdf", "png"): fig.savefig(f"fig_si_measure_alternatives.{e}", dpi=170)
