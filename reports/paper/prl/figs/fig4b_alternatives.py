#!/usr/bin/env python3
"""Alternatives for PRL Fig. 4(b) (2026-10-22, user: draw alternatives, decide later). Data: shapes/shape_runs.csv (675 runs: circles, ellipses, Cassini ovals;
R_eq 1-2 l_a; five wall densities rho = 2k/P; 3 seeds). Each panel is one candidate, drawn at the printed single-panel size.
(A) escape probability vs tilt coherence C, per rho (binned medians, log)   (B) vacant fraction vs C, per rho (binned medians)
(C) vacant fraction vs aspect, per rho (the current 4b)                       (D) escape probability vs aspect, per rho (log)
(E) map of vacant fraction over aspect x rho (ellipses + ovals, medians), C = 0.85 contour   (F) all runs: vacant fraction vs C, coloured by rho."""
import csv, numpy as np
from collections import defaultdict
from prlstyle import plt, PAGE, DARK, label, save, ROOT
from matplotlib.lines import Line2D
R = list(csv.DictReader(open(ROOT + "/reports/master/theory/shapes/shape_runs.csv")))
RHO = (0.04, 0.06, 0.09, 0.12, 0.14); RCOL = ("#fdae6b", "#f98c3a", "#e6550d", "#a63603", "#4d1a02")
rho = np.array([float(r["rho"]) for r in R]); lev = np.argmin(np.abs(rho[:, None] - np.array(RHO)), 1)
C = np.array([float(r["coh"]) for r in R]); f = np.array([float(r["vac"]) / (2 * int(r["k"]) - 2) for r in R])
esc = np.array([float(r["esc"]) for r in R]); asp = np.round(np.array([float(r["aspect"]) for r in R]), 2); fam = np.array([r["family"] for r in R])
CB = np.array([0.3, 0.45, 0.55, 0.65, 0.75, 0.85, 0.92, 1.0])
def binned(ax, y, logy=False):
    for i in range(5):
        xs, ys = [], []
        for a, b in zip(CB[:-1], CB[1:]):
            m = (lev == i) & (C >= a) & (C < b)
            if m.sum() >= 3: xs.append(np.median(C[m])); ys.append(np.median(y[m]))
        ax.plot(xs, ys, "-o", ms=2.8, lw=1.0, color=RCOL[i])
    ax.axvline(0.85, color="0.75", lw=0.6, ls=":", zorder=0); ax.set_xlabel("Tilt coherence $C$"); ax.set_xlim(0.3, 1.0)
    if logy: ax.set_yscale("log")
def by_aspect(ax, y):
    for i in range(5):
        for fm, ls, mk in (("ellipse", "-", "o"), ("oval", "--", "s")):
            ks = sorted(set(asp[(fam == fm) & (lev == i)]))
            xs = [1.0] + ks; ys = [np.median(y[(fam == "circle") & (lev == i)])] + [np.median(y[(fam == fm) & (lev == i) & (asp == a)]) for a in ks]
            ax.plot(xs, ys, ls=ls, marker=mk, ms=2.4, lw=0.9, color=RCOL[i])
    ax.set_xlabel("Aspect ratio"); ax.set_xticks([1.0, 1.5, 2.0, 2.5])
def rho_legend(ax, loc, ncol=3, fam_too=False):
    h = [Line2D([], [], color=RCOL[i], lw=1.4, label=f"{r:.2f}") for i, r in enumerate(RHO)]
    if fam_too: h += [Line2D([], [], color=DARK, ls="-", marker="o", ms=2.4, label="Ellipse"), Line2D([], [], color=DARK, ls="--", marker="s", ms=2.4, label="Oval")]
    ax.legend(handles=h, title=r"$\rho=2k/P$", title_fontsize=6, loc=loc, fontsize=6, ncol=ncol, handlelength=1.2, columnspacing=0.8, borderaxespad=0.1)
fig, axs = plt.subplots(2, 3, figsize=(PAGE, 4.6)); fig.subplots_adjust(wspace=0.62, hspace=0.62, left=0.08, right=0.97, top=0.93, bottom=0.1)
pe = np.where(esc > 0, esc, 1e-4)
ax = axs[0, 0]; binned(ax, pe, True); ax.set_ylabel(r"Escape probability $p_\mathrm{esc}$"); ax.set_ylim(5e-5, 2.0); rho_legend(ax, "upper left")
ax.set_title("A: escape vs coherence", fontsize=7, loc="left", pad=10)
ax = axs[0, 1]; binned(ax, f); ax.set_ylabel(r"Vacant fraction $f_\mathrm{vac}$"); ax.set_ylim(-0.02, 0.75); rho_legend(ax, "upper right")
ax.set_title("B: vacancies vs coherence", fontsize=7, loc="left", pad=10)
ax = axs[0, 2]; by_aspect(ax, f); ax.set_ylabel(r"Vacant fraction $f_\mathrm{vac}$"); ax.set_ylim(-0.02, 0.8); rho_legend(ax, "upper left")
ax.set_title("C: vacancies vs aspect (current)", fontsize=7, loc="left", pad=10)
ax = axs[1, 0]; by_aspect(ax, pe); ax.set_yscale("log"); ax.set_ylabel(r"Escape probability $p_\mathrm{esc}$"); ax.set_ylim(5e-5, 8.0); rho_legend(ax, "upper left", ncol=4, fam_too=True)
ax.set_title("D: escape vs aspect", fontsize=7, loc="left", pad=10)
ax = axs[1, 1]; AS = sorted(set(asp[fam != "circle"]) | {1.0}); AB = np.array([1.0, 1.1, 1.25, 1.33, 1.45, 1.6, 1.8, 2.2, 2.6])
M = np.full((5, len(AB) - 1), np.nan); MC = np.full_like(M, np.nan)
for i in range(5):
    for j, (a, b) in enumerate(zip(AB[:-1], AB[1:])):
        m = (lev == i) & (asp >= a) & (asp < b)
        if m.sum(): M[i, j] = np.median(f[m]); MC[i, j] = np.median(C[m])
im = ax.pcolormesh(AB, np.arange(6) - 0.5, M, cmap="Oranges", vmin=0, vmax=0.55, shading="flat")
xc = 0.5 * (AB[:-1] + AB[1:]); ax.contour(xc, np.arange(5), MC, levels=[0.85], colors=DARK, linewidths=1.0, linestyles="--")
ax.set_yticks(range(5)); ax.set_yticklabels([f"{r:.2f}" for r in RHO]); ax.set_ylabel(r"Wall density $\rho=2k/P$"); ax.set_xlabel("Aspect ratio"); ax.set_xticks([1.0, 1.5, 2.0, 2.5])
cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.06); cb.set_label(r"$f_\mathrm{vac}$", fontsize=7, labelpad=1); cb.ax.tick_params(labelsize=6)
ax.text(1.42, 3.0, "$C=0.85$\n(dashed)", fontsize=6, color=DARK)
ax.set_title("E: vacancy map", fontsize=7, loc="left", pad=10)
ax = axs[1, 2]
for i in range(5): ax.scatter(C[lev == i], f[lev == i], s=4, lw=0, color=RCOL[i], alpha=0.8, rasterized=True)
ax.axvline(0.85, color="0.75", lw=0.6, ls=":", zorder=0); ax.set_xlabel("Tilt coherence $C$"); ax.set_ylabel(r"Vacant fraction $f_\mathrm{vac}$"); ax.set_xlim(0.25, 1.0); ax.set_ylim(-0.05, 0.8)
rho_legend(ax, "upper right"); ax.set_title("F: all runs", fontsize=7, loc="left", pad=10)
for ax, s in zip(axs.flat, "abcdef"): label(ax, s, y=1.12)
save(fig, "fig4b_alternatives")
