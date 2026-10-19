#!/usr/bin/env python3
"""SI Fig. S_mirror / long-draft Fig. 10 (2026-10-19): +k (wall +1/2) vs -k (wall -1/2) in the same 675 cavities. Medians over k, rho, 3 seeds per (family, aspect, R_eq).
(a) wall defects lost vs the circle at the same R_eq; (b) bulk charge Q_b; (c) escape probability p_esc per wall defect per 20 t.u.; (d) |circulation|.
Sources: shape_runs.csv (+k), shape_runs_mirror.csv (-k). Run via reports/figstyle/run.sh. Output: theory_report/figs/fig_mirror_shapes.{pdf,png}."""
import csv, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = "/home/baloewe/Physics/Overcharged/reports/theory_report/figs/fig_mirror_shapes"
def load(f):
    G = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(f)):
        for c in ("vac", "qbulk", "esc", "circ"): G[(r["family"], round(float(r["aspect"]), 2), r["Req_la"][:3])][c].append(float(r[c]))
    return G
P, M = load("shape_runs.csv"), load("shape_runs_mirror.csv")
cols = {"1.0": "#9ecae1", "1.5": "#4a90c2", "2.0": "#08306b"}
fig, axs = plt.subplots(1, 4, figsize=(13, 3.9)); fig.subplots_adjust(wspace=0.38, left=0.05, right=0.99, bottom=0.15, top=0.78)
labs = ["Wall defects lost, $-\\Delta N_\\mathrm{wall}$", "Bulk charge $Q_b$", "Escape probability $p_\\mathrm{esc}$", "Circulation $|\\Gamma|$"]
for j, (c, yl) in enumerate(zip(("vac", "qbulk", "esc", "circ"), labs)):
    ax = axs[j]
    for G, ls, mf in ((P, "-", None), (M, ":", "w")):
        for fam, mk in (("ellipse", "o"), ("oval", "s")):
            for req in cols:
                ks = sorted((k for k in G if k[2] == req and k[0] in ("circle", fam)), key=lambda k: k[1]); c0 = np.median(G[("circle", 1.0, req)][c]) if c == "vac" else 0
                ax.plot([k[1] for k in ks], [np.median(G[k][c]) - c0 for k in ks], ls=ls, marker=mk, ms=4, lw=1.1, color=cols[req], mfc=mf or cols[req], alpha=1 if fam == "ellipse" else 0.6)
    ax.set_xlabel("Aspect ratio"); ax.set_ylabel(yl); ax.text(-0.18, 1.04, "(" + "abcd"[j] + ")", transform=ax.transAxes, fontweight="bold")
    if c == "esc": ax.set_yscale("log")
    if c == "vac": ax.axhline(0, color="0.8", lw=0.6, zorder=0); ax.set_ylim(-1.3, 7.0)
h = [Line2D([], [], color=cols[r], lw=1.5, label=fr"$R_\mathrm{{eq}}={r}\,\ell_a$") for r in cols]
h += [Line2D([], [], color="0.25", ls="-", marker="o", label="$+k$ (wall $+1/2$)"), Line2D([], [], color="0.25", ls=":", marker="o", mfc="w", label="$-k$ (wall $-1/2$)"),
      Line2D([], [], color="0.25", ls="", marker="o", label="Ellipse"), Line2D([], [], color="0.25", ls="", marker="s", alpha=0.6, label="Cassini oval")]
fig.legend(handles=h, loc="upper center", fontsize=9, ncol=7, columnspacing=1.2, handlelength=2.0, frameon=False)
fig.savefig(OUT + ".pdf"); fig.savefig(OUT + ".png", dpi=200); print("wrote", OUT)
