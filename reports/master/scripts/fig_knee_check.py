#!/usr/bin/env python3
"""
Knee-finder check for the master document (replaces the interactive knee_check.png of 2 Oct).

Six main-sweep circles, 3 seeds each. Q(r) from lobe_profiles.profile (200 shells, frames after 20%); knee from
charge_knee.knee with 5-bin smoothing (old) and 9-bin smoothing (current, 0.045 in r/R).
Per run two stacked axes: full Q(r) (top) and a zoom on the outer region (bottom) with raw Q (thin steps),
9-bin smoothed Q (thick), old knee (open circle) and current knee (filled circle).
Writes ../figs/knee_check_v2.png and prints the knees.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "runs"))
sys.path.insert(0, os.path.join(RUNS, "cassini"))
from charge_knee import knee        # noqa: E402
from lobe_profiles import profile   # noqa: E402

LA = 22.21
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
SEEDCOL = ["#2a78d6", "#eb6834", "#1baf7a"]
B = os.path.join(RUNS, "cavity", "server_main", "main_act0.07")
PANELS = [  # (rho, k, tag, short reading)
    ("0.08", 5, "failure", "old knee on a layer step"),
    ("0.09", 5, "failure", "old knee on a layer step"),
    ("0.08", 3, "failure", "old knee on a layer step; tiny cavity"),
    ("0.06", 5, "control", "both agree"),
    ("0.04", 6, "control", "both agree"),
    ("0.12", 18, "open case", "real tilted plateau at Q ≈ 1.4–1.6"),
]


def smooth(Q, win):
    ker = np.ones(win) / win
    return np.convolve(np.r_[np.full(win, Q[0]), Q, np.full(win, Q[-1])], ker, "same")[win:-win]


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7.5); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


fig = plt.figure(figsize=(12.5, 11))
gs = fig.add_gridspec(4, 3, height_ratios=[1, 1.15, 1, 1.15], hspace=0.42, wspace=0.25)
letters = "abcdef"
for n, (rho, k, tag, reading) in enumerate(PANELS):
    r0, c = 2 * (n // 3), n % 3
    axT = fig.add_subplot(gs[r0, c]); axZ = fig.add_subplot(gs[r0 + 1, c])
    R = k / (np.pi * float(rho))
    lines = []
    for s, col in zip((1, 2, 3), SEEDCOL):
        p = profile(os.path.join(B, f"rho{rho}", f"k{k:02d}", f"s0{s}"), 200, 0.2)
        x = p["edges"][1:]; Q = p["Q"]; Qs = smooth(Q, 9)
        k5, k9 = knee(p, win=5), knee(p, win=9)
        axT.plot(x, Q, color=col, lw=1.2)
        axZ.step(x, Q, where="post", color=col, lw=0.7, alpha=0.55)
        axZ.plot(x, Qs, color=col, lw=2.0)
        for ax in (axT, axZ):
            if k5: ax.plot(k5[0], k5[1], "o", ms=8, mfc="none", mec=col, mew=1.6, zorder=5)
            if k9: ax.plot(k9[0], k9[1], "o", ms=7, mfc=col, mec="white", mew=0.8, zorder=6)
        lines.append(f"s{s}: old {'—' if k5 is None else f'{k5[1]:.2f}'}, now {'—' if k9 is None else f'{k9[1]:.2f}'}")
        print(f"rho{rho} k{k} s{s}: 5-bin {k5[:2] if k5 else None}  9-bin {k9[:2] if k9 else None}")
    for ax in (axT, axZ):
        ax.axhline(1, color=INK2, lw=0.8, ls=":"); style(ax)
    axT.set_title(f"({letters[n]}) {tag.upper()}: k = {k}, R = {R:.1f} = {R / LA:.2f} ℓa (ρ = {rho})", loc="left", fontsize=8.5, color=INK)
    axT.set_ylabel("Q(r)", fontsize=8.5, color=INK); axT.set_xlim(0, 1)
    axZ.set_xlim(0.55, 1.0); axZ.set_ylim(min(-0.2, axZ.get_ylim()[0]), min(k, 3.2))
    axZ.set_xlabel("r / R", fontsize=8.5, color=INK); axZ.set_ylabel("Q(r), zoom", fontsize=8.5, color=INK)
    axZ.text(0.02, 0.97, reading + "\n" + "\n".join(lines), transform=axZ.transAxes, va="top", fontsize=6.6, color=INK2,
             bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
from matplotlib.lines import Line2D
h = [Line2D([], [], color=c, lw=2, label=f"seed {i + 1}") for i, c in enumerate(SEEDCOL)]
h += [Line2D([], [], color=INK2, lw=0.7, label="raw Q(r) (zoom: thin steps)"), Line2D([], [], color=INK2, lw=2, label="9-bin smoothed Q (zoom)"),
      Line2D([], [], ls="", marker="o", ms=8, mfc="none", mec=INK, label="old knee, 5-bin smoothing"),
      Line2D([], [], ls="", marker="o", ms=7, mfc=INK, mec="white", label="current knee, 9-bin smoothing"),
      Line2D([], [], color=INK2, ls=":", label="Q = +1")]
fig.legend(handles=h, loc="lower center", ncol=7, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.0))
fig.savefig(os.path.join(HERE, "..", "figs", "knee_check_v2.png"), dpi=170, bbox_inches="tight")
print("written knee_check_v2.png")
