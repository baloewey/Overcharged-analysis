#!/usr/bin/env python3
"""
Plot the active length l_a vs activity for the L=100 and L=200 periodic boxes, with error bars
and the power-law fit l_a = A act^-alpha (fitted to L=100).

Output: active_length_vs_activity.png / .pdf (next to this script).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))

# Reference palette (dataviz skill, light mode): categorical slots 1-2, surface and ink tokens
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de"


def load(path):
    rows = [r for r in csv.DictReader(open(path)) if r["estimator"] == "all"]
    a = np.array([float(r["act"]) for r in rows])
    l = np.array([float(r["l_active"]) for r in rows])
    e = np.array([float(r["l_active_err"]) for r in rows])
    return a, l, e


a1, l1, e1 = load(os.path.join(HERE, "L100", "active_length.csv"))
a2, l2, e2 = load(os.path.join(HERE, "L200", "active_length.csv"))

f = lambda x, A, al: A * x ** (-al)
p, c = curve_fit(f, a1, l1, p0=[5, 0.5], sigma=e1, absolute_sigma=True)
xs = np.linspace(0.047, 0.093, 200)

plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, ax = plt.subplots(figsize=(6.4, 4.4), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

ax.plot(xs, f(xs, *p), ls="--", lw=1.5, color=INK2, zorder=1,
        label=rf"fit (L=100): $\ell_a \propto \zeta^{{-{p[1]:.2f}\pm{np.sqrt(c[1,1]):.2f}}}$")
ax.errorbar(a1, l1, yerr=e1, fmt="o", ms=7, color=BLUE, mec=SURFACE, mew=1.5,
            ecolor=BLUE, elinewidth=1.5, capsize=3, zorder=3, label="L = 100 box")
ax.errorbar(a2, l2, yerr=e2, fmt="s", ms=7, color=ORANGE, mec=SURFACE, mew=1.5,
            ecolor=ORANGE, elinewidth=1.5, capsize=3, zorder=4, label="L = 200 box")

# selective direct label: the cavity reference value
i = int(np.argmin(abs(a2 - 0.07)))
ax.annotate(rf"$\ell_a(0.07) = {l2[i]:.1f} \pm {e2[i]:.1f}$", xy=(a2[i], l2[i]),
            xytext=(a2[i] + 0.004, l2[i] + 2.3), color=INK, fontsize=10,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))

ax.set_xlabel(r"activity $\zeta$ (act)")
ax.set_ylabel(r"active length $\ell_a = \sqrt{A/\langle N\rangle}$ (cells)")
ax.set_title("Active length vs activity, periodic box", loc="left", color=INK, fontsize=12)
ax.set_xlim(0.046, 0.094)
leg = ax.legend(frameon=False, loc="upper right")
for t in leg.get_texts():
    t.set_color(INK)
fig.text(0.01, 0.005, r"Error bars: $\pm1\sigma$ (time-correlation corrected, one seed per point)."
         r"  $A=(L-2)^2$.", fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"active_length_vs_activity.{ext}"), dpi=200, facecolor=SURFACE)
print("saved active_length_vs_activity.png/.pdf; fit A=%.2f alpha=%.3f" % (p[0], p[1]))
