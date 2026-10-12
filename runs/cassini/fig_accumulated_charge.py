#!/usr/bin/env python3
"""
Figure: accumulated topological charge Q(c) inside Cassini level c, circles vs ovals (k = 6, zeta = 0.07).
c = sqrt(d1 d2) (circle: c = r); c/b = 1 is the wall, Q(1) = k. A layer + bulk structure shows as a plateau
at Q = +1 just inside the wall (circles); its absence as a smooth rise (ovals).
Curves: mean over seeds, band = min..max over seeds (when > 1 seed). Writes fig_accumulated_charge.{png,pdf}.
Edit SETS below to add runs (e.g. the half-length-34 ovals of P3 when done).
"""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lobe_profiles import profile  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LA = 22.21
C = os.path.join(HERE, "..", "cavity", "server_k06")
P = HERE
# (label, runs, colour, linestyle)
SETS = [
    ("circle, R = 1.43 ℓa", sorted(glob.glob(f"{C}/rho0.06/k06/s0?")), "#86b6ef", "-"),
    ("circle, R = 2.15 ℓa", sorted(glob.glob(f"{C}/rho0.04/k06/s0?")), "#2a78d6", "-"),
    ("circle, R = 2.87 ℓa", sorted(glob.glob(f"{C}/rho0.03/k06/s0?")), "#104281", "-"),
    ("oval b/a 1.6, half-length 2.16 ℓa", [f"{P}/P1/H48_ba1.6_k06"], "#f4a27f", "-"),
    ("oval b/a 1.1, half-length 2.16 ℓa", [f"{P}/P1/H48_ba1.1_k06"], "#eb6834", "-"),
    ("oval b/a 1.02, half-length 2.16 ℓa", [f"{P}/P2/H48_ba1.02_k06"], "#c24d1c", "-"),
    ("oval b/a 1.1, half-length 2.88 ℓa", [f"{P}/P1/H64_ba1.1_k06"], "#8f3413", "--"),
]
small = sorted(glob.glob(f"{P}/P3/H34_ba*_k06"), key=lambda d: -float(d.split("_ba")[1].split("_")[0]))
purples = {"2.5": "#c7b8ec", "2.0": "#a38bdd", "1.6": "#7c5cc4", "1.3": "#4b2f8f"}
for d in small:
    runs = [s for s in sorted(glob.glob(f"{d}/s0?")) if os.path.exists(os.path.join(s, "DONE"))]
    if runs:
        ba = d.split("_ba")[1].split("_")[0]
        SETS.append((f"oval b/a {ba}, half-length 1.53 ℓa ({len(runs)} seeds)", runs, purples.get(ba, "#7c5cc4"), "-."))


def main():
    plt.rcParams.update({"font.size": 10, "axes.labelcolor": "#555a62", "xtick.color": "#555a62",
                         "ytick.color": "#555a62", "axes.edgecolor": "#555a62"})
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for lab, runs, col, ls in SETS:
        if not runs:
            continue
        Q = []
        for r in runs:
            p = profile(r, 80, 0.2)
            x = np.r_[0, p["edges"][1:]]
            Q.append(np.r_[0, p["Q"]])
        Q = np.array(Q)
        ax.plot(x, Q.mean(0), color=col, lw=2.0, ls=ls, label=lab + ("" if len(runs) == 1 or "seeds" in lab else f" ({len(runs)} seeds)"))
        if len(runs) > 1:
            ax.fill_between(x, Q.min(0), Q.max(0), color=col, alpha=0.18, lw=0)
    ax.axhline(1, color="#555a62", lw=0.8, ls=":")
    ax.text(0.02, 1.12, "Q = +1", color="#555a62", fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(-1.5, 6.3)
    ax.set_xlabel("c / b   (Cassini level: 0 = centre / foci, 1 = wall;  circle: r / R)")
    ax.set_ylabel("charge inside level c")
    ax.set_title("Accumulated topological charge, k = 6, ζ = 0.07", fontsize=11, color="#15171b")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e9ebee", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    out = os.path.join(HERE, "fig_accumulated_charge")
    fig.savefig(out + ".png", dpi=200); fig.savefig(out + ".pdf")
    print("written", out + ".png/.pdf")


if __name__ == "__main__":
    main()
