#!/usr/bin/env python3
"""Figure for local_lock.npz: (a) tilt against twice the near-wall director angle just inside the layer (2 alpha_1), per run; dashed y = x
(lock to the local field). (b) The two parts of the tilt against spacing: 2 alpha_1 and the excess tau - 2 alpha_1. (c) Measured depth against the
elastic equilibrium depth of the ring (one constant, images + neighbours + bulk +1). -> fig_local_lock.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z = np.load("local_lock.npz"); st = {"active": ("o", "C0"), "passive": ("s", "C1")}
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
for nm, (mk, c) in st.items():
    g = lambda q: Z[f"{nm}|{q}"]
    ax[0].scatter(2 * g("a1"), g("tau"), marker=mk, color=c, s=35, edgecolors="k", linewidths=0.4, label=f"{nm.capitalize()} ({len(g('tau'))} runs)")
    ax[1].scatter(g("a"), 2 * g("a1"), marker=mk, color=c, s=35, edgecolors="k", linewidths=0.4, label=rf"{nm.capitalize()}, $2\alpha_1$")
    ax[1].scatter(g("a"), g("tau") - 2 * g("a1"), marker=mk, facecolors="none", edgecolors=c, s=35, linewidths=1.2, label=rf"{nm.capitalize()}, $\tau-2\alpha_1$")
    ax[2].scatter(g("d_law"), g("d"), marker=mk, color=c, s=35, edgecolors="k", linewidths=0.4, label=nm.capitalize())
x = np.array([35, 75]); ax[0].plot(x, x, "k--", lw=1.2, label="Tilt = twice the local field"); ax[0].axhline(60, color="0.5", ls=":", lw=1)
ax[0].set_xlabel(r"Twice the Near-Wall Director Angle $2\alpha_1$ (deg)"); ax[0].set_ylabel(r"Layer Tilt $\tau$ (deg)"); ax[0].set_title("Tilt Against the Local Field", pad=10); ax[0].legend(fontsize=9, loc="lower right")
ax[1].set_xscale("log"); ax[1].set_xlabel("Wall-Defect Spacing (cells)"); ax[1].set_ylabel("Angle (deg)"); ax[1].set_title(r"Local Field ($2\alpha_1$) and Excess ($\tau-2\alpha_1$)", pad=10); ax[1].legend(fontsize=9, loc="upper center", ncol=2, columnspacing=0.8); ax[1].set_ylim(-5, 88)
y = np.array([0.8, 2.6]); ax[2].plot(y, y, "k--", lw=1.2, label="Equal"); ax[2].set_xlabel("Elastic Equilibrium Depth (cells)"); ax[2].set_ylabel("Measured Depth (cells)")
ax[2].set_title("Depth of the Wall Defects", pad=10); ax[2].legend(fontsize=10, loc="upper left")
for a, l in zip(ax, "abc"): a.text(-0.14, 1.04, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.05, r"$\alpha_1$: ring-averaged director angle from the local tangent (signed by the circulation) one cell inside the wall defects; $\tau$: coherent-layer tilt. Active: main sweep (seeds 1--3) and" + "\n" + r"switch-on reruns; passive: imposed vortex. Equilibrium depth: one-constant ring of $2k-2$ wall $+\frac{1}{2}$ with a central $+1$ in a disk of radius $R$ (images and neighbours).", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("fig_local_lock.pdf", bbox_inches="tight"); fig.savefig("fig_local_lock.png", dpi=150, bbox_inches="tight")
