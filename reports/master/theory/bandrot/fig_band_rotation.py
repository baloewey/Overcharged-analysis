#!/usr/bin/env python3
"""Profiles of the flow-induced director rotation (band_rotation.npz): (a) all azimuths, selected runs; (b) mid-gap vs at-defect for two runs.
Dashed: the band estimate's sign and size, chi gamma1 s d^2/(4K) at the layer depth. -> fig_band_rotation.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z = np.load("band_rotation.npz"); uc = Z["uc"]
sel = [("swirl_crowd/k04_s1", "k 4, s 0.043, d 2.2"), ("swirl_shear/k06_A0.010", "k 6, s 0.043, d 1.9"), ("swirl_crowd/A03_k06_s1", "k 6, s 0.077, d 1.7"),
       ("swirl_crowd/k08_s2", "k 8, s 0.043, d 1.5"), ("swirl_shear/k13_A0.100", "k 13, s 0.147, d 1.3"), ("swirl_shear/k06_A0.100", "k 6, s 0.147, d 1.7")]
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for (r, lab), c in zip(sel, [f"C{i}" for i in range(6)]):
    ax[0].plot(uc, Z[f"{r}|all"], "-o", ms=4, color=c, label=lab)
for (r, lab), c in zip(sel[1:3], ["C1", "C2"]):
    ax[1].plot(uc, Z[f"{r}|mid"], "-o", ms=4, color=c, label=f"{lab}: between defects"); ax[1].plot(uc, Z[f"{r}|at"], "--s", ms=4, color=c, label=f"{lab}: at the defects")
for a in ax:
    a.axhline(0, color="k", lw=0.8); a.axvspan(0, 1.9, color="0.9", zorder=0); a.set_xlabel("Distance From the Wall $u$ (cells)"); a.set_xlim(0, 10)
    a.set_ylabel(r"Flow-Induced Rotation $\delta\theta$ (deg)"); a.text(0.15, 5.6, "band (wall to defects)", fontsize=10, color="0.35")
    a.annotate("band estimate +1 to +4", xy=(1.0, 2.0), xytext=(3.6, 2.6), fontsize=10, arrowprops=dict(arrowstyle="->", color="0.3"), color="0.3")
ax[0].set_title("Rotation Against Distance From the Wall", pad=10); ax[0].legend(fontsize=9, loc="lower left"); ax[0].set_ylim(-38, 8)
ax[1].set_title("Between and At the Wall Defects", pad=10); ax[1].legend(fontsize=9, loc="lower left"); ax[1].set_ylim(-20, 8)
for a, l in zip(ax, "ab"): a.text(-0.13, 1.04, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.05, r"Passive MPCD circles under an imposed vortex. $\delta\theta$ = measured director minus the exact static one-constant field of the same defects (anchoring $\theta=k\phi$), cells more" + "\n" + r"than 1.5 cells from any core; positive = rotated further in the sense of the tilt. $s$: wall shear rate; $d$: depth of the wall defects. Grey: band between wall and defects.", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("fig_band_rotation.pdf", bbox_inches="tight"); fig.savefig("fig_band_rotation.png", dpi=150, bbox_inches="tight")
