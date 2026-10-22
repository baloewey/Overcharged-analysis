#!/usr/bin/env python3
"""Figure for lock_shift.npz: (a) position shift of the wall layer from the homeotropic points, as N x (angular shift) = arg Z_N, against the
layer tilt, every run; dashed y = x (prediction: flow sets the tilt, anchoring ties the shift to it); dotted: what a shift frozen at the sparse value
would give. (b) The deviation D = arg Z_N - tilt against the wall-defect spacing. -> fig_lock_shift.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z = np.load("lock_shift.npz"); sets = sorted({k.split("|")[0] for k in Z.files})
st = {"active main sweep": ("o", "C0", 0.25, 12), "passive imposed vortex": ("s", "C1", 0.9, 40), "lambda scan k 6": ("^", "C2", 0.9, 40), "lambda scan k 14": ("D", "C3", 0.9, 40)}
fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
for nm in sets:
    A, T, D, sp = (Z[f"{nm}|{q}"] for q in ("A", "T", "D", "sp")); mk, col, al, ms = st.get(nm, ("o", "k", 0.5, 20))
    ax[0].scatter(np.degrees(T), np.degrees(A), marker=mk, color=col, alpha=al, s=ms, label=f"{nm[0].upper() + nm[1:]} ({len(A)} runs)", edgecolors="none" if al < 0.5 else "k", linewidths=0.4)
    ax[1].scatter(sp, np.degrees(D), marker=mk, color=col, alpha=al, s=ms, edgecolors="none" if al < 0.5 else "k", linewidths=0.4)
x = np.array([30, 90]); ax[0].plot(x, x, "k--", lw=1.2, label="Shift follows the tilt (prediction)"); ax[0].axhline(64.5, color="0.4", ls=":", lw=1.2, label="Shift frozen at the sparse value")
ax[0].set_xlabel(r"Layer Tilt $\tau$ (deg)"); ax[0].set_ylabel(r"Position Shift $\times\,(2k-2)$, arg $Z_{2k-2}$ (deg)"); ax[0].set_xlim(30, 90); ax[0].set_ylim(30, 90)
ax[0].set_title("Position Shift Against Tilt", pad=10); ax[0].legend(fontsize=9, loc="lower right")
ax[1].axhline(0, color="k", ls="--", lw=1); ax[1].set_xscale("log"); ax[1].set_xlabel("Wall-Defect Spacing (cells)"); ax[1].set_ylabel(r"Deviation arg $Z_{2k-2}-\tau$ (deg)")
ax[1].set_title("Deviation From the Prediction", pad=10); ax[1].set_ylim(-30, 40)
for a, l in zip(ax, "ab"): a.text(-0.14, 1.04, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.04, r"Each point: one run, frames with exactly $2k-2$ wall $+\frac{1}{2}$ and coherent orientation. Shift measured from the homeotropic points (anchoring $\theta=k\phi$);" + "\n" + r"a shift of the whole layer by $\delta\phi$ along the wall gives arg $Z_{2k-2}=(2k-2)\delta\phi$. Signs set by the circulation.", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("fig_lock_shift.pdf", bbox_inches="tight"); fig.savefig("fig_lock_shift.png", dpi=150, bbox_inches="tight")
