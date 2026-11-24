#!/usr/bin/env python3
"""Collapse test (user, 2026-11-24): does the gap tilt - 2 x band depend only on how deep the wall defects' near field reaches compared with the flow-aligned band?
Data: band_windows.npy (MPCD chi 0.5 R 31.8, chi 1 R 21.6; static-subtracted; coherent frames, exactly 2k-2 wall defects).
 (a) tilt and 2 x band read just outside the near field ([edge, edge+1), edge = d + a/2pi) vs k.
 (b) gap with the edge window vs edge depth in cells.
 (c) same vs edge / l_b (l_b = depth where the band angle profile falls to half its maximum; chi 1: the band does not halve within the measured 12 cells,
     so l_b > 12 and the chi-1 points are upper limits in x, drawn with left arrows).
-> fig_collapse.pdf/.png"""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); D = np.load(os.path.join(H, "band_windows.npy"), allow_pickle=True).item()
C = {"MPCD, chi = 0.5": "#c0392b", "MPCD, chi = 1": "#7d3c98"}; LB1 = 12.0
fig, axs = plt.subplots(1, 3, figsize=(15, 5.0))
for lab, rows in D.items():
    c = C[lab]; L = lab.replace("chi", "χ"); k = np.array([r["k"] for r in rows]); tau = np.array([r["tau"] for r in rows]); ec = np.array([r["ec"] for r in rows])
    edge = np.array([r["edge"] for r in rows]); lb = np.array([r["lb"] for r in rows]); gap = tau - ec
    axs[0].plot(k, tau, "o-", color=c, lw=2, ms=7, label=f"{L}: layer tilt"); axs[0].plot(k, ec, "s--", color=c, mfc="white", mew=2, lw=2, ms=7, label=f"{L}: 2 × band just outside near field")
    axs[1].plot(edge, gap, "o-", color=c, lw=2, ms=7, label=L)
    if np.all(np.isnan(lb)):
        x = edge / LB1; axs[2].plot(x, gap, "o", color=c, ms=8, label=f"{L} (l_b > {LB1:g}: upper limits)")
        for xi, gi in zip(x, gap): axs[2].annotate("", (xi - 0.06, gi), (xi, gi), arrowprops=dict(arrowstyle="->", color=c, lw=1.5))
    else:
        axs[2].plot(edge / lb, gap, "o-", color=c, lw=2, ms=7, label=L)
axs[0].axhline(60, color="0.6", ls="--", lw=1.2); axs[0].set_xlabel("Winding k"); axs[0].set_ylabel("Angle (deg)"); axs[0].set_ylim(0, 75)
axs[0].set_title("(a) Tilt vs band just outside the near field", fontsize=12.5); axs[0].legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=1)
for ax in axs[1:]: ax.axhline(0, color="0.6", lw=1); ax.set_ylabel("Gap: tilt − 2 × band (deg)"); ax.legend(fontsize=10, frameon=False)
axs[1].set_xlabel("Near-field edge d + a/2π (cells)"); axs[1].set_title("(b) Gap vs depth of the near field", fontsize=12.5)
axs[2].set_xlabel(r"Near-field edge / band half-depth $\ell_b$"); axs[2].set_title("(c) Gap vs near field relative to the band", fontsize=12.5)
fig.text(0.5, -0.22, "Band: static-subtracted ring average of the director angle from the flow, in the 1-cell window just outside each layer's near field. "
         r"MPCD circles at $R/\ell_a$ = 1.43 (χ = 0.5: R 31.8; χ = 1: R 21.6).", ha="center", fontsize=10.5, color="0.3")
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_collapse.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_collapse.png"), dpi=180, bbox_inches="tight")
