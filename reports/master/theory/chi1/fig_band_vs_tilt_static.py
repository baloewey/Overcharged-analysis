#!/usr/bin/env python3
"""Version 2 of fig_band_vs_tilt with the STATIC-SUBTRACTED band (user, 2026-11-24): layer tilt vs 2 alpha_corr = 2 alpha_meas - 2 alpha_static + tau_static
(static_band.py; depth window 2-3 cells, coherent frames with exactly 2k-2 wall defects), with the plain ring average 2 alpha_meas in light dotted for
comparison. Six setups as v1; values are means over seeds of per-run medians (static_band_mpcd.npy, static_band_cont.npy). Panels without data say so.
-> fig_band_vs_tilt_static.pdf/.png"""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
D = {}
for w in ("mpcd", "cont"):
    f = os.path.join(H, f"static_band_{w}.npy")
    if os.path.exists(f): D.update(np.load(f, allow_pickle=True).item())
order = ["Continuum laminar, chi = 1", "Continuum laminar, chi = 0.5", "Continuum turbulent, chi = 1", "Continuum turbulent, chi = 0.5", "MPCD, chi = 1", "MPCD, chi = 0.5"]
fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.2), sharey=True)
for ax, lab, let in zip(axs.ravel(), order, "abcdef"):
    ax.axhline(60, color="0.6", ls="--", lw=1.2, label=r"$2\theta_L = \arccos(1/\lambda)$")
    rows = sorted(D.get(lab, []))
    if rows:
        k = [r[0] for r in rows]; tau = [r[1] for r in rows]; am = [r[2] for r in rows]; corr = [r[3] for r in rows]
        ax.plot(k, tau, "o-", color="#c0392b", ms=8, lw=2, label="Layer tilt")
        ax.plot(k, corr, "s--", color="#1f5fa8", mfc="white", mew=2, ms=8, lw=2, label=r"$2\times$ band, static-subtracted")
        ax.plot(k, am, ":", color="#7fa7d6", lw=1.8, label=r"$2\times$ band, plain ring average")
        ax.set_xlim(min(k) - 1, max(k) + 1)
    else:
        ax.text(0.5, 0.5, "pending (continuum runs\nwith field snapshots)", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="0.5")
    ax.set_title(f"({let}) {lab.replace('chi', chr(967))}", fontsize=13); ax.set_xlabel("Winding k"); ax.set_ylim(25, 75)
for ax in axs[:, 0]: ax.set_ylabel("Angle (deg)")
axs[1, 2].legend(loc="lower left", fontsize=10, frameon=False)
fig.text(0.5, -0.01, "Band: director angle from the flow direction over depths 2–3 cells, with the static texture of the same defects subtracted (fig_estimator_diagram). "
         r"Circles at $R/\ell_a$ = 1.43; coherent frames.", ha="center", fontsize=10.5, color="0.3")
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_band_vs_tilt_static.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_band_vs_tilt_static.png"), dpi=200, bbox_inches="tight")
