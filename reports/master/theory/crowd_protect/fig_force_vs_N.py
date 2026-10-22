#!/usr/bin/env python3
"""Elastic force on one wall +1/2 vs number of wall defects (force_vs_N.npz; analytic one constant K 44, circle R 31.83, P 200, central +1).
(a) restoring force F(delta) on one defect pulled inward from the equilibrium ring, others fixed; dotted: Coulomb pull of the hole and its image, (pi K/2)(1/delta + 1/(2 d* + delta)).
(b) equilibrium wall distance d* of the ring (force zero) vs rho = N/P, and measured mean resident depth in the sweep circles (minus the 0.61-cell tracker/wall offset of the depth law).
(c) peak restoring force within 8 cells and mean force over 0-6 cells vs rho; dashed: pi^2 K rho/2.
(d) measured escape probability per wall density level (shape_runs.csv, medians; circles and aspect >= 1.4) against the peak restoring force at that rho."""
import csv, os, numpy as np, matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); Z = np.load(H + "/force_vs_N.npz"); NS, DEL = Z["NS"], Z["DEL"]; P = 200.0; K = 44.0
rho = NS / P; ds = Z["dstar"]; cm = plt.get_cmap("Oranges")
fmax = np.array([Z[f"F_N{N}"][DEL <= 8].max() for N in NS]); fmean = np.array([np.interp(6, DEL, Z[f"E_N{N}"]) / 6 for N in NS])
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6)); fig.subplots_adjust(wspace=0.34, left=0.05, right=0.99, bottom=0.15, top=0.9)
for N in (8, 12, 18, 24, 28):
    c = cm(0.35 + 0.6 * (N - 8) / 20); ax[0].plot(DEL, Z[f"F_N{N}"], color=c, lw=1.8, label=fr"$N={N}$ ($\rho={N / P:.2f}$)")
    d0 = ds[list(NS).index(N)]; dd = DEL[DEL > 3.5]; ax[0].plot(dd, np.pi * K / 2 * (1 / dd + 1 / (2 * d0 + dd)), ":", color=c, lw=1.2)
ax[0].set_xlabel(r"Drag distance $\delta$ (cells)"); ax[0].set_ylabel(r"Restoring force $F$ ($k_BT$/cell)"); ax[0].set_ylim(0, 45); ax[0].legend(loc="upper right", fontsize=9.5, framealpha=1)
ax[0].text(0.97, 0.03, "Dotted: hole + image, Coulomb", transform=ax[0].transAxes, ha="right", fontsize=9.5)
meas = Z["meas"]; ax[1].plot(rho, ds, "-", color="k", lw=1.8, label="Elastic ring equilibrium $d^*$")
ax[1].plot(meas[:, 0], meas[:, 1] - 0.61, "o", color="#e6550d", ms=7, label="Measured, circles (minus 0.61)")
ax[1].set_xlabel(r"Wall density $\rho=N/P$"); ax[1].set_ylabel("Wall distance of the layer (cells)"); ax[1].legend(loc="upper right", fontsize=9.5)
ax[2].plot(rho, fmax, "-o", ms=4, color="#a63603", label="Peak force (within 8 cells)"); ax[2].plot(rho, fmean, "-s", ms=4, color="#fd8d3c", label="Mean force, 0–6 cells")
ax[2].plot(rho, np.pi ** 2 * K * rho / 2, "--", color="0.5", label=r"Line charge $\pi^2K\rho/2$"); ax[2].set_xlabel(r"Wall density $\rho=N/P$"); ax[2].set_ylabel(r"Restoring force ($k_BT$/cell)")
ax[2].legend(loc="upper left", fontsize=9.5)
R = list(csv.DictReader(open(H + "/../shapes/shape_runs.csv"))); r_ = np.array([float(r["rho"]) for r in R]); a_ = np.array([float(r["aspect"]) for r in R]); e_ = np.array([float(r["esc"]) for r in R])
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); lev = np.argmin(np.abs(r_[:, None] - LV), 1)
for nm, sel, c, mk in (("Circles", a_ < 1.05, "#3b6fb6", "o"), ("Aspect ≥ 1.4", a_ >= 1.4, "#a63603", "s")):
    y = [np.median(np.maximum(e_[(lev == i) & sel], 1e-5)) for i in range(5)]; x = np.interp(LV, rho, fmax)
    ax[3].plot(x, y, "-" + mk, color=c, ms=6, label=nm)
    for xi, yi, r0 in zip(x, y, LV): ax[3].annotate(f"{r0:.2f}", (xi, yi), textcoords="offset points", xytext=(5, 4), fontsize=8.5, color=c)
ax[3].set_yscale("log"); ax[3].set_xlim(13, 36); ax[3].set_xlabel(r"Peak elastic restoring force ($k_BT$/cell)"); ax[3].set_ylabel(r"Escape probability $p_\mathrm{esc}$"); ax[3].legend(loc="lower left", fontsize=9.5)
ax[3].text(0.97, 0.95, r"Labels: $\rho$", transform=ax[3].transAxes, ha="right", va="top", fontsize=9.5)
for a, s in zip(ax, "abcd"): a.text(-0.13, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
ax[0].set_title("Pull on one leaving defect", fontsize=12); ax[1].set_title("Layer depth", fontsize=12); ax[2].set_title("Restoring force vs density", fontsize=12); ax[3].set_title("Escape vs restoring force", fontsize=12)
for e in ("pdf", "png"): fig.savefig(H + f"/fig_force_vs_N.{e}", dpi=200)
