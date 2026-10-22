#!/usr/bin/env python3
"""Switch-on across crowding and shape (switchon2.npz): wall +1/2 count relative to 2k-2, tilt |tau| and coherence C against time from the switch
(passive before 0), 100 t.u. running means. Circles k 4/10/16 (rho 0.04/0.10/0.16) and ellipses aspect 1.5/2 (k 6); solid: initial director along x,
dashed: random."""
import numpy as np, matplotlib.pyplot as plt
Z = np.load("switchon2.npz"); names = sorted({k.rsplit("_", 1)[0] for k in Z.files})
col = {"circle_k04": "#fdae6b", "circle_k10": "#e6550d", "circle_k16": "#4d1a02", "ellipse_e1.5": "#6baed6", "ellipse_e2.0": "#08306b"}
lab = {"circle_k04": r"Circle $k=4$ ($\rho=0.04$)", "circle_k10": r"Circle $k=10$ ($\rho=0.10$)", "circle_k16": r"Circle $k=16$ ($\rho=0.16$)", "ellipse_e1.5": "Ellipse 1.5, $k=6$", "ellipse_e2.0": "Ellipse 2, $k=6$"}
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2)); fig.subplots_adjust(wspace=0.3, left=0.05, right=0.99, bottom=0.25, top=0.9)
def sm(x, n=50):
    y = np.convolve(np.nan_to_num(x, nan=np.nanmean(x)), np.ones(n) / n, mode="same"); y[:n // 2] = np.nan; y[-n // 2:] = np.nan; return y
for nm in names:
    A, B = Z[nm + "_passive"], Z[nm + "_active"]; base = nm.rsplit("_", 1)[0]; ls = "-" if nm.endswith("_x") else "--"; ts = B[0, 0]
    k = int(base.split("_k")[1]) if "_k" in base else 6; N = 2 * k - 2; M = np.vstack([A[A[:, 0] > A[-1, 0] - 1000], B]); t = M[:, 0] - ts
    lb = lab[base] if ls == "-" else None
    ax[0].plot(t, sm(M[:, 1] - N), ls, color=col[base], lw=1.3, label=lb); ax[1].plot(t, sm(np.degrees(np.abs(M[:, 2]))), ls, color=col[base], lw=1.3, label=lb)
    ax[2].plot(t, sm(M[:, 3]), ls, color=col[base], lw=1.3, label=lb)
for a, yl, s, tt in zip(ax, ("Wall $+\\frac{1}{2}$ minus $(2k-2)$", r"Layer tilt $|\tau|$ (deg)", "Tilt coherence $C$"), "abc", ("Wall count", "Tilt", "Coherence")):
    a.axvline(0, color="0.5", ls=":", lw=1); a.set_xlabel("Time from the switch (t.u.)"); a.set_ylabel(yl); a.set_title(tt, fontsize=12); a.text(-0.12, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=5, fontsize=10, frameon=False); ax[1].set_ylim(0, 175); ax[2].set_ylim(0, 1.05); ax[0].text(-950, ax[0].get_ylim()[0] * 0.9, "passive", fontsize=9, color="0.4")
ax[1].text(0.98, 0.97, "Solid: director along $x$; dashed: random", transform=ax[1].transAxes, ha="right", va="top", fontsize=9)
for e in ("pdf", "png"): fig.savefig(f"fig_switchon2.{e}", dpi=180)
