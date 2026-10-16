#!/usr/bin/env python3
"""Tilt of the wall layer vs crowding (user, 2026-10-15): from tilt_crowding_runs.npz (circles, lambda 2, k 2-10, 0.8 <= R/l_a <= 2.2, 10 seeds; tilt of the
coherent layer = median |frame mean tau| over frames with C > 0.9).
(a) tilt vs k, one line per cavity size band R/l_a (seed means +- s.e.); (b) tilt vs wall-defect spacing 2 pi R/(2k-2) in cells (collapse test: does the
absolute spacing set the tilt?); (c) heatmap of the tilt over (k, R/l_a); (d) heatmap over (rho = k/pi R, R/l_a).  -> fig_tilt_k.png (+pdf)"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
D = np.load(os.path.join(H, "tilt_crowding_runs.npz")); k, R, rho, T = D["k"], D["R"], D["rho"], D["tilt_coh"]; Rl = R / LA
pts = {}
for i in range(len(k)):
    if np.isfinite(T[i]): pts.setdefault((int(k[i]), round(R[i], 2)), []).append(T[i])
P = np.array([(kk, RR, np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0, len(v)) for (kk, RR), v in pts.items()])
pk, pR, pT, pE = P[:, 0], P[:, 1], P[:, 2], P[:, 3]; pRl = pR / LA; prho = pk / (np.pi * pR); sp = 2 * np.pi * pR / (2 * pk - 2)
fig, ax = plt.subplots(2, 2, figsize=(11, 8.6))
bands = [(0.8, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.25)]; cols = plt.cm.plasma(np.linspace(0.05, 0.8, len(bands)))
for (lo, hi), c in zip(bands, cols):
    m = (pRl >= lo) & (pRl < hi)
    for kk in sorted(set(pk[m])):
        pass
    ks = np.array(sorted(set(pk[m]))); mt = [np.average(pT[m & (pk == kk)]) for kk in ks]; me = [np.sqrt(np.sum(pE[m & (pk == kk)] ** 2)) / max((m & (pk == kk)).sum(), 1) for kk in ks]
    ax[0, 0].errorbar(ks, mt, me, fmt="o-", color=c, capsize=2, ms=4, label=f"$R/\\ell_a$ {lo}–{hi}")
    ax[0, 1].errorbar(sp[m], pT[m], pE[m], fmt="o", color=c, ms=4, capsize=1.5, label=f"$R/\\ell_a$ {lo}–{hi}")
for a_ in (ax[0, 0], ax[0, 1]): a_.axhline(60, color="k", ls="--", lw=0.8); a_.text(a_.get_xlim()[1] if False else 0.99, 60.8, "arccos(1/λ), λ = 2", transform=a_.get_yaxis_transform(), ha="right", fontsize=8)
ax[0, 0].set_xlabel("imposed winding k (wall defects 2k − 2)"); ax[0, 0].set_ylabel("tilt of the coherent wall layer (deg)"); ax[0, 0].set_title("(a) Tilt vs winding, by cavity size", fontsize=10)
ax[0, 0].legend(fontsize=8, frameon=False, loc="lower left")
ax[0, 1].set_xlabel("wall-defect spacing 2πR/(2k − 2) (cells)"); ax[0, 1].set_ylabel("tilt of the coherent wall layer (deg)"); ax[0, 1].set_title("(b) Tilt vs spacing (collapse test)", fontsize=10)
ax[0, 1].set_xscale("log"); ax[0, 1].set_xticks([7, 10, 15, 20, 30, 40]); ax[0, 1].set_xticklabels(["7", "10", "15", "20", "30", "40"]); ax[0, 1].minorticks_off(); ax[0, 1].legend(fontsize=8, frameon=False, loc="lower right")
def heat(a_, x, y, xl, yl, title, xb, yb):
    G = np.full((len(xb) - 1, len(yb) - 1), np.nan)
    for i in range(len(xb) - 1):
        for j in range(len(yb) - 1):
            m = (x >= xb[i]) & (x < xb[i + 1]) & (y >= yb[j]) & (y < yb[j + 1])
            if m.any(): G[i, j] = np.mean(pT[m])
    im = a_.pcolormesh(xb, yb, G.T, cmap="magma", vmin=35, vmax=70, shading="flat"); a_.set_xlabel(xl); a_.set_ylabel(yl); a_.set_title(title, fontsize=10)
    return im
im = heat(ax[1, 0], pk, pRl, "imposed winding k", "cavity radius R/ℓa", "(c) Tilt over (k, R/ℓa)", np.arange(1.5, 11), np.linspace(0.8, 2.25, 9))
heat(ax[1, 1], prho, pRl, "wall density ρ = k/πR", "cavity radius R/ℓa", "(d) Tilt over (ρ, R/ℓa)", np.array([0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.11, 0.13, 0.16]), np.linspace(0.8, 2.25, 9))
cax = fig.add_axes([0.88, 0.11, 0.015, 0.33]); cb = fig.colorbar(im, cax=cax); cb.set_label("tilt of the coherent wall layer (deg)")
for a_ in ax.flat: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
fig.subplots_adjust(hspace=0.35, wspace=0.28, right=0.85); fig.savefig(os.path.join(H, "fig_tilt_k.png"), dpi=200)
# collapse statistics: residual scatter of tilt vs spacing vs vs rho vs k
for name, x in (("spacing", np.log(sp)), ("rho", prho), ("k", pk), ("R/la", pRl)):
    c = np.polyfit(x, pT, 2); r = pT - np.polyval(c, x); print(f"quadratic fit of tilt on {name:8s}: residual rms {r.std():.2f} deg (raw sd {pT.std():.2f})")
