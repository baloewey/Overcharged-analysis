#!/usr/bin/env python3
"""Channel crossover figure (channel_crossover.txt; W 40, L 192, validated tracker) with the corrected cavity escape rates for comparison
(crowd_protect/exchange/shape_esc2.csv, aspect >= 1.4 and circles, medians per rho level). (a) escape probability per wall +1/2 per 20 t.u.; (b) clean
escape and entry hazards (channels); (c) wall occupancy N_wall/(2m) and bulk charge per wall length (channels). Shaded: predicted crossover rho_c 0.08-0.1."""
import csv, numpy as np, matplotlib.pyplot as plt
rows = [l.split("|") for l in open("channel_crossover.txt").read().splitlines()[2:]]; m = np.array([int(r[0]) for r in rows]); rho = 2 * m / 192
pe, he, hn, occ, qb = (np.array([float(r[j]) for r in rows]) for j in (3, 4, 5, 6, 7)); ok = m < 18
E2 = {r["run"]: float(r["p_esc"]) for r in csv.DictReader(open("../../crowd_protect/exchange/shape_esc2.csv")) if not r.get("err")}
S = list(csv.DictReader(open("../../shapes/shape_runs.csv"))); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
def cav(lo, hi): return [np.median([E2[r["run"]] for r in S if lo <= float(r["aspect"]) <= hi and np.argmin(abs(float(r["rho"]) - LV)) == i]) for i in range(5)]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6)); fig.subplots_adjust(wspace=0.3, left=0.06, right=0.99, bottom=0.15, top=0.9)
for a in ax: a.axvspan(0.08, 0.10, color="0.92", zorder=0)
ax[0].semilogy(rho[ok], pe[ok], "-o", color="k", label="Channels (W 40)"); ax[0].semilogy(rho[~ok], pe[~ok], "o", color="k", mfc="w", label="Channels m ≥ 18 (unreliable)")
ax[0].semilogy(LV, cav(1.39, 9), "-s", color="#a63603", label="Cavities, aspect ≥ 1.4"); ax[0].semilogy(LV, cav(0.99, 1.01), "-s", color="#3b6fb6", label="Circles")
ax[0].set_xlabel(r"Wall density $\rho$"); ax[0].set_ylabel(r"Escape probability $p_\mathrm{esc}$ (per 20 t.u.)"); ax[0].legend(fontsize=9); ax[0].set_title("Escape vs wall density", fontsize=12)
ax[1].semilogy(rho[ok], he[ok], "-o", color="#a63603", label="Escape hazard per wall defect"); ax[1].semilogy(rho[ok], hn[ok], "-o", color="#3b6fb6", label="Entry hazard per bulk $+\\frac{1}{2}$")
ax[1].set_xlabel(r"Wall density $\rho=2m/L$"); ax[1].set_ylabel("Hazard per 1000 t.u."); ax[1].legend(fontsize=9); ax[1].set_title("Channels: exits blocked, entries not", fontsize=12)
ax[2].plot(rho[ok], occ[ok], "-o", color="k", label="Wall occupancy $N_\\mathrm{wall}/2m$"); ax[2].set_ylim(0, 1.05); ax[2].set_xlabel(r"Wall density $\rho=2m/L$"); ax[2].set_ylabel("Wall occupancy")
a2 = ax[2].twinx() if False else None
ax[2].text(0.97, 0.06, "Bulk charge per length: " + ", ".join(f"{v:.4f}" for v in qb[ok][[0, 3, 7]]) + "\n(m = 3, 8, 12)", transform=ax[2].transAxes, ha="right", fontsize=9)
ax[2].set_title("Channels: the layer fills and saturates", fontsize=12)
for a, s in zip(ax, "abc"): a.text(-0.12, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
ax[0].set_ylim(8e-4, 0.4); ax[0].text(0.09, 1.0e-3, "predicted crossover", ha="center", fontsize=8.5, color="0.4")
for e in ("pdf", "png"): fig.savefig(f"fig_crossover.{e}", dpi=180)
