#!/usr/bin/env python3
"""Crowding protection, first tests (2026-10-22). (a) elastic cost of dragging one wall +1/2 inward out of the layer, relaxed one-constant field (K 44), wall cores
held at the layer tilt 60 deg (solid) or free (dashed), circle R 31.83, five wall densities (drag_energy.npz). (b) Gauss split at a drag of 6 cells (analytic,
gauss_split.py): uniform (smeared) layer, hole (-1/2 left behind + image), corrugation (discreteness), self-image + central +1, total. (c) rms inward flow vs
wall distance, circles of the shape sweep (wallflow.npz, medians over runs)."""
import csv, numpy as np, matplotlib.pyplot as plt, importlib.util, os
H = os.path.dirname(os.path.abspath(__file__)); S = H + "/../shapes"
RHO = (0.04, 0.06, 0.09, 0.12, 0.14); KS = (4, 6, 9, 12, 14); RCOL = ("#fdae6b", "#f98c3a", "#e6550d", "#a63603", "#4d1a02")
D = np.load(H + "/drag_energy.npz"); DEL = D["DEL"]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6)); fig.subplots_adjust(wspace=0.32, left=0.06, right=0.98, bottom=0.15, top=0.9)
for k, r, c in zip(KS, RHO, RCOL):
    ax[0].plot(DEL, D[f"fixed_k{k}_t60"], "-o", ms=4, color=c, label=fr"$\rho={r:.2f}$"); ax[0].plot(DEL, D[f"free_k{k}_t0"], "--", color=c, lw=1)
ax[0].set_xlabel(r"Drag distance $\delta$ (cells)"); ax[0].set_ylabel(r"Elastic cost $E(\delta)-E(0)$ ($k_BT$)"); ax[0].legend(loc="upper left", fontsize=10)
ax[0].text(0.97, 0.05, "Solid: tilt held at 60°\nDashed: free orientation", transform=ax[0].transAxes, ha="right", fontsize=10)
spec = importlib.util.spec_from_file_location("g", H + "/gauss_split.py"); src = open(H + "/gauss_split.py").read().split("\nL = [\"GAUSS")[0]
g = {}; exec(src, g); R, K = g["R"], g["K"]; parts = {"Uniform layer": [], "Hole + image": [], "Corrugation": [], "Self-image, centre": [], "Total": []}
for k in KS:
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; z0 = rs + 0j; M = 4000
    zr = rs * np.exp(2j * np.pi * np.arange(N) / N); sm = rs * np.exp(2j * np.pi * (np.arange(M) + 0.5) / M)
    def parts_at(dl):
        zp = z0 - dl; d = g["pair"](zp, zr[1:], [0.5] * (N - 1)); s = g["pair"](zp, sm, [0.5 * N / M] * M); h = g["hole"](zp, z0)
        o = -0.25 * np.log(R ** 2 - abs(zp) ** 2) + g["pair"](zp, [0j], [1.0]); return np.array([s, h, d - s - h, o, d + o])
    v = np.pi * K * (parts_at(6.0) - parts_at(0.0))
    for nm, x in zip(parts, v): parts[nm].append(x)
sty = {"Uniform layer": ("0.6", ":"), "Hole + image": ("#3b6fb6", "-"), "Corrugation": ("#6a9f58", "-"), "Self-image, centre": ("#9a6fb0", "--"), "Total": ("k", "-")}
for nm, v in parts.items(): ax[1].plot(RHO, v, marker="o", ms=4, color=sty[nm][0], ls=sty[nm][1], lw=2 if nm == "Total" else 1.3, label=nm)
ax[1].axhline(0, color="0.85", lw=0.8, zorder=0); ax[1].set_xlabel(r"Wall density $\rho=2k/P$"); ax[1].set_ylabel(r"Contribution to $E(6)-E(0)$ ($k_BT$)")
ax[1].legend(loc="upper left", ncol=2, fontsize=9.5, columnspacing=0.8); ax[1].set_ylim(-100, 300); ax[1].set_xticks(RHO)
W = np.load(H + "/wallflow.npz"); SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
rho = np.array([float(SR[r]["rho"]) for r in W["runs"]]); asp = np.array([float(SR[r]["aspect"]) for r in W["runs"]]); lev = np.argmin(np.abs(rho[:, None] - np.array(RHO)), 1)
dc = 0.5 * (W["DB"][:-1] + W["DB"][1:])
for i, (r, c) in enumerate(zip(RHO, RCOL)): ax[2].plot(dc, np.median(W["prof"][(lev == i) & (asp < 1.05), 1], 0), "-", color=c, lw=1.6, label=fr"$\rho={r:.2f}$")
ax[2].set_xlabel("Distance from the wall $d$ (cells)"); ax[2].set_ylabel(r"RMS inward flow $\langle v_\perp^2\rangle^{1/2}$"); ax[2].legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=3, handlelength=1.2, fontsize=9.5, columnspacing=0.8); ax[2].set_xlim(0, 12); ax[2].set_ylim(0, 0.17); ax[2].set_xlim(-0.3, 12)
for a, s in zip(ax, "abc"): a.text(-0.12, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
ax[0].set_title("Drag a wall +1/2 out of the layer", fontsize=12); ax[1].set_title("Gauss split at δ = 6 cells", fontsize=12); ax[2].set_title("Near-wall flow, circles", fontsize=12)
for e in ("pdf", "png"): fig.savefig(H + f"/fig_crowd_protect.{e}", dpi=200)
