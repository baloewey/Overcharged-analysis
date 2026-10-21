#!/usr/bin/env python3
"""Figure of the lock valley (lock_valley.npz): (a, b) energy per wall defect against the registry shift (as lock tilt tau_lock) and the orientation
offset delta at k = 6, one constant K = 44 (a) and K1/K3 = 43.1/21.9 (b); (c) valley profile E(tau_lock, delta = 0) for k = 4, 6, 10, 16;
(d) twist stiffness d2E/d delta^2 at delta = 0 (mean over the valley) against the wall-defect spacing. -> fig_lock_valley.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = np.load("lock_valley.npz"); TL, DL = D["TL"], np.array(D["DL"]); KS = [4, 6, 10, 16]; R = 31.830989
fig, ax = plt.subplots(2, 2, figsize=(13, 10.5)); ax = ax.ravel()
for a, key, ttl in ((ax[0], "K44_44_k6", "One Constant $K=44$, $k=6$"), (ax[1], "K43.1_21.9_k6", "$K_1/K_3=43.1/21.9$, $k=6$")):
    M = D[key]; im = a.pcolormesh(TL, DL, M.T, shading="gouraud", cmap="viridis", vmin=0, vmax=20)
    a.set_xlabel(r"Registry Shift as Lock Tilt $\tau_\mathrm{lock}$ (deg)"); a.set_ylabel(r"Orientation Offset $\delta$ (deg)"); a.set_title(ttl, pad=10)
    a.axhline(0, color="w", ls="--", lw=1); a.set_xticks([-180, -90, 0, 90, 180]); fig.colorbar(im, ax=a, label=r"Energy per Wall Defect ($k_BT$)")
for k, c in zip(KS, ["C0", "C1", "C2", "C3"]):
    for key, ls, lab in ((f"K44_44_k{k}", "-", "one constant"), (f"K43.1_21.9_k{k}", "--", "anisotropic")):
        ax[2].plot(TL, D[key][:, list(DL).index(0)], ls, color=c, marker="o", ms=4, label=f"$k={k}$, {lab}")
ax[2].set_xlabel(r"Registry Shift as Lock Tilt $\tau_\mathrm{lock}$ (deg)"); ax[2].set_ylabel(r"Energy per Wall Defect ($k_BT$)")
ax[2].set_title(r"Along the Lock ($\delta=0$)", pad=10); ax[2].set_xticks([-180, -90, 0, 90, 180]); ax[2].legend(fontsize=9, ncol=2, loc="upper center"); ax[2].set_ylim(-2, 22)
for key, mk, lab in (("K44_44", "o", "One constant $K=44$"), ("K43.1_21.9", "s", "$K_1/K_3=43.1/21.9$")):
    sp, tw = [], []
    for k in KS:
        M = D[f"{key}_k{k}"]; i0, ip, im_ = list(DL).index(0), list(DL).index(30), list(DL).index(-30)
        tw.append(np.mean((M[:, ip] + M[:, im_] - 2 * M[:, i0]) / np.radians(30) ** 2)); N = 2 * k - 2
        sp.append(2 * np.pi * R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) / N)
    ax[3].plot(sp, tw, mk + "-", ms=8, label=lab)
ax[3].set_xlabel("Wall-Defect Spacing (cells)"); ax[3].set_ylabel(r"Twist Stiffness $\partial^2E/\partial\delta^2$ ($k_BT$/rad$^2$)")
ax[3].set_title("Across the Lock", pad=10); ax[3].legend(fontsize=10); ax[3].set_ylim(0, 35)
for a, l in zip(ax, "abcd"): a.text(-0.14, 1.04, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.06, r"Circle $R=31.8$, ring of $2k-2$ wall $+\frac{1}{2}$ plus a central $+1$; rigid cores (radius 1.5), field relaxed outside." + "\n" + r"Registry shift $\alpha$ moves every defect along the wall ($\tau_\mathrm{lock}=(2k-2)\alpha$ is the tilt the lock gives);" + "\n" + r"$\delta$ rotates every core away from that tilt at fixed position; actual tilt $\tau=\tau_\mathrm{lock}+\delta$.", ha="center", fontsize=11)
fig.tight_layout(); fig.savefig("fig_lock_valley.pdf", bbox_inches="tight"); fig.savefig("fig_lock_valley.png", dpi=150, bbox_inches="tight")
