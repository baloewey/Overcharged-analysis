import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
R = json.load(open("scan.json")); fig, ax = plt.subplots(1, 3, figsize=(14, 4.8))
# (a) tilt vs lambda: model (ejected, coherent runs), MPCD passive swirl ranges, 2 theta_L
lm = [(r["lam"], r["tilt"]) for r in R if r["variant"] == "noslip" and r["coh"] > 0.9]
ax[0].plot(*zip(*lm), "o", color="#c2452d", label="Q-tensor model (no-slip flow)")
for l, lo, hi in ((0.5, 28, 36), (1.2, 57, 65), (2.0, 57, 68), (4.0, 75, 83)): ax[0].plot([l, l], [lo, hi], color="#1f5fa8", lw=4, alpha=0.6, label="MPCD passive swirl" if l == 0.5 else None)
ll = np.linspace(1.01, 4.2, 100); ax[0].plot(ll, np.degrees(np.arccos(1 / ll)), "k--", lw=0.8, label=r"$2\theta_L=\arccos(1/\lambda)$")
ax[0].set_xlabel(r"flow-alignment parameter $\lambda$"); ax[0].set_ylabel("tilt of the wall layer (deg)"); ax[0].legend(fontsize=7, frameon=False); ax[0].set_title("(a) tilt vs flow alignment", fontsize=9)
# (b) U scan lam 2
us = sorted([r for r in R if r["lam"] == 2.0 and r["variant"] == "noslip" and r["ring"] == 12], key=lambda r: r["U"]); us = list({r["U"]: r for r in us}.values())
U = [r["U"] for r in us]; ax[1].semilogx(U, [r["coh"] for r in us], "o-", color="#3a8a4f", label="tilt coherence")
ax1 = ax[1].twinx(); ax1.semilogx(U, [r["nwall"] for r in us], "s--", color="0.3", label="wall $+1/2$"); ax1.set_ylabel("wall $+1/2$"); ax1.set_ylim(9.5, 12.5)
ax[1].set_xlabel("imposed flow strength $U$"); ax[1].set_ylabel("tilt coherence"); ax[1].set_title(r"(b) $\lambda=2$: ejection at a threshold", fontsize=9)
ax[1].legend(loc="upper left", fontsize=7, frameon=False); ax1.legend(loc="lower right", fontsize=7, frameon=False)
# (c) variants: colour = flow profile (legend), groups = lambda
from matplotlib.patches import Patch
VN = {"noslip": "measured no-slip profile", "wallonly": "shear only in the wall band", "slip": "no shear in the wall band (wall slips)"}; VC = {"noslip": "#c2452d", "wallonly": "#d08a1e", "slip": "#1f5fa8"}
for gi, l in enumerate((2.0, 4.0)):
    x0 = gi * 4.5
    ax[2].axvspan(x0 - 0.6, x0 + 2.6, color="0.94", zorder=0)                       # group background
    for vi, v in enumerate(("noslip", "wallonly", "slip")):
        r = [x for x in R if x["lam"] == l and x["U"] == 1.0 and x["variant"] == v and x["ring"] == 12][0]
        ax[2].bar(x0 + vi, r["tilt"], color=VC[v], zorder=2); ax[2].text(x0 + vi, r["tilt"] + 1.5, f"{r['tilt']:.0f}°", ha="center", fontsize=7, zorder=3)
    ref = np.degrees(np.arccos(1 / l)); ax[2].plot([x0 - 0.5, x0 + 2.5], [ref, ref], "k--", lw=0.8, zorder=3)
ax[2].plot([], [], "k--", lw=0.8, label=r"$\arccos(1/\lambda)$")
ax[2].set_xticks([1, 5.5]); ax[2].set_xticklabels([r"$\lambda=2$", r"$\lambda=4$"], fontsize=10); ax[2].tick_params(axis="x", length=0)
ax[2].set_xlim(-0.8, 7.3); ax[2].set_ylabel("tilt (deg)"); ax[2].set_ylim(0, 112)
ax[2].legend(handles=[Patch(color=VC[v], label=VN[v]) for v in VN] + [plt.Line2D([], [], color="k", ls="--", lw=0.8, label=r"$\arccos(1/\lambda)$")], fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1)
ax[2].set_title("(c) flow variants at $U=1$ (above threshold)", fontsize=9)
fig.tight_layout(); fig.savefig("shear_tilt.png", dpi=200); fig.savefig("../../figs/shear_tilt.png", dpi=200)
