import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
R = [r for r in json.load(open("/home/baloewe/Physics/Overcharged/reports/master/theory/shear_tilt/flatwall_scan.json")) if "tilt" in r]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.3)); cm = plt.get_cmap("viridis")
ss = sorted(set(r["s"] for r in R if r["a_s"] == 20 and r["H"] == 30))
for i, s in enumerate(ss):
    rr = sorted([r for r in R if r["s"] == s and r["a_s"] == 20 and r["H"] == 30], key=lambda r: r["lam"])
    ax[0].plot([r["lam"] for r in rr], [r["tilt"] for r in rr], "o-", ms=3, color=cm(i / (len(ss) - 1)), label=f"$s={s:g}$")
l = np.linspace(1.02, 4.2, 200); ax[0].plot(l, np.degrees(np.arccos(1 / l)), "k--", lw=1, label=r"$\arccos(1/\lambda)$")
ax[0].set_xlabel(r"$\lambda$"); ax[0].set_ylabel("tilt (deg)"); ax[0].legend(fontsize=6.5, frameon=False); ax[0].set_title("(a) flat wall: tilt vs flow alignment", fontsize=9)
x = [np.degrees(r["depth"] * np.sqrt(r["s"]) * (r["lam"] ** 2 - 1) ** -0.25) for r in R]; y = [r["tilt"] - np.degrees(np.arccos(1 / r["lam"])) for r in R]
sc = ax[1].scatter(x, y, s=12, c=[r["lam"] for r in R], cmap="plasma"); cb = fig.colorbar(sc, ax=ax[1], pad=0.02); cb.set_label(r"$\lambda$"); m = max(max(x), max(y)) * 1.05; ax[1].plot([0, m], [0, m], "k:", lw=0.8)
ax[1].set_xlabel("predicted excess tilt (deg)\n" + r"$d\sqrt{\gamma_1 s/K}\,(\lambda^2-1)^{-1/4}$")
ax[1].set_ylabel("measured excess tilt (deg)\n" + r"tilt $-\arccos(1/\lambda)$")
ax[1].set_title("(b) boundary-layer correction, one point per run", fontsize=9)
fig.tight_layout(); fig.savefig("figs/fig_flatwall.png", dpi=200)
