import numpy as np, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.argv = ["x"]; exec(open("restoring.py").read().split("groups = ")[0])        # loads data, edges, c, prof()
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(9, 3.2)); col = {1.0: "#1f5fa8", 1.5: "#d08a1e", 2.0: "#3a8a4f"}
for e, ls, lab in ((False, "--", "circle-like"), (True, "-", "elongated")):
    for R in (1.0, 1.5, 2.0):
        mu, var, F, U, n = prof((EL == e) & (RQ == R))
        ax[0].plot(c, -F, ls, marker="o", ms=3, color=col[R], label=f"{lab}, $R_{{eq}}={R}$")
        ax[1].plot(c, U, ls, marker="o", ms=3, color=col[R])
ax[0].axhline(0, color="0.6", lw=0.7); ax[0].set_xlabel("depth $d$ (cells)"); ax[0].set_ylabel(r"$-U'(d)/k_BT_\mathrm{eff}$ (1/cell)"); ax[0].legend(fontsize=6.5, frameon=False)
ax[0].set_title("(a) effective force (+: away from the wall)", fontsize=9); ax[1].set_xlabel("depth $d$ (cells)"); ax[1].set_ylabel(r"$U(d)/k_BT_\mathrm{eff}$")
ax[1].set_title("(b) effective potential of a wall $+1/2$", fontsize=9); fig.tight_layout(); fig.savefig("/home/baloewe/Physics/Overcharged/reports/theory_report/figs/fig_restoring.png", dpi=200)
