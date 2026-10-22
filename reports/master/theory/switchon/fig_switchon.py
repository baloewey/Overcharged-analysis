#!/usr/bin/env python3
"""Time course of the passive -> active switch (switchon.npz): (a) wall +1/2 count, (b) circulation <v_phi>/<|v|> (signed by the final chirality),
(c) layer tilt (signed), (d) registry arg Z_10 and lock-corrected D (signed), against time after the switch; 8 initial conditions. 20 t.u. running
means for (a, c, d). -> fig_switchon.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z = np.load("switchon.npz"); runs = sorted({k.rsplit("_", 1)[0] for k in Z.files if k.endswith("_t")})
def sm(t, y, w=20.0):
    out = np.full(len(t), np.nan)
    for i, ti in enumerate(t):
        m = (np.abs(t - ti) <= w / 2) & np.isfinite(y)
        if m.any(): out[i] = np.mean(y[m])
    return out
fig, ax = plt.subplots(2, 2, figsize=(13, 9)); ax = ax.ravel(); T0 = 3000.0; XL = (-100, 1200)
for j, r in enumerate(runs):
    t, nw, tau, z10, d10, tf, cf = (Z[f"{r}_{q}"] for q in ("t", "nw", "tau", "z10", "d10", "tf", "circ"))
    sg = np.sign(np.nanmean(cf[tf > tf.max() - 1000])); col = f"C{j}"; sel = t < T0 + XL[1] + 50
    ax[0].plot(t[sel] - T0, sm(t[sel], nw[sel].astype(float)), color=col, lw=1.3, label=r.replace("_", " "))
    ax[1].plot(tf - T0, sg * cf, color=col, lw=1.3)
    ax[2].plot(t[sel] - T0, np.degrees(sm(t[sel], sg * tau[sel])), color=col, lw=1.3)
    ax[3].plot(t[sel] - T0, sm(t[sel], sg * z10[sel]), color=col, lw=1.3)
for a in ax: a.axvline(0, color="k", ls=":", lw=1); a.set_xlim(*XL); a.set_xlabel("Time After Activity Is Switched On (t.u.)")
ax[0].axhline(10, color="0.5", ls="--", lw=1); ax[0].axhline(12, color="0.5", ls="--", lw=1); ax[0].set_ylabel(r"Wall $+\frac{1}{2}$ Count"); ax[0].set_title("Wall Count", pad=10)
ax[0].legend(fontsize=9, ncol=2, loc="lower right"); ax[0].set_ylim(0, 14.5)
ax[1].set_ylabel(r"Circulation $\langle v_\phi\rangle/\langle|v|\rangle$ (signed)"); ax[1].set_title("Coherent Flow", pad=10); ax[1].set_ylim(-0.4, 1.05)
ax[2].axhline(np.degrees(np.arccos(0.5)), color="0.5", ls="--", lw=1); ax[2].text(1180, 40, r"$\arccos(1/\lambda)=60^\circ$", ha="right", color="0.4")
ax[2].set_ylabel(r"Layer Tilt $\tau$ (deg, signed)"); ax[2].set_title("Tilt", pad=10); ax[2].set_ylim(-60, 100)
ax[3].axhline(1.17, color="0.5", ls="--", lw=1); ax[3].text(1180, 2.0, "always-active runs: 1.17", ha="right", color="0.4")
ax[3].set_ylabel(r"Registry arg $Z_{10}$ (rad, signed)"); ax[3].set_title("Registry", pad=10); ax[3].set_ylim(-0.5, 2.6)
for a, l in zip(ax, "abcd"): a.text(-0.13, 1.04, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.03, r"Circle $k=6$, $R=31.8$: passive for 3000 t.u. (registry frozen by 8 different initial directors), then activity 0.07. Signs set by each run's final circulation." + "\n" + r"$Z_{10}=\langle e^{10i\phi_j}\rangle$ over wall $+\frac{1}{2}$; registry shown only in frames with exactly 10 wall defects; lock: arg $Z_{10}=\tau$ (rad).", ha="center", fontsize=11)
fig.tight_layout(); fig.savefig("fig_switchon.pdf", bbox_inches="tight"); fig.savefig("fig_switchon.png", dpi=150, bbox_inches="tight")
