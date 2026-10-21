#!/usr/bin/env python3
"""The Frank anisotropy of 2D MPCD is set by the flow coupling (2026-10-21). Passive periodic 64^2 boxes, 2 seeds each: production coupling
(shearSusc 0.5) at flow alignment lambda = -2, 0, 1, 2 (runs/Kcalib_lambda tum*, runs/Kcalib L64), and with the shear coupling off (shearSusc 0, lambda 2).
Constants: director-frame maximum-likelihood fit of the fluctuation spectrum, 0 < |q| < 0.8 (runs/Kcalib/frank_constant.py); data/frank_cache.pkl.
(a) K1 (splay) and K3 (bend) against lambda; band: coupling off (both constants, mean +- range). (b) K1/K3 against lambda, with the
coupling-off value. (c) Effective constant per mode, K_eff(q), coupling off vs on (lambda 2), seeds averaged.
-> figs/frank_coupling.{pdf,png}, data/frank_coupling_summary.txt"""
import pickle, numpy as np, matplotlib, matplotlib.ticker; matplotlib.use("Agg"); import matplotlib.pyplot as plt
C = pickle.load(open("data/frank_cache.pkl", "rb"))
SETS = {-2: ["Kcalib_lambda/L64_tumm2_s1", "Kcalib_lambda/L64_tumm2_s2"], 0: ["Kcalib_lambda/L64_tum0_s1", "Kcalib_lambda/L64_tum0_s2"],
        1: ["Kcalib_lambda/L64_tum1_s1", "Kcalib_lambda/L64_tum1_s2"], 2: ["Kcalib/L64_s1", "Kcalib/L64_s2"]}
OFF = ["Kcalib_lambda/L64_susc0_s1", "Kcalib_lambda/L64_susc0_s2"]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8)); out = ["setting | seed K1 K3 K1/K3 (director-frame ML, qmax 0.8)"]
off = np.array([[C[r]["K1"], C[r]["K3"]] for r in OFF]); lo, hi = off.min(), off.max()
ax[0].axhspan(lo, hi, color="0.85", zorder=0); ax[0].text(-1.55, hi + 1.2, f"Coupling off: {off.mean():.1f} $k_BT$ (both)", fontsize=11, color="0.3")
for lam, rr in SETS.items():
    for j, r in enumerate(rr):
        K1, K3 = C[r]["K1"], C[r]["K3"]; out.append(f"lambda {lam:+d}, shearSusc 0.5 | {r} {K1:.2f} {K3:.2f} {K1/K3:.3f}")
        ax[0].plot(lam - 0.06, K1, "o", color="C3", ms=8, mfc="none" if j else "C3", label="$K_1$ (splay)" if (lam, j) == (-2, 0) else None)
        ax[0].plot(lam + 0.06, K3, "s", color="C0", ms=8, mfc="none" if j else "C0", label="$K_3$ (bend)" if (lam, j) == (-2, 0) else None)
        ax[1].plot(lam, K1 / K3, "D", color="k", ms=8, mfc="none" if j else "k")
for r in OFF: out.append(f"lambda +2, shearSusc 0 | {r} {C[r]['K1']:.2f} {C[r]['K3']:.2f} {C[r]['K1']/C[r]['K3']:.3f}")
ax[1].plot([2.25, 2.25], [C[r]["K1"] / C[r]["K3"] for r in OFF], "D", color="C2", ms=8, label="Coupling off ($\\lambda=2$)")
ax[1].axhline(1, color="0.5", lw=1, ls=":")
ax[0].set_xlabel(r"Flow-Alignment Parameter $\lambda$"); ax[0].set_ylabel(r"Frank Constant from Fluctuations ($k_BT$)"); ax[0].set_ylim(15, 52)
ax[0].set_title("Splay and Bend vs Flow Alignment", pad=10); ax[0].legend(fontsize=11, loc="lower center"); ax[0].set_xticks([-2, -1, 0, 1, 2])
ax[1].set_yscale("log"); ax[1].set_yticks([0.5, 1, 2]); ax[1].set_yticklabels(["0.5", "1", "2"]); ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_ylim(0.4, 2.5); ax[1].set_xticks([-2, -1, 0, 1, 2])
ax[1].set_xlabel(r"Flow-Alignment Parameter $\lambda$"); ax[1].set_ylabel(r"Anisotropy $K_1/K_3$"); ax[1].set_title("Anisotropy", pad=10); ax[1].legend(fontsize=11, loc="upper left")
for rr, ls, nm in ((OFF, "--", "off"), (SETS[2], "-", "on")):
    for key, col in (("splay", "C3"), ("bend", "C0")):
        q = np.array([a for a, _ in C[rr[0]][key]]); ke = np.mean([[b for _, b in C[r][key]] for r in rr], axis=0)
        ax[2].plot(q, ke, ls, marker="o", ms=4, color=col, label=f"{key.capitalize()}, {nm}")
ax[2].set_xlabel(r"Wavenumber $|q|$ (cells$^{-1}$)"); ax[2].set_ylabel(r"Effective Constant $K_\mathrm{eff}(q)$ ($k_BT$)"); ax[2].set_ylim(0, 80); ax[2].set_xlim(0, 0.82)
ax[2].set_title("Effective Constant per Mode", pad=10); ax[2].legend(fontsize=10, loc="upper center", ncol=2, title=r"Coupling ($\lambda=2$)", title_fontsize=10)
for a, l in zip(ax, "abc"): a.text(-0.15, 1.05, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.03, r"Passive periodic $64^2$ boxes, density 20, $U=10\,k_BT$, 2 seeds (filled/open). Coupling: rotation of the director by velocity gradients (Jeffery), shear susceptibility $\chi_\mathrm{HI}=0.5$ (on) or 0 (off).", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("figs/frank_coupling.pdf", bbox_inches="tight"); fig.savefig("figs/frank_coupling.png", dpi=150, bbox_inches="tight")
open("data/frank_coupling_summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
