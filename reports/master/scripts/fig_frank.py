#!/usr/bin/env python3
"""Frank constants of the 2D MPCD nematic at the production coupling (tumble 2, shearSusc 0.5), from the director-fluctuation spectrum of the
passive periodic boxes runs/Kcalib (L64 s1, s2; L96 s1). Each frame's wavevectors are rotated into that frame's director frame; equipartition:
A kBT/<|theta_q|^2> = K1 q_perp^2 + K3 q_par^2 + c (q_par along the director = bend, q_perp across = splay), fitted by maximum likelihood over
0 < |q| < 0.8 (runs/Kcalib/frank_constant.py, dir_ml). Data from data/frank_cache.pkl (scripts/frank_cache.py).
(a) A kBT/<|theta_q|^2> against q^2 for modes within 10 deg of the director (bend) and of its normal (splay), with the fitted lines;
(b) K_eff(q) = A kBT/(q^2 <|theta_q|^2>) for the same modes, dashed: mean fitted K1 and K3.
-> figs/frank_spectrum.{pdf,png}, data/frank_constants.txt"""
import pickle, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
C = pickle.load(open("data/frank_cache.pkl", "rb")); RUNS = ["Kcalib/L64_s1", "Kcalib/L64_s2", "Kcalib/L96_s1"]; QMAX = 0.8
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6)); out = []; mk = dict(zip(RUNS, "os^"))
for r in RUNS:
    d = C[r]; lab = f"{d['L']}$^2$, seed {r[-1]}"
    out.append(f"runs/{r}: L {d['L']}, {d['frames']} frames, <S> {d['S']:.3f}; director frame, ML: K1 (splay) = {d['K1']:.2f} kBT, K3 (bend) = {d['K3']:.2f} kBT, "
               f"K1/K3 = {d['K1']/d['K3']:.2f}; c = {d['c']:.3g}; mean director angle {d['th0_mean']:+.1f} deg (sd {d['th0_sd']:.1f})")
    for key, K, col in (("bend", d["K3"], "C0"), ("splay", d["K1"], "C3")):
        q = np.array([a for a, _ in d[key]]); ke = np.array([b for _, b in d[key]])
        ax[0].plot(q**2, ke * q * q, mk[r], color=col, mfc="none" if d["L"] == 64 else col, ms=6, label=f"{key.capitalize()}, {lab}")
        ax[1].plot(q, ke, mk[r], color=col, mfc="none" if d["L"] == 64 else col, ms=6)
        qq = np.linspace(0, QMAX, 50); ax[0].plot(qq**2, K * qq**2 + d["c"], "-", color=col, lw=1, alpha=0.6)
K1m, K3m = np.mean([C[r]["K1"] for r in RUNS]), np.mean([C[r]["K3"] for r in RUNS])
for K, col, nm, dy in ((K1m, "C3", "$K_1$ (splay)", 5), (K3m, "C0", "$K_3$ (bend)", -9)):
    ax[1].axhline(K, color=col, lw=1.2, ls="--"); ax[1].text(0.79, K + dy, f"{nm} = {K:.1f} $k_BT$", color=col, ha="right", fontsize=12)
ax[0].set_xlabel(r"Squared Wavenumber $q^2$ (cells$^{-2}$)"); ax[0].set_ylabel(r"Inverse Spectrum $A\,k_BT/\langle|\theta_{\mathbf{q}}|^2\rangle$ ($k_BT$)")
ax[0].set_title("Fluctuation Spectrum, Director Frame", pad=10); ax[0].legend(fontsize=9, loc="upper left"); ax[0].set_xlim(0, 0.66)
ax[1].set_xlabel(r"Wavenumber $|q|$ (cells$^{-1}$)"); ax[1].set_ylabel(r"Effective Constant $K_\mathrm{eff}(q)$ ($k_BT$)")
ax[1].set_title("Effective Constant per Mode", pad=10); ax[1].set_ylim(0, 70); ax[1].set_xlim(0, 0.82)
for a, l in zip(ax, "ab"): a.text(-0.16, 1.05, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.02, r"Passive periodic boxes, density 20, $U=10\,k_BT$, flow alignment $\lambda=2$, shear coupling 0.5. $\theta$: director angle from the frame's mean director; $A=L^2$.", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("figs/frank_spectrum.pdf", bbox_inches="tight"); fig.savefig("figs/frank_spectrum.png", dpi=150, bbox_inches="tight")
open("data/frank_constants.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
