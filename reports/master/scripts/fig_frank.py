#!/usr/bin/env python3
"""Frank constants of the 2D MPCD nematic: director-fluctuation spectrum of the passive periodic boxes (runs/Kcalib L64_s1, L64_s2, L96_s1).
Equipartition for a director along x: A kBT/<|theta_q|^2> = K1 q_y^2 + K3 q_x^2 + c (c: white cell-level noise), fitted over 0 < |q| < 0.8.
(a) A kBT/<|theta_q|^2> against q^2 along the bend (q_x) and splay (q_y) axes, with the fitted lines; (b) K_eff(q) = A kBT/(q^2 <|theta_q|^2>) along the axes (raw, as in the data file; one mode per point, so noisy at small q).
-> figs/frank_spectrum.{pdf,png}, data/frank_constants.txt"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/Kcalib"); import frank_constant as FC
RUNS = ["L64_s1", "L64_s2", "L96_s1"]; QMAX = 0.8
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4)); out = []
mk = {"L64_s1": "o", "L64_s2": "s", "L96_s1": "^"}
for r in RUNS:
    d = FC.analyse(f"/home/baloewe/Physics/Overcharged/runs/Kcalib/{r}", QMAX, 200.0)
    out.append(f"runs/Kcalib/{r}: L {d['L']}, {d['frames']} frames, <S> {d['S']:.3f};  K1 (splay) = {d['K1']:.2f} kBT,  K3 (bend) = {d['K3']:.2f} kBT,  "
               f"mean {0.5*(d['K1']+d['K3']):.2f} kBT;  noise constant {d['c']:.3g}")
    out.append("   K_eff along q_x (bend): " + " ".join(f"{q:.2f}:{k:.2f}" for q, k in d["bend"][:6]))
    out.append("   K_eff along q_y (splay): " + " ".join(f"{q:.2f}:{k:.2f}" for q, k in d["splay"][:6]))
    lab = f"{d['L']}$^2$, seed {r[-1]}"
    for key, K, col in (("bend", d["K3"], "C0"), ("splay", d["K1"], "C3")):
        q = np.array([a for a, _ in d[key]]); ke = np.array([b for _, b in d[key]]); y = ke * q * q
        ax[0].plot(q**2, y, mk[r], color=col, mfc="none" if d["L"] == 64 else col, ms=6, label=f"{key.capitalize()}, {lab}")
        ax[1].plot(q, ke, mk[r], color=col, mfc="none" if d["L"] == 64 else col, ms=6)
        qq = np.linspace(0, QMAX, 50)
        ax[0].plot(qq**2, K * qq**2 + d["c"], "-", color=col, lw=1, alpha=0.6)
for K, col, nm in ((42.0, "C3", "$K_1$ (splay) = 42"), (23.0, "C0", "$K_3$ (bend) = 23")):
    ax[1].axhline(K, color=col, lw=1.2, ls="--"); ax[1].text(0.79, K + 6.5, nm + r" $k_BT$", color=col, ha="right", fontsize=12)
ax[0].set_xlabel(r"Squared Wavenumber $q^2$ (cells$^{-2}$)"); ax[0].set_ylabel(r"$A\,k_BT/\langle|\theta_{\mathbf{q}}|^2\rangle$ ($k_BT$)")
ax[0].set_title("Fluctuation Spectrum Along the Axes", pad=10); ax[0].legend(fontsize=9, ncol=1, loc="upper left")
ax[1].set_xlabel(r"Wavenumber $|q|$ (cells$^{-1}$)"); ax[1].set_ylabel(r"Effective Constant $K_\mathrm{eff}(q)$ ($k_BT$)")
ax[1].set_title(r"Effective Constant per Mode", pad=10); ax[1].set_ylim(0, 65); ax[1].set_xlim(0, 0.82)
for a, l in zip(ax, "ab"): a.text(-0.16, 1.05, f"({l})", transform=a.transAxes, fontsize=15, fontweight="bold")
fig.tight_layout(); fig.savefig("figs/frank_spectrum.pdf"); fig.savefig("figs/frank_spectrum.png", dpi=150)
open("data/frank_constants.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
