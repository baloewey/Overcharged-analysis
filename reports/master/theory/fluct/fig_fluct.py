#!/usr/bin/env python3
"""Figure: passive wall-lattice fluctuations, harmonic elastic theory (fluct_theory.json) vs measurement (fluct_measured.json).
(a) rms tangential (internal) and normal displacement vs N in circles R = 31.83; theory with K = 32.5 kBT (dashed) and with the fitted
K_eff (solid, one parameter for all circles and ellipses, tangential only). (b) tangential spectrum per Fourier mode m, circle N = 12.
(c) ellipses: internal tangential rms, theory vs measured. Writes fluct.png/pdf and fluct_fit.txt."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
T = json.load(open(os.path.join(H, "fluct_theory.json"))); M = json.load(open(os.path.join(H, "fluct_measured.json")))
labs = [f"circle_R31.83_k{k:02d}" for k in range(2, 11)]
ells = ["ellipse_eps1.2_k06", "ellipse_eps1.5_k06", "ellipse_eps2.0_k06"]
def tint(l): return T[l]["sig_t_int"] if "sig_t_int" in T[l] else T[l]["sig_t"]
def mint(l): return np.sqrt(np.mean(np.square(M[l]["sig_t_internal"])))
allL = labs + ["circle_R27.28_k06"] + ells
r = np.array([mint(l) / tint(l) for l in allL])
fac = np.exp(np.mean(np.log(r))); Keff = 32.5 / fac ** 2
txt = [f"tangential (internal) rms, measured/theory at K = 32.5: " + ", ".join(f"{l}: {x:.2f}" for l, x in zip(allL, r)),
       f"geometric-mean factor {fac:.3f} (spread {r.std() / r.mean():.1%}) -> K_eff = {Keff:.1f} kBT",
       "normal rms, measured vs theory(K_eff): " + ", ".join(f"{l}: {np.sqrt(np.mean(np.square(M[l]['sig_n']))):.2f} vs {T[l].get('sig_n_int', T[l]['sig_n']) * fac:.2f}" for l in allL),
       "collective (sliding) rms measured: " + ", ".join(f"{l}: {M[l]['sig_coll']:.2f}" for l in allL),
       "ellipse sliding mode, theory: " + ", ".join(f"{l}: eigenvalue {T[l]['sliding_eig']:.2e} (internal modes from {T[l]['second_eig']:.2e}), rms {T[l]['sig_coll']:.1f} cells" for l in ells)]
open(os.path.join(H, "fluct_fit.txt"), "w").write("\n".join(txt) + "\n"); print("\n".join(txt))
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
N = [T[l]["N"] for l in labs]
a = ax[0]
a.plot(N, [mint(l) for l in labs], "o", color="C0", label="tangential, measured")
a.plot(N, [tint(l) * fac for l in labs], "-", color="C0", label=f"tangential, theory ($K_\\mathrm{{eff}}$ = {Keff:.0f} $k_BT$)")
a.plot(N, [tint(l) for l in labs], "--", color="C0", lw=0.8, label="tangential, theory ($K$ = 32.5 $k_BT$)")
a.plot(N, [np.sqrt(np.mean(np.square(M[l]["sig_n"]))) for l in labs], "s", color="C3", label="normal, measured")
a.plot(N, [T[l]["sig_n"] * fac for l in labs], "-", color="C3", label="normal, theory ($K_\\mathrm{eff}$)")
a.set_xlabel("wall defects $N = 2k$"); a.set_ylabel("rms displacement (cells)"); a.set_title("(a) circles $R = 31.83$, passive", fontsize=10)
a.legend(fontsize=7); a.set_ylim(0, None)
a = ax[1]; l = "circle_R31.83_k06"; m = np.arange(len(M[l]["spec_t"]))
a.plot(m[1:], np.array(M[l]["spec_t"])[1:], "o", label="tangential, measured")
a.plot(m[1:], np.array(T[l]["spec_t"])[1:] * fac ** 2, "-", color="C0", label="tangential, theory ($K_\\mathrm{eff}$)")
a.plot(m[1:], np.array(M[l]["spec_n"])[1:], "s", color="C3", label="normal, measured")
a.plot(m[1:], np.array(T[l]["spec_n"])[1:] * fac ** 2, "-", color="C3", label="normal, theory ($K_\\mathrm{eff}$)")
a.axhline(M[l]["spec_t"][0], color="C0", ls=":", lw=0.8); a.text(2.2, M[l]["spec_t"][0] * 0.93, "m = 0 (sliding), measured; zero mode in theory", fontsize=7)
a.set_xlabel("Fourier mode $m$"); a.set_ylabel(r"$\langle|u_m|^2\rangle$ (cells$^2$)"); a.set_title("(b) circle $N = 12$: spectrum", fontsize=10)
a.legend(fontsize=7)
a = ax[2]; eps = [1.0, 1.2, 1.5, 2.0]; L4 = ["circle_R27.28_k06"] + ells
a.plot(eps, [mint(l) for l in L4], "o", label="tangential internal, measured")
a.plot(eps, [tint(l) * fac for l in L4], "-", color="C0", label="tangential internal, theory ($K_\\mathrm{eff}$)")
a.plot(eps, [M[l]["sig_coll"] for l in L4], "^", color="k", label="collective sliding, measured")
a.plot(eps[1:], [T[l]["sig_coll"] for l in ells], "v--", color="gray", label="collective sliding, theory")
a.set_yscale("log"); a.set_xlabel(r"aspect ratio $\varepsilon$ (circle: $R = 27.28$)"); a.set_ylabel("rms (cells)")
a.set_title("(c) $N = 12$: ellipses", fontsize=10); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(H, "fluct.png"), dpi=150); fig.savefig(os.path.join(H, "fluct.pdf"))
