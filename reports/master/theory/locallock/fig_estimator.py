#!/usr/bin/env python3
"""What made the near-wall-field analyses disagree (2026-10-22). (a) Per run: near-wall field from all frames vs coherent frames only (same window,
signed by the circulation): they agree, so incoherent frames are not the cause. (b) Profiles of the ring-averaged field 2 alpha(u) in 0.25-cell rings,
coherent frames, sparse and crowded walls: in sparse walls the value jumps by 10-20 deg between rings half a cell apart (the ring samples the layer's
periodic texture and the discrete cell positions), so 'the field one cell inside the layer' depends on where the ring is placed. (c) The measured
director minus the static texture of the same defects (bandrot): smooth, because the pattern cancels. -> fig_estimator.{pdf,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = np.load("coherence_artifact.npz"); P = np.load("coherence_profiles.npz"); B = np.load("../bandrot/band_rotation.npz")
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
sc = ax[0].scatter(A["f_coh"], A["f_all"], c=np.log10(A["a"]), cmap="viridis", s=45, edgecolors="k", linewidths=0.4)
x = np.array([30, 65]); ax[0].plot(x, x, "k--", lw=1.2, label="Equal"); ax[0].legend(fontsize=10, loc="upper left")
cb = fig.colorbar(sc, ax=ax[0]); cb.set_label(r"$\log_{10}$ Wall-Defect Spacing (cells)")
ax[0].set_xlabel(r"Near-Wall Field $2\alpha_1$, Coherent Frames (deg)"); ax[0].set_ylabel(r"Near-Wall Field $2\alpha_1$, All Frames (deg)")
ax[0].set_title("Coherent vs All Frames", pad=10)
uc, a, d = P["uc"], P["a"], P["d"]
for (lo, hi), c, lab in (((0, 8), "C3", "crowded (spacing < 8)"), ((20, 80), "C0", "sparse (spacing > 20)")):
    m = (a >= lo) & (a < hi); prof = np.degrees(np.angle(np.nanmean(P["P_coh"][m], 0))); dm = np.median(d[m])
    ax[1].plot(uc, prof, "-o", ms=4, color=c, label=f"{lab}, depth {dm:.1f}"); ax[1].axvspan(dm + 0.5, dm + 1.5, color=c, alpha=0.12)
ax[1].set_xlabel("Distance From the Wall $u$ (cells)"); ax[1].set_ylabel(r"Ring-Averaged Field $2\alpha(u)$ (deg)"); ax[1].set_xlim(0.5, 6); ax[1].set_ylim(20, 110)
ax[1].set_title("Field Depends on the Ring", pad=10); ax[1].legend(fontsize=9, loc="upper right")
ub = B["uc"]
for r, c, lab in (("swirl_crowd/k04_s1", "C0", "k 4 (sparse)"), ("swirl_crowd/k08_s2", "C3", "k 8 (crowded)"), ("swirl_shear/k06_A0.010", "C2", "k 6")):
    ax[2].plot(ub, B[f"{r}|all"], "-o", ms=4, color=c, label=f"Passive {lab}")
ax[2].axhline(0, color="k", lw=0.8); ax[2].set_xlim(0, 6); ax[2].set_xlabel("Distance From the Wall $u$ (cells)"); ax[2].set_ylabel(r"Measured $-$ Static Texture, $\delta\theta$ (deg)")
ax[2].set_title("Measured Minus Static: Smooth", pad=10); ax[2].legend(fontsize=9, loc="lower left")
for q, l in zip(ax, "abc"): q.text(-0.14, 1.04, f"({l})", transform=q.transAxes, fontsize=15, fontweight="bold")
fig.text(0.5, -0.05, r"Active MPCD main sweep (39 runs, director every 100 t.u., second half). $2\alpha$: angle of the order-weighted ring average of $e^{2i(\theta-\phi-\pi/2)}$, signed by the circulation." + "\n" + r"Shaded in (b): the window $[d+0.5, d+1.5)$ used for $\alpha_1$. (c) From the band-rotation analysis (passive imposed vortex): measured director minus the exact static field of the same defects.", ha="center", fontsize=10)
fig.tight_layout(); fig.savefig("fig_estimator.pdf", bbox_inches="tight"); fig.savefig("fig_estimator.png", dpi=150, bbox_inches="tight")
