#!/usr/bin/env python3
"""Diagram of the sparse-layer overtilt (user, 2026-11-27).
(a) Ideal: the defect's surroundings at the Leslie angle theta_L -> tail at 2 theta_L. Exact one-constant flat-wall field with winding anchoring
    (as PRL Fig. 3a): theta = 1/2 arg sin(pi (z - z0)/a) + 1/2 arg sin(pi (z - conj z0)/a) + kappa x + pi/2, kappa = pi/a; far field kappa x0.
(b) Overtilt (mechanism = HYPOTHESIS): the defect's own neighbourhood lies inside the sheared layer. There the director takes all orientations, so the
    aligning part of the shear torque averages out and the vorticity rotates the neighbourhood further by delta, held back elastically; the rotation
    relaxes to the flow-aligned far field over ~l_E. Drawn as the lock-consistent field for theta_L + delta minus delta (1 - exp(-r / l_E)):
    far field theta_L, core phase theta_L + delta, tail at 2 (theta_L + delta). delta exaggerated (8 deg) for visibility.
    Measured support (tilt report, Test 6): the band between wall and defects is NOT rotated forward; in MPCD the core leads its 2-3 cell ring by
    +1 to +5 deg and its 3-4.5 cell ring by +5 to +14 deg (passive imposed vortex), i.e. a rotated neighbourhood around each defect.
(c) Test 2 (swirl_thin/test2.txt, k = 4, passive, lambda 2): the overtilt grows with how deep the shear reaches, not with its strength.
-> fig_overtilt_diagram.pdf/.png"""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch
H = os.path.dirname(os.path.abspath(__file__)); a = 1.0; kap = np.pi / a; h = 0.24; thL = np.radians(30.0); LE = 0.32
RED, BLUE, GREY, ORANGE = "#c0392b", "#1f5fa8", "0.55", "#e67e22"


def field(X, Y, x0, extra=0.0):
    z = X + 1j * Y; z0 = x0 + 1j * h
    th = 0.5 * np.angle(np.sin(np.pi * (z - z0) / a)) + 0.5 * np.angle(np.sin(np.pi * (z - np.conj(z0)) / a)) + kap * X + np.pi / 2
    if extra: th = th - extra * (1 - np.exp(-np.abs(z - z0) / LE))
    return th


def tail_dir(x0, extra):
    ph = np.linspace(-np.pi, np.pi, 20000); r = 0.002; z = x0 + r * np.cos(ph) + 1j * (h + r * np.sin(ph))
    th = field(z.real, z.imag, x0, extra); m = np.angle(np.exp(2j * (th - ph)))                  # director along the radius: th = ph (mod pi)
    return ph[np.argmin(np.abs(m))]


def panel(ax, tilt_core, extra, title, shade=False):
    x0 = tilt_core / kap                                                     # lock: core phase = kappa x0
    xs = np.arange(-0.62, 0.66, 0.062); ys = np.arange(0.04, 0.80, 0.062); X, Y = np.meshgrid(xs, ys); TH = field(X, Y, x0, extra)
    far = np.abs((X - x0) + 1j * (Y - h)) > 0.075
    if shade:
        ax.add_patch(Circle((x0, h), 1.0 * LE, color=ORANGE, alpha=0.16, lw=0, zorder=0))
        ax.text(x0, h + LE + 0.02, "Rotated neighbourhood ($\\sim\\ell_E$)", color=ORANGE, fontsize=11, ha="center", va="bottom", bbox=dict(fc="w", ec="none", alpha=0.85, pad=1))
    L = 0.024
    for x, y, t, f in zip(X.ravel(), Y.ravel(), TH.ravel(), far.ravel()):
        if f: ax.plot([x - L * np.cos(t), x + L * np.cos(t)], [y - L * np.sin(t), y + L * np.sin(t)], color="0.25", lw=1.1, solid_capstyle="round")
    ax.fill_between([-0.7, 0.7], -0.07, 0, color="0.85", zorder=2); ax.plot([-0.7, 0.7], [0, 0], color="k", lw=2)
    psi = tail_dir(x0, extra)                                                # tail direction (from the core)
    ax.plot([x0, x0 + 0.15 * np.cos(psi)], [h, h + 0.15 * np.sin(psi)], color=RED, lw=3.2, solid_capstyle="round", zorder=5)
    ax.plot(x0, h, "o", color=RED, ms=9, mec="w", mew=1, zorder=6)
    ax.plot([x0, x0], [h, h - 0.2], color=GREY, ls="--", lw=1, zorder=4)     # inward normal (towards the wall)
    tau = np.degrees(np.angle(np.exp(1j * (psi + np.pi / 2))))
    ax.add_patch(Arc((x0, h), 0.2, 0.2, theta1=-90 + min(0, tau), theta2=-90 + max(0, tau), color=RED, lw=1.4, zorder=5))
    ax.text(x0 + 0.13, h - 0.17, f"$\\tau$ = {abs(tau):.0f}°", color=RED, fontsize=11, ha="left")
    for yy in (0.25, 0.45, 0.65):                                            # shear: flow grows away from the wall
        ax.add_patch(FancyArrowPatch((0.56, yy), (0.56 + 0.16 * yy / 0.65, yy), arrowstyle="-|>", mutation_scale=9, color=BLUE, lw=1.2))
    ax.text(0.6, 0.72, "Shear", color=BLUE, fontsize=10, ha="center")
    ax.plot([-0.55, -0.55 + 0.11 * np.cos(thL)], [0.7, 0.7 + 0.11 * np.sin(thL)], color=BLUE, lw=2.2)
    ax.text(-0.42, 0.74, "Far field at $\\theta_L$", color=BLUE, fontsize=10, va="center")
    ax.set_xlim(-0.66, 0.8); ax.set_ylim(-0.06, 0.86); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=12.5)


fig = plt.figure(figsize=(20, 7.2)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.7], wspace=0.22)
a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(3))
panel(a1, thL, 0.0, "(a) Ideal: tail at $2\\theta_L$")
d = np.radians(8.0); panel(a2, thL + d, d, "(b) Overtilt: a rotated neighbourhood", shade=True)
a2.text(0.0, -0.02, "HYPOTHESIS: near the core the director takes all angles,\nso the aligning torque averages out and the vorticity\nrotates the neighbourhood further by $\\delta$; elasticity\nties it to the far field over $\\sim\\ell_E$. Tail: $2(\\theta_L+\\delta)$.\n($\\delta$ = 8° here, exaggerated)", transform=a2.transAxes, fontsize=9.6, va="top", ha="left", color="0.25")
# (c) data
labs = ["Outer\n3 cells\n(0.155)", "Outer\n~7 cells\n(0.141)", "Whole\ncavity\n(0.072)"]; tilt = [61.3, 70.2, 77.3]; ws = [0.155, 0.141, 0.072]
cols = ["#f5b041", "#e67e22", "#a04000"]
a3.bar(range(3), tilt, color=cols, width=0.62)
for i, (t, w) in enumerate(zip(tilt, ws)):
    a3.text(i, t + 0.8, f"{t:.0f}°", ha="center", fontsize=13)
a3.axhline(60, color=BLUE, ls="--", lw=1.5); a3.text(-0.4, 80.5, "Dashed: $2\\theta_L$ = 60°", color=BLUE, fontsize=11, ha="left", va="top")
a3.set_xticks(range(3)); a3.set_xticklabels(labs, fontsize=11.5); a3.set_xlabel("How deep the shear reaches\n(wall shear rate in brackets)"); a3.set_ylim(45, 82); a3.set_ylabel("Wall-defect tilt $\\tau$ (deg)")
a3.set_title("(c) Measured: overtilt vs shear depth", fontsize=12.5)
a3.text(0.0, -0.25, "Passive MPCD, imposed vortex, k = 4 (sparse), λ = 2.\nThe shallowest shear has the strongest wall shear\nbut the smallest overtilt. Active sparse layers: 63–66°.\nCrowded layers: no overtilt (they lock to the band).",
        transform=a3.transAxes, fontsize=9.4, va="top", color="0.25")
for s in ("top", "right"): a3.spines[s].set_visible(False)
fig.savefig(os.path.join(H, "fig_overtilt_diagram.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_overtilt_diagram.png"), dpi=170, bbox_inches="tight")
print("wrote fig_overtilt_diagram")
