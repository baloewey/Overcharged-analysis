#!/usr/bin/env python3
"""Where does the tail of a tilted wall +1/2 meet the wall? (user, 2026-11-27). Exact one-constant flat-wall field with winding anchoring
(theta = 1/2 arg sin(pi (z - z0)/a) + 1/2 arg sin(pi (z - conj z0)/a) + kappa x + pi/2, kappa = pi/a; one +1/2 per site spacing a, at depth h).
The tail is traced as the director field line leaving the core along the tail direction psi (where the director is radial).
(a) untilted defect on a homeotropic site: the tail line runs straight down. (b) tilted, h/a = 0.07 (MPCD k = 6: depth 1.5 cells, spacing ~20):
the tail leaves the core at tau and bends; the anchoring directly under the core is rotated by tau/2 (the lock: tau = 2 kappa D), at the contact point by
less than tau. (c) anchoring at the contact point / tau vs tau for several depths, with the shallow limit 1/2 and the "tail matches anchoring" guess 1.
-> fig_tail_line_diagram.pdf/.png"""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Arc
H = os.path.dirname(os.path.abspath(__file__)); a = 1.0; kap = np.pi / a
RED, BLUE, GREEN, GREY = "#c0392b", "#1f5fa8", "#1e8449", "0.55"


def theta(x, y, x0, h):
    z = x + 1j * y; z0 = x0 + 1j * h
    return 0.5 * np.angle(np.sin(np.pi * (z - z0) / a)) + 0.5 * np.angle(np.sin(np.pi * (z - np.conj(z0)) / a)) + kap * x + np.pi / 2


def tail_dir(x0, h):
    ph = np.linspace(-np.pi, np.pi, 40000); r = 1e-4; th = theta(x0 + r * np.cos(ph), h + r * np.sin(ph), x0, h)
    return ph[np.argmin(np.abs(np.angle(np.exp(2j * (th - ph)))))]


def trace(x0, h, psi, ds=1e-4, nmax=200000):
    x, y = x0 + 1e-3 * np.cos(psi), h + 1e-3 * np.sin(psi); d = np.array([np.cos(psi), np.sin(psi)]); P = [(x0, h), (x, y)]
    for _ in range(nmax):
        t = theta(x, y, x0, h); v = np.array([np.cos(t), np.sin(t)]); v = v if v @ d > 0 else -v; d = v; x += ds * v[0]; y += ds * v[1]; P.append((x, y))
        if y <= 0 or abs(x - x0) > 2 * a: break
    return np.array(P), y <= 0


def panel(ax, tau_deg, h, title, notes=True):
    D = np.radians(tau_deg) / (2 * kap); psi = tail_dir(D, h); tau = psi + np.pi / 2
    x0v, x1v, y1v = D - 0.13, D + 0.21, 0.17; xs = np.linspace(x0v + 0.006, x1v - 0.006, 22); ys = np.linspace(0.008, y1v - 0.006, 11)
    X, Y = np.meshgrid(xs, ys); TH = theta(X, Y, D, h); L = 0.0062
    for x, y, t in zip(X.ravel(), Y.ravel(), TH.ravel()):
        if np.hypot(x - D, y - h) > 0.014: ax.plot([x - L * np.cos(t), x + L * np.cos(t)], [y - L * np.sin(t), y + L * np.sin(t)], color="0.6", lw=1.0, solid_capstyle="round")
    ax.fill_between([x0v, x1v], -0.075, 0, color="0.85"); ax.plot([x0v, x1v], [0, 0], "k", lw=2)
    for xw in np.linspace(x0v + 0.01, x1v - 0.01, 22):                              # anchoring drawn on the wall
        t = np.pi / 2 + kap * xw; ax.plot([xw - L * np.cos(t), xw + L * np.cos(t)], [-0.013 - L * np.sin(t), -0.013 + L * np.sin(t)], color="k", lw=1.1)
    ax.plot(0, 0, "o", ms=8, mfc="w", mec=RED, mew=1.6, zorder=6) if x0v < 0 < x1v else None
    P, hit = trace(D, h, psi); m = P[:, 1] >= 0
    ax.plot(P[m, 0], P[m, 1], color=RED, lw=2.6, zorder=5, label="Tail (field line from the core)")
    ax.plot([D, D + 0.095 * np.cos(psi)], [h, h + 0.095 * np.sin(psi)], color=RED, lw=1.1, ls="--", zorder=4, label="Straight continuation of the tail")
    ax.plot(D, h, "o", color=RED, ms=10, mec="w", mew=1.2, zorder=7)
    ax.plot([D, D], [h, 0], color=GREY, ls=":", lw=1.2, zorder=3)
    if tau_deg > 0:
        ax.add_patch(Arc((D, h), 0.09, 0.09, theta1=-90, theta2=-90 + np.degrees(tau), color=RED, lw=1.5, zorder=6))
        ax.text(D + 0.035, h + 0.004, f"$\\tau$ = {np.degrees(tau):.0f}°", color=RED, fontsize=12, bbox=dict(fc="w", ec="none", alpha=0.85, pad=1))
        xc = P[m][-1, 0] if hit else None; anch_below = np.degrees(kap * D)
        ax.plot(D, 0, "s", color=BLUE, ms=8, zorder=7); ax.text(D - 0.01, -0.017, f"Under the core:\nanchoring {anch_below:.0f}° = $\\tau/2$", color=BLUE, fontsize=10.5, ha="right", va="top")
        if hit:
            ax.plot(xc, 0, "D", color=GREEN, ms=8, zorder=7); ax.text(xc + 0.01, -0.017, f"Tail meets\nthe wall:\nanchoring {np.degrees(kap * xc):.0f}°\n(not $\\tau$)", color=GREEN, fontsize=10.5, va="top", ha="center")
    else:
        ax.text(D + 0.012, 0.035, "Tail straight down;\nanchoring 0° at the site", color=RED, fontsize=11, bbox=dict(fc="w", ec="none", alpha=0.85, pad=1))
    ax.text(0.0 if x0v < 0 else x0v + 0.01, -0.017, "homeotropic site", color=RED, fontsize=8.5, ha="center", va="top") if x0v < 0 < x1v else None
    ax.set_xlim(x0v, x1v); ax.set_ylim(-0.075, y1v); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=12.5)


fig = plt.figure(figsize=(18, 6.2)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.28)
a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(3))
panel(a1, 0, 0.07, "(a) On a homeotropic site: untilted")
panel(a2, 60, 0.07, "(b) Slid by D: tilted by $\\tau = 2\\kappa D$")
a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=10, frameon=False, ncol=1)
# (c) contact-point anchoring / tau vs tau for several depths
for h, c in ((0.03, "#f5b041"), (0.07, "#e67e22"), (0.12, "#a04000")):
    taus = np.arange(5, 81, 5.0); ratio = []
    for td in taus:
        D = np.radians(td) / (2 * kap); psi = tail_dir(D, h); P, hit = trace(D, h, psi, ds=2e-4)
        ratio.append(np.degrees(kap * P[-1, 0]) / np.degrees(psi + np.pi / 2) if hit else np.nan)
    a3.plot(taus, ratio, "o-", color=c, lw=2, ms=4, label=f"depth h = {h:.2f} a")
a3.axhline(0.5, color=BLUE, ls="--", lw=1.5); a3.text(79, 0.47, "under the core: exactly 1/2 (the lock)", color=BLUE, ha="right", va="top", fontsize=10)
a3.axhline(1.0, color=GREY, ls=":", lw=1.5); a3.text(79, 1.02, "'tail matches anchoring' guess", color="0.4", ha="right", va="bottom", fontsize=10)
a3.set_xlabel("Tail tilt $\\tau$ (deg)"); a3.set_ylabel("Anchoring at the tail's wall contact / $\\tau$"); a3.set_ylim(0.3, 1.25); a3.set_xlim(0, 82)
a3.set_title("(c) Anchoring where the tail meets the wall", fontsize=12); a3.legend(fontsize=9.5, frameon=False, loc="upper left")
for s in ("top", "right"): a3.spines[s].set_visible(False)
fig.text(0.5, -0.03, "Exact one-constant field at a flat winding wall (one +½ per site spacing a; MPCD k = 6 sits at h ≈ 0.07 a). Missing points in (c): the tail line does not reach the nearby wall.",
         ha="center", fontsize=10.5, color="0.3")
fig.savefig(os.path.join(H, "fig_tail_line_diagram.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_tail_line_diagram.png"), dpi=170, bbox_inches="tight")
print("wrote fig_tail_line_diagram")
