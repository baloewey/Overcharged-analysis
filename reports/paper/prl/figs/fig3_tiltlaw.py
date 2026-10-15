#!/usr/bin/env python3
"""PRL Fig. 3 (tilt law), full page width (2026-10-15).
(a) Flat-wall schematic, drawn from the exact one-constant solution of a wall with winding anchoring theta(x, 0) = pi/2 + kappa x (kappa = pi/a_s, one +1/2 per
    period a_s at depth h, same-sign images):  theta = (1/2) arg sin(pi (z - z0)/a_s) + (1/2) arg sin(pi (z - conj z0)/a_s) + kappa x + pi/2.
    Its far field is theta_inf = kappa x0 (zero mode: sliding the defect by D rotates the far field by kappa D and the defect by 2 kappa D). The defect is placed where the
    far field sits at the Leslie angle of the shear, theta_L = (1/2) arccos(1/lambda), lambda = 2 (D = theta_L / kappa); its tilt from the wall normal is then tau = 2 theta_L.
(b) Tilt |tau| of the wall +1/2 vs lambda in active MPCD circles (runs/lambda_circle; master/theory/energetics/tiltlaw.txt, mean +- s.e. over seeds) with arccos(1/lambda);
    inset: wall +1/2 count minus 2k - 2 (0 = full overcharged ring).
(c) Passive Q-tensor cavity (k = 6) under an imposed vortex (U = 1, above the ejection threshold), three flow profiles (master/theory/shear_tilt/scan.json, ring 12).
Run: figs/run.sh fig3_tiltlaw.py -> fig3_tiltlaw.pdf/.png"""
import os, json, numpy as np
from prlstyle import plt, label, save, PAGE, KCOL, RED, BLUE, GREY, DARK, ROOT
from matplotlib.patches import Arc, Patch, FancyArrowPatch

TH = ROOT + "/reports/master/theory"


def read_tiltlaw():
    """k, lambda -> (wall count, mean tilt, s.e. of the tilt) from tiltlaw.txt"""
    out = {}
    for ln in open(TH + "/energetics/tiltlaw.txt").read().splitlines()[1:]:
        c = [s.strip() for s in ln.split("|")]
        k, lam = int(c[0].split()[0]), float(c[0].split()[1])
        nw = float(c[2].split("±")[0]); mt, se = c[4].split("/")[1].split("±")
        out[(k, lam)] = (nw, float(mt), float(se))
    return out


def theta(x, y, x0, h, a):
    z, z0 = x + 1j * y, x0 + 1j * h
    return 0.5 * np.angle(np.sin(np.pi * (z - z0) / a)) + 0.5 * np.angle(np.sin(np.pi * (z - np.conj(z0)) / a)) + np.pi * x / a + np.pi / 2


def schematic(ax, lam=2.0):
    a, h = 1.0, 0.22; thL = 0.5 * np.arccos(1 / lam); x0 = thL * a / np.pi; tau = 2 * thL
    Y1 = 1.0; xs = np.linspace(-0.95, 0.95, 17); ys = np.linspace(0.08, Y1 - 0.04, 10); L = 0.045
    for x in xs:
        for y in ys:
            if np.hypot(x - x0, y - h) < 0.07 or np.hypot(x - x0 + a, y - h) < 0.07: continue
            t = theta(x, y, x0, h, a); ax.plot([x - L * np.cos(t), x + L * np.cos(t)], [y - L * np.sin(t), y + L * np.sin(t)], color="0.45", lw=0.8, solid_capstyle="round")
    # wall and its anchoring (fixed, winds by pi per period)
    ax.fill_between([-1.05, 1.05], -0.07, 0, color="0.82", lw=0); ax.plot([-1.05, 1.05], [0, 0], color="k", lw=1.2)
    for x in np.linspace(-1.0, 1.0, 21):
        t = np.pi / 2 + np.pi * x / a; ax.plot([x - L * np.cos(t), x + L * np.cos(t)], [-0.035 - L * np.sin(t), -0.035 + L * np.sin(t)], color="k", lw=0.9)
    # homeotropic sites (x = n a) and the defects displaced by D
    for n in (0, -1):
        ax.plot(n * a, 0.0, "o", ms=4, mfc="w", mec=RED, mew=1.0, zorder=5)
    psi = np.pi / 2 + tau                       # comet axis: normal (pi/2) turned by tau = 2 kappa D
    for xd in (x0, x0 - a):
        ax.plot([xd, xd + 0.13 * np.cos(psi)], [h, h + 0.13 * np.sin(psi)], color=RED, lw=2.0, solid_capstyle="round", zorder=6)
        ax.plot(xd, h, "o", ms=4.5, color=RED, mec="w", mew=0.5, zorder=7)
    # tilt angle tau between the wall normal and the comet axis (right defect)
    ax.plot([x0, x0], [h, h + 0.2], color=DARK, lw=0.6, ls=":", zorder=6)
    ax.add_patch(Arc((x0, h), 0.34, 0.34, theta1=90, theta2=90 + np.degrees(tau), color=DARK, lw=0.7, zorder=6))
    ax.text(x0 + 0.03, h + 0.24, r"$\tau=2\theta_L$", fontsize=7, ha="center", va="bottom", color=DARK, bbox=dict(fc="w", ec="none", pad=0.4), zorder=8)
    # displacement D from the homeotropic site
    ax.add_patch(FancyArrowPatch((0, 0.06), (x0, 0.06), arrowstyle="-|>", mutation_scale=6, color=DARK, lw=0.7))
    ax.text(x0 / 2, 0.075, r"$D$", fontsize=7, ha="center", va="bottom", color=DARK, bbox=dict(fc="w", ec="none", pad=0.3), zorder=8)
    # far field at the Leslie angle; shear flow arrows on the right
    ax.plot([0.25, 0.25 + 0.22 * np.cos(thL)], [Y1 + 0.06, Y1 + 0.06 + 0.22 * np.sin(thL)], color=BLUE, lw=1.6, solid_capstyle="round")
    ax.plot([0.25, 0.52], [Y1 + 0.06, Y1 + 0.06], color=BLUE, lw=0.6, ls=":")
    ax.add_patch(Arc((0.25, Y1 + 0.06), 0.36, 0.36, theta1=0, theta2=np.degrees(thL), color=BLUE, lw=0.7))
    ax.text(0.56, Y1 + 0.09, r"$\theta_L$", fontsize=7, color=BLUE, va="center")
    ax.text(-0.95, Y1 + 0.08, "Far field at the\nLeslie angle", fontsize=6.5, color=BLUE, va="center")
    for y in (0.3, 0.55, 0.8):
        ax.add_patch(FancyArrowPatch((1.08, y), (1.08 + 0.32 * y, y), arrowstyle="-|>", mutation_scale=5, color=GREY, lw=0.8))
    ax.plot([1.08, 1.08], [0.0, 0.95], color=GREY, lw=0.6)
    ax.text(1.12, 0.98, "Shear", fontsize=6.5, color=GREY, va="bottom")
    ax.text(0.0, -0.1, "No-slip wall, winding anchoring", fontsize=6.5, ha="center", va="top")
    ax.set_xlim(-1.05, 1.45); ax.set_ylim(-0.22, Y1 + 0.22); ax.set_aspect("equal"); ax.axis("off")


def tiltlaw(ax):
    T = read_tiltlaw(); lx = np.linspace(1.0, 6.4, 300)
    ax.plot(lx, np.degrees(np.arccos(1 / lx)), "k-", lw=1.0, zorder=1)
    ax.text(2.9, 88, r"$\arccos(1/\lambda)$", fontsize=6.5, ha="left", va="center")
    mk = {4: "s", 6: "o", 8: "^"}
    for k in (4, 6, 8):
        ls = sorted(l for (k_, l) in T if k_ == k)
        ax.errorbar(ls, [T[(k, l)][1] for l in ls], [T[(k, l)][2] for l in ls], fmt=mk[k] + "-", color=KCOL[k], ms=3.2, lw=0.9, capsize=1.5, label=f"$k={k}$", zorder=3)
    ax.set_xlabel(r"Flow-alignment parameter $\lambda$"); ax.set_ylabel(r"Wall-defect tilt $|\tau|$ (deg)")
    ax.set_xlim(0.6, 6.4); ax.set_ylim(0, 94); ax.set_yticks([0, 30, 60, 90])
    ax.legend(loc="lower left", bbox_to_anchor=(0.14, 0.06), fontsize=6.3, handlelength=1.4, title="Active MPCD", title_fontsize=6.3, alignment="left")
    ins = ax.inset_axes([0.64, 0.11, 0.35, 0.26])
    for k in (4, 6, 8):
        ls = sorted(l for (k_, l) in T if k_ == k)
        ins.plot(ls, [T[(k, l)][0] - (2 * k - 2) for l in ls], mk[k] + "-", color=KCOL[k], ms=2.2, lw=0.8)
    ins.axhline(0, color=GREY, lw=0.5); ins.set_xlim(0.6, 6.4); ins.set_ylim(-2.8, 0.4); ins.set_yticks([-2, -1, 0]); ins.set_xticks([2, 4, 6])
    ins.tick_params(labelsize=5.5, length=1.5, pad=1); ins.set_xlabel(r"$\lambda$", fontsize=6, labelpad=0)
    ins.set_title("Wall $+\\frac{1}{2}$ count\n$-\\,(2k-2)$", fontsize=5.5, pad=2); ins.patch.set_alpha(0.0)


def variants(ax):
    R = json.load(open(TH + "/shear_tilt/scan.json"))
    VN = {"noslip": "Measured no-slip flow", "wallonly": "Shear only in the wall band", "slip": "No shear in the wall band"}
    VC = {"noslip": RED, "wallonly": "#d08a1e", "slip": BLUE}
    for gi, l in enumerate((2.0, 4.0)):
        x0 = gi * 4.0
        for vi, v in enumerate(VN):
            r = [x for x in R if x["lam"] == l and x["U"] == 1.0 and x["variant"] == v and x["ring"] == 12][0]
            ax.bar(x0 + vi, r["tilt"], width=0.82, color=VC[v], zorder=2)
            ax.text(x0 + vi, r["tilt"] + 2, f"{r['tilt']:.0f}°", ha="center", va="bottom", fontsize=6)
        ref = np.degrees(np.arccos(1 / l)); ax.plot([x0 - 0.55, x0 + 2.55], [ref, ref], "k--", lw=0.8, zorder=3)
    ax.set_xticks([1, 5]); ax.set_xticklabels([r"$\lambda=2$", r"$\lambda=4$"]); ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.7, 6.7); ax.set_ylim(0, 150); ax.set_yticks([0, 30, 60, 90, 120]); ax.set_ylabel("Tilt of the wall layer (deg)")
    h = [Patch(color=VC[v], label=VN[v]) for v in VN] + [plt.Line2D([], [], color="k", ls="--", lw=0.8, label=r"$\arccos(1/\lambda)$")]
    ax.legend(handles=h, loc="upper left", fontsize=6, ncol=1, handlelength=1.2, borderaxespad=0.1)
    ax.text(1.0, 1.0, "Passive Q-tensor cavity,\nimposed vortex", transform=ax.transAxes, fontsize=6, ha="right", va="bottom", color=DARK)


if __name__ == "__main__":
    fig = plt.figure(figsize=(PAGE, 2.25))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.9], wspace=0.38, left=0.0, right=1.0, bottom=0.17, top=0.93)
    a0, a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    schematic(a0); tiltlaw(a1); variants(a2)
    label(a1, "b"); label(a2, "c")
    fig.canvas.draw(); y_lab = a1.get_position().y1 + 0.02 * a1.get_position().height
    fig.text(a0.get_position().x0 + 0.01, y_lab, "(a)", fontsize=9, fontweight="bold", va="bottom", ha="left")
    save(fig, "fig3_tiltlaw")
