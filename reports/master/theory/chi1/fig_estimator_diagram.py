#!/usr/bin/env python3
"""Diagram of the band estimator (user, 2026-11-24). Exact one-constant static texture of a tilted 2k-2 ring + central +1 (ejection/tilted_structures.fields;
anchoring theta = k phi), circle R 31.83, layer tilt 60 deg, defects at the depth-law radius, shown UNROLLED (x = arclength along the wall, y = depth;
director drawn at its angle from the local tangent = flow direction; a uniform zero-mode field would be sticks at tau/2 = 30 deg).
 (a) sparse (k 4), (b) crowded (k 16): shaded = near-field zone y < d + a/2pi (field varies strongly along the wall); dashed = band window y 2-3 cells.
 (c) bias of the ring-average estimator on REAL MPCD frames: static texture of the measured defects ring-averaged over lattice cells in 1-cell depth bins (as
     static_band.py), minus the static tilt; median over coherent frames. A perfectly regular ring gives exactly zero at every depth (checked numerically).
-> fig_estimator_diagram.pdf/.png"""
import glob, json, os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import importlib.util
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"; MAIN = 'if __name__ == ' + '"__main__":'
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split("\n" + MAIN)[0], "ts", "exec"), TS.__dict__)
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); R = 31.830989; TAU = np.radians(60)
RED, BLUE, GREY = "#c0392b", "#1f5fa8", "#8a9099"


def ring(k):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + TAU / N))) + [0j], [0.5] * N + [1.0], N, rs


def strip(ax, k, title, win=64.0):
    zs, qs, N, rs = ring(k); d = R - rs; a = 2 * np.pi * rs / N
    s = np.arange(-win / 2, win / 2 + 0.01, 1.0); yy = np.arange(0.5, 10.01, 0.75); S, YY = np.meshgrid(s, yy); ph = S / R
    r = R - YY; th2, _, _ = TS.fields(r * np.cos(ph), r * np.sin(ph), zs, qs, R); al = 0.5 * np.angle(np.exp(1j * (th2 - 2 * ph - np.pi)))
    ax.axhspan(0, d + a / (2 * np.pi), color="#f6c9a8", alpha=0.55, lw=0, zorder=0)
    ax.axhline(2, color=BLUE, ls="--", lw=1.8, zorder=1); ax.axhline(3, color=BLUE, ls="--", lw=1.8, zorder=1)
    L = 0.42
    ax.plot(np.array([S - L * np.cos(al), S + L * np.cos(al)]).reshape(2, -1), np.array([YY - L * np.sin(al), YY + L * np.sin(al)]).reshape(2, -1),
            "-", color=GREY, lw=1.6, solid_capstyle="round", zorder=2)          # depth axis points inwards: alpha > 0 slopes down-right
    for z in zs[:N]:
        for shift in (0, -2 * np.pi * R, 2 * np.pi * R):
            sj = R * np.angle(z) + shift
            if -win / 2 - 2 < sj < win / 2 + 2:
                yj = R - abs(z)
                ax.plot([sj, sj + 1.6 * np.sin(TAU)], [yj, yj - 1.6 * np.cos(TAU)], "-", color=RED, lw=4, solid_capstyle="round", zorder=4)
                ax.plot(sj, yj, "o", color=RED, ms=9, mec="white", mew=1.3, zorder=5)
    ax.fill_between([-win / 2, win / 2], -1.2, 0, color="0.8", zorder=3); ax.axhline(0, color="k", lw=3, zorder=3)
    ax.text(win / 2 - 0.5, 9.9, f"Shaded: near-field zone, depth < d + a/2π = {d + a / (2 * np.pi):.1f} cells", ha="right", va="bottom", fontsize=11.5, color="#a0522d", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    ax.text(-win / 2 + 0.5, 9.9, "Dashed: band window (2–3 cells)", ha="left", va="bottom", fontsize=11.5, color=BLUE, zorder=6, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    ax.set_xlim(-win / 2, win / 2); ax.set_ylim(10.4, -1.2); ax.set_aspect("equal"); ax.set_xlabel("Distance along the wall (cells)"); ax.set_ylabel("Depth (cells)")
    ax.set_title(f"{title}: k = {k}, defect spacing a = {a:.1f} cells", fontsize=13.5, loc="left")


def real_bias(k, nmax=40):
    runs = sorted(os.path.dirname(p) for p in glob.glob(f"/home/baloewe/Physics/Overcharged/runs/cavity/server_main/main_act0.07/rho*/k{k:02d}/s0[1-3]/directorfield.dat.gz")
                  if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - R) < 0.5)
    bc = json.load(open(runs[0] + "/input.json"))["BC"][0]; Q = np.array(bc["Q"][:2], float); t, fr = read_with_angle(runs[0], Q, 0.5)
    n = 67; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(n) + 0.5 - Q[1], indexing="ij"); ph = np.arctan2(Y, X); u = R - np.hypot(X, Y)
    ys = np.arange(0.0, 10.0, 1.0); out = []; N = 2 * k - 2
    for f in fr[::5]:
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        if w.sum() != N or abs(f[:, 2].sum() - k) > 1e-6: continue
        ez = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))
        if abs(ez) < 0.9: continue
        sg = np.sign(np.angle(ez)) or 1.0; zs = list(f[:, 0] + 1j * f[:, 1]); qs = list(f[:, 2]); th2s, _, _ = TS.fields(X, Y, zs, qs, R)
        nt = TS.natural_tilts(zs, qs, R, len(zs)); ts_ = sg * np.angle(np.mean(np.exp(1j * nt[np.flatnonzero(w)])))
        out.append([np.angle(np.sum(np.exp(1j * sg * (th2s[(u >= y) & (u < y + 1)] - 2 * ph[(u >= y) & (u < y + 1)] - np.pi))) * np.exp(-1j * ts_)) for y in ys])
        if len(out) >= nmax: break
    zr, qr, Nn, rs = ring(k)
    return ys + 0.5, np.degrees(np.median(np.array(out), 0)), R - rs + rs / Nn, len(out)


def bias_curve(ax):
    for k, c, lab in ((4, "#f6a04d", "sparse"), (16, "#4a1c0c", "crowded")):
        y, b, edge, nf = real_bias(k); ax.plot(y, b, "o-", color=c, lw=2.4, ms=6, label=f"k = {k} ({lab}): real frames (n = {nf})")
        ax.axvline(edge, color=c, ls=":", lw=1.6)
    ax.axhline(0, color="0.4", lw=2, ls="--", label="Perfectly regular ring: exactly 0")
    ax.axvspan(2, 3, color=BLUE, alpha=0.12, lw=0); ax.text(2.5, 9.0, "band\nwindow", ha="center", va="top", fontsize=11, color=BLUE)
    ax.set_xlabel("Depth of the ring (cells)"); ax.set_ylabel("Bias of the ring average (deg)")
    ax.set_title("(c) Bias from the defects' own field (dotted: d + a/2π)", fontsize=13.5, loc="left"); ax.legend(fontsize=10.5, frameon=False, loc="lower right")
    ax.set_ylim(-12, 12); ax.set_xlim(0, 10)


if __name__ == "__main__":
    fig = plt.figure(figsize=(14, 13)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], hspace=0.42, wspace=0.25)
    strip(fig.add_subplot(gs[0, :]), 4, "(a) Sparse layer"); strip(fig.add_subplot(gs[1, :]), 16, "(b) Crowded layer"); bias_curve(fig.add_subplot(gs[2, 0]))
    ax = fig.add_subplot(gs[2, 1]); ax.axis("off")
    ax.text(0.0, 1.0, "How the band angle is measured\n\n"
            "1. Ring average: every cell at depth 2–3 cells,\n    complex mean of the director angle from the\n    flow direction (horizontal in a, b).\n\n"
            "2. A perfectly regular ring: the defects' own field\n    averages out exactly, at any depth.\n\n"
            "3. Real layers are irregular (spacings, depths)\n    and sampled on lattice cells: inside the\n    near-field zone the average picks up a few\n    degrees of bias, more for sparse walls.\n\n"
            "4. Static subtraction: static texture of the SAME\n    defects, averaged on the SAME cells, removes it:\n"
            r"    $2\alpha_{\rm corr} = 2\alpha_{\rm meas} - 2\alpha_{\rm static} + \tau_{\rm static}$",
            va="top", ha="left", fontsize=12.5, color="0.15", transform=ax.transAxes)
    fig.text(0.02, 0.0, "(a, b): static one-constant texture of a regular layer tilted by τ = 60°, unrolled (horizontal = flow direction). A uniform field would be sticks at τ/2 = 30° everywhere.",
             fontsize=11, color="0.3", va="top")
    fig.savefig(os.path.join(H, "fig_estimator_diagram.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(H, "fig_estimator_diagram.png"), dpi=150, bbox_inches="tight")
