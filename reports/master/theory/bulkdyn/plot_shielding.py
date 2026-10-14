#!/usr/bin/env python3
"""Figure: elastic shielding by the wall layer (Japan-poster construction) and the measured potential of mean force (user, 2026-10-14).
Top: potential V(z) felt by one bulk +1/2 from the wall layer alone (layer charges, their wall images and the boundary winding; one-constant, pi K units -> kBT
with K = 32.5 kBT), from shield_wl/probe_maps.wl: (a) circle with a measured layer (3D surface, as in the poster), (b) aspect-2 ellipse with an idealised full layer,
(c) aspect-2 ellipse with a measured (depleted) layer. Bottom: (d, e) potential of mean force -ln(rho_+/<rho_+>) of +1/2 from the simulations (shielding.pmf),
(f) shielding measure: range of V inside the 0.6 and 0.8 levels (fraction of the way to the wall), per layer.
-> fig_shielding_layer.png (+pdf), shielding_layer.txt"""
import os, sys, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from shielding import pmf, RUNS  # noqa
PIK = np.pi * 32.5
J = json.load(open(os.path.join(H, "shield_wl", "layers.json"))); A2, B2, R1 = J["A2"], J["B2"], J["R1"]


def load(name):
    a = np.loadtxt(os.path.join(H, "shield_wl", f"maps_{name}.csv"), delimiter=","); out = []
    for l in np.unique(a[:, 0]):
        b = a[a[:, 0] == l]; out.append((b[:, 1], b[:, 2], b[:, 3] * PIK, (b[:, 3] - b[:, 4]) * PIK))   # total V, and the layer part (V minus the probe's own wall term)
    return out


def level(x, y, shape):
    return np.hypot(x, y) / R1 if shape == "circle" else np.hypot(x / A2, y / B2)


def ranges(maps, shape):
    r6, r8 = [], []
    for x, y, v, vl in maps:
        lv = level(x, y, shape); r6.append(np.ptp(vl[lv < 0.6])); r8.append(np.ptp(vl[lv < 0.8]))
    return np.array(r6), np.array(r8)


def grid(x, y, v):
    xs, ys = np.unique(x), np.unique(y); G = np.full((len(xs), len(ys)), np.nan)
    G[np.searchsorted(xs, x), np.searchsorted(ys, y)] = v; return xs, ys, G


if __name__ == "__main__":
    M = {k: load(k) for k in ("circle_ideal", "circle_meas", "ellipse_ideal", "ellipse_meas")}
    L = ["SHIELDING by the wall layer: potential of one bulk +1/2 from the layer alone (charges + images + boundary winding), one-constant K = 32.5 kBT.",
         "case | layers | wall defects | range of the LAYER part of V (V minus the probe's own wall term) inside the 0.6 level (kBT): median [min-max] | inside the 0.8 level"]
    R = {}
    for k, shape in (("circle_ideal", "circle"), ("circle_meas", "circle"), ("ellipse_ideal", "ellipse"), ("ellipse_meas", "ellipse")):
        r6, r8 = ranges(M[k], shape); R[k] = (r6, r8); nw = [len(l) for l in J[k]]
        L.append(f"{k:14s} | {len(r6):2d} | {min(nw)}-{max(nw)} | {np.median(r6):6.1f} [{r6.min():.1f}-{r6.max():.1f}] | {np.median(r8):6.1f} [{r8.min():.1f}-{r8.max():.1f}]")
    fig = plt.figure(figsize=(16, 9.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.45)
    vmax = 100
    # (a) circle, measured layer: 3D surface (poster style)
    ax = fig.add_subplot(gs[0, 0], projection="3d"); x, y, v, _ = M["circle_meas"][0]; xs, ys, G = grid(x, y, v - np.nanmin(v))
    XX, YY = np.meshgrid(xs, ys, indexing="ij"); ax.plot_surface(XX, YY, np.minimum(G, vmax), cmap="copper", rstride=1, cstride=1, linewidth=0, antialiased=False)
    ax.set_zlim(0, vmax); ax.set_xticks([]); ax.set_yticks([]); ax.set_zlabel("$V$ ($k_BT$)"); ax.view_init(28, -60)
    ax.set_title(f"(a) Circle, measured layer ({len(J['circle_meas'][0])} wall $+1/2$)", fontsize=10)
    th = np.linspace(0, 2 * np.pi, 400)
    def flat(ax, maps, lay, A, B, title, k=0):
        x, y, v, _ = maps[k]; xs, ys, G = grid(x, y, v - np.nanmin(v))
        from matplotlib.colors import PowerNorm
        im = ax.pcolormesh(xs, ys, np.minimum(G, 300).T, cmap="magma_r", shading="nearest", norm=PowerNorm(0.5, vmin=0, vmax=300))
        ax.contour(xs, ys, G.T, levels=np.arange(20, 300, 20), colors="w", linewidths=0.4, alpha=0.7)
        ax.plot(A * np.cos(th), B * np.sin(th), "k-", lw=1); z = np.array(lay[k]); ax.plot(z[:, 0], z[:, 1], "o", ms=5, color="#1f5fa8", mec="w", mew=0.6)
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=10); return im
    flat(fig.add_subplot(gs[0, 1]), M["ellipse_ideal"], J["ellipse_ideal"], A2, B2, "(b) Ellipse 2, idealised full layer (10)")
    k = int(np.argsort(R["ellipse_meas"][0])[len(R["ellipse_meas"][0]) // 2])
    im = flat(fig.add_subplot(gs[0, 2]), M["ellipse_meas"], J["ellipse_meas"], A2, B2, f"(c) Ellipse 2, measured layer ({len(J['ellipse_meas'][k])} wall $+1/2$)", k)
    cb = fig.colorbar(im, ax=fig.axes[1:3], fraction=0.025, pad=0.01, location="bottom"); cb.set_label("potential of one bulk $+1/2$ from the layer, $V$ ($k_BT$; square-root scale, contours every 20)")
    for i, (key, ttl) in enumerate((("circle", "(d) Circle: measured $-\\ln\\rho_+$"), ("ellipse", "(e) Ellipse 2: measured $-\\ln\\rho_+$"))):
        Xp, Yp, Up, A_, B_ = pmf(RUNS[key], key); ax = fig.add_subplot(gs[1, i]); im2 = ax.pcolormesh(Xp, Yp, np.clip(Up, -3, 3), cmap="RdBu_r", shading="nearest", vmin=-3, vmax=3)
        ax.plot(A_ * np.cos(th), B_ * np.sin(th), "k-", lw=1); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(ttl, fontsize=10)
    cb2 = fig.colorbar(im2, ax=fig.axes[-2:], fraction=0.03, pad=0.01, location="bottom"); cb2.set_label(r"potential of mean force $-\ln(\rho_+/\langle\rho_+\rangle)$ ($k_BT_\mathrm{eff}$)")
    ax = fig.add_subplot(gs[1, 2]); labs = ["circle\nideal", "circle\nmeasured", "ellipse 2\nideal", "ellipse 2\nmeasured"]
    for j, kk in enumerate(("circle_ideal", "circle_meas", "ellipse_ideal", "ellipse_meas")):
        r6, r8 = R[kk]; ax.bar(j - 0.18, np.median(r6), 0.36, color="#1f5fa8", label="inside 0.6 level" if j == 0 else None)
        ax.bar(j + 0.18, np.median(r8), 0.36, color="#c2452d", label="inside 0.8 level" if j == 0 else None)
        if len(r6) > 1: ax.plot([j - 0.18] * len(r6), r6, ".", color="0.2", ms=3); ax.plot([j + 0.18] * len(r8), r8, ".", color="0.2", ms=3)
    ax.set_yscale("log"); ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=9); ax.set_ylabel("range of the layer's field ($k_BT$)")
    ax.set_title("(f) Leftover field: how well the layer shields", fontsize=10); ax.legend(fontsize=8, frameon=False); [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.savefig(os.path.join(H, "fig_shielding_layer.png"), dpi=200, bbox_inches="tight")
    print("\n".join(L)); open(os.path.join(H, "shielding_layer.txt"), "w").write("\n".join(L) + "\n")
