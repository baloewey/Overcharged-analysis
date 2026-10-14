#!/usr/bin/env python3
"""Shielding by the wall layer, positive vs mirror (negative) winding (user, 2026-10-14). Potential felt by one bulk +1/2 from the wall layer alone (layer charges,
same-sign wall images, boundary winding; one-constant K = 32.5 kBT), from shield_wl/probe_maps.wl (k = +6: wall +1/2) and mirror_maps.wl (k = -6: wall -1/2, which
attract the self-propelled bulk +1/2). Each panel: V minus its minimum over the interior (inside the 0.8 level), own colour scale (diverging: blue = deeper than any
interior point, i.e. the wells at the wall sites), contours every 10 kBT, dashed: the 0.8 level. Right column: range of the layer's field (probe's own wall term
removed) inside the 0.6 and 0.8 levels for ideal and measured layers.
Data: k = +6 fine-output runs (fineflow, area 2500); k = -6 seeded runs (neglayer, area 2500, active). -> fig_shielding_mirror.png (+pdf), shielding_mirror.txt"""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "shield_wl"); PIK = np.pi * 32.5
J = json.load(open(os.path.join(D, "layers.json"))); A2, B2, R1 = J["A2"], J["B2"], J["R1"]


def load(name):
    a = np.loadtxt(os.path.join(D, f"maps_{name}.csv"), delimiter=","); return [(a[a[:, 0] == l, 1], a[a[:, 0] == l, 2], a[a[:, 0] == l, 3] * PIK, (a[a[:, 0] == l, 3] - a[a[:, 0] == l, 4]) * PIK) for l in np.unique(a[:, 0])]


lev = lambda x, y, shape: np.hypot(x, y) / R1 if shape == "circle" else np.hypot(x / A2, y / B2)


def ranges(M, shape):
    return np.array([[np.ptp(vl[lev(x, y, shape) < 0.6]), np.ptp(vl[lev(x, y, shape) < 0.8])] for x, y, v, vl in M])


def grid(x, y, v):
    xs, ys = np.unique(x), np.unique(y); G = np.full((len(xs), len(ys)), np.nan); G[np.searchsorted(xs, x), np.searchsorted(ys, y)] = v; return xs, ys, G


CASES = [("circle_meas", "circle", "circle_ideal"), ("ellipse_ideal", "ellipse", None), ("ellipse_meas", "ellipse", "ellipse_ideal"),
         ("mcircle_meas", "circle", "mcircle_ideal"), ("mellipse_ideal", "ellipse", None), ("mellipse_meas", "ellipse", "mellipse_ideal")]
if __name__ == "__main__":
    M = {k: load(k) for k in ("circle_ideal", "circle_meas", "ellipse_ideal", "ellipse_meas", "mcircle_ideal", "mcircle_meas", "mellipse_ideal", "mellipse_meas")}
    shape_of = lambda k: "circle" if "circle" in k else "ellipse"
    R = {k: ranges(M[k], shape_of(k)) for k in M}
    L = ["SHIELDING, positive vs mirror winding: range of the layer's field felt by a bulk +1/2 (kBT), inside the 0.6 / 0.8 levels; median [min-max] over measured frames.",
         "case | layers | wall defects | inside 0.6 | inside 0.8"]
    for k in M:
        nw = [len(l) for l in J[k]]
        L.append(f"{k:15s} | {len(R[k]):2d} | {min(nw)}-{max(nw)} | {np.median(R[k][:,0]):6.1f} [{R[k][:,0].min():.1f}-{R[k][:,0].max():.1f}] | {np.median(R[k][:,1]):6.1f} [{R[k][:,1].min():.1f}-{R[k][:,1].max():.1f}]")
    fig, axs = plt.subplots(2, 4, figsize=(19, 8.2), gridspec_kw=dict(width_ratios=[1, 1.35, 1.35, 1.1])); th = np.linspace(0, 2 * np.pi, 400)
    titles = ["(a) $k=+6$, circle, measured layer", "(b) $k=+6$, ellipse 2, ideal layer (10)", "(c) $k=+6$, ellipse 2, measured layer",
              "(e) $k=-6$, circle, measured layer", "(f) $k=-6$, ellipse 2, ideal layer (14)", "(g) $k=-6$, ellipse 2, measured layer"]
    for i, (k, shape, _) in enumerate(CASES):
        ax = axs[i // 3, i % 3]; rr = R[k][:, 0]; j = int(np.argsort(rr)[len(rr) // 2]); x, y, v, vl = M[k][j]
        inner = lev(x, y, shape) < 0.8; v = v - np.nanmin(v[inner]); xs, ys, G = grid(x, y, v)
        hi = max(np.nanpercentile(v[inner], 98), 5); lo = min(np.nanmin(v), -1)
        im = ax.pcolormesh(xs, ys, G.T, cmap="RdBu_r", shading="nearest", norm=TwoSlopeNorm(0, max(lo, -hi), hi))
        ax.contour(xs, ys, G.T, levels=np.arange(10, hi, 10), colors="k", linewidths=0.35, alpha=0.6)
        A_, B_ = (R1, R1) if shape == "circle" else (A2, B2)
        ax.plot(A_ * np.cos(th), B_ * np.sin(th), "k-", lw=1); ax.plot(0.8 * A_ * np.cos(th), 0.8 * B_ * np.sin(th), "k--", lw=0.6)
        z = np.array(J[k][j]); ax.plot(z[:, 0], z[:, 1], "o", ms=4, mec="w", mew=0.5, color=("#1f5fa8" if z[0, 2] > 0 else "#e08a3c"))
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(titles[i] + (f" ({len(z)})" if "meas" in k else ""), fontsize=9.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, shrink=0.8); cb.set_label("$V-V_\\mathrm{min,interior}$ ($k_BT$)", fontsize=8)
    for row, (lab, keys) in enumerate((("$k=+6$ (wall $+1/2$)", ["circle_ideal", "circle_meas", "ellipse_ideal", "ellipse_meas"]),
                                        ("$k=-6$ (wall $-1/2$)", ["mcircle_ideal", "mcircle_meas", "mellipse_ideal", "mellipse_meas"]))):
        ax = axs[row, 3]
        for jj, kk in enumerate(keys):
            ax.bar(jj - 0.18, np.median(R[kk][:, 0]), 0.36, color="#1f5fa8", label="inside 0.6 level" if jj == 0 else None)
            ax.bar(jj + 0.18, np.median(R[kk][:, 1]), 0.36, color="#c2452d", label="inside 0.8 level" if jj == 0 else None)
            if len(R[kk]) > 1: ax.plot([jj - 0.18] * len(R[kk]), R[kk][:, 0], ".", color="0.2", ms=3); ax.plot([jj + 0.18] * len(R[kk]), R[kk][:, 1], ".", color="0.2", ms=3)
        ax.set_yscale("log"); ax.set_ylim(0.5, 1000); ax.set_xticks(range(4)); ax.set_xticklabels(["circle\nideal", "circle\nmeas.", "ellipse\nideal", "ellipse\nmeas."], fontsize=8.5)
        ax.set_ylabel("range of the layer's field ($k_BT$)"); ax.set_title(f"({'dh'[row]}) Leftover field, {lab}", fontsize=9.5); ax.legend(fontsize=7.5, frameon=False, loc="upper left")
        [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_shielding_mirror.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "shielding_mirror.txt"), "w").write("\n".join(L) + "\n")
