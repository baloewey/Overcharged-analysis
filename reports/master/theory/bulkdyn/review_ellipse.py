#!/usr/bin/env python3
"""Ellipse shielding review (TO REVIEW, RESEARCH.md 2026-10-14): (i) sensitivity to the winding coefficient of the uniform-flux term (charges present vs the imposed
winding k); (ii) are the tips better shielded? Leftover field of the LAYER term (probe's own wall term removed) in tip regions (|x| > 0.6 A) vs the centre/flank region
(|x| < 0.6 A), both inside the 0.8 level; layer spacing at the tips vs flanks in the measured layers; (iii) maps to the 0.98 level, layer term only.
Input: shield_wl/review_*.csv (ellipse_review.wl). -> fig_review_ellipse.png (+pdf), review_ellipse.txt"""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "shield_wl"); PIK = np.pi * 32.5
J = json.load(open(os.path.join(D, "layers.json"))); A, B = J["A2"], J["B2"]
th = np.linspace(0, 2 * np.pi, 20001); Wl = A * np.cos(th) + 1j * B * np.sin(th); sarc = np.r_[0, np.cumsum(np.abs(np.diff(Wl)))]; P = sarc[-1]


def load(name):
    a = np.loadtxt(os.path.join(D, f"review_{name}.csv"), delimiter=","); out = []
    for l in np.unique(a[:, 0]):
        b = a[a[:, 0] == l]; out.append(dict(x=b[:, 1], y=b[:, 2], lay=(b[:, 3] - b[:, 4]) * PIK, layk=(b[:, 5] - b[:, 4]) * PIK))
    return out


def regions(m, key):
    lv = np.hypot(m["x"] / A, m["y"] / B); v = m[key]; inn = lv < 0.8; tip = inn & (np.abs(m["x"]) > 0.6 * A); mid = inn & (np.abs(m["x"]) <= 0.6 * A)
    return np.ptp(v[lv < 0.6]), np.ptp(v[tip & (m["x"] > 0)]), np.ptp(v[tip & (m["x"] < 0)]), np.ptp(v[mid])


if __name__ == "__main__":
    L = ["ELLIPSE SHIELDING REVIEW (aspect 2, area 2500; layer term = field of the layer felt by a bulk +1/2, probe's own wall term removed; kBT).",
         "case | flux coefficient | range inside 0.6 level | right tip region | left tip region | centre/flank region (|x| < 0.6 A), all inside the 0.8 level  [medians over frames]"]
    M = {k: load(k) for k in ("ellipse_ideal", "ellipse_meas", "mellipse_ideal", "mellipse_meas")}
    for k, m in M.items():
        for key, lab in (("lay", "charges present"), ("layk", "imposed winding")):
            R = np.array([regions(x, key) for x in m]); md = np.median(R, 0)
            L.append(f"{k:15s} | {lab:16s} | {md[0]:6.1f} | {md[1]:6.1f} | {md[2]:6.1f} | {md[3]:6.1f}")
    # layer spacing at tips vs flanks in measured layers (arclength gaps in units of the ideal spacing)
    for k, nideal in (("ellipse_meas", 10), ("mellipse_meas", 14)):
        tipg, midg, tipocc, midocc = [], [], [], []
        for lay in J[k]:
            z = np.array([complex(a, b) for a, b, q in lay]); sa = np.sort([sarc[np.argmin(np.abs(Wl - x))] for x in z]); g = np.diff(np.r_[sa, sa[0] + P]) / (P / nideal)
            mids = (sa + np.diff(np.r_[sa, sa[0] + P]) / 2) % P; xm = np.array([Wl[np.argmin(np.abs(sarc - s_))].real for s_ in mids])
            tipg += list(g[np.abs(xm) > 0.6 * A]); midg += list(g[np.abs(xm) <= 0.6 * A])
            xz = np.array([x.real for x in z]); tipocc.append(np.sum(np.abs(xz) > 0.6 * A)); midocc.append(np.sum(np.abs(xz) <= 0.6 * A))
        ftip = np.mean(np.abs(Wl.real) > 0.6 * A)   # fraction of the arclength in the tip regions (by parameter count ~ arclength-weighted below)
        arc_tip = np.sum(np.abs(np.diff(Wl))[np.abs(Wl[:-1].real) > 0.6 * A]) / P
        L.append(f"{k}: gaps (ideal spacing = 1) median tips {np.median(tipg):.2f}, centre/flanks {np.median(midg):.2f}; mean wall defects in tip regions {np.mean(tipocc):.2f} "
                 f"(ideal {nideal*arc_tip:.2f}), centre/flanks {np.mean(midocc):.2f} (ideal {nideal*(1-arc_tip):.2f})")
    fig, axs = plt.subplots(2, 2, figsize=(14, 6.4)); tt = np.linspace(0, 2 * np.pi, 300)
    for ax, (k, ttl) in zip(axs.flat, (("ellipse_ideal", "(a) $k=+6$, ideal layer (10)"), ("ellipse_meas", "(b) $k=+6$, measured layer"),
                                        ("mellipse_ideal", "(c) $k=-6$, ideal layer (14)"), ("mellipse_meas", "(d) $k=-6$, measured layer"))):
        m = M[k]; j = int(np.argsort([regions(x, "lay")[0] for x in m])[len(m) // 2]); x = m[j]; lv = np.hypot(x["x"] / A, x["y"] / B)
        v = x["lay"] - np.min(x["lay"][lv < 0.8]); xs, ys = np.unique(x["x"]), np.unique(x["y"]); G = np.full((len(xs), len(ys)), np.nan)
        G[np.searchsorted(xs, x["x"]), np.searchsorted(ys, x["y"])] = v; hi = np.nanpercentile(np.abs(v), 97)
        im = ax.pcolormesh(xs, ys, G.T, cmap="RdBu_r", shading="nearest", norm=TwoSlopeNorm(0, -hi, hi))
        ax.contour(xs, ys, G.T, levels=np.arange(10, hi, 10), colors="k", linewidths=0.3, alpha=0.5)
        ax.plot(A * np.cos(tt), B * np.sin(tt), "k-", lw=1); ax.plot(0.8 * A * np.cos(tt), 0.8 * B * np.sin(tt), "k--", lw=0.5)
        for xx in (-0.6 * A, 0.6 * A): ax.plot([xx, xx], [-B, B], color="0.5", lw=0.6, ls=":")
        z = np.array(J[k][j]); ax.plot(z[:, 0], z[:, 1], "o", ms=4, mec="w", mew=0.5, color="#1f5fa8" if z[0, 2] > 0 else "#e08a3c")
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(ttl + f" ({len(z)} wall defects)" if "meas" in k else ttl, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01).set_label("layer term, $V-V_\\mathrm{min,0.8}$ ($k_BT$)", fontsize=8)
    fig.text(0.5, 0.01, "Maps to the 0.98 level, layer term only (probe's own wall repulsion removed). Dotted verticals: |x| = 0.6 A (tip regions outside).", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(H, "fig_review_ellipse.png"), dpi=170)
    print("\n".join(L)); open(os.path.join(H, "review_ellipse.txt"), "w").write("\n".join(L) + "\n")
