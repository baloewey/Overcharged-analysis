#!/usr/bin/env python3
"""
Figures for the spatial-distribution / edge-site report (analysis_spatial/report/).

All positions in units of the cavity radius: r/R for the radial position (0 centre, 1 wall) and
r_pair/R for pair distances. Uses defect_spatial.py (maps, radial profiles, speeds, g(r | s)) and
edge_sites.py (pinned sites, classification).

  fig1_maps.png         time-averaged density of +1/2 and -1/2 defects, sites marked
  fig2_radial.png       radial number density n_±(r/R) and mean speed v_±(r/R)
  fig3_sitedist.png     distance of each +1/2 to the nearest site, per k (3 seeds pooled)
  fig4_bulkcharge.png   bulk charge and P(q_b = +1) vs k: site rule vs density cut (3 seeds)
  fig5_g_same.png       g_++ and g_-- vs r_pair/R, reference in the outer bulk / centre (R-scan pooled)
  fig6_g_pm.png         g_+- vs r_pair/R, same zones
  fig7_g_act.png        g_++ and g_-- (bulk reference) at activity 0.05 / 0.07 / 0.09 (same R/l_a)
"""
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from defect_spatial import analyse, load, ring_gr  # noqa: E402
from edge_sites import classify  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402

PLUS, MINUS, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
OUT = "analysis_spatial/report"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": INK2, "text.color": INK})


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def zone_g(runs, zones, keys=("++", "--", "+-")):
    """Pool g(r | s) over runs; zones are (name, r/R min, r/R max) for the reference defect."""
    O = {z[0]: {k: 0 for k in keys} for z in zones}
    E = {z[0]: {k: 0 for k in keys} for z in zones}
    x = None
    for run in runs:
        d = load(run, 0.2)
        R = d["R"]
        r = analyse(d, 0.25, 3.0, 4.2)
        ss, rc, obs, exp, _ = ring_gr(d["frames"], R, r["sc"], r["prof"], raw=True)
        rr = 1 - ss / R                                      # ring position in r/R
        xb = np.linspace(0, 2, 51)                            # common r_pair/R bins (0.04 R > 1 cell)
        xi = np.clip(np.digitize(rc / R, xb) - 1, 0, len(xb) - 2)
        for zn, a, b in zones:
            m = (rr >= a) & (rr < b)
            for k in keys:
                o = np.bincount(xi, obs[k][m].sum(0), len(xb) - 1)
                e = np.bincount(xi, exp[k][m].sum(0), len(xb) - 1)
                O[zn][k] = O[zn][k] + o
                E[zn][k] = E[zn][k] + e
        x = 0.5 * (xb[1:] + xb[:-1])
    return x, O, E


def plot_g(ax, x, o, e, col, label):
    g = np.where(e > 20, o / np.maximum(e, 1e-12), np.nan)
    err = np.where(o > 0, g / np.sqrt(np.maximum(o, 1)), np.nan)
    ax.plot(x, g, color=col, lw=1.9, label=label)
    ax.fill_between(x, g - err, g + err, color=col, alpha=0.2, lw=0)


def main():
    os.makedirs(OUT, exist_ok=True)
    ex = [("Rscan_R31.83_passive/k05", "passive, k = 5"), ("Rscan_R31.83_seeds/s02/k02", "active, k = 2"),
          ("Rscan_R31.83_seeds/s02/k05", "active, k = 5"), ("Rscan_R31.83_seeds/s02/k10", "active, k = 10")]
    D = [load(r, 0.2) for r, _ in ex]
    A = [analyse(d, 0.25, 3.0, 4.2) for d in D]
    C = [classify(d, 1.0, 0.02, 0.0, 4.2) for d in D]

    # fig 1: maps
    cm = {"+": LinearSegmentedColormap.from_list("p", ["#ffffff", "#cde2fb", "#6da7ec", "#2a78d6", "#104281"]),
          "-": LinearSegmentedColormap.from_list("m", ["#ffffff", "#fde0d2", "#f4a27f", "#eb6834", "#8f3413"])}
    fig, axs = plt.subplots(2, 4, figsize=(10, 5.3))
    th = np.linspace(0, 2 * np.pi, 400)
    for c, (d, r, cl, (_, lab)) in enumerate(zip(D, A, C, ex)):
        R = d["R"]
        for row, key in enumerate(("+", "-")):
            ax = axs[row, c]
            H = r["maps"][key].T
            vmax = max(H.max(), 1e-6)
            ext = [r["edges"][0] / R, r["edges"][-1] / R] * 2
            ax.imshow(np.where(H > 0, H, np.nan), origin="lower", extent=ext, cmap=cm[key],
                      norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), interpolation="nearest")
            ax.plot(np.cos(th), np.sin(th), color=INK2, lw=0.8)
            if key == "+" and cl["nsite"]:
                ax.scatter(cl["sites"][:, 0] / R, cl["sites"][:, 1] / R, s=26, facecolor="none", edgecolor=INK, lw=0.8)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if row == 0:
                ax.set_title(f"{lab}\n{cl['nsite']} sites", fontsize=9)
        axs[0, 0].set_ylabel("+½ defects", color=INK); axs[1, 0].set_ylabel("−½ defects", color=INK)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig1_maps.png", dpi=160); plt.close(fig)

    # fig 2: radial density and speed vs r/R
    fig, axs = plt.subplots(2, 4, figsize=(10.5, 5.2), sharex=True)
    for c, (d, r, (_, lab)) in enumerate(zip(D, A, ex)):
        R = d["R"]; x = 1 - r["sc"] / R
        ax = axs[0, c]; style(ax)
        ax.plot(x, r["prof"]["+"], color=PLUS, lw=1.7, label="+½")
        ax.plot(x, r["prof"]["-"], color=MINUS, lw=1.7, label="−½")
        ax.set_yscale("log"); ax.set_ylim(1e-5, 0.3); ax.set_title(lab)
        ax2 = axs[1, c]; style(ax2)
        ax2.plot(x, r["speed"]["+"], color=PLUS, lw=1.7)
        ax2.plot(x, r["speed"]["-"], color=MINUS, lw=1.7)
        ax2.set_ylim(0, 1.1); ax2.set_xlabel("radial position r/R")
        for a_ in (ax, ax2):
            a_.axvline(1 - 4.2 / R, color=INK2, lw=0.8, ls=":")
            a_.set_xlim(0, 1)
    axs[0, 0].set_ylabel("defects per cell² per frame"); axs[1, 0].set_ylabel("mean speed (cells per time unit)")
    axs[0, 0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig2_radial.png", dpi=160); plt.close(fig)

    # R-scan: all runs
    rs = sorted(glob.glob("Rscan_R31.83/k[01][0-9]")) + sorted(glob.glob("Rscan_R31.83_seeds/s0[23]/k[01][0-9]"))
    dens = {}
    for row in csv.DictReader(open("analysis_lunch/layer_rscan.csv")):
        dens[row["run"].rstrip("/")] = row
    rows, dists = [], {}
    for run in rs:
        d = load(run, 0.2)
        cl = classify(d, 1.0, 0.02, 0.0, 4.2)
        rows.append((d["k"], cl, dens.get(run)))
        dists.setdefault(d["k"], []).append((cl["dist"] / d["R"], cl["dsite"] / d["R"]))

    # fig 3: distance to nearest site, per k
    fig, ax = plt.subplots(figsize=(6.6, 3.7)); style(ax)
    ks = sorted(dists); cmb = plt.get_cmap("Blues")
    for i, k in enumerate(ks):
        dd = np.concatenate([x for x, _ in dists[k]]); dd = dd[np.isfinite(dd)]
        h, e = np.histogram(dd, np.linspace(0, 0.6, 97), density=True)
        col = cmb(0.35 + 0.65 * i / (len(ks) - 1))
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color=col, lw=1.5, label=f"k = {k}")
        ax.plot(np.mean([s for _, s in dists[k]]), 0.02, marker="v", color=col, ms=6, ls="none")
    ax.set_yscale("log"); ax.set_ylim(1e-2, None)
    ax.set_xlabel("distance of a +½ defect to the nearest site / R"); ax.set_ylabel("probability density")
    ax.legend(frameon=False, fontsize=7, ncol=3)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_sitedist.png", dpi=160); plt.close(fig)

    # fig 4: bulk charge and P1 vs k, site rule vs density cut
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax, key_s, key_d, yl in ((axs[0], "qb", "qb_density", "mean bulk charge ⟨q_b⟩"),
                                 (axs[1], "P1", "P1_density", "fraction of frames with q_b = +1")):
        style(ax)
        for off, lab, col, get in ((-0.12, "pinned-site rule", PLUS, lambda cl, dr: cl[key_s]),
                                   (0.12, "density-minimum cut", MINUS, lambda cl, dr: float(dr[key_d]))):
            m, s = [], []
            for k in ks:
                v = [get(cl, dr) for kk, cl, dr in rows if kk == k and dr is not None]
                m.append(np.mean(v)); s.append(np.std(v, ddof=1) / np.sqrt(len(v)))
            ax.errorbar(np.array(ks) + off, m, s, color=col, marker="o", ms=5, lw=0, elinewidth=1.4, capsize=0, label=lab)
        ax.set_xlabel("imposed winding k"); ax.set_ylabel(yl); ax.set_xticks(ks)
    axs[0].axhline(1, color=INK2, lw=0.7, ls="--"); axs[0].set_ylim(0.9, 1.2)
    axs[1].set_ylim(0.75, 1.0); axs[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig4_bulkcharge.png", dpi=160); plt.close(fig)

    # fig 5 / 6: pooled g vs r_pair/R, zones by the reference's r/R
    R0 = 31.830989
    zones = [("outer bulk", 0.5, 1 - 4.2 / R0), ("centre", 0.0, 0.5)]
    x, O, E = zone_g(rs, zones)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for ax, (zn, a, b) in zip(axs, zones):
        style(ax)
        plot_g(ax, x, O[zn]["++"], E[zn]["++"], PLUS, "+½ / +½")
        plot_g(ax, x, O[zn]["--"], E[zn]["--"], MINUS, "−½ / −½")
        ax.axhline(1, color=INK2, lw=0.7, ls="--")
        ax.set_title(f"reference defect in the {zn} ({a:.2f} < r/R < {b:.2f})"); ax.set_xlabel("pair distance / R")
        ax.set_xlim(0, 1.3)
    axs[0].set_ylabel("g"); axs[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig5_g_same.png", dpi=160); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for ax, (zn, a, b) in zip(axs, zones):
        style(ax)
        plot_g(ax, x, O[zn]["+-"], E[zn]["+-"], AQUA, "+½ / −½")
        ax.axhline(1, color=INK2, lw=0.7, ls="--")
        ax.set_title(f"reference +½ in the {zn} ({a:.2f} < r/R < {b:.2f})"); ax.set_xlabel("pair distance / R")
        ax.set_xlim(0, 1.3)
    axs[0].set_ylabel("g"); axs[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig6_g_pm.png", dpi=160); plt.close(fig)

    # fig 7: activity comparison, bulk reference (r/R < 1 - 4.2/R), pooled over k
    sets = [("ζ = 0.05", sorted(glob.glob("H3pilot/act0.05/k*")) + sorted(glob.glob("H3pilot_seeds/act0.05/s0*/k*/DONE")), "#9ec5f4"),
            ("ζ = 0.07", rs, "#2a78d6"), ("ζ = 0.09", sorted(glob.glob("H3pilot/act0.09/k*")), "#104281")]
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for lab, runs, col in sets:
        runs = [os.path.dirname(r) if r.endswith("DONE") else r for r in runs]
        Rr = load(runs[0], 0.2)["R"]
        x7, O7, E7 = zone_g(runs, [("bulk", 0.0, 1 - 4.2 / Rr)], keys=("++", "--"))
        for ax, key in zip(axs, ("++", "--")):
            plot_g(ax, x7, O7["bulk"][key], E7["bulk"][key], col, f"{lab} ({len(runs)} runs, R = {Rr:.1f})")
    for ax, t in zip(axs, ("+½ / +½", "−½ / −½")):
        style(ax); ax.axhline(1, color=INK2, lw=0.7, ls="--"); ax.set_xlim(0, 1.3)
        ax.set_title(f"{t}, reference in the bulk"); ax.set_xlabel("pair distance / R")
    axs[0].set_ylabel("g"); axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig7_g_act.png", dpi=160); plt.close(fig)
    print("written", OUT)


if __name__ == "__main__":
    main()
